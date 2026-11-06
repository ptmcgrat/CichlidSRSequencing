"""Per-run processing: download -> FASTQ to uBAM -> align -> mark duplicates -> binned depth -> keep only the depth track.

The alignment and duplicate-marking commands are the LAB PIPELINE'S OWN (helper_modules.seq_pipeline.align_command / markdup_command, installed by lab_repo_patch_v6): the same minimap2
preset, the same MergeBamAlignment options (including overlapping-read clipping, which changes depth) and the same MarkDuplicates settings as the mbuna samples, so the blocks found here
and the B genotype computed later on the mbuna see the same mapping behaviour.

Every stage writes its output under a .partial name and renames it when the command succeeded; a marker file per stage makes a re-run skip what is done. Large files are deleted as soon as
the next stage has consumed them. A SpaceLedger keeps the runs that are in flight from filling the disk."""
import json, os, shutil, threading, time
from . import ena, depth as D

STAGES = ['fetch', 'ubam', 'align', 'markdup', 'slice', 'depth', 'finish']


class StageError(RuntimeError):
    pass


def dir_size(path):
    tot = 0
    for r, _, fs in os.walk(path):
        for f in fs:
            try:
                tot += os.path.getsize(os.path.join(r, f))
            except OSError:
                pass
    return tot


class Executor:
    """Runs a command list; raises StageError with the tail of stderr. Tests replace this class with a fake that fabricates the outputs."""
    def run(self, cmd, label, cwd=None):
        import subprocess
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
        if r.returncode:
            raise StageError(f'{label} failed (exit {r.returncode}): {" ".join(map(str, cmd))[:300]}\n{(r.stderr or "")[-1500:]}')
        return r


class LabCommands:
    """The command builders of the lab pipeline. Imported lazily: a clear message if lab_repo_patch_v6 is not installed. The repo root (the folder that contains helper_modules/) is found
    by walking up from this package, as yhsd.common.Project does."""
    def __init__(self):
        from cichlid_common.common_light import ensure_lab_modules_on_path
        ensure_lab_modules_on_path()
        sp = None
        for name in ('align_pipeline.seq_pipeline', 'helper_modules.seq_pipeline'):             # module 03 after the move, the original place before it
            try:
                sp = __import__(name, fromlist=['x']); break
            except ImportError as e:
                err = e
        if sp is None:
            raise SystemExit('bchrom needs the lab alignment pipeline (seq_pipeline.py from lab_repo_patch_v6): install that patch first. (' + str(err) + ')')
        self.align, self.markdup, self.java_threads = sp.align_command, sp.markdup_command, getattr(sp, 'JAVA_THREADS', '')


def fastq_to_sam_command(r1, r2, out_bam, sample, tmp_dir, heap_gb=8, java_threads=''):
    return ['gatk', '--java-options', f'-Xmx{heap_gb}g {java_threads}'.strip(), 'FastqToSam', '--FASTQ', r1, '--FASTQ2', r2, '--OUTPUT', out_bam, '--SAMPLE_NAME', sample,
            '--READ_GROUP_NAME', sample, '--LIBRARY_NAME', sample, '--PLATFORM', 'ILLUMINA', '--TMP_DIR', tmp_dir]


class SpaceLedger:
    """Admission control for the disk. A run may start when the peak space of all runs in flight (each an estimate) fits into what is free now plus what they already occupy, minus a margin.
    free_fn(path) and usage_fn(path) are injectable for tests."""
    def __init__(self, base, min_free_gb, free_fn=None, usage_fn=None, sleep=time.sleep, poll=30.0, timeout=24 * 3600, log=print):
        self.base, self.min_free = base, min_free_gb * 1e9
        self.free_fn = free_fn or (lambda p: shutil.disk_usage(p).free)
        self.usage_fn = usage_fn or dir_size
        self.sleep, self.poll, self.timeout, self.log = sleep, poll, timeout, log
        self.inflight, self.lock = {}, threading.Lock()

    def try_admit(self, run, need):
        with self.lock:
            used = sum(self.usage_fn(p) for p, _ in self.inflight.values())
            budget = self.free_fn(self.base) + used - self.min_free
            if sum(n for _, n in self.inflight.values()) + need <= budget:
                self.inflight[run] = (os.path.join(self.base, run), need)
                return True
            return False

    def acquire(self, run, need):
        t0 = time.time()
        while not self.try_admit(run, need):
            if time.time() - t0 > self.timeout:
                raise StageError(f'{run}: waited {self.timeout / 3600:.0f} h for {need / 1e9:.0f} GB of disk (min free {self.min_free / 1e9:.0f} GB)')
            self.log(f'  [{run}] waiting for disk: need {need / 1e9:.0f} GB, free {self.free_fn(self.base) / 1e9:.0f} GB')
            self.sleep(self.poll)

    def release(self, run):
        with self.lock:
            self.inflight.pop(run, None)


def sources_of(cfg):
    """The sources tried in turn for every run: ena = ENA's files (https, md5-verified); sra = NCBI's copy on AWS through module 02 (aria2c, vdb-validate, fasterq-dump, read-count and quality checks);
    auto = SRA first, ENA as the fallback. ENA is ~0.3 MB/s from most networks and SRA ~100 MB/s, hence the default."""
    return {'ena': ['ena'], 'sra': ['sra'], 'auto': ['sra', 'ena']}[cfg['download'].get('source', 'ena')]


class RunPaths:
    def __init__(self, work, run):
        self.run, self.dir = run, os.path.join(work, run) + '/'
        self.fastq, self.tmp, self.stages = self.dir + 'fastq/', self.dir + 'tmp/', self.dir + 'stages/'
        self.ubam, self.sorted_bam, self.bam = self.dir + f'{run}.u.bam', self.dir + f'{run}.sorted.bam', self.dir + f'{run}.all.bam'
        self.bai = self.dir + f'{run}.all.bai'
        self.metrics, self.depth_prefix = self.dir + f'{run}.dup_metrics.txt', self.dir + f'{run}.mosdepth'
        self.slice_bam = os.path.join(work, 'slices', run + '.slice.bam')
        self.track = os.path.join(work, 'depth', run)            # <track>.npz / .json: kept across runs until `bchrom clean`
        self.summary = os.path.join(work, 'depth', run + '.summary.json')


class RunProcessor:
    def __init__(self, spec, rec, work, cfg, genome, executor=None, cmds=None, ledger=None, log=print, downloader=None, keep_bams=False, keep_fastq=False, sra_fetcher=None, sleep=time.sleep):
        self.spec, self.rec, self.cfg, self.g = spec, rec, cfg, genome            # genome: dict(reference, index, lengths_total)
        self.p = RunPaths(work, spec['id'])
        self.ex, self.cmds, self.ledger, self.log = executor or Executor(), cmds, ledger, log
        self.downloader = downloader or (lambda r, d: ena.download_file(r, d, scheme=cfg['download'].get('scheme', 'https'), attempts=cfg['download']['attempts'], backoff=cfg['download']['backoff_s'], log=log))
        self.keep_bams, self.keep_fastq = keep_bams, keep_fastq
        self.sra_fetcher, self.sleep = sra_fetcher, sleep
        self.id = spec['id']

    def say(self, msg):
        self.log(f'  [{time.strftime("%H:%M:%S")}] [{self.id}] {msg}')

    # ---- checkpoints
    def done(self, stage):
        return os.path.exists(self.p.stages + stage + '.done')

    def mark(self, stage):
        os.makedirs(self.p.stages, exist_ok=True)
        json.dump(dict(stage=stage, time=time.time()), open(self.p.stages + stage + '.done.tmp', 'w'))
        os.replace(self.p.stages + stage + '.done.tmp', self.p.stages + stage + '.done')

    def stage(self, name, fn):
        if self.done(name):
            self.say(f'{name}: already done, skipped'); return
        t0 = time.time(); self.say(f'{name}: start')
        fn()
        self.mark(name)
        self.say(f'{name}: done in {(time.time() - t0) / 60:.1f} min')

    def finished(self):
        if self.g.get('slice_regions') and not os.path.exists(self.p.slice_bam):          # slices were asked for and this run has none yet: not finished
            return False
        return os.path.exists(self.p.summary) and os.path.exists(self.p.track + '.npz')

    # ---- stages
    def _fetch_ena(self):
        r1, r2 = ena.paired_files(self.rec)
        self.files = [self.downloader(r, self.p.fastq) for r in (r1, r2)]

    def _fetch_sra(self):
        """NCBI's copy through module 02's SRA source. Transient errors are retried with a back-off; `source_unavailable` (object missing, counts or qualities differ from ENA) is final for this source."""
        from download_reads.sources import sra as SF
        from download_reads.sources.ena import EnaError
        d = self.cfg['download']; fetcher = self.sra_fetcher or SF.fetch_sra
        r1, _ = ena.paired_files(self.rec)
        spec = dict(run=self.id, ena_read_count=self.rec.get('read_count'), paired=str(self.rec.get('layout') or 'PAIRED').upper() == 'PAIRED', ena_url_r1=ena.full_url(r1['url'], d.get('scheme', 'https')),
                    instrument='', fq_threads=d.get('sra_fq_threads', 8), sra_connections=d.get('sra_connections', 8))
        attempts = int(d.get('attempts', 4))
        for k in range(attempts):
            try:
                paths, stats = fetcher(spec, self.p.fastq)
            except EnaError as e:
                if e.kind in ('transient', 'throttle') and k + 1 < attempts:
                    self.say(f'SRA {e.kind} error ({str(e)[:160]}); retrying'); self.sleep(float(d.get('backoff_s', 30)) * (k + 1)); continue
                raise ena.DownloadError(f'SRA ({e.kind}): {e}')
            self.say(f"SRA: {stats.get('mb_per_s', '?')} MB/s, quality check: {stats.get('quality_check')}, read count: {stats.get('count_check') or 'ok'}")
            self.files = paths
            return

    def _fetch(self):
        os.makedirs(self.p.fastq, exist_ok=True)
        srcs, errors = sources_of(self.cfg), []
        for k, src in enumerate(srcs):
            try:
                (self._fetch_sra if src == 'sra' else self._fetch_ena)()
                errors = None
                break
            except ena.DownloadError as e:
                errors.append(f'{src}: {str(e)[:300]}')
                if k + 1 < len(srcs):                                                      # the next source starts clean; with a single source a partial download is kept (it resumes)
                    self.say(f'{src} failed ({str(e)[:200]}); trying {srcs[k + 1]}')
                    shutil.rmtree(self.p.fastq, ignore_errors=True); os.makedirs(self.p.fastq, exist_ok=True)
        if errors:
            raise ena.DownloadError(' | '.join(errors))
        un = ena.unpaired_files(self.rec)
        if un:
            self.say(f'ignoring {len(un)} unpaired FASTQ file(s): {[f["name"] for f in un]}')

    def _disk_factor(self):
        """Peak disk of a run as a multiple of its ENA FASTQ size: the SRA route holds the .sra file, uncompressed FASTQ (about 4x the gz) and fasterq-dump's scratch at the same time."""
        f = float(self.cfg['resources']['disk_factor'])
        return max(f, 9.0) if 'sra' in sources_of(self.cfg) else f

    def _fastq_pair(self):
        """The two FASTQ files of the run: the ENA names (.fastq.gz) or, after the SRA Toolkit fallback, <run>_1.fastq / <run>_2.fastq."""
        ena_names = [f['name'] for f in ena.paired_files(self.rec)]
        for names in (ena_names, [f'{self.id}_1.fastq', f'{self.id}_2.fastq']):
            if all(os.path.exists(self.p.fastq + n) for n in names):
                return [self.p.fastq + n for n in names]
        raise StageError(f'{self.id}: no FASTQ pair found in {self.p.fastq}')

    def _ubam(self):
        r1, r2 = self._fastq_pair()
        part = self.p.ubam + '.partial.bam'
        os.makedirs(self.p.tmp, exist_ok=True)
        if os.path.exists(part):
            os.remove(part)
        self.ex.run(fastq_to_sam_command(r1, r2, part, self.id, self.p.tmp, self.cfg['resources']['fastqtosam_heap_gb'], self.cmds.java_threads), 'FastqToSam')
        os.replace(part, self.p.ubam)
        if not self.keep_fastq:
            shutil.rmtree(self.p.fastq, ignore_errors=True)

    def _align(self):
        res = self.cfg['resources']
        part = self.p.sorted_bam + '.partial.bam'
        if os.path.exists(part):
            os.remove(part)
        os.makedirs(self.p.tmp, exist_ok=True)
        logs = [f'{self.p.tmp}Alignment_errors_{k}.txt' for k in ('StF', 'MM', 'Merge')]
        cmd = self.cmds.align(self.p.ubam, self.g['index'], self.g['reference'], self.p.tmp, part, res['align_threads'], res['stf_heap_gb'], res['mba_heap_gb'], logs)
        self.ex.run(cmd, 'align (SamToFastq | minimap2 | MergeBamAlignment)')
        self.ex.run(['samtools', 'quickcheck', '-v', part], 'quickcheck aligned bam')
        os.replace(part, self.p.sorted_bam)

    def _markdup(self):
        res = self.cfg['resources']
        self.ex.run(self.cmds.markdup(self.p.sorted_bam, self.p.bam, self.p.metrics, self.p.tmp, res['markdup_heap_gb']), 'MarkDuplicates')
        self.ex.run(['samtools', 'quickcheck', '-v', self.p.bam], 'quickcheck final bam')
        if os.path.exists(self.p.bai) and not os.path.exists(self.p.bam + '.bai'):      # Picard writes <name>.bai; make <name>.bam.bai visible too
            os.symlink(os.path.basename(self.p.bai), self.p.bam + '.bai')
        for f in (self.p.sorted_bam, self.p.ubam):
            if os.path.exists(f):
                os.remove(f)

    def _slice(self):
        """Small BAM of the configured regions (config slices.regions), cut from the duplicate-marked BAM before it is deleted."""
        regs = self.g.get('slice_regions') or []
        if not regs:
            return
        os.makedirs(os.path.dirname(self.p.slice_bam), exist_ok=True)
        part = self.p.slice_bam + '.partial.bam'
        for f in (part, part + '.bai'):
            if os.path.exists(f):
                os.remove(f)
        self.ex.run(['samtools', 'view', '-b', '-o', part, self.p.bam] + regs, 'slice')
        self.ex.run(['samtools', 'index', part], 'index slice')
        os.replace(part, self.p.slice_bam); os.replace(part + '.bai', self.p.slice_bam + '.bai')

    def _depth(self):
        b = self.cfg['bin_size']
        for f in (self.p.depth_prefix + '.regions.bed.gz',):
            if os.path.exists(f):
                os.remove(f)
        self.ex.run(D.mosdepth_command(self.p.depth_prefix, self.p.bam, b, self.cfg['resources']['mosdepth_threads']), 'mosdepth')
        arrays, lengths = D.parse_regions(self.p.depth_prefix + '.regions.bed.gz', b)
        stats = D.depth_stats(arrays, b, self.cfg['scale_min_chrom_bp'])
        os.makedirs(os.path.dirname(self.p.track), exist_ok=True)
        extra = dict(run=self.id, role=self.spec['role'], label=self.spec.get('label'), baseline=self.spec.get('baseline'), reference=self.g['reference'])
        D.save_track(self.p.track, arrays, lengths, stats, b, extra)
        self.stats, self.lengths = stats, lengths

    def _finish(self):
        """Summary with the checks that catch a wrong accession or a failed mapping, then delete everything large."""
        st = getattr(self, 'stats', None) or json.load(open(self.p.track + '.json'))['stats']
        lengths = getattr(self, 'lengths', None) or json.load(open(self.p.track + '.json'))['lengths']
        genome_len = sum(lengths.values())
        expected = self.rec['base_count'] / genome_len if self.rec.get('base_count') else None
        ratio = st['mean_all_sequences'] / expected if expected else None
        flags = []
        if self.spec.get('bases_gb') and self.rec.get('base_count'):
            if abs(self.rec['base_count'] / 1e9 - self.spec['bases_gb']) > 0.01 * self.spec['bases_gb'] + 0.01:
                flags.append(f"ENA base count {self.rec['base_count'] / 1e9:.2f} Gb differs from the configured {self.spec['bases_gb']} Gb: is this the right run?")
        lo, hi = self.cfg['qc']['depth_ratio_range']
        if ratio is not None and not (lo <= ratio <= hi):
            flags.append(f'mean depth is {ratio:.2f} x the depth expected from the number of bases (expected {lo}-{hi}): check the mapping')
        if st['zero_fraction'] > self.cfg['qc']['max_zero_fraction']:
            flags.append(f"{st['zero_fraction']:.0%} of the bins on the placed chromosomes have no coverage")
        summ = dict(run=self.id, role=self.spec['role'], label=self.spec.get('label'), alias=self.rec.get('alias'), baseline=self.spec.get('baseline'), bases=self.rec.get('base_count'),
                    fastq_md5={f['name']: f['md5'] for f in self.rec['files']}, observed_mean_depth=st['mean_all_sequences'], expected_mean_depth=expected, depth_ratio=ratio,
                    median_depth=st['median'], mean_depth_placed=st['mean'], zero_fraction=st['zero_fraction'], flags=flags, bin_size=self.cfg['bin_size'], reference=self.g['reference'],
                    dup_metrics=self._dup_line(), finished=time.strftime('%Y-%m-%dT%H:%M:%S'))
        json.dump(summ, open(self.p.summary + '.tmp', 'w'), indent=1); os.replace(self.p.summary + '.tmp', self.p.summary)
        for fl in flags:
            self.say('WARNING ' + fl)

    def _dup_line(self):
        try:
            lines = [l for l in open(self.p.metrics) if l.strip() and not l.startswith('#')]
            return dict(zip(lines[0].rstrip('\n').split('\t'), lines[1].rstrip('\n').split('\t')))
        except Exception:
            return None

    def cleanup(self):
        if self.keep_bams:
            return
        shutil.rmtree(self.p.dir, ignore_errors=True)

    def process(self):
        if self.finished():
            self.say('complete (depth track and summary present): skipped'); return json.load(open(self.p.summary))
        try:
            if self.ledger:                                                        # admission once per run, also when it resumes half way
                self.ledger.acquire(self.id, self._disk_factor() * sum(f['bytes'] for f in self.rec['files']))
            for name, fn in (('fetch', self._fetch), ('ubam', self._ubam), ('align', self._align), ('markdup', self._markdup), ('slice', self._slice), ('depth', self._depth), ('finish', self._finish)):
                self.stage(name, fn)
            self.cleanup()
            return json.load(open(self.p.summary))
        finally:
            if self.ledger:
                self.ledger.release(self.id)
