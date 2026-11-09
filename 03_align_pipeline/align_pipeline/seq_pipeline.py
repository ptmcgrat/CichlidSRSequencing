"""A pipelined replacement for the stage-by-stage alignFastQReads.py (AlignmentWorker).

What the old design did: every stage was a barrier over ALL samples (download everything, then align everything one sample at a time, then mark duplicates for everything, ...).
Nothing was uploaded or written to the database until the very end, a failure late on lost the whole batch, `max_processes` equalled the core count although most jobs use more than
one core, and a crash in `SamToFastq` or `minimap2` could go unnoticed because only the last process of the pipe was checked.

What this does instead, per sample, independently of the other samples:
    download -> align (+merge) -> mark duplicates -> [split | stats | slices | panel | count uBAM reads] -> qc -> verified upload -> database row -> cleanup
  * a shared CORE and RAM budget (ResourceBudget): every command states what it needs and waits for it, so the machine is never oversubscribed;
  * several samples are in flight at once; downloads run only a few samples ahead of the aligner (disk window), not for the whole batch;
  * every stage leaves a checkpoint; a rerun continues where a sample stopped, a failure of one sample never stops the others;
  * every command is checked, including EVERY process of the alignment pipe (bash -o pipefail);
  * the database row is written as soon as the sample's data is verified on Dropbox;
  * no GVCFs (not needed for genotyping known variants from BAMs); only the chimeric split is made (decision 41: the panel uses no other); region slices (LG2, LG10, LG13) are cut from the final BAM;
  * the panel stage (panel.py, decision 40): module 04's measurements run on the final BAM while it is on disk; a failure of it never fails the sample;
  * every command that talks to the network has a time limit (v6 had none: a hung rclone would have hung the sample, as one did for 88 days on utaka);
  * downloads and uploads overlap with compute: a continuous, SPACE-AWARE download manager (disk_ledger.py, download_manager.py) runs ahead of the aligner but never ahead of the disk;
    finished samples are uploaded by a separate pool, so a slow upload never holds a compute slot; the uBAMs are deleted as soon as the duplicate-marked BAM is verified;
  * later stages are served before new alignments are started (priorities in the core budget), so finished work drains and uploads start early.
The alignment, duplicate-marking and read-routing recipes are unchanged (same GATK/minimap2 commands and options as the old worker)."""
import collections
import itertools
import json
import os
import re
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace

from . import __version__
from . import alignment_db as adb
from . import bam_slices as sl
from . import bam_stats as bs
from . import panel as pn
from .disk_ledger import DiskLedger, dir_size
from .download_manager import DownloadManager
from .procs import bam_check

HERE = os.path.dirname(os.path.abspath(__file__))
SPLIT_SCRIPT = os.path.join(HERE, 'bam_split.py')                  # run by path: it needs nothing but pysam
# A JVM sizes its garbage-collector pool to every core of the machine (92 threads per HaplotypeCaller JVM were seen on a 96-core host); each of our jobs needs ~1-2 cores.
JAVA_THREADS = '-XX:ParallelGCThreads=2'
SPLIT_TYPES = ('chimeric',)                                        # decision 41: the panel uses no other split


def short_reason(message, max_len=220):
    """One line for the sheet and the log out of a StageError message: 'exit code 3' (not the long command that follows it) and the last two lines the tool wrote."""
    lines = [l for l in message.splitlines() if l.strip()]
    if not lines:
        return ''
    m = re.match(r'(exit code -?\d+):', lines[0])
    head = m.group(1) if m else lines[0][:120]
    tail = ' / '.join(l.strip() for l in lines[1:][-2:])
    return (head + (': ' + tail if tail else ''))[:max_len]


def sh_join(cmd):
    import shlex
    return ' '.join(shlex.quote(str(c)) for c in cmd)


class StageError(Exception):
    def __init__(self, sample, stage, message):
        super().__init__(f'{sample} / {stage}: {message}')
        self.sample, self.stage, self.message = sample, stage, message


class StopRequested(Exception):
    pass


class CommandTimeout(Exception):
    pass


class PanelFailed(Exception):
    """The panel measurement of a sample failed. Not a StageError: the sample goes on (decision 40)."""
    def __init__(self, sample, message):
        super().__init__(f'{sample} / panel: {message}')
        self.sample, self.message = sample, message


# ============================================================================================ resources
class ResourceBudget:
    """Cores and RAM (GB) shared by every running command. The request at the head of the queue is served first; the head is the request with the LOWEST priority number, and among equal
    priorities the oldest. Later stages have lower numbers than alignment, so finished work drains (and uploads start) before new alignments take cores, and a big alignment request
    is still never starved by a stream of requests of its own priority."""

    def __init__(self, cores, ram_gb, stop=None):
        self.cores, self.ram = int(cores), float(ram_gb)
        self.free_cores, self.free_ram = self.cores, self.ram
        self.cv, self.queue, self.stop = threading.Condition(), [], stop or threading.Event()
        self._seq = itertools.count()
        self.peak_cores = self.peak_ram = 0

    def acquire(self, cores, ram_gb, priority=1):
        cores, ram_gb = min(int(cores), self.cores), min(float(ram_gb), self.ram)
        ticket = SimpleNamespace(prio=priority, seq=next(self._seq))
        with self.cv:
            self.queue.append(ticket)
            try:
                while True:
                    head = min(self.queue, key=lambda x: (x.prio, x.seq))
                    if head is ticket and self.free_cores >= cores and self.free_ram >= ram_gb:
                        break
                    if self.stop.is_set():
                        raise StopRequested()
                    self.cv.wait(timeout=1.0)
            except BaseException:
                self.queue.remove(ticket)
                self.cv.notify_all()
                raise
            self.queue.remove(ticket)
            self.free_cores -= cores
            self.free_ram -= ram_gb
            self.peak_cores = max(self.peak_cores, self.cores - self.free_cores)
            self.peak_ram = max(self.peak_ram, self.ram - self.free_ram)
            self.cv.notify_all()
        return cores, ram_gb

    def release(self, held):
        with self.cv:
            self.free_cores += held[0]
            self.free_ram += held[1]
            self.cv.notify_all()

    def in_use(self):
        with self.cv:
            return self.cores - self.free_cores


def kill_group(p, grace=5):
    try:
        os.killpg(p.pid, signal.SIGTERM)
    except Exception:
        pass
    try:
        p.wait(timeout=grace)
    except Exception:
        try:
            os.killpg(p.pid, signal.SIGKILL)
        except Exception:
            pass


def real_exec(cmd, stdout_path, stderr_path, produces, stop, timeout_s=None):
    out = open(stdout_path, 'w') if stdout_path else subprocess.DEVNULL
    err = open(stderr_path, 'w') if stderr_path else subprocess.DEVNULL
    try:
        p = subprocess.Popen(cmd, stdout=out, stderr=err, start_new_session=True)      # own process group: everything it starts can be stopped with it
        started = time.time()
        while True:
            try:
                return p.wait(timeout=1.0)
            except subprocess.TimeoutExpired:
                if stop.is_set():
                    kill_group(p)
                    raise StopRequested()
                if timeout_s and time.time() - started > timeout_s:
                    kill_group(p)
                    raise CommandTimeout(timeout_s)
    finally:
        for f in (out, err):
            if f is not subprocess.DEVNULL:
                f.close()


class Runner:
    def __init__(self, budget, nice=0, stop=None, exec_fn=None, log=None):
        self.budget, self.nice, self.stop = budget, nice, stop or budget.stop
        self.exec_fn, self.log = exec_fn or real_exec, log or (lambda s, m: None)

    @staticmethod
    def tail(path, n=12):
        try:
            with open(path, errors='replace') as f:
                return ''.join(f.readlines()[-n:]).strip()
        except OSError:
            return ''

    def run(self, cmd, cores, ram_gb, label, sample='', stdout=None, stderr=None, produces=(), extra_logs=(), priority=1, timeout_s=None):
        full = (['nice', '-n', str(self.nice)] if self.nice else []) + list(cmd)
        held = self.budget.acquire(cores, ram_gb, priority)
        try:
            try:
                rc = self.exec_fn(full, stdout, stderr, list(produces), self.stop, **({'timeout_s': timeout_s} if timeout_s else {}))
            except FileNotFoundError as e:
                raise StageError(sample, label, f'program not found: {e}')
            except CommandTimeout:
                raise StageError(sample, label, f'timed out after {timeout_s:.0f} s and was killed: {sh_join(cmd)[:300]}')
        finally:
            self.budget.release(held)
        if rc != 0:
            tails = '\n'.join(t for t in ([self.tail(stderr)] if stderr else []) + [self.tail(x) for x in extra_logs] if t)
            raise StageError(sample, label, f'exit code {rc}: {sh_join(cmd)[:300]}\n{tails}')
        missing = [p for p in produces if not os.path.exists(p)]
        if missing:
            raise StageError(sample, label, f'finished but did not create: {missing}')


# ============================================================================================ paths and checkpoints
class SamplePaths:
    def __init__(self, sample, ubams, bam_dir, tmp_dir, scratch_dir=None, need_bytes=None):
        self.sample, self.ubams = sample, list(ubams)
        self.bam_dir, self.tmp_dir = bam_dir.rstrip('/') + '/', tmp_dir.rstrip('/') + '/'
        self.scratch_dir = (scratch_dir or tmp_dir).rstrip('/') + '/'                 # GATK/Picard scratch files (--TMP_DIR); a fast disk if --tmp-dir is given
        self.need_bytes = float(need_bytes) if need_bytes else 20e9 * max(1, len(self.ubams))     # size of the uBAMs, for the disk ledger (20 GB per file if unknown)
        self.all_bam = f'{self.bam_dir}{sample}.all.bam'
        self.all_bai = f'{self.bam_dir}{sample}.all.bai'
        self.metrics = self.all_bam + '.duplication_metrics.txt'
        self.sorted_bam = f'{self.tmp_dir}{sample}.sorted.bam'
        self.split_bams = {t: f'{self.bam_dir}{sample}.{t}.bam' for t in SPLIT_TYPES}
        self.split_prefix = f'{self.bam_dir}{sample}'
        self.slice_dir = f'{self.bam_dir}slices/'
        self.stages_dir = f'{self.tmp_dir}stages/'
        self.stats_txt, self.idxstats_txt = self.all_bam + '.stats.txt', self.all_bam + '.idxstats.txt'
        self.ubam_reads_txt = f'{self.tmp_dir}ubam_reads.txt'
        self.qc_json = f'{self.bam_dir}{sample}.qc.json'
        self.panel_json = f'{self.bam_dir}{sample}.panel.json'                    # next to the BAM, uploaded with the sample's folder
        self.panel_status, self.panel_error, self.panel_note = 'off', '', ''
        self.holds_slot = False

    @staticmethod
    def bai(bam):
        return bam[:-4] + '.bai'

    @classmethod
    def from_fm(cls, fm, sample, scratch_base=None, need_bytes=None):
        """Snapshot the paths FileManager would use. createSampleFiles mutates shared state, so it is only called here, one sample at a time, before any thread starts."""
        fm.createSampleFiles(sample)
        sp = cls(sample, fm.localRawBamFiles, fm.localSampleBamDir, fm.localSampleTempDir, scratch_dir=(scratch_base.rstrip('/') + '/' + sample) if scratch_base else None, need_bytes=need_bytes)
        for mine, theirs in ((sp.all_bam, fm.localBamFile), (sp.sorted_bam, fm.localTempSortedBamFile), (sp.split_bams['chimeric'], fm.localChimericBamFile)):
            if mine != theirs:
                raise RuntimeError(f'file naming differs from FileManager ({mine} vs {theirs}); refusing to run')
        return sp


def stage_marker(sp, stage):
    return f'{sp.stages_dir}{stage}.json'


# Stages whose outputs are deliberately deleted by the stage that consumes them: once the consumer is done, the earlier stage still counts as done
# (without this a resume would repeat the alignment, the most expensive stage, because the sorted BAM is removed after duplicate marking).
CONSUMED_BY = {'align': 'markdup'}


def _marker_ok(sp, stage, need_outputs):
    try:
        with open(stage_marker(sp, stage)) as f:
            rec = json.load(f)
    except (OSError, ValueError):
        return False
    return (not need_outputs) or all(os.path.exists(p) and os.path.getsize(p) == size for p, size in rec.get('outputs', {}).items())


def stage_done(sp, stage):
    """A stage is done when its marker exists and every output it recorded still exists with the recorded size (or a later stage has consumed them and is itself done)."""
    if _marker_ok(sp, stage, True):
        return True
    consumer = CONSUMED_BY.get(stage)
    return bool(consumer) and _marker_ok(sp, stage, False) and stage_done(sp, consumer)


def mark_stage(sp, stage, outputs):
    os.makedirs(sp.stages_dir, exist_ok=True)
    tmp = stage_marker(sp, stage) + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(dict(stage=stage, time=time.time(), outputs={p: os.path.getsize(p) for p in outputs}), f)
    os.replace(tmp, stage_marker(sp, stage))


# ============================================================================================ command builders (pure, unit tested)
def align_command(ubam, index, reference, tmp_dir, out_bam, threads, stf_heap_gb, mba_heap_gb, logs):
    """The old pipeline, unchanged, as ONE bash pipe with pipefail so that a failure of ANY of the three processes fails the command (the old code only looked at the last one:
    MergeBamAlignment writes reads that never got an alignment as unmapped records, so a minimap2 killed half way could have looked like a success)."""
    jo = lambda heap: f'-Xmx{heap}g {JAVA_THREADS}'
    stf = ['gatk', '--java-options', jo(stf_heap_gb), 'SamToFastq', '-I', ubam, '--FASTQ', '/dev/stdout', '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2',
           '--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', tmp_dir]
    mm = ['minimap2', '-x', 'sr', '-a', '-t', str(threads), index, '-']
    mba = ['gatk', '--java-options', jo(mba_heap_gb), 'MergeBamAlignment', '-R', reference, '--UNMAPPED_BAM', ubam, '--ALIGNED_BAM', '/dev/stdin', '-O', out_bam,
           '--ADD_MATE_CIGAR', 'true', '--CLIP_ADAPTERS', 'false', '--CLIP_OVERLAPPING_READS', 'true', '--INCLUDE_SECONDARY_ALIGNMENTS', 'true',
           '--MAX_INSERTIONS_OR_DELETIONS', '-1', '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant', '--ATTRIBUTES_TO_RETAIN', 'XS', '--TMP_DIR', tmp_dir]
    pipe = ' | '.join(f'{sh_join(c)} 2> {sh_join([l])}' for c, l in zip((stf, mm, mba), logs))
    return ['bash', '-o', 'pipefail', '-c', pipe]


def merge_command(inputs, out_bam, tmp_dir, heap_gb):
    cmd = ['gatk', '--java-options', f'-Xmx{heap_gb}g {JAVA_THREADS}', 'MergeSamFiles', '--TMP_DIR', tmp_dir]
    for i in inputs:
        cmd += ['-I', i]
    return cmd + ['-O', out_bam]


def markdup_command(sorted_bam, out_bam, metrics, tmp_dir, heap_gb):
    return ['gatk', '--java-options', f'-Xmx{heap_gb}g {JAVA_THREADS}', 'MarkDuplicates', '-I', sorted_bam, '-O', out_bam, '-M', metrics, '--TMP_DIR', tmp_dir,
            '--SORTING_COLLECTION_SIZE_RATIO', '.2', '--CREATE_INDEX']


# ============================================================================================ configuration
class Config:
    PRIORITY = {'align': 2, 'merge': 1, 'markdup': 1, 'split': 1, 'stats': 1, 'count': 1, 'slice': 1, 'panel': 1, 'index': 0, 'small': 0, 'upload': 0}      # lower = served first

    def __init__(self, max_cores=None, max_ram_gb=None, samples_in_flight=6, align_threads=16, prefetch=2, max_downloads=4, uploads=2, keep_local=False,
                 slices_file=None, nice=5, verify_counts=True, upload=True, align_stf_heap_gb=4, align_mba_heap_gb=8, markdup_heap_gb=16, merge_heap_gb=8, python=None,
                 max_pending_uploads=None, min_free_gb=200, disk_reserve_factor=5.5, scratch_base=None, status_s=600, check_workers=2, download_attempts=4, download_backoff=30.0,
                 panel=None, pipeline_version=None, upload_min_mb_s=5.0, small_timeout_s=3600):
        self.max_cores = int(max_cores or max(1, (os.cpu_count() or 8) - 8))
        self.max_ram_gb = float(max_ram_gb or self.default_ram())
        self.samples_in_flight, self.align_threads, self.prefetch, self.max_downloads, self.uploads = samples_in_flight, align_threads, prefetch, max_downloads, uploads
        self.keep_local, self.slices_file, self.nice, self.verify_counts, self.upload = keep_local, slices_file, nice, verify_counts, upload
        self.align_stf_heap_gb, self.align_mba_heap_gb, self.markdup_heap_gb, self.merge_heap_gb = align_stf_heap_gb, align_mba_heap_gb, markdup_heap_gb, merge_heap_gb
        self.python = python or sys.executable
        self.max_pending_uploads = int(max_pending_uploads or (uploads + 2))          # finished samples that may wait locally for an upload slot; more and the compute side waits
        self.min_free_gb, self.disk_reserve_factor, self.scratch_base, self.status_s = min_free_gb, disk_reserve_factor, scratch_base, status_s
        self.check_workers, self.download_attempts, self.download_backoff = check_workers, download_attempts, download_backoff
        self.panel = panel                                                              # a panel.PanelConfig, or None: no panel stage
        self.pipeline_version = pipeline_version or f'align_pipeline {__version__}'
        self.upload_min_mb_s, self.small_timeout_s = upload_min_mb_s, small_timeout_s    # the slowest upload speed that is still believable; the limit of the quick commands (quickcheck, idxstats, rclone check)

    @staticmethod
    def default_ram():
        try:
            with open('/proc/meminfo') as f:
                for line in f:
                    if line.startswith('MemAvailable:'):
                        return 0.7 * int(line.split()[1]) / 1e6                          # kB -> GB
        except (OSError, ValueError):
            pass
        try:
            import psutil
            return 0.7 * psutil.virtual_memory().available / 1e9
        except Exception:
            return 64.0

    @property
    def disk_slots(self):
        return self.samples_in_flight + self.prefetch

    def cost(self, what):
        """(cores, RAM GB) a command of this kind needs."""
        return {'align': (self.align_threads + 2, 14 + self.align_stf_heap_gb + self.align_mba_heap_gb), 'merge': (2, self.merge_heap_gb + 2), 'markdup': (2, self.markdup_heap_gb + 2),
                'split': (3, 3), 'index': (2, 1), 'stats': (4, 2), 'count': (4, 1), 'slice': (2, 1), 'panel': (2, 2), 'small': (1, 0.5), 'upload': (1, 1)}[what]

    def priority(self, what):
        return self.PRIORITY.get(what, 1)

    def timeout_for(self, what, nbytes=0):
        """Seconds after which a command of this kind is killed, or None for no limit. Alignment, duplicate marking and the like legitimately run for hours and have none; a transfer must make progress
        at `upload_min_mb_s` at least; the quick commands and the panel are given generous fixed limits."""
        if what == 'upload':
            return max(1800.0, nbytes / (self.upload_min_mb_s * 1e6))
        if what == 'small':
            return float(self.small_timeout_s)
        if what == 'panel':
            return float(self.panel.timeout_s) if self.panel else None
        return None


# ============================================================================================ the pipeline
class SeqPipeline:
    def __init__(self, cfg, runner, reference, minimap_index, reference_size, download_start_fn, upload_cmds_fn, db_add_fn, run_ids=None, read_lengths=None, versions=None,
                 genome_version='', log=None, check_fn=bam_check, slices=None, clock=time.time, ledger=None, download_kwargs=None, sample_meta=None, pipeline_version=None):
        self.cfg, self.runner, self.reference, self.index, self.ref_size = cfg, runner, reference, minimap_index, reference_size
        self.download_start_fn, self.upload_cmds_fn, self.db_add_fn = download_start_fn, upload_cmds_fn, db_add_fn
        self.run_ids, self.read_lengths, self.versions, self.genome_version = run_ids or {}, read_lengths or {}, versions or {}, genome_version
        self.check_fn, self.slices, self.clock = check_fn, slices or [], clock
        self.sample_meta = sample_meta or {}                                        # per sample: Species, Instrument, ReadLength, RunIDs, NRuns, InputFingerprint (what the ledger row needs)
        self.pipeline_version = pipeline_version or cfg.pipeline_version
        self.stop = runner.stop
        self._log_lock = threading.Lock()
        self.log_fn = log or (lambda line: print(line, flush=True))
        self.disk_slots = threading.Semaphore(cfg.disk_slots)                         # samples whose uBAMs are (being) downloaded or still needed
        self.upload_backlog = threading.Semaphore(cfg.max_pending_uploads)           # finished samples waiting locally for an upload slot (bounded: no pile-up of BAMs)
        self.upload_pool = None
        self.ledger = ledger or DiskLedger('/', 0, 0, usage_fn=lambda p: SimpleNamespace(free=float('inf')))        # no ledger given: no disk gate (tests)
        self.ledger.size_fn = self.sample_bytes
        self.download_kwargs = download_kwargs or {}
        self.paths = {}
        self.ready, self.download_failures, self.results = {}, {}, {}
        self.stage_of, self.peak_ratio = {}, {}
        self.bytes_down = self.bytes_up = 0
        self.t_down = self.t_up = 0.0                                                   # seconds spent transferring (summed over transfers)
        self.timings = collections.defaultdict(list)
        self._finish_futures = []
        self._lock = threading.Lock()
        self._batch_over = threading.Event()

    # ---------------------------------------------------------------- logging and bookkeeping
    def log(self, sample, msg):
        with self._log_lock:
            self.log_fn(f'  [{time.strftime("%H:%M:%S")}] [{sample}] {msg}')

    def _set(self, sp, state):
        with self._lock:
            self.stage_of[sp.sample] = state

    def sample_bytes(self, sample):
        """Bytes this sample occupies on disk right now: its uBAMs, its BAM folder, its temp and scratch folders."""
        sp = self.paths.get(sample)
        if sp is None:
            return 0
        total = dir_size(sp.bam_dir) + dir_size(sp.tmp_dir) + sum(os.path.getsize(u) for u in sp.ubams if os.path.exists(u))
        if sp.scratch_dir != sp.tmp_dir:
            total += dir_size(sp.scratch_dir)
        return total

    # ---------------------------------------------------------------- one stage with its checkpoint
    def stage(self, sp, name, fn):
        if stage_done(sp, name):
            self.log(sp.sample, f'{name}: already done, skipped')
            return
        t0 = self.clock()
        self.log(sp.sample, f'{name}: start')
        outputs = fn()
        mark_stage(sp, name, outputs or [])
        minutes = (self.clock() - t0) / 60
        with self._lock:
            self.timings[name].append(minutes)
        self.log(sp.sample, f'{name}: done in {minutes:.1f} min')

    def run(self, sp, cmd, what, label, **kw):
        cores, ram = self.cfg.cost(what)
        nbytes = dir_size(sp.bam_dir) if what == 'upload' else 0
        self.runner.run(cmd, cores, ram, label, sample=sp.sample, priority=self.cfg.priority(what), timeout_s=self.cfg.timeout_for(what, nbytes), **kw)

    # ---------------------------------------------------------------- stages
    @staticmethod
    def _read(path):
        with open(path) as f:
            return f.read()

    def _quickcheck(self, sp, bam, label):
        self.run(sp, ['samtools', 'quickcheck', '-v', bam], 'small', label)

    def _align(self, sp):
        for d in (sp.tmp_dir, sp.bam_dir, sp.scratch_dir):
            os.makedirs(d, exist_ok=True)
        if os.path.exists(sp.sorted_bam):                                   # a crash between the rename and the checkpoint
            try:
                self._quickcheck(sp, sp.sorted_bam, 'quickcheck sorted bam')
                return [sp.sorted_bam]
            except StageError:
                os.remove(sp.sorted_bam)
        parts = []
        for i, ubam in enumerate(sp.ubams):
            t_bam = f'{sp.tmp_dir}{sp.sample}.{i}.sorted.bam'
            if not os.path.exists(t_bam):
                part = t_bam + '.partial.bam'
                logs = [f'{sp.tmp_dir}Alignment_errors_{k}_{i}.txt' for k in ('StF', 'MM', 'Merge')]
                cmd = align_command(ubam, self.index, self.reference, sp.scratch_dir, part, self.cfg.align_threads, self.cfg.align_stf_heap_gb, self.cfg.align_mba_heap_gb, logs)
                self.run(sp, cmd, 'align', f'align uBAM {i + 1}/{len(sp.ubams)}', produces=[part], extra_logs=logs)
                self._quickcheck(sp, part, 'quickcheck aligned bam')
                os.replace(part, t_bam)
            parts.append(t_bam)
        if len(parts) == 1:
            os.replace(parts[0], sp.sorted_bam)
        else:
            part = sp.sorted_bam + '.partial.bam'
            self.run(sp, merge_command(parts, part, sp.scratch_dir, self.cfg.merge_heap_gb), 'merge', 'merge uBAM alignments', stderr=f'{sp.tmp_dir}MergeSamFiles_errors.txt', produces=[part])
            self._quickcheck(sp, part, 'quickcheck merged bam')
            os.replace(part, sp.sorted_bam)
            for p in parts:
                os.remove(p)
        return [sp.sorted_bam]

    def _markdup(self, sp):
        cmd = markdup_command(sp.sorted_bam, sp.all_bam, sp.metrics, sp.scratch_dir, self.cfg.markdup_heap_gb)
        self.run(sp, cmd, 'markdup', 'MarkDuplicates', stderr=f'{sp.tmp_dir}MarkDuplicates_errors.txt', produces=[sp.all_bam, sp.all_bai, sp.metrics])
        self._quickcheck(sp, sp.all_bam, 'quickcheck final bam')
        if os.path.exists(sp.sorted_bam):
            os.remove(sp.sorted_bam)
        return [sp.all_bam, sp.all_bai, sp.metrics]

    def _split(self, sp):
        cmd = [self.cfg.python, SPLIT_SCRIPT, sp.all_bam, sp.split_prefix, '--threads', '2']
        summary = sp.split_prefix + '.split_summary.json'
        self.run(sp, cmd, 'split', 'split BAM (one pass)', stderr=f'{sp.tmp_dir}split_errors.txt', produces=list(sp.split_bams.values()) + [summary])
        outs = list(sp.split_bams.values()) + [summary]
        for bam in sp.split_bams.values():
            self.run(sp, sl.index_command(bam, sp.bai(bam)), 'index', 'index split BAM', produces=[sp.bai(bam)])
            outs.append(sp.bai(bam))
        return outs

    def _stats_all(self, sp):
        self.run(sp, ['samtools', 'stats', '-@', '4', sp.all_bam], 'stats', 'samtools stats', stdout=sp.stats_txt, produces=[sp.stats_txt])
        self.run(sp, ['samtools', 'idxstats', sp.all_bam], 'small', 'samtools idxstats', stdout=sp.idxstats_txt, produces=[sp.idxstats_txt])
        return [sp.stats_txt, sp.idxstats_txt]

    def _slices(self, sp):
        os.makedirs(sp.slice_dir, exist_ok=True)
        outs = []
        for s in self.slices:
            out = sl.slice_path(sp.slice_dir, sp.sample, s)
            part = out + '.partial.bam'
            self.run(sp, sl.view_command(sp.all_bam, part, s), 'slice', f'slice {s.name}', produces=[part])
            os.replace(part, out)
            self.run(sp, sl.index_command(out, sp.bai(out)), 'index', f'index slice {s.name}', produces=[sp.bai(out)])
            outs += [out, sp.bai(out)]
        return outs

    def _panel(self, sp):
        """Module 04's measurements of the final BAM (decision 40). A failure is a PanelFailed, never a StageError: the sample goes on without a JSON."""
        cmd = pn.measure_command(self.cfg.panel, sp.all_bam, sp.sample, sp.panel_json)
        try:
            self.run(sp, cmd, 'panel', 'panel measure', stderr=f'{sp.tmp_dir}panel_errors.txt', produces=[sp.panel_json], extra_logs=[f'{sp.tmp_dir}panel_errors.txt'])
        except StageError as e:
            raise PanelFailed(sp.sample, e.message)
        status, why = pn.identity_status(sp.panel_json, sp.all_bam)
        if status in ('different', 'bad'):                                       # a JSON that names another BAM must never be uploaded as this sample's measurement
            try:
                os.remove(sp.panel_json)
            except OSError:
                pass
            raise PanelFailed(sp.sample, why)
        return [sp.panel_json]

    def _count_ubams(self, sp):
        total = 0
        for i, ubam in enumerate(sp.ubams):
            if not os.path.exists(ubam):
                self.log(sp.sample, 'count_ubam: uBAM no longer on disk, read-count check skipped')
                with open(sp.ubam_reads_txt, 'w') as f:
                    f.write('NA\n')
                return [sp.ubam_reads_txt]
            out = f'{sp.tmp_dir}ubam_count_{i}.txt'
            self.run(sp, ['samtools', 'view', '-c', '-@', '4', ubam], 'count', 'count uBAM reads', stdout=out, produces=[out])
            total += int(self._read(out).strip())
        with open(sp.ubam_reads_txt, 'w') as f:
            f.write(f'{total}\n')
        return [sp.ubam_reads_txt]

    def _post(self, sp):
        """Everything that only needs the final BAM runs side by side; the shared budget decides how many actually run at once."""
        tasks = {'split': self._split, 'stats': self._stats_all}
        if self.slices:
            tasks['slices'] = self._slices
        if self.cfg.panel:
            tasks['panel'] = self._panel
            sp.panel_status, sp.panel_error = 'pending', ''
        with ThreadPoolExecutor(max_workers=len(tasks)) as ex:
            futures = {name: ex.submit(self.stage, sp, name, lambda fn=fn: fn(sp)) for name, fn in tasks.items()}
            errors = []
            for name, fut in futures.items():
                try:
                    fut.result()
                except StopRequested:
                    raise
                except PanelFailed as e:                                          # the panel never fails the sample: record it, go on (the stage is not marked done, so a rerun retries it)
                    sp.panel_status, sp.panel_error = 'failed', e.message
                    self.log(sp.sample, 'PANEL FAILED, the sample goes on without a panel JSON: ' + short_reason(e.message))
                except Exception as e:
                    errors.append(e)
        if errors:
            raise errors[0]
        if self.cfg.panel and sp.panel_status != 'failed':
            sp.panel_status = 'ok'
            status, why = pn.identity_status(sp.panel_json, sp.all_bam)
            sp.panel_note = why if status == 'unknown' else ''                    # accepted, but the row says so

    def _qc(self, sp):
        counts, idx_all = {}, None
        for t in bs.LEGACY_BAM_TYPES:
            bam = sp.all_bam if t == 'all' else sp.split_bams[t]
            out = sp.idxstats_txt if t == 'all' else bam + '.idxstats.txt'
            if t != 'all':
                self.run(sp, ['samtools', 'idxstats', bam], 'small', f'idxstats {t}', stdout=out, produces=[out])
            totals = bs.idxstats_totals(bs.parse_idxstats(self._read(out)))
            counts[t] = totals['records']
            if t == 'all':
                idx_all = totals
        sn = bs.parse_samtools_stats(self._read(sp.stats_txt))
        picard = bs.parse_picard_metrics(self._read(sp.metrics))
        raw = self._read(sp.ubam_reads_txt).strip() if os.path.exists(sp.ubam_reads_txt) else 'NA'
        ubam_reads = int(raw) if raw not in ('NA', '') else None
        cols = {**bs.legacy_columns(counts, self.read_lengths.get(sp.sample), self.ref_size), **bs.derive_columns(sn, idx_all, picard, self.ref_size, ubam_reads, pipeline_version=self.pipeline_version)}
        with open(sp.qc_json, 'w') as f:
            json.dump(cols, f, indent=1)
        if cols['QCWarnings']:
            self.log(sp.sample, 'QC WARNING: ' + cols['QCWarnings'])
        return [sp.qc_json]

    def _upload(self, sp):
        copy_cmd, check_cmd = self.upload_cmds_fn(sp)
        self.run(sp, copy_cmd, 'upload', 'upload to Dropbox', stderr=f'{sp.tmp_dir}upload_errors.txt')
        self.run(sp, check_cmd, 'small', 'verify upload (rclone check)', stderr=f'{sp.tmp_dir}upload_check_errors.txt')
        return []

    def _db(self, sp):
        """The ledger row: the statistics of qc.json plus what is known about the sample and the run (alignment_db.ledger_row decides which of it goes to the sheet). The panel's status is added
        here, not to qc.json: a sample whose panel failed once and succeeded on a resume must not keep an old 'failed' in a QC file whose stage is already done."""
        with open(sp.qc_json) as f:
            qc = json.load(f)
        meta = dict(self.sample_meta.get(sp.sample, {}))
        meta.setdefault('RunIDs', self.run_ids.get(sp.sample, ''))
        warnings = qc.get('QCWarnings', '')
        if sp.panel_status == 'failed':
            warnings = '; '.join(w for w in (warnings, 'panel measurement failed: ' + short_reason(sp.panel_error)) if w)
        if sp.panel_note:
            warnings = '; '.join(w for w in (warnings, sp.panel_note) if w)
        meta.update(SampleID=sp.sample, GenomeVersion=self.genome_version, BamSize=os.path.getsize(sp.all_bam), minimap2_version=self.versions.get('minimap2', ''),
                    gatk_version=self.versions.get('gatk', ''), pysam_version=self.versions.get('pysam', ''), AlignedOn=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                    Machine=socket.gethostname(), PipelineVersion=self.pipeline_version, PanelStatus=sp.panel_status, QCWarnings=warnings)
        self.db_add_fn(adb.ledger_row(qc, meta))
        return []

    def _delete_ubams(self, sp):
        for u in sp.ubams:
            if os.path.exists(u):
                os.remove(u)
        self._release_slot(sp)

    def _release_slot(self, sp):
        with self._lock:
            if sp.holds_slot:
                sp.holds_slot = False
                self.disk_slots.release()

    def _cleanup(self, sp):
        for d in (sp.bam_dir, sp.tmp_dir) + ((sp.scratch_dir,) if sp.scratch_dir != sp.tmp_dir else ()):
            if os.path.basename(d.rstrip('/')) != sp.sample:            # never delete anything that is not this sample's own folder
                raise RuntimeError(f'refusing to delete {d}')
            shutil.rmtree(d, ignore_errors=True)

    # ---------------------------------------------------------------- compute side of one sample: download -> align -> mark duplicates -> split / stats / slices -> qc
    def process(self, sp):
        current, count_pool, count_future = 'download', None, None
        try:
            if not stage_done(sp, 'align'):
                self._set(sp, 'waiting for download')
                self.ready[sp.sample].wait()
                if self.stop.is_set():
                    raise StopRequested()
                if sp.sample in self.download_failures:
                    raise StageError(sp.sample, 'download', '; '.join(self.download_failures[sp.sample]))
            if self.cfg.verify_counts and not stage_done(sp, 'count_ubam'):          # counting the uBAM reads only READS the uBAMs: it runs alongside the alignment
                count_pool = ThreadPoolExecutor(1)
                count_future = count_pool.submit(self.stage, sp, 'count_ubam', lambda: self._count_ubams(sp))
            for current, fn in (('align', self._align), ('markdup', self._markdup)):
                self._set(sp, current)
                self.stage(sp, current, lambda fn=fn: fn(sp))
            if count_future is not None:
                current = 'count_ubam'
                count_future.result()
            self._delete_ubams(sp)                                                    # the uBAMs are no longer needed: free the disk (and the disk slot) now
            current = 'post'
            self._set(sp, 'split/stats/slices')
            self._post(sp)
            current = 'qc'
            self._set(sp, 'qc')
            self.stage(sp, 'qc', lambda: self._qc(sp))
            self.ledger.settle(sp.sample)                                            # the large writes are over
            if not self.cfg.upload:                                                   # a test run: keep everything local, touch neither Dropbox nor the database
                self.results[sp.sample] = ('done', 'local only (--no-upload)')
                self._set(sp, 'done (local)')
                self.ledger.release(sp.sample)
                self.log(sp.sample, 'COMPLETE (local only: nothing uploaded, no database row)')
                return
            current = 'upload'
            while not self.upload_backlog.acquire(timeout=1.0):                      # bounded: finished samples never pile up on disk waiting for an upload slot
                if self.stop.is_set():
                    raise StopRequested()
            self._set(sp, 'queued for upload')
            fut = self.upload_pool.submit(self._finish, sp)                          # the compute worker is free for the next sample from here on
            with self._lock:
                self._finish_futures.append(fut)
        except StopRequested:
            self.results[sp.sample] = ('interrupted', current)
            self.ledger.release(sp.sample)
        except Exception as e:
            msg = e.message if isinstance(e, StageError) else f'{type(e).__name__}: {e}'
            self.results[sp.sample] = ('failed', f'{current}: {msg}')
            self._set(sp, 'failed')
            self.log(sp.sample, f'FAILED at {current}: {msg}')
            self._write_failed(sp, current, msg)
            self.ledger.release(sp.sample)
        finally:
            if count_future is not None:                                              # never leave a command running behind a failed sample
                try:
                    count_future.result()
                except BaseException:
                    pass
                count_pool.shutdown(wait=True)
            self._release_slot(sp)

    @staticmethod
    def _write_failed(sp, current, msg):
        try:
            os.makedirs(sp.tmp_dir, exist_ok=True)
            with open(sp.tmp_dir + 'FAILED.txt', 'w') as f:
                f.write(f'{current}\n{msg}\n')
        except OSError:
            pass

    # ---------------------------------------------------------------- transfer side of one sample (own pool: a slow upload never holds a compute slot)
    def _finish(self, sp):
        current = 'upload'
        try:
            self._set(sp, 'uploading')
            size, t0 = dir_size(sp.bam_dir), time.time()
            self.stage(sp, 'upload', lambda: self._upload(sp))
            secs = max(time.time() - t0, 1e-6)
            with self._lock:
                self.bytes_up += size
                self.t_up += secs
            self.log(sp.sample, f'uploaded {size / 1e9:.1f} GB in {secs / 60:.1f} min ({size / 1e6 / secs:.0f} MB/s)')
            current = 'database'
            self._set(sp, 'database')
            self.stage(sp, 'db', lambda: self._db(sp))
            if not self.cfg.keep_local:
                self._cleanup(sp)
            self.results[sp.sample] = ('done', '')
            self._set(sp, 'done')
            self.log(sp.sample, 'COMPLETE')
        except StopRequested:
            self.results[sp.sample] = ('interrupted', current)
        except Exception as e:
            msg = e.message if isinstance(e, StageError) else f'{type(e).__name__}: {e}'
            self.results[sp.sample] = ('failed', f'{current}: {msg}')
            self._set(sp, 'failed')
            self.log(sp.sample, f'FAILED at {current}: {msg}')
            self._write_failed(sp, current, msg)
        finally:
            self.ledger.release(sp.sample)
            self.upload_backlog.release()

    # ---------------------------------------------------------------- downloads: continuous, and never ahead of the disk
    def _admit(self, sp):
        """Called by the download manager before a sample's files are queued: a disk slot AND room on the disk (ledger)."""
        if not self.disk_slots.acquire(blocking=False):
            return False, f'{self.cfg.disk_slots} samples already have uBAMs on disk'
        ok, why = self.ledger.admit(sp.sample, sp.need_bytes)
        if not ok:
            self.disk_slots.release()
            return False, why
        sp.holds_slot = True
        sp.t_admit = time.time()
        self._set(sp, 'downloading')
        return True, ''

    def _on_ready(self, sp):
        size = sum(os.path.getsize(u) for u in sp.ubams if os.path.exists(u))
        secs = max(time.time() - getattr(sp, 't_admit', time.time()), 1e-6)
        with self._lock:
            self.bytes_down += size
            self.t_down += secs
        self.log(sp.sample, f'downloaded {size / 1e9:.1f} GB in {secs / 60:.1f} min ({size / 1e6 / secs:.0f} MB/s for this sample)')
        self._set(sp, 'downloaded')
        self.ready[sp.sample].set()

    def _on_failed(self, sp, reasons):
        self.download_failures[sp.sample] = reasons
        self._release_slot(sp)
        self.ledger.release(sp.sample)
        self.ready[sp.sample].set()

    def _download_thread(self, samples):
        todo = []
        for sp in samples:
            if stage_done(sp, 'align'):
                self.ready[sp.sample].set()                                  # nothing to download: the aligned BAM already exists
            else:
                todo.append(sp)
                os.makedirs(sp.tmp_dir, exist_ok=True)
        kw = dict(max_parallel=self.cfg.max_downloads, check_workers=self.cfg.check_workers, attempts=self.cfg.download_attempts, backoff=self.cfg.download_backoff)
        kw.update(self.download_kwargs)
        try:
            DownloadManager(todo, self.download_start_fn, self._admit, self._on_ready, self._on_failed, self.stop, check_fn=self.check_fn, idle_fn=self.ledger.idle,
                            log=lambda m: self.log('download', m), **kw).run()
        except Exception as e:                                               # a crashed manager must never leave a worker waiting for ever
            for sp in todo:
                if not self.ready[sp.sample].is_set():
                    self.download_failures.setdefault(sp.sample, [f'download manager crashed: {type(e).__name__}: {e}'])
                    self._release_slot(sp)
                    self.ledger.release(sp.sample)
            self.log('download', f'DOWNLOAD MANAGER CRASHED: {type(e).__name__}: {e}')
        finally:
            for sp in samples:
                self.ready[sp.sample].set()

    # ---------------------------------------------------------------- periodic status: what every sample is doing, the transfers, the disk, the cores
    def status_line(self):
        with self._lock:
            stage_items = list(self.stage_of.items())
            states = collections.Counter(v for _, v in stage_items)
            down, up, t_down, t_up = self.bytes_down, self.bytes_up, self.t_down, self.t_up
        snap = self.ledger.snapshot()
        for s, st in stage_items:                                              # calibrate the disk reservation: the largest size a sample reaches, in uBAM sizes
            sp = self.paths.get(s)
            if sp is not None and st not in ('done', 'done (local)', 'failed', 'waiting for download', 'queued'):
                self.peak_ratio[s] = max(self.peak_ratio.get(s, 0.0), self.sample_bytes(s) / max(sp.need_bytes, 1.0))
        parts = ', '.join(f'{n} {k}' for k, n in sorted(states.items()))
        return (f'STATUS: {parts or "starting"} | downloaded {down / 1e9:.0f} GB, uploaded {up / 1e9:.0f} GB | free disk {snap["free_gb"]:.0f} GB ({snap["reserved_samples"]} samples reserved) | '
                f'cores {self.runner.budget.in_use()}/{self.runner.budget.cores}')

    def _status_loop(self):
        while not self._batch_over.wait(self.cfg.status_s) and not self.stop.is_set():
            self.log('pipeline', self.status_line())
            if self.ledger.snapshot()['free_gb'] < self.cfg.min_free_gb / 2:
                self.log('pipeline', f'LOW DISK: {self.ledger.snapshot()["free_gb"]:.0f} GB free, below half of the {self.cfg.min_free_gb} GB margin')

    # ---------------------------------------------------------------- the whole batch
    def run_batch(self, samples):
        self.paths = {sp.sample: sp for sp in samples}
        for sp in samples:
            self.ready[sp.sample] = threading.Event()
            self.stage_of[sp.sample] = 'queued'
        self.upload_pool = ThreadPoolExecutor(max_workers=self.cfg.uploads)
        feeder = threading.Thread(target=self._download_thread, args=(samples,), daemon=True)
        feeder.start()
        if self.cfg.status_s:
            threading.Thread(target=self._status_loop, daemon=True).start()
        t0 = self.clock()
        ex = ThreadPoolExecutor(max_workers=self.cfg.samples_in_flight)
        futures = [ex.submit(self.process, sp) for sp in samples]
        try:
            for f in futures:
                f.result()
            for f in list(self._finish_futures):                              # compute is over: wait for the uploads that are still running
                f.result()
        except KeyboardInterrupt:
            self.log('pipeline', 'interrupted: stopping all jobs (finished stages are kept; rerun the same command to continue)')
            self.stop.set()
            raise
        finally:
            ex.shutdown(wait=True)
            self.upload_pool.shutdown(wait=True)
            self._batch_over.set()
            feeder.join(timeout=10)
        return self.summary(samples, (self.clock() - t0) / 3600)

    def summary(self, samples, hours):
        done = [s for s, (st, _) in self.results.items() if st == 'done']
        failed = {s: why for s, (st, why) in self.results.items() if st == 'failed'}
        lines = [f'{len(done)} of {len(samples)} samples complete in {hours:.1f} h; {len(failed)} failed']
        for s, why in sorted(failed.items()):
            lines.append(f'  FAILED {s}: {why.splitlines()[0][:200]}')
        for name, v in self.timings.items():
            lines.append(f'  {name:11s} mean {sum(v) / len(v):6.1f} min over {len(v)} samples')
        if self.t_down:
            lines.append(f'  downloaded {self.bytes_down / 1e9:.0f} GB, upload {self.bytes_up / 1e9:.0f} GB; mean per-sample transfer rate: down {self.bytes_down / 1e6 / self.t_down:.0f} MB/s, up {self.bytes_up / 1e6 / max(self.t_up, 1e-6):.0f} MB/s')
        if self.peak_ratio:
            r = list(self.peak_ratio.values())
            lines.append(f'  observed peak local disk use per sample (sampled at status intervals): max {max(r):.1f}x, mean {sum(r) / len(r):.1f}x the uBAM size (the ledger reserves {self.ledger.factor:.1f}x)')
        lines.append(f'  peak use: {self.runner.budget.peak_cores} of {self.runner.budget.cores} cores, {self.runner.budget.peak_ram:.0f} of {self.runner.budget.ram:.0f} GB RAM')
        return '\n'.join(lines), done, failed

    # ---------------------------------------------------------------- what would be run (for review; nothing is executed)
    def plan(self, sp):
        logs = [f'{sp.tmp_dir}Alignment_errors_{k}_0.txt' for k in ('StF', 'MM', 'Merge')]
        c = self.cfg
        steps = [('align', align_command(sp.ubams[0], self.index, self.reference, sp.scratch_dir, sp.tmp_dir + sp.sample + '.0.sorted.bam.partial.bam', c.align_threads,
                                         c.align_stf_heap_gb, c.align_mba_heap_gb, logs)),
                 ('markdup', markdup_command(sp.sorted_bam, sp.all_bam, sp.metrics, sp.scratch_dir, c.markdup_heap_gb)),
                 ('split', [c.python, SPLIT_SCRIPT, sp.all_bam, sp.split_prefix, '--threads', '2']),
                 ('stats', ['samtools', 'stats', '-@', '4', sp.all_bam])]
        steps += [('slice ' + s.name, sl.view_command(sp.all_bam, sl.slice_path(sp.slice_dir, sp.sample, s), s)) for s in self.slices]
        if c.panel:
            steps.append(('panel', pn.measure_command(c.panel, sp.all_bam, sp.sample, sp.panel_json)))
        copy_cmd, check_cmd = self.upload_cmds_fn(sp)
        return steps + [('upload', copy_cmd), ('verify upload', check_cmd)]
