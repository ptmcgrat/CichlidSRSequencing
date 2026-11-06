"""Synthetic tests of the bchrom analysis (called from selftest). No network and no external tools: the downloader, curl, GATK, minimap2 and mosdepth are replaced by fakes that
fabricate the files the real tools would write, so the ORCHESTRATION (stages, checkpoints, resume, cleanup, disk admission, failure isolation) is tested for real, and the
ALGORITHMS (block calling, core set, leave-one-out) are tested on planted signal."""
import gzip, hashlib, json, os, shutil, tempfile, types, zlib
import numpy as np
import pandas as pd


def _ok(name, cond, extra=''):
    print(('PASS ' if cond else 'FAIL ') + name, extra)
    if not cond:
        raise SystemExit(1)


def _raises(fn, exc=Exception):
    try:
        fn()
    except exc:
        return True
    return False


# ------------------------------------------------------------------------------------------------ synthetic data
def filereport(rid, nbytes=2000, bases=2_000_000_000, alias='x', layout='PAIRED', files=None):
    """An ENA filereport as the portal writes it (fastq_ftp without scheme, three ';'-separated lists)."""
    d = f'ftp.sra.ebi.ac.uk/vol1/fastq/{rid[:6]}/00{rid[-1]}/{rid}/'
    files = files or [f'{rid}_1.fastq.gz', f'{rid}_2.fastq.gz']
    ftp = ';'.join(d + f for f in files); md5 = ';'.join(hashlib.md5(f.encode()).hexdigest() for f in files); sz = ';'.join(str(nbytes) for _ in files)
    return ('run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\tbase_count\tread_count\tlibrary_layout\tsample_alias\n'
            f'{rid}\t{ftp}\t{md5}\t{sz}\t{bases}\t{bases // 100}\t{layout}\t{alias}\n')


class FakeCmds:
    java_threads = ''
    align = staticmethod(lambda ubam, index, reference, tmp, out, threads, stf, mba, logs: ['fake-align', ubam, '-O', out])
    markdup = staticmethod(lambda sorted_bam, out, metrics, tmp, heap: ['fake-markdup', sorted_bam, '-O', out, '-M', metrics])


class FakeExecutor:
    """Fabricates what each tool would write. depth_of(run) -> {sequence: per-bin depth array}. fail_once: label that fails the first time it is called for a given run."""
    def __init__(self, depth_of, bin_size, fail_once=None, lengths=None):
        self.depth_of, self.bin, self.fail_once, self.calls, self.failed = depth_of, bin_size, fail_once, [], set()

    def run(self, cmd, label, cwd=None):
        from sd_genotypes.bchrom.workflow import StageError
        run = next((t for t in ' '.join(map(str, cmd)).replace('/', ' ').replace('.', ' ').split() if t.startswith('SRR')), None)
        self.calls.append((label, run))
        if self.fail_once and label == self.fail_once[0] and run == self.fail_once[1] and (label, run) not in self.failed:
            self.failed.add((label, run)); raise StageError(f'{label} failed (fake)')
        if label == 'sra-toolkit':
            if cmd[0] == 'fasterq-dump':
                out = cmd[cmd.index('-O') + 1]
                for k in (1, 2):
                    open(f'{out}/{run}_{k}.fastq', 'w').write('@r\nA\n+\nI\n')
        elif label == 'FastqToSam':
            self.fq = (cmd[cmd.index('--FASTQ') + 1], cmd[cmd.index('--FASTQ2') + 1])
            open(cmd[cmd.index('--OUTPUT') + 1], 'wb').write(b'ubam')
        elif label.startswith('align'):
            open(cmd[cmd.index('-O') + 1], 'wb').write(b'bam')
        elif label == 'MarkDuplicates':
            out = cmd[cmd.index('-O') + 1]
            open(out, 'wb').write(b'bam'); open(out[:-4] + '.bai', 'wb').write(b'bai')
            open(cmd[cmd.index('-M') + 1], 'w').write('## METRICS CLASS\tx\nLIBRARY\tPERCENT_DUPLICATION\nlib\t0.11\n')
        elif label == 'slice':
            self.slice_cmd = cmd
            open(cmd[cmd.index('-o') + 1], 'wb').write(b'slice')
        elif label == 'index slice':
            open(cmd[-1] + '.bai', 'wb').write(b'bai')
        elif label == 'mosdepth':
            prefix = cmd[-2]
            with gzip.open(prefix + '.regions.bed.gz', 'wt') as f:
                for c, arr in self.depth_of(run).items():
                    for k, v in enumerate(arr):
                        f.write(f'{c}\t{k * self.bin}\t{(k + 1) * self.bin}\t{v:.2f}\n')


def make_proj(tmp, upload=False):
    """A stand-in for yhsd.common.Project: only what the bchrom steps use."""
    from cichlid_common.runner import Run
    root = tmp + 'P/G/'
    os.makedirs(root + 'tmp/', exist_ok=True)
    proj = types.SimpleNamespace(fm=types.SimpleNamespace(localMasterDir=tmp, cloudMasterDir=tmp + 'cloud/'), root=root, genome='G', cfg={'runner': {'upload': upload, 'archive_inputs': False}},
                                 dirs={'tmp': root + 'tmp/'}, genome_fa='/ref/G.fa', lg={'chrA': 'LG1', 'chrB': 'LG2'}, chrom_for_lg=lambda lg: {'LG1': 'chrA', 'LG2': 'chrB', 'LG5': 'chrB'}[lg])
    proj.run = lambda analysis, step, params=None, tar=(): Run(proj, analysis, step, params, tar)
    return proj


BIN = 50
LENS = {'chrA': 24000, 'chrB': 16000}                       # bins: 1.2 Mb and 0.8 Mb
CORE = [('chrA', 100_000, 102_000), ('chrA', 400_000, 401_500), ('chrA', 900_000, 903_000), ('chrB', 200_000, 201_200), ('chrB', 600_000, 602_500)]       # in every carrier
VARIABLE = ('chrA', 1_100_000, 1_102_000)                       # in 6 of 12 carriers only
CORE_LOW = ('chrB', 300_000, 301_500)                            # the male pools have only 0.3x here (poorly mapped A sequence): carriers 5x; B-free controls 1x -> SCR 3.3 in the controls too
SHARED_AMP = ('chrB', 700_000, 702_000)                          # an A-genome duplication: in EVERY sample (carriers, pools, controls) -> must not be called


def synth_depth(rng, role, cid, base_depth, core_copies=8):
    d = {c: rng.poisson(base_depth, n).astype(np.float32) for c, n in LENS.items()}
    def amp(reg, factor):
        c, s, e = reg; d[c][s // BIN: e // BIN] = rng.poisson(base_depth * factor, (e - s) // BIN)
    amp(SHARED_AMP, 3.0)
    if role == 'baseline':
        amp(CORE_LOW, 0.3)
    if role == 'carrier':
        for reg in CORE + [CORE_LOW]:
            amp(reg, (2 + core_copies) / 2)
        if cid < 6:
            amp(VARIABLE, 6.0)
    return d


def run():
    from sd_genotypes.bchrom import ena, depth as D, calling as C, workflow as W, steps as S, compare as K
    rng = np.random.default_rng(11)

    # =========================================================================================== ENA
    txt = filereport('SRR9672203', nbytes=7_700_000_000, bases=17293545422, alias='M_lombardoi_2014_1021')
    rec = ena.parse_filereport(txt)['SRR9672203']
    r1, r2 = ena.paired_files(rec)
    _ok('ENA filereport parsed: paired files, md5 and sizes line up', r1['name'] == 'SRR9672203_1.fastq.gz' and r2['name'].endswith('_2.fastq.gz') and rec['base_count'] == 17293545422
        and rec['alias'] == 'M_lombardoi_2014_1021' and r1['bytes'] == 7_700_000_000)
    _ok('ftp paths are served over https', ena.full_url(r1['url']).startswith('https://ftp.sra.ebi.ac.uk/') and ena.full_url('ftp://x/y') == 'ftp://x/y')
    rec3 = ena.parse_filereport(filereport('SRR1', files=['SRR1_1.fastq.gz', 'SRR1_2.fastq.gz', 'SRR1.fastq.gz']))['SRR1']
    _ok('a third (unpaired) FASTQ is ignored for pairing but reported', len(ena.unpaired_files(rec3)) == 1 and ena.paired_files(rec3)[0]['name'] == 'SRR1_1.fastq.gz')
    _ok('a single-end run is refused', _raises(lambda: ena.paired_files(ena.parse_filereport(filereport('SRR2', files=['SRR2.fastq.gz']))['SRR2']), ValueError))
    bad = txt.replace(';', '', 1)
    _ok('lists of different lengths are refused', _raises(lambda: ena.parse_filereport(bad), ValueError))
    _ok('filereport URL asks for the run list with the fields the code needs', 'accession=SRR9672203' in ena.filereport_url('SRR9672203') and 'fastq_md5' in ena.filereport_url('SRR9672203'))

    # ------------- download: resume, corrupted file, failure
    tmp = tempfile.mkdtemp() + '/'
    payload = os.urandom(3000)
    frec = dict(name='A_1.fastq.gz', url='ftp.example/A_1.fastq.gz', md5=hashlib.md5(payload).hexdigest(), bytes=len(payload))

    def curl_factory(script):
        state = {'n': 0}
        def run_(cmd, capture_output=True, text=True):
            dest = cmd[cmd.index('-o') + 1]; state['n'] += 1
            action = script[min(state['n'] - 1, len(script) - 1)]
            have = os.path.getsize(dest) if os.path.exists(dest) else 0
            if action == 'half':
                open(dest, 'ab').write(payload[have:len(payload) // 2]); return types.SimpleNamespace(returncode=18, stderr='partial')
            if action == 'rest':
                open(dest, 'ab').write(payload[have:]); return types.SimpleNamespace(returncode=0, stderr='')
            if action == 'corrupt':
                open(dest, 'wb').write(b'x' * len(payload)); return types.SimpleNamespace(returncode=0, stderr='')
            if action == 'fail':
                return types.SimpleNamespace(returncode=22, stderr='404')
        run_.state = state
        return run_
    quiet = lambda *a, **k: None
    c1 = curl_factory(['half', 'rest'])
    p = ena.download_file(frec, tmp + 'd1/', run=c1, sleep=quiet, log=quiet)
    _ok('download: a partial file is resumed and the md5 is verified', open(p, 'rb').read() == payload and c1.state['n'] == 2 and not os.path.exists(p + '.partial'))
    c2 = curl_factory(['corrupt', 'rest'])
    p = ena.download_file(frec, tmp + 'd2/', run=c2, sleep=quiet, log=quiet)
    _ok('download: a complete but corrupt file is deleted and fetched again', open(p, 'rb').read() == payload and c2.state['n'] == 2)
    c3 = curl_factory(['rest'])
    ena.download_file(frec, tmp + 'd1/', run=c3, sleep=quiet, log=quiet)
    _ok('download: a verified file is not downloaded again', c3.state['n'] == 0)
    _ok('download: gives up after the attempts, with the reason', _raises(lambda: ena.download_file(frec, tmp + 'd3/', attempts=3, run=curl_factory(['fail']), sleep=quiet, log=quiet), ena.DownloadError))
    waits = []
    try:
        ena.download_file(frec, tmp + 'd4/', attempts=3, backoff=10, run=curl_factory(['fail']), sleep=waits.append, log=quiet)
    except ena.DownloadError:
        pass
    _ok('download: back-off grows between attempts', waits == [10, 20], str(waits))
    cmd = W.fastq_to_sam_command('/a_1.fq.gz', '/a_2.fq.gz', '/o.bam', 'SRR1', '/t', 8, '-XX:ParallelGCThreads=2')
    _ok('FastqToSam command: both mates, sample and read group, uBAM out', cmd[:3] == ['gatk', '--java-options', '-Xmx8g -XX:ParallelGCThreads=2'] and cmd[cmd.index('--FASTQ2') + 1] == '/a_2.fq.gz'
        and cmd[cmd.index('--SAMPLE_NAME') + 1] == 'SRR1' and cmd[cmd.index('--PLATFORM') + 1] == 'ILLUMINA')
    try:
        lc = W.LabCommands()
    except SystemExit:
        print('SKIP lab pipeline signatures: helper_modules.seq_pipeline is not importable here (checked on a machine that has lab_repo_patch_v6)')
    else:
        import inspect
        _ok('lab pipeline: align_command / markdup_command still take the parameters bchrom passes, in that order',
            list(inspect.signature(lc.align).parameters) == ['ubam', 'index', 'reference', 'tmp_dir', 'out_bam', 'threads', 'stf_heap_gb', 'mba_heap_gb', 'logs']
            and list(inspect.signature(lc.markdup).parameters) == ['sorted_bam', 'out_bam', 'metrics', 'tmp_dir', 'heap_gb'])

    # =========================================================================================== disk admission
    free = {'v': 500e9}; used = {'a': 0.0, 'b': 0.0}
    led = W.SpaceLedger('/w', min_free_gb=100, free_fn=lambda p: free['v'], usage_fn=lambda p: used[os.path.basename(p)])
    _ok('ledger: first run admitted', led.try_admit('a', 250e9))
    _ok('ledger: a second run that would not fit next to it waits', not led.try_admit('b', 250e9))
    used['a'] = 200e9; free['v'] = 300e9                                    # run a has written 200 GB: free fell by the same, the budget is unchanged
    _ok('ledger: what a run has already written is not counted twice', not led.try_admit('b', 250e9))
    led.release('a'); free['v'] = 500e9; used['a'] = 0.0                  # a run is released after its scratch is deleted: the space is free again
    _ok('ledger: after release the waiting run is admitted', led.try_admit('b', 250e9))

    # =========================================================================================== depth tracks
    arrs = {c: rng.poisson(20, n).astype(np.float32) for c, n in LENS.items()}
    f = tempfile.mkdtemp() + '/x.regions.bed.gz'
    with gzip.open(f, 'wt') as fh:
        for c, a in arrs.items():
            for k, v in enumerate(a):
                fh.write(f'{c}\t{k * BIN}\t{(k + 1) * BIN}\t{v:.2f}\n')
    a2, l2 = D.parse_regions(f, BIN)
    _ok('mosdepth regions parsed into per-sequence arrays', set(a2) == set(LENS) and np.allclose(a2['chrA'], arrs['chrA']) and l2['chrA'] == LENS['chrA'] * BIN)
    with gzip.open(f, 'wt') as fh:
        fh.write('chrA\t0\t50\t1.0\nchrA\t100\t150\t1.0\n')
    _ok('non-contiguous bins are refused (mosdepth run with another --by)', _raises(lambda: D.parse_regions(f, BIN), ValueError))
    st = D.depth_stats(a2, BIN, 500_000)
    _ok('depth statistics: median and mean on the placed chromosomes', abs(st['median'] - 20) < 1 and abs(st['mean'] - 20) < 0.5 and st['scale_chroms'] == 2)
    _ok('scale factor: median or mean, nothing else', D.scale_factor(st, 'mean') == st['mean'] and _raises(lambda: D.scale_factor(st, 'mode'), ValueError))
    base = tempfile.mkdtemp() + '/t'; D.save_track(base, a2, l2, st, BIN, dict(run='x'))
    a3, m3 = D.load_track(base)
    _ok('track round trip', np.array_equal(a3['chrB'], a2['chrB']) and m3['bin_size'] == BIN and m3['run'] == 'x' and m3['lengths']['chrB'] == l2['chrB'])

    # =========================================================================================== calling
    from scipy.stats import binom
    x = np.array([0, 3, 20, 40]); n = np.array([10, 50, 60, 80])
    _ok('binomial tail equals scipy', np.allclose(C.binom_upper_p(x, n, 0.3)[1:], [binom.sf(2, 50, .3), binom.sf(19, 60, .3), binom.sf(39, 80, .3)]) and C.binom_upper_p([0], [10], .3)[0] == 1.0
        and C.binom_upper_p([5], [0], .3)[0] == 1.0)
    _ok('runs of True and the merge distance (300 bp = 6 bins)', C.runs_of_true(np.array([0, 1, 1, 0, 1], bool)) == [(1, 3), (4, 5)]
        and C.merge_runs([(0, 2), (8, 10)], 50, 300) == [[0, 10]] and C.merge_runs([(0, 2), (9, 10)], 50, 300) == [[0, 2], [9, 10]])
    thr = dict(C.DEFAULTS)
    n = 4000
    dm = rng.poisson(25, n).astype(np.float32); db = rng.poisson(15, n).astype(np.float32)
    def put(s_bp, e_bp, f_, arr=db, base=15):
        arr[s_bp // BIN: e_bp // BIN] = rng.poisson(base * f_, (e_bp - s_bp) // BIN)
    put(10_000, 12_000, 5)                                                # 2.0 kb, real
    put(50_000, 50_500, 5)                                                # exactly 500 bp: dropped (<= 500)
    put(60_000, 60_550, 5)                                                # 550 bp: kept
    put(100_000, 100_250, 5); put(100_500, 100_750, 5)                    # two pieces 250 bp apart: merged into one 750 bp block
    put(150_000, 150_300, 5); put(150_650, 150_950, 5)                    # two pieces 350 bp apart: NOT merged, each < 500 bp: gone
    put(160_000, 162_000, 3, dm, 25); put(160_000, 162_000, 3, db, 15)    # A duplication in both: SCR 1
    put(170_000, 172_000, 0.0, dm, 25); put(170_000, 172_000, 5, db, 15)  # baseline has no reads there: no SCR
    sfb, sfm = float(np.median(db)), float(np.median(dm)); pexp = db.mean() / (db.mean() + dm.mean())
    bl = C.call_blocks({'c': db}, {'c': dm}, sfb, sfm, pexp, BIN, thr, lengths={'c': n * BIN})
    got = [(r.start, r.end) for r in bl.itertuples()]
    _ok('planted blocks recovered at bin resolution', (10_000, 12_000) in got and (60_000, 60_550) in got and (100_000, 100_750) in got, str(got))
    _ok('a block of exactly 500 bp is dropped, one of 550 bp kept', (50_000, 50_500) not in got and (60_000, 60_550) in got)
    _ok('pieces 350 bp apart are not merged', not any(150_000 <= s < 151_000 for s, e in got))
    _ok('an A-genome duplication present in both samples is not called', not any(160_000 <= s < 162_000 for s, e in got))
    _ok('no ratio where the baseline has no reads (no calls from an empty baseline)', not any(170_000 <= s < 172_000 for s, e in got))
    _ok('block statistics: SCR ~ (2+c)/2 and the copy-number estimate', abs(bl[bl.start == 10_000].SCR.iloc[0] - 5) < 1.0 and bl[bl.start == 10_000].est_copies.iloc[0] > 6)
    mk = {'c': np.zeros(n, bool)}; mk['c'][10_000 // BIN: 12_000 // BIN] = True
    bl2 = C.call_blocks({'c': db}, {'c': dm}, sfb, sfm, pexp, BIN, thr, mask=mk)
    _ok('regions of known misassembly are not called', not any(10_000 <= r.start < 12_000 for r in bl2.itertuples()) and len(bl2) == len(bl) - 1)
    dbs = db.copy(); lo, hi = 180_000 // BIN, 182_000 // BIN
    dbs[lo:hi] = rng.poisson(15, hi - lo); dbs[lo:hi:4] = 200                       # every fourth bin of a 2 kb stretch passes (gaps of 150 bp merge): 25 % of the bins pass
    kept = C.call_blocks({'c': dbs}, {'c': dm}, sfb, sfm, pexp, BIN, thr, lengths={'c': n * BIN})
    k = kept[(kept.start >= 179_000) & (kept.start < 182_000)]
    dropped = C.call_blocks({'c': dbs}, {'c': dm}, sfb, sfm, pexp, BIN, {**thr, 'min_fraction_passing': 0.5}, lengths={'c': n * BIN})
    _ok('fraction filter: a block in which 25 % of the bins pass is kept at the default (> 10 %) and dropped at 50 %', len(k) == 1 and abs(k.frac_passing.iloc[0] - 0.25) < 0.05
        and not ((dropped.start >= 179_000) & (dropped.start < 182_000)).any(), str(k.frac_passing.tolist()))
    _ok('copy-number equations (Clark Eq. 3 and 4)', C.est_copies(5.0, 1.0, 5.0) == 8.0 and C.est_copies(5.0, 0.8, 4.0) == 8.0 and C.est_copies(5.0, 0.99, 5.0) == 10.0)
    a_, b_ = C.b_space_mb(pd.DataFrame(dict(length=[1000, 2000], est_copies=[10.0, 4.0])))
    _ok('A and B space totals', abs(a_ - 0.003) < 1e-9 and abs(b_ - 0.018) < 1e-9)

    # ------------- core set
    def mkb(rows):
        return pd.DataFrame(rows, columns=['chrom', 'start', 'end'])
    carriers = {f'c{i}': mkb([('chrA', 1000, 3000)] + ([('chrA', 5000, 6000)] if i < 10 else []) + ([('chrA', 8000, 9000)] if i < 11 else []) + ([('chrB', 0, 500)] if i < 3 else [])) for i in range(12)}
    core = C.core_regions(carriers, BIN, 11, 100)
    _ok('core: positions covered by >= 11 of 12 carriers only', [(r.chrom, r.start, r.end) for r in core.itertuples()] == [('chrA', 1000, 3000), ('chrA', 8000, 9000)] and core.n_carriers.tolist() == [12, 11])
    _ok('core: 10 of 12 is not enough; the carriers that overlap each piece are listed', len(C.core_regions(carriers, BIN, 10, 100)) == 3 and len(core.carriers.iloc[1].split(';')) == 11)
    _ok('core: pieces shorter than min_len are dropped', len(C.core_regions({f'c{i}': mkb([('chrA', 1000, 1050)]) for i in range(12)}, BIN, 11, 100)) == 0)
    ragged = {f'c{i}': mkb([('chrA', 1000 + 50 * (i % 3), 3000 - 50 * (i % 2))]) for i in range(12)}
    cr = C.core_regions(ragged, BIN, 11, 100)
    _ok('core: partly overlapping blocks give their intersection (>= 11 carriers)', len(cr) == 1 and cr.start.iloc[0] == 1050 or cr.start.iloc[0] in (1000, 1050, 1100))
    _ok('overlap_bp', C.overlap_bp(mkb([('c', 0, 100), ('c', 500, 600)]), mkb([('c', 50, 550)])) == 100)

    # =========================================================================================== workflow with fake tools
    cfg = S.cfg_of(); cfg['download']['source'] = 'ena'                     # this self-test exercises the ENA route with a fake downloader; the sources are tested in test_bchrom_sources.py
    cfg['bin_size'] = BIN; cfg['scale_min_chrom_bp'] = 500_000; cfg['download']['backoff_s'] = 0
    cfg['qc']['depth_ratio_range'] = [0.5, 1.15]
    cfg['resources']['min_free_gb'] = 0                                   # the real margin (150 GB) is for the real disk; the synthetic runs are a few kB
    by = S.validate_config(cfg)
    ids = list(by)
    carrier_ids = [i for i in ids if by[i]['role'] == 'carrier']
    depth_cache = {}
    def depth_of(rid):
        if rid not in depth_cache:
            role = by[rid]['role']
            depth_cache[rid] = synth_depth(np.random.default_rng(zlib.crc32(rid.encode())), role, carrier_ids.index(rid) if role == 'carrier' else 0, 15 if role == 'carrier' else (25 if role == 'baseline' else 25))
        return depth_cache[rid]
    glen = sum(LENS.values()) * BIN                                       # synthetic genome length; the ENA base count is set so that the expected depth equals the depth of the fake track
    reports = {i: filereport(i, nbytes=1000, bases=int((15 if by[i]['role'] == 'carrier' else 25) * glen), alias=by[i].get('label', i)) for i in ids}
    # the configured bases_gb are for the real data: neutralise the identity check for the synthetic runs
    for s in cfg['runs']:
        s.pop('bases_gb', None)

    fetched = []
    def fetch(i):
        fetched.append(i); return reports[i]
    downloads = []
    def downloader(rec, d):
        os.makedirs(d, exist_ok=True); downloads.append(rec['name'])
        open(d + rec['name'], 'wb').write(b'fq'); return d + rec['name']

    tmp = tempfile.mkdtemp() + '/'
    proj = make_proj(tmp)
    ex = FakeExecutor(depth_of, BIN)
    # one run: stage order, files deleted, track kept
    one = ids[0]
    work = S.work_dir(proj, cfg)
    recs = S.get_records(work, [by[one]], fetch)
    rp = W.RunProcessor(by[one], recs[one], work, cfg, dict(reference='/ref.fa', index='/ref.mmi'), executor=ex, cmds=FakeCmds, downloader=downloader)
    summ = rp.process()
    labels = [c[0] for c in ex.calls]
    _ok('workflow: stages run in order (uBAM, align, quickcheck, markdup, quickcheck, mosdepth)', labels == ['FastqToSam', 'align (SamToFastq | minimap2 | MergeBamAlignment)', 'quickcheck aligned bam', 'MarkDuplicates', 'quickcheck final bam', 'mosdepth'], str(labels))   # this g has no slice regions, so no slice stage
    _ok('workflow: reads, BAMs and scratch are deleted; only the depth track and the summary remain', not os.path.exists(rp.p.dir) and os.path.exists(rp.p.track + '.npz') and os.path.exists(rp.p.summary))
    _ok('workflow: summary carries md5s, the depth check and the duplicate rate', summ['depth_ratio'] is not None and 0.5 < summ['depth_ratio'] < 1.15 and summ['dup_metrics']['PERCENT_DUPLICATION'] == '0.11'
        and len(summ['fastq_md5']) == 2 and not summ['flags'], str(summ['depth_ratio']))
    n_calls = len(ex.calls)
    rp2 = W.RunProcessor(by[one], recs[one], work, cfg, dict(reference='/ref.fa', index='/ref.mmi'), executor=ex, cmds=FakeCmds, downloader=downloader)
    rp2.process()
    _ok('workflow: a finished run is skipped completely on a second invocation', len(ex.calls) == n_calls and len(downloads) == 2)
    # resume after a crash in the align stage: the download and the uBAM are not redone
    two = ids[1]
    ex_f = FakeExecutor(depth_of, BIN, fail_once=('align (SamToFastq | minimap2 | MergeBamAlignment)', two))
    recs2 = S.get_records(work, [by[two]], fetch)
    n_dl = len(downloads)
    rpf = W.RunProcessor(by[two], recs2[two], work, cfg, dict(reference='/ref.fa', index='/ref.mmi'), executor=ex_f, cmds=FakeCmds, downloader=downloader)
    _ok('workflow: a failing stage raises', _raises(rpf.process, W.StageError))
    _ok('workflow: after the crash the finished stages are marked and the scratch is kept', rpf.done('fetch') and rpf.done('ubam') and not rpf.done('align') and os.path.exists(rpf.p.ubam))
    rpf2 = W.RunProcessor(by[two], recs2[two], work, cfg, dict(reference='/ref.fa', index='/ref.mmi'), executor=ex_f, cmds=FakeCmds, downloader=downloader)
    rpf2.process()
    done_labels = [c[0] for c in ex_f.calls if c[1] == two]
    _ok('workflow: the re-run resumes at the failed stage (no new download, no second FastqToSam)', len(downloads) == n_dl + 2 and done_labels.count('FastqToSam') == 1 and rpf2.finished(), str(done_labels))
    four = ids[3]
    def down_fail(rec, d):
        raise ena.DownloadError('ENA is down')
    _ok('an ENA failure is an error when ENA is the only source', _raises(lambda: W.RunProcessor(by[ids[4]], S.get_records(work, [by[ids[4]]], fetch)[ids[4]], work, cfg, dict(reference='/ref.fa', index='/ref.mmi'), executor=FakeExecutor(depth_of, BIN), cmds=FakeCmds, downloader=down_fail).process(), ena.DownloadError))
    # wrong base count -> flag
    three = ids[2]
    bad_rep = filereport(three, nbytes=1000, bases=int(40 * glen), alias='x')
    spec_b = dict(by[three], bases_gb=1.0)
    rb = W.RunProcessor(spec_b, ena.parse_filereport(bad_rep)[three], work, cfg, dict(reference='/ref.fa', index='/ref.mmi'), executor=ex, cmds=FakeCmds, downloader=downloader).process()
    _ok('workflow: a run whose bases do not match the configured accession, or whose depth is off, is flagged', len(rb['flags']) >= 1, str(rb['flags']))
    shutil.rmtree(work + 'depth/'); os.makedirs(work + 'depth/')

    # =========================================================================================== prepare end to end (21 runs, fakes), including failure isolation
    ex2 = FakeExecutor(depth_of, BIN)
    fails = {'n': 0}
    def downloader2(rec, d):
        if rec['name'].startswith(carrier_ids[3] + '_1'):
            fails['n'] += 1; raise ena.DownloadError('simulated outage')
        return downloader(rec, d)
    proj = make_proj(tempfile.mkdtemp() + '/')
    raised = False
    try:
        S.prepare_step(proj, executor=ex2, cmds=FakeCmds, fetch=fetch, downloader=downloader2, cfg=cfg, runs_in_flight=3)
    except SystemExit as e:
        raised = e.code == 1
    work = S.work_dir(proj, cfg)
    done = [i for i in ids if os.path.exists(f'{work}depth/{i}.npz')]
    _ok('prepare: one failed download does not stop the other 20 runs; the exit code reports it', raised and len(done) == 20 and carrier_ids[3] not in done, f'{len(done)} done')
    runs_dirs = os.listdir(proj.root + 'bchrom/prepare/')
    _ok('prepare: provenance run written even though one run failed (runs.tsv with md5s, summaries, ENA listings)', any(os.path.exists(proj.root + 'bchrom/prepare/' + d + '/runs.tsv') for d in runs_dirs))
    S.status_step(proj, cfg=cfg)
    downloader3 = downloader
    S.prepare_step(proj, executor=ex2, cmds=FakeCmds, fetch=fetch, downloader=downloader3, cfg=cfg, runs_in_flight=3)
    _ok('prepare: the re-run completes only the missing run', all(os.path.exists(f'{work}depth/{i}.npz') for i in ids) and sum(1 for c in ex2.calls if c[0] == 'FastqToSam') == 21)
    sl = work + 'slices/'
    _ok('slices: every run keeps a slice and its index although its BAM was deleted', all(os.path.exists(f'{sl}{i}.slice.bam') and os.path.exists(f'{sl}{i}.slice.bam.bai') for i in ids))
    _ok('slices: samtools is asked for both LG5 breakpoint windows (1-based, inclusive), in the duplicate-marked BAM',
        'chrB:19400001-19910000' in ex2.slice_cmd and 'chrB:26470001-26980000' in ex2.slice_cmd and ex2.slice_cmd[5].endswith('.all.bam'), str(ex2.slice_cmd[-3:]))
    pr = proj.root + 'bchrom/prepare/'; last = pr + open(pr + 'LATEST').read().strip() + '/'
    _ok('slices: copied into the prepare run (provenance) with the region definitions, so `clean` does not lose them', sum(f.endswith('.slice.bam') for f in os.listdir(last + 'slices/')) == 21
        and json.load(open(last + 'slices/regions.json'))['LG5_breakpoint_left'] == 'chrB:19400001-19910000')
    rp_s = W.RunProcessor(by[ids[0]], recs[ids[0]] if ids[0] in recs else S.get_records(work, [by[ids[0]]], fetch)[ids[0]], work, cfg, dict(reference='/r', index='/i', slice_regions=['chrB:1-10']), executor=ex2, cmds=FakeCmds)
    os.rename(sl + ids[0] + '.slice.bam', sl + ids[0] + '.hidden')
    _ok('slices: a run whose slice is missing is not "finished" (slices were requested after it ran)', not rp_s.finished())
    os.rename(sl + ids[0] + '.hidden', sl + ids[0] + '.slice.bam')
    cfg_ns = json.loads(json.dumps(cfg)); cfg_ns['slices']['enabled'] = False
    _ok('slices: disabled in the config = no region strings, no slice stage', S.slice_regions(proj, cfg_ns) == [] and len(S.slice_regions(proj, cfg)) == 2)
    _ok('prepare: nothing large is left in the work folder (21 tracks + listings only)', W.dir_size(work) < 30e6 and not any(os.path.isdir(work + i) for i in ids), f'{W.dir_size(work) / 1e6:.1f} MB')
    _ok('preflight: names every missing tool and file before anything is downloaded', S.preflight(dict(reference='/no/ref.fa', index='/no/ref.mmi'), which=lambda t: None if t in ('gatk', 'mosdepth') else '/bin/' + t)
        == ['`gatk` is not on PATH', '`mosdepth` is not on PATH', 'reference file not found: /no/ref.fa', 'index file not found: /no/ref.mmi'])
    _ok('preflight: nothing to report when all is there', S.preflight(dict(reference='/r', index='/i'), which=lambda t: '/bin/' + t, exists=lambda p: True) == [])
    _ok('prepare: unknown --only ids are refused', _raises(lambda: S.prepare_step(proj, only=['SRR0'], executor=ex2, cmds=FakeCmds, fetch=fetch, downloader=downloader, cfg=cfg), SystemExit))

    # =========================================================================================== call on the 21 synthetic tracks
    core_df, summ = S.call_step(proj, cfg=cfg)
    found = [(r.chrom, r.start, r.end) for r in core_df.itertuples()]
    def covered(reg):
        return any(c == reg[0] and s <= reg[1] + 100 and e >= reg[2] - 100 for c, s, e in found)
    _ok('call: every planted core block (in all 12 carriers) is in the core set', all(covered(r) for r in CORE + [CORE_LOW]), str(found))
    _ok('call: the block present in only 6 carriers is not core, the A-genome amplification (in every sample) is not core',
        not any(c == VARIABLE[0] and s < VARIABLE[2] and e > VARIABLE[1] for c, s, e in found) and not any(c == SHARED_AMP[0] and s < SHARED_AMP[2] and e > SHARED_AMP[1] for c, s, e in found))
    tiers = dict(zip(zip(core_df.chrom, core_df.start), core_df.tier))
    low_tier = [v for (c, s0), v in tiers.items() if c == CORE_LOW[0] and abs(s0 - CORE_LOW[1]) <= 100]
    _ok('call: a block that the male pools barely cover is tier B, the others are tier A', low_tier == ['B'] and (core_df.tier == 'A').sum() == len(CORE), str(core_df.tier.tolist()))
    ca = summ[summ.role == 'carrier']; co = summ[summ.role == 'control']
    _ok('call: the B-free controls show the low-baseline block (a false positive of a single B-free sample, as in Clark\'s controls) and nothing else', all(len(S_blocks) == 1 for S_blocks in [pd.read_csv(proj.root + 'bchrom/call/' + open(proj.root + 'bchrom/call/LATEST').read().strip() + '/' + f, sep='\t') for f in os.listdir(proj.root + 'bchrom/call/' + open(proj.root + 'bchrom/call/LATEST').read().strip() + '/') if 'Female' in f]))
    _ok('call: tier A signal is reported separately (carriers high, controls ~1)', ca.tierA_median_scaled.min() > 3 and co.tierA_median_scaled.max() < 1.5 and ca.tierA_n_blocks.min() == len(CORE))
    _ok('call: carriers have high signal over the core (in sample and leave-one-out), controls ~1', ca.median_scaled.min() > 3 and ca.loo_median_scaled.min() > 3 and co.median_scaled.max() < 1.5 and ca.loo_frac_ge2.min() > 0.9
        and co.frac_ge2.max() < 0.2, f'carriers {ca.median_scaled.min():.1f}, LOO {ca.loo_median_scaled.min():.1f}, controls {co.median_scaled.max():.2f}')
    _ok('call: leave-one-out defines the core without the carrier (10 of the other 11) and still finds the planted blocks', ca.loo_core_blocks.min() >= len(CORE))
    out = proj.root + 'bchrom/call/'; rid = open(out + 'LATEST').read().strip()
    files = os.listdir(out + rid + '/')
    _ok('call: outputs written (BED, tables, per-carrier blocks, summary)', 'core_blocks_G.bed' in files and 'core_blocks.tsv' in files and 'sample_summary.tsv' in files and 'block_signal.tsv' in files
        and 'summary.txt' in files and sum(f.startswith('blocks_') for f in files) == 14, str(sorted(files)))
    txt = open(out + rid + '/summary.txt').read()
    _ok('call: the summary states the checks and compares with Clark', 'OK' in txt and 'Clark' in txt and 'leave-one-out' in txt and 'CHECK' not in txt, txt.split('\n')[0])
    bed = pd.read_csv(out + rid + '/core_blocks_G.bed', sep='\t', header=None)
    _ok('call: BED has 5 columns, core ids and the carrier count', bed.shape[1] == 5 and bed[3].str.startswith('core').all())
    # exclusion of core pieces that are also called in a control
    cfg_x = json.loads(json.dumps(cfg)); cfg_x['core']['exclude_blocks_in_controls'] = True
    core_x, _ = S.call_step(proj, cfg=cfg_x)
    _ok('call: the optional control filter drops the block that a B-free control also shows (the low-baseline block) and keeps the rest', len(core_x) == len(core_df) - 1 and (core_x.tier == 'A').all())
    for i in carrier_ids[:2]:
        os.remove(f'{work}depth/{i}.npz')
    _ok('call: refuses to run with fewer carriers than core.min_individuals (10 < 11)', _raises(lambda: S.call_step(proj, cfg=cfg), SystemExit))

    # =========================================================================================== compare with Clark's blocks
    from cichlid_common import paf as P
    def pline(q, qlen, qs, qe, st, t, tlen, ts, te, cg):
        n = min(qe - qs, te - ts); return '\t'.join(map(str, [q, qlen, qs, qe, st, t, tlen, ts, te, n, max(qe - qs, te - ts), 60])) + '\tcg:Z:' + cg
    # Clark's lg1 = chrA shifted by +1000; lg2 = chrB reversed; one block on an unplaced contig that is not aligned
    pf = tempfile.mkdtemp() + '/x.paf'
    open(pf, 'w').write(pline('lg1', 1_300_000, 0, 1_200_000, '+', 'chrA', 1_200_000, 0, 1_200_000, '1200000M') + '\n' + pline('lg2', 800_000, 0, 800_000, '-', 'chrB', 800_000, 0, 800_000, '800000M') + '\n')
    cb = tempfile.mkdtemp() + '/clark.bed'
    open(cb, 'w').write('lg1 100100 101900 1800 1 1 1 1 1\nlg2 %d %d 1200 1 1 1 1 1\nunplaced1 0 800 800 1 1 1 1 1\nlg1 50000 50900 900 1 1 1 1 1\n' % (800_000 - 201_200, 800_000 - 200_000))
    S.compare_step(proj, cb, paf=pf, cfg=cfg)
    cmp_dir = proj.root + 'bchrom/compare/' + open(proj.root + 'bchrom/compare/LATEST').read().strip() + '/'
    cj = json.load(open(cmp_dir + 'comparison.json'))
    _ok('compare: Clark blocks on chromosomes lift (also on the minus strand) and are recovered; unplaced and non-core ones are reported separately',
        cj['all']['n_clark_blocks'] == 4 and cj['anchored']['n_lifted'] == 3 and cj['anchored']['n_recovered'] == 2 and cj['unplaced']['n_lifted'] == 0 and cj['unplaced']['n_clark_blocks'] == 1, json.dumps(cj['all']))

    # =========================================================================================== clean
    n_before = W.dir_size(work)
    S.clean_step(proj, cfg=cfg)
    _ok('clean: deletes the work folder once a `call` run exists', not os.path.exists(work) and n_before > 0)
    other = tempfile.mkdtemp() + '/'
    cfg_o = json.loads(json.dumps(cfg)); cfg_o['work_dir'] = other
    os.makedirs(other + 'depth/'); os.remove(other + '.bchrom_work') if os.path.exists(other + '.bchrom_work') else None
    open(other + 'precious.txt', 'w').write('x')
    shutil.rmtree(proj.root + 'bchrom/call')
    S.work_dir(proj, cfg)
    open(proj.root + 'bchrom/call/LATEST', 'w').write('x') if os.path.isdir(proj.root + 'bchrom/call') else None
    os.makedirs(proj.root + 'bchrom/call/', exist_ok=True); open(proj.root + 'bchrom/call/LATEST', 'w').write('x')
    _ok('clean: refuses to delete a folder that bchrom did not create (a work_dir pointing at shared data)', _raises(lambda: S.clean_step(proj, cfg=cfg_o), SystemExit) and os.path.exists(other + 'precious.txt'))
    shutil.rmtree(proj.root + 'bchrom/call')
    _ok('clean: refuses to delete the only copy of the tracks when no `call` run exists', _raises(lambda: S.clean_step(proj, cfg=cfg), SystemExit))

    # =========================================================================================== the configuration itself
    real = S.cfg_of(); rb = S.validate_config(real)
    ena_bases = {'SRR9672203': 17293545422, 'SRR9672201': 14720395288, 'SRR9672202': 15391497262, 'SRR9672205': 14738619930, 'SRR9672206': 16373908910, 'SRR9672207': 15169642076,
                 'SRR9672212': 16360590242, 'SRR9672208': 14681886412, 'SRR9672204': 11866558882, 'SRR9672211': 13405532646, 'SRR9672209': 14932552858, 'SRR9672210': 12608853130,
                 'SRR18083689': 24812117530, 'SRR18083719': 37009940454, 'SRR18083711': 13308147032, 'SRR18083714': 12784761396, 'SRR18083735': 29997863146}
    _ok('config: the base counts of the 12 carriers and 5 pools equal ENA\'s (transcribed from the listings of PRJNA552724 / PRJNA808791)',
        all(abs(rb[i]['bases_gb'] - b / 1e9) < 0.006 for i, b in ena_bases.items()), str([i for i, b in ena_bases.items() if abs(rb[i]['bases_gb'] - b / 1e9) >= 0.006]))
    alias = {'SRR9672203': 'M_lombardoi_2014_1021', 'SRR9672201': 'M_mbenji_2012_3997', 'SRR9672202': 'M_zeb_boadzulu_2005_0976', 'SRR9672205': 'M_greshakei_2012_3493', 'SRR9672206': 'M_lombardoi_2014_1018',
             'SRR9672207': 'L_trewavasae_2005_1306', 'SRR9672212': 'M_zebra_NkhBay_2012_5347', 'SRR9672208': 'M_auratus_2008_1601', 'SRR9672204': 'M_lombardoi_2014_1108',
             'SRR9672211': 'M_zebra_NkhBay_2012_5340', 'SRR9672209': 'M_zeb_boadzulu_2005_0983', 'SRR9672210': 'M_zeb_boadzulu_2005_0986'}
    _ok('config: carrier labels equal the ENA sample aliases', all(rb[i]['label'] == a for i, a in alias.items()))
    pair = {'M_lombardoi': 'SRR18083689', 'M_zeb_boadzulu': 'SRR18083689', 'M_mbenji': 'SRR18083735', 'M_greshakei': 'SRR18083723', 'L_trewavasae': 'SRR18083714', 'M_auratus': 'SRR18083711', 'M_zebra_NkhBay': 'SRR18083731'}
    _ok('config: every carrier is paired with the male pool of its population (Clark Table 1; M. lombardoi uses the Boadzulu pool)', all(rb[i]['baseline'] == pair[next(k for k in pair if a.startswith(k))] for i, a in alias.items()))
    _ok('config: 12 carriers, 7 baselines, 2 controls; the two controls use the Maison male pool', [sum(r['role'] == k for r in real['runs']) for k in ('carrier', 'baseline', 'control')] == [12, 7, 2]
        and all(rb[i]['baseline'] == 'SRR18083719' for i in rb if rb[i]['role'] == 'control'))
    _ok('config: published A and B space of Clark Table 3 are recorded for the carriers', abs(sum(rb[i]['published_A_mb'] for i in alias) - 95.23) < 0.05)
    print('bchrom selftests passed')
