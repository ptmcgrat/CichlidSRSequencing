"""The orchestration, tested with a fake command executor: scheduling, resource budget, checkpoints, failure isolation, download window, database rows and the exact commands.
No GATK, minimap2, samtools, rclone or Google access is needed; what is real here is every line of seq_pipeline.py except the subprocess call itself."""
import json, os, shutil, sys, tempfile, threading, time, unittest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from helper_modules import seq_pipeline as sp
from helper_modules import bam_slices as sl
from test_bam_stats import STATS, IDXSTATS, PICARD


class FakeProc:
    pid = 1
    def poll(self):
        return 0


class Fake:
    """Stands in for subprocess: creates the files a command promises, writes stdout for the samtools queries, can fail on demand, records what ran and when."""
    def __init__(self, nubams, fail=None, delay=0.003, delays=None, hook=None):
        self.nubams, self.fail, self.delay = nubams, fail or (lambda kind, sample, text: False), delay
        self.delays, self.hook, self.ends = delays or {}, hook, []
        self.calls, self.lock, self.t0 = [], threading.Lock(), time.time()

    @staticmethod
    def kind(cmd):
        c = [x for x in cmd if x != 'nice']
        t = ' '.join(c)
        if c[0] == 'bash':
            return 'align'
        for k, key in (('markdup', 'MarkDuplicates'), ('merge', 'MergeSamFiles'), ('split', 'split_bamfile_single_pass'), ('stats', 'samtools stats'), ('idxstats', 'idxstats'),
                       ('index', 'samtools index'), ('quickcheck', 'quickcheck'), ('copy', 'rclone copy'), ('check', 'rclone check')):
            if key in t:
                return k
        if t.startswith('samtools view -c'):
            return 'count'
        if t.startswith('samtools view -b'):
            return 'slice'
        return 'other'

    def sample_of(self, text):
        return next((s for s in self.nubams if s in text), '?')

    def __call__(self, cmd, stdout, stderr, produces, stop):
        text, kind = ' '.join(cmd), self.kind(cmd)
        sample = self.sample_of(text)
        with self.lock:
            self.calls.append((time.time() - self.t0, kind, sample, text))
        if self.hook:
            self.hook(kind, sample)
        end = time.time() + self.delays.get(kind, self.delay)
        while time.time() < end:
            if stop.is_set():
                raise sp.StopRequested()
            time.sleep(0.001)
        if self.fail(kind, sample, text):
            import re
            for path in ([stderr] if stderr else []) + re.findall(r'2> (\S+)', text):          # the alignment pipe writes its three error logs via 2>
                with open(path, 'w') as f:
                    f.write('simulated failure\nsecond line\n')
            return 1
        for p in produces:
            os.makedirs(os.path.dirname(p), exist_ok=True)
            with open(p, 'w') as f:
                f.write('x')
        if kind == 'markdup':
            with open(cmd[cmd.index('-M') + 1], 'w') as f:
                f.write(PICARD)
        if stdout:
            body = {'stats': STATS, 'idxstats': IDXSTATS, 'count': f'{100_000_000 // self.nubams[sample]}\n'}.get(kind, '')
            with open(stdout, 'w') as f:
                f.write(body)
        with self.lock:
            self.ends.append((time.time() - self.t0, kind, sample))
        return 0

    def kinds(self, sample=None):
        return [k for _, k, s, _ in self.calls if sample in (None, s)]


def build(tmp, nubams, fake=None, prefetch=2, in_flight=3, max_cores=24, align_threads=8, keep_local=False, slices=None, delay=0.003, downloads=None, fail=None, nice=0, dl_kwargs=None,
          uploads=2, max_pending=None, ledger=None, need_bytes=None, scratch=False, delays=None, hook=None):
    fake = fake or Fake(nubams, fail=fail, delay=delay, delays=delays, hook=hook)
    stop = threading.Event()
    cfg = sp.Config(max_cores=max_cores, max_ram_gb=400, samples_in_flight=in_flight, align_threads=align_threads, prefetch=prefetch, max_downloads=4, uploads=uploads, keep_local=keep_local, nice=nice, max_pending_uploads=max_pending, status_s=0)
    budget = sp.ResourceBudget(cfg.max_cores, cfg.max_ram_gb, stop)
    runner = sp.Runner(budget, nice=nice, stop=stop, exec_fn=fake)
    paths = [sp.SamplePaths(s, [f'{tmp}/reads/{s}_{i}.bam' for i in range(n)], f'{tmp}/bam/{s}', f'{tmp}/temp/{s}', scratch_dir=(f'{tmp}/scratch/{s}' if scratch else None), need_bytes=need_bytes) for s, n in nubams.items()]
    for p in paths:
        os.makedirs(os.path.dirname(paths[0].ubams[0]), exist_ok=True)
    rows, events = [], []
    on_disk = {'max': 0}

    def start(job):
        with open(job.path, 'w') as f:
            f.write('bam')
        occ = sum(1 for p in paths if any(os.path.exists(u) for u in p.ubams))
        on_disk['max'] = max(on_disk['max'], occ)
        events.append(('download', job.key, time.time()))
        return FakeProc()

    def upload_cmds(p):
        return ['rclone', 'copy', p.bam_dir, 'cloud/' + p.sample], ['rclone', 'check', p.bam_dir, 'cloud/' + p.sample, '--one-way', '--size-only']

    def db(row):
        rows.append(row)
        events.append(('db', row['SampleID'], time.time()))

    pipe = sp.SeqPipeline(cfg, runner, '/ref.fna', '/ref.mmi', 962_000_000, downloads or start, upload_cmds, db, run_ids={s: f'run_{s}' for s in nubams}, read_lengths={s: 300 for s in nubams},
                          versions={'minimap2': '2.28', 'gatk': '4.5', 'pysam': '0.22'}, genome_version='Mzebra_GT3_NCBI', log=lambda l: None, check_fn=lambda p: (os.path.exists(p), ''), slices=slices, ledger=ledger,
                          download_kwargs=dict(kill_fn=lambda proc: None, poll=0.002, admit_poll=0.01, backoff=0.0, **(dl_kwargs or {})))
    return pipe, paths, fake, rows, events, on_disk


class TestPipeline(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_two_samples_end_to_end(self):
        slices = [sl.Slice('LG2_bp_left', 'NC_135168.1', 19550000, 20050000)]
        pipe, paths, fake, rows, events, _ = build(self.tmp, {'SMP_1': 1, 'SMP_2': 2}, slices=slices)
        text, done, failed = pipe.run_batch(paths)
        self.assertEqual((sorted(done), failed), (['SMP_1', 'SMP_2'], {}), text)
        self.assertEqual(sorted(r['SampleID'] for r in rows), ['SMP_1', 'SMP_2'])
        r = next(r for r in rows if r['SampleID'] == 'SMP_1')
        for key in ('Coverage', 'TotalReads', 'UnmappedReads', 'DuplicationReads', 'PrimaryReads', 'UnmappedFraction', 'DuplicateFraction', 'QCWarnings', 'RunIDs', 'BamSize', 'minimap2_version'):
            self.assertIn(key, r)
        self.assertEqual(r['PrimaryReads'], 100_000_000)
        self.assertAlmostEqual(r['DuplicateFraction'], 0.09)
        self.assertAlmostEqual(r['UnmappedFraction'], 0.03)
        self.assertEqual(r['QCWarnings'], '', r['QCWarnings'])
        self.assertEqual(r['RunIDs'], 'run_SMP_1')
        k1, k2 = fake.kinds('SMP_1'), fake.kinds('SMP_2')
        self.assertLess(k1.index('align'), k1.index('markdup'))
        self.assertLess(k1.index('markdup'), k1.index('split')); self.assertLess(k1.index('markdup'), k1.index('stats'))
        self.assertLess(max(i for i, k in enumerate(k1) if k in ('split', 'stats', 'slice', 'count', 'idxstats')), k1.index('copy'))
        self.assertLess(k1.index('copy'), k1.index('check'))
        self.assertIn('slice', k1)
        self.assertEqual(k2.count('align'), 2, 'one alignment per uBAM')
        self.assertEqual(k2.count('merge'), 1, 'two uBAM alignments are merged')
        self.assertEqual(k1.count('merge'), 0)
        ev = {(e, s): t for e, s, t in events}
        calls = {(k, s): t for t, k, s, _ in reversed(fake.calls)}                      # time of the first call of each kind
        self.assertGreater(ev[('db', 'SMP_1')], calls[('check', 'SMP_1')], 'the database row is written only after the upload was verified')
        for p in paths:                                                                   # keep_local=False: everything local is gone after the verified upload
            self.assertFalse(os.path.exists(p.bam_dir) or os.path.exists(p.tmp_dir) or any(os.path.exists(u) for u in p.ubams), p.sample)

    def test_multi_uBAM_sample_qc_has_no_false_warning(self):
        pipe, paths, fake, rows, *_ = build(self.tmp, {'SMP_2': 2})
        pipe.run_batch(paths)
        self.assertEqual(rows[0]['QCWarnings'], '', 'two uBAMs of 50M reads must add up to the 100M primary reads')

    def test_keep_local_keeps_bams_but_not_ubams(self):
        pipe, paths, *_ = build(self.tmp, {'SMP_1': 1}, keep_local=True)
        pipe.run_batch(paths)
        self.assertTrue(os.path.exists(paths[0].all_bam) and os.path.exists(paths[0].qc_json))
        self.assertFalse(os.path.exists(paths[0].ubams[0]))

    def test_failure_of_one_sample_does_not_stop_the_others_and_frees_its_disk_slot(self):
        bad = lambda kind, sample, text: kind == 'align' and sample == 'SMP_B'
        pipe, paths, fake, rows, *_ = build(self.tmp, {'SMP_A': 1, 'SMP_B': 1, 'SMP_C': 1, 'SMP_D': 1}, in_flight=1, prefetch=0, fail=bad)      # ONE disk slot: a leaked slot would deadlock C and D
        text, done, failed = pipe.run_batch(paths)
        self.assertEqual(sorted(done), ['SMP_A', 'SMP_C', 'SMP_D'], text)
        self.assertEqual(list(failed), ['SMP_B'])
        self.assertIn('align', failed['SMP_B']); self.assertIn('exit code 1', failed['SMP_B']); self.assertIn('simulated failure', failed['SMP_B'], 'the tail of the error log is reported')
        self.assertEqual(sorted(r['SampleID'] for r in rows), ['SMP_A', 'SMP_C', 'SMP_D'], 'a failed sample never gets a database row')
        self.assertTrue(os.path.exists(paths[1].tmp_dir + 'FAILED.txt'))
        self.assertEqual(fake.kinds('SMP_B').count('markdup'), 0, 'nothing runs after a failed stage')

    def test_a_late_failure_is_resumed_without_redoing_the_early_stages(self):
        fail_upload = lambda kind, sample, text: kind == 'copy'
        pipe, paths, fake, rows, *_ = build(self.tmp, {'SMP_1': 1}, fail=fail_upload)
        _, done, failed = pipe.run_batch(paths)
        self.assertEqual((done, list(failed)), ([], ['SMP_1']))
        first = fake.kinds('SMP_1')
        self.assertTrue(os.path.exists(paths[0].all_bam), 'the final BAM survives a failed upload')
        pipe2, paths2, fake2, rows2, events2, _ = build(self.tmp, {'SMP_1': 1})              # a new run: same folders
        _, done2, failed2 = pipe2.run_batch(paths2)
        self.assertEqual((done2, failed2), (['SMP_1'], {}))
        redone = set(fake2.kinds('SMP_1'))
        self.assertTrue(redone.isdisjoint({'align', 'markdup', 'split', 'stats', 'slice', 'count'}), f'finished stages were repeated: {redone}')
        self.assertEqual([e for e in events2 if e[0] == 'download'], [], 'nothing is downloaded again once the alignment exists')
        self.assertEqual(len(rows2), 1)

    def test_core_and_ram_budget_is_never_exceeded_and_samples_overlap(self):
        pipe, paths, fake, rows, *_ = build(self.tmp, {f'S{i}': 1 for i in range(8)}, in_flight=6, max_cores=24, align_threads=8, delay=0.01)
        _, done, failed = pipe.run_batch(paths)
        self.assertEqual((len(done), failed), (8, {}))
        b = pipe.runner.budget
        self.assertLessEqual(b.peak_cores, 24)
        self.assertGreaterEqual(b.peak_cores, 16, 'two alignments (10 cores each) should have run together')

    def test_downloads_run_only_a_window_ahead_and_alignment_starts_before_the_last_download(self):
        n = 10
        pipe, paths, fake, rows, events, on_disk = build(self.tmp, {f'S{i:02d}': 1 for i in range(n)}, in_flight=2, prefetch=1, delay=0.01)
        _, done, failed = pipe.run_batch(paths)
        self.assertEqual(len(done), n)
        self.assertLessEqual(on_disk['max'], pipe.cfg.disk_slots + 1, f'{on_disk["max"]} samples had uBAMs on disk at once')
        first_align = min(t for t, k, s, _ in fake.calls if k == 'align')
        last_download = max(t for e, s, t in events if e == 'download') - fake.t0
        self.assertLess(first_align, last_download, 'alignment must start while later samples are still downloading')

    def test_a_download_that_never_works_fails_only_that_sample(self):
        def start(job):
            if job.key == 'SMP_B':
                raise FileNotFoundError('Cant find file for download: SMP_B')
            with open(job.path, 'w') as f:
                f.write('bam')
            return FakeProc()
        pipe, paths, fake, rows, *_ = build(self.tmp, {'SMP_A': 1, 'SMP_B': 1, 'SMP_C': 1}, downloads=start, dl_kwargs=dict(attempts=2))
        _, done, failed = pipe.run_batch(paths)
        self.assertEqual((sorted(done), list(failed)), (['SMP_A', 'SMP_C'], ['SMP_B']))
        self.assertIn('download', failed['SMP_B'])

    def test_stop_request_ends_the_run_quickly_and_leaves_nothing_hanging(self):
        pipe, paths, fake, rows, *_ = build(self.tmp, {f'S{i}': 1 for i in range(6)}, delay=0.3)
        threading.Timer(0.4, pipe.stop.set).start()
        t0 = time.time()
        _, done, failed = pipe.run_batch(paths)
        self.assertLess(time.time() - t0, 10)
        self.assertLess(len(done), 6)
        self.assertEqual(failed, {}, 'an interrupt is not a failure')


    def test_no_upload_mode_keeps_everything_local_and_touches_neither_dropbox_nor_the_database(self):
        pipe, paths, fake, rows, *_ = build(self.tmp, {'SMP_1': 1})
        pipe.cfg.upload = False
        _, done, failed = pipe.run_batch(paths)
        self.assertEqual((done, failed, rows), (['SMP_1'], {}, []))
        self.assertTrue(set(fake.kinds()).isdisjoint({'copy', 'check'}))
        self.assertTrue(os.path.exists(paths[0].all_bam) and os.path.exists(paths[0].qc_json))

    def test_read_count_check_is_skipped_not_failed_when_the_ubams_are_gone(self):
        pipe, paths, fake, rows, *_ = build(self.tmp, {'SMP_1': 1}, fail=lambda kind, sample, text: kind == 'copy')
        pipe.run_batch(paths)                                                           # fails at the upload; uBAMs are already deleted
        os.remove(paths[0].tmp_dir + 'stages/count_ubam.json'); os.remove(paths[0].ubam_reads_txt)
        pipe2, paths2, fake2, rows2, *_ = build(self.tmp, {'SMP_1': 1})
        _, done, failed = pipe2.run_batch(paths2)
        self.assertEqual((done, failed), (['SMP_1'], {}))
        self.assertIn('UBamReads', rows2[0])


class TestOverlapAndSpace(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_a_slow_upload_does_not_hold_a_compute_slot(self):
        """TWO compute workers, four samples, every upload takes 1 s while a sample's compute takes ~0.3 s. If a worker stayed with its sample through the upload (the old design),
        samples 2 and 3 could only start computing after the uploads of samples 0 and 1 had finished. Here they must start long before."""
        pipe, paths, fake, rows, *_ = build(self.tmp, {f'SMP_{i}': 1 for i in range(4)}, in_flight=2, uploads=2, max_pending=4, delays={'copy': 1.0}, delay=0.01)
        _, done, failed = pipe.run_batch(paths)
        self.assertEqual((len(done), failed), (4, {}))
        copy_end = {smp: t for t, k, smp in fake.ends if k == 'copy'}
        first_align = {smp: min(t for t, k, sm, _ in fake.calls if sm == smp and k == 'align') for smp in copy_end}
        self.assertLess(first_align['SMP_2'], copy_end['SMP_0'] - 0.4, 'sample 2 waited for the upload of sample 0: a compute slot was held during the upload')
        self.assertLess(first_align['SMP_3'], copy_end['SMP_1'] - 0.4, 'sample 3 waited for the upload of sample 1: a compute slot was held during the upload')

    def test_the_number_of_simultaneous_uploads_is_the_upload_pool_size(self):
        for pool in (1, 2):
            tmp = tempfile.mkdtemp()
            try:
                pipe, paths, fake, rows, *_ = build(tmp, {f'SMP_{i}': 1 for i in range(6)}, in_flight=3, uploads=pool, max_pending=6, delays={'copy': 0.15}, delay=0.005)
                _, done, failed = pipe.run_batch(paths)
                self.assertEqual((len(done), failed), (6, {}))
                spans = [(t0, t1) for (t0, k, s, _) in fake.calls if k == 'copy' for (t1, k2, s2) in fake.ends if k2 == 'copy' and s2 == s]
                most = max(sum(1 for (a0, a1) in spans if a0 <= t0 + 1e-6 < a1 - 1e-6) for (t0, _) in spans)
                self.assertEqual(most, pool, f'with {pool} upload slot(s) the most uploads running at once was {most}')
            finally:
                shutil.rmtree(tmp, ignore_errors=True)

    def test_finished_samples_never_pile_up_waiting_for_an_upload_slot(self):
        """Uploads are slow (0.4 s) and the compute side is fast; at most `max_pending` finished samples may be waiting or uploading."""
        pipe, paths, fake, rows, *_ = build(self.tmp, {f'SMP_{i}': 1 for i in range(8)}, in_flight=4, uploads=1, max_pending=2, delays={'copy': 0.4}, delay=0.005)
        peak, stop = {'n': 0}, threading.Event()

        def monitor():
            while not stop.is_set():
                with pipe._lock:
                    n = sum(1 for v in pipe.stage_of.values() if v in ('queued for upload', 'uploading', 'database'))
                peak['n'] = max(peak['n'], n); time.sleep(0.003)
        m = threading.Thread(target=monitor); m.start()
        _, done, failed = pipe.run_batch(paths)
        stop.set(); m.join()
        self.assertEqual((len(done), failed), (8, {}))
        self.assertLessEqual(peak['n'], 2, f'{peak["n"]} finished samples were waiting for upload at once')

    def test_ubams_are_deleted_before_the_post_stages_and_the_read_count_runs_alongside_the_alignment(self):
        seen = {}

        def hook(kind, sample):
            sp = next(p for p in paths if p.sample == sample) if sample != '?' else None
            if sp is not None:
                seen.setdefault((kind, sample), any(os.path.exists(u) for u in sp.ubams))
        pipe, paths, fake, rows, *_ = build(self.tmp, {'SMP_1': 1, 'SMP_2': 2}, hook=hook)
        pipe.run_batch(paths)
        for s in ('SMP_1', 'SMP_2'):
            self.assertTrue(seen[('align', s)], 'the uBAM is needed for the alignment')
            self.assertTrue(seen[('count', s)], 'the read count reads the uBAM')
            self.assertFalse(seen[('split', s)] or seen[('stats', s)] or seen[('slice', s)] if ('slice', s) in seen else seen[('split', s)] or seen[('stats', s)], 'the uBAMs must be gone before the split and the statistics')
            first = lambda kind: min(t for t, k, sm, _ in fake.calls if k == kind and sm == s)
            self.assertLess(first('count'), first('markdup'), 'the read count must run alongside the alignment, not after it')

    def test_drain_first_a_later_stage_is_served_before_a_new_alignment_that_asked_earlier(self):
        b = sp.ResourceBudget(10, 100)
        held, order = b.acquire(10, 1), []

        def align():
            h = b.acquire(10, 1, priority=2); order.append('align'); b.release(h)

        def markdup():
            h = b.acquire(2, 1, priority=1); order.append('markdup'); b.release(h)
        t1 = threading.Thread(target=align); t1.start(); time.sleep(0.1)
        t2 = threading.Thread(target=markdup); t2.start(); time.sleep(0.1)
        b.release(held); t1.join(5); t2.join(5)
        self.assertEqual(order, ['markdup', 'align'])

    def test_the_disk_gate_inside_the_pipeline_limits_how_many_samples_are_downloaded_at_once(self):
        from types import SimpleNamespace
        from helper_modules.disk_ledger import DiskLedger
        ledger = DiskLedger('/data', 3.5, 200, usage_fn=lambda p: SimpleNamespace(free=400e9))        # (400 - 200) / (3.5 x 15) = 3 samples at a time
        pipe, paths, fake, rows, events, on_disk = build(self.tmp, {f'SMP_{i}': 1 for i in range(10)}, in_flight=6, prefetch=2, ledger=ledger, need_bytes=15e9, delay=0.01)
        _, done, failed = pipe.run_batch(paths)
        self.assertEqual((len(done), failed), (10, {}), 'the samples must all get through, one batch of three after another')
        self.assertLessEqual(on_disk['max'], 3, f'{on_disk["max"]} samples had uBAMs on disk at once; the ledger allows 3')
        self.assertTrue(ledger.idle(), 'every reservation must be released at the end')

    def test_a_sample_that_can_never_fit_fails_fast_with_a_clear_message_instead_of_hanging(self):
        from types import SimpleNamespace
        from helper_modules.disk_ledger import DiskLedger
        ledger = DiskLedger('/data', 3.5, 200, usage_fn=lambda p: SimpleNamespace(free=100e9))         # 100 GB free: nothing can be admitted
        pipe, paths, fake, rows, *_ = build(self.tmp, {f'SMP_{i}': 1 for i in range(3)}, ledger=ledger, need_bytes=15e9, dl_kwargs=dict(idle_wait=0.05))
        t0 = time.time()
        _, done, failed = pipe.run_batch(paths)
        self.assertLess(time.time() - t0, 15)
        self.assertEqual((done, sorted(failed)), ([], ['SMP_0', 'SMP_1', 'SMP_2']))
        self.assertIn('not enough free disk', failed['SMP_0']); self.assertEqual(fake.kinds(), [], 'no work may start for a sample that was never admitted')

    def test_scratch_files_go_to_the_tmp_dir_and_are_removed_with_the_sample(self):
        pipe, paths, fake, rows, *_ = build(self.tmp, {'SMP_1': 1, 'SMP_2': 2}, scratch=True)
        pipe.run_batch(paths)
        for sample in ('SMP_1', 'SMP_2'):
            texts = [t for _, k, s, t in fake.calls if s == sample and k in ('align', 'markdup', 'merge')]
            self.assertTrue(texts and all(f'{self.tmp}/scratch/{sample}/' in t for t in texts), 'GATK --TMP_DIR must point at the scratch folder')
            self.assertFalse(os.path.exists(f'{self.tmp}/scratch/{sample}'), 'the scratch folder is deleted with the sample')

    def test_the_status_line_and_the_peak_disk_ratio(self):
        pipe, paths, fake, rows, *_ = build(self.tmp, {'SMP_1': 1}, need_bytes=1000.0)
        pipe.paths = {p.sample: p for p in paths}
        pipe.stage_of['SMP_1'] = 'markdup'
        os.makedirs(paths[0].bam_dir); os.makedirs(paths[0].tmp_dir)
        with open(paths[0].bam_dir + 'x.bam', 'wb') as f:
            f.write(b'1' * 2500)                                                          # 2.5 uBAM sizes on disk at this moment
        line = pipe.status_line()
        self.assertIn('1 markdup', line); self.assertIn('free disk', line); self.assertIn('cores 0/24', line)
        self.assertAlmostEqual(pipe.peak_ratio['SMP_1'], 2.5, places=2)
        pipe.results['SMP_1'] = ('done', '')
        self.assertIn('observed peak local disk use per sample', pipe.summary(paths, 0.1)[0])

    def test_markdup_heap_is_configurable(self):
        c = sp.markdup_command('a', 'b', 'm', 't', 32)
        self.assertIn('-Xmx32g', c[2])


class TestCommands(unittest.TestCase):
    def test_alignment_is_one_pipe_with_pipefail_and_the_old_options(self):
        cmd = sp.align_command('/r/u.bam', '/g/ref.mmi', '/g/ref.fna', '/t/', '/t/o.bam', 16, 4, 8, ['/t/a', '/t/b', '/t/c'])
        self.assertEqual(cmd[:3], ['bash', '-o', 'pipefail'])
        pipe = cmd[-1]
        self.assertEqual(pipe.count(' | '), 2)
        for opt in ('SamToFastq', '--CLIPPING_ATTRIBUTE XT', '--CLIPPING_ACTION 2', '--INTERLEAVE true', '--NON_PF true', 'minimap2 -x sr -a -t 16 /g/ref.mmi -', 'MergeBamAlignment',
                    '--UNMAPPED_BAM /r/u.bam', '--ALIGNED_BAM /dev/stdin', '--ADD_MATE_CIGAR true', '--CLIP_ADAPTERS false', '--CLIP_OVERLAPPING_READS true',
                    '--INCLUDE_SECONDARY_ALIGNMENTS true', '--MAX_INSERTIONS_OR_DELETIONS -1', '--PRIMARY_ALIGNMENT_STRATEGY MostDistant', '--ATTRIBUTES_TO_RETAIN XS', '--FASTQ /dev/stdout'):
            self.assertIn(opt, pipe)
        for log in ('2> /t/a', '2> /t/b', '2> /t/c'):
            self.assertIn(log, pipe)

    def test_markdup_matches_the_old_command(self):
        old = ['MarkDuplicates', '-I', 'sorted', '-O', 'all', '-M', 'all.metrics', '--TMP_DIR', 'tmp', '--SORTING_COLLECTION_SIZE_RATIO', '.2', '--CREATE_INDEX']
        new = sp.markdup_command('sorted', 'all', 'all.metrics', 'tmp', 16)
        self.assertEqual(new[new.index('MarkDuplicates'):], old)
        self.assertIn('-Xmx16g', new[2])
        self.assertIn('-XX:ParallelGCThreads=2', new[2])

    def test_no_gvcf_anywhere_in_the_plan(self):
        pipe, paths, *_ = build(tempfile.mkdtemp(), {'SMP_1': 1}, slices=[sl.Slice('x', 'c', 1, 100)])
        text = ' '.join(sp.sh_join(c) for _, c in pipe.plan(paths[0]))
        for word in ('HaplotypeCaller', 'GatherVcfs', 'g.vcf', 'CountReads'):
            self.assertNotIn(word, text)
        self.assertIn('samtools view -b', text)

    def test_naming_must_match_filemanager(self):
        class FM:                       # the attribute formulas of FileManager.createSampleFiles
            localRawBamFiles = ['/r/x.bam']
            def createSampleFiles(self, s):
                self.localSampleBamDir = f'/d/Bamfiles/G/{s}/'; self.localSampleTempDir = f'/d/Temp/{s}/'
                self.localTempSortedBamFile = self.localSampleTempDir + s + '.sorted.bam'
                self.localBamFile = self.localSampleBamDir + s + '.all.bam'
                self.localUnmappedBamFile = self.localSampleBamDir + s + '.unmapped.bam'
                self.localChimericBamFile = self.localSampleBamDir + s + '.chimeric.bam'
        p = sp.SamplePaths.from_fm(FM(), 'SAMEA1')
        self.assertEqual((p.all_bam, p.sorted_bam), ('/d/Bamfiles/G/SAMEA1/SAMEA1.all.bam', '/d/Temp/SAMEA1/SAMEA1.sorted.bam'))
        class Drifted(FM):
            def createSampleFiles(self, s):
                super().createSampleFiles(s); self.localBamFile = self.localSampleBamDir + s + '.bam'
        with self.assertRaises(RuntimeError):
            sp.SamplePaths.from_fm(Drifted(), 'SAMEA1')


class TestBudget(unittest.TestCase):
    def test_first_come_first_served_so_a_big_job_is_not_starved(self):
        b = sp.ResourceBudget(10, 100)
        small = b.acquire(2, 1)
        order = []
        def big():
            h = b.acquire(10, 1); order.append('big'); b.release(h)
        def later():
            time.sleep(0.1); h = b.acquire(1, 1); order.append('small'); b.release(h)
        t1, t2 = threading.Thread(target=big), threading.Thread(target=later)
        t1.start(); time.sleep(0.05); t2.start(); time.sleep(0.2)
        self.assertEqual(order, [], 'neither may start while the big request waits at the head')
        b.release(small)
        t1.join(5); t2.join(5)
        self.assertEqual(order, ['big', 'small'])


class TestDatabase(unittest.TestCase):
    def setUp(self):
        import pandas as pd, types
        self.pd = pd
        self.calls = []
        mod = types.ModuleType('gspread_dataframe')
        mod.set_with_dataframe = lambda ws, df, **kw: self.calls.append((ws, len(df), kw))
        sys.modules['gspread_dataframe'] = mod
        from helper_modules.alignment_db import DatabaseWriter, LEGACY_KEYS
        self.LEGACY_KEYS = LEGACY_KEYS
        class Sheet:
            def worksheet(s, name): return ('worksheet', name)
        class GC:
            def open_by_key(s, key): return Sheet()
        class FM:
            gc, a_ID, genome_version = GC(), 'id', 'Mzebra_GT3_NCBI'
            alignment_dt = pd.DataFrame([dict(SampleID='OLD', TotalReads=5)])
            def readAlignmentDatabase(s): pass
        self.fm = FM()
        self.writer = DatabaseWriter(self.fm, attempts=2, backoff=0, sleep=lambda s: None)

    def tearDown(self):
        sys.modules.pop('gspread_dataframe', None)

    def row(self, sample='NEW', **extra):
        return dict({k: 1 for k in self.LEGACY_KEYS}, SampleID=sample, **extra)

    def test_added_with_new_columns_and_the_worksheet_is_resized(self):
        self.assertEqual(self.writer.add_row(self.row(PrimaryReads=7)), 'added')
        ws, n, kw = self.calls[0]
        self.assertEqual((ws, n, kw), (('worksheet', 'Mzebra_GT3_NCBI'), 2, {'resize': True}))
        self.assertIn('PrimaryReads', self.fm.alignment_dt.columns)

    def test_existing_sample_is_left_alone_and_missing_legacy_columns_are_refused(self):
        self.assertEqual(self.writer.add_row(self.row('OLD')), 'exists')
        self.assertEqual(self.calls, [])
        with self.assertRaises(ValueError):
            self.writer.add_row({'SampleID': 'X'})

    def test_a_failing_write_raises_instead_of_being_ignored(self):
        sys.modules['gspread_dataframe'].set_with_dataframe = lambda *a, **k: (_ for _ in ()).throw(RuntimeError('quota exceeded'))
        with self.assertRaises(RuntimeError) as cm:
            self.writer.add_row(self.row())
        self.assertIn('quota exceeded', str(cm.exception))


if __name__ == '__main__':
    unittest.main()
