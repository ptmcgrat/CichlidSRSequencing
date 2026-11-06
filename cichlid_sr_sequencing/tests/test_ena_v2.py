"""Tests for downloadShortReadData_v2 / grab_ena_v2 / ena_fetch / ena_scheduler. No network, no GATK, no ascp: those are replaced by fakes that fabricate what the real tools would write.
python3 -m unittest discover -s tests -v   (from cichlid_sr_sequencing/)"""
import gzip, hashlib, io, json, os, sys, tempfile, types, unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
from helper_modules import ena_fetch as E, ena_scheduler as SCH, grab_ena_v2 as W
import downloadShortReadData_v2 as D

HERE = os.path.dirname(os.path.abspath(__file__))


def rec(run, files=('_1', '_2'), layout='PAIRED', size=1000, orphan_first=False, aspera=True, **kw):
    names = [f'{run}{s}.fastq.gz' for s in files]
    if orphan_first:
        names = [f'{run}.fastq.gz'] + names
    base = f'ftp.sra.ebi.ac.uk/vol1/fastq/{run[:6]}/00{run[-1]}/{run}/'
    d = dict(run_accession=run, sample_accession='SAMEA1', secondary_sample_accession='', experiment_accession='ERX1', study_accession='PRJEB1', library_name='lib', library_layout=layout,
             library_source='GENOMIC', library_selection='RANDOM', instrument_platform='ILLUMINA', instrument_model='HiSeq', scientific_name='Fish', read_count='100', base_count='20000',
             fastq_ftp=';'.join(base + n for n in names), fastq_aspera=';'.join(base.replace('ftp.sra.ebi.ac.uk/', 'fasp.sra.ebi.ac.uk:/') + n for n in names) if aspera else '',
             fastq_md5=';'.join(hashlib.md5(n.encode()).hexdigest() for n in names), fastq_bytes=';'.join(str(size) for _ in names), submitted_ftp='', submitted_format='')
    d.update(kw)
    return d


class TestSelect(unittest.TestCase):
    def test_orphan_listed_first_is_not_taken_as_mate_1(self):          # the positional-selection bug of the old script
        s = E.select_fastqs(rec('ERR1', orphan_first=True))
        self.assertEqual(s['status'], 'ok'); self.assertEqual(s['r1']['name'], 'ERR1_1.fastq.gz'); self.assertEqual(s['r2']['name'], 'ERR1_2.fastq.gz')
        self.assertEqual([o['name'] for o in s['orphans']], ['ERR1.fastq.gz']); self.assertIn('ignored', s['note'])

    def test_classification(self):
        self.assertEqual(E.select_fastqs(rec('ERR2', files=('',), layout='SINGLE'))['status'], 'single_end')
        n = E.select_fastqs(rec('ERR3', fastq_ftp='', fastq_md5='', fastq_bytes='', submitted_ftp='ftp.sra.ebi.ac.uk/x/y.cram', submitted_format='CRAM'))
        self.assertEqual(n['status'], 'no_fastq'); self.assertIn('CRAM', n['note'])
        self.assertEqual(E.select_fastqs(rec('ERR4', files=('_1',)))['status'], 'missing_mate')
        self.assertEqual(E.select_fastqs(rec('ERR5', fastq_md5='abc'))['status'], 'bad_listing')

    def test_aspera_path_from_ena_or_derived(self):
        a = E.select_fastqs(rec('ERR6'))['r1']['aspera']; b = E.select_fastqs(rec('ERR6', aspera=False))['r1']['aspera']
        self.assertTrue(a.startswith('fasp.sra.ebi.ac.uk:/vol1/fastq/')); self.assertEqual(a, b)

    def test_runinfo_columns(self):
        ri = E.runinfo_from_ena(pd.DataFrame([rec('ERR7', read_count='100', base_count='25000')]))
        self.assertEqual(list(ri.columns), D.RUNINFO_COLUMNS); self.assertEqual(int(ri.AvgSpotLen[0]), 250)


class TestDownload(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp() + '/'; self.payload = os.urandom(3000)
        self.f = dict(name='A_1.fastq.gz', url='ftp.example/A_1.fastq.gz', aspera='fasp.sra.ebi.ac.uk:/vol1/A_1.fastq.gz', md5=hashlib.md5(self.payload).hexdigest(), bytes=3000)
        self.q = lambda *a, **k: None

    def curl(self, script):
        st = {'n': 0}
        def run(cmd, capture_output=True, text=True):
            dest = cmd[cmd.index('-o') + 1]; st['n'] += 1; act = script[min(st['n'] - 1, len(script) - 1)]
            have = os.path.getsize(dest) if os.path.exists(dest) else 0
            if act == 'half':
                open(dest, 'ab').write(self.payload[have:1500]); return types.SimpleNamespace(returncode=18, stderr='partial')
            if act == 'rest':
                open(dest, 'ab').write(self.payload[have:]); return types.SimpleNamespace(returncode=0, stderr='')
            if act == 'corrupt':
                open(dest, 'wb').write(b'x' * 3000); return types.SimpleNamespace(returncode=0, stderr='')
            return types.SimpleNamespace(returncode=22, stderr='404')
        run.st = st
        return run

    def test_https_resume_and_corrupt(self):
        c = self.curl(['half', 'rest']); p = E.https_download(self.f, self.tmp, run=c, sleep=self.q, log=self.q)
        self.assertEqual(open(p, 'rb').read(), self.payload); self.assertEqual(c.st['n'], 2)
        c2 = self.curl(['corrupt', 'rest']); p2 = E.https_download(dict(self.f, name='B.gz'), self.tmp, run=c2, sleep=self.q, log=self.q)
        self.assertEqual(c2.st['n'], 2)
        with self.assertRaises(E.EnaError):
            E.https_download(dict(self.f, name='C.gz'), self.tmp, attempts=2, run=self.curl(['fail']), sleep=self.q, log=self.q)

    def test_ascp_command(self):
        a = E.select_fastqs(rec('ERR9'))
        cmd = E.ascp_command([a['r1'], a['r2']], '/dest/', '/bin/ascp', '/key.pem', rate='300m')
        self.assertEqual(cmd[:9], ['/bin/ascp', '-QT', '-l', '300m', '-P', '33001', '-k', '1', '-i'])
        self.assertEqual(cmd[-3:-1], ['era-fasp@' + a['r1']['aspera'], 'era-fasp@' + a['r2']['aspera']]); self.assertEqual(cmd[-1], '/dest/')

    def test_ascp_error_kinds(self):
        self.assertEqual(E.classify_ascp_error('Error: No such file or directory'), 'permanent')
        self.assertEqual(E.classify_ascp_error('Session Stop (Error: Server aborted session: Too many connections)'), 'throttle')
        self.assertEqual(E.classify_ascp_error('Completed: 0 files'), 'transient')

    def fake_ascp(self, write, rc=0, err=''):
        def run(cmd, capture_output=True, timeout=None):
            if write:
                for p in cmd[9 + 1:-1]:
                    pass
                dest = cmd[-1]
                for src in [c for c in cmd if c.startswith('era-fasp@')]:
                    open(os.path.join(dest, os.path.basename(src)), 'wb').write(self.payload)
            return types.SimpleNamespace(returncode=rc, stderr=err.encode())
        return run

    def sel(self):
        mk = lambda n: dict(name=n, url='ftp.example/' + n, aspera='fasp.sra.ebi.ac.uk:/v/' + n, md5=self.f['md5'], bytes=3000)
        return dict(status='ok', r1=mk('R_1.fastq.gz'), r2=mk('R_2.fastq.gz'), orphans=[])

    def test_download_run_routes(self):
        ascp = os.path.join(self.tmp, 'ascp'); open(ascp, 'w').write('x')
        calls = []
        def dispatcher(asp, https):                                           # ascp calls go to `asp`, curl calls to `https`
            def run(cmd, **k):
                calls.append(cmd[0])
                return asp(cmd, **k) if cmd[0] == ascp else https(cmd)
            return run
        throttled = self.fake_ascp(False, 1, 'Session Stop: Too many connections')
        # 1. aspera works: one session for both mates
        _, _, st = E.download_run(self.sel(), self.tmp + 'd2/', ascp=ascp, key='k', log=self.q, run=dispatcher(self.fake_ascp(True), None))
        self.assertEqual((st['method'], calls.count(ascp), 'curl' in calls), ('aspera', 1, False))
        # 2. aspera throttled -> https fallback delivers
        _, _, st = E.download_run(self.sel(), self.tmp + 'd3/', ascp=ascp, key='k', log=self.q, sleep=self.q, run=dispatcher(throttled, self.curl(['rest'])))
        self.assertEqual(st['method'], 'https')
        # 3. aspera throttled AND https failing -> the error says throttle (so that the scheduler slows down)
        with self.assertRaises(E.EnaError) as cm:
            E.download_run(self.sel(), self.tmp + 'd4/', ascp=ascp, key='k', log=self.q, sleep=self.q, attempts=1, run=dispatcher(throttled, self.curl(['fail'])))
        self.assertEqual(cm.exception.kind, 'throttle')
        # 4. permanent ascp error ('no such file'): no https fallback
        calls.clear()
        with self.assertRaises(E.EnaError) as cm:
            E.download_run(self.sel(), self.tmp + 'd5/', ascp=ascp, key='k', log=self.q, sleep=self.q, run=dispatcher(self.fake_ascp(False, 1, 'No such file or directory'), self.curl(['rest'])))
        self.assertEqual((cm.exception.kind, 'curl' in calls), ('permanent', False))
        # 5. both mates in separate sessions when asked
        calls.clear()
        E.download_run(self.sel(), self.tmp + 'd6/', ascp=ascp, key='k', log=self.q, both_in_one_session=False, run=dispatcher(self.fake_ascp(True), None))
        self.assertEqual(calls.count(ascp), 2)
        # 6. no ascp at all -> https
        _, _, st = E.download_run(self.sel(), self.tmp + 'd7/', ascp=None, log=self.q, run=dispatcher(None, self.curl(['rest'])))
        self.assertEqual(st['method'], 'https')

    def test_wrong_md5_is_deleted_not_resumed(self):
        ascp = os.path.join(self.tmp, 'ascp'); open(ascp, 'w').write('x'); s = self.sel(); s['r1']['md5'] = 'bad'
        ok, kind, msg = E.aspera_download([s['r1']], self.tmp + 'w/', ascp, 'k', run=self.fake_ascp(True))
        self.assertFalse(ok); self.assertFalse(os.path.exists(self.tmp + 'w/R_1.fastq.gz'))


class TestGovernor(unittest.TestCase):
    def test_aimd(self):
        t = [0.0]; g = SCH.Governor(lo=1, hi=4, up_after=2, cooldown=60, now=lambda: t[0])
        self.assertEqual(g.allowed, 4); g.on_throttle(); self.assertEqual(g.allowed, 2); self.assertFalse(g.may_start(0))
        t[0] = 61; self.assertTrue(g.may_start(1)); self.assertFalse(g.may_start(2))
        g.on_success(); g.on_success(); self.assertEqual(g.allowed, 3); g.on_throttle(); g.on_throttle(); self.assertEqual(g.allowed, 1)
        t[0] = 500; g.on_success(); g.on_success(); self.assertEqual(g.allowed, 2)


class TestSanitize(unittest.TestCase):
    def gz(self, path, text):
        with gzip.open(path, 'wt') as f:
            f.write(text)

    def test_drops_pairs_with_an_empty_mate_and_renames(self):
        d = tempfile.mkdtemp() + '/'
        self.gz(d + 'a.gz', '@r1/1 x\nACGT\n+\nIIII\n@r2/1\n\n+\n\n@r3/1\nGG\n+\nII\n'); self.gz(d + 'b.gz', '@r1/2\nTTTT\n+\nIIII\n@r2/2\nAC\n+\nII\n@r3/2\nCC\n+\nII\n')
        n, kept = E.sanitize_pair(d + 'a.gz', d + 'b.gz', d + 'o1', d + 'o2')
        self.assertEqual((n, kept), (3, 2)); self.assertEqual(open(d + 'o1').read().split('\n')[0], '@r1 1:N:0:2'); self.assertEqual(open(d + 'o2').read().split('\n')[0], '@r1 2:N:0:2')

    def test_unequal_mate_files_raise_instead_of_being_zipped(self):
        d = tempfile.mkdtemp() + '/'
        self.gz(d + 'a.gz', '@r1\nA\n+\nI\n@r2\nA\n+\nI\n'); self.gz(d + 'b.gz', '@r1\nA\n+\nI\n')
        with self.assertRaises(E.EnaError) as cm:
            E.sanitize_pair(d + 'a.gz', d + 'b.gz', d + 'o1', d + 'o2')
        self.assertEqual(cm.exception.kind, 'permanent')


class FakeFM:
    def __init__(self, tmp, in_db=()):
        self.localReadDownloadDir, self.localTempDir, self.localReadsDir = tmp + 'dl/', tmp + 'temp/', tmp + 'reads/'
        for d in (self.localReadDownloadDir, self.localTempDir, self.localReadsDir):
            os.makedirs(d, exist_ok=True)
        self.reads_dt = pd.DataFrame(dict(RunID=list(in_db))); self.cloud, self.rows, self.uploaded, self.drop_upload = {}, [], [], set()

    def readSampleDatabase(self): pass
    def downloadData(self, path): pass
    def uploadData(self, path):
        self.uploaded.append(path)
        if os.path.basename(path) not in self.drop_upload:
            self.cloud.setdefault(os.path.dirname(path) + '/', set()).add(os.path.basename(path))
    def returnCloudFiles(self, d): return sorted(self.cloud.get(d, set()))
    def returnFileSize(self, p): return 123
    def addDNAReadRow(self, read_data, sample_data): self.rows.append((read_data, sample_data))


def fake_gatk(fail_for=(), empty_error_for=(), log=None):
    def run(cmd, capture_output=True):
        sub = cmd[cmd.index('FastqToSam') if 'FastqToSam' in cmd else cmd.index('MarkIlluminaAdapters')]
        run_id = next(c for c in cmd if 'ERR' in c)[-30:]
        rid = next((r for r in list(fail_for) + list(empty_error_for) if r in ' '.join(cmd)), None)
        if log is not None:
            log.append((sub, cmd))
        if sub == 'FastqToSam' and rid in empty_error_for and not any('fixed_' in c for c in cmd):
            return types.SimpleNamespace(returncode=3, stderr=b'Exception: Sequence and quality lengths differ for read')
        if rid in fail_for:
            return types.SimpleNamespace(returncode=1, stderr=b'Exception: boom in ' + sub.encode())
        out = cmd[cmd.index('--OUTPUT') + 1] if sub == 'FastqToSam' else cmd[cmd.index('-O') + 1]
        open(out, 'wb').write(b'bam')
        if sub == 'MarkIlluminaAdapters':
            open(cmd[cmd.index('-M') + 1], 'w').write('m')
        return types.SimpleNamespace(returncode=0, stderr=b'')
    return run


def make_spec(tmp, run_id, platform='ILLUMINA', sel=None):
    return dict(run=run_id, sample='SAMEA1', library='lib', platform=platform, sel=sel or E.select_fastqs(rec(run_id)), temp_dir=tmp + 'temp/', output_bam=tmp + f'reads/PRJEB1/{run_id}.unmapped_marked_adapters.bam',
                status_json=tmp + f'temp/{run_id}.status.json')


def fake_download(raise_for=None):
    def dl(sel, dest, **kw):
        if raise_for and sel['r1']['name'].startswith(raise_for[0]) and raise_for[1]:
            e = raise_for[1]; raise_for_copy = e
            raise E.EnaError('simulated', e)
        os.makedirs(dest, exist_ok=True)
        for k in ('r1', 'r2'):
            open(os.path.join(dest, sel[k]['name']), 'wb').write(b'fq')
        return os.path.join(dest, sel['r1']['name']), os.path.join(dest, sel['r2']['name']), dict(method='aspera', seconds=1.0, bytes=2000, mb_per_s=2.0)
    return dl


class TestWorker(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp() + '/'; self.fm = FakeFM(self.tmp); os.makedirs(self.fm.localTempDir, exist_ok=True)

    def left(self):
        return [f for f in os.listdir(self.fm.localTempDir) if not f.endswith('.json')]

    def test_success_fast_path_uploads_verifies_and_cleans(self):
        log = []; sanit = []
        st = W.process_run(make_spec(self.tmp, 'ERR100'), self.fm, run=fake_gatk(log=log), download=fake_download(), sanitize=lambda *a: sanit.append(1), log=lambda *a: None)
        self.assertEqual((st['state'], st['kind']), ('done', '')); self.assertEqual(sanit, [])
        self.assertEqual([c[0] for c in log], ['FastqToSam', 'MarkIlluminaAdapters']); self.assertEqual(len(self.fm.uploaded), 2); self.assertEqual(self.left(), [])
        self.assertEqual(json.load(open(self.tmp + 'temp/ERR100.status.json'))['state'], 'done')

    def test_rejected_fastq_goes_through_the_sanitiser_fallback(self):
        san = []
        def sanitize(a, b, c, d): san.append((c, d)); open(c, 'w').write('x'); open(d, 'w').write('x'); return 10, 9
        st = W.process_run(make_spec(self.tmp, 'ERR101'), self.fm, run=fake_gatk(empty_error_for=['ERR101']), download=fake_download(), sanitize=sanitize, log=lambda *a: None)
        self.assertEqual(st['state'], 'done'); self.assertEqual(st['sanitized'], dict(pairs=10, kept=9)); self.assertEqual(len(san), 1); self.assertEqual(self.left(), [])

    def test_failed_conversion_is_a_failure_not_a_success(self):      # grabENA.py exited 0 here
        for stage, rid in (('FastqToSam', 'ERR102'), ('MarkIlluminaAdapters', 'ERR103')):
            fm = FakeFM(tempfile.mkdtemp() + '/')
            run = fake_gatk(fail_for=[rid]) if stage == 'FastqToSam' else (lambda cmd, capture_output=True, _g=fake_gatk(): types.SimpleNamespace(returncode=1, stderr=b'adapter boom') if 'MarkIlluminaAdapters' in cmd else _g(cmd))
            spec = make_spec(fm.localTempDir.replace('temp/', ''), rid)
            st = W.process_run(spec, fm, run=run, download=fake_download(), log=lambda *a: None)
            self.assertEqual((st['state'], st['stage'], st['kind']), ('failed', stage, 'permanent')); self.assertEqual(fm.uploaded, [])
            self.assertEqual([f for f in os.listdir(fm.localTempDir) if not f.endswith('.json')], [])    # no leftovers

    def test_exit_codes(self):
        self.assertEqual(W.EXIT, {'permanent': 10, 'transient': 11, 'throttle': 12})

    def test_non_illumina_platform_skips_adapter_marking_but_uploads(self):
        log = []
        st = W.process_run(make_spec(self.tmp, 'ERR104', platform='DNBSEQ'), self.fm, run=fake_gatk(log=log), download=fake_download(), log=lambda *a: None)
        self.assertEqual(st['state'], 'done'); self.assertEqual([c[0] for c in log], ['FastqToSam']); self.assertEqual(len(self.fm.uploaded), 2)

    def test_upload_that_is_not_on_the_cloud_afterwards_fails_transiently(self):
        self.fm.drop_upload.add('ERR105.unmapped_marked_adapters.bam')
        st = W.process_run(make_spec(self.tmp, 'ERR105'), self.fm, run=fake_gatk(), download=fake_download(), log=lambda *a: None)
        self.assertEqual((st['state'], st['stage'], st['kind']), ('failed', 'upload', 'transient'))

    def test_fetch_failures_carry_their_kind_and_keep_or_drop_partials(self):
        for kind in ('throttle', 'permanent', 'transient'):
            spec = make_spec(self.tmp, 'ERR106'); open(self.fm.localTempDir + 'ERR106_1.fastq.gz.partial', 'w').write('p')
            st = W.fetch_stage(spec, download=fake_download(('ERR106', kind)), log=lambda *a: None)
            self.assertEqual((st['state'], st['kind']), ('failed', kind))
            self.assertEqual(os.path.exists(self.fm.localTempDir + 'ERR106_1.fastq.gz.partial'), kind != 'permanent')

    def test_convert_without_fastqs_says_so(self):
        st = W.convert_stage(make_spec(self.tmp, 'ERR107'), self.fm, run=fake_gatk(), log=lambda *a: None)
        self.assertEqual((st['state'], st['kind']), ('failed', 'transient'))


class FakeProc:
    """A worker that 'runs' for n polls, then writes its status file and returns its exit code."""
    def __init__(self, spec_path, polls, outcome):
        self.spec = json.load(open(spec_path)); self.polls, self.outcome, self.n = polls, outcome, 0

    def poll(self):
        self.n += 1
        if self.n <= self.polls:
            return None
        st = dict(run=self.spec['run'], state=self.outcome['state'], stage=self.outcome.get('stage', 'x'), kind=self.outcome.get('kind', ''), error=self.outcome.get('error', ''),
                  download=dict(method='aspera', mb_per_s=50.0))
        if self.outcome.get('nostatus'):
            return 1
        json.dump(st, open(self.spec['status_json'], 'w'))
        return {'': 0, 'permanent': 10, 'transient': 11, 'throttle': 12}[st['kind']]


class TestScheduler(unittest.TestCase):
    def setup(self, n=10, **over):
        self.tmp = tempfile.mkdtemp() + '/'; os.makedirs(self.tmp + 'temp/')
        self.cfg = dict(temp_dir=self.tmp + 'temp/', dl_max=3, workers=2, max_ready=2, stagger_s=5, max_attempts=3, backoff_s=30, min_free_gb=0, disk_factor=4, status_every_s=1e9, poll_s=1); self.cfg.update(over)
        self.tasks = [dict(run=f'ERR{i:03d}', sample='S', fastq_gb=1.0, spec=dict(run=f'ERR{i:03d}')) for i in range(n)]
        self.t = [0.0]; self.now = lambda: self.t[0]
        def sleep(dt): self.t[0] += dt
        self.sleep = sleep; self.peak = dict(fetch=0, convert=0, ready=0); self.starts = []; self.running = dict(fetch=0, convert=0)

    def launcher(self, outcomes=None, fetch_polls=3, convert_polls=4):
        outcomes = outcomes or {}
        active = {'fetch': [], 'convert': []}
        def launch(kind, spec_path, logf):
            spec = json.load(open(spec_path)); run = spec['run']
            if kind == 'fetch':
                self.starts.append(self.t[0])
            key = (kind, run); o = outcomes.get(key)
            if isinstance(o, list):
                o = o.pop(0) if len(o) > 1 else o[0]
            o = o or dict(state='fetched' if kind == 'fetch' else 'done')
            p = FakeProc(spec_path, fetch_polls if kind == 'fetch' else convert_polls, o)
            active[kind].append(p); self.peak[kind] = max(self.peak[kind], sum(1 for q in active[kind] if q.n <= q.polls))
            return p
        return launch

    def run_(self, launch, free=lambda p: 1e12):
        self.done = []
        return SCH.run_tasks(self.tasks, self.cfg, launch, lambda t, st: self.done.append(t['run']), sleep=self.sleep, now=self.now, free_fn=free, log=lambda *a: None)

    def test_all_runs_complete_within_the_pool_limits(self):
        self.setup(10); res = self.run_(self.launcher())
        self.assertEqual(sorted(self.done), [t['run'] for t in self.tasks]); self.assertTrue(all(v['state'] == 'done' for v in res.values()))
        self.assertLessEqual(self.peak['fetch'], 3); self.assertLessEqual(self.peak['convert'], 2)
        self.assertTrue(all(b - a >= 5 - 1e-9 for a, b in zip(self.starts, self.starts[1:])))     # downloads are staggered

    def test_throttle_halves_the_downloads_and_the_run_is_retried(self):
        self.setup(6); res = self.run_(self.launcher({('fetch', 'ERR000'): [dict(state='failed', stage='download', kind='throttle', error='Too many'), dict(state='fetched')]}))
        self.assertEqual(sorted(self.done), [t['run'] for t in self.tasks]); self.assertEqual(res['ERR000']['attempts'], 2)

    def test_permanent_failures_are_not_retried_and_never_reach_the_database(self):
        self.setup(4); res = self.run_(self.launcher({('fetch', 'ERR001'): dict(state='failed', stage='download', kind='permanent', error='gone'),
                                                       ('convert', 'ERR002'): dict(state='failed', stage='FastqToSam', kind='permanent', error='bad')}))
        self.assertEqual(sorted(self.done), ['ERR000', 'ERR003']); self.assertEqual((res['ERR001']['state'], res['ERR002']['state']), ('failed', 'failed'))

    def test_transient_failures_retry_then_give_up(self):
        self.setup(2, max_attempts=3); res = self.run_(self.launcher({('fetch', 'ERR000'): dict(state='failed', stage='download', kind='transient', error='net')}))
        self.assertEqual((res['ERR000']['state'], res['ERR000']['attempts']), ('failed', 3)); self.assertEqual(self.done, ['ERR001'])

    def test_a_crashed_worker_without_a_status_file_is_recorded(self):
        self.setup(1, max_attempts=1); res = self.run_(self.launcher({('fetch', 'ERR000'): dict(state='x', nostatus=True)}))
        self.assertEqual((res['ERR000']['state'], res['ERR000']['stage']), ('failed', 'worker crashed'))

    def test_backpressure_and_disk_admission(self):
        self.setup(8, workers=1, max_ready=1, dl_max=4); res = self.run_(self.launcher(convert_polls=30))
        self.assertEqual(len(self.done), 8)
        self.setup(2); res = self.run_(self.launcher(), free=lambda p: 3e9)          # needs 4 GB + margin, only 3 GB free and nothing in flight
        self.assertTrue(all(v['stage'] == 'disk' for v in res.values())); self.assertEqual(self.done, [])

    def test_a_database_failure_after_upload_is_reported(self):
        self.setup(1)
        def bad(t, st): raise RuntimeError('sheet locked')
        res = SCH.run_tasks(self.tasks, self.cfg, self.launcher(), bad, sleep=self.sleep, now=self.now, free_fn=lambda p: 1e12, log=lambda *a: None)
        self.assertEqual((res['ERR000']['state'], res['ERR000']['stage']), ('failed', 'database'))


class TestFileManagerFailures(unittest.TestCase):
    def test_upload_exceptions_are_transient(self):
        import subprocess as sp
        fm = FakeFM(tempfile.mkdtemp() + '/')
        def boom(path): raise sp.CalledProcessError(1, ['rclone', 'check'], stderr='differences found')
        fm.uploadData = boom
        st = W.process_run(make_spec(fm.localTempDir.replace('temp/', ''), 'ERR200'), fm, run=fake_gatk(), download=fake_download(), log=lambda *a: None)
        self.assertEqual((st['state'], st['stage'], st['kind']), ('failed', 'upload', 'transient'))

    def run_one(self, fm):
        ena = pd.DataFrame([rec('ERR9100001', study_accession='PRJEB78', sample_accession='SAMEA9', library_name='')])
        def launch(kind, spec_path, logf):
            spec = json.load(open(spec_path))
            st = W.fetch_stage(spec, download=fake_download(), log=lambda *a: None) if kind == 'fetch' else W.convert_stage(spec, fm, run=fake_gatk(), log=lambda *a: None)
            return types.SimpleNamespace(poll=lambda: 0 if st['state'] in ('done', 'fetched') else W.EXIT[st['kind']])
        clock = [0.0]; logs = []
        rc = D.run(D.parse_args(['--ena', 'PRJEB78', '--stagger', '0', '--ascp', '/no/ascp']), fm, fetch=lambda acc: ena, launch=launch, sleep=lambda s: clock.__setitem__(0, clock[0] + s),
                   now=lambda: clock[0], free_fn=lambda p: 1e12, log=logs.append)
        return rc, logs

    def test_a_failing_database_write_keeps_the_row_pending(self):
        fm = FakeFM(tempfile.mkdtemp() + '/')
        def bad(rd, sd): raise AssertionError('run already in database')
        fm.addDNAReadRow = bad
        rc, logs = self.run_one(fm)
        self.assertEqual(rc, 0)                                                      # the uBAM is on the cloud: the run is not a failure
        self.assertEqual(len(open(fm.localTempDir + 'pending_db_rows.jsonl').read().splitlines()), 1)
        self.assertTrue(any('--add-pending' in str(m) for m in logs))
        st = pd.read_csv([fm.localTempDir + f for f in os.listdir(fm.localTempDir) if f.startswith('download_status_')][0], sep='\t', keep_default_na=False)
        self.assertEqual(st.stage.iloc[0], 'done (row pending)')

    def test_a_failing_size_lookup_goes_pending_and_add_pending_fills_it(self):
        fm = FakeFM(tempfile.mkdtemp() + '/')
        def nosize(p): raise ValueError('not listed yet')
        fm.returnFileSize = nosize
        rc, logs = self.run_one(fm)
        self.assertEqual((rc, fm.rows), (0, []))
        row = json.loads(open(fm.localTempDir + 'pending_db_rows.jsonl').readline()); self.assertEqual(row['read_data']['FileSize'], '')
        fm.returnFileSize = lambda p: 4242
        D.add_pending(fm, fm.localTempDir + 'pending_db_rows.jsonl', log=lambda *a: None)
        self.assertEqual(fm.rows[0][0]['FileSize'], 4242)


def csv_to_ena(csv):
    return pd.DataFrame([rec(r.Run, sample_accession=r.BioSample, study_accession=r.BioProject, experiment_accession=r.Experiment, instrument_model=r.Instrument, scientific_name=r.Organism,
                             read_count=str(int(int(r.Bases) // int(r.AvgSpotLen))), base_count=r.Bases, library_name='', size=2_000_000_000) for _, r in csv.iterrows()])


class TestTwoRunInfoFiles(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp() + '/'; self.fm = FakeFM(self.tmp)
        self.csv = pd.read_csv(os.path.join(HERE, 'MalinskyPaperData.csv'), dtype=str, keep_default_na=False).head(40)
        self.runs = list(self.csv.Run); self.ena = csv_to_ena(self.csv); self.calls = []
        def fetch(acc):
            self.calls.append(acc)
            return self.ena[(self.ena.study_accession == acc) | (self.ena.run_accession == acc)].reset_index(drop=True)
        self.fetch = fetch

    def put_cloud(self, runs):
        proj = dict(zip(self.csv.Run, self.csv.BioProject))
        for r in runs:
            self.fm.cloud.setdefault(self.fm.localReadsDir + proj[r] + '/', set()).add(r + '.unmapped_marked_adapters.bam')

    def no_launch(self, *a):
        raise AssertionError('nothing should be downloaded')

    def test_two_files_in_database_on_cloud_and_to_download(self):
        self.csv.iloc[:25].to_csv(self.fm.localReadDownloadDir + 'blumer.csv', index=False); self.csv.iloc[15:].to_csv(self.fm.localReadDownloadDir + 'malinsky.csv', index=False)   # 10 runs in both
        in_db, on_cloud = self.runs[0:3], self.runs[3:5]
        self.fm.reads_dt = pd.DataFrame(dict(RunID=in_db)); self.put_cloud(on_cloud)
        logs = []
        rc = D.run(D.parse_args(['blumer.csv', 'malinsky.csv', '--dry-run']), self.fm, fetch=self.fetch, launch=self.no_launch, log=logs.append)
        self.assertEqual(rc, 0)
        plan = pd.read_csv([self.fm.localTempDir + f for f in os.listdir(self.fm.localTempDir) if f.startswith('ena_plan_2')][0], sep='\t', keep_default_na=False)
        self.assertEqual(len(plan), 40)                                                         # the 10 shared runs appear once
        st = plan.set_index('run').status
        self.assertTrue(all(st[r] == 'in_database' for r in in_db) and all(st[r] == 'on_cloud' for r in on_cloud)); self.assertEqual(int((plan.status == 'ok').sum()), 35)
        self.assertTrue(any('to download: 35 runs' in str(m) and '2 runs have a uBAM on the cloud' in str(m) and '--adopt-cloud' in str(m) for m in logs))
        self.assertEqual(len(self.calls), len(set(self.calls)))                                 # each ENA accession asked once, although both files list the same projects
        ss = pd.read_csv([self.fm.localTempDir + f for f in os.listdir(self.fm.localTempDir) if f.startswith('ena_plan_samples_')][0], sep='\t', keep_default_na=False)
        exp_complete = {s for s, g in self.csv.groupby('BioSample') if set(g.Run) <= set(in_db)}
        self.assertEqual(set(ss[ss.state == 'complete']['sample']), exp_complete); self.assertEqual(int(ss.runs.sum()), 40)
        self.assertEqual(self.fm.rows, [])                                                      # a dry run writes nothing

    def test_adopt_cloud_adds_rows_for_ubams_already_on_dropbox_without_downloading(self):
        small = self.csv.iloc[:6]; small.to_csv(self.fm.localReadDownloadDir + 'small.csv', index=False)
        runs = list(small.Run); self.fm.reads_dt = pd.DataFrame(dict(RunID=runs[:2])); self.put_cloud(runs[2:6])
        logs = []
        D.run(D.parse_args(['small.csv']), self.fm, fetch=self.fetch, launch=self.no_launch, log=logs.append)
        self.assertEqual(self.fm.rows, []); self.assertTrue(any('--adopt-cloud' in str(m) for m in logs))      # default: only reported
        D.run(D.parse_args(['small.csv', '--adopt-cloud']), self.fm, fetch=self.fetch, launch=self.no_launch, log=lambda *a: None)
        self.assertEqual(sorted(r[0]['RunID'] for r in self.fm.rows), sorted(runs[2:6]))
        self.assertEqual(self.fm.rows[0][0]['FileLocations'], f"{self.fm.rows[0][0]['ProjectID']}/{self.fm.rows[0][0]['RunID']}.unmapped_marked_adapters.bam")


class TestAsperaDiscovery(unittest.TestCase):
    def test_order_and_missing(self):
        have = {'/home/student/bin/ascp', '/home/student/.aspera/sdk/aspera_bypass_rsa.pem'}
        ex = lambda p: p in have
        old = os.environ.get('HOME'); os.environ['HOME'] = '/home/student'
        try:
            self.assertEqual(D.find_first(None, D.ASCP_CANDIDATES, exists=ex), '/home/student/bin/ascp')
            self.assertEqual(D.find_first('/x/ascp', D.ASCP_CANDIDATES, exists=ex), '/home/student/bin/ascp')          # an explicit path that does not exist is skipped, not fatal
            self.assertEqual(D.find_first(None, D.KEY_CANDIDATES, exists=ex), '/home/student/.aspera/sdk/aspera_bypass_rsa.pem')
            self.assertIsNone(D.find_first(None, D.ASCP_CANDIDATES, exists=lambda p: False))
        finally:
            os.environ['HOME'] = old


class TestNoDb(unittest.TestCase):
    def test_no_db_keeps_rows_aside_and_add_pending_adds_only_verified_ones(self):
        tmp = tempfile.mkdtemp() + '/'; fm = FakeFM(tmp)
        runs = ['ERR9000001', 'ERR9000002', 'ERR9000003']
        ena = pd.DataFrame([rec(r, study_accession='PRJEB77', sample_accession='SAMEA' + r[-1], library_name='') for r in runs])
        def launch(kind, spec_path, logf):
            spec = json.load(open(spec_path))
            st = W.fetch_stage(spec, download=fake_download(), log=lambda *a: None) if kind == 'fetch' else W.convert_stage(spec, fm, run=fake_gatk(), log=lambda *a: None)
            return types.SimpleNamespace(poll=lambda: 0 if st['state'] in ('done', 'fetched') else W.EXIT[st['kind']])
        clock = [0.0]
        args = D.parse_args(['--ena', 'PRJEB77', '--no-db', '--stagger', '0', '--ascp', '/no/ascp'])
        rc = D.run(args, fm, fetch=lambda acc: ena, launch=launch, sleep=lambda s: clock.__setitem__(0, clock[0] + s), now=lambda: clock[0], free_fn=lambda p: 1e12, log=lambda *a: None)
        self.assertEqual((rc, fm.rows), (0, []))                                    # the sample database was not touched
        lines = [json.loads(l) for l in open(fm.localTempDir + 'pending_db_rows.jsonl')]
        self.assertEqual(sorted(l['read_data']['RunID'] for l in lines), runs)
        self.assertEqual(len(fm.cloud[fm.localReadsDir + 'PRJEB77/']), 6)           # but the uBAMs and metrics are uploaded
        # one uBAM disappears from the cloud, one run is already in the database
        fm.cloud[fm.localReadsDir + 'PRJEB77/'].discard('ERR9000002.unmapped_marked_adapters.bam'); fm.reads_dt = pd.DataFrame(dict(RunID=['ERR9000003']))
        added, skipped, missing = D.add_pending(fm, fm.localTempDir + 'pending_db_rows.jsonl', log=lambda *a: None)
        self.assertEqual((added, skipped, missing), (1, 1, ['ERR9000002'])); self.assertEqual([r[0]['RunID'] for r in fm.rows], ['ERR9000001'])
        # --add-pending through run(): no download, no plan
        fm.rows.clear(); fm.reads_dt = pd.DataFrame(dict(RunID=[]))
        self.assertEqual(D.run(D.parse_args(['--add-pending', fm.localTempDir + 'pending_db_rows.jsonl']), fm, log=lambda *a: None), 0); self.assertEqual(len(fm.rows), 2)


class TestOrchestratorOnTheExampleCsv(unittest.TestCase):
    """Your MalinskyPaperData.csv (673 runs) through the real orchestrator, scheduler and workers, with fake ENA, fake downloads and fake GATK."""
    def test_end_to_end(self):
        tmp = tempfile.mkdtemp() + '/'; fm = FakeFM(tmp, in_db=['ERR13389714'])
        csv = pd.read_csv(os.path.join(HERE, 'MalinskyPaperData.csv'), dtype=str, keep_default_na=False).head(40)
        csv.to_csv(fm.localReadDownloadDir + 'runs.csv', index=False)
        self.assertNotIn('Library Name', csv.columns)                      # the real file has none; the old script exits on that
        runs = list(csv.Run)
        special = dict(orphan=set(runs[1:6]), nofastq=set(runs[6:8]), single=set(runs[8:9]), badgatk=runs[10], throttle=runs[11], gone=runs[13])
        recs = []
        for _, r in csv.iterrows():
            kw = dict(run_accession=r.Run, sample_accession=r.BioSample, study_accession=r.BioProject, experiment_accession=r.Experiment, instrument_model=r.Instrument, scientific_name=r.Organism,
                      read_count=str(int(int(r.Bases) // int(r.AvgSpotLen))), base_count=r.Bases, library_name='')
            if r.Run in special['nofastq']:
                recs.append(rec(r.Run, fastq_ftp='', fastq_md5='', fastq_bytes='', fastq_aspera='', submitted_ftp='ftp.sra.ebi.ac.uk/x/a.cram', submitted_format='CRAM', **kw))
            elif r.Run in special['single']:
                recs.append(rec(r.Run, files=('',), layout='SINGLE', **kw))
            else:
                recs.append(rec(r.Run, orphan_first=r.Run in special['orphan'], **kw))
        ena = pd.DataFrame(recs)
        fetch = lambda acc: ena[(ena.study_accession == acc) | (ena.run_accession == acc)].reset_index(drop=True)
        state = {'throttled': False}

        def download(sel, dest, **kw):
            run = sel['r1']['name'].split('_')[0]
            if run == special['throttle'] and not state['throttled']:
                state['throttled'] = True; raise E.EnaError('Session Stop: Too many connections', 'throttle')
            if run == special['gone']:
                raise E.EnaError('No such file', 'permanent')
            assert sel['r1']['name'].endswith('_1.fastq.gz') and sel['r2']['name'].endswith('_2.fastq.gz')     # never an orphan file
            return fake_download()(sel, dest)

        def launch(kind, spec_path, logf):
            spec = json.load(open(spec_path))
            st = W.fetch_stage(spec, download=download, log=lambda *a: None) if kind == 'fetch' else W.convert_stage(spec, fm, run=fake_gatk(fail_for=[special['badgatk']]), log=lambda *a: None)
            return types.SimpleNamespace(poll=lambda: 0 if st['state'] in ('done', 'fetched') else W.EXIT[st['kind']])

        args = D.parse_args(['runs.csv', '--stagger', '0', '--backoff', '0', '--workers', '3', '--dl-workers', '3', '--ascp', '/no/ascp'])
        clock = [0.0]
        rc = D.run(args, fm, fetch=fetch, launch=launch, sleep=lambda s: clock.__setitem__(0, clock[0] + s), now=lambda: clock[0], free_fn=lambda p: 1e12, log=lambda *a: None)
        added = {r[0]['RunID'] for r in fm.rows}
        expected_ok = set(runs) - {'ERR13389714'} - special['nofastq'] - special['single'] - {special['badgatk'], special['gone']}
        self.assertEqual(added, expected_ok)                              # database rows only for runs whose uBAM was verified on the cloud (incl. the orphan-listed and the throttled one)
        self.assertEqual(rc, 1)
        status = pd.read_csv([os.path.join(fm.localTempDir, f) for f in os.listdir(fm.localTempDir) if f.startswith('download_status_')][0], sep='\t', keep_default_na=False)
        self.assertEqual(len(status), 40)                                 # every run has a row: nothing vanishes
        by = status.set_index('run')
        self.assertEqual(by.loc['ERR13389714', 'state'], 'skipped'); self.assertTrue(all(by.loc[list(special['nofastq']), 'plan'] == 'no_fastq')); self.assertTrue(all(by.loc[list(special['single']), 'plan'] == 'single_end'))
        self.assertEqual((by.loc[special['badgatk'], 'state'], by.loc[special['badgatk'], 'stage']), ('failed', 'FastqToSam'))
        self.assertEqual(by.loc[special['gone'], 'kind'], 'permanent')
        self.assertEqual(int(by.loc[special['throttle'], 'attempts']), 2)
        rows = {r[0]['RunID']: r[0] for r in fm.rows}
        self.assertEqual(rows[runs[20]]['LibraryID'], runs[20])            # no Library Name column and none on ENA: the run id stands in
        self.assertEqual(rows[runs[20]]['ReadLength'], 250 if csv.set_index('Run').loc[runs[20], 'AvgSpotLen'] == '250' else int(csv.set_index('Run').loc[runs[20], 'AvgSpotLen']))
        # --retry takes only the retryable failures (here none are retryable that failed 3 times? the permanent ones are excluded)
        a2 = D.parse_args(['runs.csv', '--retry', [os.path.join(fm.localTempDir, f) for f in os.listdir(fm.localTempDir) if f.startswith('download_status_')][0], '--dry-run'])
        plan_log = []
        D.run(a2, fm, fetch=fetch, launch=launch, sleep=lambda s: None, now=lambda: 0.0, free_fn=lambda p: 1e12, log=plan_log.append)
        self.assertTrue(any('retrying' in str(m) for m in plan_log))


if __name__ == '__main__':
    unittest.main()


class TestBench(unittest.TestCase):
    def test_bench_measures_aggregate_and_counts_dead_sessions(self):
        import subprocess as sp
        sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'unit_scripts'))
        import ena_bench
        d = tempfile.mkdtemp() + '/'
        fake = d + 'ascp'
        open(fake, 'w').write('#!/usr/bin/env python3\nimport sys,time,os\nsrc=[a for a in sys.argv if a.startswith("era-fasp@")][0]\ndest=sys.argv[-1]\n'
                              'if "dies" in src: sys.stderr.write("Session Stop: Too many connections"); sys.exit(1)\n'
                              'f=open(os.path.join(dest,os.path.basename(src)),"wb")\nfor _ in range(100):\n    f.write(b"x"*100000); f.flush(); time.sleep(0.02)\n')
        os.chmod(fake, 0o755)
        mk = lambda n: dict(name=n, aspera='fasp:/v/' + n, md5='', bytes=1)
        r = ena_bench.bench([mk('a.gz'), mk('b.gz')], 2, 0.5, '300m', fake, 'k')
        self.assertEqual((r['concurrency'], r['died_early']), (2, 0)); self.assertGreater(r['aggregate_mb_s'], 0.5)
        r2 = ena_bench.bench([mk('a.gz'), mk('dies.gz')], 2, 0.5, '300m', fake, 'k')
        self.assertEqual(r2['died_early'], 1); self.assertIn('Too many', r2['errors'][0])
