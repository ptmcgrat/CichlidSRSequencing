"""The SRA source: unit tests of helper_modules/sra_fetch.py with the real numbers of the pilot (ERR13389718), the worker stages, the scheduler's source rotation and the downloader's source choice.
aria2c, vdb-validate and fasterq-dump are replaced by fakes that write what the real tools write."""
import gzip, json, os, sys, tempfile, types, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
from test_ena_v2 import TestScheduler, FakeFM, fake_gatk, fake_download, rec, csv_to_ena, HERE
from test_manifest_mode import Harness, sheet
from helper_modules import ena_fetch as E, ena_scheduler as SCH, grab_ena_v2 as W, sra_fetch as SF
import downloadShortReadData_v2 as D

# the first read of ERR13389718 as it came out of fasterq-dump and from ENA (pilot)
SEQ = 'TCTGCAGCTGTCCATGGACCAAAGTAATGCTGTGTTTCCTTGTTATAAAAGGCTCCTTTGGTGCCCTCGTGATACAGTGGGGACCAGAGTCTAGATTTCCTGAGGTATAATCGCGCCATTTTAAA'
QUAL = 'BBBBBFFFFFFFF///BFFFBFFBBBFFBF/BFFFFFBFFF/<FFFF<F/FF/BBF/<F<<FFFFFFFBBF<FBFF<B/<FF/FF/</<FBF<<BFFFBFBFB7<B7/77B/</BFF//BFF/F/'
ALPHA = 'BF/<7#'


def reads(n, alpha=ALPHA, name='r'):
    out = []
    for i in range(n):
        seq = (SEQ[i % 7:] + SEQ[:i % 7]); q = ''.join(alpha[(i + j) % len(alpha)] for j in range(len(seq)))
        out.append((f'{name}{i}', seq, q))
    return out


def fastq(rs, header=lambda i, n: f'@{n}'):
    return ''.join(f'{header(i, n)}\n{s}\n+\n{q}\n' for i, (n, s, q) in enumerate(rs))


class TestFacts(unittest.TestCase):
    def test_ena_read_count_counts_reads_so_spots_are_half(self):
        self.assertEqual(E.ena_spots(dict(read_count='109512218', library_layout='PAIRED')), 54756109)         # measured: fasterq-dump found 54,756,109 spots
        self.assertEqual(E.ena_spots(dict(read_count='100', library_layout='SINGLE')), 100); self.assertIsNone(E.ena_spots(dict(read_count='', library_layout='PAIRED')))

    def test_fasterq_summary_is_parsed_and_matches(self):
        c = SF.parse_counts('spots read      : 54,756,109\nreads read      : 109,512,218\nreads written   : 109,512,218\n')
        self.assertEqual(c['spots'], E.ena_spots(dict(read_count='109512218', library_layout='PAIRED')))

    def test_urls_and_commands(self):
        self.assertEqual(SF.sra_url('ERR13389718'), 'https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR13389718/ERR13389718')
        c = SF.aria2_command(SF.sra_url('ERR1'), '/d', 'ERR1.sra', 8)
        self.assertEqual((c[0], c[c.index('-x') + 1], c[c.index('-s') + 1], '-c' in c, c[c.index('-o') + 1], c[-1]), ('aria2c', '8', '8', True, 'ERR1.sra', SF.sra_url('ERR1')))
        f = SF.fasterq_command('/d/ERR1.sra', '/d', '/d/t', 6)
        self.assertEqual((f[0], '--split-3' in f, f[f.index('-e') + 1], f[-1]), ('fasterq-dump', True, '6', '/d/ERR1.sra'))


class TestQualityCheck(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp() + '/'; self.mine = fastq(reads(300), lambda i, n: f'@ERR13389718.{i + 1} {i + 1} length=125'); open(self.d + 'm.fastq', 'w').write(self.mine)

    def ena(self, text):
        return lambda url: text

    def test_identical_reads_with_different_headers_are_identical(self):
        theirs = fastq(reads(300), lambda i, n: f'@ERR13389718.{i + 1} HS30_17226:6:1101:1208:2202#1/1')
        self.assertEqual(SF.quality_check(self.d + 'm.fastq', 'u', fetch=self.ena(theirs))[0], 'identical-to-ENA')

    def test_a_different_quality_string_is_a_mismatch(self):
        rs = reads(300); n, s, q = rs[150]; rs[150] = (n, s, q[:-1] + ('#' if q[-1] != '#' else 'B'))
        v, d = SF.quality_check(self.d + 'm.fastq', 'u', fetch=self.ena(fastq(rs)))
        self.assertEqual(v, 'MISMATCH'); self.assertIn('1 quality', d)

    def test_unreachable_ena_is_skipped_and_the_six_symbol_binned_run_is_NOT_rejected(self):
        v, d = SF.quality_check(self.d + 'm.fastq', 'u', instrument='Illumina HiSeq 2500', fetch=self.ena(''))
        self.assertTrue(v.startswith('skipped')); self.assertIn('6 distinct', d)                    # the HiSeq 2500 run of the pilot has 6 symbols and is genuinely binned

    def test_only_one_or_two_symbols_looks_like_sra_lite(self):
        open(self.d + 'lite.fastq', 'w').write(fastq(reads(300, alpha='?#')))
        self.assertEqual(SF.quality_check(self.d + 'lite.fastq', 'u', fetch=self.ena(''))[0], 'simplified-qualities-suspected')
        self.assertEqual(SF.quality_check(self.d + 'lite.fastq', 'u', fetch=self.ena(fastq(reads(300))))[0], 'MISMATCH')            # and when ENA answers, the comparison decides

    def test_ena_head_decodes_a_truncated_gzip_stream(self):
        raw = gzip.compress(fastq(reads(2000)).encode())[:4000]
        text = SF.gunzip_partial(raw).decode(); self.assertGreater(len(SF.records(text, 5000)), 10)
        self.assertEqual(SF.records(fastq(reads(3)), 2)[0][0], SEQ)


class TestDownloadAndFetch(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp() + '/'; self.calls = []

    def monitor(self, rc=0, err='', stalled=False, write=True, leave_ctrl=False):
        def m(cmd, watch, **kw):
            self.calls.append(cmd[0])
            if cmd[0] == 'aria2c' and write:
                out = os.path.join(cmd[cmd.index('-d') + 1], cmd[cmd.index('-o') + 1]); open(out, 'wb').write(b'sra' * 1000)
                if leave_ctrl: open(out + '.aria2', 'w').write('x')
            return rc, err, stalled
        return m

    def test_download_outcomes(self):
        p, st = SF.download_sra('ERR1', self.d, monitor=self.monitor())
        self.assertEqual((os.path.basename(p), st['source'], st['bytes']), ('ERR1.sra', 'SRA-aws', 3000))
        for kw, kind in ((dict(rc=3), 'source_unavailable'), (dict(rc=1, err='network problem'), 'transient'), (dict(stalled=True), 'throttle'), (dict(leave_ctrl=True), 'transient'), (dict(rc=1, err='HTTP/1.1 404 Not Found'), 'source_unavailable')):
            with self.assertRaises(E.EnaError) as cm:
                SF.download_sra('ERR2', self.d, monitor=self.monitor(**kw))
            self.assertEqual(cm.exception.kind, kind, str(kw))

    def spec(self, **kw):
        return dict(run='ERR13389718', reads_expected=54756109, ena_url_r1='https://ftp.sra.ebi.ac.uk/x_1.fastq.gz', instrument='Illumina HiSeq 2500', fq_threads=4) | kw

    def run_fake(self, validate_rc=0, make_pair=True, spots='54,756,109'):
        def run(cmd, capture_output=True, text=True):
            self.calls.append(cmd[0])
            if cmd[0] == 'vdb-validate':
                return types.SimpleNamespace(returncode=validate_rc, stdout='', stderr='corrupt' if validate_rc else '')
            if cmd[0] == 'fasterq-dump':
                out = cmd[cmd.index('-O') + 1]
                if make_pair:
                    for k in (1, 2): open(f'{out}/ERR13389718_{k}.fastq', 'w').write(fastq(reads(5)))
                else:
                    open(f'{out}/ERR13389718.fastq', 'w').write(fastq(reads(5)))
                return types.SimpleNamespace(returncode=0, stdout='', stderr=f'spots read      : {spots}\nreads read      : 109,512,218\nreads written   : 109,512,218\n')
        return run

    def left(self):
        return sorted(f for f in os.listdir(self.d))

    def test_success_with_the_real_counts_and_clean_up(self):
        paths, st = SF.fetch_sra(self.spec(), self.d, run=self.run_fake(), monitor=self.monitor(), qcheck=lambda *a: ('identical-to-ENA', '1000 reads compared'))
        self.assertEqual([os.path.basename(p) for p in paths], ['ERR13389718_1.fastq', 'ERR13389718_2.fastq'])
        self.assertEqual((st['source'], st['quality_check'], st['reads_fastq']), ('SRA-aws', 'identical-to-ENA', 54756109))
        self.assertEqual(self.left(), ['ERR13389718_1.fastq', 'ERR13389718_2.fastq'])                  # the .sra file, its log and the scratch folder are gone

    def test_each_failure_is_classified_and_cleans_up(self):
        q = lambda *a: ('identical-to-ENA', '')
        cases = [(dict(validate_rc=1), 'transient'), (dict(spots='109,512,218'), 'source_unavailable'), (dict(make_pair=False), 'permanent')]
        for kw, kind in cases:
            with self.assertRaises(E.EnaError) as cm:
                SF.fetch_sra(self.spec(), self.d, run=self.run_fake(**kw), monitor=self.monitor(), qcheck=q)
            self.assertEqual(cm.exception.kind, kind, str(kw)); self.assertEqual(self.left(), [], str(kw))
        with self.assertRaises(E.EnaError) as cm:
            SF.fetch_sra(self.spec(), self.d, run=self.run_fake(), monitor=self.monitor(), qcheck=lambda *a: ('MISMATCH', '3 quality strings differ'))
        self.assertEqual((cm.exception.kind, self.left()), ('source_unavailable', []))

    def test_a_skipped_check_is_recorded_or_refused_on_request(self):
        sk = lambda *a: ('skipped (ENA unreachable)', '6 distinct quality symbols')
        paths, st = SF.fetch_sra(self.spec(), self.d, run=self.run_fake(), monitor=self.monitor(), qcheck=sk)
        self.assertTrue(st['quality_check'].startswith('skipped'))
        with self.assertRaises(E.EnaError) as cm:
            SF.fetch_sra(self.spec(require_quality_check=True), self.d, run=self.run_fake(), monitor=self.monitor(), qcheck=sk)
        self.assertEqual(cm.exception.kind, 'transient')


class TestWorkerSra(unittest.TestCase):
    def test_fetch_stage_records_the_fastq_and_convert_uses_and_deletes_them(self):
        tmp = tempfile.mkdtemp() + '/'; fm = FakeFM(tmp); os.makedirs(fm.localTempDir, exist_ok=True)
        spec = dict(run='ERR9300001', sample='S', library='lib', platform='ILLUMINA', sel=E.select_fastqs(rec('ERR9300001')), temp_dir=fm.localTempDir, source='sra',
                    output_bam=fm.localReadsDir + 'P/ERR9300001.unmapped_marked_adapters.bam', status_json=fm.localTempDir + 'x.status.json')
        def fake_sra(sp, dest, log=print):
            ps = [f'{dest}/ERR9300001_1.fastq', f'{dest}/ERR9300001_2.fastq']; [open(p, 'w').write('x') for p in ps]
            open(f'{dest}/ERR9300001.sra', 'w').write('sra'); os.makedirs(f'{dest}/ERR9300001.fqtmp', exist_ok=True)
            return ps, dict(source='SRA-aws', method='aria2', seconds=9.0, bytes=3, mb_per_s=17.0, quality_check='identical-to-ENA', reads_fastq=7)
        st = W.fetch_stage(spec, sra=fake_sra, log=lambda *a: None)
        self.assertEqual((st['state'], st['download']['source'], len(st['fastq'])), ('fetched', 'SRA-aws', 2))
        seen = []
        def gatk(cmd, capture_output=True):
            if 'FastqToSam' in cmd: seen.append(cmd[cmd.index('--FASTQ') + 1])
            return fake_gatk()(cmd, capture_output)
        done = W.convert_stage(dict(spec, fastq=st['fastq'], source=None), fm, run=gatk, log=lambda *a: None)
        self.assertEqual(done['state'], 'done'); self.assertTrue(seen[0].endswith('ERR9300001_1.fastq'))                       # the uncompressed FASTQs of the SRA route
        self.assertEqual([f for f in os.listdir(fm.localTempDir) if not f.endswith('.json')], [])                              # FASTQs, .sra and the scratch folder deleted

    def test_a_failed_sra_fetch_carries_its_kind(self):
        spec = dict(run='ERR9300002', temp_dir=tempfile.mkdtemp() + '/', source='sra', sel=E.select_fastqs(rec('ERR9300002')))
        def bad(sp, dest, log=print): raise E.EnaError('no such object', 'source_unavailable')
        st = W.fetch_stage(spec, sra=bad, log=lambda *a: None)
        self.assertEqual((st['state'], st['kind']), ('failed', 'source_unavailable'))


class TestSourceRotation(TestScheduler):
    def run_sources(self, sources, outcomes):
        self.setup(3, stagger_s=0)
        for t in self.tasks: t['sources'] = sources
        seen = []
        base = self.launcher(outcomes)
        def launch(kind, spec_path, logf):
            spec = json.load(open(spec_path)); seen.append((kind, spec['run'], spec.get('source'), spec.get('fastq')))
            return base(kind, spec_path, logf)
        res = self.run_(launch); return res, seen

    def test_a_run_the_first_source_cannot_deliver_goes_to_the_next_at_once(self):
        res, seen = self.run_sources(['sra', 'ena'], {('fetch', 'ERR000'): [dict(state='failed', stage='FastqToSam', kind='source_unavailable', error='SRA Lite?'), dict(state='fetched')]})
        self.assertEqual([s[2] for s in seen if s[:2] == ('fetch', 'ERR000')], ['sra', 'ena']); self.assertEqual(res['ERR000']['attempts'], 2)
        self.assertTrue(all(s[2] == 'sra' for s in seen if s[0] == 'fetch' and s[1] != 'ERR000')); self.assertEqual(len(self.done), 3)

    def test_with_one_source_unavailable_is_final(self):
        res, seen = self.run_sources(['sra'], {('fetch', 'ERR000'): dict(state='failed', stage='x', kind='source_unavailable', error='gone')})
        self.assertEqual((res['ERR000']['state'], res['ERR000']['attempts']), ('failed', 1)); self.assertEqual(len([s for s in seen if s[:2] == ('fetch', 'ERR000')]), 1)

    def test_the_fastq_files_of_the_fetch_reach_the_convert_spec(self):
        res, seen = self.run_sources(['sra'], {('fetch', 'ERR001'): dict(state='fetched', fastq=['/t/a_1.fastq', '/t/a_2.fastq'])})
        self.assertEqual([s[3] for s in seen if s[:2] == ('convert', 'ERR001')], [['/t/a_1.fastq', '/t/a_2.fastq']])


class TestOrchestratorSources(unittest.TestCase):
    def go(self, extra, probe=(False, 0.0, 0.0), ascp_exists=False, launch_none=False, which=lambda t: '/bin/' + t):
        tmp = tempfile.mkdtemp() + '/'; fm = FakeFM(tmp)
        ena = pd.DataFrame([rec('ERR9400001', study_accession='PRJEB80', sample_accession='SAMEA8', library_name='', read_count='200', library_layout='PAIRED')])
        specs, logs = [], []
        def launch(kind, spec_path, logf):
            spec = json.load(open(spec_path)); specs.append(spec)
            st = W.fetch_stage(spec, download=fake_download(), sra=lambda sp, d, log=print: ([d + '/ERR9400001_1.fastq', d + '/ERR9400001_2.fastq'], dict(source='SRA-aws', quality_check='identical-to-ENA')) if
                               [open(d + f'/ERR9400001_{k}.fastq', 'w').write('x') for k in (1, 2)] else None, log=lambda *a: None) if kind == 'fetch' else W.convert_stage(spec, fm, run=fake_gatk(), log=lambda *a: None)
            return types.SimpleNamespace(poll=lambda: 0 if st['state'] in ('done', 'fetched') else W.EXIT[st['kind']])
        argv = ['--ena', 'PRJEB80', '--no-db', '--stagger', '0'] + extra
        if ascp_exists:
            a = tmp + 'ascp'; k = tmp + 'k.pem'; open(a, 'w').write('x'); open(k, 'w').write('k'); argv += ['--ascp', a, '--key', k]
        else:
            argv += ['--ascp', '/no/ascp']
        args = D.parse_args(argv); clock = [0.0]
        rc = D.run(args, fm, fetch=lambda acc: ena, launch=None if launch_none else launch, sleep=lambda s: clock.__setitem__(0, clock[0] + s), now=lambda: clock[0], free_fn=lambda p: 1e12, log=logs.append,
                   probe=lambda *a, **k: probe, which=which)
        return args, specs, logs, rc

    def test_auto_prefers_sra_when_aspera_moves_no_data_and_ena_when_it_does(self):
        args, specs, logs, rc = self.go([], probe=(False, 0.0, 0.0), ascp_exists=True)
        self.assertEqual(specs[0]['source'], 'sra'); self.assertTrue(any("['sra', 'ena']" in str(m) for m in logs)); self.assertEqual((rc, args.disk_factor), (0, 8.0))
        args, specs, logs, rc = self.go([], probe=(True, 400.0, 20.0), ascp_exists=True)
        self.assertEqual(specs[0]['source'], 'ena'); self.assertTrue(any("['ena', 'sra']" in str(m) for m in logs)); self.assertEqual(args.disk_factor, 8.0)      # SRA stays the fallback: its space is reserved

    def test_explicit_sources(self):
        self.assertEqual(self.go(['--source', 'sra'])[1][0]['source'], 'sra'); self.assertEqual(self.go(['--source', 'ena'])[1][0]['source'], 'ena')

    def test_the_spec_carries_what_the_sra_route_checks_against(self):
        spec = self.go(['--source', 'sra', '--fq-threads', '5', '--sra-connections', '6', '--sra-require-check'])[1][0]
        self.assertEqual((spec['reads_expected'], spec['fq_threads'], spec['sra_connections'], spec['require_quality_check']), (100, 5, 6, True))   # read_count 200 = 100 pairs
        self.assertTrue(spec['ena_url_r1'].startswith('https://ftp.sra.ebi.ac.uk/'))

    def test_missing_sra_tools_are_reported_before_anything_starts(self):
        with self.assertRaises(SystemExit) as cm:
            self.go(['--source', 'sra'], launch_none=True, which=lambda t: None)
        self.assertIn('aria2c', str(cm.exception))

    def test_max_runs_takes_the_smallest_runs_only(self):
        tmp = tempfile.mkdtemp() + '/'; fm = FakeFM(tmp)
        ena = pd.DataFrame([rec(f'ERR95000{i}', study_accession='PRJEB81', sample_accession=f'SAMEA9{i}', library_name='', size=1000 * (5 - i)) for i in range(5)])
        specs = []
        def launch(kind, spec_path, logf):
            specs.append(json.load(open(spec_path))['run']); raise SystemExit('stop')
        args = D.parse_args(['--ena', 'PRJEB81', '--source', 'ena', '--max-runs', '2', '--ascp', '/no/ascp', '--dry-run']); logs = []
        D.run(args, fm, fetch=lambda acc: ena, launch=launch, free_fn=lambda p: 1e12, log=logs.append)
        self.assertTrue(any("--max-runs 2" in str(m) and 'ERR950004' in str(m) and 'ERR950003' in str(m) for m in logs))

    def test_main_accepts_the_manifest_option(self):
        D.validate_args(D.parse_args(['--manifest', 'Blumer', 'Malinsky', '--dry-run']))
        with self.assertRaises(SystemExit):
            D.validate_args(D.parse_args(['--dry-run']))


class TestManifestWithSra(Harness):
    def fake_sra(self, spec, dest, log=print):
        if spec['run'] == self.runs[5] and not getattr(self, 'tried', False):
            self.tried = True; raise E.EnaError('SRA Lite?', 'source_unavailable')
        ps = [f"{dest}/{spec['run']}_{k}.fastq" for k in (1, 2)]; [open(p, 'w').write('x') for p in ps]
        return ps, dict(source='SRA-aws', method='aria2', seconds=10.0, bytes=5_000_000_000, mb_per_s=17.0, quality_check='identical-to-ENA', reads_fastq=1234)

    def launch(self, kind, spec_path, logf):
        spec = json.load(open(spec_path))
        st = W.fetch_stage(spec, download=self.download, sra=self.fake_sra, log=lambda *a: None) if kind == 'fetch' else W.convert_stage(spec, self.fm, run=fake_gatk(fail_for=[self.runs[4]]), log=lambda *a: None)
        return types.SimpleNamespace(poll=lambda: 0 if st['state'] in ('done', 'fetched') else W.EXIT[st['kind']])

    def test_source_quality_check_and_read_counts_land_in_the_sheet(self):
        self.go('--source', 'auto'); s = sheet(self.fm); r = self.runs
        row = s.loc[r[8]]
        self.assertEqual((row.DownloadStatus, row.Source, row.QualityCheck, row.ReadsInUBAM, row.Verified), ('in_database', 'SRA-aws', 'identical-to-ENA', '1234', 'TRUE'))
        self.assertEqual(float(row.MBPerSec), 17.0)
        fb = s.loc[r[5]]                                                    # attempt 1: SRA cannot deliver it correctly; attempt 2: ENA (throttled by the shared harness); attempt 3: SRA again, now fine
        self.assertEqual((fb.DownloadStatus, fb.Source, fb.Attempts), ('in_database', 'SRA-aws', '3')); self.assertEqual(fb.QualityCheck, 'identical-to-ENA')
        self.assertEqual(s.loc[r[4], 'DownloadStatus'], 'failed')
        self.assertTrue(s.loc[r[8], 'ReadsSource'] != '' and int(s.loc[r[8], 'ReadsSource']) > 0)


if __name__ == '__main__':
    unittest.main()
