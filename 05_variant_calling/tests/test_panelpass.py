"""variant_calling panel (panelpass.py, panelpass_cli.py): module 04's panel over existing BAMs and its calls into the SampleDatabase, with injected fakes for Dropbox (a folder), rclone, samtools, the panel
subprocess and module 04's call layer. The call layer is a stand-in written from the documented interface (decision 76: call_panel, sheet_row, a record with `calls`); the real one is tested by module 04."""
import io, json, os, shutil, subprocess, sys, tempfile, time, types, unittest
from unittest import mock
import contextlib, importlib
import pandas as pd
from cichlid_common.fake_sheet import FakeWorksheet
from variant_calling import build as B, panelpass as P, panelpass_cli as C

G = 'Mzebra_GT3_NCBI'
KEYS_TEXT = '[panel] gsdf     OK       5d71ee8d6c2b \n[panel] lg5inv   OK       80bbea57e10b \n[panel] bchrom   OK       75367d0dfba8 \n[panel] te       OK       2af781e90119 '
KEYS_STRING = 'bchrom=75367d0dfba8;gsdf=5d71ee8d6c2b;lg5inv=80bbea57e10b;te=2af781e90119'


class DefinitionMismatch(Exception):
    pass


class FakeCalls:
    """What module 04's sd_genotypes.panel.calls is documented to do, on the fake JSONs below."""
    def call_panel(self, path):
        d = json.load(open(path))
        if d.get('refuse'):
            raise DefinitionMismatch('measured under another definition than the current one')
        calls = [dict(c, value_to_write=(c['value'] if c['confident'] else None), call_tag='pilot-v1', definition='x') for c in d['fake_calls']]
        return dict(bam=d['bam'], bam_size_bytes=d['bam_size_bytes'], calls=calls)

    def sheet_row(self, r):
        values = {c['column']: c['value_to_write'] for c in r['calls'] if c['confident'] and c['value_to_write'] is not None}
        flags = {c['flags_column']: ';'.join(['pilot-v1'] + list(c['flags'])) for c in r['calls']}
        return values, flags


def fake_calls(n):
    return [dict(column='gsdf number', flags_column='gsdf flags', value=2 + (n % 2), flags=[], confident=True, details={'ratio': round(0.5 + 0.1 * n, 2)}),
            dict(column='chr7ins', flags_column='chr7ins flags', value='0/0', flags=['dilution_possible'] if n == 4 else [], confident=True, details={'ref': 20, 'alt': 0}),
            dict(column='TE1', flags_column='TE1 flags', value=None if n == 2 else '0/0', flags=['no_call'] if n == 2 else [], confident=(n != 2), details={'ref': 5 * n, 'alt': 0}),
            dict(column='TE2', flags_column='TE2 flags', value='0/1', flags=[], confident=True, details={}),
            dict(column='Inversion5', flags_column='Inversion5 flags', value=1, flags=['near_threshold'], confident=True, details={'median_site_fraction': 0.5}),
            dict(column='B score', flags_column='B flags', value=1.07, flags=[], confident=True, details={'blocks': 836})]


class PassBase(unittest.TestCase):
    ids = ['S1', 'S2', 'S3', 'S4']
    sizes = {'S1': 1000, 'S2': 2000, 'S3': 1500, 'S4': 500}
    delay = 0.0                                                 # seconds a fake download takes

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.local, self.cloud = self.tmp + '/local/', self.tmp + '/cloud/'
        for s in self.ids:
            d = f'{self.cloud}Bamfiles/{G}/{s}/'
            os.makedirs(d)
            open(d + f'{s}.all.bam', 'wb').write(b'B' * self.sizes[s]); open(d + f'{s}.all.bai', 'wb').write(b'I' * 10)
        self.env = B.Env(self.local, self.cloud, G, self.tmp + '/ref.fa')
        self.pcfg = P.load_config(); self.pcfg['default_bam_bytes'] = 1000.0
        self.panel = P.make_panel_config(self.pcfg, python='python-fake')
        self.ledger = pd.DataFrame({'SampleID': self.ids, 'GenomeVersion': G, 'BamSize': [self.sizes[s] for s in self.ids]})
        self.fail_on, self.bad_size, self.refuse, self.timeout, self.exit_after_json = set(), set(), set(), set(), set()
        self.panel_ok, self.bad_upload, self.log = True, False, []

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    # the fakes
    def fake_run(self, argv, input_text=None):
        self.log.append(list(argv))
        if argv[:2] == ['rclone', 'copy']:
            src, dst = argv[2], argv[3]
            os.makedirs(dst, exist_ok=True)
            time.sleep(self.delay)
            for n in [argv[i + 1] for i, x in enumerate(argv) if x == '--include']:
                if os.path.exists(src + n):
                    shutil.copy(src + n, dst + n)
            return 0, '', ''
        if argv[:2] == ['rclone', 'copyto']:
            os.makedirs(os.path.dirname(argv[3]), exist_ok=True)
            shutil.copy(argv[2], argv[3])
            return 0, '', ''
        if argv[:2] == ['rclone', 'lsf']:
            dest = argv[-1]
            return (0, f'{os.path.getsize(dest) + (1 if self.bad_upload else 0)}\n', '') if os.path.exists(dest) else (3, '', 'directory not found')
        if argv[:2] == ['samtools', 'quickcheck']:
            return 0, '', ''
        return 127, '', 'unknown tool ' + argv[0]

    def fake_runner(self, cmd, **kw):
        i = cmd.index('sd_genotypes')
        verb = cmd[i + 2]
        if verb == 'check':
            return types.SimpleNamespace(returncode=0 if self.panel_ok else 1, stdout=KEYS_TEXT if self.panel_ok else 'config/gsdf.yaml is missing')
        bam, sample, out = cmd[cmd.index('--bam') + 1], cmd[cmd.index('--sample') + 1], cmd[cmd.index('--out-json') + 1]
        if sample in self.timeout:
            raise subprocess.TimeoutExpired(cmd, 1)
        if sample in self.fail_on:
            return types.SimpleNamespace(returncode=1, stdout='Traceback ... site_error')
        json.dump({'schema': 2, 'bam': os.path.basename(bam), 'bam_size_bytes': os.path.getsize(bam) + (5 if sample in self.bad_size else 0), 'sample': sample,
                   'fake_calls': fake_calls(int(sample[1:])), 'refuse': sample in self.refuse}, open(out, 'w'))
        return types.SimpleNamespace(returncode=2 if sample in self.exit_after_json else 0, stdout='warning: it wrote its JSON and then died')

    def plan(self, members=None, **kw):
        keys = kw.pop('keys', P.parse_keys(KEYS_TEXT))
        return P.make_plan(self.env, self.pcfg, members or self.ids, self.ledger, keys, **kw)

    def measure(self, s, plan=None, **kw):
        plan = plan or self.plan()
        return P.measure_sample(self.env, self.pcfg, self.panel, s, plan, self.sizes[s], run=self.fake_run, runner=self.fake_runner, free_bytes=kw.pop('free_bytes', 10 ** 12), **kw)


class TestConfigAndKeys(unittest.TestCase):
    def test_the_shipped_config(self):
        c = P.load_config()
        self.assertEqual(c['expect_keys'], ['5d71ee8d6c2b', '2af781e90119', '80bbea57e10b', '75367d0dfba8']); self.assertTrue(all(isinstance(k, str) for k in c['expect_keys']))
        self.assertEqual(c['sheet']['write_columns'], ['gsdf number', 'chr7ins', 'TE1', 'TE2', 'Inversion5', 'B score'], 'module 04 lifted the hold on 2026-11-17 (decision 78)')
        self.assertEqual(c['sheet']['scope'], {'Inversion5': ['Mbuna'], 'B score': ['Mbuna']}, 'LG5 and B were validated on the mbuna only (decision 78)')
        self.assertEqual((c['sheet']['nocall']['value'], c['sheet']['nocall']['columns']), ('./.', ['chr7ins', 'TE1', 'TE2', 'Inversion5'])); self.assertIn('unvalidated_value', c['sheet']['nocall']['exclude_flags']); self.assertIn('site_error', c['sheet']['nocall']['exclude_flags'])
        self.assertTrue(c['sheet']['flags_in_sheet'])
        self.assertEqual((c['jobs'], c['default_bam_bytes'], c['min_free_factor']), (2, 12e9, 1.15)); self.assertEqual(c['sheet']['tab'], 'SampleDatabase')

    def test_bad_numbers_and_an_empty_write_list_are_refused(self):
        d = tempfile.mkdtemp()
        try:
            for text, frag in (('jobs: many\n', 'jobs must be a number'), ('default_bam_bytes: lots\n', 'default_bam_bytes'), ('sheet: {write_columns: []}\n', 'write_columns is empty')):
                open(d + '/c.yaml', 'w').write(text)
                with self.assertRaises(ValueError) as cm: P.load_config(d + '/c.yaml')
                self.assertIn(frag, str(cm.exception))
            open(d + '/c.yaml', 'w').write("expect_keys: [123456789012, '2af781e90119']\njobs: '3'\n")
            c = P.load_config(d + '/c.yaml'); self.assertEqual(c['expect_keys'], ['123456789012', '2af781e90119']); self.assertEqual(c['jobs'], 3)         # an all-digit key is a number to YAML
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_the_definition_keys_are_read_from_what_panel_check_prints(self):
        k = P.parse_keys(KEYS_TEXT)
        self.assertEqual(k, {'gsdf': '5d71ee8d6c2b', 'lg5inv': '80bbea57e10b', 'bchrom': '75367d0dfba8', 'te': '2af781e90119'}); self.assertEqual(P.keys_string(k), KEYS_STRING)
        self.assertEqual(P.parse_keys('[panel] gsdf     FAILED   no config\nnoise'), {}); self.assertEqual(P.parse_keys(None), {})


class TestPlan(PassBase):
    def log_row(self, s, **kw):
        row = dict(SampleID=s, GenomeVersion=G, Status='ok', BamBytes=self.sizes[s], DefinitionKeys=KEYS_STRING); row.update(kw)
        os.makedirs(self.plan().json_dir, exist_ok=True)
        P.append_row(self.plan().log_path, row)

    def with_json(self, s):
        os.makedirs(self.plan().json_dir, exist_ok=True)
        open(f'{self.plan().json_dir}{s}.panel.json', 'w').write('{}')

    def test_a_fresh_plan_has_everything_to_do_and_names_its_places(self):
        p = self.plan()
        self.assertEqual((p.todo, p.done, p.stale), (self.ids, {}, {})); self.assertTrue(p.json_dir.endswith('/local/Outputs/PanelJSON/Mzebra_GT3_NCBI/')); self.assertTrue(p.staged.endswith('/local/PanelPass/staged/'))
        self.assertTrue(p.log_path.endswith('/Outputs/PanelJSON/Mzebra_GT3_NCBI/panel_pass.Mzebra_GT3_NCBI.tsv')); self.assertEqual(p.sizes['S2'], 2000.0)

    def test_a_measured_sample_with_its_json_is_done(self):
        self.log_row('S1'); self.with_json('S1')
        p = self.plan(); self.assertEqual((list(p.done), p.todo, p.stale), (['S1'], ['S2', 'S3', 'S4'], {}))

    def test_the_four_reasons_for_measuring_again(self):
        self.log_row('S1'); self.with_json('S1'); self.log_row('S2', DefinitionKeys='gsdf=OLD'); self.with_json('S2')
        self.log_row('S3', BamBytes=1499); self.with_json('S3'); self.log_row('S4')                           # S4: no JSON in the folder
        p = self.plan()
        self.assertEqual(p.todo, ['S2', 'S3', 'S4']); self.assertIn('definition keys', p.stale['S2']); self.assertIn('realigned', p.stale['S3']); self.assertIn('not in the local folder', p.stale['S4'])

    def test_a_failed_row_and_redo_and_limit(self):
        self.log_row('S1', Status='failed', Error='boom'); self.with_json('S1'); self.log_row('S2'); self.with_json('S2')
        self.assertEqual(self.plan().todo, ['S1', 'S3', 'S4']); self.assertEqual(self.plan(redo=True).todo, self.ids)
        p = self.plan(limit=2); self.assertEqual((p.todo, p.todo_all), (['S1', 'S3'], ['S1', 'S3', 'S4']))
        self.log_row('S2', Status='failed'); self.assertIn('S2', self.plan().todo, 'the LAST row of a sample decides')

    def test_without_keys_the_keys_are_not_compared_and_duplicates_in_the_selection_are_dropped(self):
        self.log_row('S1', DefinitionKeys='anything'); self.with_json('S1')
        p = self.plan(['S1', 'S1', 'S2'], keys={}); self.assertEqual((list(p.done), p.members), (['S1'], ['S1', 'S2']))

    def test_a_sample_not_in_the_ledger_is_listed_and_its_size_is_unknown(self):
        p = P.make_plan(self.env, self.pcfg, ['S1', 'NEW'], self.ledger, P.parse_keys(KEYS_TEXT)); self.assertEqual(p.unknown, ['NEW']); self.assertIsNone(p.sizes['NEW'])
        p = P.make_plan(self.env, self.pcfg, ['S1'], None, {}); self.assertEqual(p.unknown, ['S1'])


class TestMeasureSample(PassBase):
    def test_the_whole_pass_for_one_sample(self):
        row = self.measure('S1'); p = self.plan()
        self.assertEqual((row['Status'], row['Error']), ('ok', '')); self.assertEqual((row['BamBytes'], row['JsonBytes'] > 0, row['DefinitionKeys']), (1000, True, KEYS_STRING))
        self.assertTrue(os.path.exists(f'{p.json_dir}S1.panel.json')); self.assertFalse(os.path.exists(f'{p.json_dir}S1.panel.json.part'))
        self.assertTrue(os.path.exists(f'{self.cloud}Bamfiles/{G}/S1/S1.panel.json'), 'the durable copy next to the BAM'); self.assertFalse(os.path.exists(f'{p.staged}S1/'), 'the BAM is deleted')
        self.assertTrue(all(row[k] != '' for k in ('DownloadSeconds', 'MeasureSeconds', 'UploadSeconds', 'Seconds', 'Machine', 'MeasuredOn')))
        cmds = [a for a in self.log if a[:2] == ['rclone', 'copy']]; self.assertEqual(len(cmds), 1); self.assertIn('S1.all.bam', cmds[0]); self.assertIn('S1.all.bai', cmds[0])
        self.assertEqual(json.load(open(f'{p.json_dir}S1.panel.json'))['bam_size_bytes'], 1000)

    def test_no_upload_and_keep_bam(self):
        row = self.measure('S1', upload=False, keep_bam=True); p = self.plan()
        self.assertEqual(row['Status'], 'ok'); self.assertFalse(os.path.exists(f'{self.cloud}Bamfiles/{G}/S1/S1.panel.json')); self.assertEqual(row['UploadSeconds'] if 'UploadSeconds' in row else '', '')
        self.assertTrue(os.path.exists(f'{p.staged}S1/S1.all.bam')); self.assertFalse([a for a in self.log if a[:2] == ['rclone', 'copyto']])

    def test_the_panel_failing_fails_the_sample_and_still_deletes_the_bam(self):
        self.fail_on.add('S1'); row = self.measure('S1'); p = self.plan()
        self.assertEqual(row['Status'], 'failed'); self.assertIn('failed (exit 1)', row['Error']); self.assertIn('site_error', row['Error'])
        self.assertFalse(os.path.exists(f'{p.staged}S1/')); self.assertFalse(os.path.exists(f'{p.json_dir}S1.panel.json')); self.assertFalse(os.path.exists(f'{self.cloud}Bamfiles/{G}/S1/S1.panel.json'))

    def test_a_non_zero_exit_is_a_failure_even_when_a_json_was_written(self):
        self.exit_after_json.add('S1'); row = self.measure('S1'); p = self.plan()
        self.assertEqual(row['Status'], 'failed'); self.assertIn('failed (exit 2)', row['Error']); self.assertIn('wrote its JSON and then died', row['Error'])
        self.assertFalse(os.path.exists(f'{p.json_dir}S1.panel.json')); self.assertFalse(os.path.exists(f'{p.json_dir}S1.panel.json.part')); self.assertFalse(os.path.exists(f'{self.cloud}Bamfiles/{G}/S1/S1.panel.json'))

    def test_a_timeout_is_a_failure(self):
        self.timeout.add('S1'); row = self.measure('S1'); self.assertEqual(row['Status'], 'failed'); self.assertIn('did not finish within', row['Error'])

    def test_a_json_that_names_another_bam_is_refused_and_nothing_is_kept_or_uploaded(self):
        self.bad_size.add('S1'); row = self.measure('S1'); p = self.plan()
        self.assertEqual(row['Status'], 'failed'); self.assertIn("does not belong to the BAM (different)", row['Error'])
        self.assertFalse(os.path.exists(f'{p.json_dir}S1.panel.json')); self.assertFalse(os.path.exists(f'{p.json_dir}S1.panel.json.part')); self.assertFalse(os.path.exists(f'{self.cloud}Bamfiles/{G}/S1/S1.panel.json'))

    def test_an_upload_whose_size_does_not_check_out_is_not_trusted(self):
        self.bad_upload = True; row = self.measure('S1'); self.assertEqual(row['Status'], 'failed'); self.assertIn('on Dropbox has', row['Error']); self.assertIn('not trusted', row['Error'])

    def test_a_bam_missing_on_dropbox_and_a_full_disk_are_failures(self):
        os.remove(f'{self.cloud}Bamfiles/{G}/S1/S1.all.bam'); row = self.measure('S1'); self.assertEqual(row['Status'], 'failed'); self.assertIn('is not in', row['Error'])
        row = self.measure('S2', free_bytes=100); self.assertEqual(row['Status'], 'failed'); self.assertIn('only 0 GB free', row['Error'])
        self.assertFalse([a for a in self.log if a[:2] == ['rclone', 'copy'] and 'S2.all.bam' in a], 'no download is started when the disk cannot take it')

    def test_a_stale_part_file_is_not_mistaken_for_a_result(self):
        p = self.plan(); os.makedirs(p.json_dir); open(f'{p.json_dir}S1.panel.json.part', 'w').write('old')
        self.fail_on.add('S1'); row = self.measure('S1'); self.assertEqual(row['Status'], 'failed'); self.assertFalse(os.path.exists(f'{p.json_dir}S1.panel.json.part'))


class TestRunAndResume(PassBase):
    def go(self, plan=None, **kw):
        plan = plan or self.plan(); lines = []
        res = P.run_pass(self.env, self.pcfg, plan, self.panel, log=lines.append, run=self.fake_run, runner=self.fake_runner, free_bytes=10 ** 12, **kw)
        return res, lines, plan

    def test_every_sample_is_measured_logged_and_its_bam_deleted(self):
        self.delay = 0.2; res, lines, p = self.go()
        self.assertEqual((res['ok'], res['failed'], res['skipped']), (4, [], 0)); rows = P.read_rows(p.log_path); self.assertEqual(sorted(r['SampleID'] for r in rows), self.ids)
        self.assertTrue(all(r['Status'] == 'ok' for r in rows)); self.assertEqual(P.list_json_samples(p.json_dir), self.ids); self.assertEqual(os.listdir(p.staged), [])
        self.assertIn('4 measured, 0 failed, 0 skipped', lines[-1]); self.assertIn('MB/s over the whole run', lines[-1]); self.assertIn('MB/s per download stream', lines[-1]); self.assertTrue(any('1/4' in l or '2/4' in l for l in lines))

    def test_a_download_that_took_no_time_does_not_break_the_run(self):
        """A BAM kept from an earlier run, or a download too fast to time, has a download time of 0: the log line divided by it and stopped the whole pass."""
        self.assertEqual(self.delay, 0.0); res, lines, p = self.go()
        self.assertEqual((res['ok'], res['failed']), (4, [])); self.assertTrue(all('no download' in l for l in lines if 'GB,' in l)); self.assertEqual(len(P.read_rows(p.log_path)), 4)
        self.assertNotIn('MB/s over the whole run', lines[-1], 'no rate is claimed without a download time')

    def test_a_problem_in_the_log_line_never_stops_the_run(self):
        calls = []
        def bad_log(m):
            calls.append(m)
            if 'GB,' in m: raise RuntimeError('the log line broke')
        plan = self.plan(); res = P.run_pass(self.env, self.pcfg, plan, self.panel, log=bad_log, run=self.fake_run, runner=self.fake_runner, free_bytes=10 ** 12)
        self.assertEqual(res['ok'], 4); self.assertEqual(len(P.read_rows(plan.log_path)), 4); self.assertTrue(any('could not be made' in m for m in calls))

    def test_one_failed_sample_does_not_stop_the_others_and_a_rerun_retries_only_it(self):
        self.fail_on.add('S2'); res, lines, p = self.go()
        self.assertEqual((res['ok'], [s for s, _ in res['failed']]), (3, ['S2'])); self.assertTrue(any('S2: FAILED' in l for l in lines))
        self.fail_on.clear(); p2 = self.plan(); self.assertEqual(p2.todo, ['S2']); res2, _, _ = self.go(p2)
        self.assertEqual((res2['ok'], res2['failed']), (1, [])); self.assertEqual(self.plan().todo, [], 'now everything is done'); self.assertEqual(len(P.read_rows(p.log_path)), 5)

    def test_a_run_with_nothing_to_do_says_so(self):
        self.go(); res, lines, _ = self.go(self.plan()); self.assertEqual((res['ok'], res['skipped']), (0, 4)); self.assertIn('0 measured, 0 failed, 4 skipped', lines[-1])

    def test_a_changed_definition_measures_everything_again(self):
        self.go(); res, _, _ = self.go(self.plan(keys={'gsdf': 'NEWKEY'})); self.assertEqual(res['ok'], 4)

    def test_limit_takes_the_first_samples_still_to_do_and_batches_resume(self):
        res, _, _ = self.go(self.plan(limit=2)); self.assertEqual(res['ok'], 2); self.assertEqual(self.plan(limit=2).todo, ['S3', 'S4']); res2, _, _ = self.go(self.plan(limit=2)); self.assertEqual(res2['ok'], 2)

    def test_two_jobs_at_a_time_give_the_same_result(self):
        p = self.plan(jobs=2); res, _, _ = self.go(p); self.assertEqual(res['ok'], 4); self.assertEqual(sorted(r['SampleID'] for r in P.read_rows(p.log_path)), self.ids)


class TestPreflight(PassBase):
    def pre(self, plan=None, ok=True, text=KEYS_TEXT, **kw):
        plan = plan or self.plan(); return P.preflight(self.env, self.pcfg, plan, ok, text, which=kw.pop('which', lambda t: '/bin/' + t), free_bytes=kw.pop('free_bytes', 10 ** 12)), plan

    def test_a_ready_run_has_no_problem_and_knows_what_it_needs(self):
        bad, p = self.pre(); self.assertEqual(bad, []); self.assertAlmostEqual(p.need_bytes, (2000 + 1500) * 1.15)

    def test_each_problem_is_named(self):
        self.assertIn('rclone is not on the PATH', self.pre(which=lambda t: None if t == 'rclone' else '/bin/' + t)[0][0])
        self.assertIn('the panel is not ready: config/gsdf.yaml is missing', self.pre(ok=False, text='config/gsdf.yaml is missing')[0][0])
        self.assertIn('no samples are selected', self.pre(P.make_plan(self.env, self.pcfg, [], self.ledger, {}))[0][0])
        bad, _ = self.pre(free_bytes=1000); self.assertIn('lower --jobs', bad[0]); self.assertIn('2 BAMs at a time', bad[0])

    def test_the_expected_keys_stop_a_run_through_module_3s_preflight(self):
        def runner(cmd, **kw):
            return types.SimpleNamespace(returncode=0, stdout='[panel] gsdf     OK       5d71ee8d6c2b \n[panel] te       OK       000000000000 ')
        ok, text, keys = P.check_panel(self.panel, runner=runner)
        self.assertFalse(ok); self.assertIn('not the expected ones', text); self.assertIn("'2af781e90119'", text); self.assertEqual(keys, {})
        ok, text, keys = P.check_panel(self.panel, runner=self.fake_runner); self.assertTrue(ok); self.assertEqual(keys['gsdf'], '5d71ee8d6c2b')

    def test_the_plan_text_tells_what_will_happen(self):
        p = self.plan(); text = P.format_plan(p, [], self.pcfg)
        for frag in ('4 samples (0 already measured, 4 to measure', KEYS_STRING, 'download: 0 GB in all, 2 BAMs at a time', 'copied to Dropbox next to its BAM', 'ready: add --run'): self.assertIn(frag, text)
        self.assertIn('NOT copied to Dropbox', P.format_plan(p, [], self.pcfg, upload=False)); self.assertIn('PROBLEM: x', P.format_plan(p, ['x'], self.pcfg)); self.assertIn('starting', P.format_plan(p, [], self.pcfg, running=True))
        self.assertIn('this run takes the first 2', P.format_plan(self.plan(limit=2), [], self.pcfg))


class CliBase(PassBase):
    def setUp(self):
        super().setUp()
        self.st, self.al = self.tmp + '/samples.tsv', self.tmp + '/align.tsv'
        pd.DataFrame({'SampleID': self.ids + ['NOBAM'], 'Species': ['Cm', 'Cm', 'Lf', 'Lf', 'Cm'], 'Ecogroup': ['Mbuna'] * 4 + ['Mbuna'], 'Sex': ['m', 'f', 'm', 'f', 'm']}).to_csv(self.st, sep='\t', index=False)
        self.ledger.to_csv(self.al, sep='\t', index=False)
        self.common = ['--local-master', self.local, '--cloud-master', self.cloud, '--genome', G, '--sample-table', self.st, '--alignment-table', self.al]

    def call(self, *args, **kw):
        out, err = io.StringIO(), io.StringIO()
        kw = dict(dict(run=self.fake_run, runner=self.fake_runner, free_bytes=10 ** 12, which=lambda t: '/bin/' + t, calls_module=FakeCalls()), **kw)
        rc = C.main(list(args) + self.common, out=out, err=err, **kw)
        return rc, out.getvalue(), err.getvalue()

    def jd(self):
        return f'{self.local}Outputs/PanelJSON/{G}/'

    def base_cfg(self, **sheet):
        """A config file: the SHIPPED config with the four original columns, no nocall block and no scope (and the given sheet keys set). The writer tests are about the writer; what the shipped config chooses (./. for a no call) is
        checked in test_panelpass_integration.py."""
        import yaml
        cfg = P.load_config(); cfg['sheet']['write_columns'] = ['gsdf number', 'chr7ins', 'TE1', 'TE2']; cfg['sheet']['nocall'] = {}; cfg['sheet']['scope'] = {}; cfg['sheet'].update(sheet)
        p = self.tmp + f'/cfg{len(os.listdir(self.tmp))}.yaml'
        yaml.safe_dump(cfg, open(p, 'w'))
        return p


class TestMeasureCli(CliBase):
    def test_without_run_it_prints_the_plan_and_stages_and_measures_nothing(self):
        rc, out, err = self.call('measure', '--ecogroup', 'Mbuna')
        self.assertEqual(rc, 0, err); self.assertIn('4 samples (0 already measured, 4 to measure', out); self.assertIn('ready: add --run', out); self.assertFalse(os.path.exists(self.jd()))
        self.assertIn('5 samples of Ecogroup Mbuna in the SampleDatabase', out); self.assertIn('1 selected samples have no row in the AlignmentDatabase of Mzebra_GT3_NCBI (no BAM there): left out (NOBAM)', out)
        self.assertFalse([a for a in self.log if a[:2] == ['rclone', 'copy']])

    def test_run_measures_and_a_second_run_does_nothing(self):
        rc, out, err = self.call('measure', '--ecogroup', 'Mbuna', '--run'); self.assertEqual(rc, 0, err + out)
        self.assertEqual(P.list_json_samples(self.jd()), self.ids); self.assertIn('4 measured, 0 failed', out)
        rc, out, _ = self.call('measure', '--ecogroup', 'Mbuna', '--run'); self.assertEqual(rc, 0); self.assertIn('0 to measure', out)

    def test_limit_and_the_samples_file_and_samples(self):
        f = self.tmp + '/ids.txt'; open(f, 'w').write('# the pilot\nS3\nS1  # first\n\n')
        rc, out, _ = self.call('measure', '--samples-file', f, '--run'); self.assertEqual((rc, P.list_json_samples(self.jd())), (0, ['S1', 'S3']))
        rc, out, _ = self.call('measure', '-s', 'S2', 'S4', 'S1', '--limit', '1', '--run'); self.assertEqual(P.list_json_samples(self.jd()), ['S1', 'S2', 'S3']); self.assertIn('takes the first 1', out)

    def test_a_failure_gives_exit_1_and_names_the_log(self):
        self.fail_on.add('S2'); rc, out, _ = self.call('measure', '--ecogroup', 'Mbuna', '--run'); self.assertEqual(rc, 1); self.assertIn('1 samples failed', out); self.assertIn('panel_pass.Mzebra_GT3_NCBI.tsv', out); self.assertIn('S2: FAILED', out)

    def test_problems_stop_it_before_anything_is_touched(self):
        self.panel_ok = False; rc, out, _ = self.call('measure', '--ecogroup', 'Mbuna', '--run'); self.assertEqual(rc, 3); self.assertIn('PROBLEM: the panel is not ready', out); self.assertFalse(os.path.exists(self.jd()))
        self.panel_ok = True; rc, out, _ = self.call('measure', '--ecogroup', 'Mbuna', free_bytes=1000); self.assertEqual(rc, 3); self.assertIn('lower --jobs', out)
        rc, out, _ = self.call('measure', '--ecogroup', 'Mbuna', which=lambda t: None); self.assertEqual(rc, 3); self.assertIn('rclone is not on the PATH', out)

    def test_the_selection_must_be_exactly_one(self):
        for args in (['measure'], ['measure', '--ecogroup', 'Mbuna', '-s', 'S1']):
            rc, _, err = self.call(*args); self.assertEqual(rc, 2); self.assertIn('choose exactly one of', err)

    def test_the_offline_options_go_together(self):
        out, err = io.StringIO(), io.StringIO(); rc = C.main(['measure', '--ecogroup', 'Mbuna', '--local-master', self.local], out=out, err=err); self.assertEqual(rc, 2); self.assertIn('go together', err.getvalue())

    def test_no_upload_and_no_ecogroup_column(self):
        rc, out, _ = self.call('measure', '-s', 'S1', '--run', '--no-upload'); self.assertEqual(rc, 0); self.assertFalse(os.path.exists(f'{self.cloud}Bamfiles/{G}/S1/S1.panel.json')); self.assertIn('NOT copied to Dropbox', out)
        pd.DataFrame({'SampleID': self.ids}).to_csv(self.st, sep='\t', index=False); rc, _, err = self.call('measure', '--ecogroup', 'Mbuna'); self.assertEqual(rc, 2); self.assertIn('needs the SampleDatabase with an Ecogroup column', err)


class WriteBase(CliBase):
    HEAD = ['SampleID', 'Species', 'gsdf number', 'TE1', 'TE2']

    def setUp(self):
        super().setUp()
        rc, out, err = self.call('measure', '--ecogroup', 'Mbuna', '--run'); self.assertEqual(rc, 0, out + err)
        self.ws = FakeWorksheet(self.HEAD, [[s, 'Cm', '', '', ''] for s in self.ids] + [['OTHER', 'Lf', '3', '0/0', '']])

    def write(self, *args, ws=None, **kw):
        out, err = io.StringIO(), io.StringIO()
        args = list(args) + ([] if '--config' in args else ['--config', self.base_cfg()])
        rc = C.main(['write'] + list(args) + self.common, out=out, err=err, worksheet_fn=lambda: (ws or self.ws), calls_module=FakeCalls(), **kw)
        return rc, out.getvalue(), err.getvalue()


class TestWrite(WriteBase):
    def test_the_default_is_a_plan_that_writes_nothing_and_names_the_new_columns(self):
        rc, out, err = self.write(); self.assertEqual(rc, 0, err)
        self.assertEqual(self.ws.writes(), []); self.assertEqual(self.ws.grid[0], self.HEAD); self.assertIn('4 samples with a JSON', out); self.assertIn('4 have something to write', out)
        self.assertIn('columns to ADD at the end of the header: gsdf flags, chr7ins, chr7ins flags, TE1 flags, TE2 flags', out); self.assertIn('nothing was written: add --apply', out)

    def test_apply_writes_the_confident_values_and_every_flags_cell_and_nothing_else(self):
        rc, out, err = self.write('--apply'); self.assertEqual(rc, 0, err + out)
        self.assertEqual(self.ws.grid[0], self.HEAD + ['gsdf flags', 'chr7ins', 'chr7ins flags', 'TE1 flags', 'TE2 flags'], 'each value column is followed by its flags column')
        self.assertEqual(self.ws.column('gsdf number'), [3, 2, 3, 2, '3']); self.assertEqual(self.ws.column('chr7ins'), ['0/0'] * 4 + ['']); self.assertEqual(self.ws.column('TE2'), ['0/1'] * 4 + [''])
        self.assertEqual(self.ws.column('TE1'), ['0/0', '', '0/0', '0/0', '0/0'], 'S2 has TE1 no_call: no value; OTHER keeps its own')
        self.assertEqual(self.ws.column('TE1 flags'), ['pilot-v1', 'pilot-v1;no_call', 'pilot-v1', 'pilot-v1', ''], 'the flags cell is written even when the value is not, and starts with the call tag')
        self.assertEqual(self.ws.column('chr7ins flags'), ['pilot-v1', 'pilot-v1', 'pilot-v1', 'pilot-v1;dilution_possible', '']); self.assertEqual(self.ws.column('gsdf flags'), ['pilot-v1'] * 4 + [''])
        self.assertEqual(self.ws.column('Species'), ['Cm'] * 4 + ['Lf'], 'no other column is touched'); self.assertIn('every written cell was read back and is in the right row', out)

    def test_the_new_columns_keep_their_order_even_when_the_first_sample_has_a_blocked_value(self):
        p = f'{self.jd()}S1.panel.json'; d = json.load(open(p))
        for c in d['fake_calls']:
            if c['column'] == 'chr7ins': c.update(value=None, confident=False, flags=['no_call'])          # S1 (alphabetically first) has a flags cell and no value for chr7ins
        json.dump(d, open(p, 'w')); self.write('--apply')
        self.assertEqual(self.ws.grid[0], self.HEAD + ['gsdf flags', 'chr7ins', 'chr7ins flags', 'TE1 flags', 'TE2 flags'])
        self.assertEqual(self.ws.column('chr7ins'), ['', '0/0', '0/0', '0/0', '']); self.assertEqual(self.ws.column('chr7ins flags')[0], 'pilot-v1;no_call')

    def test_column_order_follows_the_config_with_each_flags_column_after_its_value(self):
        flag_of = {('A', 'TE1 flags'): 'TE1', ('A', 'gsdf flags'): 'gsdf number', ('B', 'gsdf flags'): 'gsdf number', ('B', 'TE2 flags'): 'TE2'}
        self.assertEqual(P.column_order(['gsdf number', 'chr7ins', 'TE1', 'TE2'], flag_of), ['gsdf number', 'gsdf flags', 'chr7ins', 'TE1', 'TE1 flags', 'TE2', 'TE2 flags'])

    def test_the_new_columns_keep_their_order_when_a_conflict_removes_the_first_samples_flags_cell(self):
        ws = FakeWorksheet(self.HEAD, [['S1', 'Cm', '', '9', ''], ['S2', 'Cm', '', '', ''], ['S3', 'Cm', '', '', ''], ['S4', 'Cm', '', '', '']])      # S1's TE1 holds a different value: its TE1 flags cell is left alone
        self.write('--apply', ws=ws); self.assertEqual(ws.grid[0], self.HEAD + ['gsdf flags', 'chr7ins', 'chr7ins flags', 'TE1 flags', 'TE2 flags'])
        self.assertEqual(ws.column('TE1'), ['9', '', '0/0', '0/0']); self.assertEqual(ws.column('TE1 flags'), ['', 'pilot-v1;no_call', 'pilot-v1', 'pilot-v1'])

    def test_inversion5_and_the_b_score_are_never_written(self):
        self.write('--apply')
        for col in ('Inversion5', 'Inversion5 flags', 'B score', 'B flags'): self.assertNotIn(col, self.ws.grid[0])
        for c in self.ws.writes(): self.assertNotIn(1.07, c); self.assertFalse('near_threshold' in str(c))

    def test_a_second_apply_writes_nothing(self):
        self.write('--apply'); n = len(self.ws.calls); rc, out, _ = self.write('--apply'); self.assertEqual(rc, 0); self.assertIn('nothing to write', out); self.assertEqual(len(self.ws.calls), n)

    def test_a_different_value_already_there_is_reported_and_never_overwritten_and_its_flags_cell_is_left_alone_too(self):
        ws = FakeWorksheet(self.HEAD + ['gsdf flags'], [['S1', 'Cm', '4', '', '', ''], ['S2', 'Cm', '2', '', '', ''], ['S3', 'Cm', '', '', '', ''], ['S4', 'Cm', '', '', '', '']])
        rc, out, err = self.write('--apply', ws=ws); self.assertEqual(rc, 0, err + out)
        self.assertEqual(ws.column('gsdf number'), ['4', '2', 3, 2], 'S1 holds 4 against the panel\'s 3: left alone; S2 holds 2 = the panel\'s 2: same')
        self.assertEqual(ws.column('gsdf flags'), ['', 'pilot-v1', 'pilot-v1', 'pilot-v1'], 'S1: the conflicting value\'s flags cell is not written; S2 (same value): it is')
        self.assertIn('conflict 1', out); self.assertIn('a different value is already there (left alone)', out); self.assertIn('S1', out); self.assertIn("sheet '4', measured 3", out); self.assertIn('1 flags cells are left alone', out)

    def test_there_is_no_overwrite_option(self):
        with self.assertRaises(SystemExit): C.build_parser().parse_args(['write', '--overwrite'])

    def test_limit_writes_the_first_samples_only(self):
        rc, out, err = self.write('--apply', '--limit', '2'); self.assertEqual(rc, 0, err)
        self.assertEqual(self.ws.column('gsdf number')[:4], [3, 2, '', '']); self.assertIn('--limit 2: this run writes only 2 samples (S1, S2)', out)

    def test_a_refused_measurement_and_a_missing_json_are_reported(self):
        self.refuse.add('S3'); self.fail_on.discard('S3'); self.write('--apply'); jd = self.jd()
        json.dump({'schema': 2, 'bam': 'S3.all.bam', 'bam_size_bytes': 1500, 'fake_calls': fake_calls(3), 'refuse': True}, open(jd + 'S3.panel.json', 'w'))
        os.remove(jd + 'S4.panel.json')
        rc, out, _ = self.write('-s', 'S3', 'S4', 'S1'); self.assertEqual(rc, 0)
        self.assertIn('REFUSED S3: DefinitionMismatch: measured under another definition', out); self.assertIn('1 selected samples have no JSON (S4)', out)
        self.assertEqual(self.ws.column('gsdf number')[2], 3, 'the refused sample was written by the first apply and is not touched now')

    def test_the_journal_lists_the_changes(self):
        self.write('--apply'); jdir = self.jd() + 'sheet_writes/'; files = os.listdir(jdir); self.assertEqual(len(files), 1); self.assertTrue(files[0].startswith('SampleDatabase.') and files[0].endswith('.tsv'))
        text = open(jdir + files[0]).read(); self.assertIn('gsdf number', text); self.assertIn('S1', text)

    def test_a_missing_json_folder_and_nothing_to_write_are_errors(self):
        shutil.rmtree(self.jd()); rc, _, err = self.write(); self.assertEqual(rc, 2); self.assertIn('has no *.panel.json files', err)

    def test_the_tab_goes_to_the_live_sheet_and_offline_options_are_refused_there(self):
        seen = []
        class FM:
            s_ID = 'SPREADSHEET'
            class gc:
                @staticmethod
                def open_by_key(key):
                    seen.append(('key', key)); return types.SimpleNamespace(worksheet=lambda tab: (seen.append(('tab', tab)), self.ws)[1])
        out, err = io.StringIO(), io.StringIO()
        rc = C.main(['write', '--json-dir', self.jd(), '--tab', 'SampleDatabase_panel_test', '--config', self.base_cfg()], out=out, err=err, calls_module=FakeCalls(), fm_factory=FM)
        self.assertEqual(rc, 0, err.getvalue()); self.assertIn(('key', 'SPREADSHEET'), seen); self.assertIn(('tab', 'SampleDatabase_panel_test'), seen); self.assertIn("tab 'SampleDatabase_panel_test'", out.getvalue())
        out, err = io.StringIO(), io.StringIO(); rc = C.main(['write', '--json-dir', self.jd(), '--sample-table', self.st], out=out, err=err, calls_module=FakeCalls()); self.assertEqual(rc, 2); self.assertIn('writing needs the live sheet', err.getvalue())

    def test_the_real_call_layer_missing_is_a_clear_error(self):
        try:
            pkg = importlib.import_module('sd_genotypes.panel')
        except ImportError:
            pkg = None
        with contextlib.ExitStack() as st:
            st.enter_context(mock.patch.dict(sys.modules, {'sd_genotypes.panel.calls': None}))      # an import of it fails ...
            if pkg is not None:                                                                       # ... even if an earlier test left `calls` as an attribute of the package (`from package import calls` would find it)
                st.enter_context(mock.patch.object(pkg, 'calls', types.SimpleNamespace(), create=True))
            with self.assertRaises(P.PanelPassError) as cm: P.load_calls_module()
            self.assertIn('sd_genotypes.panel.calls cannot be imported', str(cm.exception)); self.assertIn('call layer (decision 76)', str(cm.exception))
            out, err = io.StringIO(), io.StringIO(); rc = C.main(['write', '--json-dir', self.jd()], out=out, err=err, worksheet_fn=lambda: self.ws)
            self.assertEqual(rc, 2); self.assertIn('sd_genotypes.panel.calls cannot be imported', err.getvalue()); self.assertEqual(self.ws.writes(), [])


class NoCallBase(CliBase):
    """Hand-made JSONs whose fake calls cover the cases: confident; no_call; low_gq (a value that is not confident); site_error; a number with no estimate; an information-only flag."""
    HEAD = ['SampleID', 'gsdf number', 'TE1', 'TE2']

    def call_rec(self, column, flags_column, value, flags, confident):
        return dict(column=column, flags_column=flags_column, value=value, flags=flags, confident=confident, details={})

    def setUp(self):
        super().setUp()
        jd = self.jd(); os.makedirs(jd, exist_ok=True)
        rec = self.call_rec
        self.cases = {
            'N1': [rec('gsdf number', 'gsdf flags', 2, [], True), rec('chr7ins', 'chr7ins flags', '0/0', [], True), rec('TE1', 'TE1 flags', '0/0', [], True), rec('TE2', 'TE2 flags', '0/1', [], True)],
            'N2': [rec('gsdf number', 'gsdf flags', None, ['low_gq'], False), rec('chr7ins', 'chr7ins flags', None, ['no_call', 'site_error'], False), rec('TE1', 'TE1 flags', None, ['no_call'], False), rec('TE2', 'TE2 flags', '0/1', ['low_gq'], False)],
            'N3': [rec('gsdf number', 'gsdf flags', 3, [], True), rec('chr7ins', 'chr7ins flags', '1/1', ['dilution_possible'], True), rec('TE1', 'TE1 flags', '0/0', [], True), rec('TE2', 'TE2 flags', '0/0', [], True)],
            'N4': [rec('gsdf number', 'gsdf flags', 2, [], True), rec('chr7ins', 'chr7ins flags', '0/0', [], True), rec('TE1', 'TE1 flags', '0/0', ['low_gq'], True), rec('TE2', 'TE2 flags', None, ['junction_uninformative'], False)],      # N4: a CONFIDENT value that carries low_gq (a config whose blocking flags differ); a withheld call with an unrelated flag
        }
        for s, calls in self.cases.items():
            json.dump(dict(bam=f'{s}.bam', bam_size_bytes=10, fake_calls=calls), open(f'{jd}{s}.panel.json', 'w'))
        self.ws = FakeWorksheet(self.HEAD, [[s, '', '', ''] for s in self.cases] + [['OTHER', '5', '0/0', '']])

    def cfg(self, **sheet):
        return self.base_cfg(**sheet)

    def write(self, *args, ws=None):
        out, err = io.StringIO(), io.StringIO()
        args = list(args) + ([] if '--config' in args else ['--config', self.base_cfg()])
        rc = C.main(['write', '--json-dir', self.jd()] + list(args) + self.common, out=out, err=err, worksheet_fn=lambda: (ws or self.ws), calls_module=FakeCalls())
        return rc, out.getvalue(), err.getvalue()

    NOCALL = dict(value='./.', columns=['chr7ins', 'TE1', 'TE2'], flags=['no_call', 'low_gq'], exclude_flags=['site_error'])


class TestNoCall(NoCallBase):
    def test_off_unless_configured_a_withheld_genotype_stays_empty(self):
        rc, out, err = self.write('--apply'); self.assertEqual(rc, 0, err + out)
        self.assertEqual(self.ws.column('TE1')[:4], ['0/0', '', '0/0', '0/0']); self.assertEqual(self.ws.column('TE2')[:4], ['0/1', '', '0/0', ''])
        self.assertNotIn("'./.'", out)

    def test_no_call_and_low_gq_become_dot_slash_dot_in_genotype_columns_only(self):
        rc, out, err = self.write('--apply', '--config', self.cfg(nocall=self.NOCALL)); self.assertEqual(rc, 0, err + out)
        self.assertEqual(self.ws.column('TE1')[:4], ['0/0', './.', '0/0', '0/0'], 'no_call -> ./.; N4 has a CONFIDENT value that carries low_gq: it stays')
        self.assertEqual(self.ws.column('TE2')[:4], ['0/1', './.', '0/0', ''], 'low_gq: the genotype exists but is below the GQ cut; an unrelated flag (N4) gives no ./.')
        self.assertEqual(self.ws.column('chr7ins')[:4], ['0/0', '', '1/1', '0/0'], 'no_call together with site_error stays empty: an error of the site, not a lack of data')
        self.assertEqual(self.ws.column('gsdf number')[:4], [2, '', 3, 2], 'a number is never ./., even with a nocall-type flag')
        self.assertIn("2 genotype cells of the calls get './.'", out)
        self.assertEqual(self.ws.grid[5][1:4], ['5', '0/0', ''], 'a row that is not a sample of the pass is untouched'); self.assertEqual(set(self.ws.grid[5][4:]), {''})

    def test_an_information_only_flag_does_not_turn_a_value_into_no_call(self):
        rc, out, err = self.write('--apply', '--config', self.cfg(nocall=self.NOCALL)); self.assertEqual(rc, 0)
        self.assertEqual(self.ws.column('chr7ins')[2], '1/1')

    def test_a_value_already_in_the_sheet_is_never_overwritten_by_dot_slash_dot(self):
        self.ws = FakeWorksheet(self.HEAD, [['N1', '', '', ''], ['N2', '', '0/1', ''], ['N3', '', '', ''], ['N4', '', '', '']])
        rc, out, err = self.write('--apply', '--config', self.cfg(nocall=self.NOCALL)); self.assertEqual(rc, 0, err + out)
        self.assertEqual(self.ws.column('TE1')[1], '0/1'); self.assertIn('conflict', out.lower())

    def test_a_second_run_finds_nothing_to_write(self):
        c = self.cfg(nocall=self.NOCALL); self.write('--apply', '--config', c)
        rc, out, err = self.write('--apply', '--config', c); self.assertEqual(rc, 0); self.assertIn('nothing to write', out)

    def test_the_config_is_checked(self):
        for bad, why in ((dict(nocall=dict(value='./.', columns=['TE9'], flags=['no_call'])), 'not in sheet.write_columns'), (dict(nocall=dict(value='./.', columns=['TE1'])), 'sheet.nocall.flags is missing'), (dict(flags_in_sheet='no'), 'true or false')):
            p = self.cfg(**bad)
            with self.assertRaises(ValueError) as cm:
                P.load_config(p)
            self.assertIn(why, str(cm.exception))

    def test_record_updates_directly(self):
        r = dict(calls=[dict(c, value_to_write=c['value'] if c['confident'] else None, call_tag='pilot-v1') for c in self.cases['N2']])
        upd, pairs = P.record_updates(r, FakeCalls(), ['gsdf number', 'chr7ins', 'TE1', 'TE2'], dict(flags_in_sheet=False, nocall=self.NOCALL))
        self.assertEqual(upd, {'TE1': './.', 'TE2': './.'}); self.assertEqual(pairs, {})


class TestFlagsInSheet(NoCallBase):
    def test_the_default_still_writes_the_flags_columns(self):
        rc, out, err = self.write('--apply'); self.assertEqual(rc, 0, err + out)
        self.assertIn('TE1 flags', self.ws.grid[0]); self.assertEqual(self.ws.column('TE1 flags')[1], 'pilot-v1;no_call')

    def test_false_writes_the_values_and_no_flags_cells_at_all(self):
        rc, out, err = self.write('--apply', '--config', self.cfg(flags_in_sheet=False)); self.assertEqual(rc, 0, err + out)
        self.assertEqual(self.ws.grid[0], self.HEAD + ['chr7ins'], 'only the value column chr7ins is new')
        self.assertFalse([h for h in self.ws.grid[0] if 'flags' in h])
        self.assertEqual(self.ws.column('TE1')[:4], ['0/0', '', '0/0', '0/0']); self.assertEqual(self.ws.column('chr7ins')[:4], ['0/0', '', '1/1', '0/0'])
        self.assertIn('the flags are not written to the sheet', out)

    def test_the_journal_records_the_call_tag_in_its_name(self):
        self.write('--apply', '--config', self.cfg(flags_in_sheet=False))
        files = os.listdir(self.jd() + 'sheet_writes/'); self.assertEqual(len(files), 1)
        self.assertTrue(files[0].startswith('SampleDatabase.') and files[0].endswith('.pilot-v1.tsv'), files[0])

    def test_the_table_keeps_every_flag_whatever_the_sheet_does(self):
        out, err = io.StringIO(), io.StringIO(); path = self.tmp + '/calls.tsv'
        rc = C.main(['table', '--json-dir', self.jd(), '--out', path] + self.common, out=out, err=err, calls_module=FakeCalls()); self.assertEqual(rc, 0, err.getvalue())
        t = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False).set_index('SampleID')
        self.assertEqual(t.loc['N2', 'TE1 flags'], 'pilot-v1;no_call'); self.assertEqual(t.loc['N2', 'chr7ins flags'], 'pilot-v1;no_call;site_error'); self.assertEqual(t.loc['N2', 'TE1'], '')


class ScopeBase(NoCallBase):
    """Six calls per sample, an Ecogroup column in the sheet, and the cases the scope and the extended nocall rule have to tell apart."""
    HEAD = ['SampleID', 'Ecogroup', 'gsdf number', 'TE1', 'TE2']

    def setUp(self):
        CliBase.setUp(self)
        jd = self.jd(); os.makedirs(jd, exist_ok=True); rec = self.call_rec
        def six(inv, binv, b, bf=(), te1=('0/0', (), True)):
            return [rec('gsdf number', 'gsdf flags', 2, [], True), rec('chr7ins', 'chr7ins flags', '0/0', [], True), rec('TE1', 'TE1 flags', te1[0], list(te1[1]), te1[2]), rec('TE2', 'TE2 flags', '0/0', [], True),
                    rec('Inversion5', 'Inversion5 flags', inv[0], list(inv[1]), inv[2]), rec('B score', 'B flags', b[0], list(b[1]), b[2])]
        self.cases = {
            'S1': ('Mbuna', six((0, (), True), None, (1.1, (), True))),
            'S2': ('Utaka', six((1, (), True), None, (5.0, (), True))),
            'S3': ('Mbuna', six((None, ('low_depth',), False), None, (None, ('low_depth',), False))),
            'S4': ('Mbuna', six((2, ('unvalidated_value',), False), None, (1.0, (), True))),
            'S5': ('', six((0, (), True), None, (1.2, (), True))),
            'S6': ('Mbuna', six((None, ('few_markers', 'low_depth'), False), None, (1.0, (), True))),
        }
        for s, (eco, calls) in self.cases.items():
            json.dump(dict(bam=f'{s}.bam', bam_size_bytes=10, fake_calls=calls), open(f'{jd}{s}.panel.json', 'w'))
        self.ws = FakeWorksheet(self.HEAD, [[s, eco, '', '', ''] for s, (eco, _) in self.cases.items()])

    def cfg6(self, **sheet):
        import yaml
        cfg = P.load_config(); cfg['sheet'].update(sheet)                          # the SHIPPED config: six columns, nocall, scope
        p = self.tmp + f'/cfg{len(os.listdir(self.tmp))}.yaml'; yaml.safe_dump(cfg, open(p, 'w')); return p

    def got(self):
        return {str(r[0]): {k: str(v) for k, v in zip(self.ws.grid[0], r)} for r in self.ws.grid[1:]}


class TestScope(ScopeBase):
    def test_inversion5_and_b_score_are_written_for_the_listed_ecogroups_only(self):
        rc, out, err = self.write('--apply', '--config', self.cfg6()); self.assertEqual(rc, 0, err + out); g = self.got()
        self.assertEqual([g['S1'][c] for c in ('Inversion5', 'Inversion5 flags', 'B score', 'B flags')], ['0', 'pilot-v1', '1.1', 'pilot-v1'])
        for s in ('S2', 'S5'):                                                     # Utaka, and a sample whose Ecogroup is empty: outside every scope
            self.assertEqual([g[s][c] for c in ('Inversion5', 'Inversion5 flags', 'B score', 'B flags')], ['', '', '', ''], s)
            self.assertEqual((g[s]['gsdf number'], g[s]['gsdf flags'], g[s]['chr7ins'], g[s]['TE1']), ('2', 'pilot-v1', '0/0', '0/0'), 'the other columns are written for every ecogroup')
        self.assertIn('Inversion5: left out for 2 samples whose Ecogroup is not Mbuna', out); self.assertIn('B score: left out for 2 samples', out)

    def test_low_depth_inversion5_is_dot_slash_dot_but_a_value_held_back_on_purpose_and_a_score_stay_empty(self):
        rc, out, err = self.write('--apply', '--config', self.cfg6()); self.assertEqual(rc, 0, err + out); g = self.got()
        self.assertEqual((g['S3']['Inversion5'], g['S6']['Inversion5']), ('./.', './.'), 'low_depth / few_markers: too little data to genotype')
        self.assertEqual((g['S3']['B score'], g['S3']['B flags']), ('', 'pilot-v1;low_depth'), 'a score is never ./.')
        self.assertEqual((g['S4']['Inversion5'], g['S4']['Inversion5 flags']), ('', 'pilot-v1;unvalidated_value'), 'an Inversion5 of 2 is held back on purpose: not a lack of data, so not ./.')
        self.assertEqual(g['S4']['B score'], '1.0')

    def test_no_scope_means_every_ecogroup(self):
        rc, out, err = self.write('--apply', '--config', self.cfg6(scope={})); self.assertEqual(rc, 0, err + out); g = self.got()
        self.assertEqual((g['S2']['Inversion5'], g['S2']['B score'], g['S5']['Inversion5']), ('1', '5.0', '0')); self.assertNotIn('left out', out)

    def test_a_tab_without_an_ecogroup_column_is_refused_and_nothing_is_written(self):
        ws = FakeWorksheet(['SampleID', 'gsdf number', 'TE1', 'TE2'], [[s, '', '', ''] for s in self.cases])
        rc, out, err = self.write('--apply', '--config', self.cfg6(), ws=ws)
        self.assertEqual(rc, 3); self.assertIn('cannot give it', err); self.assertIn("'Ecogroup'", err); self.assertEqual(list(ws.calls), [])

    def test_the_config_is_checked(self):
        for bad, why in ((dict(scope={'TE9': ['Mbuna']}), 'not in sheet.write_columns'), (dict(scope={'B score': []}), 'non-empty list'), (dict(scope={'B score': 'Mbuna'}), 'non-empty list')):
            with self.assertRaises(ValueError) as cm:
                P.load_config(self.cfg6(**bad))
            self.assertIn(why, str(cm.exception))

    def test_collect_updates_needs_the_ecogroups_when_a_scope_is_set(self):
        cfg = P.load_config()
        with self.assertRaises(ValueError):
            P.collect_updates(self.jd(), ['S1'], FakeCalls(), cfg['sheet']['write_columns'], cfg['sheet'])

    def test_record_updates_directly(self):
        eco, calls = self.cases['S3']; r = dict(calls=[dict(c, value_to_write=c['value'] if c['confident'] else None, call_tag='pilot-v1') for c in calls]); sc = P.load_config()['sheet']
        upd, pairs = P.record_updates(r, FakeCalls(), sc['write_columns'], sc, 'Mbuna'); self.assertEqual((upd['Inversion5'], 'B score' in upd, upd['B flags']), ('./.', False, 'pilot-v1;low_depth'))
        upd, pairs = P.record_updates(r, FakeCalls(), sc['write_columns'], sc, 'Utaka'); self.assertFalse({'Inversion5', 'Inversion5 flags', 'B score', 'B flags'} & set(upd)); self.assertNotIn('Inversion5 flags', pairs)


class TestTable(WriteBase):
    def test_every_call_and_its_numbers_in_one_file(self):
        out = self.tmp + '/calls.tsv'; rc, text, err = self.call('table', '--out', out); self.assertEqual(rc, 0, err)
        t = pd.read_csv(out, sep='\t', dtype=str, keep_default_na=False); self.assertEqual(list(t.SampleID), self.ids); self.assertEqual(list(t.Species), ['Cm', 'Cm', 'Lf', 'Lf'])
        for col in ('gsdf number', 'chr7ins', 'TE1', 'TE2', 'Inversion5', 'B score', 'Inversion5 flags', 'B flags', 'gsdf number.ratio', 'Inversion5.median_site_fraction', 'B score.blocks', 'TE1 confident'): self.assertIn(col, t.columns)
        self.assertEqual(list(t['TE1 flags']), ['pilot-v1', 'pilot-v1;no_call', 'pilot-v1', 'pilot-v1']); self.assertEqual(list(t['TE1 confident']), ['True', 'False', 'True', 'True']); self.assertEqual(list(t['TE1']), ['0/0', '', '0/0', '0/0'])
        self.assertIn('4 samples', text); self.assertIn('TE1: 3 confident', text); self.assertIn('no_call 1', text); self.assertIn('wrote ' + out, text)

    def test_an_unreadable_sheet_does_not_stop_the_table_it_only_loses_the_species_columns(self):
        from variant_calling import cohort_cli
        out, err = io.StringIO(), io.StringIO(); path = self.tmp + '/t.tsv'
        with mock.patch.object(cohort_cli, 'load_tables', side_effect=ImportError('No module named gspread')):            # never the real sheets in a test
            rc = C.main(['table', '--json-dir', self.jd(), '--out', path], out=out, err=err, calls_module=FakeCalls())
        self.assertEqual(rc, 0, err.getvalue()); t = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
        self.assertNotIn('Species', t.columns); self.assertEqual(len(t), 4); self.assertIn('note: the SampleDatabase could not be read, so the table has no Species or Ecogroup column (ImportError', out.getvalue())

    def test_the_default_place_and_a_selection(self):
        rc, text, _ = self.call('table', '-s', 'S1', 'S3'); self.assertEqual(rc, 0); t = pd.read_csv(self.jd() + 'panel_calls.tsv', sep='\t'); self.assertEqual(list(t.SampleID), ['S1', 'S3'])


if __name__ == '__main__':
    unittest.main()
