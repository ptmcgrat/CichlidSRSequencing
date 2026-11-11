"""`effcov` (module 05): the effective coverage of BAMs. Pure arithmetic and parsing, the definition, the output file, one sample and a whole run with FAKE samtools / bcftools / rclone
(fake_tools.py describes a BAM by a few numbers, so every expected value here is exact), resume and failure, the preflight, and the command line.
The fakes establish that the right command lines are issued and the numbers are combined correctly. How real samtools and bcftools behave is what the first real run is for."""
import contextlib, io, json, os, shutil, stat, sys, tempfile, unittest
from unittest import mock
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pandas as pd
from cohort_fixtures import *
from variant_calling import build as B, cli, effcov as E, effcov_cli
from cichlid_common.fake_sheet import FakeWorksheet

CFG = B.load_build_config()
ECFG = dict(E.load_config(), window_bp=1000, jobs=2, threads=2)
FAI = [('NC_1', 20_000_000), ('NC_2', 15_000_000), ('NC_3', 8_000_000), ('NC_MT', 16_000), ('NW_1', 900_000)]
# what each fake BAM looks like: 1000 primary reads, 291 duplicates, 700 mapped and not duplicate, 500 of them with MAPQ >= 13; the pileup emits 90% of the window positions with
# INFO/DP cycling 5, 6 (mean 5.5) and FORMAT/DP cycling 3, 4, 4, 4, 4 (mean 3.8)
GOOD = dict(primary=1000, dups=291, mapped_nondup=700, mapq_pass=500, raw=[5, 6], hq=[3, 4, 4, 4, 4], cover=0.9)


def fake_install(tmp):
    bindir = os.path.join(tmp, 'bin'); os.makedirs(bindir)
    for tool in ('bcftools', 'samtools', 'rclone'):
        p = os.path.join(bindir, tool)
        open(p, 'w').write(f"#!{sys.executable}\nimport sys\nsys.path.insert(0, {HERE!r})\nimport fake_tools\nsys.exit(fake_tools.main({tool!r}, sys.argv[1:]))\n")
        os.chmod(p, os.stat(p).st_mode | stat.S_IEXEC)
    return bindir


class EffBase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.saved = {k: os.environ.get(k) for k in ('PATH', 'FAKE_LOG')}
        os.environ['PATH'] = fake_install(self.tmp) + os.pathsep + os.environ['PATH']
        self.logfile = os.path.join(self.tmp, 'calls.log'); os.environ['FAKE_LOG'] = self.logfile
        self.local, self.cloud = os.path.join(self.tmp, 'local') + '/', os.path.join(self.tmp, 'cloud') + '/'
        os.makedirs(self.local)
        self.fa = os.path.join(self.tmp, 'ref.fa'); open(self.fa, 'w').write('>NC_1\nACGT\n')
        open(self.fa + '.fai', 'w').write(''.join(f'{c}\t{n}\t0\t60\t61\n' for c, n in FAI))
        self.env = B.Env(self.local, self.cloud, G, self.fa)
        self.out = os.path.join(self.tmp, 'out', 'eff.tsv')
        self.ids = ['S1', 'S2', 'S3', 'S4']
        self.sizes = {}
        for s in self.ids:
            self.put(s, GOOD)

    def tearDown(self):
        for k, v in self.saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        shutil.rmtree(self.tmp, ignore_errors=True)

    def put(self, s, desc, raw_text=None):
        d = f'{self.cloud}Bamfiles/{G}/{s}/'; os.makedirs(d, exist_ok=True)
        open(d + f'{s}.all.bam', 'w').write(raw_text if raw_text is not None else json.dumps(desc)); open(d + f'{s}.all.bai', 'w').write('BAI')
        self.sizes[s] = os.path.getsize(d + f'{s}.all.bam')

    def ledger(self, ids=None, cov=10.0):
        return pd.DataFrame([dict(SampleID=s, GenomeVersion=G, Coverage=cov, BamSize=self.sizes[s]) for s in (ids or self.ids)])

    def plan(self, members=None, ecfg=None, redo=False, ledger=None):
        return E.make_plan(self.env, CFG, ecfg or ECFG, self.ids if members is None else members, self.ledger() if ledger is None else ledger, self.out, redo=redo)

    def calls(self, tool):
        return [json.loads(l)['argv'] for l in open(self.logfile) if json.loads(l)['tool'] == tool] if os.path.exists(self.logfile) else []


class TestArithmetic(unittest.TestCase):
    def test_windows_one_per_chromosome_inside_it(self):
        w = E.windows_from_fai(FAI, 1000, 0.4, 5_000_000)
        self.assertEqual(w, [('NC_1', 8_000_001, 8_001_000), ('NC_2', 6_000_001, 6_001_000), ('NC_3', 3_200_001, 3_201_000)], 'NC_MT and NW_1 are not chromosomes of at least 5 Mb')
        self.assertTrue(all(e - s + 1 == 1000 for _, s, e in w))
        self.assertEqual(E.regions_arg(w), 'NC_1:8000001-8001000,NC_2:6000001-6001000,NC_3:3200001-3201000')

    def test_a_window_is_clipped_to_the_end_of_the_chromosome(self):
        self.assertEqual(E.windows_from_fai([('NC_9', 5_000_000)], 1_000_000, 0.99, 5_000_000), [('NC_9', 4_000_001, 5_000_000)])

    def test_flagstat_takes_the_primary_lines_only(self):
        text = ('63075860 + 0 in total (QC-passed reads + QC-failed reads)\n63075860 + 0 primary\n0 + 0 secondary\n641137 + 0 supplementary\n18341603 + 0 duplicates\n'
                '18341603 + 0 primary duplicates\n62336968 + 0 mapped (97.83% : N/A)\n61695831 + 0 primary mapped (97.81% : N/A)\n')
        self.assertEqual(E.parse_flagstat(text), dict(primary=63075860, primary_duplicates=18341603, primary_mapped=61695831), 'primary mapped is not confused with primary')
        no_mapped = '1000 + 0 primary\n291 + 0 primary duplicates\n'
        self.assertEqual(E.parse_flagstat(no_mapped), dict(primary=1000, primary_duplicates=291), 'older output: the count is simply absent')
        with self.assertRaises(E.MeasureError):
            E.parse_flagstat('1000 + 0 in total\n291 + 0 duplicates\n')

    def test_the_window_line(self):
        self.assertEqual(E.parse_window_line('2700\t14850\t10260\n'), (2700, 14850, 10260))
        for bad in ('', 'x\ty\tz', '1\t2'):
            with self.assertRaises(E.MeasureError):
                E.parse_window_line(bad)

    def test_the_columns_from_the_counts(self):
        d = E.derive(1000, 291, 700, 500, 2700, 14850, 10260, 3000)
        self.assertAlmostEqual(d['DuplicateFraction'], 0.291); self.assertAlmostEqual(d['MapqPassFraction'], 500 / 700)
        self.assertAlmostEqual(d['WindowRawDepth'], 5.5); self.assertAlmostEqual(d['EffectiveCoverage'], 3.8)
        self.assertAlmostEqual(d['BaseSurvivalFraction'], 3.8 / 5.5); self.assertAlmostEqual(d['WindowCoveredFraction'], 0.9); self.assertEqual(d['PrimaryReads'], 1000)

    def test_the_duplication_rate_is_over_mapped_reads_like_the_pipelines_own(self):
        d = E.derive(63075860, 18341603, 43354228, 31489962, 2700, 14850, 10260, 3000, primary_mapped=61695831)
        self.assertAlmostEqual(d['DuplicationRate'], 18341603 / 61695831); self.assertAlmostEqual(d['DuplicateFraction'], 18341603 / 63075860)
        self.assertGreater(d['DuplicationRate'], d['DuplicateFraction'], 'the all-primary fraction is lower by the unmapped share'); self.assertEqual(d['PrimaryMappedReads'], 61695831)
        self.assertNotIn('DuplicationRate', E.derive(1000, 291, 700, 500, 2700, 14850, 10260, 3000))

    def test_nothing_to_divide_by_is_an_error_not_a_nan(self):
        for args in ((0, 0, 0, 0, 10, 10, 10, 100), (1000, 5, 0, 0, 10, 10, 10, 100), (1000, 5, 700, 500, 0, 0, 0, 100)):
            with self.assertRaises(E.MeasureError):
                E.derive(*args)

    def test_the_definition_key_changes_with_every_parameter_that_changes_the_numbers(self):
        w = E.windows_from_fai(FAI, 1000, 0.4, 5_000_000)
        base = E.definition_key(CFG, ECFG, w)
        self.assertEqual(base, E.definition_key(CFG, ECFG, w)); self.assertEqual(len(base), 12)
        other_cfg = json.loads(json.dumps(CFG)); other_cfg['bcftools']['min_mapping_quality'] = 20
        other_bq = json.loads(json.dumps(CFG)); other_bq['bcftools']['min_base_quality'] = 10
        self.assertNotEqual(base, E.definition_key(other_cfg, ECFG, w)); self.assertNotEqual(base, E.definition_key(other_bq, ECFG, w))
        self.assertNotEqual(base, E.definition_key(CFG, dict(ECFG, definition_version=2), w))
        self.assertNotEqual(base, E.definition_key(CFG, ECFG, E.windows_from_fai(FAI, 1000, 0.5, 5_000_000)))
        self.assertEqual(base, E.definition_key(CFG, dict(ECFG, jobs=9, threads=1, min_depth=3), w), 'parallelism and the report threshold do not change a measurement')


class TestConfig(unittest.TestCase):
    def test_the_shipped_config_loads_with_numbers_only(self):
        c = E.load_config()
        for k in ('window_bp', 'jobs', 'threads', 'default_bam_bytes', 'min_depth', 'position_fraction', 'min_free_factor'):
            self.assertIsInstance(c[k], (int, float), k)
        self.assertEqual(c['default_bam_bytes'], 15_000_000_000)

    def test_a_number_yaml_reads_as_text_is_still_a_number(self):
        p = os.path.join(tempfile.mkdtemp(), 'c.yaml'); open(p, 'w').write('default_bam_bytes: 15.0e9\njobs: "3"\n')
        c = E.load_config(p); self.assertEqual((c['default_bam_bytes'], c['jobs']), (15e9, 3))
        open(p, 'w').write('window_bp: many\n')
        with self.assertRaises(ValueError):
            E.load_config(p)


class TestOutputFile(unittest.TestCase):
    def setUp(self):
        self.path = os.path.join(tempfile.mkdtemp(), 'sub', 'eff.tsv')

    def test_header_once_rows_appended_and_read_back(self):
        E.append_row(self.path, dict(SampleID='A', Status='ok', Error='', EffectiveCoverage=3.8, DefinitionKey='k1'))
        E.append_row(self.path, dict(SampleID='B', Status='failed', Error='bad\tthing\nhappened', DefinitionKey='k1'))
        lines = open(self.path).read().splitlines()
        self.assertEqual(len(lines), 3); self.assertEqual(lines[0].split('\t'), E.COLUMNS)
        rows = E.read_rows(self.path)
        self.assertEqual((rows[0]['SampleID'], rows[0]['EffectiveCoverage']), ('A', '3.8')); self.assertEqual(rows[1]['Error'], 'bad thing happened', 'a tab or a newline cannot break a line')
        self.assertEqual(E.read_rows(os.path.join(os.path.dirname(self.path), 'none.tsv')), [])

    def test_latest_ok_is_the_newest_success_of_the_same_definition(self):
        rows = [dict(SampleID='A', Status='ok', DefinitionKey='k1', EffectiveCoverage='1'), dict(SampleID='A', Status='ok', DefinitionKey='k1', EffectiveCoverage='2'),
                dict(SampleID='B', Status='failed', DefinitionKey='k1'), dict(SampleID='C', Status='ok', DefinitionKey='old')]
        got = E.latest_ok(rows, 'k1')
        self.assertEqual(sorted(got), ['A']); self.assertEqual(got['A']['EffectiveCoverage'], '2')


class TestMigration(unittest.TestCase):
    def test_a_file_from_an_older_version_is_migrated_before_a_longer_row_is_appended(self):
        d = tempfile.mkdtemp(); path = os.path.join(d, 'eff.tsv')
        old = E.COLUMNS[:-2]
        open(path, 'w').write('\t'.join(old) + '\n' + '\t'.join(['A', G, 'ok', '', 'k', '2026-11-11', 'host', 'tools', '1000', '0.29', '0.72', '5.5', '3.8', '0.69', '0.9', '9.9', '123', '40.0']) + '\n')
        E.append_row(path, dict(SampleID='B', Status='ok', DefinitionKey='k', EffectiveCoverage=4.0, PrimaryMappedReads=995, DuplicationRate=0.29246))
        lines = open(path).read().splitlines()
        self.assertEqual(lines[0].split('\t'), E.COLUMNS); self.assertEqual(len(lines), 3)
        rows = E.read_rows(path); a, b = rows[0], rows[1]
        self.assertEqual((a['SampleID'], a['EffectiveCoverage'], a['TableCoverage'], a['Seconds'], a['DuplicationRate']), ('A', '3.8', '9.9', '40.0', ''))
        self.assertEqual((b['SampleID'], b['EffectiveCoverage'], b['PrimaryMappedReads'], b['DuplicationRate']), ('B', '4', '995', '0.29246'))
        bak = [f for f in os.listdir(d) if '.bak-' in f]
        self.assertEqual(len(bak), 1, 'the old file is kept'); self.assertEqual(open(os.path.join(d, bak[0])).read().count('\n'), 2)

    def test_a_current_file_is_not_touched_and_an_unknown_column_is_refused(self):
        d = tempfile.mkdtemp(); path = os.path.join(d, 'eff.tsv')
        E.append_row(path, dict(SampleID='A', Status='ok')); E.append_row(path, dict(SampleID='B', Status='ok'))
        self.assertEqual(os.listdir(d), ['eff.tsv'])
        open(path, 'w').write('SampleID\tFromTheFuture\nA\t1\n')
        with self.assertRaises(RuntimeError):
            E.append_row(path, dict(SampleID='B', Status='ok'))

    def test_latest_ok_can_require_the_newest_columns(self):
        rows = [dict(SampleID='A', Status='ok', DefinitionKey='k', DuplicationRate=''), dict(SampleID='B', Status='ok', DefinitionKey='k', DuplicationRate='0.3')]
        self.assertEqual(sorted(E.latest_ok(rows, 'k')), ['A', 'B']); self.assertEqual(sorted(E.latest_ok(rows, 'k', E.REQUIRED)), ['B'])


class TestReport(unittest.TestCase):
    def rows(self):
        def r(s, eff, table, dup=0.1, mq=0.8, bs=0.9, cov=0.95):
            return dict(SampleID=s, Status='ok', DefinitionKey='k', EffectiveCoverage=eff, TableCoverage=table, DuplicateFraction=dup, MapqPassFraction=mq, BaseSurvivalFraction=bs, WindowCoveredFraction=cov)
        return [r('A', '3.8', '9.9', 0.291, 0.726, 0.68), r('B', '16.9', '19.9'), r('C', '11.9', '')]

    def test_the_table_and_what_stands_out(self):
        t = E.report_text(self.rows(), key='k', sex={'A': 'F', 'B': 'M', 'C': 'M'}, min_depth=5)
        self.assertIn('below 5x', t); self.assertIn('3 samples; effective coverage 3.80 to 16.90x (median 11.90x)', t); self.assertIn('1 below 5x', t)
        self.assertIn('effective / table: 0.38 to 0.85', t); self.assertIn('sex M: 2 samples, median effective 14.40x, 0 below 5x', t); self.assertIn('sex F: 1 samples', t)
        a = [l for l in t.splitlines() if l.startswith('A ')][0]
        self.assertIn('0.38', a); self.assertIn('29.1', a); self.assertIn('72.6', a); self.assertIn('68.0', a)

    def test_other_definitions_and_other_samples_are_left_out(self):
        rows = self.rows() + [dict(SampleID='D', Status='ok', DefinitionKey='other', EffectiveCoverage='1', TableCoverage='1', DuplicateFraction=0, MapqPassFraction=1, BaseSurvivalFraction=1, WindowCoveredFraction=1)]
        self.assertNotIn('\nD ', E.report_text(rows, key='k'))
        self.assertIn('2 samples', E.report_text(rows, key='k', members={'A', 'B'}))
        self.assertEqual(E.report_text([], key='k'), 'no successful measurement to report')


class TestOneSample(EffBase):
    def measure(self, s='S1', **kw):
        w = E.windows_from_fai(FAI, 1000, 0.4, 5_000_000)
        return E.measure_sample(self.env, CFG, ECFG, s, w, E.definition_key(CFG, ECFG, w), f'{self.local}stage/', self.sizes[s], 10.0, 'samtools 1.99; bcftools 1.99', threads=2, **kw)

    def test_the_numbers_of_a_sample(self):
        r = self.measure()
        self.assertEqual(r['Status'], 'ok', r['Error'])
        self.assertAlmostEqual(r['EffectiveCoverage'], 3.8); self.assertAlmostEqual(r['WindowRawDepth'], 5.5); self.assertAlmostEqual(r['BaseSurvivalFraction'], 3.8 / 5.5)
        self.assertAlmostEqual(r['WindowCoveredFraction'], 0.9); self.assertAlmostEqual(r['DuplicateFraction'], 0.291); self.assertAlmostEqual(r['MapqPassFraction'], 500 / 700)
        self.assertEqual((r['PrimaryReads'], r['TableCoverage'], r['GenomeVersion'], r['SampleID']), (1000, 10.0, G, 'S1')); self.assertEqual(r['BamBytes'], self.sizes['S1'])
        self.assertEqual(r['PrimaryMappedReads'], 995); self.assertAlmostEqual(r['DuplicationRate'], 291 / 995)

    def test_the_pileup_is_the_one_the_build_issues_over_the_windows(self):
        self.measure()
        mp = self.calls('bcftools')
        mp = [c for c in mp if c[0] == 'mpileup'][0]
        w = E.windows_from_fai(FAI, 1000, 0.4, 5_000_000)
        want = B.mpileup_cmd(self.fa, mp[mp.index('-b') + 1], E.regions_arg(w), CFG)[1:]
        self.assertEqual(mp, want, 'the same function builds both commands: nothing can drift')
        for flag in ('--ignore-RG', '-q', '13', '-Q', '20', '-a', 'FORMAT/AD,FORMAT/DP'):
            self.assertIn(flag, mp)
        self.assertEqual(mp[mp.index('-r') + 1], 'NC_1:8000001-8001000,NC_2:6000001-6001000,NC_3:3200001-3201000')

    def test_the_whole_bam_counts_use_the_callers_mapq_and_the_read_filter(self):
        self.measure()
        views = [c for c in self.calls('samtools') if c[0] == 'view']
        self.assertEqual(len(views), 2)
        self.assertTrue(all('-c' in v and v[v.index('-F') + 1] == '3844' for v in views))
        self.assertEqual([v[v.index('-q') + 1] for v in views if '-q' in v], ['13'])

    def test_the_bam_is_deleted_afterwards_unless_kept(self):
        self.measure('S1'); self.assertFalse(os.path.exists(f'{self.local}stage/S1'))
        self.measure('S2', keep_bam=True); self.assertTrue(os.path.exists(f'{self.local}stage/S2/S2.all.bam'))

    def test_a_failure_is_a_row_and_not_an_exception(self):
        self.put('S1', dict(GOOD, fail='flagstat'))
        r = self.measure('S1'); self.assertEqual(r['Status'], 'failed'); self.assertIn('simulated failure of flagstat', r['Error']); self.assertFalse(os.path.exists(f'{self.local}stage/S1'), 'nothing is left on the disk')
        self.put('S2', dict(GOOD, fail='mpileup'))
        r = self.measure('S2'); self.assertEqual(r['Status'], 'failed'); self.assertIn('the pileup of S2', r['Error'])
        shutil.rmtree(f'{self.cloud}Bamfiles/{G}/S3')
        r = self.measure('S3'); self.assertEqual(r['Status'], 'failed'); self.assertIn('S3', r['Error'])
        self.put('S4', GOOD, raw_text='CORRUPT')
        r = self.measure('S4'); self.assertEqual(r['Status'], 'failed'); self.assertIn('damaged', r['Error'])

    def test_a_window_with_no_covered_position_fails_with_the_reason(self):
        self.put('S1', dict(GOOD, cover=0.0))
        r = self.measure('S1'); self.assertEqual(r['Status'], 'failed'); self.assertIn('no position of the windows is covered', r['Error'])

    def test_a_downloaded_bam_the_wrong_size_is_refused(self):
        w = E.windows_from_fai(FAI, 1000, 0.4, 5_000_000)
        r = E.measure_sample(self.env, CFG, ECFG, 'S1', w, 'k', f'{self.local}stage/', self.sizes['S1'] + 1, None, '')
        self.assertEqual(r['Status'], 'failed'); self.assertIn('the ledger says', r['Error'])


class TestRun(EffBase):
    def run_all(self, plan, **kw):
        msgs = []
        res = E.run_effcov(self.env, CFG, plan, log=msgs.append, **kw)
        return res, msgs

    def test_every_sample_is_measured_and_written_as_it_finishes(self):
        plan = self.plan()
        res, msgs = self.run_all(plan)
        self.assertEqual((res['ok'], res['failed'], res['skipped']), (4, [], 0))
        rows = E.read_rows(self.out)
        self.assertEqual(sorted(r['SampleID'] for r in rows), self.ids); self.assertTrue(all(r['Status'] == 'ok' and r['DefinitionKey'] == plan.key for r in rows))
        self.assertTrue(all(abs(float(r['EffectiveCoverage']) - 3.8) < 1e-6 for r in rows)); self.assertEqual(len(open(self.out).read().splitlines()), 5)
        self.assertEqual(len([m for m in msgs if m.startswith('[effcov] ') and '/4 ' in m]), 4); self.assertIn('4 measured, 0 failed', msgs[-1])
        self.assertEqual(os.listdir(plan.staged), [], 'no BAM is left behind')

    def test_a_failed_sample_does_not_stop_the_others_and_is_retried_by_a_rerun(self):
        self.put('S3', dict(GOOD, fail='mpileup'))
        res, _ = self.run_all(self.plan())
        self.assertEqual(res['ok'], 3); self.assertEqual([s for s, _ in res['failed']], ['S3'])
        self.assertEqual(len(E.read_rows(self.out)), 4)
        self.put('S3', GOOD)                                                          # the cause is fixed
        plan = self.plan(); self.assertEqual((sorted(plan.done), plan.todo), (['S1', 'S2', 'S4'], ['S3']))
        res, _ = self.run_all(plan)
        self.assertEqual((res['ok'], res['skipped']), (1, 3)); self.assertEqual(len(E.read_rows(self.out)), 5, 'the failed row stays as history')
        self.assertEqual(sorted(E.latest_ok(E.read_rows(self.out), plan.key)), self.ids)
        self.assertEqual(len([c for c in self.calls('rclone') if c[0] == 'copy']), 5, 'four downloads, then one more for the retry: nothing is downloaded twice')

    def test_a_changed_definition_measures_everything_again_and_redo_forces_it(self):
        self.run_all(self.plan())
        again = self.plan(ecfg=dict(ECFG, window_bp=2000)); self.assertEqual((len(again.done), len(again.todo)), (0, 4)); self.assertNotEqual(again.key, self.plan().key)
        same = self.plan(); self.assertEqual((len(same.done), len(same.todo)), (4, 0))
        redo = self.plan(redo=True); self.assertEqual((len(redo.done), len(redo.todo)), (0, 4))

    def test_the_rows_keep_the_tables_coverage_for_comparison(self):
        led = self.ledger(); led['Coverage'] = led['Coverage'].astype(object); led.loc[led.SampleID == 'S2', 'Coverage'] = 9.9; led.loc[led.SampleID == 'S4', 'Coverage'] = ''
        plan = self.plan(ledger=led); self.run_all(plan)
        got = {r['SampleID']: r['TableCoverage'] for r in E.read_rows(self.out)}
        self.assertEqual((got['S1'], got['S2'], got['S4']), ('10', '9.9', ''))


class TestPlanAndPreflight(EffBase):
    def test_the_plan_counts_and_says_how_much_will_be_downloaded(self):
        led = self.ledger(['S1', 'S2'])                                              # S3 and S4 are not in the AlignmentDatabase
        p = self.plan(ledger=led)
        self.assertEqual(sorted(p.unknown), ['S3', 'S4']); self.assertIsNone(p.sizes['S3']); self.assertEqual(p.sizes['S1'], self.sizes['S1'])
        text = E.format_plan(p, [], CFG)
        self.assertIn('4 samples (0 already measured with this definition, 4 to measure)', text); self.assertIn(p.key, text); self.assertIn('2 samples are not in the AlignmentDatabase', text)
        self.assertIn('MAPQ >= 13, base quality >= 20', text); self.assertIn('3 windows of 0.001 Mb', text); self.assertIn('ready: add --run to measure', text)

    def test_a_clean_plan_has_no_problems(self):
        self.assertEqual(E.preflight(self.env, self.plan(), free_bytes=10 ** 12), [])

    def test_problems_are_named(self):
        p = self.plan()
        self.assertTrue(any('samtools is not on the PATH' in b for b in E.preflight(self.env, p, free_bytes=10 ** 12, which=lambda t: None if t == 'samtools' else '/x/' + t)))
        self.assertTrue(any('awk is not on the PATH' in b for b in E.preflight(self.env, p, free_bytes=10 ** 12, which=lambda t: None if t == 'awk' else '/x/' + t)))
        bad = E.preflight(B.Env(self.local, self.cloud, G, '/no/such.fa'), p, free_bytes=10 ** 12)
        self.assertTrue(any('reference FASTA is not at /no/such.fa' in b for b in bad))
        self.assertTrue(any('disk has 0 GB free' in b for b in E.preflight(self.env, self.plan(ecfg=dict(ECFG, default_bam_bytes=15e9)), free_bytes=10)))
        self.assertTrue(any('no samples are selected' in b for b in E.preflight(self.env, self.plan(members=[]), free_bytes=10 ** 12)))

    def test_the_disk_need_is_the_largest_jobs_bams_and_not_all_of_them(self):
        big = pd.DataFrame([dict(SampleID=f'X{i}', GenomeVersion=G, Coverage=10, BamSize=10e9) for i in range(30)])
        p = E.make_plan(self.env, CFG, dict(ECFG, jobs=3), [f'X{i}' for i in range(30)], big, self.out)
        E.preflight(self.env, p, free_bytes=10 ** 12)
        self.assertAlmostEqual(p.need_bytes, 3 * 10e9 * 1.2, delta=1)

    def test_stage_bam_keeps_a_sound_bam_and_replaces_nothing(self):
        stage = f'{self.local}st/'
        self.assertEqual(B.stage_bam(self.env, stage, 'S1', self.sizes['S1']), 'downloaded')
        self.assertEqual(B.stage_bam(self.env, stage, 'S1', self.sizes['S1']), 'kept'); self.assertTrue(os.path.lexists(f'{stage}S1/S1.all.bam.bai'))
        self.assertEqual(len([c for c in self.calls('rclone') if c[0] == 'copy']), 1)


class CliBase(EffBase):
    def setUp(self):
        super().setUp()
        rows = [fish(s, 'M' if i % 2 else 'F', eco='Deep_Benthic' if i < 3 else 'Mbuna') for i, s in enumerate(self.ids)]
        self.paths = dict(sample=os.path.join(self.tmp, 's.tsv'), ledger=os.path.join(self.tmp, 'l.tsv'), ids=os.path.join(self.tmp, 'ids.txt'), cohorts=os.path.join(self.tmp, 'c.yaml'))
        sheet(rows).to_csv(self.paths['sample'], sep='\t', index=False); self.ledger().to_csv(self.paths['ledger'], sep='\t', index=False)
        open(self.paths['ids'], 'w').write('# a comment\nS1\nS2  # trailing comment\n\nS1\n')
        open(self.paths['cohorts'], 'w').write("defaults:\n  cohort: {genome: Mzebra_GT3_NCBI, caller: bcftools, groups: per_sample}\ncallsets:\n  deep:\n    kind: sex_set\n    species: Aulonocara Nhkata Yellow Head\n    ecogroups: [Deep_Benthic]\n")
        self.common = ['--sample-table', self.paths['sample'], '--alignment-table', self.paths['ledger'], '--local-master', self.local, '--cloud-master', self.cloud, '--fasta', self.fa,
                       '--out', self.out, '--cohorts', self.paths['cohorts'], '--config', os.path.join(self.tmp, 'none.yaml')]
        open(os.path.join(self.tmp, 'none.yaml'), 'w').write('window_bp: 1000\njobs: 2\nthreads: 2\n')

    def call(self, *args):
        out, err = io.StringIO(), io.StringIO()
        rc = effcov_cli.main(list(args) + self.common, out=out, err=err, free_bytes=10 ** 12)
        return rc, out.getvalue(), err.getvalue()



class TestCli(CliBase):
    def test_the_default_is_the_plan_and_nothing_is_downloaded(self):
        rc, out, err = self.call('measure', '--samples-file', self.paths['ids'])
        self.assertEqual(rc, 0, err + out); self.assertIn('2 samples (0 already measured with this definition, 2 to measure)', out); self.assertIn('add --run', out)
        self.assertEqual(self.calls('rclone'), []); self.assertFalse(os.path.exists(self.out))

    def test_run_measures_writes_and_reports(self):
        rc, out, err = self.call('measure', '-s', 'S1', 'S2', '--run')
        self.assertEqual(rc, 0, err + out); self.assertEqual(sorted(r['SampleID'] for r in E.read_rows(self.out)), ['S1', 'S2'])
        self.assertIn('[effcov] 2 measured, 0 failed', out); self.assertIn('2 samples; effective coverage 3.80 to 3.80x', out); self.assertIn('2 below 5x', out); self.assertIn('0.38', out)
        rc, out2, _ = self.call('measure', '-s', 'S1', 'S2', '--run')
        self.assertEqual(rc, 0); self.assertIn('2 already measured', out2); self.assertEqual(len(E.read_rows(self.out)), 2, 'nothing is measured twice')

    def test_a_failure_gives_exit_code_1_and_a_message(self):
        self.put('S2', dict(GOOD, fail='mpileup'))
        rc, out, _ = self.call('measure', '-s', 'S1', 'S2', '--run')
        self.assertEqual(rc, 1); self.assertIn('FAILED', out); self.assertIn('1 samples failed', out); self.assertEqual(len(E.read_rows(self.out)), 2)

    def test_selection_by_cohort_and_by_all_with_an_ecogroup_left_out(self):
        rc, out, err = self.call('measure', '--cohort', 'deep')
        self.assertEqual(rc, 0, err + out); self.assertIn('3 samples', out)
        rc, out, err = self.call('measure', '--all', '--exclude-ecogroup', 'Mbuna')
        self.assertEqual(rc, 0, err + out); self.assertIn('3 samples', out); self.assertIn('1 samples of Ecogroup Mbuna are left out', out)
        rc, out, _ = self.call('measure', '--all'); self.assertIn('4 samples', out)

    def test_usage_errors_exit_2(self):
        self.assertEqual(self.call('measure')[0], 2)
        self.assertEqual(self.call('measure', '--all', '-s', 'S1')[0], 2)
        rc, _, err = self.call('measure', '--cohort', 'nope'); self.assertEqual(rc, 2); self.assertIn("unknown cohort 'nope'", err)
        out, err = io.StringIO(), io.StringIO()
        self.assertEqual(effcov_cli.main(['measure', '-s', 'S1', '--local-master', self.local], out=out, err=err), 2)

    def test_a_problem_in_the_plan_stops_a_run_with_exit_3(self):
        out, err = io.StringIO(), io.StringIO()
        rc = effcov_cli.main(['measure', '-s', 'S1', '--run'] + self.common, out=out, err=err, free_bytes=10)
        self.assertEqual(rc, 3); self.assertIn('PROBLEM: the disk has', out.getvalue()); self.assertEqual(self.calls('rclone'), [])

    def test_report_reads_the_file(self):
        self.call('measure', '--all', '--run')
        rc, out, err = self.call('report')
        self.assertEqual(rc, 0, err); self.assertIn('4 samples; effective coverage 3.80', out); self.assertIn('sex M', out)
        rc, out, _ = self.call('report', '--cohort', 'deep'); self.assertIn('3 samples;', out)
        os.remove(self.out); self.assertEqual(self.call('report')[0], 2)

    def test_the_command_is_reachable_through_the_package_entry_point(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = cli.main(['effcov', 'measure', '-s', 'S1'] + self.common)
        self.assertEqual(rc, 0); self.assertIn('1 samples', out.getvalue().replace('1 samples (', '1 samples ('))


class TestGroupedReport(CliBase):
    def setUp(self):
        super().setUp()
        for s, d in (('S1', dict(GOOD, dups=500)), ('S2', dict(GOOD, dups=100)), ('S3', dict(GOOD, dups=300)), ('S4', dict(GOOD))):
            self.put(s, d)
        rows = [fish('S1', 'M', species='Mchenga conophoros'), fish('S2', 'F', species='Mchenga conophoros'), fish('S3', 'M', species='Aulonocara saulosi'), fish('S4', 'F', species='Aulonocara saulosi')]
        sheet(rows).to_csv(self.paths['sample'], sep='\t', index=False)
        led = self.ledger(); led['RunIDs'] = ['run1,,run2', 'run1', 'run3', 'run3']; led.to_csv(self.paths['ledger'], sep='\t', index=False)
        self.assertEqual(self.call('measure', '--all', '--run')[0], 0)

    def test_the_medians_by_species(self):
        rc, out, err = self.call('report', '--by', 'Species', '--summary-only')
        self.assertEqual(rc, 0, err); lines = out.splitlines()
        self.assertTrue(lines[0].startswith('Species')); self.assertFalse([l for l in lines if l.startswith('S1 ')], 'summary only: no per-sample lines')
        mc = [l for l in lines if l.startswith('Mchenga')][0]; sa = [l for l in lines if l.startswith('Aulonocara')][0]
        self.assertIn('30.2', mc, 'median of S1 50.3 and S2 10.1'); self.assertIn('10.1-50.3', mc); self.assertIn('29.2-30.2', sa, 'S4 29.2 and S3 30.2')

    def test_by_a_column_of_the_alignment_database(self):
        rc, out, err = self.call('report', '--by', 'RunIDs', '--summary-only')
        self.assertEqual(rc, 0, err); self.assertIn('run1,,run2', out); self.assertIn('run3', out)

    def test_the_groups_follow_the_per_sample_table_unless_summary_only(self):
        rc, out, _ = self.call('report', '--by', 'Species')
        self.assertIn('4 samples; effective coverage', out); self.assertIn('Species', out.split('4 samples;')[1])

    def test_an_unknown_column_is_refused(self):
        rc, out, err = self.call('report', '--by', 'NoSuchColumn')
        self.assertEqual(rc, 2); self.assertIn("no column 'NoSuchColumn'", err)

    def test_a_blank_group_is_labelled_not_dropped(self):
        sheet([fish('S1', 'M', species=''), fish('S2', 'F', species='X'), fish('S3', 'M', species='X'), fish('S4', 'F', species='X')]).to_csv(self.paths['sample'], sep='\t', index=False)
        rc, out, _ = self.call('report', '--by', 'Species', '--summary-only')
        self.assertIn('(blank)', out)


class TestLabReared(CliBase):
    def setUp(self):
        super().setUp()
        self.put('S5', GOOD); self.sizes.pop('S5')
        rows = [fish('S1', 'M', lab='1', eco='Deep_Benthic'), fish('S2', 'F', lab='1', eco='Mbuna'), fish('S3', 'M', lab='0'), fish('S4', 'F', lab='TRUE', eco='Utaka'), fish('S6', 'M', lab='1')]
        sheet(rows).to_csv(self.paths['sample'], sep='\t', index=False)               # S6 is LabReared but has no row in the AlignmentDatabase

    def test_the_lab_reared_samples_that_have_a_bam(self):
        rc, out, err = self.call('measure', '--lab-reared')
        self.assertEqual(rc, 0, err + out); self.assertIn('3 samples (0 already measured', out); self.assertIn('1 samples marked LabReared have no row in the AlignmentDatabase and are skipped', out)

    def test_it_combines_with_leaving_the_mbuna_out(self):
        rc, out, _ = self.call('measure', '--lab-reared', '--exclude-ecogroup', 'Mbuna')
        self.assertIn('2 samples (0 already', out); self.assertIn('1 samples of Ecogroup Mbuna are left out', out)

    def test_it_is_one_selector_among_the_others(self):
        self.assertEqual(self.call('measure', '--lab-reared', '--all')[0], 2); self.assertEqual(self.call('measure', '--lab-reared', '-s', 'S1')[0], 2)

    def test_a_sheet_without_the_column_says_so(self):
        sheet([fish('S1', 'M')]).drop(columns=['LabReared']).to_csv(self.paths['sample'], sep='\t', index=False)
        rc, out, err = self.call('measure', '--lab-reared')
        self.assertIn('no LabReared column', out); self.assertEqual(rc, 3)


class TestWrite(CliBase):
    HEAD = ['SampleID', 'GenomeVersion', 'Coverage', 'UsableCoverage', 'NRuns', 'DuplicationRate']          # the real header has DuplicationRate (the pipeline fills it for its own rows)

    def setUp(self):
        super().setUp()
        self.ws = FakeWorksheet(self.HEAD, [[s, G, '10', '6.8', '1', ''] for s in self.ids] + [['OTHER', G, '5', '', '1', '']])
        self.measured = self.call('measure', '--all', '--run')
        self.assertEqual(self.measured[0], 0, self.measured[1] + self.measured[2])

    def write(self, *args, ws=None, **kw):
        out, err = io.StringIO(), io.StringIO()
        rc = effcov_cli.main(['write'] + list(args) + self.common, out=out, err=err, worksheet_fn=lambda: (ws or self.ws), **kw)
        return rc, out.getvalue(), err.getvalue()

    def test_the_default_is_a_plan_and_writes_nothing(self):
        rc, out, err = self.write()
        self.assertEqual(rc, 0, err); self.assertEqual(self.ws.writes(), []); self.assertIn('4 samples measured with definition', out); self.assertIn('DuplicationRate: write 4', out)
        self.assertIn('columns to ADD at the end of the header: EffectiveCoverage, MapqPassFraction, BaseSurvivalFraction', out); self.assertIn('nothing was written: add --apply', out)
        self.assertEqual(self.ws.grid[0], self.HEAD)

    def test_apply_adds_the_columns_at_the_end_and_fills_the_measured_samples(self):
        rc, out, err = self.write('--apply')
        self.assertEqual(rc, 0, err + out); self.assertEqual(self.ws.grid[0], self.HEAD + ['EffectiveCoverage', 'MapqPassFraction', 'BaseSurvivalFraction'])
        self.assertEqual(self.ws.column('EffectiveCoverage'), [3.8, 3.8, 3.8, 3.8, '']); self.assertEqual(self.ws.column('MapqPassFraction'), [0.7143] * 4 + ['']); self.assertEqual(self.ws.column('BaseSurvivalFraction'), [0.6909] * 4 + [''])
        self.assertEqual(self.ws.column('DuplicationRate'), [0.2925] * 4 + [''], '291 duplicates / 995 mapped primary reads')
        self.assertEqual(self.ws.column('UsableCoverage'), ['6.8'] * 4 + ['']); self.assertEqual(self.ws.column('NRuns'), ['1'] * 5, 'no other column is touched')
        self.assertIn('wrote 16 cells', out); self.assertIn('every written cell was read back and is in the right row', out)
        j = [f for f in os.listdir(os.path.join(os.path.dirname(self.out), 'sheet_writes'))]; self.assertEqual(len(j), 1); self.assertTrue(j[0].startswith('alignment.Mzebra_GT3_NCBI.'))

    def test_a_second_apply_writes_nothing(self):
        self.write('--apply'); n = len(self.ws.calls)
        rc, out, _ = self.write('--apply')
        self.assertEqual(rc, 0); self.assertIn('nothing to write', out); self.assertEqual(len(self.ws.calls), n)

    def test_a_cell_with_another_value_is_left_alone_unless_overwrite(self):
        self.write('--apply'); c = self.ws.grid[0].index('EffectiveCoverage'); self.ws.grid[1][c] = 99.0
        rc, out, _ = self.write('--apply'); self.assertIn('conflict 1', out); self.assertEqual(self.ws.grid[1][c], 99.0)
        rc, out, _ = self.write('--apply', '--overwrite'); self.assertEqual(self.ws.grid[1][c], 3.8); self.assertIn('wrote 1 cells', out)

    def test_limit_is_a_first_trial_on_a_few_samples(self):
        rc, out, _ = self.write('--apply', '--limit', '2')
        self.assertEqual(rc, 0); self.assertEqual(self.ws.column('EffectiveCoverage'), [3.8, 3.8, '', '', '']); self.assertIn('wrote 8 cells', out)

    def test_the_pipelines_own_duplication_rate_is_kept_and_the_difference_is_shown(self):
        for r, v in ((1, '0.2931'), (2, '0.2925')):                                                # S1: Picard's slightly different figure; S2: the same value
            self.ws.grid[r][self.ws.grid[0].index('DuplicationRate')] = v
        rc, out, _ = self.write('--apply')
        self.assertEqual(rc, 0, out); self.assertEqual(self.ws.column('DuplicationRate'), ['0.2931', '0.2925', 0.2925, 0.2925, ''])
        self.assertIn("DuplicationRate: write 2, same 1, conflict 1", out); self.assertIn("sheet '0.2931', measured 0.2925", out)

    def _as_an_older_version_wrote_it(self):
        rows = E.read_rows(self.out)
        with open(self.out, 'w') as f:
            f.write('\t'.join(E.COLUMNS[:-2]) + '\n')
            for r in rows:
                f.write('\t'.join(r[c] for c in E.COLUMNS[:-2]) + '\n')

    def test_rows_from_an_older_version_still_write_the_other_columns_and_not_a_guess(self):
        self._as_an_older_version_wrote_it()
        rc, out, _ = self.write('--apply')
        self.assertEqual(rc, 0, out); self.assertEqual(self.ws.column('EffectiveCoverage'), [3.8, 3.8, 3.8, 3.8, '']); self.assertEqual(self.ws.column('DuplicationRate'), [''] * 5, 'no mapped count, no DuplicationRate')

    def test_the_plan_for_older_rows_says_they_will_be_measured_again(self):
        self._as_an_older_version_wrote_it()
        rc, out, _ = self.call('measure', '--all')
        self.assertIn('0 already measured with this definition, 4 to measure', out); self.assertIn('4 of the samples to measure were measured by an older version (no DuplicationRate)', out)
        rc, out, _ = self.call('measure', '--all', '--run')
        self.assertEqual(rc, 0, out); rc, out, _ = self.call('measure', '--all'); self.assertIn('4 already measured', out)
        self.assertEqual(len(E.read_rows(self.out)), 8, 'the older rows stay as history'); self.assertTrue(any(f.startswith('eff.tsv.bak-') for f in os.listdir(os.path.dirname(self.out))))

    def test_the_plan_is_alphabetical_and_says_what_a_limited_run_writes(self):
        rc, out, _ = self.write('--limit', '2')
        first = [l.split()[1] for l in out.splitlines() if l.startswith('    ') and '->' in l]
        self.assertEqual(first, sorted(first)); self.assertIn('write 16', out)
        self.assertIn('--limit 2: this run writes only 2 samples (S1, S2): 8 cells', out); self.assertEqual(self.ws.writes(), [])

    def test_a_selection_restricts_what_is_written(self):
        rc, out, _ = self.write('--apply', '-s', 'S1', 'S3')
        self.assertEqual(self.ws.column('EffectiveCoverage'), [3.8, '', 3.8, '', ''])

    def test_a_failed_measurement_is_never_written(self):
        self.put('S2', dict(GOOD, fail='mpileup')); os.remove(self.out)
        self.assertEqual(self.call('measure', '--all', '--run')[0], 1)                          # S2 fails now
        ws = FakeWorksheet(self.HEAD, [[s, G, '10', '', '1', ''] for s in self.ids])
        rc, out, _ = self.write('--apply', ws=ws)
        self.assertEqual(rc, 0); self.assertEqual(ws.column('EffectiveCoverage'), [3.8, '', 3.8, 3.8]); self.assertEqual(ws.column('DuplicationRate'), [0.2925, '', 0.2925, 0.2925])

    def test_samples_missing_from_the_sheet_are_reported(self):
        ws = FakeWorksheet(self.HEAD, [[s, G, '10', '', '1', ''] for s in ['S1', 'S3']])
        rc, out, _ = self.write(ws=ws)
        self.assertIn('the sample is not in the sheet: 2 samples (S2, S4)', out)

    def test_nothing_is_written_without_a_measurement_of_the_current_definition(self):
        os.remove(self.out)
        rc, _, err = self.write(); self.assertEqual(rc, 2); self.assertIn('no rows', err)
        lines = open(os.path.join(self.tmp, 'x.tsv'), 'w'); lines.write('\t'.join(E.COLUMNS) + '\n'); lines.write('\t'.join(['S1', G, 'ok', '', 'oldkey'] + [''] * (len(E.COLUMNS) - 5)) + '\n'); lines.close()
        out, err = io.StringIO(), io.StringIO()
        common = [a for a in self.common]; common[common.index('--out') + 1] = os.path.join(self.tmp, 'x.tsv')
        rc = effcov_cli.main(['write'] + common, out=out, err=err, worksheet_fn=lambda: self.ws)
        self.assertEqual(rc, 2); self.assertIn('no successful measurement under the current definition', err.getvalue()); self.assertEqual(self.ws.writes(), [])

    def test_the_live_sheet_is_needed_and_the_offline_flags_cannot_reach_it(self):
        out, err = io.StringIO(), io.StringIO()
        rc = effcov_cli.main(['write'] + self.common, out=out, err=err)
        self.assertEqual(rc, 2); self.assertIn('writing needs the live sheet', err.getvalue())

    def test_a_sheet_that_changes_while_writing_is_a_loud_failure(self):
        self.ws.after_write = lambda w: w.insert_row_at_top(['NEW', G, '1', '', '1'])
        rc, out, _ = self.write('--apply')
        self.assertEqual(rc, 1); self.assertIn('VERIFICATION FAILED', out); self.assertIn('the sheet changed while writing', out)

    def test_two_selectors_are_refused(self):
        self.assertEqual(self.write('--all', '--lab-reared')[0], 2)

    def test_the_sheet_updates_are_rounded_numbers_and_only_the_known_columns(self):
        ok = {'A': dict(SampleID='A', EffectiveCoverage='4.012345', MapqPassFraction='0.72601', BaseSurvivalFraction='', Seconds='9'), 'B': dict(EffectiveCoverage='x', MapqPassFraction='', BaseSurvivalFraction='')}
        self.assertEqual(E.sheet_updates(ok), {'A': {'EffectiveCoverage': 4.0123, 'MapqPassFraction': 0.726}})
        self.assertEqual(E.sheet_updates(ok, members={'B'}), {})


class TestTheRealSheetsHaveNaN(EffBase):
    """The first real run (2026-11-11) measured a sample correctly and then crashed in the report: the lab sheets give NaN, a float, for a blank cell, and the test tables had empty strings."""
    def tables(self):
        nan = float('nan')
        dt = pd.DataFrame(dict(SampleID=['S1', 'S2', 'S3', 'S4', nan], Sex=['M', nan, 'F', None, 'M'], Ecogroup=[nan, 'Deep_Benthic', 'Mbuna', 'Mbuna', nan]))
        led = pd.DataFrame(dict(SampleID=['S1', 'S2', 'S3', 'S4', nan], GenomeVersion=[G] * 5, Coverage=[10.0, nan, 9.9, 10.0, nan], BamSize=[float(self.sizes[s]) for s in self.ids] + [nan]))
        return dt, led

    def test_a_cell_is_its_text_whatever_gspread_made_of_it(self):
        import numpy as np
        t = effcov_cli._text
        self.assertEqual([t(1), t(1.0), t(0), t(True), t(np.bool_(True)), t(np.int64(1)), t(' M '), t('1')], ['1', '1.0', '0', 'True', 'True', '1', 'M', '1'])
        self.assertEqual([t(None), t(float('nan')), t(pd.NA), t('')], ['', '', '', ''])

    def test_lab_reared_read_as_numbers_booleans_or_text(self):
        """The live sheet gave LabReared as numbers and the first version found no sample at all."""
        nan = float('nan')
        for name, flags in (('ints', [1, 0, 1, nan]), ('floats', [1.0, 0.0, 1.0, nan]), ('bools', [True, False, True, None]), ('text', ['1', '0', 'TRUE', '']), ('words', ['yes', 'no', 'Yes', None])):
            dt = pd.DataFrame(dict(SampleID=['S1', 'S2', 'S3', 'S4'], Sex=['M', 'F', 'M', 'F'], LabReared=flags, Ecogroup=['Deep_Benthic'] * 4))
            led = self.ledger()
            a = effcov_cli.build_parser().parse_args(['measure', '--lab-reared'])
            members, notes = effcov_cli.select_members(a, dt, led)
            self.assertEqual(members, ['S1', 'S3'], name)

    def test_nothing_selected_says_what_the_column_holds(self):
        dt = pd.DataFrame(dict(SampleID=['S1', 'S2'], LabReared=['maybe', float('nan')]))
        a = effcov_cli.build_parser().parse_args(['measure', '--lab-reared'])
        members, notes = effcov_cli.select_members(a, dt, self.ledger())
        self.assertEqual(members, []); self.assertTrue(any("no sample has a true LabReared" in n and "'maybe' x1" in n for n in notes), notes)

    def test_groups_by_a_numeric_column_are_not_blank(self):
        dt = pd.DataFrame(dict(SampleID=['S1', 'S2', 'S3'], BroodID=[12, 12.0, float('nan')]))
        self.assertEqual(effcov_cli.group_map('BroodID', dt, None), {'S1': '12.0', 'S2': '12.0', 'S3': ''}, 'a float column: the same label for the same group, and a blank stays blank')

    def test_the_whole_path_with_tables_typed_like_gspread_delivers_them(self):
        """measure --lab-reared --run, report --by and write, on tables whose LabReared, BamSize and Coverage are numbers and whose blanks are NaN."""
        nan = float('nan')
        dt = pd.DataFrame(dict(SampleID=['S1', 'S2', 'S3', 'S4', 'S9'], Sex=['M', nan, 'F', 'M', 'F'], Species=['A sp', 'A sp', 'B sp', nan, 'B sp'], LabReared=[1, 1, 0, 1, 1], BroodID=[12, 12, nan, nan, 3],
                            Ecogroup=['Deep_Benthic', 'Deep_Benthic', 'Mbuna', 'Utaka', nan]))
        led = pd.DataFrame(dict(SampleID=['S1', 'S2', 'S3', 'S4'], GenomeVersion=[G] * 4, Coverage=[10.0, 9.9, nan, 10], BamSize=[int(self.sizes[s]) for s in self.ids], RunIDs=['r1', 'r1', 'r2', 'r2']))
        base = ['--local-master', self.local, '--cloud-master', self.cloud, '--fasta', self.fa, '--out', self.out, '--config', os.path.join(self.tmp, 'none.yaml')]
        ws = FakeWorksheet(['SampleID', 'GenomeVersion', 'Coverage', 'DuplicationRate'], [[s, G, '10', ''] for s in self.ids])
        with mock.patch('variant_calling.cohort_cli.load_tables', return_value=(dt, led)):
            out, err = io.StringIO(), io.StringIO()
            rc = effcov_cli.main(['measure', '--lab-reared', '--run'] + base, out=out, err=err, free_bytes=10 ** 12)
            self.assertEqual(rc, 0, err.getvalue() + out.getvalue()); self.assertIn('3 samples (0 already', out.getvalue()); self.assertIn('1 samples marked LabReared have no row in the AlignmentDatabase', out.getvalue())
            self.assertEqual(sorted(r['SampleID'] for r in E.read_rows(self.out)), ['S1', 'S2', 'S4'])
            for by in ('Species', 'RunIDs', 'BroodID', 'Ecogroup'):
                out = io.StringIO(); rc = effcov_cli.main(['report', '--by', by, '--summary-only'] + base, out=out, err=err)
                self.assertEqual(rc, 0, f'{by}: {err.getvalue()}'); self.assertIn(by, out.getvalue().splitlines()[0])
            out = io.StringIO(); rc = effcov_cli.main(['write', '--apply'] + base, out=out, err=err, worksheet_fn=lambda: ws)
        self.assertEqual(rc, 0, err.getvalue() + out.getvalue()); self.assertEqual(ws.column('EffectiveCoverage'), [3.8, 3.8, '', 3.8]); self.assertEqual(ws.column('DuplicationRate'), [0.2925, 0.2925, '', 0.2925])

    def test_the_sex_map_keeps_only_real_values(self):
        self.assertEqual(effcov_cli.sex_map(self.tables()[0]), {'S1': 'M', 'S3': 'F'})
        self.assertEqual(effcov_cli.sex_map(pd.DataFrame(dict(SampleID=['A']))), {}); self.assertEqual(effcov_cli.sex_map(None), {})

    def test_the_report_survives_nan_and_none_in_the_sex_map(self):
        rows = [dict(SampleID='A', Status='ok', DefinitionKey='k', EffectiveCoverage='3.8', TableCoverage='9.9', DuplicateFraction=0.29, MapqPassFraction=0.72, BaseSurvivalFraction=0.7, WindowCoveredFraction=0.9)]
        t = E.report_text(rows, key='k', sex={'A': float('nan'), 'B': None, 'C': 'F'})
        self.assertIn('1 samples', t); self.assertNotIn('nan', t)

    def test_all_ignores_blank_rows_of_the_ledger(self):
        dt, led = self.tables()
        a = effcov_cli.build_parser().parse_args(['measure', '--all', '--exclude-ecogroup', 'Mbuna'])
        members, notes = effcov_cli.select_members(a, dt, led)
        self.assertEqual(members, ['S1', 'S2'], 'no sample called nan, and the Mbuna are left out'); self.assertTrue(any('2 samples of Ecogroup Mbuna' in n for n in notes))

    def test_a_whole_run_on_nan_tables_ends_in_a_report_not_a_traceback(self):
        out, err = io.StringIO(), io.StringIO()
        with mock.patch('variant_calling.cohort_cli.load_tables', return_value=self.tables()):
            rc = effcov_cli.main(['measure', '--all', '--run', '--local-master', self.local, '--cloud-master', self.cloud, '--fasta', self.fa, '--out', self.out,
                                  '--config', os.path.join(self.tmp, 'none.yaml')], out=out, err=err, free_bytes=10 ** 12)
        self.assertEqual(rc, 0, err.getvalue() + out.getvalue()); self.assertEqual(sorted(r['SampleID'] for r in E.read_rows(self.out)), self.ids)
        self.assertIn('4 samples; effective coverage', out.getvalue()); self.assertIn('sex M: 1 samples', out.getvalue()); self.assertIn('sex F: 1 samples', out.getvalue())
        with mock.patch('variant_calling.cohort_cli.load_tables', return_value=self.tables()):
            out = io.StringIO()
            rc = effcov_cli.main(['report', '--local-master', self.local, '--cloud-master', self.cloud, '--fasta', self.fa, '--out', self.out], out=out, err=err)
        self.assertEqual(rc, 0); self.assertIn('4 samples;', out.getvalue())

    def setUp(self):
        super().setUp()
        open(os.path.join(self.tmp, 'none.yaml'), 'w').write('window_bp: 1000\njobs: 2\nthreads: 2\n')


if __name__ == '__main__':
    unittest.main()
