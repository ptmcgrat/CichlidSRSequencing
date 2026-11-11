"""`effcov` (module 05): the effective coverage of BAMs. Pure arithmetic and parsing, the definition, the output file, one sample and a whole run with FAKE samtools / bcftools / rclone
(fake_tools.py describes a BAM by a few numbers, so every expected value here is exact), resume and failure, the preflight, and the command line.
The fakes establish that the right command lines are issued and the numbers are combined correctly. How real samtools and bcftools behave is what the first real run is for."""
import contextlib, io, json, os, shutil, stat, sys, tempfile, unittest
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pandas as pd
from cohort_fixtures import *
from variant_calling import build as B, cli, effcov as E, effcov_cli

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
        self.assertEqual(E.parse_flagstat(text), dict(primary=63075860, primary_duplicates=18341603))
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


class TestCli(EffBase):
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


if __name__ == '__main__':
    unittest.main()
