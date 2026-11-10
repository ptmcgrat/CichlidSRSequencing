"""`python -m variant_calling cohort ...`: list, members and check, offline from TSV exports; exit codes and the JSON report."""
import contextlib, io, json, os, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_fixtures import *
from variant_calling import cli, cohort_cli

YAML = """
defaults:
  cohort: {genome: Mzebra_GT3_NCBI, caller: bcftools, groups: per_sample}
callsets:
  good:
    kind: family
    description: a clean brood
    broods: [B1]
  small:
    kind: family
    broods: [B2]
  pool:
    kind: sex_set
    species: Aulonocara Nhkata Yellow Head
    ecogroups: [Deep_Benthic]
  broken:
    kind: trio
    broods: [B1]
"""


class TestCohortCli(unittest.TestCase):
    def setUp(self):
        d = tempfile.mkdtemp()
        rows = brood('B1') + brood('B2', sire='S2_m', dam='D2_f', n_m=2, n_f=2) + [fish(f'DB{i}', 'M' if i % 2 else 'F', eco='Deep_Benthic') for i in range(8)]
        self.paths = dict(cohorts=os.path.join(d, 'c.yaml'), sample=os.path.join(d, 's.tsv'), ledger=os.path.join(d, 'l.tsv'), json=os.path.join(d, 'out.json'))
        open(self.paths['cohorts'], 'w').write(YAML)
        sheet(rows).to_csv(self.paths['sample'], sep='\t', index=False)
        ledger([r['SampleID'] for r in rows]).to_csv(self.paths['ledger'], sep='\t', index=False)

    def call(self, *args, common=True):
        out = io.StringIO(); err = io.StringIO()
        extra = ['--cohorts', self.paths['cohorts'], '--sample-table', self.paths['sample'], '--alignment-table', self.paths['ledger']] if common else ['--cohorts', self.paths['cohorts']]
        with contextlib.redirect_stderr(err):
            rc = cohort_cli.main(list(args) + extra, out=out)
        return rc, out.getvalue(), err.getvalue()

    def test_list_shows_kind_genome_caller_and_groups_without_any_table(self):
        rc, out, _ = self.call('list', common=False)
        self.assertEqual(rc, 0); self.assertIn('good', out); self.assertIn('family', out); self.assertIn('caller=bcftools', out); self.assertIn('groups=per_sample', out); self.assertIn('a clean brood', out)

    def test_members_prints_one_id_per_line(self):
        rc, out, err = self.call('members', 'good')
        self.assertEqual(rc, 0); ids = out.split(); self.assertEqual(len(ids), 26); self.assertIn('S_m', ids)
        rc, out, err = self.call('members', 'pool'); self.assertEqual(len(out.split()), 8)

    def test_members_of_an_unusable_definition_exits_2(self):
        rc, out, err = self.call('members', 'broken')
        self.assertEqual(rc, 2); self.assertIn('kind', err)

    def test_check_ok_exits_0_and_blocked_exits_1(self):
        rc, out, _ = self.call('check', 'good'); self.assertEqual(rc, 0); self.assertIn('RESULT: ok to build', out)
        rc, out, _ = self.call('check', 'small'); self.assertEqual(rc, 1); self.assertIn('RESULT: BLOCKED', out); self.assertIn('too_few_per_sex', out)

    def test_a_definition_error_is_a_blocked_check_not_a_crash(self):
        rc, out, _ = self.call('check', 'broken'); self.assertEqual(rc, 1); self.assertIn('bad_definition', out)

    def test_json_report_is_written(self):
        rc, out, _ = self.call('check', 'small', '--json', self.paths['json'])
        d = json.load(open(self.paths['json'])); self.assertEqual((d['cohort'], d['kind'], d['blocked']), ('small', 'family', True)); self.assertGreaterEqual(d['n_block'], 1)

    def test_an_unknown_cohort_exits_2_and_lists_the_names(self):
        rc, out, err = self.call('check', 'nope'); self.assertEqual(rc, 2); self.assertIn('good', err)
        rc, out, err = self.call('check'); self.assertEqual(rc, 2)

    def test_offline_mode_needs_both_tables(self):
        with self.assertRaises(SystemExit):
            cohort_cli.main(['check', 'good', '--cohorts', self.paths['cohorts'], '--sample-table', self.paths['sample']], out=io.StringIO())

    def test_the_genome_can_be_overridden(self):
        rc, out, _ = self.call('check', 'good', '--genome', 'Other_Genome')
        self.assertEqual(rc, 1); self.assertIn('genome_version_mismatch', out); self.assertIn('Other_Genome', out)

    def test_the_module_command_line_dispatches_cohort_without_the_shared_parser(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = cli.main(['cohort', 'list', '--cohorts', self.paths['cohorts']])
        self.assertEqual(rc, 0); self.assertIn('good', out.getvalue())

    def test_the_shipped_cohorts_list_without_error(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = cli.main(['cohort', 'list'])
        self.assertEqual(rc, 0)
        for n in ('bc_family', 'yh_broods', 'mc_broods'):
            self.assertIn(n, out.getvalue())


if __name__ == '__main__':
    unittest.main()
