"""panel.py: the command lines are the contract with module 04; the preflight stops a run that cannot measure; the identity check says whether a panel JSON belongs to a BAM."""
import json, os, subprocess, sys, tempfile, unittest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from align_pipeline import panel as pn


class Done:
    def __init__(self, rc, out=''):
        self.returncode, self.stdout = rc, out


class TestCommands(unittest.TestCase):
    def test_the_measure_and_check_commands_are_exactly_the_contract(self):
        cfg = pn.PanelConfig(python='/env/bin/python3')
        self.assertEqual(pn.check_command(cfg), ['/env/bin/python3', '-m', 'sd_genotypes', 'panel', 'check', '--modules', 'all'])
        self.assertEqual(pn.measure_command(cfg, '/b/S.all.bam', 'S', '/b/S.panel.json'),
                         ['/env/bin/python3', '-m', 'sd_genotypes', 'panel', 'measure', '--bam', '/b/S.all.bam', '--sample', 'S', '--out-json', '/b/S.panel.json', '--modules', 'all'])
        self.assertEqual(pn.check_command(pn.PanelConfig(python='p', modules='gsdf,te'))[-2:], ['--modules', 'gsdf,te'])

    def test_a_checkout_that_is_not_installed_is_put_on_pythonpath_in_front_of_the_existing_one(self):
        old = os.environ.get('PYTHONPATH')
        try:
            os.environ.pop('PYTHONPATH', None)
            self.assertEqual(pn.PanelConfig(python='p', pythonpath='/repo/04').prefix(), ['env', 'PYTHONPATH=/repo/04'])
            os.environ['PYTHONPATH'] = '/other'
            self.assertEqual(pn.PanelConfig(python='p', pythonpath='/repo/04').prefix(), ['env', 'PYTHONPATH=/repo/04' + os.pathsep + '/other'])
        finally:
            os.environ.pop('PYTHONPATH', None)
            if old is not None:
                os.environ['PYTHONPATH'] = old
        self.assertEqual(pn.PanelConfig(python='p').prefix(), [])
        cmd = pn.check_command(pn.PanelConfig(python='p', pythonpath='/repo/04'))
        self.assertEqual((cmd[0], cmd[2]), ('env', 'p')); self.assertTrue(cmd[1].startswith('PYTHONPATH=/repo/04'))

    def test_the_default_interpreter_is_the_one_running_the_pipeline(self):
        self.assertEqual(pn.PanelConfig().python, sys.executable)


class TestPreflight(unittest.TestCase):
    cfg = pn.PanelConfig(python='p')

    def test_ok_returns_what_the_check_printed(self):
        self.assertEqual(pn.preflight(self.cfg, run=lambda cmd, **kw: Done(0, 'gsdf 5d71ee8d6c2b\nte 2af781e90119\n')), (True, 'gsdf 5d71ee8d6c2b\nte 2af781e90119'))

    def test_a_non_zero_exit_a_missing_interpreter_and_a_hang_all_stop_the_run(self):
        ok, text = pn.preflight(self.cfg, run=lambda cmd, **kw: Done(1, 'te: definition file missing; build it with ...'))
        self.assertFalse(ok); self.assertIn('definition file missing', text)
        def no_python(cmd, **kw):
            raise FileNotFoundError(2, 'No such file', cmd[0])
        ok, text = pn.preflight(self.cfg, run=no_python); self.assertFalse(ok); self.assertIn('could not start p', text)
        def hang(cmd, **kw):
            raise subprocess.TimeoutExpired(cmd, kw['timeout'])
        ok, text = pn.preflight(self.cfg, run=hang, timeout=5); self.assertFalse(ok); self.assertIn('did not answer within 5 s', text)

    def test_the_pilots_frozen_keys_are_checked_when_they_are_expected(self):
        printed = 'gsdf 5d71ee8d6c2b\nte 2af781e90119\nlg5inv 80bbea57e10b\nbchrom 75367d0dfba8'
        cfg = pn.PanelConfig(python='p', expect=['5d71ee8d6c2b', '2af781e90119', '80bbea57e10b', '75367d0dfba8'])
        self.assertEqual(pn.preflight(cfg, run=lambda cmd, **kw: Done(0, printed)), (True, printed))
        ok, text = pn.preflight(cfg, run=lambda cmd, **kw: Done(0, printed.replace('2af781e90119', 'ffffffffffff')))
        self.assertFalse(ok); self.assertIn('not the expected ones', text); self.assertIn("['2af781e90119']", text); self.assertIn('ffffffffffff', text, 'what was printed is shown')
        self.assertTrue(pn.preflight(pn.PanelConfig(python='p'), run=lambda cmd, **kw: Done(0, 'anything'))[0], 'without expected keys nothing is compared')
        self.assertFalse(pn.preflight(cfg, run=lambda cmd, **kw: Done(1, printed))[0], 'a failing check stays a failure whatever it printed')

    def test_the_check_is_run_with_the_contract_command_and_a_time_limit(self):
        seen = {}
        def run(cmd, **kw):
            seen.update(cmd=cmd, kw=kw); return Done(0, '')
        pn.preflight(self.cfg, run=run)
        self.assertEqual(seen['cmd'], pn.check_command(self.cfg)); self.assertEqual(seen['kw']['timeout'], 300)


class TestIdentity(unittest.TestCase):
    """The panel JSON of sd_genotypes schema 2 records `bam` and `bam_size_bytes`; the rule is the one of `python -m sd_genotypes.panel.identity`."""
    def setUp(self):
        self.dir = tempfile.mkdtemp() + '/'
        self.bam = self.dir + 'S1.all.bam'
        open(self.bam, 'wb').write(b'x' * 100)
        self.json = self.dir + 'S1.panel.json'

    def write(self, rec):
        json.dump(rec, open(self.json, 'w'))

    def test_the_same_name_and_size_is_the_same_bam(self):
        self.write({'schema': 2, 'bam': 'S1.all.bam', 'bam_size_bytes': 100})
        self.assertEqual(pn.identity_status(self.json, self.bam), ('same', ''))
        self.assertEqual(pn.identity_status({'schema': 2, 'bam': '/other/machine/S1.all.bam', 'bam_size_bytes': 100}, self.bam)[0], 'same', 'only the file name counts, not the machine\'s path')

    def test_a_realigned_sample_keeps_its_name_but_not_its_size(self):
        self.write({'schema': 2, 'bam': 'S1.all.bam', 'bam_size_bytes': 100})
        open(self.bam, 'ab').write(b'more reads')
        status, why = pn.identity_status(self.json, self.bam)
        self.assertEqual(status, 'different'); self.assertIn('measured on a BAM of 100 bytes, this one has 110', why)

    def test_another_file_name_is_a_different_bam(self):
        self.write({'schema': 2, 'bam': 'S2.all.bam', 'bam_size_bytes': 100})
        self.assertEqual(pn.identity_status(self.json, self.bam)[0], 'different')

    def test_no_recorded_size_is_unknown_and_schema_1_is_still_read(self):
        self.assertEqual(pn.identity_status({'schema': 1, 'sample': 'S1'}, self.bam)[0], 'unknown')
        self.assertEqual(pn.identity_status({'schema': 2, 'bam': 'S1.all.bam', 'bam_size_bytes': 'unknown'}, self.bam)[0], 'unknown')
        self.assertEqual(pn.identity_status({'schema': 2, 'bam': 'S1.all.bam', 'bam_size_bytes': 'many'}, self.bam)[0], 'unknown')
        self.assertEqual(pn.identity_status({'schema': 1, 'bam': '/x/S1.all.bam', 'bam_size': 100}, self.bam)[0], 'same', 'schema 1 files that have bam and bam_size')
        self.assertEqual(pn.identity_status({'schema': 1, 'bam': '/x/S1.all.bam', 'bam_size': 7}, self.bam)[0], 'different')

    def test_a_missing_or_damaged_json_is_bad(self):
        self.assertEqual(pn.identity_status(self.dir + 'nothing.json', self.bam)[0], 'bad')
        open(self.json, 'w').write('{not json'); status, why = pn.identity_status(self.json, self.bam)
        self.assertEqual(status, 'bad'); self.assertIn('cannot be read', why)


if __name__ == '__main__':
    unittest.main()
