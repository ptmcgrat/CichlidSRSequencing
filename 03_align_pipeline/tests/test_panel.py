"""panel.py: the command lines are the contract with module 04; the preflight stops a run that cannot measure; the sidecar says which BAM a panel JSON belongs to."""
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

    def test_the_check_is_run_with_the_contract_command_and_a_time_limit(self):
        seen = {}
        def run(cmd, **kw):
            seen.update(cmd=cmd, kw=kw); return Done(0, '')
        pn.preflight(self.cfg, run=run)
        self.assertEqual(seen['cmd'], pn.check_command(self.cfg)); self.assertEqual(seen['kw']['timeout'], 300)


class TestSidecar(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp() + '/'
        self.bam = self.dir + 'S1.all.bam'
        open(self.bam, 'wb').write(b'x' * 100)
        self.meta = self.dir + 'S1.panel.meta.json'
        self.cfg = pn.PanelConfig(python='p')

    def test_it_records_the_bam_and_the_runs_and_is_current_for_exactly_those(self):
        rec = pn.write_sidecar(self.meta, self.bam, 'S1', 'fp1', self.cfg, now=0)
        self.assertEqual((rec['bam'], rec['bam_size'], rec['input_fingerprint'], rec['measured_on'], rec['sample']), ('S1.all.bam', 100, 'fp1', '1970-01-01T00:00:00Z', 'S1'))
        self.assertEqual(json.load(open(self.meta)), rec); self.assertFalse(os.path.exists(self.meta + '.tmp'), 'written atomically')
        self.assertTrue(pn.sidecar_is_current(self.meta, self.bam, 'fp1'))

    def test_a_different_bam_other_runs_or_a_damaged_sidecar_are_not_current(self):
        pn.write_sidecar(self.meta, self.bam, 'S1', 'fp1', self.cfg)
        self.assertFalse(pn.sidecar_is_current(self.meta, self.bam, 'fp2'), 'the runs changed')
        open(self.bam, 'ab').write(b'more'); self.assertFalse(pn.sidecar_is_current(self.meta, self.bam, 'fp1'), 'the BAM was re-made: its size differs')
        self.assertFalse(pn.sidecar_is_current(self.dir + 'nothing.json', self.bam, 'fp1'), 'no sidecar: not current')
        open(self.meta, 'w').write('{not json'); self.assertFalse(pn.sidecar_is_current(self.meta, self.bam, 'fp1'))

    def test_an_unknown_fingerprint_matches_an_unknown_fingerprint(self):
        pn.write_sidecar(self.meta, self.bam, 'S1', None, self.cfg)
        self.assertTrue(pn.sidecar_is_current(self.meta, self.bam, ''))


if __name__ == '__main__':
    unittest.main()
