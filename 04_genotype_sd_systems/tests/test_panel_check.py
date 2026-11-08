"""The panel's pipeline-facing API: module resolution ('all'), the preflight check, and the command line exit codes."""
import os, subprocess, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import yaml
from cichlid_common.common_light import config_path, module_dirs, REPO_ROOT
from sd_genotypes.panel import measure as pm
import test_te_panel as T


class TestResolve(unittest.TestCase):
    def test_all_a_list_a_string_and_unknown_names(self):
        self.assertEqual(pm.resolve_modules('all'), list(pm.MODULES)); self.assertEqual(pm.resolve_modules(['all']), list(pm.MODULES))
        self.assertEqual(pm.resolve_modules('gsdf, te'), ['gsdf', 'te']); self.assertEqual(pm.resolve_modules(('lg5inv',)), ['lg5inv'])
        self.assertTrue({'gsdf', 'te', 'lg5inv', 'bchrom'} <= set(pm.MODULES))
        with self.assertRaises(ValueError):
            pm.resolve_modules('gsdf,nope')


class TestCheck(unittest.TestCase):
    def test_a_module_with_its_definition_is_ok_and_one_without_is_missing_with_the_build_command(self):
        cfg, _ = T.make_def()
        rows = {r['module']: r for r in pm.check(['gsdf', 'te'], cfgs={'te': cfg})}
        self.assertEqual((rows['gsdf']['status'], len(rows['gsdf']['key'])), ('OK', 12))
        self.assertEqual((rows['te']['status'], len(rows['te']['key'])), ('OK', 12))
        empty = dict(cfg); empty['definition_dir'] = tempfile.mkdtemp()
        r = pm.check(['te'], cfgs={'te': empty})[0]
        self.assertEqual(r['status'], 'MISSING'); self.assertIn('sd_genotypes.te.definition', r['message'])

    def test_a_broken_config_is_an_error_not_a_crash(self):
        cfg, _ = T.make_def(); broken = {k: v for k, v in cfg.items() if k != 'flanking'}        # a valid definition, a config without a required key
        r = pm.check(['te'], cfgs={'te': broken})[0]
        self.assertEqual(r['status'], 'ERROR'); self.assertIn('KeyError', r['message'])

    def test_the_command_line_exit_codes(self):
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(module_dirs()))
        run = lambda *a: subprocess.run([sys.executable, '-m', 'sd_genotypes'] + list(a), capture_output=True, text=True, cwd=REPO_ROOT, env=env)
        ok = run('panel', 'check', '--modules', 'gsdf')
        self.assertEqual(ok.returncode, 0, ok.stderr[-300:]); self.assertIn('gsdf', ok.stdout); self.assertIn('OK', ok.stdout)
        bad = run('panel', 'check', '--modules', 'nope')
        self.assertNotEqual(bad.returncode, 0); self.assertIn('unknown panel modules', bad.stdout + bad.stderr)
        use = run('panel', 'measure')
        self.assertNotEqual(use.returncode, 0); self.assertIn('--modules all', use.stdout + use.stderr)


if __name__ == '__main__':
    unittest.main()
