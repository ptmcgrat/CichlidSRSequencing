'''tools/migrate_lab_files.sh, run in a throwaway git repository that imitates the old layout.'''
import os, shutil, subprocess, sys, tempfile, unittest
from cichlid_common.common_light import REPO_ROOT

OLD_FM = ("import platform, os\n\nclass FileManager():\n\tdef __init__(self):\n\t\tself.s_ID = 'S'\n\t\tself.a_ID = 'A' # Alignment Database\n\t\tself.dl_ID = 'X'\n\n"
          "\t# ---- Download database\n\tdef downloadManifest(self, tab, who = None, flush_seconds = 15.0):\n\t\ttry:\n\t\t\tfrom helper_modules import download_manifest as dm\n\t\texcept ImportError:\n\t\t\timport download_manifest as dm\n"
          "\t\twho = who or 'x'\n\t\treturn None\n\n\tdef _authenticateGS(self):\n\t\tpass\n")
OLD_FILES = ['helper_modules/ena_fetch.py', 'helper_modules/sra_fetch.py', 'helper_modules/ena_scheduler.py', 'helper_modules/grab_ena_v2.py', 'helper_modules/download_manifest.py', 'helper_modules/seq_pipeline.py',
             'downloadShortReadData_v2.py', 'downloadManifest.py', 'downloadShortReadData.py', 'grabENA.py', 'alignFastQReads_v2.py', 'unit_scripts/ena_bench.py', 'unit_scripts/patch_file_manager.py',
             'tests/test_ena_v2.py', 'tests/test_sra.py', 'YH_SD_Paper/yhsd/cli.py']


def sh(cmd, cwd, **kw):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, **kw)


@unittest.skipUnless(shutil.which('git') and shutil.which('bash'), 'needs git and bash')
class T(unittest.TestCase):
    def setUp(self):
        d = tempfile.mkdtemp(); self.r = d + '/repo'
        shutil.copytree(REPO_ROOT, self.r, ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '.git'))
        old = self.r + '/cichlid_sr_sequencing'
        for f in OLD_FILES:
            os.makedirs(os.path.dirname(f'{old}/{f}'), exist_ok=True); open(f'{old}/{f}', 'w').write('x\n')
        open(old + '/helper_modules/file_manager.py', 'w').write(OLD_FM)
        for c in (['git', 'init', '-q'], ['git', 'config', 'user.email', 't@t'], ['git', 'config', 'user.name', 't'], ['git', 'add', '-A'], ['git', 'commit', '-qm', 'base']):
            self.assertEqual(sh(c, self.r).returncode, 0)
        self.env = dict(os.environ, MIGRATE_NO_BUSY_CHECK='1')

    def mig(self, *args):
        return subprocess.run(['bash', self.r + '/tools/migrate_lab_files.sh'] + list(args), cwd=self.r, capture_output=True, text=True, env=self.env)

    def exists(self, p):
        return os.path.exists(self.r + '/cichlid_sr_sequencing/' + p)

    def test_check_and_dry_run_change_nothing(self):
        r = self.mig('check'); self.assertEqual(r.returncode, 0); self.assertIn('helper_modules/ena_fetch.py', r.stdout); self.assertIn('YH_SD_Paper', r.stdout)
        r = self.mig('retire-old-download'); self.assertEqual(r.returncode, 0); self.assertIn('would: git', r.stdout); self.assertTrue(self.exists('helper_modules/ena_fetch.py'))
        self.assertEqual(sh(['git', 'status', '--porcelain'], self.r).stdout, '')

    def test_execute_retires_the_download_code_and_leaves_the_alignment_code(self):
        self.assertEqual(self.mig('retire-old-download', '--execute').returncode, 0)
        for gone in ('helper_modules/ena_fetch.py', 'helper_modules/download_manifest.py', 'downloadShortReadData_v2.py', 'unit_scripts/ena_bench.py', 'tests/test_sra.py'):
            self.assertFalse(self.exists(gone), gone)
        for kept in ('helper_modules/seq_pipeline.py', 'helper_modules/file_manager.py', 'alignFastQReads_v2.py', 'downloadShortReadData.py'):
            self.assertTrue(self.exists(kept), kept)                                                        # module 03, the file manager and the old scripts are not touched by this command
        self.assertEqual(self.mig('retire-old-download', '--execute').returncode, 4)                        # the staged removals are uncommitted: refused until committed, as the instructions say
        self.assertEqual(sh(['git', 'commit', '-qm', 'retire'], self.r).returncode, 0)
        self.assertEqual(self.mig('retire-old-download', '--execute').returncode, 0)                        # idempotent

    def test_yhsd_and_the_old_scripts_move_to_legacy(self):
        self.assertEqual(self.mig('retire-yhsd', '--execute').returncode, 0); self.assertEqual(sh(['git', 'commit', '-qm', 'yhsd'], self.r).returncode, 0)
        self.assertEqual(self.mig('retire-old-scripts', '--execute').returncode, 0)
        self.assertTrue(os.path.exists(self.r + '/legacy/YH_SD_Paper_v34/yhsd/cli.py')); self.assertTrue(os.path.exists(self.r + '/legacy/grabENA.py'))
        self.assertFalse(self.exists('YH_SD_Paper')); self.assertFalse(self.exists('grabENA.py'))

    def test_a_dirty_tree_is_refused(self):
        open(self.r + '/cichlid_sr_sequencing/alignFastQReads_v2.py', 'a').write('edit\n')
        r = self.mig('retire-old-download', '--execute'); self.assertEqual(r.returncode, 4); self.assertTrue(self.exists('helper_modules/ena_fetch.py'))

    def test_running_jobs_are_refused_unless_bypassed(self):
        self.env.pop('MIGRATE_NO_BUSY_CHECK'); p = subprocess.Popen(['bash', '-c', 'exec -a "python3 downloadShortReadData_v2.py --manifest Blumer" sleep 30'])
        try:
            r = self.mig('retire-old-download', '--execute'); self.assertEqual(r.returncode, 3); self.assertIn('running', r.stdout)
        finally:
            p.kill()

    def test_the_file_manager_upgrade_and_the_shim(self):
        fm = self.r + '/cichlid_sr_sequencing/helper_modules/file_manager.py'
        r = self.mig('upgrade-file-manager'); self.assertEqual(r.returncode, 0); self.assertIn('+\t\t\tfrom download_reads import manifest as dm', r.stdout); self.assertEqual(open(fm).read(), OLD_FM)       # dry run: the diff, nothing written
        r = self.mig('upgrade-file-manager', '--execute'); self.assertEqual(r.returncode, 0, r.stderr); self.assertIn('methods upgraded', r.stdout)
        new = open(fm).read(); self.assertIn('from download_reads import manifest as dm', new); self.assertIn('from helper_modules import download_manifest as dm', new)
        self.assertTrue([f for f in os.listdir(os.path.dirname(fm)) if f.startswith('file_manager.py.bak-')])
        self.assertEqual(sh(['git', 'commit', '-qam', 'fm'], self.r).returncode, 0)
        r2 = self.mig('upgrade-file-manager', '--execute'); self.assertEqual(r2.returncode, 0); self.assertIn('already patched', r2.stdout)                  # idempotent
        code = 'from cichlid_common.file_manager import FileManagerClass; print(FileManagerClass().__name__)'
        r = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, cwd='/', env=dict(os.environ, PYTHONPATH=self.r + '/common'))
        self.assertEqual((r.returncode, r.stdout.strip()), (0, 'FileManager'), r.stderr)                      # the shim finds the (upgraded) file manager in the old place


if __name__ == '__main__':
    unittest.main()
