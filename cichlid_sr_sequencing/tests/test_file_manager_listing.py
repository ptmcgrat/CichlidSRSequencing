"""Tests for the patched FileManager.downloadData / _cloudObjects with a fake rclone (no Dropbox, no Google). python3 tests/test_file_manager_listing.py"""
import importlib.util, os, subprocess, sys, tempfile, types, unittest
from unittest import mock
HERE = os.path.dirname(os.path.abspath(__file__))
for name in ('gspread', 'gspread_dataframe', 'google', 'google.oauth2', 'google.oauth2.service_account'):
    sys.modules.setdefault(name, types.ModuleType(name))
sys.modules['gspread_dataframe'].get_as_dataframe = lambda *a, **k: None
sys.modules['gspread_dataframe'].set_with_dataframe = lambda *a, **k: None
sys.modules['google.oauth2.service_account'].Credentials = object
spec = importlib.util.spec_from_file_location('fm_patched', os.path.join(HERE, '..', 'helper_modules', 'file_manager.py'))
fmmod = importlib.util.module_from_spec(spec); spec.loader.exec_module(fmmod)

CLOUD = 'ptm_dropbox:/COS/x/'
LOCAL = '/tmp/lm/'


class FakeRclone:
    """Scripted rclone: lsf results are popped from `lsf_results` [(returncode, stdout, stderr), ...]; the last one repeats."""
    def __init__(self, lsf_results):
        self.lsf_results = list(lsf_results); self.lsf_calls = 0; self.popen_args = []; self.copy_calls = []

    def run(self, args, **k):
        if args[:2] == ['rclone', 'lsf']:
            self.lsf_calls += 1
            rc, out, err = self.lsf_results.pop(0) if len(self.lsf_results) > 1 else self.lsf_results[0]
            return types.SimpleNamespace(returncode=rc, stdout=out, stderr=err)
        if args[:2] == ['rclone', 'copy']:
            self.copy_calls.append(args)
            os.makedirs(args[3], exist_ok=True); open(os.path.join(args[3], os.path.basename(args[2].rstrip('/'))), 'w').write('x')
            return types.SimpleNamespace(returncode=0, stdout='', stderr='')
        raise AssertionError(args)

    def popen(self, args, **k):
        self.popen_args.append(args); return types.SimpleNamespace(args=args)


class T(unittest.TestCase):
    def mk(self, lsf_results, local=None):
        fm = fmmod.FileManager.__new__(fmmod.FileManager)            # skip __init__ (needs Google credentials)
        fm.localMasterDir = local or LOCAL; fm.cloudMasterDir = CLOUD
        self.rc = FakeRclone(lsf_results)
        patches = [mock.patch.object(fmmod.subprocess, 'run', self.rc.run), mock.patch.object(fmmod.subprocess, 'Popen', self.rc.popen), mock.patch.object(fmmod.time, 'sleep', lambda s: None)]
        for p in patches: p.start(); self.addCleanup(p.stop)
        return fm

    def test_transient_listing_failure_is_retried_not_reported_as_missing(self):
        fm = self.mk([(1, '', 'CRITICAL: get current account failed'), (1, '', 'Error (500)'), (0, 'ERR1.bam\nERR2.bam\n', '')])
        p = fm.downloadData(LOCAL + 'Reads/P/ERR2.bam', parallel=True)
        self.assertEqual(self.rc.lsf_calls, 3); self.assertEqual(p.args, ['rclone', 'copy', CLOUD + 'Reads/P/ERR2.bam', LOCAL + 'Reads/P/'])

    def test_persistent_listing_failure_is_a_runtime_error_with_the_rclone_message(self):
        fm = self.mk([(1, '', 'dial tcp: network is unreachable')])
        with self.assertRaises(RuntimeError) as c:
            fm.downloadData(LOCAL + 'Reads/P/ERR2.bam', parallel=True)
        self.assertIn('network is unreachable', str(c.exception)); self.assertNotIsInstance(c.exception, FileNotFoundError); self.assertEqual(self.rc.lsf_calls, 3)

    def test_file_really_absent_from_a_successful_listing_is_still_file_not_found(self):
        fm = self.mk([(0, 'ERR1.bam\n', '')])
        with self.assertRaises(FileNotFoundError) as c:
            fm.downloadData(LOCAL + 'Reads/P/ERR9.bam', parallel=True)
        self.assertIn('Cant find file for download', str(c.exception)); self.assertEqual(self.rc.lsf_calls, 1)       # fresh listing: no pointless second one

    def test_directory_not_found_exit_code_means_absent(self):
        fm = self.mk([(3, '', 'directory not found')])
        with self.assertRaises(FileNotFoundError):
            fm.downloadData(LOCAL + 'Reads/NOPE/ERR1.bam', parallel=True)

    def test_a_folder_is_listed_once_for_many_files(self):
        fm = self.mk([(0, ''.join(f'ERR{i}.bam\n' for i in range(50)), '')])
        for i in range(50):
            fm.downloadData(LOCAL + f'Reads/P/ERR{i}.bam', parallel=True)
        self.assertEqual(self.rc.lsf_calls, 1); self.assertEqual(len(self.rc.popen_args), 50)

    def test_stale_cache_is_refreshed_before_giving_up(self):
        fm = self.mk([(0, 'ERR1.bam\n', ''), (0, 'ERR1.bam\nERR2.bam\n', '')])
        fm.downloadData(LOCAL + 'Reads/P/ERR1.bam', parallel=True)
        p = fm.downloadData(LOCAL + 'Reads/P/ERR2.bam', parallel=True)         # uploaded after the first listing
        self.assertEqual(self.rc.lsf_calls, 2); self.assertEqual(p.args[2], CLOUD + 'Reads/P/ERR2.bam')

    def test_names_with_spaces_and_directories_are_listed_correctly(self):
        fm = self.mk([(0, 'my file.bam\nsubdir/\n', '')])
        self.assertEqual(fm._cloudObjects(CLOUD + 'x/'), ['my file.bam', 'subdir/'])
        fm.downloadData(LOCAL + 'x/my file.bam', parallel=True)                  # used to be split on whitespace and reported missing
        self.assertEqual(self.rc.popen_args[0][2], CLOUD + 'x/my file.bam')

    def test_non_parallel_download_and_directory_download_still_work(self):
        d = tempfile.mkdtemp() + '/'
        fm = self.mk([(0, 'a.txt\nsub/\n', '')], local=d)
        fm.downloadData(d + 'Reads/a.txt')
        self.assertTrue(os.path.exists(d + 'Reads/a.txt'))
        fm.downloadData(d + 'Reads/sub/')
        self.assertTrue(os.path.isdir(d + 'Reads/sub'))


if __name__ == '__main__':
    unittest.main(verbosity=2)
