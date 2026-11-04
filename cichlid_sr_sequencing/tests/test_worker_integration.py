"""Runs the PATCHED AlignmentWorker.downloadReadData against a stand-in FileManager that launches real subprocesses.
pysam is stubbed (not needed for downloading). python3 tests/test_worker_integration.py"""
import os, sys, tempfile, types, subprocess, unittest, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
tmp = tempfile.mkdtemp()
os.makedirs(tmp + '/helper_modules')
open(tmp + '/helper_modules/__init__.py', 'w').write('')
open(tmp + '/helper_modules/file_manager.py', 'w').write('class FileManager: pass\n')
open(tmp + '/helper_modules/Timer.py', 'w').write('class Timer:\n    def start(self, *a): pass\n    def stop(self): pass\n')
for f in ('download_queue.py', 'alignment_worker.py'):
    open(f'{tmp}/helper_modules/{f}', 'w').write(open(os.path.join(ROOT, 'helper_modules', f)).read())
sys.modules['pysam'] = types.ModuleType('pysam')
sys.path.insert(0, tmp)
from helper_modules.alignment_worker import AlignmentWorker


class StubFM:
    """Just enough FileManager for AlignmentWorker.__init__ and downloadReadData. 'Downloads' by running a real python subprocess."""
    def __init__(self, d):
        self.d = d; self.samples = ['S1', 'S2', 'S3', 'S4']; self.launches = {}; self.localRawBamFiles = []
    def readSampleDatabase(self): pass
    def readAlignmentDatabase(self): pass
    def createSampleFiles(self, s): self.localRawBamFiles = [f'{self.d}{s}.run{i}.bam' for i in range(2 if s == 'S2' else 1)]
    def downloadData(self, path, parallel=False):
        self.launches[path] = self.launches.get(path, 0) + 1
        n = self.launches[path]
        if 'S3' in path and n < 3:
            raise FileNotFoundError('Cant find file for download: ' + path)            # transient listing failure
        if 'S4' in path:
            raise FileNotFoundError('Cant find file for download: ' + path)            # really missing
        code = 'import sys,time; time.sleep(0.05); open(sys.argv[1],"wb").write(b"BAMDATA")'
        if 'S2.run0' in path and n == 1:
            code = 'import sys; sys.exit(1)'                                              # first run of S2 fails once, then succeeds
        return subprocess.Popen([sys.executable, '-c', code, path])


class T(unittest.TestCase):
    def test_patched_download_read_data(self):
        d = tempfile.mkdtemp() + '/'
        fm = StubFM(d)
        aw = AlignmentWorker('G', fm, check_size=False)
        t0 = time.time()
        import helper_modules.download_queue as dq
        dq.bam_check = lambda p: (os.path.exists(p) and os.path.getsize(p) > 0, 'bad')       # no samtools here: size check only
        orig = dq.DownloadQueue.__init__
        def init(self, *a, **k): k.setdefault('check_fn', dq.bam_check); k['poll'] = 0.02; orig(self, *a, **k)
        dq.DownloadQueue.__init__ = init
        bad = aw.downloadReadData(max_parallel=2, attempts=3, backoff=0.05)
        self.assertEqual(bad, ['S4'])                                                       # only the genuinely missing sample is reported, nothing raised
        for p in ('S1.run0.bam', 'S2.run0.bam', 'S2.run1.bam', 'S3.run0.bam'):
            self.assertEqual(open(d + p, 'rb').read(), b'BAMDATA', p)
        self.assertEqual(fm.launches[d + 'S3.run0.bam'], 3)                                 # listing failed twice, third try worked
        self.assertEqual(fm.launches[d + 'S2.run0.bam'], 2)                                 # exit code 1 once, retried
        self.assertIn('Cant find file', aw.download_failures['S4'][0])
        # a rerun downloads nothing that is already complete
        fm.launches.clear()
        bad2 = aw.downloadReadData(max_parallel=2, attempts=1, backoff=0.01)
        self.assertEqual(bad2, ['S4']); self.assertEqual([k for k in fm.launches if 'S4' not in k], [])


if __name__ == '__main__':
    unittest.main(verbosity=2)
