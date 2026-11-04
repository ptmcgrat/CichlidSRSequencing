"""AlignmentWorker.monitorProcesses: finishes, reports failed samples, removes the error file of successful ones, and does not spin a core while waiting.
python3 tests/test_monitor_processes.py"""
import os, sys, tempfile, time, types, unittest
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
tmp = tempfile.mkdtemp(); os.makedirs(tmp + '/helper_modules')
open(tmp + '/helper_modules/__init__.py', 'w').write('')
open(tmp + '/helper_modules/file_manager.py', 'w').write('class FileManager: pass\n')
open(tmp + '/helper_modules/Timer.py', 'w').write('class Timer:\n    def start(self, *a): pass\n    def stop(self): pass\n')
for f in ('download_queue.py', 'alignment_worker.py'):
    open(f'{tmp}/helper_modules/{f}', 'w').write(open(os.path.join(ROOT, 'helper_modules', f)).read())
sys.modules['pysam'] = types.ModuleType('pysam'); sys.path.insert(0, tmp)
from helper_modules.alignment_worker import AlignmentWorker


class StubFM:
    samples = []
    def readSampleDatabase(self): pass
    def readAlignmentDatabase(self): pass


class T(unittest.TestCase):
    def test_monitor(self):
        aw = AlignmentWorker('G', StubFM(), check_size=False)
        d = tempfile.mkdtemp() + '/'
        mk = lambda s, code: types.SimpleNamespace(sampleID=s, error_file=d + s + '.err', command=[sys.executable, '-c', code])
        cmds = [mk('good1', 'import time; time.sleep(1.3)'), mk('bad', 'import sys,time; time.sleep(1.3); sys.stderr.write("boom"); sys.exit(2)'), mk('good2', 'import time; time.sleep(0.2)'), mk('good3', 'pass')]
        cpu0, t0 = time.process_time(), time.time()
        bad = aw.monitorProcesses(cmds, num_parallel=2, verbose=False)
        cpu, wall = time.process_time() - cpu0, time.time() - t0
        self.assertEqual(bad, ['bad'])
        self.assertTrue(os.path.exists(d + 'bad.err') and 'boom' in open(d + 'bad.err').read())
        self.assertFalse(any(os.path.exists(d + s + '.err') for s in ('good1', 'good2', 'good3')))
        self.assertLess(cpu, 0.5 * wall, f'busy loop: {cpu:.2f}s CPU in {wall:.2f}s')


if __name__ == '__main__':
    unittest.main(verbosity=2)
