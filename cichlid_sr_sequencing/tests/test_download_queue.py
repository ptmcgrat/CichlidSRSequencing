"""Tests for helper_modules/download_queue.py with fake downloads (no Dropbox, no rclone).  python3 tests/test_download_queue.py"""
import os, subprocess, sys, tempfile, time, unittest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'helper_modules'))
import download_queue as dq


class FakeProc:
    """Pretends to be a Popen. Finishes after `duration` seconds; `rc` is its exit code; writes `content` to `path` when it finishes OK (or a corrupt partial if rc != 0)."""
    active = 0
    max_active = 0
    terminated = []

    def __init__(self, path, rc=0, duration=0.05, content=b'GOOD', partial=None):
        self.path, self.rc, self.content, self.partial = path, rc, content, partial
        self.end = time.time() + duration; self.pid = os.getpid() + 100000; self._done = False
        FakeProc.active += 1; FakeProc.max_active = max(FakeProc.max_active, FakeProc.active)

    def poll(self):
        if time.time() < self.end:
            return None
        if not self._done:
            self._done = True; FakeProc.active -= 1
            if self.rc == 0: open(self.path, 'wb').write(self.content)
            elif self.partial is not None: open(self.path, 'wb').write(self.partial)
        return self.rc

    def wait(self, timeout=None): return self.rc


def ok_if_good(path):
    if not os.path.exists(path): return False, 'file missing after download'
    return (open(path, 'rb').read() == b'GOOD'), 'content is not GOOD (corrupt)'


class T(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp() + '/'
        FakeProc.active = 0; FakeProc.max_active = 0; FakeProc.terminated = []
        self.logs = []
        self._orig_kill = dq.kill_tree
        dq.kill_tree = lambda p, grace=5: FakeProc.terminated.append(p.path)

    def tearDown(self):
        dq.kill_tree = self._orig_kill

    def q(self, start, **k):
        return dq.DownloadQueue(start, check_fn=ok_if_good, log=self.logs.append, poll=0.01, backoff=0.01, **k)

    def test_concurrency_is_bounded(self):
        q = self.q(lambda j: FakeProc(j.path, duration=0.05), max_parallel=3)
        jobs = [(f's{i}', self.d + f'f{i}.bam') for i in range(12)]
        self.assertEqual(q.run(jobs), {})
        self.assertTrue(all(os.path.exists(p) for _, p in jobs))
        self.assertLessEqual(FakeProc.max_active, 3); self.assertGreaterEqual(FakeProc.max_active, 2)

    def test_transient_failure_is_retried_and_corrupt_partial_removed(self):
        calls = {}
        def start(j):
            calls[j.path] = calls.get(j.path, 0) + 1
            return FakeProc(j.path, rc=1 if calls[j.path] == 1 else 0, partial=b'BAD')   # first attempt leaves a corrupt partial file
        self.assertEqual(self.q(start, max_parallel=2, attempts=3).run([('s', self.d + 'a.bam')]), {})
        self.assertEqual(calls[self.d + 'a.bam'], 2); self.assertEqual(open(self.d + 'a.bam', 'rb').read(), b'GOOD')

    def test_corrupt_download_with_exit_code_zero_is_caught_and_redone(self):
        calls = {}
        def start(j):
            calls[j.path] = calls.get(j.path, 0) + 1
            return FakeProc(j.path, content=b'TRUNCATED' if calls[j.path] == 1 else b'GOOD')
        self.assertEqual(self.q(start, attempts=3).run([('s', self.d + 'a.bam')]), {})
        self.assertEqual(calls[self.d + 'a.bam'], 2)

    def test_file_not_found_at_launch_is_retried_then_reported_without_crashing(self):
        calls = {}
        def start(j):
            calls[j.path] = calls.get(j.path, 0) + 1
            if 'missing' in j.path: raise FileNotFoundError('Cant find file for download: x')
            if 'flaky' in j.path and calls[j.path] < 3: raise FileNotFoundError('Cant find file for download: y')   # listing failed twice (network / throttling)
            return FakeProc(j.path)
        jobs = [('A', self.d + 'ok.bam'), ('B', self.d + 'flaky.bam'), ('C', self.d + 'missing.bam')]
        fails = self.q(start, attempts=3).run(jobs)
        self.assertEqual(list(fails), ['C']); self.assertIn('Cant find file', fails['C'][0]); self.assertIn('3 attempts', fails['C'][0])
        self.assertTrue(os.path.exists(self.d + 'ok.bam') and os.path.exists(self.d + 'flaky.bam') and calls[self.d + 'flaky.bam'] == 3)

    def test_every_file_of_a_sample_is_tracked(self):
        """Regression: processes[sample] = ... kept only the LAST file of a sample, so a failure of an earlier run's uBam went unnoticed."""
        start = lambda j: FakeProc(j.path, rc=1 if 'run1' in j.path else 0)
        fails = self.q(start, attempts=2).run([('S', self.d + 'run1.bam'), ('S', self.d + 'run2.bam')])
        self.assertEqual(list(fails), ['S']); self.assertEqual(len(fails['S']), 1); self.assertIn('run1.bam', fails['S'][0])
        self.assertTrue(os.path.exists(self.d + 'run2.bam'))

    def test_finished_files_are_skipped_and_incomplete_ones_replaced(self):
        open(self.d + 'done.bam', 'wb').write(b'GOOD'); open(self.d + 'half.bam', 'wb').write(b'HALF')
        started = []
        def start(j): started.append(os.path.basename(j.path)); return FakeProc(j.path)
        self.assertEqual(self.q(start).run([('a', self.d + 'done.bam'), ('b', self.d + 'half.bam')]), {})
        self.assertEqual(started, ['half.bam']); self.assertEqual(open(self.d + 'half.bam', 'rb').read(), b'GOOD')

    def test_interrupt_terminates_everything_that_is_running(self):
        n = {'k': 0}
        def start(j):
            n['k'] += 1
            if n['k'] == 4: raise KeyboardInterrupt
            return FakeProc(j.path, duration=5)
        with self.assertRaises(KeyboardInterrupt):
            self.q(start, max_parallel=6).run([(f's{i}', self.d + f'f{i}.bam') for i in range(6)])
        self.assertEqual(len(FakeProc.terminated), 3)

    def test_timeout_kills_a_stuck_download(self):
        q = dq.DownloadQueue(lambda j: FakeProc(j.path, duration=60), check_fn=ok_if_good, log=self.logs.append, poll=0.01, backoff=0.01, attempts=1, timeout=0.05)
        fails = q.run([('s', self.d + 'a.bam')])
        self.assertIn('timed out', fails['s'][0]); self.assertEqual(len(FakeProc.terminated), 1)

    def test_rclone_env_defaults_do_not_override_the_user(self):
        for k in dq.RCLONE_DEFAULTS: os.environ.pop(k, None)
        os.environ['RCLONE_TRANSFERS'] = '5'
        eff = dq.configure_rclone_env()
        self.assertEqual(eff['RCLONE_TRANSFERS'], '5'); self.assertEqual(eff['RCLONE_MULTI_THREAD_STREAMS'], '1')

    def test_bam_check_on_real_bgzf_files(self):
        """Real BGZF files built by hand (no samtools needed). A uBAM has no @SQ lines: it must still pass."""
        import struct, zlib
        def block(data):
            c = zlib.compressobj(6, zlib.DEFLATED, -15); cd = c.compress(data) + c.flush()
            return (b'\x1f\x8b\x08\x04' + b'\x00' * 4 + b'\x00\xff' + struct.pack('<H', 6) + b'BC' + struct.pack('<H', 2) + struct.pack('<H', len(cd) + 25)
                    + cd + struct.pack('<II', zlib.crc32(data) & 0xffffffff, len(data)))
        self.assertEqual(block(b''), dq.BGZF_EOF)
        good = block(os.urandom(20000)) + block(os.urandom(20000)) + block(b'')          # incompressible payloads: a byte offset inside a block is unambiguous
        p = self.d + 'u.bam'
        def put(data):
            with open(p, 'wb') as f: f.write(data)
        self.assertFalse(dq.bam_check(p)[0])                                         # missing
        put(b''); self.assertEqual(dq.bam_check(p), (False, 'empty file'))
        put(good); self.assertEqual(dq.bam_check(p), (True, ''))                      # complete BAM passes (no @SQ lines needed)
        put(good[:len(good) // 2]); ok, why = dq.bam_check(p)                         # truncated mid-file
        self.assertFalse(ok); self.assertIn('end-of-file marker', why)
        put(good[:-28]); self.assertFalse(dq.bam_check(p)[0])                         # cut exactly at a block boundary: only the EOF marker tells
        bad = bytearray(good); bad[10000] ^= 0xFF                                     # a flipped byte inside a block's payload: only gzip -t tells
        put(bytes(bad)); ok, why = dq.bam_check(p)
        self.assertFalse(ok); self.assertIn('gzip -t', why)
        self.assertTrue(dq.bam_check(p, full=False)[0])                               # (the fast check alone cannot see it)
        put(b'x' * 100); self.assertFalse(dq.bam_check(p)[0])                         # not a BAM at all


class RealProcessTree(unittest.TestCase):
    def test_kill_tree_kills_children(self):
        import psutil
        p = subprocess.Popen(['bash', '-c', 'sleep 60 & sleep 60; wait'])
        time.sleep(0.3)
        kids = psutil.Process(p.pid).children(recursive=True)
        self.assertTrue(len(kids) >= 1)
        dq.kill_tree(p, grace=3)
        time.sleep(0.2)
        self.assertTrue(all(not k.is_running() or k.status() == psutil.STATUS_ZOMBIE for k in kids) and p.poll() is not None)


if __name__ == '__main__':
    unittest.main(verbosity=2)
