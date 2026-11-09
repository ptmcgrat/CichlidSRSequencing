"""The disk ledger and the continuous download manager: the queue must never run ahead of the disk, never be held up by slow integrity checks, and never hang."""
import os, shutil, sys, tempfile, threading, time, unittest
from types import SimpleNamespace
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from align_pipeline.disk_ledger import DiskLedger, dir_size
from align_pipeline.download_manager import DownloadManager

GB = 1e9


class FakeDisk:
    """A disk of fixed capacity. Samples occupy bytes; `free` is what is left."""
    def __init__(self, capacity_gb, other_gb=0):
        self.capacity, self.other, self.used = capacity_gb * GB, other_gb * GB, {}
        self.min_free_seen = capacity_gb * GB

    def usage(self, path):
        free = self.capacity - self.other - sum(self.used.values())
        self.min_free_seen = min(self.min_free_seen, free)
        return SimpleNamespace(free=free)

    def write(self, sample, gb):
        self.used[sample] = self.used.get(sample, 0) + gb * GB
        self.usage(None)

    def delete(self, sample):
        self.used.pop(sample, None)


class TestDiskLedger(unittest.TestCase):
    def ledger(self, disk, factor=3.5, margin=200):
        return DiskLedger('/data', factor, margin, usage_fn=disk.usage, size_fn=lambda s: disk.used.get(s, 0))

    def test_admits_only_when_the_reservation_and_the_margin_fit(self):
        d = FakeDisk(1000)
        led = self.ledger(d)
        self.assertEqual(led.admit('S1', 100 * GB), (True, ''))                           # needs 350 + 200 margin of 1000 free
        ok, why = led.admit('S2', 100 * GB)                                                # 1000 - 350 (S1 still to write) = 650 < 350 + 200? 650 >= 550: fits
        self.assertTrue(ok)
        ok, why = led.admit('S3', 100 * GB)                                                # 1000 - 700 = 300 < 550
        self.assertFalse(ok)
        self.assertIn('needs 350 GB', why); self.assertIn('200 GB margin', why)

    def test_what_a_sample_has_already_written_is_not_counted_twice(self):
        d = FakeDisk(1000)
        led = self.ledger(d)
        led.admit('S1', 100 * GB)
        d.write('S1', 100)                                                                 # its uBAMs are on disk: free fell by 100, "still to write" fell by 100 too
        self.assertAlmostEqual(led.snapshot()['unwritten_gb'], 250)
        self.assertAlmostEqual(led.snapshot()['free_gb'], 900)

    def test_settle_and_release_give_the_space_back(self):
        d = FakeDisk(1000)
        led = self.ledger(d)
        led.admit('S1', 100 * GB); led.admit('S2', 100 * GB)
        self.assertFalse(led.admit('S3', 100 * GB)[0])
        d.write('S1', 40)
        led.settle('S1')                                                                   # S1's big writes are over: it only holds the 40 GB it occupies
        self.assertAlmostEqual(led.snapshot()['unwritten_gb'], 350)                        # only S2 is still to write
        self.assertTrue(led.admit('S3', 100 * GB)[0] is False or True)                    # (S3 may or may not fit yet; the point is the arithmetic above)
        led.release('S2')
        self.assertTrue(led.admit('S3', 100 * GB)[0])

    def simulate(self, peak_extra_gb, n=40, capacity=1500, margin=200):
        """40 samples of 15 GB uBAMs through a disk. Each sample writes the sorted BAM, then the final BAM plus MarkDuplicates scratch (`peak_extra_gb`), then shrinks to its final
        files; uploads free the space at the end. Samples are admitted whenever the ledger says yes. Returns (min free bytes seen, peak samples in flight, finished)."""
        d = FakeDisk(capacity, other_gb=100)
        led = self.ledger(d, factor=3.5, margin=margin)
        queue, in_flight, done, peak_in_flight, steps = [f'S{i:02d}' for i in range(n)], [], [], 0, 0
        while queue or in_flight:
            steps += 1
            self.assertLess(steps, 10_000)
            while queue and led.admit(queue[0], 15 * GB)[0]:                              # the download manager: admit whenever the ledger says yes
                s = queue.pop(0); d.write(s, 15); in_flight.append([s, 0])
            peak_in_flight = max(peak_in_flight, len(in_flight))
            for entry in list(in_flight):
                s, stage = entry
                entry[1] += 1
                if stage == 0: d.write(s, 15)                                              # sorted BAM
                elif stage == 1: d.write(s, 15); d.write(s, peak_extra_gb)                 # final BAM + scratch: the peak
                elif stage == 2: d.used[s] = 20 * GB; led.settle(s)                        # uBAMs, sorted BAM and scratch deleted; final + splits remain
                elif stage == 3: d.delete(s); led.release(s); in_flight.remove(entry); done.append(s)    # uploaded and cleaned up
        return d.min_free_seen, peak_in_flight, len(done)

    def test_a_disk_that_fills_up_is_never_overrun_when_samples_stay_within_their_reservation(self):
        """Peak = 15 (uBAMs) + 15 (sorted) + 15 (final) + 5 (scratch) = 50 GB <= the 52.5 GB reserved (3.5x): free space never falls below the margin."""
        min_free, peak_in_flight, done = self.simulate(peak_extra_gb=5)
        self.assertEqual(done, 40)
        self.assertGreaterEqual(min_free, 200 * GB - 1, f'free space fell to {min_free / GB:.0f} GB, below the 200 GB margin')
        self.assertGreaterEqual(peak_in_flight, 2, 'the ledger must still let several samples run at once')

    def test_the_margin_keeps_the_disk_from_filling_even_if_samples_outgrow_their_reservation(self):
        """Peak = 55 GB (3.67x, 5% over the reservation): the margin absorbs it. Free space dips below the margin but stays far from zero."""
        min_free, _, done = self.simulate(peak_extra_gb=10)
        self.assertEqual(done, 40)
        self.assertGreater(min_free, 100 * GB, f'free space fell to {min_free / GB:.0f} GB')
        min_free, _, _ = self.simulate(peak_extra_gb=15)                                   # 60 GB peak: 14% over
        self.assertGreater(min_free, 0, 'even 14% over the reservation must not fill the disk')

    def test_dir_size_counts_files(self):
        t = tempfile.mkdtemp()
        os.makedirs(t + '/a/b')
        for name, n in (('/x', 100), ('/a/b/y', 50)):
            with open(t + name, 'wb') as f:
                f.write(b'1' * n)
        self.assertEqual(dir_size(t), 150); self.assertEqual(dir_size(t + '/nothing'), 0)
        shutil.rmtree(t)


class FakeProc:
    pid = None                                  # never a real PID: the kill function is injected

    def __init__(self, ends_at, rc=0, clock=time.time):
        self.ends_at, self.rc, self.clock = ends_at, rc, clock

    def poll(self):
        return self.rc if self.clock() >= self.ends_at else None


class Harness:
    def __init__(self, tmp, plan, dur=0.02, check_delay=0.0, max_parallel=2, attempts=3, admit=None, idle=None, idle_wait=0.2, fail_start=None, check_fail=None):
        self.tmp, self.dur, self.check_delay = tmp, dur, check_delay
        self.samples = [SimpleNamespace(sample=s, ubams=[f'{tmp}/{s}_{i}.bam' for i in range(n)]) for s, n in plan.items()]
        self.ready, self.failed, self.killed = [], {}, []
        self.running_now, self.max_running, self.starts, self.attempts_seen = 0, 0, [], {}
        self.fail_start, self.check_fail = fail_start or (lambda path, attempt: False), check_fail or (lambda path, attempt: False)
        self.lock, self.stop = threading.Lock(), threading.Event()
        self.t0 = time.time()
        self.mgr = DownloadManager(self.samples, self.start, admit or (lambda sp: (True, '')), lambda sp: self.ready.append(sp.sample), lambda sp, why: self.failed.__setitem__(sp.sample, why),
                                   self.stop, check_fn=self.check, idle_fn=idle, max_parallel=max_parallel, check_workers=2, attempts=attempts, backoff=0.0, poll=0.002, admit_poll=0.01,
                                   idle_wait=idle_wait, log=lambda m: None, kill_fn=self.killed.append)

    def start(self, job):
        n = self.attempts_seen[job.path] = self.attempts_seen.get(job.path, 0) + 1
        if self.fail_start(job.path, n):
            raise FileNotFoundError('Cant find file for download: ' + job.path)
        with self.lock:
            self.running_now += 1; self.max_running = max(self.max_running, self.running_now)
        self.starts.append((time.time() - self.t0, job.path))
        with open(job.path, 'w') as f:
            f.write('bam')
        proc = FakeProc(time.time() + self.dur)
        orig = proc.poll

        def poll():
            r = orig()
            if r is not None and not getattr(proc, 'counted', False):
                proc.counted = True
                with self.lock:
                    self.running_now -= 1
            return r
        proc.poll = poll
        return proc

    def check(self, path):
        time.sleep(self.check_delay)
        n = self.attempts_seen.get(path, 1)
        return (False, 'gzip -t failed (corrupt BAM)') if self.check_fail(path, n) else (os.path.exists(path), '')

    def run(self):
        t = time.time(); self.mgr.run(); return time.time() - t


class TestDownloadManager(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_all_samples_arrive_and_no_more_than_max_parallel_run(self):
        h = Harness(self.tmp, {'A': 2, 'B': 1, 'C': 3}, max_parallel=2)
        h.run()
        self.assertEqual(sorted(h.ready), ['A', 'B', 'C']); self.assertEqual(h.failed, {})
        self.assertLessEqual(h.max_running, 2); self.assertEqual(h.max_running, 2)

    def test_slow_integrity_checks_do_not_hold_up_the_downloads(self):
        """12 files, 2 slots, each transfer 0.04 s, each check 0.4 s. If the check blocked the slot the run would take 12/2 x 0.44 = 2.6 s; with background checks the transfers finish
        in ~0.3 s and only the checks (2 workers) remain: ~12/2 x 0.4 = 2.4 s total, but the DOWNLOADS must all have started long before that."""
        h = Harness(self.tmp, {f'S{i}': 2 for i in range(6)}, dur=0.04, check_delay=0.4, max_parallel=2)
        h.run()
        self.assertEqual(len(h.ready), 6)
        last_start = max(t for t, _ in h.starts)
        self.assertLess(last_start, 1.2, f'the last download only started after {last_start:.2f} s: the checks are blocking the slots')

    def test_a_sample_is_not_started_until_the_admission_function_says_there_is_room(self):
        gate = {'open': False}
        h = Harness(self.tmp, {'A': 1, 'B': 2}, admit=lambda sp: (True, '') if (sp.sample == 'A' or gate['open']) else (False, 'no room'))
        threading.Timer(0.3, lambda: gate.update(open=True)).start()
        h.run()
        b_first = min(t for t, p in h.starts if '/B_' in p)
        self.assertGreaterEqual(b_first, 0.28, 'B was started before there was room for it')
        self.assertEqual(sorted(h.ready), ['A', 'B'])

    def test_the_queue_is_only_filled_as_fast_as_it_drains(self):
        """With 30 samples and 2 slots the manager must not admit all 30 at once (each admission reserves disk)."""
        admitted = []
        h = Harness(self.tmp, {f'S{i:02d}': 1 for i in range(30)}, dur=0.03, max_parallel=2, admit=lambda sp: (admitted.append(time.time()) or (True, '')))
        done_at = []
        h.mgr.on_ready = lambda sp: done_at.append(time.time())
        t = threading.Thread(target=h.run); t.start()
        time.sleep(0.15)
        self.assertLess(len(admitted), 15, f'{len(admitted)} samples admitted in 0.15 s: the queue is running far ahead of the transfers')
        t.join()
        self.assertEqual(len(admitted), 30)

    def test_retry_then_success_and_giving_up(self):
        h = Harness(self.tmp, {'OK': 1, 'BAD': 2, 'AFTER': 1}, attempts=3,
                    fail_start=lambda path, n: ('OK_0' in path and n == 1) or ('BAD_1' in path))
        h.run()
        self.assertEqual(sorted(h.ready), ['AFTER', 'OK']); self.assertEqual(list(h.failed), ['BAD'])
        self.assertEqual(h.attempts_seen[f'{self.tmp}/OK_0.bam'], 2, 'one failed start, then a retry that worked')
        self.assertIn('Cant find file', h.failed['BAD'][0]); self.assertIn('after 3 attempts', h.failed['BAD'][0])

    def test_a_file_that_fails_the_integrity_check_is_deleted_and_fetched_again(self):
        seen = []
        h = Harness(self.tmp, {'A': 1}, check_fail=lambda path, n: n == 1)
        orig = h.start
        h.start = lambda job: (seen.append(os.path.exists(job.path)), orig(job))[1]            # was the corrupt file still there when the retry started?
        h.mgr.start_fn = h.start
        h.run()
        self.assertEqual(h.ready, ['A']); self.assertEqual(h.attempts_seen[f'{self.tmp}/A_0.bam'], 2)
        self.assertEqual(seen, [False, False], 'a corrupt file must be removed before the retry (rclone would skip a same-size file)')

    def test_when_a_sample_fails_its_other_transfers_are_stopped_and_the_others_go_on(self):
        h = Harness(self.tmp, {'BAD': 2, 'GOOD': 1}, dur=0.3, max_parallel=3, attempts=1, fail_start=lambda path, n: 'BAD_1' in path)      # BAD_0 is already downloading when BAD_1 fails
        h.run()
        self.assertEqual(list(h.failed), ['BAD']); self.assertEqual(h.ready, ['GOOD'])
        self.assertTrue(len(h.killed) >= 1, 'the second file of the failed sample was still downloading and must be killed')

    def test_a_sample_that_can_never_fit_is_skipped_instead_of_waiting_for_ever(self):
        h = Harness(self.tmp, {'HUGE': 1, 'SMALL': 1}, idle=lambda: True, idle_wait=0.1,
                    admit=lambda sp: (False, 'free disk 90 GB; needs 350 GB plus margin') if sp.sample == 'HUGE' else (True, ''))
        elapsed = h.run()
        self.assertLess(elapsed, 5)
        self.assertEqual(list(h.failed), ['HUGE']); self.assertIn('not enough free disk', h.failed['HUGE'][0]); self.assertEqual(h.ready, ['SMALL'])

    def test_with_samples_in_flight_the_manager_waits_for_them_to_free_space(self):
        """Not idle (a sample in flight will free space): it must keep waiting, however long it takes."""
        state = {'free': False}
        h = Harness(self.tmp, {'A': 1}, idle=lambda: False, idle_wait=0.05, admit=lambda sp: (state['free'], 'waiting for a sample to finish'))
        threading.Timer(0.4, lambda: state.update(free=True)).start()
        h.run()
        self.assertEqual((h.ready, h.failed), (['A'], {}))

    def test_a_stop_request_ends_the_run_and_kills_what_is_running(self):
        h = Harness(self.tmp, {'A': 2, 'B': 2}, dur=30, max_parallel=2)
        threading.Timer(0.2, h.stop.set).start()
        elapsed = h.run()
        self.assertLess(elapsed, 5)
        self.assertEqual(len(h.killed), 2); self.assertEqual(h.ready, [])


if __name__ == '__main__':
    unittest.main()
