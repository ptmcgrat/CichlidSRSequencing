"""Continuous, space-aware downloads of the uBAMs of many samples.

Compared with DownloadQueue (download_queue.py, one blocking call per batch of files):
  * continuous: as soon as a download slot is free the next file starts, across samples; there is no "wait for the whole chunk" gap;
  * admission: a sample's files are only queued when `admit_fn(sample)` says there is room (a disk-slot limit AND the free-disk ledger, see disk_ledger.py), so the queue can never
    run ahead of the disk. When there is no room it waits and tries again; it never queues more than it can hold;
  * integrity checks (BGZF end marker + gzip -t, which read the whole file) run in background threads: the download slot is freed the moment the transfer ends;
  * every file is retried with exponential backoff (a failed start such as 'Cant find file', a non-zero exit, a file that fails the check); a sample whose file never arrives is reported
    failed and its other files are dropped; the other samples go on;
  * a sample for which there is never going to be room (nothing in flight that could free space) fails with a clear message after `idle_wait` seconds instead of waiting for ever;
  * a stop request kills every running transfer (and everything it started).
Uses kill_tree and bam_check from procs.py."""
import glob
import os
import random
import time
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace

from .procs import bam_check, kill_tree


def _log(msg):
    print('  [download] ' + time.strftime('%H:%M:%S') + ' ' + msg, flush=True)


class DownloadManager:
    def __init__(self, samples, start_fn, admit_fn, on_ready, on_failed, stop, check_fn=bam_check, idle_fn=None, max_parallel=4, check_workers=2, attempts=4, backoff=30.0,
                 poll=1.0, admit_poll=15.0, idle_wait=1800.0, log=_log, sleep=time.sleep, clock=time.time, kill_fn=kill_tree):
        """samples: objects with .sample and .ubams (list of local paths). start_fn(job) -> Popen-like (poll(), pid), job has .key (sample) and .path;
        admit_fn(sample_obj) -> (ok, reason); on_ready(sample_obj) / on_failed(sample_obj, [reasons]); idle_fn() -> True when nothing in flight could free space."""
        self.samples, self.start_fn, self.admit_fn, self.on_ready, self.on_failed, self.stop = samples, start_fn, admit_fn, on_ready, on_failed, stop
        self.check_fn, self.idle_fn = check_fn, idle_fn or (lambda: False)
        self.max_parallel, self.check_workers, self.attempts, self.backoff = max(1, int(max_parallel)), max(1, int(check_workers)), max(1, int(attempts)), backoff
        self.poll, self.admit_poll, self.idle_wait = poll, admit_poll, idle_wait
        self.log, self.sleep, self.clock, self.kill_fn = log, sleep, clock, kill_fn            # kill_fn(proc): tests inject a recorder so no real PID is ever signalled
        self.stats = dict(files_ok=0, retries=0, samples_ready=0, samples_failed=0)

    # ---------------------------------------------------------------- one file failed an attempt: clean up, retry later or give up on the sample
    def _fail(self, job, reason, pending):
        st = job.state
        for f in [job.path] + glob.glob(job.path + '.*partial'):
            if os.path.exists(f):
                os.remove(f)                                         # rclone would skip a corrupt file of the same size on the retry
        if st['failed']:
            return
        if job.attempt < self.attempts:
            job.not_before = self.clock() + self.backoff * 2 ** (job.attempt - 1) * (1 + 0.25 * random.random())
            self.stats['retries'] += 1
            pending.append(job)
            self.log(f'{os.path.basename(job.path)}: {reason}; attempt {job.attempt}/{self.attempts}, retrying in {job.not_before - self.clock():.0f}s')
            return
        st['failed'] = True
        st['reasons'].append(f'{os.path.basename(job.path)}: {reason} (after {job.attempt} attempts)')
        self.log(f'{os.path.basename(job.path)}: GIVING UP: {reason}')
        self.stats['samples_failed'] += 1
        self.on_failed(st['sample'], list(st['reasons']))

    def _ok(self, job):
        st = job.state
        st['n_ok'] += 1
        self.stats['files_ok'] += 1
        if st['failed']:
            return
        self.log(f'{os.path.basename(job.path)} ok ({st["sample"].sample}: {st["n_ok"]}/{st["n_total"]})')
        if st['n_ok'] == st['n_total']:
            self.stats['samples_ready'] += 1
            self.on_ready(st['sample'])

    # ---------------------------------------------------------------- the loop
    def run(self):
        todo = deque(self.samples)
        pending, running, checking = deque(), [], []
        checker = ThreadPoolExecutor(self.check_workers)
        last_admit, waiting_since, last_reason_log = -1e18, None, -1e18
        try:
            while (todo or pending or running or checking) and not self.stop.is_set():
                now = self.clock()
                for job in list(running):                              # a sample that has failed: stop its other transfers at once
                    if job.state['failed'] and job.proc.poll() is None:
                        self.kill_fn(job.proc)
                for job in list(running):                              # transfers that ended
                    rc = job.proc.poll()
                    if rc is None:
                        continue
                    running.remove(job)
                    if job.state['failed']:
                        self._fail(job, 'sample already failed', pending)
                    elif rc != 0:
                        self._fail(job, f'download process exited with code {rc}', pending)
                    else:
                        job.future = checker.submit(self.check_fn, job.path)      # the slot is free at once; the check reads the file in the background
                        checking.append(job)
                for job in list(checking):                             # checks that ended
                    if not job.future.done():
                        continue
                    checking.remove(job)
                    try:
                        ok, why = job.future.result()
                    except Exception as e:
                        ok, why = False, f'check crashed: {type(e).__name__}: {e}'
                    if ok:
                        self._ok(job)
                    else:
                        self._fail(job, why, pending)
                # admit the next sample only when the queue is about to run dry, and only when there is room for it
                if todo and len(pending) + len(running) < self.max_parallel and now - last_admit >= (self.admit_poll if waiting_since is not None else 0):
                    sp = todo[0]
                    last_admit = now
                    ok, why = self.admit_fn(sp)
                    if ok:
                        todo.popleft()
                        waiting_since = None
                        st = dict(sample=sp, n_total=len(sp.ubams), n_ok=0, failed=False, reasons=[])
                        if not sp.ubams:
                            st['failed'] = True
                            self.stats['samples_failed'] += 1
                            self.on_failed(sp, ['no uBAM files listed for this sample'])
                        for path in sp.ubams:
                            pending.append(SimpleNamespace(key=sp.sample, path=path, attempt=0, not_before=0.0, proc=None, future=None, state=st))
                    else:
                        if waiting_since is None:
                            waiting_since = now
                        if now - last_reason_log > 600:
                            self.log(f'waiting for room before starting {sp.sample}: {why}')
                            last_reason_log = now
                        if self.idle_fn() and now - waiting_since > self.idle_wait and not running and not pending and not checking:
                            todo.popleft()                             # nothing in flight can ever free the space: do not wait for ever
                            waiting_since = None
                            self.stats['samples_failed'] += 1
                            self.log(f'{sp.sample}: NOT ENOUGH FREE DISK and nothing in flight that could free any; skipped ({why})')
                            self.on_failed(sp, [f'not enough free disk: {why}'])
                for job in list(pending):                              # start what is allowed
                    if len(running) >= self.max_parallel:
                        break
                    if job.state['failed']:
                        pending.remove(job)
                        continue
                    if job.not_before > now:
                        continue
                    pending.remove(job)
                    job.attempt += 1
                    try:
                        job.proc = self.start_fn(job)
                        running.append(job)
                    except Exception as e:                             # e.g. FileNotFoundError('Cant find file for download') raised while starting
                        self._fail(job, f'could not start: {type(e).__name__}: {e}', pending)
                if todo or pending or running or checking:
                    self.sleep(self.poll)
        finally:
            for job in running:                                        # a stop request or an unexpected error: leave no transfer behind
                self.kill_fn(job.proc)
            checker.shutdown(wait=False)
