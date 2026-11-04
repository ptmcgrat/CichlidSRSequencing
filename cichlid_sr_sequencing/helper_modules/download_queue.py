"""Bounded, retrying, verifying parallel downloads for the alignment pipeline.

Why: AlignmentWorker.downloadReadData used to start one rclone process per file for ALL samples at once, kept only the LAST process of each sample, let a single
FileNotFoundError kill the run (leaving the already started rclone processes behind), and never checked what it downloaded. DownloadQueue replaces that:
  * at most `max_parallel` downloads run at a time (the rest wait in a queue);
  * every file is tracked (a sample with several uBams is no longer judged by its last file only);
  * a failed start ('Cant find file' raised while launching), a non-zero exit or a file that fails the check is retried with exponential backoff;
  * a downloaded file must end with the BGZF EOF marker and pass `gzip -t` (truncated or corrupt BAMs are caught; samtools quickcheck is not used because it rejects
    unmapped BAMs); a bad file is deleted before the retry, because rclone would otherwise skip a same-size corrupt file;
  * files that are already complete are not downloaded again (a rerun only fetches what is missing);
  * on Ctrl-C or any error every running rclone (and its children) is terminated;
  * failures are returned per sample instead of raised."""
import glob, os, random, subprocess, time
from collections import deque
from types import SimpleNamespace

try:
    import psutil
except ImportError:                      # pragma: no cover
    psutil = None

# Gentle defaults for every rclone child (an exported value always wins): one stream per file, so a few files at a time do not become hundreds of connections.
RCLONE_DEFAULTS = {'RCLONE_MULTI_THREAD_STREAMS': '1', 'RCLONE_TRANSFERS': '1', 'RCLONE_CHECKERS': '2'}


def configure_rclone_env(extra=None):
    for k, v in {**RCLONE_DEFAULTS, **(extra or {})}.items():
        os.environ.setdefault(k, v)
    return {k: os.environ[k] for k in RCLONE_DEFAULTS}


def kill_tree(p, grace=5):
    """Terminate a process and everything it started."""
    procs = []
    try:
        if psutil is not None:
            parent = psutil.Process(p.pid)
            procs = parent.children(recursive=True) + [parent]
    except Exception:                    # already gone
        return
    for q in procs:
        try:
            q.terminate()
        except Exception:
            pass
    if psutil is not None:
        _, alive = psutil.wait_procs(procs, timeout=grace)
        for q in alive:
            try:
                q.kill()
            except Exception:
                pass
    else:                                # pragma: no cover
        p.terminate()
    try:
        p.wait(timeout=grace)
    except Exception:
        pass


# A BAM file is BGZF: a series of gzip members that must END with this fixed 28-byte empty block. A truncated download has no (or a cut) marker.
BGZF_EOF = bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')


def bam_check(path, full=True):
    """(ok, reason). Works for unmapped BAMs too: `samtools quickcheck` is NOT used because it demands at least one @SQ reference line in the header and so rejects every uBAM.
    1) exists and is not empty;  2) ends with the BGZF EOF marker (catches truncation, also at a block boundary);
    3) if full: `gzip -t` reads the whole file and verifies every block's CRC (catches corruption). Skipped, with ok=True and a warning reason, if gzip is not installed."""
    if not os.path.exists(path):
        return False, 'file missing after download'
    size = os.path.getsize(path)
    if size == 0:
        return False, 'empty file'
    if size < len(BGZF_EOF):
        return False, 'file too short to be a BAM'
    with open(path, 'rb') as f:
        f.seek(-len(BGZF_EOF), 2)
        if f.read(len(BGZF_EOF)) != BGZF_EOF:
            return False, 'BGZF end-of-file marker missing (truncated download)'
    if not full:
        return True, ''
    try:
        r = subprocess.run(['gzip', '-t', path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        return True, 'gzip not available: file contents not verified (end-of-file marker only)'
    return (r.returncode == 0), ('' if r.returncode == 0 else 'gzip -t failed (corrupt BAM)')


def _log(msg):
    print('  [download] ' + time.strftime('%H:%M:%S') + ' ' + msg, flush=True)


class DownloadQueue:
    def __init__(self, start_fn, max_parallel=6, attempts=4, backoff=30.0, poll=2.0, check_fn=bam_check, log=_log, sleep=time.sleep, clock=time.time, timeout=None):
        """start_fn(job) -> subprocess.Popen-like object (poll(), pid) or raises. job has .key (sample) and .path (local file)."""
        self.start_fn, self.max_parallel, self.attempts, self.backoff, self.poll = start_fn, max(1, int(max_parallel)), max(1, int(attempts)), backoff, poll
        self.check_fn, self.log, self.sleep, self.clock, self.timeout = check_fn, log, sleep, clock, timeout
        self._warned = False

    # ---- one file failed an attempt: clean up, then retry later or give up
    def _fail(self, job, reason, pending, failures):
        for f in [job.path] + glob.glob(job.path + '.*partial'):
            if os.path.exists(f):
                os.remove(f)                                     # a corrupt same-size file would be skipped by rclone on the retry
        if job.attempt < self.attempts:
            delay = self.backoff * 2 ** (job.attempt - 1) * (1 + 0.25 * random.random())
            job.not_before = self.clock() + delay
            pending.append(job)
            self.log(f'{os.path.basename(job.path)}: {reason}; attempt {job.attempt}/{self.attempts}, retrying in {delay:.0f}s')
        else:
            failures.setdefault(job.key, []).append(f'{os.path.basename(job.path)}: {reason} (after {job.attempt} attempts)')
            self.log(f'{os.path.basename(job.path)}: GIVING UP: {reason}')

    def run(self, jobs):
        """jobs: iterable of (key, local_path). Returns {key: [failure messages]} (empty dict = everything downloaded and verified)."""
        pending, running, failures, n_skipped = deque(), [], {}, 0
        for key, path in jobs:
            ok, _ = self.check_fn(path) if os.path.exists(path) else (False, '')
            if ok:
                n_skipped += 1
                continue
            if os.path.exists(path):
                os.remove(path)                                  # incomplete earlier download
            pending.append(SimpleNamespace(key=key, path=path, attempt=0, not_before=0.0, proc=None, started=None))
        if n_skipped:
            self.log(f'{n_skipped} file(s) already downloaded and verified; skipped')
        total = len(pending)
        done = 0
        try:
            while pending or running:
                now = self.clock()
                for job in list(pending):                         # start what is allowed
                    if len(running) >= self.max_parallel:
                        break
                    if job.not_before > now:
                        continue
                    pending.remove(job)
                    job.attempt += 1
                    try:
                        job.proc = self.start_fn(job); job.started = now; running.append(job)
                    except KeyboardInterrupt:
                        raise
                    except Exception as e:                       # e.g. FileNotFoundError('Cant find file for download') raised while launching
                        self._fail(job, f'could not start: {type(e).__name__}: {e}', pending, failures)
                for job in list(running):                         # collect what finished
                    rc = job.proc.poll()
                    if rc is None:
                        if self.timeout and now - job.started > self.timeout:
                            kill_tree(job.proc); running.remove(job); self._fail(job, f'timed out after {self.timeout:.0f}s', pending, failures)
                        continue
                    running.remove(job)
                    if rc != 0:
                        self._fail(job, f'download process exited with code {rc}', pending, failures)
                        continue
                    ok, why = self.check_fn(job.path)
                    if ok and why and not self._warned:
                        self._warned = True; self.log('WARNING: ' + why)
                    if ok:
                        done += 1
                        self.log(f'{os.path.basename(job.path)} ok ({done}/{total}); {len(running)} running, {len(pending)} waiting')
                    else:
                        self._fail(job, why, pending, failures)
                if pending or running:
                    self.sleep(self.poll)
        finally:
            for job in running:                                  # Ctrl-C or an unexpected error: leave no rclone behind
                kill_tree(job.proc)
        return failures
