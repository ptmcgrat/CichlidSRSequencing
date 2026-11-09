"""Small process and file helpers that the pipeline needs, with no dependency on the old download queue and none on psutil (used when it happens to be installed)."""
import os
import subprocess

try:
    import psutil
except ImportError:                                                   # optional: without it a transfer's children are stopped through its process group or by terminate()
    psutil = None

BGZF_EOF = bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')      # the 28-byte end-of-file block of a complete BAM
RCLONE_DEFAULTS = {'RCLONE_MULTI_THREAD_STREAMS': '1', 'RCLONE_TRANSFERS': '1', 'RCLONE_CHECKERS': '2'}


def configure_rclone_env(extra=None):
    for k, v in {**RCLONE_DEFAULTS, **(extra or {})}.items():
        os.environ.setdefault(k, v)
    return {k: os.environ[k] for k in RCLONE_DEFAULTS}


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


def kill_tree(p, grace=5):
    """Terminate a process and everything it started (psutil when installed, otherwise the process itself)."""
    procs = []
    try:
        if psutil is not None:
            parent = psutil.Process(p.pid)
            procs = parent.children(recursive=True) + [parent]
    except Exception:                                                  # already gone
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
    else:
        p.terminate()
    try:
        p.wait(timeout=grace)
    except Exception:
        try:
            p.kill()
        except Exception:
            pass
