"""Dependency-free helpers (no FileManager import) so selftest and genome.py work anywhere."""
import os, subprocess

PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # .../YH_SD_Paper (wherever it is checked out)


def config_path(name):
    return os.path.join(PKG_ROOT, 'config', name)


def sh(cmd, check=True, **kw):
    """Run a command; raise with the tail of stderr on failure."""
    r = subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True, **kw)
    if check and r.returncode:
        raise RuntimeError(f'FAILED: {cmd}\n{r.stderr[-2000:]}')
    return r
