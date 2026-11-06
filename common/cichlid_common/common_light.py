"""Dependency-free helpers (no FileManager import) so selftest and genome.py work anywhere."""
import os, subprocess
import glob, sys

# ---- where things are. REPO_ROOT is the folder that holds SPEC.md, common/ and the numbered module folders.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PKG_ROOT = REPO_ROOT                                                    # old name, kept so that code written for the single package keeps working


def module_dirs():
    """Folders that must be on sys.path (PYTHONPATH of subprocesses): common/ and every numbered module folder."""
    out = [os.path.join(REPO_ROOT, 'common')]
    for d in sorted(glob.glob(os.path.join(REPO_ROOT, '[0-9][0-9]_*'))):
        if os.path.isdir(d):
            out.append(d)
    return out


def subprocess_env(extra=None):
    """The environment for a subprocess that runs `python -m <package>` of any module."""
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join(module_dirs() + ([env['PYTHONPATH']] if env.get('PYTHONPATH') else []))
    env.update(extra or {})
    return env


def config_path(name):
    """config/<name> of whichever module owns it (names are unique across modules): common/config, then <module>/config."""
    cands = [os.path.join(REPO_ROOT, 'common', 'config', name)] + sorted(glob.glob(os.path.join(REPO_ROOT, '[0-9][0-9]_*', 'config', name))) + [os.path.join(REPO_ROOT, 'config', name)]
    for c in cands:
        if os.path.exists(c):
            return c
    return cands[0]


def repo_version():
    try:
        for l in open(os.path.join(REPO_ROOT, 'VERSION.txt')):
            if l.startswith('version'):
                return l.split(':', 1)[1].strip()
    except OSError:
        pass
    return 'unknown'


def ensure_lab_modules_on_path():
    """The lab's FileManager and alignment helpers live in a folder called helper_modules/ (today <repo>/cichlid_sr_sequencing/helper_modules; after the move, inside the numbered folders).
    Put its parent on sys.path so that `import helper_modules.x` works, whichever layout is checked out."""
    seen = set(sys.path)
    for pat in ('helper_modules', '*/helper_modules', '*/*/helper_modules'):
        for hm in glob.glob(os.path.join(REPO_ROOT, pat)):
            parent = os.path.dirname(hm)
            if os.path.isdir(hm) and parent not in seen:
                sys.path.insert(0, parent); seen.add(parent)


def sh(cmd, check=True, **kw):
    """Run a command; raise with the tail of stderr on failure."""
    r = subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True, **kw)
    if check and r.returncode:
        raise RuntimeError(f'FAILED: {cmd}\n{r.stderr[-2000:]}')
    return r
