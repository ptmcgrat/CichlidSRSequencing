"""Run records + Dropbox as the system of record.

Every analysis step runs inside `with proj.run(analysis, step, params) as run:` and writes into `run.dir`:
    <localMasterDir>/YH_SD_Paper/<genome>/<analysis>/<step>/<run_id>/        (run_id = YYYYmmdd-HHMMSS_<git short hash>[-dirty])
On success the runner writes manifest.json (command, git commit, config, tool versions, input/output files with sizes+md5),
uploads the run folder to the same path under cloudMasterDir (Dropbox via rclone), VERIFIES it with `rclone check`, and writes a
LATEST pointer. Runs are never overwritten. Large directories of many small files can be tarred before upload.
Steps that keep a working directory (the BC genotyping pipeline) use `publish()` to upload a snapshot of chosen outputs the same way."""
import datetime, getpass, hashlib, json, os, shutil, socket, sys, tarfile
import yaml
from .common_light import sh

TOOLS = {'samtools': 'samtools --version', 'bcftools': 'bcftools --version', 'gatk': 'gatk --version',
         'bedtools': 'bedtools --version', 'mosdepth': 'mosdepth --version', 'python': 'python --version'}


def git_state(path):
    """(commit, dirty) for the run record.
    1. DEPLOYED_FROM.txt (written by scripts/deploy_to_utaka.sh on the laptop) is authoritative: commit + dirty flag at deploy time.
    2. Otherwise, if the package files are tracked by the git repo that contains them, use that repo's HEAD and dirty state.
    3. Otherwise (e.g. unzipped into a clone that does not track it) report the package version from VERSION.txt, never an unrelated HEAD."""
    from .common_light import PKG_ROOT
    dep = os.path.join(PKG_ROOT, 'DEPLOYED_FROM.txt')
    if os.path.exists(dep):
        kv = dict(l.strip().split(': ', 1) for l in open(dep) if ': ' in l)
        return kv.get('commit', 'nogit'), (kv.get('dirty', '0') == '1')
    try:
        tracked = sh(['git', '-C', PKG_ROOT, 'ls-files', '--', 'config/paper.yaml']).stdout.strip()
        if tracked:
            commit = sh(['git', '-C', PKG_ROOT, 'rev-parse', '--short', 'HEAD']).stdout.strip()
            dirty = bool(sh(['git', '-C', PKG_ROOT, 'status', '--porcelain', '--', PKG_ROOT]).stdout.strip())
            return commit, dirty
    except Exception:
        pass
    ver = os.path.join(PKG_ROOT, 'VERSION.txt')
    if os.path.exists(ver):
        for l in open(ver):
            if l.startswith('version'):
                return 'pkg-' + l.split(':', 1)[1].strip(), None
    return 'nogit', None


def tool_versions():
    out = {}
    for k, cmd in TOOLS.items():
        r = sh(f'{cmd} 2>&1 | head -2 | tr "\\n" " "', check=False)
        out[k] = r.stdout.strip() or 'n/a'
    return out


def file_record(path, md5_limit=300e6):
    size = os.path.getsize(path)
    rec = {'path': path, 'size': size}
    if size <= md5_limit:
        h = hashlib.md5()
        with open(path, 'rb') as f:
            for chunk in iter(lambda: f.read(1 << 20), b''):
                h.update(chunk)
        rec['md5'] = h.hexdigest()
    return rec


def to_cloud(proj, local):
    return local.replace(proj.fm.localMasterDir, proj.fm.cloudMasterDir)


def sync(proj, local, verify=True):
    """rclone copy local file/dir to the matching Dropbox path, then rclone check. Raises on any difference."""
    cloud = to_cloud(proj, local)
    if os.path.isdir(local):
        sh(['rclone', 'copy', local, cloud])
    else:
        sh(['rclone', 'copy', local, os.path.dirname(cloud) + '/'])
    if verify:
        r = sh(['rclone', 'check', local, cloud if os.path.isdir(local) else os.path.dirname(cloud) + '/', '--one-way'], check=False)
        if r.returncode:
            raise RuntimeError(f'UPLOAD VERIFICATION FAILED for {local}\n{r.stderr[-800:]}')
    return cloud


def tar_dir(path):
    """Tar a directory of many small files next to itself; returns the tar path."""
    t = path.rstrip('/') + '.tar'
    with tarfile.open(t, 'w') as tf:
        tf.add(path, arcname=os.path.basename(path.rstrip('/')))
    return t


class Run:
    def __init__(self, proj, analysis, step, params=None, tar=()):
        self.proj, self.analysis, self.step, self.params, self.tar = proj, analysis, step, params or {}, tuple(tar)
        self.commit, self.dirty = git_state(None)
        stamp = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        self.run_id = f'{stamp}_{self.commit}' + ('-dirty' if self.dirty else '')
        self.base = f'{proj.root}{analysis}/{step}/'
        self.dir = f'{self.base}{self.run_id}/'
        self.inputs = []

    def add_input(self, path):
        if os.path.exists(path):
            self.inputs.append(file_record(path))
        else:
            self.inputs.append({'path': path, 'missing': True})

    def __enter__(self):
        os.makedirs(self.dir, exist_ok=True)
        self.t0 = datetime.datetime.now()
        return self

    def __exit__(self, et, ev, tb):
        if et is not None:
            open(self.dir + 'FAILED.txt', 'w').write(f'{et.__name__}: {ev}\n')
            return False
        self.finalize()
        return False

    def finalize(self):
        for d in self.tar:                                   # tar directories of many small files
            p = self.dir + d
            if os.path.isdir(p):
                tar_dir(p); shutil.rmtree(p)
        outs = []
        for root, _, files in os.walk(self.dir):
            for fn in files:
                if fn != 'manifest.json':
                    outs.append(file_record(os.path.join(root, fn)))
        man = dict(run_id=self.run_id, analysis=self.analysis, step=self.step, genome=self.proj.genome,
                   command=' '.join(sys.argv), user=getpass.getuser(), host=socket.gethostname(),
                   started=self.t0.isoformat(), finished=datetime.datetime.now().isoformat(),
                   git_commit=self.commit, git_dirty=self.dirty, params=self.params,
                   config=self.proj.cfg, tools=tool_versions(), inputs=self.inputs, outputs=outs)
        json.dump(man, open(self.dir + 'manifest.json', 'w'), indent=1, default=str)
        if self.proj.cfg.get('runner', {}).get('upload', True):
            cloud = sync(self.proj, self.dir)
            open(self.base + 'LATEST', 'w').write(self.run_id + '\n')
            sync(self.proj, self.base + 'LATEST', verify=False)
            print(f'[runner] {self.analysis}/{self.step} {self.run_id}: uploaded and verified -> {cloud}')
        else:
            open(self.base + 'LATEST', 'w').write(self.run_id + '\n')
            print(f'[runner] {self.analysis}/{self.step} {self.run_id}: upload disabled (local only)')


def latest_dir(proj, analysis, step):
    """Most recent run folder; pulls it from Dropbox if it is not on local disk."""
    base = f'{proj.root}{analysis}/{step}/'
    lp = base + 'LATEST'
    if not os.path.exists(lp):
        os.makedirs(base, exist_ok=True)
        sh(['rclone', 'copy', to_cloud(proj, lp), base], check=False)
    if not os.path.exists(lp):
        raise FileNotFoundError(f'no LATEST run for {analysis}/{step} locally or on Dropbox')
    rid = open(lp).read().strip()
    d = f'{base}{rid}/'
    if not os.path.isdir(d):
        sh(['rclone', 'copy', to_cloud(proj, d), d])
    return d


def publish(proj, analysis, step, paths, params=None, tar=()):
    """Snapshot selected outputs of a working-directory step as an immutable, verified Dropbox run (manifest included)."""
    with Run(proj, analysis, step, params, tar) as run:
        for p in paths:
            if not os.path.exists(p):
                print(f'[runner] skip missing {p}'); continue
            dst = run.dir + os.path.basename(p.rstrip('/'))
            if os.path.isdir(p):
                shutil.copytree(p, dst)                       # copy so the working directory stays usable
            else:
                shutil.copy2(p, dst)
    return run.dir
