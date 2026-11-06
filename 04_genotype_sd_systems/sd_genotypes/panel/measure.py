"""Run the enabled modules on one BAM and write <sample>.panel.json. Needs only samtools on PATH; no Project, no Dropbox, no FileManager."""
import datetime, importlib, json, os
import yaml
from cichlid_common.common_light import config_path, sh, repo_version

MODULES = {'gsdf': ('sd_genotypes.gsdf.measure', 'gsdf.yaml')}              # LG5 inversion and B chromosome are added here when their definitions are fixed
SCHEMA = 1


def module_cfg(name):
    return yaml.safe_load(open(config_path(MODULES[name][1])))


def package_version():
    return repo_version()


def bam_index_present(bam):
    return any(os.path.exists(p) for p in (bam + '.bai', bam[:-4] + '.bai', bam + '.csi'))


def measure(bam, sample, modules=('gsdf',), run=None, cfgs=None, now=None):
    """Returns the JSON-able dict. run(cmd) defaults to common_light.sh. cfgs: {module: config} overrides (tests)."""
    run = run or sh
    if not os.path.exists(bam):
        raise FileNotFoundError(bam)
    if not bam_index_present(bam):
        raise FileNotFoundError(f'{bam} has no index (.bai / .csi next to it): the panel reads regions and needs one')
    unknown = [m for m in modules if m not in MODULES]
    if unknown:
        raise ValueError(f'unknown panel modules {unknown}; available: {sorted(MODULES)}')
    out = dict(schema=SCHEMA, sample=sample, bam=os.path.abspath(bam), bam_size=os.path.getsize(bam), panel_version=package_version(),
               measured=(now or datetime.datetime.now()).isoformat(timespec='seconds'), samtools=run(['samtools', '--version']).stdout.splitlines()[0], modules={})
    for m in modules:
        mod = importlib.import_module(MODULES[m][0])
        cfg = (cfgs or {}).get(m) or module_cfg(m)
        out['modules'][m] = mod.measure_bam(bam, cfg, run)
    return out


def write_json(path, obj):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + '.tmp'
    json.dump(obj, open(tmp, 'w'), indent=1)
    os.replace(tmp, path)


def read_panel_json(path):
    d = json.load(open(path))
    if d.get('schema') != SCHEMA:
        raise ValueError(f"{path}: schema {d.get('schema')} (this code reads schema {SCHEMA})")
    return d


def collect(path, module='gsdf'):
    """Raw rows of one module from a folder of *.panel.json (or a single file): [{SampleID, bins, bins_q, junction_depth, junction_depth_q, definition, ...}]."""
    files = [path] if os.path.isfile(path) else sorted(os.path.join(path, f) for f in os.listdir(path) if f.endswith('.panel.json'))
    rows, skipped = [], []
    for f in files:
        d = read_panel_json(f)
        if module not in d['modules']:
            skipped.append(os.path.basename(f)); continue
        rows.append(dict(SampleID=d['sample'], bam_size=d['bam_size'], measured=d['measured'], panel_version=d['panel_version'], **d['modules'][module]))
    return rows, skipped
