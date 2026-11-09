"""Run the enabled modules on one BAM and write <sample>.panel.json. Needs only samtools on PATH; no Project, no Dropbox, no FileManager."""
import datetime, importlib, json, os
import yaml
from cichlid_common.common_light import config_path, sh, repo_version
from . import identity as ID

MODULES = {'gsdf': ('sd_genotypes.gsdf.measure', 'gsdf.yaml'),
           'lg5inv': ('sd_genotypes.lg5inv.measure', 'lg5panel.yaml'),      # LG5 inversion: needs the definition files in sd_genotypes/lg5inv/data/
           'bchrom': ('sd_genotypes.bchrom.measure', 'bchrompanel.yaml'),  # B chromosome (disk only): needs the definition files in sd_genotypes/bchrom/data/
           'te': ('sd_genotypes.te.measure', 'tepanel.yaml')}              # transposon insertions TE1, TE2 (LG10) and chr7ins (LG7); needs sd_genotypes/te/data/te_definition.json. Modules with a regions(cfg) function can be measured retrospectively (sd_genotypes.retro)
SCHEMA = 2                                       # 2: `bam` is the file name and `bam_size_bytes` its size (the identity of the BAM measured); 1 had an absolute path and `bam_size`
SCHEMAS_READ = (1, 2)


def module_cfg(name):
    return yaml.safe_load(open(config_path(MODULES[name][1])))


def package_version():
    return repo_version()


def bam_index_present(bam):
    return any(os.path.exists(p) for p in (bam + '.bai', bam[:-4] + '.bai', bam + '.csi'))


def resolve_modules(spec):
    """'all' (every registered module), a comma-separated string, or a list/tuple of names -> list of module names (checked against the registry)."""
    names = list(MODULES) if spec in ('all', ['all'], ('all',)) else (spec.split(',') if isinstance(spec, str) else list(spec))
    names = [n.strip() for n in names if n.strip()]
    unknown = [n for n in names if n not in MODULES]
    if unknown:
        raise ValueError(f'unknown panel modules {unknown}; available: {sorted(MODULES)}')
    return names


def check(modules='all', cfgs=None):
    """Preflight for a pipeline: can every module measure? Loads its config and definition files and returns [{module, status ('OK' | 'MISSING' | 'ERROR'), key (the definition key), message}].
    Run it once before a long job: a missing definition then stops the run before anything is aligned."""
    out = []
    for m in resolve_modules(modules):
        try:
            mod = importlib.import_module(MODULES[m][0]); cfg = (cfgs or {}).get(m) or module_cfg(m)
            key = mod.definition_key(cfg) if hasattr(mod, 'definition_key') else ''
            out.append(dict(module=m, status='OK', key=key, message=''))
        except FileNotFoundError as e:
            out.append(dict(module=m, status='MISSING', key='', message=str(e)[:200]))
        except Exception as e:
            out.append(dict(module=m, status='ERROR', key='', message=f'{type(e).__name__}: {str(e)[:180]}'))
    return out


def measure(bam, sample, modules=('gsdf',), run=None, cfgs=None, now=None, identity=None):
    """Returns the JSON-able dict. run(cmd) defaults to common_light.sh. cfgs: {module: config} overrides (tests). `identity` ({'bam': name, 'bam_size_bytes': n}) replaces the identity taken from the file
    measured: used when the file is a temporary regional BAM streamed from a larger one (sd_genotypes.retro)."""
    run = run or sh
    if not os.path.exists(bam):
        raise FileNotFoundError(bam)
    if not bam_index_present(bam):
        raise FileNotFoundError(f'{bam} has no index (.bai / .csi next to it): the panel reads regions and needs one')
    modules = resolve_modules(modules)
    ident = identity or ID.bam_identity(bam)
    out = dict(schema=SCHEMA, sample=sample, bam=ident['bam'], bam_size_bytes=ident['bam_size_bytes'], panel_version=package_version(),
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
    if d.get('schema') not in SCHEMAS_READ:
        raise ValueError(f"{path}: schema {d.get('schema')} (this code reads schemas {SCHEMAS_READ})")
    return d


def collect(path, module='gsdf'):
    """Raw rows of one module from a folder of *.panel.json (or a single file): [{SampleID, bam, bam_size_bytes, measured, panel_version, <the module's raw measurements>, definition}]. Schema 1 files are read too."""
    files = [path] if os.path.isfile(path) else sorted(os.path.join(path, f) for f in os.listdir(path) if f.endswith('.panel.json'))
    rows, skipped = [], []
    for f in files:
        d = read_panel_json(f)
        if module not in d['modules']:
            skipped.append(os.path.basename(f)); continue
        name, size = ID.normalize(d)
        rows.append(dict(SampleID=d['sample'], bam=name, bam_size_bytes=size, measured=d['measured'], panel_version=d['panel_version'], **d['modules'][module]))
    return rows, skipped
