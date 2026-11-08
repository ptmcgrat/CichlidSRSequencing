"""Run the enabled modules on one BAM and write <sample>.panel.json. Needs only samtools on PATH; no Project, no Dropbox, no FileManager."""
import datetime, importlib, json, os
import yaml
from cichlid_common.common_light import config_path, sh, repo_version

MODULES = {'gsdf': ('sd_genotypes.gsdf.measure', 'gsdf.yaml'),
           'lg5inv': ('sd_genotypes.lg5inv.measure', 'lg5panel.yaml'),      # LG5 inversion: needs the definition files in sd_genotypes/lg5inv/data/
           'bchrom': ('sd_genotypes.bchrom.measure', 'bchrompanel.yaml'),  # B chromosome (disk only): needs the definition files in sd_genotypes/bchrom/data/
           'te': ('sd_genotypes.te.measure', 'tepanel.yaml')}              # transposon insertions TE1, TE2 (LG10) and chr7ins (LG7); needs sd_genotypes/te/data/te_definition.json. Modules with a regions(cfg) function can be measured retrospectively (sd_genotypes.retro)
SCHEMA = 1


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


def measure(bam, sample, modules=('gsdf',), run=None, cfgs=None, now=None):
    """Returns the JSON-able dict. run(cmd) defaults to common_light.sh. cfgs: {module: config} overrides (tests)."""
    run = run or sh
    if not os.path.exists(bam):
        raise FileNotFoundError(bam)
    if not bam_index_present(bam):
        raise FileNotFoundError(f'{bam} has no index (.bai / .csi next to it): the panel reads regions and needs one')
    modules = resolve_modules(modules)
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
