"""The panel stage: module 04's measurements (gsdf, te, lg5inv, bchrom) made in the same pass as the alignment, while the BAM is on local disk (decision 40).

The contract with module 04 (docs/decisions.md 27 and 28, SPEC.md module 04) and what this module does with it:
  * preflight, once per run, before anything is aligned:  python -m sd_genotypes panel check --modules all
    It loads every module's config and definition, prints the definition keys and exits 1 if any is missing. A run does not start when it fails (`preflight`).
  * per sample, after duplicate marking, on the coordinate-sorted, duplicate-marked BAM with its index next to it:
        python -m sd_genotypes panel measure --bam <sample>.all.bam --sample <ID> --out-json <ID>.panel.json --modules all
    one atomic JSON, or a non-zero exit and no JSON. No panel logic is copied here: 04 owns the definitions, the pipeline only calls it as a subprocess.
  * a failure never fails the sample (the BAM is the primary product): the pipeline records PanelStatus and a QC warning and `genotype-only` fills the JSON in later.
The panel JSON does not say which BAM it was measured on. Until module 04 records that (requested, SPEC.md module 03 item 6), the pipeline writes a sidecar next to it (<sample>.panel.meta.json) with the
BAM's name and size and the fingerprint of the runs it was made from, so that a re-alignment makes the old measurement detectable as stale (`sidecar_is_current`)."""
import json
import os
import subprocess
import sys
import time


class PanelConfig:
    def __init__(self, python=None, modules='all', pythonpath=None, timeout_s=1800):
        """python: the interpreter that can `import sd_genotypes` (default: the one running the pipeline). pythonpath: where sd_genotypes is checked out, when it is not installed
        (for example <repo>/04_genotype_sd_systems). timeout_s: a measurement that takes longer is killed (it normally takes seconds)."""
        self.python, self.modules, self.pythonpath, self.timeout_s = python or sys.executable, modules, pythonpath, timeout_s

    def prefix(self):
        if not self.pythonpath:
            return []
        old = os.environ.get('PYTHONPATH')
        return ['env', 'PYTHONPATH=' + self.pythonpath + ((os.pathsep + old) if old else '')]


def check_command(cfg):
    return cfg.prefix() + [cfg.python, '-m', 'sd_genotypes', 'panel', 'check', '--modules', cfg.modules]


def measure_command(cfg, bam, sample, out_json):
    return cfg.prefix() + [cfg.python, '-m', 'sd_genotypes', 'panel', 'measure', '--bam', bam, '--sample', sample, '--out-json', out_json, '--modules', cfg.modules]


def preflight(cfg, run=subprocess.run, timeout=300):
    """(ok, text): `panel check` run once. ok is False when it exits non-zero, cannot start, or does not answer within `timeout` seconds; text is what it printed (the definition keys when ok)."""
    cmd = check_command(cfg)
    try:
        r = run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=timeout)
    except FileNotFoundError as e:
        return False, f'could not start {cmd[0]}: {e}'
    except subprocess.TimeoutExpired:
        return False, f'`panel check` did not answer within {timeout} s'
    return r.returncode == 0, (r.stdout or '').strip()


def write_sidecar(path, bam, sample, fingerprint, cfg, now=None):
    """Which BAM this panel JSON belongs to; written atomically, only after a successful measurement."""
    rec = dict(schema=1, sample=sample, bam=os.path.basename(bam), bam_size=os.path.getsize(bam), input_fingerprint=fingerprint or '',
               measured_on=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(now if now is not None else time.time())), modules=cfg.modules)
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(rec, f, indent=1)
    os.replace(tmp, path)
    return rec


def sidecar_is_current(path, bam, fingerprint):
    """True when the sidecar exists and was written for THIS BAM (same name and size) made from THESE runs (same fingerprint); False otherwise, including an unreadable sidecar."""
    try:
        with open(path) as f:
            rec = json.load(f)
        return rec.get('bam') == os.path.basename(bam) and rec.get('bam_size') == os.path.getsize(bam) and rec.get('input_fingerprint', '') == (fingerprint or '')
    except (OSError, ValueError):
        return False
