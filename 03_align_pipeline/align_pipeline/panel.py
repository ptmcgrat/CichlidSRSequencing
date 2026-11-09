"""The panel stage: module 04's measurements (gsdf, te, lg5inv, bchrom) made in the same pass as the alignment, while the BAM is on local disk (decision 40).

The contract with module 04 (docs/decisions.md 27 and 28, SPEC.md module 04) and what this module does with it:
  * preflight, once per run, before anything is aligned:  python -m sd_genotypes panel check --modules all
    It loads every module's config and definition, prints the definition keys and exits 1 if any is missing. A run does not start when it fails (`preflight`).
  * per sample, after duplicate marking, on the coordinate-sorted, duplicate-marked BAM with its index next to it:
        python -m sd_genotypes panel measure --bam <sample>.all.bam --sample <ID> --out-json <ID>.panel.json --modules all
    one atomic JSON, or a non-zero exit and no JSON. No panel logic is copied here: 04 owns the definitions, the pipeline only calls it as a subprocess.
  * a failure never fails the sample (the BAM is the primary product): the pipeline records PanelStatus and a QC warning, and module 04's own commands (`panel measure`, `retro`) fill the JSON in later.
Since 2026-11-09 (decision 45) the panel JSON (schema 2) records the BAM it was measured on: `bam` (the file name) and `bam_size_bytes`. After every measurement the pipeline checks that the JSON it just
got names THIS BAM (`identity_status`, the same rule as `python -m sd_genotypes.panel.identity`): a JSON that names another BAM is a failed panel, one without an identity (schema 1, an old sd_genotypes)
is accepted with a QC warning. A name and a size are not a checksum, but a realigned sample keeps its file name and not its size, which is what the rerun case needs.
For the pilot the definition keys are frozen (decision 28): with `expect` the preflight also stops the run when the keys `panel check` prints are not the expected ones."""
import json
import os
import subprocess
import sys
import time


class PanelConfig:
    def __init__(self, python=None, modules='all', pythonpath=None, timeout_s=1800, expect=None):
        """expect: definition keys that `panel check` must print (the pilot's frozen keys); none = no check. python: the interpreter that can `import sd_genotypes` (default: the one running the pipeline). pythonpath: where sd_genotypes is checked out, when it is not installed
        (for example <repo>/04_genotype_sd_systems). timeout_s: a measurement that takes longer is killed (it normally takes seconds)."""
        self.python, self.modules, self.pythonpath, self.timeout_s = python or sys.executable, modules, pythonpath, timeout_s
        self.expect = list(expect or [])

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
    text = (r.stdout or '').strip()
    if r.returncode != 0:
        return False, text
    missing = [k for k in cfg.expect if k not in text]
    if missing:
        return False, f'the definition keys are not the expected ones (not printed by `panel check`: {missing}); a definition or a measurement setting changed since the keys were frozen: find out why before measuring.\n' + text
    return True, text


def identity_status(path_or_record, bam):
    """('same' | 'different' | 'unknown' | 'bad', reason): does this panel JSON belong to this BAM? Same rule as sd_genotypes.panel.identity: the file name AND the size must match. 'unknown' when no size is
    recorded (schema 1 files of streamed samples); schema 1 files that have `bam` and `bam_size` are read as before. 'bad' when the JSON cannot be read."""
    try:
        if isinstance(path_or_record, dict):
            rec = path_or_record
        else:
            with open(path_or_record) as f:
                rec = json.load(f)
    except (OSError, ValueError) as e:
        return 'bad', f'the panel JSON cannot be read: {type(e).__name__}: {e}'
    name = rec.get('bam')
    size = rec.get('bam_size_bytes', rec.get('bam_size'))
    if not name or size in (None, '', 'unknown'):
        return 'unknown', 'the panel JSON records no BAM identity (schema 1: an old sd_genotypes)'
    if os.path.basename(str(name)) != os.path.basename(bam):
        return 'different', f'it was measured on {os.path.basename(str(name))}, not {os.path.basename(bam)}'
    try:
        recorded = int(size)
    except (TypeError, ValueError):
        return 'unknown', f'the recorded size {size!r} is not a number'
    actual = os.path.getsize(bam)
    if recorded != actual:
        return 'different', f'it was measured on a BAM of {recorded} bytes, this one has {actual}'
    return 'same', ''
