"""Effective coverage (module 05): the depth the variant caller actually sees in a BAM, measured per sample.

Why: the AlignmentDatabase's `Coverage` is raw yield (every read x its length / the genome). The caller then drops duplicates, reads with MAPQ below 13, low-quality bases and the second copy of
every base in overlapping mates, and what is left differs a lot between libraries (on 2026-11-11: 41% of the table's coverage for one saulosi library, 81 to 85% for two others).
Sexloci needs depth >= 5 at a genotype, so the number that matters is the one the caller sees.

What is measured, per sample, from the sample's `all.bam` (downloaded, measured, deleted):
  * `EffectiveCoverage`        the mean FORMAT/DP of `bcftools mpileup` run EXACTLY as `cohort build` runs it (the same function builds the command: MAPQ >= 13, base quality >= 20, BAQ and
                               mate-overlap handling as bcftools does them) over one 1 Mb window per chromosome (22 windows, about 2% of the genome), at the positions the caller emits
                               (a position with no read passing the filters is not emitted: this is the depth at COVERED sites, which is what a genotype is made from)
  * `WindowRawDepth`           the mean INFO/DP over the same positions: the reads that passed the read-level filters (flags, duplicates, MAPQ), before the base-level ones
  * `BaseSurvivalFraction`     EffectiveCoverage / WindowRawDepth: what the base-level filters (base quality, BAQ, mate overlap) keep
  * `WindowCoveredFraction`    the share of the window positions the caller emits
  * `MapqPassFraction`         over the WHOLE BAM: reads that are mapped, primary, not duplicates and have MAPQ >= 13, divided by those without the MAPQ condition (`samtools view -c -F 3844 [-q 13]`)
  * `DuplicateFraction`        over the whole BAM: primary reads flagged as duplicates / primary reads (`samtools flagstat`)
  * `PrimaryReads`, `TableCoverage` (the AlignmentDatabase's own number, for comparison), `BamBytes`, `Seconds`, `Machine`, `Tools`.
A row carries a `DefinitionKey` (a hash of every parameter that changes the numbers); a row measured under another key is measured again.

Rows are appended to ONE tab-separated file as each sample finishes (nothing is lost if a run stops), a failed sample is recorded with its reason and does not stop the others, and a rerun
skips what is already measured with the same definition. Nothing here writes to a Google sheet: that is a separate, later step that reads this file.
"""
import csv
import datetime
import hashlib
import json
import math
import os
import shlex
import shutil
import socket
import statistics
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import yaml

from cichlid_common.common_light import config_path
from variant_calling import build as B

DEFAULTS = dict(window_bp=1_000_000, position_fraction=0.4, min_depth=5, jobs=4, threads=8, default_bam_bytes=15e9, min_free_factor=1.2, definition_version=1,
                output_subdir='Outputs/EffectiveCoverage/', work_subdir='Cohorts/effcov_staged/')
VIEW_EXCLUDE = 3844                      # samtools view -F: unmapped, secondary, QC fail, duplicate, supplementary
COLUMNS = ['SampleID', 'GenomeVersion', 'Status', 'Error', 'DefinitionKey', 'MeasuredOn', 'Machine', 'Tools', 'PrimaryReads', 'DuplicateFraction', 'MapqPassFraction', 'WindowRawDepth',
           'EffectiveCoverage', 'BaseSurvivalFraction', 'WindowCoveredFraction', 'TableCoverage', 'BamBytes', 'Seconds']
TOOLS = ('samtools', 'bcftools', 'rclone', 'awk')


class MeasureError(RuntimeError):
    pass


def load_config(path=None):
    cfg = dict(DEFAULTS)
    p = path or config_path('effcov.yaml')
    if os.path.exists(p):
        cfg.update(yaml.safe_load(open(p)) or {})
    for k in ('window_bp', 'position_fraction', 'min_depth', 'jobs', 'threads', 'default_bam_bytes', 'min_free_factor', 'definition_version'):
        try:
            cfg[k] = float(cfg[k]) if k in ('position_fraction', 'min_depth', 'default_bam_bytes', 'min_free_factor') else int(float(cfg[k]))   # YAML reads 15.0e9 as a string
        except (TypeError, ValueError):
            raise ValueError(f'effcov config: {k} must be a number, not {cfg[k]!r}')
    return cfg


# ------------------------------------------------------------------ the definition
def windows_from_fai(fai, window_bp, fraction, min_bp):
    """[(contig, start, end)], 1-based inclusive: one window of window_bp on every chromosome (the NC_ contigs of at least min_bp, as the build chooses them), starting `fraction` of the way along it."""
    lengths = dict(fai)
    out = []
    for contig, _, _ in B.choose_contigs(fai, 'chromosomes', min_bp):
        n = lengths[contig]
        end = min(n, int(n * fraction) + window_bp)
        out.append((contig, max(1, end - window_bp + 1), end))
    return out


def regions_arg(windows):
    return ','.join(f'{c}:{s}-{e}' for c, s, e in windows)


def definition_key(cfg, ecfg, windows):
    """A hash of everything that changes the numbers: the caller's own parameters (from cohort_build.yaml), the read filter of the whole-BAM counts, the windows, the definition version."""
    b = cfg['bcftools']
    d = dict(version=ecfg['definition_version'], mq=b['min_mapping_quality'], bq=b['min_base_quality'], annotate=b.get('mpileup_annotate'), indel=[b.get('indel_min_reads'), b.get('indel_min_fraction')],
             view_exclude=VIEW_EXCLUDE, windows=[list(w) for w in windows])
    return hashlib.md5(json.dumps(d, sort_keys=True).encode()).hexdigest()[:12]


def window_script(env, bamlist, windows, cfg):
    """The pileup exactly as the build issues it (B.mpileup_cmd), over the windows, reduced to `n_positions<TAB>sum INFO/DP<TAB>sum FORMAT/DP` by awk (22 million positions are not parsed in python)."""
    cmd = B.mpileup_cmd(env.genome_fa, bamlist, regions_arg(windows), cfg)
    return ('set -o pipefail; ' + shlex.join(cmd) + ' | ' + shlex.join(['bcftools', 'query', '-f', '%INFO/DP\t[%DP]\n']) +
            " | awk '{r+=$1; h+=$2; n++} END{printf \"%d\\t%.0f\\t%.0f\\n\", n, r, h}'")


# ------------------------------------------------------------------ parsing and arithmetic (pure)
def parse_flagstat(text):
    """{'primary': n, 'primary_duplicates': n} from `samtools flagstat` (the QC-passed count of each line)."""
    out = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 4 and parts[1] == '+' and parts[3:] == ['primary']:
            out['primary'] = int(parts[0])
        elif len(parts) >= 5 and parts[1] == '+' and parts[3:5] == ['primary', 'duplicates']:
            out['primary_duplicates'] = int(parts[0])
    missing = [k for k in ('primary', 'primary_duplicates') if k not in out]
    if missing:
        raise MeasureError(f'samtools flagstat has no line for {", ".join(missing)}: this samtools is too old (the `primary` lines were added in 1.10)')
    return out


def parse_window_line(text):
    """(n_positions, sum_info_dp, sum_format_dp) from the awk line."""
    try:
        n, r, h = text.strip().splitlines()[-1].split('\t')
        return int(n), int(r), int(h)
    except (ValueError, IndexError):
        raise MeasureError(f'the pileup summary is not three numbers: {text.strip()[-200:]!r}')


def derive(primary, duplicates, mapped_nd, mapq_pass, n, raw_sum, hq_sum, window_bp_total):
    """The measured columns from the counts."""
    if primary <= 0 or mapped_nd <= 0:
        raise MeasureError(f'no usable reads in the BAM (primary {primary}, mapped and not duplicate {mapped_nd})')
    if n <= 0 or raw_sum <= 0:
        raise MeasureError('no position of the windows is covered after the caller\'s filters (the BAM does not match the reference, or the windows hold no reads)')
    return dict(PrimaryReads=primary, DuplicateFraction=duplicates / primary, MapqPassFraction=mapq_pass / mapped_nd, WindowRawDepth=raw_sum / n, EffectiveCoverage=hq_sum / n,
                BaseSurvivalFraction=hq_sum / raw_sum, WindowCoveredFraction=n / window_bp_total)


def _fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ''
    if isinstance(v, float):
        return f'{v:.6g}'
    return str(v).replace('\t', ' ').replace('\n', ' ')


# ------------------------------------------------------------------ the output file
_lock = threading.Lock()


def append_row(path, row):
    """One line, written at once under a lock (the header first when the file is new)."""
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with _lock:
        new = not os.path.exists(path) or os.path.getsize(path) == 0
        with open(path, 'a') as f:
            if new:
                f.write('\t'.join(COLUMNS) + '\n')
            f.write('\t'.join(_fmt(row.get(c)) for c in COLUMNS) + '\n')


def read_rows(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def latest_ok(rows, key):
    """{SampleID: row} of the newest successful row of each sample measured under `key`."""
    out = {}
    for r in rows:
        if r.get('Status') == 'ok' and r.get('DefinitionKey') == key:
            out[r['SampleID']] = r
    return out


# ------------------------------------------------------------------ one sample
def measure_sample(env, cfg, ecfg, sample, windows, key, staged, want_bytes=None, table_coverage=None, tools='', keep_bam=False, threads=8, run=B.run_cmd, script=B.run_script):
    """The row of one sample: download, measure, delete. A failure is returned as a row with Status 'failed' and the reason; it is never raised."""
    t0 = time.time()
    row = dict(SampleID=sample, GenomeVersion=env.genome, DefinitionKey=key, MeasuredOn=datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), Machine=socket.gethostname(), Tools=tools,
               TableCoverage=table_coverage, BamBytes=want_bytes)
    d = f'{staged}{sample}/'
    try:
        B.stage_bam(env, staged, sample, want_bytes, run)
        bam = f'{d}{sample}.all.bam'
        row['BamBytes'] = os.path.getsize(bam)
        bamlist = f'{d}bams.list'
        open(bamlist, 'w').write(bam + '\n')
        mq = str(cfg['bcftools']['min_mapping_quality'])
        flag = parse_flagstat(B.must(run(['samtools', 'flagstat', '-@', str(threads), bam]), f'samtools flagstat of {sample}'))
        mapped_nd = int(B.must(run(['samtools', 'view', '-@', str(threads), '-c', '-F', str(VIEW_EXCLUDE), bam]), f'counting the reads of {sample}').strip())
        mapq_pass = int(B.must(run(['samtools', 'view', '-@', str(threads), '-c', '-F', str(VIEW_EXCLUDE), '-q', mq, bam]), f'counting the MAPQ >= {mq} reads of {sample}').strip())
        n, raw, hq = parse_window_line(B.must(script(window_script(env, bamlist, windows, cfg)), f'the pileup of {sample} over the windows'))
        row.update(derive(flag['primary'], flag['primary_duplicates'], mapped_nd, mapq_pass, n, raw, hq, sum(e - s + 1 for _, s, e in windows)))
        row['Status'], row['Error'] = 'ok', ''
    except Exception as e:                                          # noqa: one sample must never stop the others
        row['Status'], row['Error'] = 'failed', f'{type(e).__name__}: {e}'[:400]
    finally:
        if not keep_bam:
            shutil.rmtree(d, ignore_errors=True)
    row['Seconds'] = round(time.time() - t0, 1)
    return row


# ------------------------------------------------------------------ the run
class Plan:
    pass


def make_plan(env, cfg, ecfg, members, ledger, out_path, redo=False, jobs=None, threads=None, genome_fai=None):
    p = Plan()
    p.members = list(dict.fromkeys(members))
    p.ecfg = dict(ecfg)
    p.jobs, p.threads = jobs or ecfg['jobs'], threads or ecfg['threads']
    p.fai = genome_fai if genome_fai is not None else (B.read_fai(env.genome_fa + '.fai') if os.path.exists(env.genome_fa + '.fai') else None)
    p.windows = windows_from_fai(p.fai, ecfg['window_bp'], ecfg['position_fraction'], cfg['min_contig_bp']) if p.fai else []
    p.key = definition_key(cfg, ecfg, p.windows) if p.windows else ''
    p.out = out_path
    led = ledger.dropna(how='all').drop_duplicates('SampleID', keep='last').set_index('SampleID') if ledger is not None and len(ledger) else None
    p.sizes = {s: (B._bytes(led.at[s, 'BamSize']) if led is not None and s in led.index and 'BamSize' in led.columns else None) for s in p.members}
    p.table_cov = {s: (B._bytes(led.at[s, 'Coverage']) if led is not None and s in led.index and 'Coverage' in led.columns else None) for s in p.members}
    p.unknown = [s for s in p.members if led is None or s not in led.index]
    p.done = {} if redo else latest_ok(read_rows(out_path), p.key)
    p.todo = [s for s in p.members if s not in p.done]
    p.staged = f"{env.local_master}{ecfg['work_subdir']}"
    return p


def preflight(env, plan, free_bytes=None, which=shutil.which, run=B.run_cmd):
    """Problems that would stop a run, one sentence each; empty when it can start."""
    bad = []
    for t in TOOLS:
        if not which(t):
            bad.append(f'{t} is not on the PATH (activate the variant_calling environment, then `source tools/env.sh`)')
    if not bad:
        _, out, err = run(['bcftools', 'mpileup'])
        if '--ignore-RG' not in out + err:
            bad.append('this bcftools has no `mpileup --ignore-RG`: the same one is needed as for the cohort build')
    if not os.path.exists(env.genome_fa):
        bad.append(f'the reference FASTA is not at {env.genome_fa}')
    elif not plan.fai:
        bad.append(f'{env.genome_fa}.fai is missing (`samtools faidx {env.genome_fa}`)')
    if not plan.members:
        bad.append('no samples are selected')
    sizes = sorted([(plan.sizes.get(s) or plan.ecfg['default_bam_bytes']) for s in plan.todo], reverse=True)[:plan.jobs]
    plan.need_bytes = sum(sizes) * plan.ecfg['min_free_factor']
    if free_bytes is None:
        probe = plan.staged
        while probe and not os.path.exists(probe):
            probe = os.path.dirname(probe.rstrip('/'))
        free_bytes = shutil.disk_usage(probe or '/').free
    plan.free_bytes = free_bytes
    if plan.need_bytes > free_bytes:
        bad.append(f'the disk has {free_bytes / 1e9:.0f} GB free but {plan.jobs} BAMs at a time need about {plan.need_bytes / 1e9:.0f} GB (largest BAMs x {plan.ecfg["min_free_factor"]:g}): lower --jobs')
    return bad


def format_plan(plan, problems, cfg, running=False):
    b = cfg['bcftools']
    total = sum((plan.sizes.get(s) or plan.ecfg['default_bam_bytes']) for s in plan.todo)
    lines = [f"effective coverage: {len(plan.members)} samples ({len(plan.done)} already measured with this definition, {len(plan.todo)} to measure), definition key {plan.key or '?'}",
             f"  definition: the caller's own pileup (MAPQ >= {b['min_mapping_quality']}, base quality >= {b['min_base_quality']}, from cohort_build.yaml) over {len(plan.windows)} windows of {plan.ecfg['window_bp'] / 1e6:g} Mb "
             f"({plan.ecfg['position_fraction'] * 100:.0f}% along each chromosome); whole-BAM duplicate and MAPQ fractions",
             f"  download: {total / 1e9:.0f} GB in all, {plan.jobs} samples at a time ({plan.threads} samtools threads each), each BAM deleted after its measurement"
             + (f'; {len(plan.unknown)} samples are not in the AlignmentDatabase (their size is assumed {plan.ecfg["default_bam_bytes"] / 1e9:.0f} GB and the download is not verified)' if plan.unknown else ''),
             f'  output: {plan.out}']
    for p in problems:
        lines.append(f'  PROBLEM: {p}')
    if not problems:
        lines.append('  starting' if running else '  ready: add --run to measure')
    return '\n'.join(lines)


def tool_versions(run=B.run_cmd):
    out = []
    for t in ('samtools', 'bcftools'):
        try:
            rc, o, _ = run([t, '--version'])
            out.append(o.splitlines()[0].strip() if rc == 0 and o else f'{t}?')
        except OSError:
            out.append(f'{t}?')
    return '; '.join(out)


def run_effcov(env, cfg, plan, keep_bams=False, log=print, run=B.run_cmd, script=B.run_script):
    """Measure plan.todo, appending each row to plan.out as it finishes. Returns {'ok': n, 'failed': [(sample, reason)], 'skipped': n}."""
    os.makedirs(plan.staged, exist_ok=True)
    tools = tool_versions(run)
    total, count, t0 = len(plan.todo), [0], time.time()
    failed = []

    def one(s):
        row = measure_sample(env, cfg, plan.ecfg, s, plan.windows, plan.key, plan.staged, plan.sizes.get(s), plan.table_cov.get(s), tools, keep_bams, plan.threads, run, script)
        append_row(plan.out, row)
        count[0] += 1
        if row['Status'] == 'ok':
            tc = f" (table {row['TableCoverage']:.1f}x)" if row.get('TableCoverage') else ''
            log(f"[effcov] {count[0]}/{total} {s}: effective {row['EffectiveCoverage']:.2f}x{tc}, duplicates {row['DuplicateFraction']:.1%}, MAPQ>=13 {row['MapqPassFraction']:.1%}, {row['Seconds'] / 60:.1f} min")
        else:
            failed.append((s, row['Error']))
            log(f"[effcov] {count[0]}/{total} {s}: FAILED {row['Error']}")
        return row

    with ThreadPoolExecutor(max_workers=plan.jobs) as ex:
        list(ex.map(one, plan.todo))
    log(f'[effcov] {total - len(failed)} measured, {len(failed)} failed, {len(plan.done)} skipped (already measured), {(time.time() - t0) / 60:.1f} min')
    return dict(ok=total - len(failed), failed=failed, skipped=len(plan.done))


# ------------------------------------------------------------------ the report
def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def report_text(rows, key=None, sex=None, min_depth=5.0, members=None):
    """The table of the newest successful row of every sample (of one definition when `key` is given), and what stands out. `sex` is {SampleID: 'M'|'F'|...}, `members` restricts the samples."""
    best = {}
    for r in rows:
        if r.get('Status') == 'ok' and (key is None or r.get('DefinitionKey') == key) and (members is None or r['SampleID'] in members):
            best[r['SampleID']] = r
    if not best:
        return 'no successful measurement to report'
    sex = sex or {}
    head = f"{'sample':<18}{'sex':>4}{'table':>8}{'effective':>10}{'eff/table':>10}{'dup%':>7}{'mapq13%':>9}{'bases kept%':>12}{'covered%':>10}  "
    lines = [head + 'note', '-' * (len(head) + 8)]
    ratios, eff = [], {}
    for s in sorted(best):
        r = best[s]
        e, t = _num(r['EffectiveCoverage']), _num(r['TableCoverage'])
        eff[s] = e
        ratio = e / t if (e is not None and t) else None
        if ratio is not None:
            ratios.append(ratio)
        note = f'below {min_depth:g}x' if (e is not None and e < min_depth) else ''
        lines.append(f"{s:<18}{sex.get(s, ''):>4}{(f'{t:.1f}' if t else ''):>8}{e:>10.2f}{(f'{ratio:.2f}' if ratio is not None else ''):>10}{100 * _num(r['DuplicateFraction']):>7.1f}"
                     f"{100 * _num(r['MapqPassFraction']):>9.1f}{100 * _num(r['BaseSurvivalFraction']):>12.1f}{100 * _num(r['WindowCoveredFraction']):>10.1f}  {note}")
    low = [s for s, e in eff.items() if e is not None and e < min_depth]
    lines.append('')
    lines.append(f'{len(best)} samples; effective coverage {min(eff.values()):.2f} to {max(eff.values()):.2f}x (median {statistics.median(eff.values()):.2f}x)')
    if ratios:
        lines.append(f'effective / table: {min(ratios):.2f} to {max(ratios):.2f} (median {statistics.median(ratios):.2f})')
    lines.append(f'{len(low)} below {min_depth:g}x (sexloci needs depth >= {min_depth:g} at a genotype)' + (': ' + ', '.join(sorted(low)) if low else ''))
    for sx in sorted({v for v in sex.values() if v}):
        vals = [eff[s] for s in eff if sex.get(s) == sx and eff[s] is not None]
        if vals:
            lines.append(f'  sex {sx}: {len(vals)} samples, median effective {statistics.median(vals):.2f}x, {sum(1 for v in vals if v < min_depth)} below {min_depth:g}x')
    return '\n'.join(lines)
