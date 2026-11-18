"""Module 04's panel over EXISTING BAMs, and its calls into the SampleDatabase (module 05).

Why here: the alignment pipeline (03) makes the panel's measurements for the samples it aligns, in the same pass, while the BAM is on local disk. The BAMs aligned before the panel existed (the 573 mbuna,
5.8 TB on Dropbox) need the same measurement made later, one BAM at a time: that is what `effcov measure` already does for another measurement (download, measure, delete, a log that resumes), so it is done
here. No panel logic is copied: the commands, the preflight that stops a run when the definition keys are not the frozen ones, and the rule that says whether a JSON belongs to a BAM are module 03's
`align_pipeline.panel`; the measurement is `python -m sd_genotypes panel measure` run as a subprocess in whatever interpreter has module 04's dependencies (`--panel-python`); the calls are module 04's
`sd_genotypes.panel.calls` (decision 76), applied to the stored JSON and never to a BAM.

  measure  per sample: free-disk check, download `<sample>.all.bam` + `.all.bai`, `panel measure --modules all` into `<json dir>/<sample>.panel.json`, check that the JSON names THIS BAM, copy the JSON to
           Dropbox next to the BAM (the durable record: Temp is not a record) and check its size there, delete the BAM. A log (TSV, one row per sample as it finishes) makes a stopped run resume; a sample is
           measured again when the definition keys changed, the BAM changed (a realigned sample keeps its name and not its size) or its JSON is gone. One failed sample never stops the others.
  write    the calls of the JSONs into the SampleDatabase (or a copy of the tab) through the guarded sheet writer: empty cells only, a call with a blocking flag writes its flags cell and not its value, the
           flags cell starts with the call tag (`pilot-v1`), a different value already in a cell is reported and never overwritten (and its flags cell is then left alone too: a panel tag next to a
           hand-entered value would claim a call that is not in the sheet). Only the columns of `sheet.write_columns` are written; Inversion5 and B score stay in the JSONs.
  table    every call of the JSONs (all six columns, the flags and the numbers behind them) in one TSV, for reading the distributions.
"""
import csv
import datetime
import json
import os
import re
import shutil
import socket
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import yaml

from cichlid_common.common_light import config_path
from . import build as B

DEFAULT_GENOME = 'Mzebra_GT3_NCBI'
DEFAULTS = dict(jobs=2, min_free_factor=1.15, default_bam_bytes=12.0e9, work_subdir='PanelPass/staged/', json_subdir='Outputs/PanelJSON/', log_name='panel_pass.{genome}.tsv', modules='all', timeout_s=1800,
                expect_keys=[], sheet=dict(tab='SampleDatabase', key='SampleID', write_columns=['gsdf number', 'chr7ins', 'TE1', 'TE2'], flags_in_sheet=True, nocall={}, scope={}))
LOG_COLUMNS = ['SampleID', 'GenomeVersion', 'Status', 'Error', 'BamName', 'BamBytes', 'JsonBytes', 'DefinitionKeys', 'DownloadSeconds', 'MeasureSeconds', 'UploadSeconds', 'Seconds', 'MeasuredOn', 'Machine']
KEY_LINE = re.compile(r'^\[panel\]\s+(\w+)\s+OK\s+([0-9a-f]{12})\b', re.M)
_LOCK = threading.Lock()


class PanelPassError(RuntimeError):
    pass


def load_config(path=None):
    cfg = {k: (dict(v) if isinstance(v, dict) else (list(v) if isinstance(v, list) else v)) for k, v in DEFAULTS.items()}
    p = path or config_path('panelpass.yaml')
    if os.path.exists(p):
        for k, v in (yaml.safe_load(open(p)) or {}).items():
            if isinstance(v, dict) and isinstance(cfg.get(k), dict):
                cfg[k].update(v)
            else:
                cfg[k] = v
    for k in ('jobs', 'timeout_s'):
        try:
            cfg[k] = int(float(cfg[k]))
        except (TypeError, ValueError):
            raise ValueError(f'panelpass config: {k} must be a number, not {cfg[k]!r}')
    for k in ('min_free_factor', 'default_bam_bytes'):
        try:
            cfg[k] = float(cfg[k])                                  # YAML reads 12.0e9 as a float, 12e9 as a string
        except (TypeError, ValueError):
            raise ValueError(f'panelpass config: {k} must be a number, not {cfg[k]!r}')
    cfg['expect_keys'] = [str(x) for x in (cfg['expect_keys'] or [])]
    if not cfg['sheet'].get('write_columns'):
        raise ValueError('panelpass config: sheet.write_columns is empty: nothing would ever be written')
    sh = cfg['sheet']
    if not isinstance(sh.get('flags_in_sheet', True), bool):
        raise ValueError(f"panelpass config: sheet.flags_in_sheet must be true or false, not {sh['flags_in_sheet']!r}")
    nc = sh.get('nocall') or {}
    if nc:
        for k in ('value', 'columns', 'flags'):
            if not nc.get(k):
                raise ValueError(f'panelpass config: sheet.nocall.{k} is missing')
        bad = [c for c in nc['columns'] if c not in sh['write_columns']]
        if bad:
            raise ValueError(f'panelpass config: sheet.nocall.columns {bad} are not in sheet.write_columns')
        nc = dict(value=str(nc['value']), columns=list(nc['columns']), flags=list(nc['flags']), exclude_flags=list(nc.get('exclude_flags', ['site_error'])))
    sh['nocall'] = nc
    sc = sh.get('scope') or {}
    for col, groups in sc.items():
        if col not in sh['write_columns']:
            raise ValueError(f'panelpass config: sheet.scope names {col!r}, which is not in sheet.write_columns')
        if not isinstance(groups, (list, tuple)) or not groups:
            raise ValueError(f'panelpass config: sheet.scope.{col} must be a non-empty list of Ecogroups')
    sh['scope'] = {c: [str(g) for g in gs] for c, gs in sc.items()}
    return cfg


# ------------------------------------------------------------------ the panel (module 04 through module 03's wrapper)
def panel_module():
    try:
        from align_pipeline import panel as pn
    except ImportError as e:
        raise PanelPassError(f'module 03 (align_pipeline.panel) cannot be imported: `source tools/env.sh` puts every module on PYTHONPATH ({e})')
    return pn


def make_panel_config(pcfg, python=None, pythonpath=None):
    return panel_module().PanelConfig(python=python, modules=pcfg['modules'], pythonpath=pythonpath, timeout_s=pcfg['timeout_s'], expect=pcfg['expect_keys'])


def parse_keys(text):
    """{module: definition key} from what `panel check` printed ('[panel] gsdf     OK       5d71ee8d6c2b')."""
    return {m: k for m, k in KEY_LINE.findall(text or '')}


def keys_string(keys):
    return ';'.join(f'{m}={k}' for m, k in sorted(keys.items()))


def check_panel(panel, runner=subprocess.run):
    """(ok, text, keys): `panel check` once, with the expected-keys rule of module 03."""
    ok, text = panel_module().preflight(panel, run=runner)
    return ok, text, (parse_keys(text) if ok else {})


# ------------------------------------------------------------------ the log: one row per sample as it finishes
def read_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def append_row(path, row):
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with _LOCK:
        new = not os.path.exists(path) or os.path.getsize(path) == 0
        with open(path, 'a', newline='') as f:
            w = csv.DictWriter(f, fieldnames=LOG_COLUMNS, delimiter='\t', extrasaction='ignore')
            if new:
                w.writeheader()
            w.writerow({k: ('' if row.get(k) is None else row[k]) for k in LOG_COLUMNS})


# ------------------------------------------------------------------ the plan
class Plan:
    pass


def make_plan(env, pcfg, members, ledger, keys=None, redo=False, jobs=None, limit=None, json_dir=None):
    p = Plan()
    p.genome, p.members, p.jobs = env.genome, list(dict.fromkeys(members)), jobs or pcfg['jobs']
    p.keys = dict(keys or {})
    p.keys_text = keys_string(p.keys)
    p.json_dir = json_dir or f"{env.local_master}{pcfg['json_subdir']}{env.genome}/"
    p.log_path = f"{p.json_dir}{pcfg['log_name'].format(genome=env.genome)}"
    p.staged = f"{env.local_master}{pcfg['work_subdir']}"
    led = ledger.dropna(how='all').drop_duplicates('SampleID', keep='last').set_index('SampleID') if ledger is not None and len(ledger) else None
    p.sizes = {s: (B._bytes(led.at[s, 'BamSize']) if led is not None and s in led.index and 'BamSize' in led.columns else None) for s in p.members}
    p.unknown = [s for s in p.members if led is None or s not in led.index]
    last = {r['SampleID']: r for r in read_rows(p.log_path)}
    p.done, p.stale = {}, {}
    for s in p.members:
        r = last.get(s)
        if redo or r is None or r.get('Status') != 'ok':
            continue
        why = None
        if p.keys_text and r.get('DefinitionKeys') != p.keys_text:
            why = 'measured under other definition keys'
        elif p.sizes.get(s) and r.get('BamBytes') not in ('', None) and int(float(r['BamBytes'])) != int(p.sizes[s]):
            why = 'the BAM has changed since (realigned?)'
        elif not os.path.exists(f'{p.json_dir}{s}.panel.json'):
            why = 'its JSON is not in the local folder'
        if why:
            p.stale[s] = why
        else:
            p.done[s] = r
    p.todo_all = [s for s in p.members if s not in p.done]
    p.limit = limit
    p.todo = p.todo_all[:limit] if limit is not None else list(p.todo_all)
    return p


def _free(path, free_bytes=None):
    if free_bytes is not None:
        return free_bytes
    probe = path
    while probe and not os.path.exists(probe):
        probe = os.path.dirname(probe.rstrip('/'))
    return shutil.disk_usage(probe or '/').free


def preflight(env, pcfg, plan, panel_ok, panel_text, which=shutil.which, free_bytes=None, upload=True):
    """Problems that would stop a run, one sentence each; empty when it can start."""
    bad = []
    for t in ('rclone', 'samtools'):
        if not which(t):
            bad.append(f'{t} is not on the PATH (a BAM is downloaded and checked here; the panel itself runs in --panel-python)')
    if not panel_ok:
        bad.append('the panel is not ready: ' + (panel_text or '').strip()[-700:])
    if not plan.members:
        bad.append('no samples are selected')
    sizes = sorted([(plan.sizes.get(s) or pcfg['default_bam_bytes']) for s in plan.todo], reverse=True)[:plan.jobs]
    plan.need_bytes = sum(sizes) * pcfg['min_free_factor']
    plan.free_bytes = _free(plan.staged, free_bytes)
    if plan.need_bytes > plan.free_bytes:
        bad.append(f'the disk of {plan.staged} has {plan.free_bytes / 1e9:.0f} GB free but {plan.jobs} BAMs at a time need about {plan.need_bytes / 1e9:.0f} GB (largest BAMs x {pcfg["min_free_factor"]:g}): lower --jobs')
    return bad


def format_plan(plan, problems, pcfg, running=False, upload=True):
    total = sum((plan.sizes.get(s) or pcfg['default_bam_bytes']) for s in plan.todo)
    lines = [f"panel pass, {plan.genome}: {len(plan.members)} samples ({len(plan.done)} already measured, {len(plan.todo_all)} to measure"
             + (f'; this run takes the first {len(plan.todo)}' if plan.limit is not None and len(plan.todo) < len(plan.todo_all) else '') + ')',
             f"  definition keys: {plan.keys_text or '? (the panel is not ready)'}",
             f"  download: {total / 1e9:.0f} GB in all, {plan.jobs} BAMs at a time, each deleted as soon as it is measured; staged under {plan.staged}"
             + (f'; {len(plan.unknown)} samples are not in the AlignmentDatabase (size assumed {pcfg["default_bam_bytes"] / 1e9:.0f} GB, download not verified)' if plan.unknown else ''),
             f"  JSONs: {plan.json_dir}" + ('; each is also copied to Dropbox next to its BAM' if upload else '; NOT copied to Dropbox (--no-upload)'), f'  log: {plan.log_path}']
    for s, why in list(plan.stale.items())[:6]:
        lines.append(f'  measured again: {s} ({why})')
    if len(plan.stale) > 6:
        lines.append(f'  ... and {len(plan.stale) - 6} more measured again')
    for p in problems:
        lines.append(f'  PROBLEM: {p}')
    if not problems:
        lines.append('  starting' if running else '  ready: add --run to measure')
    return '\n'.join(lines)


# ------------------------------------------------------------------ one sample, and the run
def measure_sample(env, pcfg, panel, sample, plan, want_bytes=None, upload=True, keep_bam=False, run=B.run_cmd, runner=subprocess.run, free_bytes=None):
    """The log row of one sample: download, measure, check, copy to Dropbox, delete the BAM. A failure is returned as a row with Status 'failed' and the reason; it is never raised."""
    pn = panel_module()
    t0 = time.time()
    row = dict(SampleID=sample, GenomeVersion=env.genome, BamName=f'{sample}.all.bam', BamBytes=want_bytes, DefinitionKeys=plan.keys_text, MeasuredOn=datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
               Machine=socket.gethostname())
    d = f'{plan.staged}{sample}/'
    part = f'{plan.json_dir}{sample}.panel.json.part'
    final = f'{plan.json_dir}{sample}.panel.json'
    try:
        os.makedirs(plan.json_dir, exist_ok=True)
        need = (want_bytes or pcfg['default_bam_bytes']) * pcfg['min_free_factor']
        free = _free(plan.staged, free_bytes)
        if free < need:
            raise PanelPassError(f'only {free / 1e9:.0f} GB free on {plan.staged}: {sample} needs about {need / 1e9:.0f} GB')
        t = time.time()
        B.stage_bam(env, plan.staged, sample, want_bytes, run)
        row['DownloadSeconds'] = round(time.time() - t, 1)
        bam = f'{d}{sample}.all.bam'
        row['BamBytes'] = os.path.getsize(bam)
        if os.path.exists(part):
            os.remove(part)
        t = time.time()
        try:
            r = runner(pn.measure_command(panel, bam, sample, part), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=panel.timeout_s)
        except subprocess.TimeoutExpired:
            raise PanelPassError(f'`panel measure` of {sample} did not finish within {panel.timeout_s} s')
        row['MeasureSeconds'] = round(time.time() - t, 1)
        if r.returncode != 0 or not os.path.exists(part):
            raise PanelPassError(f'`panel measure` of {sample} failed (exit {r.returncode}): {(r.stdout or "").strip()[-700:]}')
        status, why = pn.identity_status(part, bam)
        if status != 'same':
            raise PanelPassError(f'the JSON does not belong to the BAM ({status}): {why}')
        os.replace(part, final)
        row['JsonBytes'] = os.path.getsize(final)
        if upload:
            t = time.time()
            dest = f'{env.cloud_sample_dir(sample)}{sample}.panel.json'
            B.must(run(['rclone', 'copyto', final, dest]), f'copying the panel JSON of {sample} to Dropbox')
            rc, out, err = run(['rclone', 'lsf', '--format', 's', dest])
            got = out.strip().splitlines()[0].strip() if rc == 0 and out.strip() else ''
            if got != str(row['JsonBytes']):
                raise PanelPassError(f'the JSON of {sample} on Dropbox has {got or "no"} bytes, the local one {row["JsonBytes"]}: the copy is not trusted')
            row['UploadSeconds'] = round(time.time() - t, 1)
        row['Status'], row['Error'] = 'ok', ''
    except Exception as e:                                          # noqa: one sample must never stop the others
        row['Status'], row['Error'] = 'failed', f'{type(e).__name__}: {e}'[:600]
        for p in (part,):
            if os.path.exists(p):
                os.remove(p)
    finally:
        if not keep_bam:
            shutil.rmtree(d, ignore_errors=True)
    row['Seconds'] = round(time.time() - t0, 1)
    return row


def run_pass(env, pcfg, plan, panel, upload=True, keep_bams=False, log=print, run=B.run_cmd, runner=subprocess.run, free_bytes=None):
    """Measure plan.todo, appending each row to the log as it finishes. Returns {'ok': n, 'failed': [(sample, reason)], 'skipped': n, 'seconds': wall seconds, 'bytes': downloaded}."""
    os.makedirs(plan.staged, exist_ok=True)
    os.makedirs(plan.json_dir, exist_ok=True)
    total, count, t0, failed, nbytes, dl_seconds = len(plan.todo), [0], time.time(), [], [0], [0.0]

    def one(s):
        row = measure_sample(env, pcfg, panel, s, plan, plan.sizes.get(s), upload, keep_bams, run, runner, free_bytes)
        append_row(plan.log_path, row)
        count[0] += 1
        try:                                                        # a problem in the log line must never stop the run: the row is already in the log
            if row['Status'] == 'ok':
                gb = float(row['BamBytes'] or 0) / 1e9
                dl = float(row.get('DownloadSeconds') or 0)
                nbytes[0] += float(row['BamBytes'] or 0)
                dl_seconds[0] += dl
                rate = f'{gb * 1000 / dl:.0f} MB/s' if dl >= 0.5 else 'no download'     # a BAM kept from an earlier run, or one too fast to time
                log(f"[panel] {count[0]}/{total} {s}: {gb:.1f} GB, download {dl:.0f} s ({rate}), measure {float(row.get('MeasureSeconds') or 0):.0f} s, "
                    f"upload {float(row.get('UploadSeconds') or 0):.0f} s, {row['Seconds'] / 60:.1f} min")
            else:
                failed.append((s, row['Error']))
                log(f"[panel] {count[0]}/{total} {s}: FAILED {row['Error']}")
        except Exception as e:                                      # noqa
            if row['Status'] != 'ok':
                failed.append((s, row['Error']))
            log(f'[panel] {count[0]}/{total} {s}: {row["Status"]} (the log line could not be made: {type(e).__name__}: {e})')
        return row

    with ThreadPoolExecutor(max_workers=plan.jobs) as ex:
        list(ex.map(one, plan.todo))
    wall = time.time() - t0
    if nbytes[0] and dl_seconds[0] > 0:
        log(f"[panel] {total - len(failed)} measured, {len(failed)} failed, {len(plan.done)} skipped (already measured), {wall / 60:.1f} min; {nbytes[0] / 1e9:.0f} GB downloaded: "
            f"{nbytes[0] / 1e6 / max(wall, 1e-9):.0f} MB/s over the whole run, {nbytes[0] / 1e6 / dl_seconds[0]:.0f} MB/s per download stream")
    else:
        log(f'[panel] {total - len(failed)} measured, {len(failed)} failed, {len(plan.done)} skipped (already measured), {wall / 60:.1f} min')
    return dict(ok=total - len(failed), failed=failed, skipped=len(plan.done), seconds=wall, bytes=nbytes[0])


# ------------------------------------------------------------------ the calls (module 04, stored measurements only)
def load_calls_module():
    try:
        import importlib
        calls = importlib.import_module('sd_genotypes.panel.calls')
    except ImportError as e:
        raise PanelPassError('sd_genotypes.panel.calls cannot be imported: `source tools/env.sh` puts module 04 on PYTHONPATH, and the call layer (decision 76) must be in this checkout '
                             f'(it needs numpy, pandas and pyyaml): {e}')
    return calls


def list_json_samples(json_dir):
    if not os.path.isdir(json_dir):
        return []
    return sorted(f[:-len('.panel.json')] for f in os.listdir(json_dir) if f.endswith('.panel.json'))


def record_updates(r, calls_mod, write_columns, sheet_cfg=None, ecogroup=None):
    """({column: value} to write, {flags column: its value column}) for one call record. Values are only the CONFIDENT ones (calls.sheet_row). Options of the sheet config:
    `nocall`: a genotype column whose call has a flag of nocall.flags (no call, low GQ, low depth: too little data) and none of nocall.exclude_flags (a site_error is an error of the site, and an
    unvalidated_value is a value that is held back on purpose: neither is a lack of data) gets nocall.value ('./.') instead of staying empty;
    `flags_in_sheet: false`: no flags cells (they stay in the JSON and in `panel table`);
    `scope`: {column: [Ecogroups]}: a column is written only for samples of those Ecogroups (neither its value nor its flags cell otherwise); a sample whose Ecogroup is unknown is outside every scope."""
    sc = sheet_cfg or {}
    scope = sc.get('scope') or {}
    cols = [c for c in write_columns if c not in scope or ecogroup in scope[c]]
    values, flags = calls_mod.sheet_row(r)
    upd = {k: v for k, v in values.items() if k in cols}
    nc = sc.get('nocall') or {}
    for c in (r['calls'] if nc else []):
        col, fl = c['column'], set(c.get('flags') or [])
        if col in cols and col in nc['columns'] and col not in upd and (fl & set(nc['flags'])) and not (fl & set(nc.get('exclude_flags', []))):
            upd[col] = nc['value']
    pairs = {}
    if sc.get('flags_in_sheet', True):
        pairs = {c['flags_column']: c['column'] for c in r['calls'] if c['column'] in cols}
        upd.update({k: v for k, v in flags.items() if k in pairs})
    return upd, pairs


def collect_updates(json_dir, samples, calls_mod, write_columns, sheet_cfg=None, eco_of=None):
    """(updates {sample: {column: value}}, flag_of {(sample, flags column): value column}, info {'no_json': [...], 'refused': [(sample, why)], 'tags': [the call tags of the calls written],
    'scoped_out': {column: number of samples left out of it by sheet.scope}}). eco_of: {SampleID: Ecogroup}, needed only when sheet.scope is set."""
    updates, flag_of, info, tags = {}, {}, dict(no_json=[], refused=[], scoped_out={}), set()
    scope = (sheet_cfg or {}).get('scope') or {}
    if scope and eco_of is None:
        raise ValueError('sheet.scope is set but no Ecogroup of the samples was given')
    for s in samples:
        path = f'{json_dir}{s}.panel.json'
        if not os.path.exists(path):
            info['no_json'].append(s)
            continue
        eco = (eco_of or {}).get(s)
        try:
            r = calls_mod.call_panel(path)
            upd, pairs = record_updates(r, calls_mod, write_columns, sheet_cfg, eco)
        except Exception as e:                                      # noqa: a measurement made under another definition is refused (DefinitionMismatch); say so and go on
            info['refused'].append((s, f'{type(e).__name__}: {e}'[:300]))
            continue
        for col, groups in scope.items():
            if eco not in groups:
                info['scoped_out'][col] = info['scoped_out'].get(col, 0) + 1
        if upd:
            updates[s] = upd
            tags.update(c['call_tag'] for c in r['calls'] if c['column'] in write_columns and c.get('call_tag'))
            for fc, vc in pairs.items():
                flag_of[(s, fc)] = vc
    info['tags'] = sorted(tags)
    return updates, flag_of, info


def column_order(write_columns, flag_of):
    """The order new columns are added to the header: each value column followed by its flags column (gsdf number, gsdf flags, chr7ins, chr7ins flags, ...)."""
    order = []
    for vc in write_columns:
        order.append(vc)
        order += sorted({fc for (_, fc), v in flag_of.items() if v == vc})
    return order


def drop_flags_of_conflicts(cu, updates, flag_of):
    """Plan `updates`; for every value cell that CONFLICTS (a different value is already there) leave its flags cell alone too (a panel tag next to a hand-entered, different value would claim a call that is not
    in the sheet); return (the final plan, the number of flags cells dropped). Mutates `updates`."""
    plan = cu.plan(updates)
    conflicted = {(c['sample'], c['column']) for c in plan.changes if c['status'] == 'conflict'}
    drop = [(s, fc) for (s, fc), vc in flag_of.items() if (s, vc) in conflicted and fc in updates.get(s, {})]
    for s, fc in drop:
        del updates[s][fc]
    return (cu.plan(updates) if drop else plan), len(drop)


def flatten_details(details, prefix):
    out = {}
    if isinstance(details, dict):
        for k, v in details.items():
            if isinstance(v, (int, float, str, bool)) and not isinstance(v, bool) or v is None:
                out[f'{prefix}.{k}'] = '' if v is None else v
    return out


def table_rows(json_dir, samples, calls_mod, sample_info=None):
    """([row dicts], [columns]): every call of every JSON, the flags cells as `sheet_row` makes them, and the numbers of each call's details."""
    rows, cols = [], ['SampleID']
    for s in samples:
        path = f'{json_dir}{s}.panel.json'
        if not os.path.exists(path):
            continue
        row = {'SampleID': s}
        row.update((sample_info or {}).get(s, {}))
        try:
            r = calls_mod.call_panel(path)
            values, flags = calls_mod.sheet_row(r)
            for c in r['calls']:
                row[c['column']] = '' if c['value'] is None else c['value']
                row[f"{c['column']} confident"] = bool(c['confident'])
                row[c['flags_column']] = flags.get(c['flags_column'], '')
                row.update(flatten_details(c.get('details'), c['column']))
        except Exception as e:                                      # noqa
            row['Error'] = f'{type(e).__name__}: {e}'[:300]
        rows.append(row)
        for k in row:
            if k not in cols:
                cols.append(k)
    return rows, cols


def summary_text(rows, columns):
    """Per value column: the values of the confident calls and the commonest flags."""
    value_cols = [c for c in columns if f'{c} confident' in columns]
    lines = [f'{len(rows)} samples']
    for c in value_cols:
        conf = [r for r in rows if r.get(f'{c} confident') is True or r.get(f'{c} confident') == 'True']
        vals = {}
        for r in conf:
            vals[str(r.get(c))] = vals.get(str(r.get(c)), 0) + 1
        flag_col = next((k for k in columns if k.endswith('flags') and k.startswith(c.split()[0])), None)
        fl = {}
        for r in rows:
            for t in str(r.get(flag_col, '')).split(';')[1:] if flag_col else []:
                if t:
                    fl[t] = fl.get(t, 0) + 1
        top = ', '.join(f'{k} {v}' for k, v in sorted(fl.items(), key=lambda kv: -kv[1])[:5])
        if c in ('B score',) or all(_isnum(r.get(c)) for r in conf if r.get(c) != ''):
            nums = sorted(float(r[c]) for r in conf if _isnum(r.get(c)))
            body = f'min {nums[0]:.3g}, median {nums[len(nums) // 2]:.3g}, max {nums[-1]:.3g}' if nums and len(set(map(str, nums))) > 8 else ', '.join(f'{k}: {v}' for k, v in sorted(vals.items()))
        else:
            body = ', '.join(f'{k}: {v}' for k, v in sorted(vals.items()))
        lines.append(f'  {c}: {len(conf)} confident ({body})' + (f'; flags: {top}' if top else ''))
    return '\n'.join(lines)


def _isnum(x):
    try:
        float(x)
        return x != ''
    except (TypeError, ValueError):
        return False


def ecogroups_of(ws, key='SampleID', column='Ecogroup'):
    """{SampleID: Ecogroup} of the SampleDatabase tab, from the same worksheet the cells are written to."""
    vals = ws.get_all_values()
    head = [str(h).strip() for h in vals[0]]
    if key not in head or column not in head:
        raise ValueError(f'the tab has no {key!r} and {column!r} columns')
    ik, ie = head.index(key), head.index(column)
    return {str(r[ik]).strip(): (str(r[ie]).strip() if len(r) > ie else '') for r in vals[1:] if len(r) > ik and str(r[ik]).strip()}
