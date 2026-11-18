"""Read-back of the SampleDatabase against the stored panel JSONs, for the writer with `nocall` and `scope` (module 05, config/panelpass.yaml). Independent of the writer: it recomputes every call from the
JSON and derives the expected cell from the call and from the config (not from the writer's code).
   PANEL=<folder of *.panel.json>  SHEET=<export of the SampleDatabase tab, strings as displayed>  [CONFIG=<panelpass.yaml>]   python readback_v2.py      (after: source tools/env.sh)"""
import glob, os, pandas as pd, yaml
from sd_genotypes.panel import calls
cfg = yaml.safe_load(open(os.environ.get('CONFIG') or os.path.join(os.environ['CICHLID_REPO'], '05_variant_calling/config/panelpass.yaml')))['sheet']
write_cols, nc, scope = cfg['write_columns'], cfg.get('nocall') or {}, cfg.get('scope') or {}
sheet = pd.read_csv(os.environ['SHEET'], sep='\t', dtype=str, keep_default_na=False).set_index('SampleID')
if scope and 'Ecogroup' not in sheet.columns: raise SystemExit('the export has no Ecogroup column, which sheet.scope needs')
stat, bad = {}, []
def count(col, key): stat.setdefault(col, {}).setdefault(key, 0); stat[col][key] += 1
for f in sorted(glob.glob(os.environ['PANEL'] + '/*.panel.json')):
    r = calls.call_panel(f); s = r['sample']
    if s not in sheet.index: continue
    eco = sheet.at[s, 'Ecogroup'] if 'Ecogroup' in sheet.columns else ''
    for c in r['calls']:
        col, fcol, want, fl = c['column'], c['flags_column'], c['value_to_write'], set(c['flags'])
        if col not in write_cols: continue
        out_of_scope = col in scope and eco not in scope[col]
        if out_of_scope: exp_val, exp_flags, kind = '', '', 'out_of_scope'
        elif want is not None: exp_val, exp_flags, kind = want, ';'.join([c['call_tag']] + c['flags']), 'value'
        elif col in nc.get('columns', []) and (fl & set(nc['flags'])) and not (fl & set(nc.get('exclude_flags', []))): exp_val, exp_flags, kind = nc['value'], ';'.join([c['call_tag']] + c['flags']), 'nocall'
        else: exp_val, exp_flags, kind = '', ';'.join([c['call_tag']] + c['flags']), 'empty'
        count(col, 'expected_' + kind)
        cell = sheet.at[s, col] if col in sheet.columns else ''
        try: same = (abs(float(cell) - float(exp_val)) < 1e-3) if isinstance(exp_val, float) else (cell == str(exp_val))
        except ValueError: same = False
        if same: count(col, 'ok')
        else: count(col, 'wrong'); bad.append((s, col, cell, f'expected {exp_val!r} ({kind})'))
        fc = sheet.at[s, fcol] if fcol in sheet.columns else ''
        if fc == exp_flags: count(col, 'flags_ok')
        else: count(col, 'flags_wrong'); bad.append((s, fcol, fc, f'expected {exp_flags!r}'))
t = pd.DataFrame(stat).T.fillna(0).astype(int); print(t.to_string()); print(len(bad), 'differences'); [print(' ', b) for b in bad[:20]]
print('PASS' if not bad else 'FAIL')
