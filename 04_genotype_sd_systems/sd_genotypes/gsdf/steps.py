"""gsdf steps. Offline-capable: `call` reads the panel JSONs, `validate` reads a calls table and a sample table. With a Project the outputs go to an immutable run folder (Dropbox); with
--out-dir they are plain files."""
import os
import numpy as np
import pandas as pd
import yaml
from cichlid_common.common_light import config_path
from . import call as C, measure as M, validate as V
from ..panel import measure as PM


def cfg_of():
    return yaml.safe_load(open(config_path('gsdf.yaml')))


def _out(proj, analysis, step, out_dir, params):
    """Context manager: a Run folder with a Project, else a plain directory."""
    if proj is not None:
        return proj.run(analysis, step, params=params)
    if not out_dir:
        raise SystemExit('give --out-dir (offline) or run with a Project')
    os.makedirs(out_dir, exist_ok=True)
    class _Dir:
        dir = out_dir.rstrip('/') + '/'
        def __enter__(self): return self
        def __exit__(self, *a): return False
    return _Dir()


def call_step(measurements, out_dir=None, proj=None, cfg=None):
    cfg = cfg or cfg_of()
    rows, skipped = PM.collect(measurements, 'gsdf')
    if not rows:
        raise SystemExit(f'no gsdf measurements in {measurements}')
    key = M.definition_key(cfg)
    wrong = [r['SampleID'] for r in rows if r.get('definition') != key]
    if wrong:
        raise SystemExit(f'{len(wrong)} measurements were made under a different gsdf definition than config/gsdf.yaml (region, bins, windows or MAPQ changed): {wrong[:5]} ...; measure them again')
    dup = pd.Series([r['SampleID'] for r in rows]); dup = dup[dup.duplicated()].unique().tolist()
    if dup:
        raise SystemExit(f'sample measured more than once in {measurements}: {dup[:5]}')
    calls = C.call_table(rows, cfg)
    calls = calls.merge(pd.DataFrame([dict(SampleID=r['SampleID'], bam=r['bam'], bam_size_bytes=r['bam_size_bytes'], measured=r['measured'], panel_version=r['panel_version']) for r in rows]), on='SampleID')
    clus = C.cluster_check(calls, cfg)
    fl = calls['call_flags'].str.split(';').explode(); fl = fl[fl != '']
    text = [f'gsdf copy number: {len(calls)} samples measured ({len(skipped)} panel files without a gsdf measurement skipped)',
            'copies (rounded): ' + ', '.join(f'{int(k)}: {v}' for k, v in calls.gsdf_number.value_counts().sort_index().items()),
            f"unflagged: {int((calls['call_flags'] == '').sum())}; flagged: {int((calls['call_flags'] != '').sum())}" + ('  (' + ', '.join(f'{k}: {v}' for k, v in fl.value_counts().items()) + ')' if len(fl) else ''),
            f"junction cross-check: {int((calls['check_flags'] != '').sum())} calls with a unique-read junction in a sample called with two copies (check_flags; these calls still count as unflagged)",
            '', 'where the unflagged samples sit against the calibration (drift check for a new batch; ratio and junction ratio by call):', clus.round(3).to_string(index=False), '']
    with _out(proj, 'gsdf', 'call', out_dir, dict(cfg=cfg, n=len(calls))) as run:
        calls.to_csv(run.dir + 'gsdf_calls.tsv', sep='\t', index=False)
        calls.assign(**{'gsdf number': calls.gsdf_number, 'gsdf flags': calls['call_flags']})[['SampleID', 'gsdf number', 'gsdf flags']].to_csv(run.dir + 'sheet_update.tsv', sep='\t', index=False)
        open(run.dir + 'summary.txt', 'w').write('\n'.join(text) + '\n')
    print('\n'.join(text))
    return calls


def validate_step(calls_path, table_path, out_dir=None, proj=None, cfg=None):
    """calls_path may be None: the sheet's own `gsdf number` is then used as the copy number (sex x species table of the existing values)."""
    cfg = cfg or cfg_of(); vc = cfg['validate']
    sheet = V.sheet_frame(pd.read_csv(table_path, sep='\t', low_memory=False))
    text = []
    if calls_path:
        calls = pd.read_csv(calls_path, sep='\t', keep_default_na=False, na_values=['']).fillna({'call_flags': ''})
        summ, mism = V.compare_with_sheet(calls, sheet)
        text += [f"agreement with the existing `gsdf number`: {summ['n_agree']} of {summ['n_compared']} unflagged calls with a sheet value agree ({summ['frac_agree']:.1%}); "
                 f"{summ['n_with_sheet_value']} samples have both; {summ['n_flagged_among_with_sheet']} of those are flagged and left out", '']
        use = calls.merge(sheet[['SampleID', 'Sex', 'Species', 'Ecogroup']], on='SampleID', how='left').rename(columns={'gsdf_number': 'copies'})
        use = use[use['call_flags'] == '']
    else:
        mism = pd.DataFrame()
        use = sheet.rename(columns={'sheet_gsdf': 'copies'})
        text += ['no calls given: sex x copy number of the EXISTING `gsdf number` values', '']
    tab, disc = V.sex_table(use[['SampleID', 'Sex', 'Species', 'Ecogroup', 'copies']], vc['min_per_sex'], vc['male_min_copies'])
    cnt = tab.pattern.value_counts().to_dict()
    text += [f"species with >= {vc['min_per_sex']} sexed males and females: {len(tab)}; " + '; '.join(f'{k}: {v}' for k, v in cnt.items()), '',
             tab.to_string(index=False, float_format=lambda x: f'{x:.2f}'), '', f'{len(disc)} samples are exceptions in species that look male-specific (a female with the male-specific copy is also what a B-carrying XY female would show)']
    with _out(proj, 'gsdf', 'validate', out_dir, dict(cfg=cfg, calls=bool(calls_path))) as run:
        tab.to_csv(run.dir + 'sex_by_species.tsv', sep='\t', index=False)
        disc.to_csv(run.dir + 'sex_exceptions.tsv', sep='\t', index=False)
        if calls_path:
            mism.to_csv(run.dir + 'mismatches_vs_sheet.tsv', sep='\t', index=False)
        open(run.dir + 'summary.txt', 'w').write('\n'.join(text) + '\n')
    print('\n'.join(text))
    return tab
