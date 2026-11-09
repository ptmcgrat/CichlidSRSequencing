"""Insertion genotypes from the stored pair counts (measure.py): per site, and one combined call for sites that are the same locus at several places of the reference (chr7ins at both gsdf copies).
Column names follow sexloci's sv_genotypes.tsv (<site>_GT, _ref, _alt, _GQ, _reads, _status) so that the two can be compared and fed to the same readers."""
import argparse
import numpy as np
import pandas as pd
import yaml
from . import measure as M
from .longindel import call_from_counts

GTS = ('0/0', '0/1', '1/1')


def _call(n_ref, n_alt, min_pairs):
    gt, gq, _ = call_from_counts(int(n_ref), int(n_alt), min_informative_pairs=min_pairs)
    return gt, (round(gq, 1) if gq == gq else np.nan)


def call_table(rows, cfg, copies=None):
    """rows: dicts with SampleID, definition and sites (as stored by measure.py, collected by panel.collect). copies: optional Series SampleID -> total gsdf copies (unflagged gsdf calls), used for
    the dilution flag of chr7ins. Returns one row per sample. Refuses measurements made under another definition than the current one."""
    key = M.definition_key(cfg); d = M.load_definition(cfg); call = cfg['call']
    bad = sorted({r['definition'] for r in rows} - {key})
    if bad:
        raise ValueError(f'measurements made under definition(s) {bad}, the current definition is {key}: re-measure, or rebuild the definition that matches them')
    dup = pd.Series([r['SampleID'] for r in rows]); dup = dup[dup.duplicated()].unique().tolist()
    if dup:
        raise ValueError(f'sample measured more than once: {dup[:5]}')
    names = [s['name'] for s in d['sites']]; groups = d['groups']; out = []
    for r in rows:
        row = dict(SampleID=r['SampleID'], bam=r.get('bam', ''), bam_size_bytes=r.get('bam_size_bytes'), definition=r['definition'], backend=r.get('backend', ''))
        for n in names:
            m = r['sites'].get(n)
            if m is None:
                raise ValueError(f"{r['SampleID']}: no measurement for site {n}")
            gt, gq = _call(m['n_ref'], m['n_alt'], call['min_informative_pairs']) if m['status'] == 'OK' else ('', np.nan)
            row.update({f'{n}_GT': gt, f'{n}_ref': m['n_ref'], f'{n}_alt': m['n_alt'], f'{n}_equal': m['n_equal'], f'{n}_uninf': m['n_uninf'], f'{n}_GQ': gq, f'{n}_reads': m['n_reads'], f'{n}_status': m['status']})
        for g, members in groups.items():
            ms = [r['sites'][x] for x in members]; ok = any(m['status'] == 'OK' for m in ms)
            n_ref, n_alt = sum(m['n_ref'] for m in ms), sum(m['n_alt'] for m in ms)
            gt, gq = _call(n_ref, n_alt, call['min_informative_pairs']) if ok else ('', np.nan)
            row.update({f'{g}_GT': gt, f'{g}_ref': n_ref if ok else np.nan, f'{g}_alt': n_alt if ok else np.nan, f'{g}_GQ': gq, f'{g}_reads': sum(m['n_reads'] for m in ms) if ok else np.nan,
                        f'{g}_status': 'OK' if ok else next((m['status'] for m in ms if m['status'] != 'OK'), '')})
        out.append(row)
    t = pd.DataFrame(out)
    for n in names + list(groups):                                                       # the value the Sample Database gets: confident calls only
        t[f'{n}_call'] = np.where(t[f'{n}_GT'].isin(GTS) & (pd.to_numeric(t[f'{n}_GQ'], errors='coerce') >= call['min_gq']), t[f'{n}_GT'], '')
        t[f'{n}_flags'] = ''
    if copies is not None:
        cp = t.SampleID.map(copies)
        for g in groups:
            dil = (cp >= call['dilution_min_copies']) & (t[f'{g}_GT'] == '0/0') & (pd.to_numeric(t[f'{g}_alt'], errors='coerce') >= 1)
            t.loc[dil, f'{g}_flags'] = 'dilution_possible'
    return t


def summary(t, cfg):
    d = M.load_definition(cfg); lines = []
    for n in [s['name'] for s in d['sites']] + list(d['groups']):
        conf = t[t[f'{n}_call'] != '']
        lines.append(f"{n}: {len(conf)} confident calls (GQ >= {cfg['call']['min_gq']:g}) of {int((t[f'{n}_status'] == 'OK').sum())} measured: " + ', '.join(f'{k} {v}' for k, v in conf[f'{n}_call'].value_counts().sort_index().items()))
    return '\n'.join(lines)


def main(argv=None):
    from cichlid_common.common_light import config_path
    from sd_genotypes.panel import measure as pm
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.te.call'); ap.add_argument('--panel', required=True, help='folder of *.panel.json, or one file'); ap.add_argument('--out', required=True)
    ap.add_argument('--gsdf-calls', help='gsdf_calls.tsv (copy numbers, for the chr7ins dilution flag)')
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(config_path('tepanel.yaml')))
    rows, skipped = pm.collect(a.panel, 'te')
    if skipped:
        print(f'[te] {len(skipped)} panel files without a te measurement skipped')
    copies = None
    if a.gsdf_calls:
        g = pd.read_csv(a.gsdf_calls, sep='\t', keep_default_na=False, na_values=['']).fillna({'call_flags': ''})
        copies = g[g.call_flags == ''].set_index('SampleID').gsdf_number
    t = call_table(rows, cfg, copies); t.to_csv(a.out, sep='\t', index=False)
    print(summary(t, cfg)); print(f'-> {a.out}')


if __name__ == '__main__':
    main()
