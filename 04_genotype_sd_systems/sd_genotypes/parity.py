"""Parity of the panel's gsdf and transposon calls with sexloci's `sv_genotypes.tsv` on the SAME samples: the check that the vendored code and the rules reproduce the reference implementation
(a port check, not a validation of the method: use independent labels for that).
  python -m sd_genotypes.parity --sexloci sv_genotypes.tsv [--te-calls te_calls.tsv] [--gsdf-calls gsdf_calls.tsv] [--out-dir DIR]"""
import argparse, os
import numpy as np
import pandas as pd

GTS = ('0/0', '0/1', '1/1')


def _read(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False, na_values=[''])


def compare_te(ours, theirs, min_gq=10.0):
    """Per site (and combined group) present in both tables: genotype agreement among samples both call confidently (GQ >= min_gq), and identity of the REF/ALT pair counts."""
    m = ours.merge(theirs, on='SampleID', suffixes=('', '_sx'))
    rows, mism = [], []
    sites = sorted({c[:-3] for c in ours.columns if c.endswith('_GT') and c in theirs.columns})
    for n in sites:
        ok_o = m[f'{n}_GT'].isin(GTS) & (pd.to_numeric(m[f'{n}_GQ'], errors='coerce') >= min_gq)
        ok_t = m[f'{n}_GT_sx'].isin(GTS) & (pd.to_numeric(m[f'{n}_GQ_sx'], errors='coerce') >= min_gq)
        both = m[ok_o & ok_t]; agree = both[f'{n}_GT'] == both[f'{n}_GT_sx']
        has = m[(m[f'{n}_status'] == 'OK') & (m[f'{n}_status_sx'] == 'OK')]
        same = (pd.to_numeric(has[f'{n}_ref'], errors='coerce') == pd.to_numeric(has[f'{n}_ref_sx'], errors='coerce')) & (pd.to_numeric(has[f'{n}_alt'], errors='coerce') == pd.to_numeric(has[f'{n}_alt_sx'], errors='coerce'))
        dmax = float(max((pd.to_numeric(has[f'{n}_ref'], errors='coerce') - pd.to_numeric(has[f'{n}_ref_sx'], errors='coerce')).abs().max(),
                         (pd.to_numeric(has[f'{n}_alt'], errors='coerce') - pd.to_numeric(has[f'{n}_alt_sx'], errors='coerce')).abs().max())) if len(has) else np.nan
        rows.append(dict(site=n, n_both_measured=len(has), n_counts_identical=int(same.sum()), max_count_difference=dmax, n_both_confident=len(both), n_genotype_agree=int(agree.sum())))
        for r in both[~agree].itertuples():
            mism.append(dict(kind='genotype', site=n, SampleID=r.SampleID, ours=getattr(r, f'{n}_GT'), sexloci=getattr(r, f'{n}_GT_sx')))
        for r in has[~same].itertuples():
            mism.append(dict(kind='counts', site=n, SampleID=r.SampleID, ours=f"{getattr(r, f'{n}_ref')}:{getattr(r, f'{n}_alt')}", sexloci=f"{getattr(r, f'{n}_ref_sx')}:{getattr(r, f'{n}_alt_sx')}"))
    return pd.DataFrame(rows), pd.DataFrame(mism, columns=['kind', 'site', 'SampleID', 'ours', 'sexloci'])


def compare_gsdf(ours, theirs):
    """gsdf: depth ratio, flank depth and junction ratio (differences) and copy number (agreement among samples unflagged in both)."""
    m = ours.merge(theirs, on='SampleID')
    m = m[m.ratio.notna() & pd.to_numeric(m.DepthRatio, errors='coerce').notna()].copy()
    for c in ('DepthRatio', 'FlankDepth', 'JunctionJ', 'CopyNumber'):
        m[c] = pd.to_numeric(m[c], errors='coerce')
    dr = (m.ratio - m.DepthRatio).abs(); fl = (m.flank_depth - m.FlankDepth).abs(); jj = (m.junction_ratio - m.JunctionJ).abs()
    un = m[(m.call_flags == '') & (m.Flag.fillna('') == '')]
    agree = un.gsdf_number == un.CopyNumber
    summ = dict(n_compared=len(m), max_ratio_difference=float(dr.max()) if len(m) else np.nan, median_ratio_difference=float(dr.median()) if len(m) else np.nan,
                max_flank_depth_difference=float(fl.max()) if len(m) else np.nan, max_junction_ratio_difference=float(jj.max()) if len(m) else np.nan,
                n_unflagged_in_both=len(un), n_copy_number_agree=int(agree.sum()), n_flagged_in_one_only=int(((m.call_flags == '') != (m.Flag.fillna('') == '')).sum()))
    mism = un[~agree][['SampleID', 'gsdf_number', 'CopyNumber', 'ratio', 'DepthRatio']].rename(columns={'gsdf_number': 'ours', 'CopyNumber': 'sexloci', 'ratio': 'ratio_ours', 'DepthRatio': 'ratio_sexloci'})
    return summ, mism


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.parity'); ap.add_argument('--sexloci', required=True); ap.add_argument('--te-calls'); ap.add_argument('--gsdf-calls'); ap.add_argument('--out-dir', default='.')
    ap.add_argument('--min-gq', type=float, default=10.0)
    a = ap.parse_args(argv); sx = _read(a.sexloci); os.makedirs(a.out_dir, exist_ok=True); text = []
    if a.te_calls:
        s, mm = compare_te(_read(a.te_calls), sx, a.min_gq); text += ['transposon insertions (pair counts and genotypes):', s.to_string(index=False), f'{len(mm)} differences']; mm.to_csv(os.path.join(a.out_dir, 'parity_te_differences.tsv'), sep='\t', index=False)
    if a.gsdf_calls:
        g = _read(a.gsdf_calls); g['call_flags'] = g['call_flags'].fillna('')
        s, mm = compare_gsdf(g, sx); text += ['gsdf:'] + [f'  {k}: {v}' for k, v in s.items()]; mm.to_csv(os.path.join(a.out_dir, 'parity_gsdf_differences.tsv'), sep='\t', index=False)
    print('\n'.join(text)); open(os.path.join(a.out_dir, 'parity_summary.txt'), 'w').write('\n'.join(text) + '\n')


if __name__ == '__main__':
    main()
