"""Fisher exact tests of marker x sex within each clade (or species): where is there evidence that an allele or inversion acts as a sex determiner?
For every marker the table is  carrier / non-carrier  x  expected-sex / other-sex, where the expected sex is the sex that should carry the allele
if it is a sex determiner (male for Y-linked alleles, female for a W). Output per marker x group: counts, odds ratio, two-sided Fisher p, BH q,
concordance, and a STATUS that separates 'no variation' and 'too few individuals' from real tests, and flags tests that could not reach significance
even with perfect association ('underpowered'). Pure functions so they can be tested and run offline."""
import operator, os
import numpy as np
import pandas as pd
import yaml
from scipy import stats
from cichlid_common.common_light import config_path

OPS = {'>=': operator.ge, '<=': operator.le, '==': operator.eq, '!=': operator.ne, '>': operator.gt, '<': operator.lt}


def load_config(path=None):
    return yaml.safe_load(open(path or config_path('sex_association.yaml')))


# ------------------------------------------------------------------ genotype parsing and groups
def to_number(series, coding):
    """'number' -> float; 'gt' -> alt-allele dosage from strings like 0/0, 0/1, 1/1 (missing -> NaN)."""
    if coding == 'number':
        return pd.to_numeric(series, errors='coerce')
    def dose(x):
        if pd.isna(x):
            return np.nan
        s = str(x).replace('|', '/').strip()
        parts = s.split('/')
        if len(parts) != 2 or any(p in ('.', '') for p in parts):
            return np.nan
        try:
            return float(sum(int(p) > 0 for p in parts))
        except ValueError:
            return np.nan
    return series.map(dose)


def cond_mask(values, conds):
    m = values.notna()
    for c in conds:
        for op in ('>=', '<=', '==', '!=', '>', '<'):
            if c.startswith(op):
                m &= OPS[op](values, float(c[len(op):])); break
        else:
            raise ValueError(f'bad condition {c!r}')
    return m


def assign_groups(df, cfg):
    """Canonical clade per sample from Subgroup (preferred) or Ecogroup; unmapped values become NaN (skipped and counted)."""
    sub = df['Subgroup'].map(cfg['subgroup_map']) if 'Subgroup' in df else pd.Series(np.nan, index=df.index)
    eco = df['Ecogroup'].map(cfg['ecogroup_map']) if 'Ecogroup' in df else pd.Series(np.nan, index=df.index)
    return sub.fillna(eco)


def filter_samples(df, cfg, marker_cfg=None):
    """Returns (filtered df, flow table of how many samples were dropped and why)."""
    flow = []
    d = df.copy()
    def step(name, keep):
        nonlocal d
        flow.append((name, int((~keep).sum()))); d = d[keep]
    ok_sex = ['M', 'F'] + (['M?', 'F?'] if cfg.get('include_uncertain_sex') else [])
    step('sex not M/F', d['Sex'].isin(ok_sex))
    d = d.assign(Sex=d['Sex'].str.replace('?', '', regex=False))
    if cfg.get('exclude_offspring_with_brood') and 'BroodID' in d:
        step('lab-reared offspring (BroodID set)', d['BroodID'].isna())
    if marker_cfg is not None and str(marker_cfg['column']).startswith('Inversion') and 'InversionSource' in d:
        pat = '|'.join(cfg.get('exclude_source_patterns', []))
        if pat:
            step('inversion genotype not observed (pedigree/BC inferred)', ~d['InversionSource'].fillna('').str.contains(pat))
    return d, flow


# ------------------------------------------------------------------ the test
def min_possible_p(n_car, n_non, n_exp, n_oth):
    """Smallest two-sided Fisher p reachable with these margins (perfect association as far as the margins allow)."""
    a = min(n_car, n_exp); b = n_car - a; c = n_exp - a; d = n_non - c
    if min(a, b, c, d) < 0:
        return np.nan
    return stats.fisher_exact([[a, b], [c, d]])[1]


def test_group(g, marker, mcfg, cfg):
    """g: samples of ONE group with columns Sex and the marker column. Returns a result dict."""
    val = to_number(g[mcfg['column']], mcfg['coding'])
    car = cond_mask(val, mcfg['carrier']); non = cond_mask(val, mcfg['noncarrier'])
    used = car | non
    sex = g['Sex'][used]; exp = mcfg['carrier_sex']; oth = 'F' if exp == 'M' else 'M'
    n_car, n_non = int(car.sum()), int(non.sum())
    a = int((car & (g['Sex'] == exp)).sum()); b = int((car & (g['Sex'] == oth)).sum())
    c = int((non & (g['Sex'] == exp)).sum()); d = int((non & (g['Sex'] == oth)).sum())
    n_M = int((sex == 'M').sum()); n_F = int((sex == 'F').sum())
    res = dict(marker=marker, label=mcfg['label'], n_used=int(used.sum()), n_M=n_M, n_F=n_F, n_carriers=n_car, n_noncarriers=n_non,
               carrier_sex_expected=exp, car_expected=a, car_other=b, non_expected=c, non_other=d,
               concordance=(a + d) / max(a + b + c + d, 1), odds_ratio=np.nan, p=np.nan, min_possible_p=np.nan)
    if n_car == 0 or n_non == 0:
        res['status'] = 'not testable: no variation'; return res
    if min(n_M, n_F) < cfg['min_per_sex'] or min(n_car, n_non) < cfg['min_per_class']:
        res['status'] = 'not testable: too few individuals'; return res
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    res.update(odds_ratio=orr, p=p, min_possible_p=min_possible_p(n_car, n_non, a + c, b + d))
    res['status'] = 'underpowered' if res['min_possible_p'] > cfg['fdr'] else 'tested'
    return res


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    ranked = p[o] * n / (np.arange(n) + 1)
    q[o] = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.minimum(q, 1.0)


def run_tests(df, cfg):
    """All markers x groups. Returns (results DataFrame, sample-flow DataFrame)."""
    df = df.copy()
    df['group'] = df['Species'] if cfg['group_by'] == 'species' else assign_groups(df, cfg)
    if cfg.get('split_lab_wild') and 'LabReared' in df:                       # lab cohorts have unknown relatedness: report wild and lab separately
        lab = pd.to_numeric(df['LabReared'], errors='coerce') == 1
        df['group'] = df['group'].where(df['group'].isna(), df['group'] + np.where(lab, ' (lab)', ' (wild)'))
    rows, flows = [], []
    for marker, mcfg in cfg['markers'].items():
        if mcfg['column'] not in df.columns:
            flows.append((marker, 'column missing in table: ' + mcfg['column'], len(df)))
            continue
        d, flow = filter_samples(df, cfg, mcfg)
        flows += [(marker, n, k) for n, k in flow]
        n_nogroup = int(d['group'].isna().sum()); flows.append((marker, 'clade not assigned', n_nogroup)); d = d[d['group'].notna()]
        for grp, g in d.groupby('group'):
            r = test_group(g, marker, mcfg, cfg); r['group'] = grp; rows.append(r)
    res = pd.DataFrame(rows)
    flow = pd.DataFrame(flows, columns=['marker', 'dropped_because', 'n_samples'])
    if len(res):
        t = res.status == 'tested'
        res['q'] = np.nan
        if t.any():
            res.loc[t, 'q'] = bh(res.loc[t, 'p'].values)
        def ev(r):
            if r.status != 'tested': return ''
            if r.q < cfg['fdr']: return 'supports sex determiner' if r.odds_ratio > 1 else 'opposite direction'
            return 'no evidence (nominal p<0.05)' if r.p < 0.05 else 'no evidence'
        res['evidence'] = res.apply(ev, axis=1)
        order = {g: i for i, g in enumerate(group_rows(cfg))}
        res['_o'] = res['group'].map(order).fillna(99)
        res = res.sort_values(['marker', '_o', 'group']).drop(columns='_o').reset_index(drop=True)
    return res, flow


# ------------------------------------------------------------------ summary matrix + figure
def group_rows(cfg):
    base = cfg['group_order']
    return [f'{g}{s}' for g in base for s in ((' (wild)', ' (lab)') if cfg.get('split_lab_wild') and cfg['group_by'] == 'clade' else ('',))]


def summary_matrix(res, cfg):
    cols = list(cfg['markers']); rows = [g for g in group_rows(cfg) if g in set(res.group)] + sorted(set(res.group) - set(group_rows(cfg)))
    M = pd.DataFrame('', index=rows, columns=cols)
    for r in res.itertuples():
        M.loc[r.group, r.marker] = r.evidence if r.status == 'tested' else r.status
    return M


def plot_matrix(res, cfg, path_base):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.size': 6.5})
    cols = list(cfg['markers']); rows = [g for g in group_rows(cfg) if g in set(res.group)] + sorted(set(res.group) - set(group_rows(cfg)))
    fig, ax = plt.subplots(figsize=(1.2 + 0.78 * len(cols), 0.9 + 0.3 * len(rows)))
    ax.set_xlim(0, len(cols)); ax.set_ylim(len(rows), 0); ax.axis('off')
    for j, c in enumerate(cols):
        ax.text(j + 0.5, -0.25, cfg['markers'][c]['label'].replace(' (', '\n('), ha='center', va='bottom', fontsize=6)
        if c not in set(res.marker):
            ax.text(j + 0.5, 0.45, 'column not\nin table', ha='center', va='top', fontsize=4.8, color='#aa4444')
    for i, g in enumerate(rows):
        ax.text(-0.05, i + 0.5, g, ha='right', va='center', fontsize=6.5)
        for j, c in enumerate(cols):
            r = res[(res.group == g) & (res.marker == c)]
            fc, txt, tc = '#f4f4f4', '', '#999999'
            if len(r):
                r = r.iloc[0]
                if r.status == 'tested':
                    s = min(-np.log10(max(r.q, 1e-12)), 6) / 6
                    if r.evidence.startswith('supports'): fc, tc = plt.cm.Oranges(0.35 + 0.6 * s), 'black'
                    elif r.evidence.startswith('opposite'): fc, tc = plt.cm.Blues(0.35 + 0.6 * s), 'black'
                    else: fc, tc = '#dcdcdc', '#444444'
                    txt = f'{int(r.car_expected)}/{int(r.car_other)}|{int(r.non_expected)}/{int(r.non_other)}'
                elif r.status == 'underpowered': fc, tc, txt = '#ececec', '#777777', f'n={int(r.n_used)} (low power)'
                else: txt = 'n/a'
            ax.add_patch(Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, fc=fc, ec='white'))
            ax.text(j + 0.5, i + 0.5, txt, ha='center', va='center', fontsize=4.6, color=tc)
    ax.text(0, len(rows) + 0.55, 'cell text: carriers expected-sex/other-sex | non-carriers expected-sex/other-sex.  Orange: supports a sex determiner (FDR); blue: opposite direction; grey: no evidence; n/a: not testable.',
            fontsize=4.8, color='#666666', va='top')
    for ext in ('svg', 'png'):
        fig.savefig(f'{path_base}.{ext}', dpi=250, bbox_inches='tight')
    plt.close(fig)


def write_outputs(df, cfg, outdir):
    os.makedirs(outdir, exist_ok=True)
    res, flow = run_tests(df, cfg)
    res.to_csv(f'{outdir}/fisher_results.tsv', sep='\t', index=False)
    flow.to_csv(f'{outdir}/sample_flow.tsv', sep='\t', index=False)
    if len(res):
        summary_matrix(res, cfg).to_csv(f'{outdir}/summary_matrix.tsv', sep='\t')
        if cfg['group_by'] == 'clade':
            plot_matrix(res, cfg, f'{outdir}/sex_association_matrix')
    return res, flow


# ------------------------------------------------------------------ CLI steps
def fisher_offline(table_path, outdir, cfg_path=None, group_by=None):
    cfg = load_config(cfg_path)
    if group_by: cfg['group_by'] = group_by
    df = pd.read_csv(table_path, sep='\t')
    res, flow = write_outputs(df, cfg, outdir)
    print_summary(res)
    return res


def fisher_step(proj, group_by=None):
    cfg = load_config()
    if group_by: cfg['group_by'] = group_by
    with proj.run('sex_association', 'fisher', params=dict(cfg=cfg)) as run:
        proj.fm.readSampleDatabase()
        df = proj.fm.sample_dt.copy()
        df.to_csv(run.dir + 'sample_table_snapshot.tsv', sep='\t', index=False)       # exact input stored with the result
        res, flow = write_outputs(df, cfg, run.dir)
        print_summary(res)


def print_summary(res):
    if not len(res):
        print('no results (check the marker column names in config/sex_association.yaml)'); return
    print(res.groupby(['marker', 'status']).size().unstack(fill_value=0).to_string())
    hit = res[res.evidence == 'supports sex determiner']
    print('\nEvidence for a sex determiner (FDR):'); print(hit[['marker', 'group', 'car_expected', 'car_other', 'non_expected', 'non_other', 'odds_ratio', 'p', 'q']].to_string(index=False) if len(hit) else '  none')
