"""LG5 inversion genotype from the stored measurements (measure.py), with the flags of the method (config/lg5panel.yaml). Thresholds are PROVISIONAL (midpoints of the expected fractions)."""
import argparse, math, os
import numpy as np
import pandas as pd
import yaml

EXPECTED = {0: 0.0, 1: 0.5, 2: 1.0}


def metrics(raw, cfg):
    th = cfg['thresholds']
    ref = np.asarray(raw['ref'], dtype=float); alt = np.asarray(raw['alt'], dtype=float); d = ref + alt
    use = d >= th['min_site_depth']
    f = alt[use] / d[use]
    jc = raw['junction']; inv = (jc['A'] + jc['B']) / 2; tot = jc['A'] + jc['B'] + jc['C']
    return dict(n_markers=len(d), n_usable=int(use.sum()), mean_depth=float(d.mean()) if len(d) else float('nan'),
                median_site_fraction=float(np.median(f)) if len(f) else float('nan'), aggregate_fraction=float(alt[use].sum() / d[use].sum()) if use.any() else float('nan'),
                junction_A=jc['A'], junction_B=jc['B'], junction_C=jc['C'], junction_reads=tot, junction_fraction=float(inv / (inv + jc['C'])) if inv + jc['C'] > 0 else float('nan'))


def klass(x, cuts):
    return 0 if x < cuts[0] else (1 if x < cuts[1] else 2)


def call_one(m, cfg):
    th = cfg['thresholds']; flags = []
    if m['n_usable'] < th['min_usable_markers']:
        flags.append('few_markers')
    if not m['mean_depth'] >= th['min_mean_depth']:
        flags.append('low_depth')
    x = m['median_site_fraction']
    if x != x:
        return dict(n_inverted=np.nan, genotype='', flags=flags + ['no_estimate'])
    k = klass(x, th['class_cuts'])
    if any(abs(x - c) < th['near_cut'] for c in th['class_cuts']):
        flags.append('near_threshold')
    if abs(m['aggregate_fraction'] - EXPECTED[k]) > th['aggregate_max_distance']:
        flags.append('aggregate_disagrees')
    if m['junction_reads'] < th['min_junction_reads'] or m['junction_fraction'] != m['junction_fraction']:
        flags.append('junction_uninformative')
    elif abs(m['junction_fraction'] - EXPECTED[k]) > th['junction_max_distance']:
        flags.append('junction_disagrees')
    return dict(n_inverted=k, genotype=cfg['labels'][k], flags=flags)


def call_table(rows, cfg):
    """rows: dicts with SampleID, definition, ref, alt, junction (as stored by measure.py / collected by panel.collect). Refuses to mix definitions."""
    defs = {r['definition'] for r in rows}
    if len(defs) > 1:
        raise ValueError(f'measurements made under different definitions: {sorted(defs)}; re-measure or call them separately')
    out = []
    for raw in rows:
        m = metrics(raw, cfg); c = call_one(m, cfg)
        out.append(dict(SampleID=raw['SampleID'], definition=raw['definition'], **m, n_inverted=c['n_inverted'], genotype=c['genotype'], call_flags=';'.join(c['flags'])))
    return pd.DataFrame(out)


def cluster_check(calls, cfg):
    """Where the unflagged samples sit in each class against the expected fractions (the drift check for a new batch; should form three clusters)."""
    ok = calls[calls.call_flags == '']; rows = []
    for k in (0, 1, 2):
        g = ok[ok.n_inverted == k]
        rows.append(dict(n_inverted=k, n=len(g), median_site_fraction=float(g.median_site_fraction.median()) if len(g) else np.nan, expected=EXPECTED[k],
                         median_junction_fraction=float(g.junction_fraction.median()) if len(g) else np.nan))
    return pd.DataFrame(rows)


def main(argv=None):
    from cichlid_common.common_light import config_path
    from sd_genotypes.panel import measure as pm
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.lg5inv.call'); ap.add_argument('--panel', required=True, help='folder of *.panel.json, or one file'); ap.add_argument('--out', required=True)
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(config_path('lg5panel.yaml')))
    rows, skipped = pm.collect(a.panel, 'lg5inv')
    if skipped:
        print(f'[lg5inv] {len(skipped)} panel files without an lg5inv measurement skipped')
    t = call_table(rows, cfg); t.to_csv(a.out, sep='\t', index=False)
    print(t.genotype.value_counts().to_string()); print(t[t.call_flags != ''].call_flags.value_counts().to_string() if (t.call_flags != '').any() else 'no flags')
    print(cluster_check(t, cfg).to_string(index=False)); print(f'-> {a.out}')


if __name__ == '__main__':
    main()
