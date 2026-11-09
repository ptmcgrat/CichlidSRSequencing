"""B chromosome presence from the stored measurements (measure.py), with the flags of the method (config/bchrompanel.yaml). Cut-offs are PROVISIONAL (the middle of the carrier / B-free gap of the discovery)."""
import argparse
import numpy as np
import pandas as pd
import yaml
from . import measure as M


def select_blocks(blocks, cfg):
    sel = blocks.tier.isin(cfg['blocks']['tiers'])
    return sel & blocks.clean.astype(bool) if cfg['blocks'].get('clean_only') else sel


def metrics(raw, cfg, blocks):
    bd = np.asarray(raw['block_depth'], dtype=float); sb = np.asarray(raw['scale_bins'], dtype=float)
    nz = sb[sb > 0]; scale = float(np.median(nz)) if len(nz) else float('nan')
    v = bd[select_blocks(blocks, cfg).to_numpy()] / scale if scale == scale and scale > 0 else np.array([])
    v = v[~np.isnan(v)]
    return dict(scale_depth=scale, n_scale_bins_nonzero=int(len(nz)), n_blocks_used=int(len(v)), score=float(np.median(v)) if len(v) else float('nan'), mean_scaled=float(v.mean()) if len(v) else float('nan'),
                p10_scaled=float(np.quantile(v, 0.1)) if len(v) else float('nan'), frac_ge2=float((v >= 2).mean()) if len(v) else float('nan'), frac_ge3=float((v >= 3).mean()) if len(v) else float('nan'))


def call_one(m, cfg):
    th, lab = cfg['thresholds'], cfg['labels']; flags = []
    if not m['scale_depth'] >= th['min_scale_depth']:
        flags.append('low_depth')
    if m['n_scale_bins_nonzero'] < th['min_scale_bins']:
        flags.append('few_scale_bins')
    if m['n_blocks_used'] < th['min_blocks_used']:
        flags.append('few_blocks')
    s = m['score']
    if s != s:
        return dict(status='', flags=flags + ['no_estimate'])
    status = 'present' if s >= th['present_min'] else ('absent' if s <= th['absent_max'] else 'uncertain')
    if (status == 'absent' and m['frac_ge2'] > th['frac_ge2_absent_max']) or (status == 'present' and m['frac_ge2'] < th['frac_ge2_present_min']):
        flags.append('fraction_disagrees')
    return dict(status=lab[status], flags=flags)


def call_table(rows, cfg):
    """rows: dicts with SampleID, definition, block_depth, scale_bins (as stored by measure.py / collected by panel.collect). Refuses measurements made under another definition than the one in the config."""
    key = M.definition_key(cfg); blocks, _, _ = M.load_definition(cfg)
    bad = sorted({r['definition'] for r in rows} - {key})
    if bad:
        raise ValueError(f'measurements made under definition(s) {bad}, the current definition is {key}: re-measure, or rebuild the definition that matches them')
    out = []
    for raw in rows:
        if len(raw['block_depth']) != len(blocks):
            raise ValueError(f"{raw['SampleID']}: {len(raw['block_depth'])} block depths, the definition has {len(blocks)} blocks")
        m = metrics(raw, cfg, blocks); c = call_one(m, cfg)
        out.append(dict(SampleID=raw['SampleID'], bam=raw.get('bam', ''), bam_size_bytes=raw.get('bam_size_bytes'), definition=raw['definition'], **m, bchrom=c['status'], call_flags=';'.join(c['flags'])))
    return pd.DataFrame(out)


def cluster_check(calls):
    """Scores of the unflagged samples by call: should form two groups far apart (B- near 1, B+ above 4)."""
    ok = calls[calls.call_flags == '']; rows = []
    for s in sorted(ok.bchrom.unique()):
        g = ok[ok.bchrom == s]
        rows.append(dict(bchrom=s, n=len(g), median_score=float(g.score.median()), min_score=float(g.score.min()), max_score=float(g.score.max())))
    return pd.DataFrame(rows)


def main(argv=None):
    from cichlid_common.common_light import config_path
    from sd_genotypes.panel import measure as pm
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.bchrom.score'); ap.add_argument('--panel', required=True, help='folder of *.panel.json, or one file'); ap.add_argument('--out', required=True)
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(config_path('bchrompanel.yaml')))
    rows, skipped = pm.collect(a.panel, 'bchrom')
    if skipped:
        print(f'[bchrom] {len(skipped)} panel files without a bchrom measurement skipped')
    t = call_table(rows, cfg); t.to_csv(a.out, sep='\t', index=False)
    print(t.bchrom.value_counts().to_string()); print(t[t.call_flags != ''].call_flags.value_counts().to_string() if (t.call_flags != '').any() else 'no flags')
    print(cluster_check(t).to_string(index=False)); print(f'-> {a.out}')


if __name__ == '__main__':
    main()
