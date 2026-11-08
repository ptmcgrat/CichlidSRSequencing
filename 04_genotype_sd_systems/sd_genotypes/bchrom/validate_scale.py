"""Does the sampled-bin depth estimate reproduce the discovery? For every depth track of the discovery: the score with the genome-wide median the discovery used (exact: what `bchrom call` and sample_summary.tsv
report) against the score with the scale taken from the sampled bins of the definition only (what the panel measurement can do from a BAM). Run on the 21 tracks before trusting the measurement; run again on the
first pilot BAMs against a mosdepth track of the same BAM.
  python -m sd_genotypes.bchrom.validate_scale --tracks DIR [--summary sample_summary.tsv]

On a pilot BAM: the panel measurement against a mosdepth run (`mosdepth -n --by 50 PREFIX S.all.bam`) of the same BAM: block by block and in the final score.
  python -m sd_genotypes.bchrom.validate_scale --panel S.panel.json --mosdepth PREFIX.regions.bed.gz"""
import argparse, os
import numpy as np
import pandas as pd
import yaml
from . import measure as M, score as S
from .depth import load_track


def track_scores(arrays, meta, blocks, bins, cfg):
    b = meta['bin_size']; exact = float(meta['stats']['median'])
    sampled = np.array([arrays[c][s // b] for c, s in zip(bins.chrom, bins.start)]); nz = sampled[sampled > 0]
    sampled_scale = float(np.median(nz))
    bd = np.array([float(arrays[r.chrom][r.start // b: -(-r.end // b)].mean()) for r in blocks.itertuples()])
    sel = S.select_blocks(blocks, cfg).to_numpy()
    return dict(exact_scale=exact, sampled_scale=sampled_scale, exact_score=float(np.median(bd[sel] / exact)), sampled_score=float(np.median(bd[sel] / sampled_scale)))


def compare(tracks_dir, cfg, summary=None):
    blocks, bins, meta = M.load_definition(cfg); rows = []
    for f in sorted(os.listdir(tracks_dir)):
        if f.endswith('.npz'):
            run = f[:-4]; arrays, tm = load_track(os.path.join(tracks_dir, run))
            r = dict(run=run, **track_scores(arrays, tm, blocks, bins, cfg)); r['scale_ratio'] = r['sampled_scale'] / r['exact_scale']; r['score_ratio'] = r['sampled_score'] / r['exact_score']
            rows.append(r)
    t = pd.DataFrame(rows)
    if summary is not None and len(t):
        t = t.merge(summary[['run', 'tierA_median_scaled']].rename(columns={'tierA_median_scaled': 'reported_score'}), on='run', how='left')
        t['exact_vs_reported'] = t.exact_score / t.reported_score
    return t


def compare_bam(panel_json, mosdepth_regions, cfg, bin_bp=50):
    """The panel measurement of a BAM (stored JSON) against mosdepth 50 bp bins of the same BAM: scale, score and block depths. Differences in read filters or in mate-overlap handling between
    `samtools bedcov` and mosdepth would show up here as a block-depth ratio away from 1."""
    import json
    from .depth import parse_regions, depth_stats
    blocks, bins, meta = M.load_definition(cfg)
    raw = json.load(open(panel_json))['modules']['bchrom']
    if raw['definition'] != M.definition_key(cfg):
        raise ValueError(f"{panel_json}: measured under definition {raw['definition']}, the current definition is {M.definition_key(cfg)}")
    arrays, lengths = parse_regions(mosdepth_regions, bin_bp)
    stats = depth_stats(arrays, bin_bp, meta['min_chrom_bp'])
    mo = track_scores(arrays, dict(bin_size=bin_bp, stats=stats), blocks, bins, cfg)
    pm_ = S.metrics(dict(block_depth=raw['block_depth'], scale_bins=raw['scale_bins']), cfg, blocks)
    mb = np.array([float(arrays[r.chrom][r.start // bin_bp: -(-r.end // bin_bp)].mean()) for r in blocks.itertuples()])
    ratio = np.asarray(raw['block_depth'], dtype=float) / np.where(mb > 0, mb, np.nan)
    return dict(panel_scale=pm_['scale_depth'], mosdepth_scale=mo['exact_scale'], panel_score=pm_['score'], mosdepth_score=mo['exact_score'], score_ratio=pm_['score'] / mo['exact_score'],
                block_ratio_median=float(np.nanmedian(ratio)), block_ratio_within_5pct=float(np.nanmean(np.abs(ratio - 1) <= 0.05)))


def main(argv=None):
    from cichlid_common.common_light import config_path
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.bchrom.validate_scale'); ap.add_argument('--tracks'); ap.add_argument('--summary'); ap.add_argument('--max-deviation', type=float, default=0.02)
    ap.add_argument('--panel'); ap.add_argument('--mosdepth')
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(config_path('bchrompanel.yaml')))
    if a.panel or a.mosdepth:
        if not (a.panel and a.mosdepth):
            raise SystemExit('--panel and --mosdepth go together')
        r = compare_bam(a.panel, a.mosdepth, cfg)
        print('\n'.join(f'  {k:24s} {v:.4f}' for k, v in r.items()))
        ok = abs(r['score_ratio'] - 1) <= a.max_deviation
        print('OK' if ok else 'NOT OK: the panel measurement of this BAM does not reproduce the mosdepth score; look at the block ratio (read filters, mate overlap) before using the measurement')
        return
    if not a.tracks:
        raise SystemExit('give --tracks DIR, or --panel and --mosdepth')
    t = compare(a.tracks, cfg, pd.read_csv(a.summary, sep='\t') if a.summary else None)
    pd.set_option('display.width', 200); print(t.round(4).to_string(index=False))
    dev = float((t.score_ratio - 1).abs().max()); sdev = float((t.scale_ratio - 1).abs().max())
    print(f'\nlargest deviation of the sampled-scale score from the exact score: {dev:.4f}; of the scale itself: {sdev:.4f} (limit {a.max_deviation})')
    if 'exact_vs_reported' in t:
        print(f'exact score against the discovery report (tierA_median_scaled): largest deviation {float((t.exact_vs_reported - 1).abs().max()):.2e} (float32 rounding, ~1e-7, is expected)')
    print('OK' if dev <= a.max_deviation else 'NOT OK: the sampled scale does not reproduce the discovery; do not use the measurement until this is understood')


if __name__ == '__main__':
    main()
