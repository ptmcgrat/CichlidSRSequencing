"""B-block calling: the coverage-ratio analysis of Clark, Conte & Kocher (2018, Genes 9:610), on binned depth.

For a B+ female and a B-free baseline (male pool of the same population):
    scaled coverage = depth / genome-wide depth of the sample;   SCR = scaled(B female) / scaled(baseline)
    a bin PASSES when SCR >= scr_min (3: at least four B-located copies; a plain A-chromosome duplication gives 2) and a one-sided binomial test of the female's reads against the
    baseline's, with the expected share taken from the two genome-wide depths, has p <= p_max
    passing bins within merge_bp are merged; blocks <= min_block_bp are dropped, and so are blocks in which <= min_fraction of the bins pass  -> "B blocks"
    core blocks = positions covered by blocks of at least `min_individuals` carriers
Differences from the paper, all deliberate and all configurable (config/bchrom.yaml): bins of bin_size bp instead of single bases (mean depth of the bin stands in for the per-base
depth, so the binomial test keeps its per-base meaning); a bin needs scaled baseline coverage >= min_baseline_scaled (SCR is meaningless where the baseline has no reads); the scale
factor is the median depth (robust to collapsed repeats) unless scale_by: mean; regions of the reference known to be misassembled are not called."""
import numpy as np
import pandas as pd
from scipy.stats import binom

DEFAULTS = dict(scr_min=3.0, p_max=0.001, merge_bp=300, min_block_bp=500, min_fraction_passing=0.10, min_baseline_scaled=0.2)


def binom_upper_p(x, n, p):
    """P(X >= x) for X ~ Binomial(n, p), vectorised; 1.0 where x <= 0, and where n == 0."""
    x = np.asarray(x, dtype=np.int64); n = np.asarray(n, dtype=np.int64)
    out = np.ones(x.shape, dtype=float)
    ok = (x > 0) & (n > 0)
    out[ok] = binom.sf(x[ok] - 1, n[ok], p)
    return out


def passing_bins(depth_b, depth_m, sf_b, sf_m, p_exp, thr, mask_bins=None):
    """Boolean array: which bins pass SCR and the binomial test. Also returns the SCR array (nan where the baseline is too low)."""
    sb = depth_b / sf_b; sm = depth_m / sf_m
    valid = sm >= thr['min_baseline_scaled']
    scr = np.full(depth_b.shape, np.nan, dtype=np.float64)
    scr[valid] = sb[valid] / sm[valid]
    cand = np.where(scr >= thr['scr_min'])[0]
    ok = np.zeros(depth_b.shape, dtype=bool)
    if len(cand):
        x = np.rint(depth_b[cand]).astype(np.int64); n = x + np.rint(depth_m[cand]).astype(np.int64)
        ok[cand] = binom_upper_p(x, n, p_exp) <= thr['p_max']
    if mask_bins is not None:
        ok &= ~mask_bins
    return ok, scr


def runs_of_true(mask):
    """[(start, end_exclusive), ...] of the consecutive True runs."""
    if not mask.any():
        return []
    d = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    return list(zip(np.where(d == 1)[0].tolist(), np.where(d == -1)[0].tolist()))


def merge_runs(runs, bin_size, merge_bp):
    """Merge runs whose gap (bp) is <= merge_bp (bedtools merge -d)."""
    out = []
    for s, e in runs:
        if out and (s - out[-1][1]) * bin_size <= merge_bp:
            out[-1][1] = e
        else:
            out.append([s, e])
    return out


def est_copies(scr, baseline_scaled, female_scaled):
    """Clark's copy-number estimate of a block on the B: Eq. 3, 2*SCR - 2 (a haploid B against a diploid A, minus the two A copies) when the baseline's scaled coverage is >= 1;
    otherwise (baseline < 1: poor alignment, SCR would over-estimate) Eq. 4, 2 * the female's scaled coverage, as printed in the paper (the text there is ambiguous about
    whether the A copies are subtracted). The B-space totals are therefore only as good as Eq. 4."""
    return (2 * scr - 2) if baseline_scaled >= 1 else 2 * female_scaled


def call_blocks(depth_b, depth_m, sf_b, sf_m, p_exp, bin_size, thr=None, mask=None, lengths=None):
    """B blocks of one carrier against one baseline. depth_b / depth_m: {sequence: per-bin depth}. mask: {sequence: boolean per-bin array} of bins that must not be called.
    Returns a DataFrame: chrom, start, end, length, n_bins, frac_passing, raw_cov_B, scaled_cov_B, SCR, scaled_cov_baseline, est_copies, (all means over the block's bins)."""
    thr = {**DEFAULTS, **(thr or {})}
    rows = []
    for c, db in depth_b.items():
        dm = depth_m.get(c)
        if dm is None or len(dm) != len(db):
            continue
        ok, scr = passing_bins(db, dm, sf_b, sf_m, p_exp, thr, None if mask is None else mask.get(c))
        for s, e in merge_runs(runs_of_true(ok), bin_size, thr['merge_bp']):
            n = e - s; length = n * bin_size
            if lengths is not None:
                length = min(e * bin_size, lengths[c]) - s * bin_size
            frac = float(ok[s:e].sum()) / n
            if length <= thr['min_block_bp'] or frac <= thr['min_fraction_passing']:
                continue
            sb = float(np.mean(db[s:e])) / sf_b; sm = float(np.mean(dm[s:e])) / sf_m
            sc = float(np.nanmean(scr[s:e])) if np.isfinite(scr[s:e]).any() else np.nan
            rows.append(dict(chrom=c, start=s * bin_size, end=s * bin_size + length, length=length, n_bins=n, frac_passing=frac, raw_cov_B=float(np.mean(db[s:e])),
                             scaled_cov_B=sb, SCR=sc, scaled_cov_baseline=sm, est_copies=max(0.0, est_copies(sc, sm, sb)) if np.isfinite(sc) else np.nan))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'length', 'n_bins', 'frac_passing', 'raw_cov_B', 'scaled_cov_B', 'SCR', 'scaled_cov_baseline', 'est_copies'])


def b_space_mb(blocks):
    """(A-space Mb, B-space Mb): total block length, and total length x estimated B copies (Clark's Table 3 columns)."""
    if not len(blocks):
        return 0.0, 0.0
    a = blocks.length.sum() / 1e6
    b = float((blocks.length * blocks.est_copies.fillna(0)).sum()) / 1e6
    return float(a), b


def core_regions(blocks_by_carrier, bin_size, min_n, min_len_bp=0):
    """Regions covered by the blocks of at least min_n carriers. Returns chrom, start, end, n_carriers (max over the region), carriers (';'-joined ids that overlap it).
    Coordinates are bin-aligned, as the blocks are."""
    ids = list(blocks_by_carrier)
    chroms = sorted({c for b in blocks_by_carrier.values() for c in b.chrom.unique()})
    rows = []
    for c in chroms:
        top = 0
        for b in blocks_by_carrier.values():
            sub = b[b.chrom == c]
            if len(sub):
                top = max(top, int(sub.end.max()))
        nb = top // bin_size + 2
        cnt = np.zeros(nb, dtype=np.int32)
        for i in ids:
            sub = blocks_by_carrier[i]; sub = sub[sub.chrom == c]
            if not len(sub):
                continue
            diff = np.zeros(nb + 1, dtype=np.int32)                      # a carrier's blocks never overlap each other (they are merged), so each adds at most 1
            np.add.at(diff, (sub.start.values // bin_size).astype(int), 1)
            np.add.at(diff, (-(-sub.end.values // bin_size)).astype(int), -1)
            cnt += np.cumsum(diff)[:nb]
        for s, e in runs_of_true(cnt >= min_n):
            st, en = s * bin_size, e * bin_size
            if en - st < min_len_bp:
                continue
            who = [i for i in ids if ((blocks_by_carrier[i].chrom == c) & (blocks_by_carrier[i].start < en) & (blocks_by_carrier[i].end > st)).any()]
            rows.append(dict(chrom=c, start=st, end=en, n_carriers=int(cnt[s:e].max()), carriers=';'.join(who)))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'n_carriers', 'carriers'])


def block_signal(core, depths, sfs, bin_size):
    """Mean scaled coverage of every core block in every sample: DataFrame (one row per block, one column per sample id)."""
    out = {}
    for sid, d in depths.items():
        v = np.full(len(core), np.nan)
        for k, r in enumerate(core.itertuples()):
            a = d.get(r.chrom)
            if a is None:
                continue
            seg = a[r.start // bin_size: -(-r.end // bin_size)]
            if len(seg):
                v[k] = float(seg.mean()) / sfs[sid]
        out[sid] = v
    return pd.DataFrame(out, index=core.index)


def signal_summary(sig):
    """Per-sample statistics of a block-signal matrix: median, mean, 10th percentile, and the fractions of blocks at >= 2x and >= 3x the genome-wide depth."""
    rows = []
    for sid in sig.columns:
        v = sig[sid].dropna()
        rows.append(dict(sample=sid, n_blocks=len(v), median_scaled=float(v.median()) if len(v) else np.nan, mean_scaled=float(v.mean()) if len(v) else np.nan,
                         p10_scaled=float(v.quantile(0.1)) if len(v) else np.nan, frac_ge2=float((v >= 2).mean()) if len(v) else np.nan, frac_ge3=float((v >= 3).mean()) if len(v) else np.nan))
    return pd.DataFrame(rows)


def overlap_bp(a, b):
    """Total bp of the intervals of `a` that lie inside the union of `b` (same sequence names)."""
    tot = 0
    for c, ga in a.groupby('chrom'):
        gb = b[b.chrom == c].sort_values('start')
        if not len(gb):
            continue
        bs, be = gb.start.values, gb.end.values
        for s, e in zip(ga.start.values, ga.end.values):
            m = (bs < e) & (be > s)
            if m.any():
                tot += int(np.minimum(be[m], e).sum() - np.maximum(bs[m], s).sum())
    return tot
