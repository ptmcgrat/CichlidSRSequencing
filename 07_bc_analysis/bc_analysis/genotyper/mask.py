"""Loci whose signal is the same in many fish are a systematic artifact, not ancestry (a backcross fish is heterozygous at a given locus in about half of the fish, never in most of them at a hundred-kilobase
scale, and a real crossover hotspot does not hold the crossovers of a tenth of all meioses). A bin that holds a short interior segment in a sizeable fraction of the fish, or a breakpoint in a larger fraction,
is masked: the HMM step leaves its counts out (no recounting is needed)."""
import collections
import math

import numpy as np
import pandas as pd

COLUMNS = ['chrom', 'bin', 'n_segment_fish', 'n_breakpoint_fish', 'reason']


def threshold(fraction, n_fish):
    """The number of fish a fraction stands for, rounded up: 0.07 of 100 fish is 7 (0.07 * 100 is 7.000000000000001 and would round up to 8)."""
    return max(1, math.ceil(round(fraction * n_fish, 6)))


def derive(shorts, bps, n_fish, mp, bin_bp):
    """DataFrame (chrom, bin, n_segment_fish, n_breakpoint_fish, reason) of the bins to mask, from the short interior segments and the breakpoints of one run.
    reason: segments / breakpoints / both (the bin itself passes the threshold) or pad (a neighbour of one that does). A fish counts once per bin."""
    seg, brk = collections.Counter(), collections.Counter()
    for r in shorts.itertuples():
        for b in range(int(r.start // bin_bp), int(r.end // bin_bp)):
            seg[(r.chrom, b)] += 1
    for r in bps.itertuples():
        brk[(r.chrom, int(r.bin_hi // bin_bp))] += 1                      # the first bin of the new state
    s_min, b_min = threshold(mp['segment_fish_fraction'], n_fish), threshold(mp['breakpoint_fish_fraction'], n_fish)
    reason = {k: 'segments' for k, c in seg.items() if c >= s_min}
    for k, c in brk.items():
        if c >= b_min:
            reason[k] = 'both' if k in reason else 'breakpoints'
    for ch, b in list(reason):
        for j in range(b - int(mp['pad_bins']), b + int(mp['pad_bins']) + 1):
            if j >= 0:
                reason.setdefault((ch, j), 'pad')
    rows = [dict(chrom=ch, bin=b, n_segment_fish=seg.get((ch, b), 0), n_breakpoint_fish=brk.get((ch, b), 0), reason=w) for (ch, b), w in sorted(reason.items())]
    return pd.DataFrame(rows, columns=COLUMNS)


def intervals(df, bin_bp):
    """Contiguous runs of masked bins: dicts chrom, start, end (bp), bins, segment_fish and breakpoint_fish (the largest count in the run)."""
    out = []
    for ch, g in df.groupby('chrom', sort=False):
        bins = g.bin.to_numpy()
        start = 0
        for i in range(1, len(bins) + 1):
            if i == len(bins) or bins[i] != bins[i - 1] + 1:
                seg = g.iloc[start:i]
                out.append(dict(chrom=ch, start=int(seg.bin.iloc[0] * bin_bp), end=int((seg.bin.iloc[-1] + 1) * bin_bp), bins=len(seg),
                                segment_fish=int(seg.n_segment_fish.max()), breakpoint_fish=int(seg.n_breakpoint_fish.max())))
                start = i
    return out


def load(path, bin_bp):
    """{chrom: sorted array of masked bin indices} from a mask file written by `genotype mask`."""
    import os
    if not os.path.exists(path):
        raise SystemExit(f'the mask file {path} is missing: run `genotype mask` first')
    df = pd.read_csv(path, sep='\t', dtype={'chrom': str})
    if len(df) and int((df.end - df.start).iloc[0]) != bin_bp:
        raise SystemExit(f'the mask was made with bins of {int((df.end - df.start).iloc[0])} bp but the HMM uses {bin_bp} bp: make the mask again')
    return {ch: np.sort(g.start.to_numpy() // bin_bp) for ch, g in df.groupby('chrom')}


def apply(counts, masked, bin_bp):
    """The counts {chrom: (pos, n_mc, n_yh, yh_is_ref)} without the sites in the masked bins, and the number of sites dropped."""
    out, dropped = {}, 0
    for ch, (pos, mc, yh, ref) in counts.items():
        m = masked.get(ch)
        if m is None or not len(m):
            out[ch] = (pos, mc, yh, ref)
            continue
        keep = ~np.isin(pos // bin_bp, m)
        dropped += int((~keep).sum())
        out[ch] = (pos[keep], mc[keep], yh[keep], ref[keep])
    return out, dropped
