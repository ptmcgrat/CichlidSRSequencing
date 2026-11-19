"""Recombination from the hard calls: Marey maps, rates in windows, a region against its chromosome and against the genome, the hottest window of every chromosome.
A crossover is a change of state between two adjacent 100 kb bins in a fish; its position is the boundary between the bins. cM = 100 x crossovers per fish. The tips of a chromosome (the first and last
`terminal_bp`) and masked bins are left out of every comparison: they are assembly artifacts rather than recombination (see the README)."""
import numpy as np
import pandas as pd


class Steps:
    """Per chromosome: D (bins-1 x fish, True where the state changes), loc (boundary position, bp) and the length."""

    def __init__(self, calls):
        self.calls = calls
        self.D, self.loc, self.len, self.lg = {}, {}, {}, {}
        for c in calls.chroms:
            ix = calls.idx(c)
            Xc = calls.X[ix]
            self.D[c] = Xc[1:] != Xc[:-1]
            self.loc[c] = calls.pos(c)[1:]                                  # the boundary between two bins = the start of the second (bins are labelled by their start)
            self.len[c] = calls.length(c)
            self.lg[c] = calls.lg_of(c)
        self.n = calls.X.shape[1]

    def counts(self, c, a, b):
        """Crossovers per fish in [a, b) bp of chromosome c (float array, one per fish)."""
        m = (self.loc[c] >= a) & (self.loc[c] < b)
        return self.D[c][m].sum(0).astype(float)


def marey_bins(steps):
    rows = []
    for c in steps.calls.chroms:
        step = steps.D[c].mean(1) * 100
        rows.append(pd.DataFrame({'chrom': c, 'LG': steps.lg[c], 'Mb': steps.loc[c] / 1e6, 'cM_step': step, 'cM': np.cumsum(step)}))
    return pd.concat(rows, ignore_index=True)


def masked_windows(c, a, b, masked, bin_bp=100_000):
    return any((c, k) in masked for k in range(int(a // bin_bp), int(b // bin_bp) + 1))


def map_lengths(steps, terminal_bp=2_000_000, masked=()):
    """Per chromosome: length, total cM, cM without the tips, cM per Mb, and the number of masked bins."""
    rows = []
    for c in steps.calls.chroms:
        L = steps.len[c]
        total = 100 * steps.D[c].sum(0).mean()
        inner = 100 * steps.counts(c, terminal_bp, L - terminal_bp).mean() if L > 2 * terminal_bp else float('nan')
        rows.append(dict(chrom=c, LG=steps.lg[c], Mb=round(L / 1e6, 2), cM=round(total, 1), cM_without_tips=round(inner, 1), cM_per_Mb_interior=round(inner / ((L - 2 * terminal_bp) / 1e6), 2) if L > 2 * terminal_bp else float('nan'),
                         masked_bins=sum(1 for (cc, _) in masked if cc == c)))
    return pd.DataFrame(rows)


def window_rates(steps, width_bp=2_500_000, step_bp=None, terminal_bp=2_000_000, masked=(), skip=()):
    """Rate (cM per Mb) of every window of the given width, tips and masked windows left out; `skip` = chromosomes (accession or LG label) left out entirely."""
    step_bp = step_bp or width_bp
    rows = []
    for c in steps.calls.chroms:
        if c in skip or steps.lg[c] in skip:
            continue
        L = steps.len[c]
        for a in range(int(terminal_bp), int(L - terminal_bp - width_bp) + 1, int(step_bp)):
            if masked_windows(c, a, a + width_bp, masked, steps.calls.bin_bp):
                continue
            rows.append(dict(chrom=c, LG=steps.lg[c], start_Mb=a / 1e6, rate=100 * steps.counts(c, a, a + width_bp).mean() / (width_bp / 1e6), dist_end_Mb=min(a, L - a - width_bp) / 1e6))
    return pd.DataFrame(rows)


def best_windows(steps, width_bp=2_500_000, step_bp=500_000, terminal_bp=2_000_000, masked=(), skip=()):
    """The hottest window of each chromosome, ranked; x_chr_mean = its rate over the chromosome's own mean window rate."""
    w = window_rates(steps, width_bp, step_bp, terminal_bp, masked, skip)
    if not len(w):
        return w
    mean_c = w.groupby('chrom').rate.mean()
    best = w.loc[w.groupby('chrom').rate.idxmax()].sort_values('rate', ascending=False).reset_index(drop=True)
    best['rank'] = best.index + 1
    best['x_chr_mean'] = (best.rate / best.chrom.map(mean_c)).round(2)
    return best


def region_test(steps, chrom, a, b, terminal_bp=2_000_000, n_boot=5000, seed=1, masked=(), null_windows=None):
    """One region [a, b) bp of `chrom`: its rate, the rate of the rest of the chromosome (tips and the region left out), their ratio with a 95% bootstrap interval over FISH, P(ratio <= 1) in the bootstrap, and
    -- when `null_windows` (window_rates of the other chromosomes at the same width) is given -- how many of those windows are as hot or hotter."""
    k = steps.counts(chrom, a, b)
    L = steps.len[chrom]
    loc = steps.loc[chrom]
    bg = (loc >= terminal_bp) & (loc < L - terminal_bp) & ~((loc >= a) & (loc < b))
    kb = steps.D[chrom][bg].sum(0).astype(float)
    mb_bg = bg.sum() * steps.calls.bin_bp / 1e6
    mb_reg = (b - a) / 1e6
    rate, rate_bg = 100 * k.mean() / mb_reg, 100 * kb.mean() / mb_bg
    boot = np.random.default_rng(seed).integers(0, steps.n, (n_boot, steps.n))
    ratios = (k[boot].mean(1) / mb_reg) / np.maximum(kb[boot].mean(1) / mb_bg, 1e-12)
    out = dict(chrom=chrom, LG=steps.lg[chrom], start_Mb=a / 1e6, end_Mb=b / 1e6, rate_cM_per_Mb=rate, background_cM_per_Mb=rate_bg, ratio=rate / rate_bg if rate_bg else float('nan'),
               ci_lo=float(np.percentile(ratios, 2.5)), ci_hi=float(np.percentile(ratios, 97.5)), p_bootstrap_ratio_le_1=float((ratios <= 1).mean()), fish_with_a_crossover=int((k >= 1).sum()), crossovers=int(k.sum()))
    if null_windows is not None and len(null_windows):
        n_ge = int((null_windows.rate.to_numpy() >= rate).sum())
        out.update(null_windows=len(null_windows), null_as_hot_or_hotter=n_ge, p_genome=(1 + n_ge) / (1 + len(null_windows)), null_p95=float(null_windows.rate.quantile(0.95)), null_p99=float(null_windows.rate.quantile(0.99)), null_max=float(null_windows.rate.max()))
    return out


def rate_profile(steps, chrom, width_bp=1_000_000, n_boot=2000, seed=2):
    """Rate in non-overlapping windows along a chromosome, with a 95% bootstrap interval over fish."""
    rng = np.random.default_rng(seed)
    boot = rng.integers(0, steps.n, (n_boot, steps.n))
    rows = []
    for a in range(0, steps.len[chrom] - int(width_bp) + 1, int(width_bp)):
        k = steps.counts(chrom, a, a + width_bp)
        b = 100 * k[boot].mean(1) / (width_bp / 1e6)
        rows.append(dict(chrom=chrom, start_Mb=a / 1e6, rate=100 * k.mean() / (width_bp / 1e6), ci_lo=float(np.percentile(b, 2.5)), ci_hi=float(np.percentile(b, 97.5))))
    return pd.DataFrame(rows)


def by_group(steps, chrom, a, b, groups):
    """Crossovers per fish in the region summarised by a grouping (sire, brood, dam): fish, crossovers per fish, rate, fish with a crossover."""
    k = steps.counts(chrom, a, b)
    g = pd.Series(groups.to_numpy(), index=range(len(k)))
    rows = []
    for v in sorted(pd.unique(g.dropna())):
        m = (g == v).to_numpy()
        rows.append(dict(group=v, fish=int(m.sum()), crossovers_per_fish=float(k[m].mean()), rate_cM_per_Mb=float(100 * k[m].mean() / ((b - a) / 1e6)), fish_with_a_crossover=int((k[m] >= 1).sum())))
    return pd.DataFrame(rows)


def windows_table(steps, width_bp=2_500_000, terminal_bp=2_000_000, masked=()):
    """Every window along every chromosome with its rate, flagged tip (inside the first or last terminal_bp) and masked."""
    rows = []
    for c in steps.calls.chroms:
        L = steps.len[c]
        for a in range(0, int(L - width_bp) + 1, int(width_bp)):
            rows.append(dict(chrom=c, LG=steps.lg[c], start_Mb=a / 1e6, rate=100 * steps.counts(c, a, a + width_bp).mean() / (width_bp / 1e6), tip=bool(a < terminal_bp or a + width_bp > L - terminal_bp),
                             masked=masked_windows(c, a, a + width_bp, masked, steps.calls.bin_bp)))
    return pd.DataFrame(rows)
