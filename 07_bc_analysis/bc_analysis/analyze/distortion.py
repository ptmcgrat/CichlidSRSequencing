"""Transmission ratio distortion and linkage between segments.
trd_scan: at every bin the share of fish with the YH-carrying state is tested against 0.5 (z = (2k - n) / sqrt(n)). Genome-wide significance comes from a SIGN FLIP: each fish's whole chromosome is
complemented with probability 1/2, which keeps its crossovers and the linkage and removes any preference for one parent; the maximum of |z| over the bins of such a data set is the null.
linkscan: for one segment, how strongly its state in every fish is inherited together with every other bin (0 = independent, 1 = always together or always opposite): it shows which chromosome a segment belongs to."""
import numpy as np
import pandas as pd
from scipy import stats


def trd_scan(calls, n_perm=2000, seed=1):
    n = calls.X.shape[1]
    S = 2 * calls.X.astype(np.float32) - 1
    k = calls.X.sum(1)
    rng = np.random.default_rng(seed)
    maxc = np.empty((len(calls.chroms), n_perm), dtype=np.float32)
    for j, c in enumerate(calls.chroms):
        eps = rng.choice(np.array([-1, 1], dtype=np.float32), size=(n, n_perm))
        maxc[j] = (np.abs(S[calls.idx(c)] @ eps) / np.sqrt(n)).max(0)
    return dict(n=n, f=k / n, z=(2 * k - n) / np.sqrt(n), maxc=maxc, G=maxc.max(0))


def p_genomewide(sc, zabs):
    """How often the maximum |z| of a sign-flipped data set reaches zabs (with the observed data counted once)."""
    return (1 + (sc['G'] >= zabs).sum()) / (1 + len(sc['G']))


def trd_regions(calls, sc, zmin=3.0, gap=3):
    """Runs of bins with |z| >= zmin of the same sign (up to `gap` bins may interrupt a run): chrom, LG, start, end, bins, peak, share at the peak, mean share, z at the peak, genome-wide p of the peak."""
    out = []
    pos_all = calls.meta.pos.to_numpy()
    for c in calls.chroms:
        ix = calls.idx(c)
        z, f, pos = sc['z'][ix], sc['f'][ix], pos_all[ix]
        hit = np.where(np.abs(z) >= zmin)[0]
        if not len(hit):
            continue
        for g in np.split(hit, np.where(np.diff(hit) > gap + 1)[0] + 1):
            sign = np.sign(z[g].sum())
            j = g[np.argmax(z[g] * sign)]
            out.append(dict(chrom=c, LG=calls.lg_of(c), start=int(pos[g[0]]), end=int(pos[g[-1]] + calls.bin_bp), bins=len(g), peak_Mb=pos[j] / 1e6, share_at_peak=float(f[j]), mean_share=float(f[g].mean()), z_peak=float(z[j]),
                            p_genomewide=float(p_genomewide(sc, abs(z[j])))))
    cols = ['chrom', 'LG', 'start', 'end', 'bins', 'peak_Mb', 'share_at_peak', 'mean_share', 'z_peak', 'p_genomewide']
    return pd.DataFrame(out, columns=cols).sort_values('z_peak', key=abs, ascending=False).reset_index(drop=True)


def region_state(calls, chrom, start, end):
    """Each fish's state over a region (the mean of its bins, rounded): 0/1 array."""
    pos = calls.pos(chrom)
    m = (pos >= start) & (pos < end)
    return (calls.X[calls.idx(chrom)][m].mean(0) > 0.5).astype(int)


def group_checks(state, groups, label):
    """Per group (sire, brood, dam, sex): fish, fish with the YH state, share, binomial p against 0.5; and the p of a chi-square test of heterogeneity across the groups."""
    g = pd.Series(groups.to_numpy(), index=range(len(state)))
    rows = []
    for v in sorted(pd.unique(g.dropna())):
        m = (g == v).to_numpy()
        k, n = int(state[m].sum()), int(m.sum())
        rows.append(dict(by=label, group=v, fish=n, YH=k, share=k / n, p_vs_half=float(stats.binomtest(k, n, 0.5).pvalue)))
    t = pd.DataFrame(rows, columns=['by', 'group', 'fish', 'YH', 'share', 'p_vs_half'])
    het = float('nan')
    if len(t) > 1:
        tab = pd.crosstab(g.dropna().to_numpy(), state[g.notna().to_numpy()])
        if tab.shape[0] > 1 and tab.shape[1] > 1:
            het = float(stats.chi2_contingency(tab)[1])
    t['heterogeneity_p'] = het
    return t


def linkscan(calls, chrom, start, end, smooth=10):
    """(best link per chromosome, the profile of the segment's own chromosome by 5 Mb). link = |agreement - 0.5| x 2, smoothed over `smooth` bins; the best position of each chromosome is reported."""
    s = region_state(calls, chrom, start, end)
    agree = (calls.X == s[None, :]).mean(1)
    link = np.abs(agree - 0.5) * 2
    rows = []
    for c in calls.chroms:
        ix = calls.idx(c)
        sm = pd.Series(link[ix]).rolling(smooth, center=True, min_periods=1).mean().to_numpy()
        j = int(sm.argmax())
        rows.append(dict(chrom=c, LG=calls.lg_of(c), best_Mb=calls.pos(c)[j] / 1e6, link=float(sm[j]), agree=float(agree[ix][j])))
    best = pd.DataFrame(rows).sort_values('link', ascending=False).reset_index(drop=True)
    ix = calls.idx(chrom)
    own = pd.DataFrame({'Mb': (calls.pos(chrom) // 5_000_000) * 5, 'link': pd.Series(link[ix]).rolling(smooth, center=True, min_periods=1).mean().to_numpy()}).groupby('Mb').link.mean().reset_index()
    return best, own
