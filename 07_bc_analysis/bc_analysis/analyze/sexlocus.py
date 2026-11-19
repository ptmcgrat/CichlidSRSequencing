"""Where is the sex determiner: the association of every bin's state with the recorded sex (sexscan), and its limits from the recombinant fish (finemap).
A male (1) carries the YH-derived state at the determiner and a female (0) the MC state. A fish with a crossover near the locus shows on which side of the crossover the determiner must lie: the side whose state
matches the fish's sex. The limits are the nearest such crossovers on either side."""
import numpy as np
import pandas as pd
from scipy import stats


def sexscan(calls, sex, chroms=None):
    """Fisher's exact test of state against sex in every bin: columns chrom, LG, pos, M_YH, F_YH, M_MC, F_MC, nlogp."""
    ok = ~np.isnan(sex)
    sx = sex[ok].astype(int)
    rows = []
    for c in (chroms or calls.chroms):
        ix = calls.idx(c)
        for i, p in zip(ix, calls.pos(c)):
            s = calls.X[i, ok]
            a, b, d, e = int(((s == 1) & (sx == 1)).sum()), int(((s == 1) & (sx == 0)).sum()), int(((s == 0) & (sx == 1)).sum()), int(((s == 0) & (sx == 0)).sum())
            rows.append((c, calls.lg_of(c), int(p), a, b, d, e, -np.log10(max(stats.fisher_exact([[a, b], [d, e]])[1], 1e-300))))
    return pd.DataFrame(rows, columns=['chrom', 'LG', 'pos', 'M_YH', 'F_YH', 'M_MC', 'F_MC', 'nlogp'])


def scan_summary(scan, drop_log=(2, 6)):
    """Per chromosome: the best bin, its -log10 p, the counts there and the interval within 2 and 6 log units of the best."""
    rows = []
    for c, g in scan.groupby('chrom', sort=False):
        j = g.nlogp.idxmax()
        r = dict(chrom=c, LG=g.LG.iloc[0], best_Mb=g.pos[j] / 1e6, nlogp=g.nlogp[j], M_with_YH=int(g.M_YH[j]), M_with_MC=int(g.M_MC[j]), F_with_YH=int(g.F_YH[j]), F_with_MC=int(g.F_MC[j]))
        for d in drop_log:
            w = g[g.nlogp >= g.nlogp.max() - d]
            r[f'within_{d}_from_Mb'], r[f'within_{d}_to_Mb'] = w.pos.min() / 1e6, (w.pos.max() + 100_000) / 1e6
        rows.append(r)
    return pd.DataFrame(rows).sort_values('nlogp', ascending=False).reset_index(drop=True)


def finemap(calls, sex, chrom, locus_bp=None, window_bp=5_000_000, exclude=(), max_k=3):
    """The limits of the determiner on `chrom`. locus_bp: where to look (default: the middle of the run of bins with the fewest fish whose state disagrees with their sex); window_bp: how far from it a crossover
    counts as nearby; exclude: fish to ignore (their recorded sex is doubted). A fish that disagrees with its sex AT the locus is ignored too (it is mis-sexed or carries a second factor) and listed.
    For every other fish, going outward from the locus on each side, the first bin where its state stops matching its sex marks its nearest crossover: the determiner lies right of the one on the left
    ('SD right of') and left of the one on the right ('SD left of'). This holds however many other crossovers the fish has farther out. The limits are the nearest such crossovers over the fish (k_wrong_per_side
    = 0), or the k-th nearest on each side if the k nearest fish were wrong (k = 1 ... max_k).
    Returns a dict: locus_bp, discordant (fish at the locus), excluded (all ignored), recombinants (DataFrame: fish, sex, boundary_bp, bound, used), limits (DataFrame: k_wrong_per_side, left_bp, right_bp,
    width_bp) and profile (DataFrame: pos, discordant)."""
    ok = ~np.isnan(sex)
    ix = calls.idx(chrom)
    pos = calls.pos(chrom)
    Xs = calls.X[ix][:, ok]
    sx = sex[ok].astype(int)
    names = np.array(calls.fish)[ok]
    disc = (Xs != sx[None, :]).sum(1)
    if locus_bp is None:
        best = np.where(disc == disc.min())[0]
        runs = np.split(best, np.where(np.diff(best) > 1)[0] + 1)
        run = max(runs, key=len)
        j0 = int(run[len(run) // 2])
    else:
        j0 = int(np.argmin(np.abs(pos - locus_bp)))
    locus_bp = int(pos[j0])
    discordant = [str(n) for n in names[Xs[j0] != sx]]
    skip = set(exclude) | set(discordant)
    rows = []
    for f in range(Xs.shape[1]):
        if str(names[f]) in discordant:
            continue
        s, sex_f = Xs[:, f], 'M' if sx[f] == 1 else 'F'
        i = j0 - 1
        while i >= 0 and s[i] == sx[f]:
            i -= 1
        if i >= 0 and locus_bp - pos[i + 1] <= window_bp:                                    # the boundary between bins i and i+1 is the crossover
            rows.append(dict(fish=str(names[f]), sex=sex_f, boundary_bp=int(pos[i + 1]), bound='SD right of', used=str(names[f]) not in skip))
        i = j0 + 1
        while i < len(s) and s[i] == sx[f]:
            i += 1
        if i < len(s) and pos[i] - locus_bp <= window_bp:
            rows.append(dict(fish=str(names[f]), sex=sex_f, boundary_bp=int(pos[i]), bound='SD left of', used=str(names[f]) not in skip))
    rec = pd.DataFrame(rows, columns=['fish', 'sex', 'boundary_bp', 'bound', 'used']).sort_values('boundary_bp').reset_index(drop=True)
    r = rec[rec.used]
    lo = sorted(r[r.bound == 'SD right of'].boundary_bp, reverse=True)
    up = sorted(r[r.bound == 'SD left of'].boundary_bp)
    lim = pd.DataFrame([dict(k_wrong_per_side=k, left_bp=int(lo[k]), right_bp=int(up[k]), width_bp=int(up[k] - lo[k])) for k in range(max_k + 1) if k < len(lo) and k < len(up)],
                       columns=['k_wrong_per_side', 'left_bp', 'right_bp', 'width_bp'])
    return dict(locus_bp=locus_bp, discordant=discordant, excluded=sorted(skip), recombinants=rec, limits=lim, profile=pd.DataFrame({'pos': pos, 'discordant': disc}))
