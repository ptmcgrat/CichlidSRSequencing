"""The position of a crossover from the READS of one fish (not from the 100 kb calls): at the informative sites that have reads, the fish is in one state (MC: YH reads at the error rate eps) on one side of
the crossover and in the other (heterozygous: YH reads at about q) on the other side. The log likelihood of every split between two read-bearing sites gives the best split and the span that stays within
1 LOD (about 95%) and 2 LOD of it. Combining the nearest left and right crossovers gives the limits of the determiner."""
import os

import numpy as np
import pandas as pd

EPS, Q = 0.02, 0.45


def _ll(n, y, p):
    return y * np.log(p) + (n - y) * np.log(1 - p)


def read_sites(counts_dir, sample, chrom, lo, hi):
    """(pos, reads, YH reads) of the sites with reads in [lo, hi) of one fish's counts file."""
    path = os.path.join(counts_dir, f'{sample}.counts.tsv.gz')
    d = pd.read_csv(path, sep='\t')
    d = d[(d.chrom == chrom) & (d.pos >= lo) & (d.pos < hi) & ((d.n_mc + d.n_yh) > 0)].sort_values('pos')
    return d.pos.to_numpy(), (d.n_mc + d.n_yh).to_numpy(), d.n_yh.to_numpy()


def localize(pos, n, y, eps=EPS, q=Q):
    """The best split of the sites into two states and its support. state = '0->1' (MC then heterozygous) or '1->0'. gap = distance between the two read-bearing sites around the best split;
    edge = the best split is before the first or after the last site: the crossover may lie outside the window."""
    if len(pos) < 2:
        raise ValueError('fewer than two sites with reads: nothing to localize')
    l0, l1 = _ll(n, y, eps), _ll(n, y, q)
    best = None
    for a, la, lb in ((0, l0, l1), (1, l1, l0)):
        tot = np.concatenate([[0], np.cumsum(la)]) + np.concatenate([np.cumsum(lb[::-1])[::-1], [0]])
        if best is None or tot.max() > best[0].max():
            best = (tot, a)
    tot, a = best
    lod = (tot - tot.max()) / np.log(10)
    k = int(tot.argmax())

    def span(th):
        ks = np.where(lod >= -th)[0]
        return int(pos[max(ks.min() - 1, 0)]), int(pos[min(ks.max(), len(pos) - 1)])
    lo, hi = int(pos[max(k - 1, 0)]), int(pos[min(k, len(pos) - 1)])
    return dict(state=f'{a}->{1 - a}', best=(lo + hi) / 2, gap=hi - lo, span1=span(1), span2=span(2), sites=len(pos), reads=int(n.sum()), edge=bool(k == 0 or k == len(pos)), split_index=k)


def fitted_rates(pos, n, y, split_index, state):
    """eps and q fitted from the flanks of the split (the YH read share on the MC side and on the heterozygous side), for the sensitivity check."""
    k = split_index
    left, right = (n[:k].sum(), y[:k].sum()), (n[k:].sum(), y[k:].sum())
    fl, fr = left[1] / max(left[0], 1), right[1] / max(right[0], 1)
    return (max(fl, 0.002), fr) if state == '0->1' else (max(fr, 0.002), fl)


def sweep(pos, n, y, eps_list=(0.01, 0.03), q_list=(0.35, 0.45, 0.55)):
    """The best split under a grid of read rates: a robust crossover does not move."""
    return {(e, q): localize(pos, n, y, e, q)['best'] for e in eps_list for q in q_list}


def interval(res, roles):
    """(best-estimate interval, conservative interval) of the determiner from the crossovers `res` {fish: localize dict} with roles {fish: 'left' | 'right'}: it lies right of the left fish's crossover and left
    of the right fish's. Conservative = the ends of the 1-LOD spans (the crossover could be anywhere in its span)."""
    left = [res[f] for f in res if roles[f] == 'left']
    right = [res[f] for f in res if roles[f] == 'right']
    if not left or not right:
        raise ValueError('the crossovers need at least one left and one right fish')
    best = (max(r['best'] for r in left), min(r['best'] for r in right))
    cons = (max(r['span1'][0] for r in left), min(r['span1'][1] for r in right))
    return best, cons


def roles_from_finemap(recombinants, per_side=2):
    """{fish: role} of the `per_side` nearest USED crossovers on each side of the limits (a fish bounding on the left is a 'left' fish)."""
    r = recombinants[recombinants.used]
    left = r[r.bound == 'SD right of'].sort_values('boundary_bp', ascending=False).drop_duplicates('fish').head(per_side)
    right = r[r.bound == 'SD left of'].sort_values('boundary_bp').drop_duplicates('fish').head(per_side)
    roles = {f: 'left' for f in left.fish}
    roles.update({f: 'right' for f in right.fish})
    return roles
