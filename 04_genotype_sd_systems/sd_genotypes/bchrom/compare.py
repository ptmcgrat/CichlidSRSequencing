"""Compare the core set found here with Clark et al.'s published core blocks (Appendix A, Core_blocks/*.bed, coordinates on the Mazinzi Reef assembly GCA_000238955.4).
His blocks are lifted to the reference through the cigars of an assembly-to-assembly alignment (the lift routine of the inversion-breakpoint module); a block lifts only if both ends land on
the same chromosome and strand with a plausible span. The report separates blocks on assembled chromosomes (lg*) from blocks on unplaced contigs: most of his core sequence (56% of the bp)
sits on 151 unplaced contigs, which lift least reliably, and the discovery on the reference exists precisely so that the result does not depend on them."""
import numpy as np
import pandas as pd
from cichlid_common import paf as P
from . import calling as C

CLARK_COLS = ['chrom', 'start', 'end', 'length', 'raw_cov_B', 'scaled_cov_B', 'SCR', 'scaled_cov_baseline', 'B_space_len']


def read_clark_bed(path):
    d = pd.read_csv(path, sep=r'\s+', header=None, names=CLARK_COLS)
    d['anchored'] = d.chrom.str.match(r'^lg\d+$')
    return d


def read_name_map(path):
    """Two columns, tab-separated: Clark's sequence name <tab> the name of the same sequence in the FASTA that was aligned."""
    m = pd.read_csv(path, sep='\t', header=None, comment='#')
    return dict(zip(m[0].astype(str), m[1].astype(str)))


def lift_blocks(blocks, paf, name_map=None, max_len_ratio=3.0):
    """One row per block: status ('ok', 'start not aligned', 'end not aligned', 'no consistent lift', 'ambiguous', 'sequence not in the alignment') and, if ok, the lifted interval."""
    have = set(paf.qname.unique())
    rows = []
    for b in blocks.itertuples():
        q = (name_map or {}).get(b.chrom, b.chrom)
        rec = dict(chrom=b.chrom, start=b.start, end=b.end, anchored=bool(getattr(b, 'anchored', True)), query=q)
        if q not in have:
            rows.append({**rec, 'status': 'sequence not in the alignment'}); continue
        hs = [h for h in P.lift_point(paf, q, b.start) if h[3] == 'aligned']
        he = [h for h in P.lift_point(paf, q, b.end - 1) if h[3] == 'aligned']
        if not hs:
            rows.append({**rec, 'status': 'start not aligned'}); continue
        if not he:
            rows.append({**rec, 'status': 'end not aligned'}); continue
        good = []
        for a in hs:
            for z in he:
                if a[0] == z[0] and a[2] == z[2]:
                    span = abs(z[1] - a[1]) + 1; ratio = span / max(1, b.end - b.start)
                    if 1 / max_len_ratio <= ratio <= max_len_ratio:
                        good.append((a[0], min(a[1], z[1]), max(a[1], z[1]) + 1, a[2]))
        good = sorted(set(good))
        if not good:
            rows.append({**rec, 'status': 'no consistent lift'}); continue
        if len(good) > 1:
            rows.append({**rec, 'status': 'ambiguous'}); continue
        t, s, e, st = good[0]
        rows.append({**rec, 'status': 'ok', 'tname': t, 'tstart': s, 'tend': e, 'strand': st})
    return pd.DataFrame(rows)


def compare(lifted, ours, slop=0):
    """Recovery of the lifted blocks by our core set. Returns (per-block table, summary dict)."""
    ok = lifted[lifted.status == 'ok'].copy()
    rec = []
    for r in ok.itertuples():
        g = ours[ours.chrom == r.tname]
        hit = bool(len(g) and ((g.start.values - slop < r.tend) & (g.end.values + slop > r.tstart)).any())
        rec.append(hit)
    ok['recovered'] = rec
    lifted = lifted.merge(ok[['chrom', 'start', 'end', 'recovered']], on=['chrom', 'start', 'end'], how='left')
    lifted['recovered'] = lifted.recovered.fillna(False).astype(bool)
    summ = {}
    for name, g in (('all', lifted), ('anchored', lifted[lifted.anchored]), ('unplaced', lifted[~lifted.anchored])):
        n = len(g); nl = int((g.status == 'ok').sum()); nr = int(g.recovered.sum())
        bp = (g.end - g.start).sum()
        summ[name] = dict(n_clark_blocks=n, bp=int(bp), n_lifted=nl, frac_lifted=nl / n if n else np.nan, n_recovered=nr, frac_recovered_of_lifted=nr / nl if nl else np.nan,
                          frac_recovered_of_all=nr / n if n else np.nan)
    lt = ok.rename(columns={'tname': 'chrom_', 'tstart': 'start_', 'tend': 'end_'})
    if len(lt) and len(ours):
        a = pd.DataFrame(dict(chrom=lt.chrom_, start=lt.start_, end=lt.end_))
        summ['ours_overlapping_lifted_clark'] = dict(n_ours=len(ours), n_overlapping=int(sum(((a.chrom == r.chrom) & (a.start < r.end) & (a.end > r.start)).any() for r in ours.itertuples())))
    return lifted, summ
