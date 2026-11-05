"""Inversion breakpoints from assembly-to-reference alignments (minimap2 PAF with cigar). Pure functions, unit-tested on synthetic alignments.
Three complementary methods, all in REFERENCE coordinates:
  A  strand switch inside one contig: flank / inverted segment / flank, whichever strand the contig has (+ - + or - + -) -> intervals = gaps between the flanks and the inverted run
     (A-half: the contig ends at one breakpoint; assumes the reference orientation, flagged strand_dependent)
  B  orientation discordance between haplotigs: reference windows covered by contigs on BOTH strands (heterozygous assembly) -> inversion span
  C  extent of an all-'-' contig that has no flanks (low confidence: the contig may simply end at the breakpoint repeats; strand_dependent)
Polarity: only the strand relative to the reference is observed; whether '-' means inverted depends on the reference's own orientation (config reference_is_standard).
Minus-strand logic: along the query contig the inverted run has DEcreasing reference coordinates, so the LEFT breakpoint (low reference position) is where
the LAST '-' block ends (min tstart) and the RIGHT breakpoint where the FIRST '-' block starts (max tend)."""
import os, re
import numpy as np
import pandas as pd

PAF_COLS = ['qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', 'tlen', 'tstart', 'tend', 'nmatch', 'alen', 'mapq']
INTCOLS = ['qlen', 'qstart', 'qend', 'tlen', 'tstart', 'tend', 'nmatch', 'alen', 'mapq']


def parse_paf(src, tnames=None):
    """src: path or iterable of lines. tnames: keep only these target sequences. cigar kept in column 'cigar' (tag cg)."""
    lines = open(src) if isinstance(src, str) and os.path.exists(src) else src
    rows = []
    for line in lines:
        f = line.rstrip('\n').split('\t')
        if len(f) < 12:
            continue
        if tnames is not None and f[5] not in tnames:
            continue
        rec = dict(zip(PAF_COLS, f[:12]))
        for k in INTCOLS:
            rec[k] = int(rec[k])
        tags = {t[:2]: t[5:] for t in f[12:] if len(t) > 5 and t[2] == ':'}
        rec['cigar'] = tags.get('cg')
        rows.append(rec)
    return pd.DataFrame(rows, columns=PAF_COLS + ['cigar'])


def drop_blocks_in_regions(df, regions, max_overlap=0.5):
    """Remove alignment blocks whose reference span lies mostly (> max_overlap) inside a region the reference is known to have wrong (e.g. a misplaced repeat chunk)."""
    if not regions or not len(df):
        return df
    keep = np.ones(len(df), bool)
    ts, te = df.tstart.values, df.tend.values
    for s, e in regions:
        ov = (np.minimum(te, e) - np.maximum(ts, s)).clip(min=0)
        keep &= ~(ov / np.maximum(te - ts, 1) > max_overlap)
    return df[keep].reset_index(drop=True)


# ------------------------------------------------------------------ method A: strand switches inside one contig
def main_target(df, qname):
    g = df[df.qname == qname]
    return g.groupby('tname').alen.sum().idxmax() if len(g) else None


def blocks_for_contig(df, qname, tname, min_block):
    g = df[(df.qname == qname) & (df.tname == tname) & ((df.qend - df.qstart) >= min_block)].sort_values('qstart')
    return g.reset_index(drop=True)


def strand_runs(g, despeckle=0):
    """Consecutive blocks (in query order) with the same strand -> list of dicts. despeckle > 0 absorbs a run shorter than that (reference span) that sits between two runs
    of the SAME other strand (the tiny opposite-strand blocks that repeats scatter along a diagonal would otherwise split one long run into many)."""
    runs = []
    for r in g.itertuples():
        if runs and runs[-1]['strand'] == r.strand:
            x = runs[-1]; x['tmin'] = min(x['tmin'], r.tstart); x['tmax'] = max(x['tmax'], r.tend); x['qend'] = r.qend; x['n'] += 1
        else:
            runs.append(dict(strand=r.strand, tmin=r.tstart, tmax=r.tend, qstart=r.qstart, qend=r.qend, n=1))
    if despeckle:
        changed = True
        while changed:
            changed = False
            for i in range(1, len(runs) - 1):
                a, b, c = runs[i - 1], runs[i], runs[i + 1]
                if b['tmax'] - b['tmin'] < despeckle and a['strand'] == c['strand'] != b['strand']:
                    runs[i - 1:i + 2] = [dict(strand=a['strand'], tmin=min(a['tmin'], b['tmin'], c['tmin']), tmax=max(a['tmax'], b['tmax'], c['tmax']),
                                              qstart=a['qstart'], qend=c['qend'], n=a['n'] + b['n'] + c['n'])]
                    changed = True
                    break
    return runs


def inversions_in_contig(g, min_inversion=500_000, max_minus_gap=None, adjacent_tol=500_000, despeckle=300_000, max_interval=1_000_000):
    """Inversions inside one contig, independent of which strand the contig has relative to the reference.
    A (full): the MIDDLE run of a triplet whose two neighbours share the opposite strand and sit on either side of it in the reference
              ( + - + for a contig in the reference orientation, - + - for a reverse-complemented contig ): breakpoint intervals from the flank adjacencies.
    A-half / C: a contig that ends at a breakpoint (two runs) or is only the inverted segment (one run). With no flank there is nothing to compare against, so these
              ASSUME the reference orientation: the '-' run is the inversion. They are flagged strand_dependent and are low confidence."""
    runs = strand_runs(g, despeckle)
    out = []

    def row(method, left, right, r, support, dependent):
        return dict(method=method, left_lo=left[0], left_hi=left[1], right_lo=right[0], right_hi=right[1], left_support=support[0], right_support=support[1],
                    inv_ref_len=r['tmax'] - r['tmin'], inv_query_len=r['qend'] - r['qstart'], inv_q_start=r['qstart'], inv_q_end=r['qend'], n_blocks=r['n'],
                    inv_strand=r['strand'], strand_dependent=dependent)

    for i in range(1, len(runs) - 1):
        r, p, n = runs[i], runs[i - 1], runs[i + 1]
        if r['tmax'] - r['tmin'] < min_inversion or p['strand'] == r['strand'] or n['strand'] == r['strand']:
            continue
        # the flank strand must agree with where the flanks lie: '+' flanks run left -> right in the reference, '-' flanks right -> left. A '+' run between two '-' runs that are
        # in the wrong places is a translocation or an assembly join, not a simple inversion.
        if p['strand'] == '+' and p['tmax'] <= r['tmin'] + adjacent_tol and n['tmin'] >= r['tmax'] - adjacent_tol:        # flanks '+', inversion '-'
            low, high = p, n
        elif p['strand'] == '-' and n['tmax'] <= r['tmin'] + adjacent_tol and p['tmin'] >= r['tmax'] - adjacent_tol:      # contig reverse-complemented: flanks '-', inversion '+'
            low, high = n, p
        else:
            continue
        left, right = tuple(sorted((low['tmax'], r['tmin']))), tuple(sorted((r['tmax'], high['tmin'])))
        if left[1] - left[0] > max_interval or right[1] - right[0] > max_interval:                   # a 'breakpoint' that is megabases wide is not a breakpoint
            continue
        out.append(row('A', left, right, r, ('flank', 'flank'), False))
    if out:
        return out
    if len(runs) == 2:
        for r, f in ((runs[0], runs[1]), (runs[1], runs[0])):
            if r['strand'] != '-' or r['tmax'] - r['tmin'] < min_inversion:
                continue
            if f['tmax'] <= r['tmin'] + adjacent_tol:
                out.append(row('A-half', tuple(sorted((f['tmax'], r['tmin']))), (r['tmax'], r['tmax']), r, ('flank', 'extent'), True))
            elif f['tmin'] >= r['tmax'] - adjacent_tol:
                out.append(row('A-half', (r['tmin'], r['tmin']), tuple(sorted((r['tmax'], f['tmin']))), r, ('extent', 'flank'), True))
    elif len(runs) == 1 and runs[0]['strand'] == '-' and runs[0]['tmax'] - runs[0]['tmin'] >= min_inversion:
        r = runs[0]; out.append(row('C', (r['tmin'], r['tmin']), (r['tmax'], r['tmax']), r, ('extent', 'extent'), True))
    return out


def call_method_AC(df, tname, min_block=20_000, min_inversion=500_000, max_minus_gap=None, despeckle=300_000, max_interval=1_000_000):
    rows = []
    for q in df[df.tname == tname].qname.unique():
        if main_target(df, q) != tname:
            continue
        g = blocks_for_contig(df, q, tname, min_block)
        for inv in inversions_in_contig(g, min_inversion, None, 500_000, despeckle, max_interval):
            rows.append(dict(qname=q, **inv))
    return pd.DataFrame(rows)


def runs_table(df, tname, contigs, min_block=20_000, despeckle=300_000):
    """The strand runs of the given contigs after despeckling (what the caller actually sees): one row per run, in query order."""
    rows = []
    for q in contigs:
        g = blocks_for_contig(df, q, tname, min_block)
        for i, r in enumerate(strand_runs(g, despeckle)):
            rows.append(dict(contig=q, run=i + 1, strand=r['strand'], ref_start=r['tmin'], ref_end=r['tmax'], ref_span=r['tmax'] - r['tmin'], query_start=r['qstart'], query_end=r['qend'], n_blocks=r['n']))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ method B: orientation discordance between contigs
def window_strand_coverage(df, tname, tlen, window, min_block):
    """Matrix of aligned bp per (window, strand) for each contig: returns dict qname -> (plus[n_win], minus[n_win])."""
    n = tlen // window + 1
    cov = {}
    for r in df[(df.tname == tname) & ((df.tend - df.tstart) >= min_block)].itertuples():
        arr = cov.setdefault(r.qname, (np.zeros(n), np.zeros(n)))[0 if r.strand == '+' else 1]
        w0, w1 = r.tstart // window, (r.tend - 1) // window
        for w in range(w0, w1 + 1):
            arr[w] += max(0, min(r.tend, (w + 1) * window) - max(r.tstart, w * window))
    return cov


def discordance_regions(df, tname, tlen, window=10_000, min_window_bp=2_000, min_block=20_000, min_len=500_000, max_gap_windows=20):
    """Windows covered on '+' by at least one contig AND on '-' by at least one contig; longest run (gaps <= max_gap_windows tolerated)."""
    cov = window_strand_coverage(df, tname, tlen, window, min_block)
    if not cov:
        return pd.DataFrame()
    plus = np.max([v[0] for v in cov.values()], axis=0); minus = np.max([v[1] for v in cov.values()], axis=0)
    disc = (plus >= min_window_bp) & (minus >= min_window_bp)
    idx = np.where(disc)[0]
    regs = []
    if len(idx):
        s = p = idx[0]
        for w in list(idx[1:]) + [None]:
            if w is None or w - p > max_gap_windows + 1:
                if (p - s + 1) * window >= min_len:
                    regs.append((int(s * window), int((p + 1) * window)))
                if w is not None:
                    s = w
            if w is not None:
                p = w
    rows = []
    for lo, hi in regs:
        rows.append(dict(method='B', left_lo=max(0, lo - window), left_hi=lo + window, right_lo=hi - window, right_hi=hi + window,
                         left_support='discordance', right_support='discordance', inv_ref_len=hi - lo, inv_query_len=np.nan, qname='haplotig-pair', n_blocks=np.nan))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ consensus across assemblies / methods
METHOD_RANK = {'A': 0, 'B': 1, 'A-half': 2, 'C': 3}     # a call with both flanks beats a haplotig comparison, which beats calls that assume the reference orientation


def consensus(calls, tol=150_000):
    """calls: DataFrame with columns assembly, method, left_lo/hi, right_lo/hi. Cluster left and right midpoints within `tol`; returns one row per cluster pair,
    ranked by the best method in the cluster (A > B > A-half > C), then by the number of supporting calls, then by length. confidence: high (has an A call), medium (B only), low."""
    if not len(calls):
        return pd.DataFrame()
    c = calls.copy(); c['left_mid'] = (c.left_lo + c.left_hi) / 2; c['right_mid'] = (c.right_lo + c.right_hi) / 2
    c = c.sort_values('left_mid').reset_index(drop=True)
    cl, cur = [], [0]
    for i in range(1, len(c)):
        if abs(c.left_mid[i] - c.left_mid[cur[-1]]) <= tol and abs(c.right_mid[i] - c.right_mid[cur[-1]]) <= tol:
            cur.append(i)
        else:
            cl.append(cur); cur = [i]
    cl.append(cur)
    rows = []
    for ix in cl:
        s = c.loc[ix]
        methods = set(s.method)
        best = min(METHOD_RANK.get(m, 9) for m in methods)
        conf = 'high' if 'A' in methods else ('medium' if 'B' in methods else 'low')       # low = only calls that assume the reference orientation (A-half, C)
        rows.append(dict(left_breakpoint=int(s.left_mid.median()), left_min=int(s.left_lo.min()), left_max=int(s.left_hi.max()),
                         right_breakpoint=int(s.right_mid.median()), right_min=int(s.right_lo.min()), right_max=int(s.right_hi.max()),
                         inversion_length=int(s.right_mid.median() - s.left_mid.median()), n_calls=len(ix), confidence=conf, _rank=best,
                         sources=';'.join(sorted(set(s.assembly + ':' + s.method)))))
    out = pd.DataFrame(rows).sort_values(['_rank', 'n_calls', 'inversion_length'], ascending=[True, False, False]).drop(columns='_rank').reset_index(drop=True)
    return out


# ------------------------------------------------------------------ liftover through the PAF cigar
_CIG = re.compile(r'(\d+)([MIDNSHP=X])')


def lift_point(df, qname, qpos):
    """Map a query position (e.g. a YH coordinate) to the target through the alignment cigars. Returns list of (tname, tpos, strand, status)."""
    out = []
    for r in df[(df.qname == qname) & (df.qstart <= qpos) & (df.qend > qpos)].itertuples():
        if r.cigar is None:
            continue
        t = r.tstart; q = r.qstart if r.strand == '+' else r.qend
        res = None
        for n, op in ((int(a), b) for a, b in _CIG.findall(r.cigar)):
            if op in 'M=X':
                if r.strand == '+' and q <= qpos < q + n:
                    res = (r.tname, t + (qpos - q), '+', 'aligned'); break
                if r.strand == '-' and q - n <= qpos < q:
                    res = (r.tname, t + (q - 1 - qpos), '-', 'aligned'); break
                t += n; q += n if r.strand == '+' else -n
            elif op == 'I':
                if r.strand == '+' and q <= qpos < q + n:
                    res = (r.tname, t, '+', 'in insertion'); break
                if r.strand == '-' and q - n <= qpos < q:
                    res = (r.tname, t, '-', 'in insertion'); break
                q += n if r.strand == '+' else -n
            elif op in 'DN':
                t += n
        if res:
            out.append(res)
    return out


# ------------------------------------------------------------------ polarity and interval liftover
def strand_polarity(reference_is_standard):
    """What the inversion means. The strand of the inverted segment is NOT meaningful (a contig can be reverse-complemented as a whole); what is meaningful is that the assembly
    differs from the reference between the breakpoints. reference_is_standard: True / False / None (unknown)."""
    if reference_is_standard is True:
        return "reference = standard orientation: the assembly's orientation between the breakpoints is the INVERTED (derived) state"
    if reference_is_standard is False:
        return "reference carries the INVERTED orientation: the assembly's orientation between the breakpoints is the standard state"
    return "orientation of the reference at this inversion is not established: the assembly differs from the reference between the breakpoints"


def lift_interval(df, qname, start, end, step=25_000, max_scan=3_000_000):
    """Lift a query interval (e.g. a centromere found in the YH assembly) to the target. Centromeric satellite arrays usually do not align, so the interval is
    bounded by the nearest ALIGNABLE positions outside it: scan left from `start` and right from `end` until a position lifts. Returns dict or None."""
    def scan(pos, direction):
        for k in range(0, max_scan + 1, step):
            p = pos + direction * k
            if p < 0:
                return None
            hits = [h for h in lift_point(df, qname, p) if h[3] == 'aligned']
            if hits:
                return hits[0], p
        return None
    left, right = scan(start, -1), scan(end, +1)
    if left is None or right is None:
        return None
    (tl, pl, sl, _), (tr, pr, sr, _) = left[0], right[0]
    if tl != tr:
        return dict(status='flanks map to different chromosomes', tname_left=tl, tname_right=tr)
    lo, hi = sorted((pl, pr))
    return dict(status='ok' if sl == sr else 'flanks on opposite strands (interval straddles an inversion breakpoint)', tname=tl, tstart=lo, tend=hi, strand=sl,
                left_flank_distance=start - left[1], right_flank_distance=right[1] - end)
