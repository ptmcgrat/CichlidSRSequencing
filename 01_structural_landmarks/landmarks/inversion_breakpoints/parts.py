"""The parts of a tandem (double) inversion: where each part lies, and which assembly carries which (module 01; LG11 = 11a + 11b, LG20 = 20a + 20b).

The caller reports a tandem inversion A B' C' D as ONE inversion with two segments and an internal breakpoint (paf.inversion_segments): right for an assembly that carries both parts.
An assembly that carries only the first part (Rhamphochromis on LG20, Kumar et al. 2025: it lacks the 4.1 Mb second part) shows an ordinary single inversion whose RIGHT breakpoint is the
internal breakpoint of the double. The consensus asks that both breakpoints agree, so it lists that assembly as a different inversion, and nothing links the two. This module does:

  * the structure comes from the assemblies that resolve two segments (the best two-segment consensus row): the outer breakpoints, and the internal breakpoint as the gap between the segments;
  * a single-part carrier is recognised by one outer breakpoint plus the internal one, and its (narrow) breakpoint interval SHARPENS the internal breakpoint, which the segment gap alone leaves
    hundreds of kb wide;
  * every assembly gets a state: both parts, the first only, the second only, both with the parts not resolved, another inversion, or none; optionally compared with what the literature says (`expect`).

Parts are named left to right in REFERENCE coordinates (parts: [20a, 20b]); only two parts are supported. Nothing here replaces the consensus: it adds `<LG>_parts.tsv`, `<LG>_assembly_states.tsv`
and `<LG>_parts_breakpoint_regions.bed` next to it."""
import numpy as np
import pandas as pd

from cichlid_common.paf import METHOD_RANK, _gap

PARTS_COLUMNS = ['part', 'left_breakpoint', 'left_min', 'left_max', 'right_breakpoint', 'right_min', 'right_max', 'length', 'carriers', 'n_carriers', 'evidence']
STATE_COLUMNS = ['assembly', 'state', 'parts_carried', 'method', 'detail', 'expected', 'matches_expectation']


def _num(x, default=np.nan):
    try:
        v = float(x)
        return v if v == v else default
    except (TypeError, ValueError):
        return default


def _carried(state, parts):
    """The set of parts a state says an assembly carries."""
    p1, p2 = parts
    return {'both': {p1, p2}, 'both (parts not resolved)': {p1, p2}, f'{p1} only': {p1}, f'{p2} only': {p2}, 'none': set()}.get(state)


def resolve_parts(calls, cons, parts, expect=None, tol=150_000, gap=100_000, assemblies=None):
    """(parts_table, assembly_states, notes). calls: the per-assembly calls (assembly, method, left_lo/hi, right_lo/hi, n_segments, segments); cons: their consensus; parts: the two part names,
    left to right in reference coordinates; expect: {assembly: [parts it should carry]} from the literature; tol/gap: the consensus tolerances (config calling)."""
    if len(parts) != 2:
        raise ValueError(f'a double inversion has two parts, not {len(parts)}: {parts}')
    p1, p2 = parts
    notes = []
    calls = calls.copy() if len(calls) else pd.DataFrame(columns=['assembly', 'method', 'left_lo', 'left_hi', 'right_lo', 'right_hi'])
    for col in ('n_segments', 'left_lo', 'left_hi', 'right_lo', 'right_hi'):
        calls[col] = pd.to_numeric(calls[col], errors='coerce') if col in calls.columns else np.nan
    names = list(dict.fromkeys(list(assemblies or []) + (calls.assembly.tolist() if len(calls) else [])))
    # ---------------------------------------------------------------- the structure, from the assemblies that resolve two segments
    two = cons[pd.to_numeric(cons.n_segments, errors='coerce') == 2] if len(cons) and 'n_segments' in cons.columns else pd.DataFrame()
    structure = None
    if len(two):
        c = two.iloc[0]
        segs = [tuple(int(v) for v in s.split('-')) for s in str(c.segments).split(';') if s]
        if len(segs) == 2:
            lo, hi = sorted((segs[0][1], segs[1][0]))
            structure = dict(L=(int(c.left_min), int(c.left_max)), Lbp=int(c.left_breakpoint), R=(int(c.right_min), int(c.right_max)), Rbp=int(c.right_breakpoint), int_lo=lo, int_hi=hi,
                             sources=str(c.sources), segs=segs)
    if structure is None:
        notes.append(f'no assembly resolves two segments: the two parts {p1} and {p2} cannot be placed; the inversions that were called are in the consensus table')
        rows = []
        for a in names:
            sub = calls[calls.assembly == a]
            rows.append(dict(assembly=a, state='another inversion' if len(sub) else 'none', parts_carried='' if not len(sub) else '?', method=sub.method.iloc[0] if len(sub) else '',
                             detail='a single inversion was called, but the parts are not resolved' if len(sub) else 'no inversion called on this chromosome', expected='', matches_expectation=''))
        states = pd.DataFrame(rows, columns=STATE_COLUMNS)
        return pd.DataFrame(columns=PARTS_COLUMNS), _compare(states, expect, parts), notes
    L, R, lo, hi = structure['L'], structure['R'], structure['int_lo'], structure['int_hi']
    # ---------------------------------------------------------------- classify every call against the structure
    def match(row):
        near_L = _gap(row.left_lo, row.left_hi, L[0], L[1]) <= gap
        near_R = _gap(row.right_lo, row.right_hi, R[0], R[1]) <= gap
        int_as_right = _gap(row.right_lo, row.right_hi, lo - tol, hi + tol) == 0
        int_as_left = _gap(row.left_lo, row.left_hi, lo - tol, hi + tol) == 0
        n = row.n_segments
        if near_L and near_R and n == 2:
            return 'both'
        if near_L and int_as_right and not near_R:
            return f'{p1} only'
        if int_as_left and near_R and not near_L:
            return f'{p2} only'
        if near_L and near_R:
            return 'both (parts not resolved)'
        return 'other'

    calls['_match'] = [match(r) for r in calls.itertuples()] if len(calls) else []
    calls['_rank'] = calls.method.map(lambda m: METHOD_RANK.get(m, 9)) if len(calls) else []
    # the internal breakpoint: the gap between the two segments, sharpened by the single-part carriers (their breakpoint interval is the breakpoint itself)
    sharp = []
    for r in calls[calls._match == f'{p1} only'].itertuples():
        sharp.append((r.right_lo, r.right_hi, r.assembly))
    for r in calls[calls._match == f'{p2} only'].itertuples():
        sharp.append((r.left_lo, r.left_hi, r.assembly))
    int_lo, int_hi, how = lo, hi, f'the gap between the two segments of {structure["sources"]}'
    if sharp:
        ilo, ihi = max(s[0] for s in sharp), min(s[1] for s in sharp)
        who = ', '.join(sorted({s[2] for s in sharp}))
        if ilo <= ihi:
            if ilo >= lo - tol and ihi <= hi + tol:
                int_lo, int_hi, how = int(ilo), int(ihi), f'the breakpoint of the single-part carrier(s) {who}, inside the segment gap'
            else:
                int_lo, int_hi, how = int(ilo), int(ihi), f'the breakpoint of the single-part carrier(s) {who}'
                notes.append(f'the internal breakpoint of the single-part carrier(s) {who} ({int(ilo)}-{int(ihi)}) lies outside the segment gap of the double carriers ({lo}-{hi}): check the alignments')
        else:
            mid = int(np.median([(s[0] + s[1]) / 2 for s in sharp]))
            int_lo, int_hi, how = mid, mid, f'the median of single-part carriers that disagree ({who})'
            notes.append(f'single-part carriers disagree on the internal breakpoint: {sorted((int(s[0]), int(s[1]), s[2]) for s in sharp)}')
    ibp = int((int_lo + int_hi) / 2)
    # ---------------------------------------------------------------- who carries what
    state_rows, carriers = [], {p1: [], p2: []}
    for a in names:
        sub = calls[calls.assembly == a]
        pick = None
        for st in ('both', f'{p1} only', f'{p2} only', 'both (parts not resolved)'):
            m = sub[sub._match == st]
            if len(m):
                pick = (st, m.sort_values('_rank').iloc[0]); break
        if pick:
            st, r = pick
            carried = _carried(st, parts)
            detail = f'{r.method}: left {int(r.left_lo)}-{int(r.left_hi)}, right {int(r.right_lo)}-{int(r.right_hi)}' + (f', {int(r.n_segments)} segment(s)' if r.n_segments == r.n_segments else '')
            state_rows.append(dict(assembly=a, state=st, parts_carried='+'.join(sorted(carried, key=parts.index)), method=r.method, detail=detail, expected='', matches_expectation=''))
            for p in carried:
                carriers[p].append(a)
        elif len(sub):
            r = sub.sort_values('_rank').iloc[0]
            state_rows.append(dict(assembly=a, state='another inversion', parts_carried='?', method=r.method,
                                   detail=f'{r.method}: left {int(r.left_lo)}-{int(r.left_hi)}, right {int(r.right_lo)}-{int(r.right_hi)}: not the breakpoints of {p1} or {p2}', expected='', matches_expectation=''))
        else:
            state_rows.append(dict(assembly=a, state='none', parts_carried='', method='', detail='no inversion called on this chromosome (collinear with the reference, or no alignment)', expected='', matches_expectation=''))
    states = _compare(pd.DataFrame(state_rows, columns=STATE_COLUMNS), expect, parts)
    # ---------------------------------------------------------------- the parts table
    ev = lambda side, extra='': f'{side}: {structure["sources"]}{extra}'
    rows = [dict(part=p1, left_breakpoint=structure['Lbp'], left_min=L[0], left_max=L[1], right_breakpoint=ibp, right_min=int_lo, right_max=int_hi, evidence=f'outer {ev("left")}; internal: {how}'),
            dict(part=p2, left_breakpoint=ibp, left_min=int_lo, left_max=int_hi, right_breakpoint=structure['Rbp'], right_min=R[0], right_max=R[1], evidence=f'internal: {how}; outer {ev("right")}')]
    for r in rows:
        r['length'] = r['right_breakpoint'] - r['left_breakpoint']
        r['carriers'] = ';'.join(sorted(carriers[r['part']]))
        r['n_carriers'] = len(carriers[r['part']])
    pt = pd.DataFrame(rows, columns=PARTS_COLUMNS)
    if not any(s == f'{p1} only' or s == f'{p2} only' for s in states.state):
        notes.append(f'no assembly carries a single part: the internal breakpoint rests on the gap between the two segments ({lo}-{hi}, {hi - lo} bp wide)')
    return pt, states, notes


def _compare(states, expect, parts):
    states = states.copy()
    states['matches_expectation'] = states['matches_expectation'].astype(object)           # True, False or '' (no expectation): a strict string column would refuse the booleans
    for i, r in states.iterrows():
        want = (expect or {}).get(r.assembly)
        if want is None:
            continue
        carried = _carried(r.state, parts)
        states.at[i, 'expected'] = '+'.join(want) if want else 'none'
        states.at[i, 'matches_expectation'] = '' if carried is None else bool(carried == set(want))
    return states


def parts_bed(parts_table, tname, lg):
    """BED lines of the breakpoint intervals of every part: <chrom> <min> <max> <LG>_<part>_<left|right>."""
    out = []
    for r in parts_table.itertuples():
        out.append(f'{tname}\t{int(r.left_min)}\t{int(r.left_max)}\t{lg}_{r.part}_left')
        out.append(f'{tname}\t{int(r.right_min)}\t{int(r.right_max)}\t{lg}_{r.part}_right')
    return out


def format_parts(parts_table, states, notes):
    out = []
    if len(parts_table):
        out.append(parts_table[['part', 'left_breakpoint', 'right_breakpoint', 'length', 'n_carriers', 'carriers']].to_string(index=False))
    out.append(states[['assembly', 'state', 'expected', 'matches_expectation']].to_string(index=False))
    out += [f'note: {n}' for n in notes]
    bad = states[states.matches_expectation == False]       # noqa: E712 (the column holds True, False or '')
    if len(bad):
        out.append('DIFFERS FROM THE LITERATURE: ' + ', '.join(f'{r.assembly} is {r.state!r}, expected {r.expected}' for r in bad.itertuples()))
    return '\n'.join(out)
