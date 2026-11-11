"""The inversion breakpoints of this project: the published intervals of Kumar et al. 2025 (Supplementary File 4) for LG2, 9, 10, 11 (11a, 11b), 13 and 20 (20a, 20b), and this pipeline's own call for LG5.

SPEC decision 70: the breakpoints are taken from the literature and not refined further here. The table is `data/published_inversions_<genome>.tsv`; this module reads it, validates it and derives
what the other modules need. Nothing else in the repository should carry a copy of these numbers.

    from landmarks.inversions_table import load_inversions, breakpoint_windows
    inv = load_inversions()                    # refuses a table with a `pending` row (LG5 until its numbers are in); allow_pending=True returns it with empty coordinates
    win = breakpoint_windows(inv, pad=10_000)  # chrom, start, end, name: one interval per breakpoint, a shared tandem breakpoint once"""
import os
import re

import numpy as np
import pandas as pd

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
DEFAULT_GENOME = 'Mzebra_GT3_NCBI'
COLUMNS = ['lg', 'part', 'chrom', 'bp1_start', 'bp1_end', 'bp2_start', 'bp2_end', 'reference_state', 'status', 'source', 'note']
COORDS = ['bp1_start', 'bp1_end', 'bp2_start', 'bp2_end']
STATUSES = ('published', 'pipeline', 'pending')
ACCESSION_BASE = 135166                       # Mzebra_GT3_NCBI: LG n is NC_(135166 + n).1


def table_path(genome=DEFAULT_GENOME):
    return os.path.join(DATA, f'published_inversions_{genome}.tsv')


def validate(d, genome=DEFAULT_GENOME):
    """The problems of a raw table (every column a string), one sentence each; empty when it is sound."""
    bad = []
    missing = [c for c in COLUMNS if c not in d.columns]
    if missing:
        return [f'missing columns: {", ".join(missing)}']
    if d.part.duplicated().any():
        bad.append('duplicate part names: ' + ', '.join(sorted(d.part[d.part.duplicated()])))
    for r in d.itertuples():
        tag = f'{r.part} ({r.lg})'
        if r.status not in STATUSES:
            bad.append(f'{tag}: status {r.status!r} is not one of {", ".join(STATUSES)}')
        if not re.fullmatch(r'LG\d+', r.lg):
            bad.append(f'{tag}: lg must look like LG11')
            continue
        if not re.fullmatch(r'NC_\d+\.\d+', r.chrom):
            bad.append(f'{tag}: chrom {r.chrom!r} is not an NC_ accession')
        elif genome == DEFAULT_GENOME and int(r.chrom[3:].split('.')[0]) - ACCESSION_BASE != int(r.lg[2:]):
            bad.append(f'{tag}: {r.chrom} is not the accession of {r.lg} in {genome}')
        if not r.source.strip():
            bad.append(f'{tag}: no source')
        vals = [getattr(r, c) for c in COORDS]
        if r.status == 'pending':
            if any(v != '' for v in vals):
                bad.append(f'{tag}: a pending row must have empty coordinates')
            continue
        if any(not re.fullmatch(r'\d+', v) for v in vals):
            bad.append(f'{tag}: the four coordinates must be whole numbers')
            continue
        a, b, c, e = (int(v) for v in vals)
        if a > b or c > e:
            bad.append(f'{tag}: an interval starts after it ends')
        if b > c:
            bad.append(f'{tag}: the left breakpoint interval ends after the right one starts')
    # a tandem pair (11a, 11b): the right breakpoint of the first part is the left breakpoint of the second
    for _, g in d.groupby('lg'):
        parts = {p: row for p, row in zip(g.part, g.itertuples())}
        for p, row in parts.items():
            m = re.fullmatch(r'(\d+)([a-z])', p)
            nxt = parts.get(f'{m.group(1)}{chr(ord(m.group(2)) + 1)}') if m else None
            if nxt is not None and (row.bp2_start, row.bp2_end, row.chrom) != (nxt.bp1_start, nxt.bp1_end, nxt.chrom):
                bad.append(f'{p} and {nxt.part} are a tandem pair but share no breakpoint: the right interval of {p} must equal the left interval of {nxt.part}')
    return bad


def load_inversions(genome=DEFAULT_GENOME, path=None, allow_pending=False):
    """The table as a DataFrame with integer coordinates and: bp1_mid, bp2_mid (the middle of each interval), length_min, length_max (the inversion's length from the interval ends:
    bp2_start - bp1_end, bp2_end - bp1_start). Raises ValueError for an unsound table, and for a `pending` row unless allow_pending."""
    d = pd.read_csv(path or table_path(genome), sep='\t', comment='#', dtype=str, keep_default_na=False)
    problems = validate(d, genome)
    if problems:
        raise ValueError('the inversion table is not sound: ' + '; '.join(problems))
    pending = d[d.status == 'pending']
    if len(pending) and not allow_pending:
        raise ValueError(f"{len(pending)} inversion(s) have no breakpoints yet ({', '.join(pending.part)}): copy their numbers into the table, or pass allow_pending=True")
    for c in COORDS:
        d[c] = pd.to_numeric(d[c].replace('', np.nan)).astype('Int64')
    d['bp1_mid'] = (d.bp1_start + d.bp1_end) // 2
    d['bp2_mid'] = (d.bp2_start + d.bp2_end) // 2
    d['length_min'] = d.bp2_start - d.bp1_end
    d['length_max'] = d.bp2_end - d.bp1_start
    return d


def breakpoint_windows(inversions, pad=0):
    """The intervals the read-level analyses look at: chrom, start, end, name for every breakpoint, each interval widened by `pad` (never below 0). The shared breakpoint of a tandem pair
    is listed once, as inv<LG number>_internal; the others are inv<part>_bp1 and inv<part>_bp2. Rows without coordinates (pending) are left out."""
    rows, seen = [], set()
    d = inversions[inversions.status != 'pending']
    tandem = {(r.lg, r.part) for r in d.itertuples() if re.fullmatch(r'\d+[a-z]', r.part)}
    for r in d.itertuples():
        for side, (s, e) in (('bp1', (r.bp1_start, r.bp1_end)), ('bp2', (r.bp2_start, r.bp2_end))):
            letter = re.fullmatch(r'\d+([a-z])', r.part)
            internal = bool(letter) and ((side == 'bp2' and (r.lg, r.part[:-1] + chr(ord(letter.group(1)) + 1)) in tandem) or (side == 'bp1' and (r.lg, r.part[:-1] + chr(ord(letter.group(1)) - 1)) in tandem))
            name = f'inv{r.lg[2:]}_internal' if internal else f'inv{r.part}_{side}'
            if name in seen:
                continue
            seen.add(name)
            rows.append(dict(chrom=r.chrom, start=max(0, int(s) - pad), end=int(e) + pad, name=name))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'name'])
