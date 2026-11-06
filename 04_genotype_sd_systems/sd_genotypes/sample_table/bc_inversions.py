"""Infer inversion genotypes of YH x MC backcross fish from the ancestry (haplotype) calls, for inversions where the F1 sires are heterozygous.

Logic: in an inversion heterozygote recombination inside the inversion is suppressed, so each BC fish carries ONE F1-sire haplotype across the whole
block. State 1 = the YH-derived haplotype, state 0 = the MC haplotype. If MC is homozygous non-inverted (Kumar 2025 Table 1 and Blumer 2025 Table S7:
Mchenga is non-inverted on chr2 and chr13) and the F1 sire is heterozygous (shown by the missing recombination), the YH haplotype is the inverted one:
BC genotype = 1 (heterozygous) if the fish carries the YH haplotype, 0 (homozygous non-inverted) if it carries the MC haplotype. Coding as Blumer Table S7.
The block boundaries are detected from the data as the longest run of grid intervals with (almost) no crossovers in any fish; they can include
adjacent low-recombination sequence (e.g. a centromere), so they are outer bounds of the inversion, not breakpoints."""
import numpy as np
import pandas as pd


def load_grid(path):
    g = pd.read_csv(path)
    g = g.drop(columns=[c for c in g.columns if str(c).startswith('Unnamed')])
    assert {'Chromosome', 'Position'} <= set(g.columns), 'expected columns Chromosome, Position, <one column per fish>'
    return g


def detect_block(pos, R, max_rate=0.03, min_len=10e6, slack=2):
    """Longest stretch of consecutive grid intervals whose crossover count (switches summed over all fish) is a small fraction of normal.
    A stretch [a, b] qualifies when  sum(switches) <= max_rate * n_intervals + slack.  With N fish and ~1 cM/Mb a normal 100 kb interval has
    ~N/1000 crossovers, so max_rate should be ~10% of that (N=266: normal ~0.27, default 0.03). Returns (lo, hi, n_inside, n_total) or None."""
    sw = (np.diff(R, axis=0) != 0).sum(1)
    cs = np.concatenate([[0], np.cumsum(sw)])
    n = len(sw)
    best = (0, 0, -1)
    for a in range(n):
        for b in range(n - 1, a + best[0] - 1, -1):               # only lengths that could beat the current best
            if cs[b + 1] - cs[a] <= max_rate * (b - a + 1) + slack:
                best = (b - a + 1, a, b)
                break
    ln, a, b = best
    if ln == 0 or pos[b + 1] - pos[a] < min_len:
        return None
    return int(pos[a]), int(pos[b + 1]), int(cs[b + 1] - cs[a]), int(sw.sum())


def infer(grid, chrom_by_lg, mc_inverted, max_rate=0.03, min_len=10e6, min_consistency=0.95):
    """grid: DataFrame(Chromosome, Position, fish...). chrom_by_lg: {'LG2': 'NC_...'}; mc_inverted: {'LG2': False, ...} (MC dam and MC haplotype state).
    Returns (per-fish table, per-LG block table). Fish with consistency < min_consistency across the block get NaN (ambiguous)."""
    fish = [c for c in grid.columns if c not in ('Chromosome', 'Position')]
    calls, blocks = pd.DataFrame(index=fish), []
    for lg, ch in chrom_by_lg.items():
        sub = grid[grid.Chromosome == ch].sort_values('Position')
        pos, R = sub.Position.values, np.round(sub[fish].values).astype(int)
        blk = detect_block(pos, R, max_rate, min_len)
        if blk is None:
            blocks.append(dict(LG=lg, chrom=ch, found=False)); continue
        lo, hi, n_in, n_all = blk
        inside = R[(pos >= lo) & (pos <= hi)]
        frac_yh = inside.mean(0)
        cons = np.maximum(frac_yh, 1 - frac_yh)
        haplo = (frac_yh >= 0.5).astype(int)                       # 1 = YH-derived haplotype
        mc_inv = int(bool(mc_inverted[lg]))
        yh_inv = 1 - mc_inv                                         # the F1 sire is heterozygous, so YH carries the other state
        copies = mc_inv + np.where(haplo == 1, yh_inv, mc_inv)      # dam contributes mc_inv; sire contributes the haplotype's state
        copies = np.where(cons >= min_consistency, copies, np.nan)
        n_num = lg.replace('LG', '')
        calls[f'Inversion{n_num}'] = copies
        calls[f'haplotype_{lg}'] = np.where(haplo == 1, 'YH', 'MC')
        calls[f'consistency_{lg}'] = cons
        blocks.append(dict(LG=lg, chrom=ch, found=True, block_start=lo, block_end=hi, block_Mb=(hi - lo) / 1e6, switches_inside=n_in,
                           switches_chromosome=n_all, n_fish=len(fish), n_YH_haplotype=int(haplo.sum()), n_ambiguous=int((cons < min_consistency).sum())))
    return calls, pd.DataFrame(blocks)


def apply_to_table(table, calls, lgs):
    """Fill Inversion2/13 (etc.) for BC fish that have no Blumer genotype, and mark the source. Returns the new table."""
    t = table.copy()
    inv_cols = [f'Inversion{lg.replace("LG", "")}' for lg in lgs if f'Inversion{lg.replace("LG", "")}' in calls.columns]
    for c in inv_cols:
        m = t.SampleID.isin(calls.index) & t[c].isna()
        t.loc[m, c] = t.loc[m, 'SampleID'].map(calls[c])
    filled = t.SampleID.isin(calls.index) & t[inv_cols].notna().any(axis=1) & (t.InversionSource == '')
    t.loc[filled, 'InversionSource'] = 'inferred_BC_haplotype:' + ','.join(c.replace('Inversion', 'chr') for c in inv_cols)
    return t


def check_no_block(grid, chrom_by_lg, max_rate=0.03, min_len=10e6):
    """For chromosomes where the F1 sires are expected NOT to be heterozygous: True if the data show no recombination-suppressed block."""
    fish = [c for c in grid.columns if c not in ('Chromosome', 'Position')]
    out = {}
    for lg, ch in chrom_by_lg.items():
        sub = grid[grid.Chromosome == ch].sort_values('Position')
        out[lg] = detect_block(sub.Position.values, np.round(sub[fish].values).astype(int), max_rate, min_len) is None
    return out


def is_bc_row(table, bc_species='YHxMC x MC BC', prefix='MCYHBC1-'):
    return (table.Species == bc_species) | table.SampleID.astype(str).str.startswith(prefix)


def apply_expected(table, expected, no_block, lg_of, bc_species='YHxMC x MC BC'):
    """Fill pedigree-expected genotypes for BC fish where blank. expected: {'Inversion9': 2, ...}; no_block: {'LG9': True, ...} (data support);
    lg_of: {'Inversion9': 'LG9', ...}. InversionSource gets '...;pedigree_expected:chr9,chr10,chr11' for the rows that received values.
    Returns (table, list of columns filled)."""
    t = table.copy()
    bc = is_bc_row(t, bc_species)
    filled, per_row = [], {}
    for col, val in expected.items():
        if not no_block.get(lg_of[col], False):
            print(f'[bc_inversions] {col}: a recombination-suppressed block was found on {lg_of[col]}; NOT filling the pedigree expectation')
            continue
        m = bc & t[col].isna()
        t.loc[m, col] = val
        for i in t.index[m]:
            per_row.setdefault(i, []).append(col.replace('Inversion', 'chr'))
        filled.append(col)
    for i, cols in per_row.items():
        cur = t.at[i, 'InversionSource']
        t.at[i, 'InversionSource'] = (cur + ';' if cur else '') + 'pedigree_expected:' + ','.join(cols)
    return t, filled
