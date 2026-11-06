"""Copy number from the depth ratio, with the flags of the method (config/gsdf.yaml)."""
import math
import numpy as np
import pandas as pd


def estimate_copies(ratio, cal):
    return 2 + (ratio - cal['ratio_at_two_copies']) / cal['ratio_per_copy']


def call_one(m, cfg):
    """m = measure.metrics(...). Returns dict: copies_estimate, gsdf_number (rounded), genotype label, flags (list), junction fields."""
    th, cal = cfg['thresholds'], cfg['calibration']
    flags = []
    if not (m['flank_depth'] == m['flank_depth']) or m['flank_depth'] < th['min_flank_depth']:
        flags.append('low_flank_depth')
    est = estimate_copies(m['ratio'], cal) if m['ratio'] == m['ratio'] else float('nan')
    if est != est:
        return dict(copies_estimate=np.nan, gsdf_number=np.nan, genotype='', flags=flags + ['no_estimate'], junction_implied_copies=np.nan)
    n = int(math.floor(est + 0.5))
    if abs(est - n) > th['max_distance_from_integer']:
        flags.append('non_integer')
    lo, hi = th['copies_range']
    if n < lo or n > hi:
        flags.append('outside_calibrated_range')
    jr = m['junction_ratio']
    n_dup_j = 2 * jr if jr == jr else float('nan')                           # duplicated chromosomes the unique junction reads imply: 0, 1 or 2
    implied = 2 + min(2, max(0, round(n_dup_j))) if n_dup_j == n_dup_j else np.nan
    if n_dup_j == n_dup_j and lo <= n <= hi and abs(n_dup_j - (n - 2)) > th['junction_max_distance']:
        flags.append('junction_disagrees')
    return dict(copies_estimate=est, gsdf_number=n, genotype=cfg['labels'].get(n, ''), flags=flags, junction_implied_copies=implied)


def call_table(raw_rows, cfg):
    """raw_rows: iterable of dicts with SampleID, bins, bins_q, junction_depth, junction_depth_q. Returns one row per sample."""
    from . import measure
    out = []
    for raw in raw_rows:
        m = measure.metrics(raw, cfg); c = call_one(m, cfg)
        out.append(dict(SampleID=raw['SampleID'], **{k: m[k] for k in ('unit_depth', 'flank_depth', 'ratio', 'flank_depth_unique', 'junction_depth_unique', 'junction_ratio')},
                        copies_estimate=c['copies_estimate'], gsdf_number=c['gsdf_number'], genotype=c['genotype'], junction_implied_copies=c['junction_implied_copies'], call_flags=';'.join(c['flags'])))
    return pd.DataFrame(out)


def cluster_check(calls, cfg):
    """Where the unflagged 2-, 3- and 4-copy samples actually sit, against the calibration (ratio 0.530 at two copies, +0.268 per copy): the drift check for a new batch."""
    cal = cfg['calibration']; ok = calls[calls['call_flags'] == '']      # NOT calls.flags: that is a pandas attribute
    rows = []
    for n in (2, 3, 4):
        g = ok[ok.gsdf_number == n]
        rows.append(dict(copies=n, n=len(g), median_ratio=float(g.ratio.median()) if len(g) else np.nan, expected_ratio=cal['ratio_at_two_copies'] + (n - 2) * cal['ratio_per_copy'],
                         median_junction_ratio=float(g.junction_ratio.median()) if len(g) else np.nan, expected_junction=(n - 2) / 2))
    return pd.DataFrame(rows)
