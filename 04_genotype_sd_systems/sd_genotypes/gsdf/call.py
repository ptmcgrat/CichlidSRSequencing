"""Copy number from the depth ratio, with the flags of the method (config/gsdf.yaml)."""
import math
import numpy as np
import pandas as pd


def estimate_copies(ratio, cal):
    return 2 + (ratio - cal['ratio_at_two_copies']) / cal['ratio_per_copy']


def call_one(m, cfg):
    """m = measure.metrics(...). Returns dict: copies_estimate, gsdf_number (rounded, clipped), genotype label, flags (the method's: low depth, non-integer, outside the calibrated range) and
    check_flags (the one-sided junction cross-check: a junction in a sample called with two copies; it does not change whether the call counts as unflagged)."""
    th, cal = cfg['thresholds'], cfg['calibration']
    flags, checks = [], []
    if not (m['flank_depth'] == m['flank_depth']) or m['flank_depth'] < th['min_flank_depth']:
        flags.append('low_flank_depth')
    est = estimate_copies(m['ratio'], cal) if m['ratio'] == m['ratio'] else float('nan')
    if est != est:
        return dict(copies_estimate=np.nan, gsdf_number=np.nan, genotype='', flags=flags + ['no_estimate'], check_flags=[])
    c_lo, c_hi = th['clip_copies']; e_lo, e_hi = th['estimate_range']
    n = int(min(max(round(est), c_lo), c_hi))                                   # round(): half to even, as pandas does in sexloci
    if abs(est - round(est)) > th['max_distance_from_integer']:
        flags.append('non_integer')
    if est < e_lo or est > e_hi:
        flags.append('outside_calibrated_range')
    jr = m['junction_ratio']                                                # one-sided: reads across the junction between the copies exist only on chromosomes that carry the tandem duplication
    if jr == jr and n == 2 and e_lo <= est <= e_hi and jr > th['junction_max_two_copy']:
        checks.append('junction_in_two_copy')
    return dict(copies_estimate=est, gsdf_number=n, genotype=cfg['labels'].get(n, f'CN{n}'), flags=flags, check_flags=checks)


def fit_anchors(ratios, copies):
    """Calibration from samples of KNOWN copy number (sexloci's --anchors): ratio = a + step * copies; returns {'ratio_at_two_copies': a + 2 * step, 'ratio_per_copy': step}. Needs at least two different copy numbers."""
    import pandas as pd
    r = pd.Series(ratios).astype(float); c = pd.Series(copies).astype(float)
    ok = r.notna() & c.notna()
    if c[ok].nunique() < 2:
        raise ValueError('anchors need samples with at least two different copy numbers')
    step, icpt = np.polyfit(c[ok].values, r[ok].values, 1)
    return dict(ratio_at_two_copies=float(icpt + 2 * step), ratio_per_copy=float(step))


THEORY_CALIBRATION = dict(ratio_at_two_copies=0.5, ratio_per_copy=0.25)         # ratio = copies / 4 (sexloci's --theory)


def call_table(raw_rows, cfg, calibration=None):
    """raw_rows: iterable of dicts with SampleID, bins, bins_q, junction_depth, junction_depth_q. Returns one row per sample. `calibration` (e.g. from fit_anchors) replaces the one in the config."""
    from . import measure
    if calibration is not None:
        cfg = dict(cfg, calibration=calibration)
    out = []
    for raw in raw_rows:
        m = measure.metrics(raw, cfg); c = call_one(m, cfg)
        out.append(dict(SampleID=raw['SampleID'], **{k: m[k] for k in ('unit_depth', 'flank_depth', 'ratio', 'flank_depth_unique', 'junction_depth_unique', 'junction_ratio')},
                        copies_estimate=c['copies_estimate'], gsdf_number=c['gsdf_number'], genotype=c['genotype'], call_flags=';'.join(c['flags']), check_flags=';'.join(c['check_flags'])))
    return pd.DataFrame(out)


def cluster_check(calls, cfg):
    """Where the unflagged 2-, 3- and 4-copy samples actually sit, against the calibration (ratio 0.530 at two copies, +0.268 per copy): the drift check for a new batch."""
    cal = cfg['calibration']; ok = calls[calls['call_flags'] == '']      # NOT calls.flags: that is a pandas attribute
    rows = []
    for n in (2, 3, 4):
        g = ok[ok.gsdf_number == n]
        rows.append(dict(copies=n, n=len(g), median_ratio=float(g.ratio.median()) if len(g) else np.nan, expected_ratio=cal['ratio_at_two_copies'] + (n - 2) * cal['ratio_per_copy'],
                         median_junction_ratio=float(g.junction_ratio.median()) if len(g) else np.nan))
    return pd.DataFrame(rows)
