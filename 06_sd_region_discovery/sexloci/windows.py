"""Helpers for grouping sorted markers into fixed-size genomic windows."""
from __future__ import annotations

import numpy as np


def window_groups(chrom: np.ndarray, pos: np.ndarray, size: int):
    """Return (starts, ends, keys) for runs of consecutive markers sharing (chrom, window index).

    keys is a list of (chrom, widx) tuples. Input must be position-sorted within contigs
    (as in any valid VCF); non-contiguous repeats of a key are harmless because callers
    merge by key.
    """
    m = len(pos)
    if m == 0:
        return np.array([], int), np.array([], int), []
    widx = (pos - 1) // size
    change = np.ones(m, bool)
    change[1:] = (chrom[1:] != chrom[:-1]) | (widx[1:] != widx[:-1])
    starts = np.flatnonzero(change)
    ends = np.append(starts[1:], m)
    keys = [(str(chrom[s]), int(widx[s])) for s in starts]
    return starts, ends, keys


def robust_z(x) -> np.ndarray:
    """(x - median) / (1.4826 * MAD); if MAD is 0 (many identical values) fall back to the mean absolute
    deviation so a single outlier among identical values is still detected."""
    x = np.asarray(x, float)
    med = np.nanmedian(x)
    dev = np.abs(x - med)
    scale = np.nanmedian(dev) * 1.4826
    if not np.isfinite(scale) or scale == 0:
        scale = np.nanmean(dev) * 1.2533
    if not np.isfinite(scale) or scale == 0:
        return np.zeros_like(x)
    return (x - med) / scale
