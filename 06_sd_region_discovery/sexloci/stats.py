"""Shared statistics: standardisation, stratified permutations, max-Z family-wise correction."""
from __future__ import annotations

import numpy as np


def standardise(X: np.ndarray, called: np.ndarray) -> np.ndarray:
    """Row-wise mean-impute, centre and scale to unit L2 norm (constant rows -> 0)."""
    cnt = called.sum(1, keepdims=True)
    mean = (X * called).sum(1, keepdims=True) / np.maximum(cnt, 1)
    Xc = np.where(called, X - mean, 0.0).astype(np.float32)
    norm = np.sqrt((Xc ** 2).sum(1, keepdims=True))
    return np.divide(Xc, norm, out=np.zeros_like(Xc), where=norm > 1e-6)


def make_permuted_labels(y, strata, n_perm: int, rng) -> np.ndarray:
    """(n, n_perm+1) float32; column 0 observed, others permuted within strata."""
    y = np.asarray(y, np.float32)
    Y = np.empty((len(y), n_perm + 1), np.float32)
    Y[:, 0] = y
    strata = np.asarray(strata)
    for lab in np.unique(strata):
        g = np.flatnonzero(strata == lab)
        Y[g, 1:] = rng.permuted(np.tile(y[g], (n_perm, 1)), axis=1).T
    return Y


def centre_unit(Y: np.ndarray) -> np.ndarray:
    Yc = Y - Y.mean(0, keepdims=True)
    norm = np.sqrt((Yc ** 2).sum(0, keepdims=True))
    return np.divide(Yc, norm, out=np.zeros_like(Yc), where=norm > 1e-9)


def maxz_adjust(T: np.ndarray):
    """T: (tests, P+1), column 0 observed. Westfall-Young max-T on per-test standardised statistics.

    Returns (p_unadjusted, p_adjusted) for the observed column. Unlike a min-p scheme this does not
    saturate at 1/(P+1) when many windows reach the top rank.
    """
    P = T.shape[1] - 1
    perm = T[:, 1:]
    mu = perm.mean(1, keepdims=True)
    sd = perm.std(1, keepdims=True)
    sd = np.where(sd < 1e-9, np.inf, sd)
    Z = (T - mu) / sd
    max_perm = Z[:, 1:].max(axis=0)
    p_adj = (1 + (max_perm[None, :] >= Z[:, [0]]).sum(1)) / (P + 1)
    p_un = (1 + (perm >= T[:, [0]]).sum(1)) / (P + 1)
    return p_un, p_adj


def perm_chi2_pvalue(labels_a, labels_b, n_perm: int = 2000, seed: int = 1) -> float:
    """Permutation p-value for association between two categorical vectors (chi-square statistic)."""
    a, b = np.asarray(labels_a), np.asarray(labels_b)
    ua, ia = np.unique(a, return_inverse=True)
    ub, ib = np.unique(b, return_inverse=True)

    def stat(ia_, ib_):
        tab = np.zeros((len(ua), len(ub)))
        np.add.at(tab, (ia_, ib_), 1)
        exp = tab.sum(1, keepdims=True) * tab.sum(0, keepdims=True) / tab.sum()
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.nansum((tab - exp) ** 2 / np.where(exp > 0, exp, np.nan))

    obs = stat(ia, ib)
    rng = np.random.default_rng(seed)
    ge = sum(stat(ia, rng.permutation(ib)) >= obs - 1e-9 for _ in range(n_perm))
    return (1 + ge) / (n_perm + 1)
