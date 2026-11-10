"""Population-genetic contrasts between sexes: windowed Weir-Cockerham Fst, heterozygosity difference,
sex-diagnostic SNP counts, and male:female read-depth ratio.

These are the signals used in the Aulonocara analysis (Fst + delta-heterozygosity) plus two extras that
help recognise highly diverged or partly deleted sex chromosomes.
  delta-het = het(male) - het(female):  > 0 suggests XY, < 0 suggests ZW
  coverage log2(M/F):                   > 0 -> X-like / Z-like dosage differences or male-specific sequence
                                        (interpret with care: derived from variant-site depth)
"""
from __future__ import annotations

import warnings
from typing import Dict, Tuple

import numpy as np
import pandas as pd

from .windows import window_groups


def wc_fst_components(D1: np.ndarray, D2: np.ndarray):
    """Weir & Cockerham (1984) a, b, c for two populations from dosage matrices (m x n_i, -1 = missing).

    Returns (a, b, c, n1, n2, p1, p2, h1, h2); rows with <2 called individuals in a group are NaN.
    """
    c1, c2 = D1 >= 0, D2 >= 0
    n1, n2 = c1.sum(1).astype(float), c2.sum(1).astype(float)
    p1 = np.where(c1, D1, 0).sum(1) / np.maximum(2 * n1, 1)
    p2 = np.where(c2, D2, 0).sum(1) / np.maximum(2 * n2, 1)
    h1 = (D1 == 1).sum(1) / np.maximum(n1, 1)
    h2 = (D2 == 1).sum(1) / np.maximum(n2, 1)
    nbar = (n1 + n2) / 2
    ok = (n1 >= 2) & (n2 >= 2)
    nbar_s = np.where(ok, nbar, 2.0)
    nc = 2 * nbar_s - (n1 ** 2 + n2 ** 2) / (2 * nbar_s)
    pbar = (n1 * p1 + n2 * p2) / (2 * nbar_s)
    s2 = (n1 * (p1 - pbar) ** 2 + n2 * (p2 - pbar) ** 2) / nbar_s
    hbar = (n1 * h1 + n2 * h2) / (2 * nbar_s)
    a = nbar_s / nc * (s2 - (pbar * (1 - pbar) - s2 / 2 - hbar / 4) / (nbar_s - 1))
    b = nbar_s / (nbar_s - 1) * (pbar * (1 - pbar) - s2 / 2 - (2 * nbar_s - 1) / (4 * nbar_s) * hbar)
    c = hbar / 2
    for arr in (a, b, c):
        arr[~ok] = np.nan
    return a, b, c, n1, n2, p1, p2, h1, h2


class PopgenAccumulator:
    name = "popgen"

    def __init__(self, unit, cfg):
        self.cfg = cfg
        ci = unit.contrast_idx
        self.m, self.f = ci[unit.sex[ci] == 1], ci[unit.sex[ci] == 0]
        self.ci = ci
        self.rows: Dict[Tuple[str, int], np.ndarray] = {}   # sums: a, abc, hM, hF, n, diagXY, diagZW
        self.cov_sum: Dict[Tuple[str, int], np.ndarray] = {}
        self.cov_cnt: Dict[Tuple[str, int], np.ndarray] = {}
        self.n_samples = len(unit.names)
        self.have_dp = False

    def update(self, chunk):
        cfg = self.cfg
        # ---- coverage (raw depth, all samples)
        if chunk.dp is not None:
            self.have_dp = True
            dp = chunk.dp
            valid = ~np.isnan(dp)
            starts, ends, keys = window_groups(chunk.chrom, chunk.pos, cfg.window_size)
            S = np.add.reduceat(np.where(valid, dp, 0.0), starts, axis=0)
            C = np.add.reduceat(valid.astype(np.int32), starts, axis=0)
            for i, k in enumerate(keys):
                self.cov_sum[k] = self.cov_sum.get(k, 0) + S[i]
                self.cov_cnt[k] = self.cov_cnt.get(k, 0) + C[i]

        # ---- genotype contrasts
        Dm, Df = chunk.dos[:, self.m], chunk.dos[:, self.f]
        a, b, c, n1, n2, p1, p2, h1, h2 = wc_fst_components(Dm, Df)
        ok = (n1 >= cfg.min_per_sex_called) & (n2 >= cfg.min_per_sex_called)
        if not ok.any():
            return
        tol = cfg.diag_tol
        maf2 = np.minimum(p2, 1 - p2)
        maf1 = np.minimum(p1, 1 - p1)
        diag_xy = (h1 >= 1 - tol) & (h2 <= tol) & (maf2 <= tol / 2)
        diag_zw = (h2 >= 1 - tol) & (h1 <= tol) & (maf1 <= tol / 2)
        abc = a + b + c
        vals = np.column_stack([
            np.nan_to_num(a), np.nan_to_num(abc), h1, h2, np.ones(len(a)),
            diag_xy.astype(float), diag_zw.astype(float)])
        vals = vals[ok]
        chrom, pos = chunk.chrom[ok], chunk.pos[ok]
        starts, ends, keys = window_groups(chrom, pos, cfg.window_size)
        R = np.add.reduceat(vals, starts, axis=0)
        for i, k in enumerate(keys):
            self.rows[k] = self.rows.get(k, 0) + R[i]

    def finalize(self):
        cfg = self.cfg
        if not self.rows:
            return pd.DataFrame()
        keys = sorted(self.rows, key=lambda k: (k[0], k[1]))
        R = np.vstack([self.rows[k] for k in keys])
        df = pd.DataFrame({"chrom": [k[0] for k in keys], "widx": [k[1] for k in keys]})
        df["pg_n_markers"] = R[:, 4].astype(int)
        with np.errstate(invalid="ignore", divide="ignore"):
            df["fst"] = np.where(R[:, 1] > 0, R[:, 0] / R[:, 1], np.nan)
            df["het_male"] = R[:, 2] / R[:, 4]
            df["het_female"] = R[:, 3] / R[:, 4]
        df["dhet"] = df["het_male"] - df["het_female"]
        df["n_diag_xy"] = R[:, 5].astype(int)
        df["n_diag_zw"] = R[:, 6].astype(int)
        small = df["pg_n_markers"] < cfg.min_window_markers
        df.loc[small, ["fst", "dhet"]] = np.nan

        if self.have_dp and self.cov_sum:
            ck = sorted(self.cov_sum, key=lambda k: (k[0], k[1]))
            S = np.vstack([np.asarray(self.cov_sum[k]) for k in ck]).astype(float)
            C = np.vstack([np.asarray(self.cov_cnt[k]) for k in ck]).astype(float)
            with np.errstate(invalid="ignore", divide="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)       # windows with no depth data are all-NaN rows
                mean = np.where(C >= 3, S / C, np.nan)
                mean = mean / np.nanmedian(mean, axis=0, keepdims=True)
                lm = np.nanmean(mean[:, self.m], axis=1)
                lf = np.nanmean(mean[:, self.f], axis=1)
                ratio = np.log2(lm / lf)
            cov = pd.DataFrame({"chrom": [k[0] for k in ck], "widx": [k[1] for k in ck],
                                "cov_log2_mf": ratio})
            df = df.merge(cov, on=["chrom", "widx"], how="left")
        return df
