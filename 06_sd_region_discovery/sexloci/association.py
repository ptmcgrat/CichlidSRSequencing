"""Marker-level sex association, computed per stratum and combined, with permutation FWER control.

Strata are the independent units of replication: each full-sib family (offspring only) and, if present,
one pool of non-family animals. Within a stratum two statistics are computed per SNP:
  * additive dosage correlation with sex        (fixed / differentiated alleles)
  * heterozygote-indicator correlation with sex  (one sex heterozygous, the other homozygous: XY / ZW)
stratum chi2 = n_s * max(r_dosage^2, r_het^2). Squared correlations are used because the allele linked
to sex differs between families (phase), so signed effects must not be pooled. The combined statistic is
the sum over strata, so a SNP supported in several independent strata outranks one supported in a single
stratum. Sex labels are permuted within each stratum, preserving kinship structure.
"""
from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
import pandas as pd

from .stats import centre_unit, make_permuted_labels, maxz_adjust, standardise
from .windows import window_groups

BLOCK_MARKERS = 2000


class AssociationAccumulator:
    name = "assoc"

    def __init__(self, unit, cfg):
        self.cfg = cfg
        self.idx = unit.contrast_idx
        lab = unit.strata[self.idx]
        sex = unit.sex[self.idx]
        rng = np.random.default_rng(cfg.seed)
        self.strata, self.skipped = [], []
        for s in np.unique(lab):
            loc = np.flatnonzero(lab == s)
            y = sex[loc]
            if len(loc) < cfg.min_stratum_size or min((y == 1).sum(), (y == 0).sum()) < 2:
                self.skipped.append(s)
                continue
            Y = centre_unit(make_permuted_labels(y, np.zeros(len(y)), cfg.n_perm, rng))
            self.strata.append(dict(label=s, loc=loc, n=len(loc), Y=Y, male=(y == 1)))
        if not self.strata:
            raise ValueError("No stratum has enough animals of both sexes for the association test")
        self.K = len(self.strata)
        self.win: Dict[Tuple[str, int], np.ndarray] = {}
        self.win_n: Dict[Tuple[str, int], int] = {}
        self.win_sig: Dict[Tuple[str, int], int] = {}
        self.hits = []
        self.n_used = sum(s["n"] for s in self.strata)

    def update(self, chunk):
        D_all = chunk.dos[:, self.idx]
        called_all = D_all >= 0
        nc = called_all.sum(1)
        alt = np.where(called_all, D_all, 0).sum(1)
        keep = np.minimum(alt, 2 * nc - alt) >= self.cfg.min_mac
        if not keep.any():
            return
        D_all, chrom, pos = D_all[keep], chunk.chrom[keep], chunk.pos[keep]
        for b in range(0, len(pos), BLOCK_MARKERS):
            self._block(D_all[b:b + BLOCK_MARKERS], chrom[b:b + BLOCK_MARKERS], pos[b:b + BLOCK_MARKERS])

    def _block(self, D, chrom, pos):
        cfg = self.cfg
        m, P1 = len(pos), cfg.n_perm + 1
        S_tot = np.zeros((m, P1), np.float32)
        S_k, obs_r = [], []
        for st in self.strata:
            Ds = D[:, st["loc"]]
            called = Ds >= 0
            R1 = standardise(Ds.astype(np.float32), called) @ st["Y"]
            R2 = standardise((Ds == 1).astype(np.float32), called) @ st["Y"]
            S = (st["n"] * np.maximum(R1 ** 2, R2 ** 2)).astype(np.float32)
            S_k.append(S)
            S_tot += S
        starts, ends, keys = window_groups(chrom, pos, cfg.window_size)
        Mt = np.maximum.reduceat(S_tot, starts, axis=0)
        Mk = [np.maximum.reduceat(S, starts, axis=0) for S in S_k]
        for i, key in enumerate(keys):
            arr = np.stack([Mt[i]] + [mk[i] for mk in Mk])
            cur = self.win.get(key)
            self.win[key] = arr if cur is None else np.maximum(cur, arr)
            self.win_n[key] = self.win_n.get(key, 0) + int(ends[i] - starts[i])
        sig = np.flatnonzero(S_tot[:, 0] >= cfg.store_min_chi2)
        for k, s, e in zip(keys, starts, ends):
            c = int((S_tot[s:e, 0] >= cfg.store_min_chi2).sum())
            if c:
                self.win_sig[k] = self.win_sig.get(k, 0) + c
        male_all = None
        for j in sig:
            row = D[j]
            nom = int(sum(S[j, 0] >= 3.84 for S in S_k))
            hm = hf = fm = ff = np.nan
            ms = np.concatenate([st["loc"][st["male"]] for st in self.strata])
            fs = np.concatenate([st["loc"][~st["male"]] for st in self.strata])
            cm, cf = row[ms] >= 0, row[fs] >= 0
            if cm.any():
                hm, fm = float((row[ms][cm] == 1).mean()), float(row[ms][cm].mean() / 2)
            if cf.any():
                hf, ff = float((row[fs][cf] == 1).mean()), float(row[fs][cf].mean() / 2)
            self.hits.append((chrom[j], int(pos[j]), float(S_tot[j, 0]), nom, hm, hf, fm, ff))

    def finalize(self):
        cfg = self.cfg
        if not self.win:
            return pd.DataFrame(), pd.DataFrame(), {}, {}
        keys = sorted(self.win, key=lambda k: (k[0], k[1]))
        W = np.stack([self.win[k] for k in keys])            # (windows, K+1, P+1)
        tot = W[:, 0, :]
        perm_max = tot[:, 1:].max(axis=0)
        padj = (1 + (perm_max[None, :] >= tot[:, [0]]).sum(1)) / (cfg.n_perm + 1)
        thr = float(np.quantile(perm_max, 1 - cfg.alpha))
        win = pd.DataFrame({
            "chrom": [k[0] for k in keys], "widx": [k[1] for k in keys],
            "assoc_n_markers": [self.win_n[k] for k in keys],
            "assoc_n_stored": [self.win_sig.get(k, 0) for k in keys],
            "assoc_max_chi2": tot[:, 0], "assoc_padj": padj})
        tables = {}
        for k, st in enumerate(self.strata):
            p_un, p_adj = maxz_adjust(W[:, k + 1, :])
            tables[st["label"]] = pd.DataFrame({
                "chrom": win["chrom"], "widx": win["widx"], "chi2": W[:, k + 1, 0], "p": p_un, "padj": p_adj})
        hits = pd.DataFrame(self.hits, columns=["chrom", "pos", "chi2", "n_strata_nominal", "het_frac_male",
                                                "het_frac_female", "alt_freq_male", "alt_freq_female"])
        if len(hits):
            hits["padj"] = (1 + (perm_max[None, :] >= hits["chi2"].to_numpy()[:, None]).sum(1)) / (cfg.n_perm + 1)
            hits = hits.sort_values("chi2", ascending=False).reset_index(drop=True)
        info = {"threshold_chi2": thr, "strata_used": [s["label"] for s in self.strata],
                "strata_skipped": self.skipped, "n_samples": self.n_used}
        return win, hits, info, tables
