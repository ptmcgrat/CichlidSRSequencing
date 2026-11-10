"""Pedigree-transmission mapping of sex (generalises the Mendelian-phasing / haplotype-painting approach).

For each full-sib family and each parent role we use *informative markers*: parent X heterozygous and the
other parent homozygous. At such a marker the allele an offspring received from parent X is directly
observable (offspring dosage minus the contribution of the homozygous parent), so no phasing is needed.

Within a window the inheritance vectors of the informative markers all track one thing - which of the
parent's two haplotypes each offspring received - up to an arbitrary per-marker sign (we do not know
which allele sits on which haplotype). The leading principal component of the marker x offspring matrix
is sign-invariant and gives one consensus 'haplotype received' score per offspring. We then ask whether
this score predicts offspring sex:   T = sum over families of  n * r^2.

  * sire role significant, dam not  -> male heterogamety (XY-like): sons inherit one sire haplotype
  * dam role significant, sire not  -> female heterogamety (ZW-like)

Significance uses within-family permutation of offspring sex with max-T (Westfall-Young) correction across all windows
and both roles, so inversions are neither required nor assumed.
"""
from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
import pandas as pd

from .stats import centre_unit, make_permuted_labels, maxz_adjust
from .windows import window_groups


def _leading_score(X: np.ndarray, min_markers: int, min_valid_frac: float):
    """Sign-invariant consensus score per offspring from a (markers x offspring) 0/1/NaN matrix."""
    valid = ~np.isnan(X)
    X = X[valid.mean(1) >= min_valid_frac]
    if X.shape[0] < min_markers:
        return None, 0
    valid = ~np.isnan(X)
    mean = np.where(valid, X, 0).sum(1, keepdims=True) / valid.sum(1, keepdims=True)
    Xc = np.where(valid, X - mean, 0.0)
    keep = (Xc ** 2).sum(1) > 0
    Xc = Xc[keep]
    if Xc.shape[0] < min_markers:
        return None, 0
    _, s, vt = np.linalg.svd(Xc, full_matrices=False)
    return vt[0], Xc.shape[0]


class PedigreeAccumulator:
    name = "pedigree"

    def __init__(self, unit, cfg):
        self.cfg = cfg
        self.unit = unit
        self.tracks = {}            # (family id, role) -> {chrom: (window indices, scores (windows x offspring))}
        self.family_meta = {}       # family id -> dict(names, sex)
        self.rng = np.random.default_rng(cfg.seed + 1)
        self.fams = []
        self.Y = []
        self.skipped = []
        for f in unit.families:
            self.fams.append(f)
            Y = make_permuted_labels(f.sex.astype(np.float32), np.zeros(len(f.sex)), cfg.n_perm, self.rng)
            self.Y.append(centre_unit(Y))
        self.store: Dict[Tuple[int, str, str, int], list] = {}
        self.count: Dict[Tuple[int, str, str, int], int] = {}
        self.n_informative = {}

    @property
    def active(self):
        return len(self.fams) > 0

    def _append(self, key, block):
        cap = self.cfg.ped_max_markers_per_window
        lst = self.store.setdefault(key, [])
        lst.append(block)
        self.count[key] = self.count.get(key, 0) + block.shape[0]
        if self.count[key] > 2 * cap:
            X = np.concatenate(lst, axis=0)
            sel = np.sort(self.rng.choice(X.shape[0], cap, replace=False))
            self.store[key] = [X[sel]]
            self.count[key] = cap

    def update(self, chunk):
        cfg = self.cfg
        G = chunk.dos_role if chunk.dos_role is not None else chunk.dos
        for fi, f in enumerate(self.fams):
            s, d = G[:, f.sire], G[:, f.dam]
            O = G[:, f.offspring]
            for role, het, hom in (("sire", s, d), ("dam", d, s)):
                mask = (het == 1) & ((hom == 0) | (hom == 2))
                if not mask.any():
                    continue
                Om = O[mask]
                x = np.where(Om >= 0, Om.astype(np.float32) - (hom[mask] // 2)[:, None], np.nan)
                x[~((x == 0) | (x == 1))] = np.nan      # Mendelian errors -> missing
                starts, ends, keys = window_groups(chunk.chrom[mask], chunk.pos[mask], cfg.window_size)
                for st, en, (chrom, w) in zip(starts, ends, keys):
                    self._append((fi, role, chrom, w), x[st:en])
                self.n_informative[(fi, role)] = self.n_informative.get((fi, role), 0) + int(mask.sum())

    def finalize(self):
        cfg = self.cfg
        acc: Dict[Tuple[str, str, int], dict] = {}
        fam_rows = {}
        raw_tracks = {}
        for (fi, role, chrom, w), blocks in self.store.items():
            X = np.concatenate(blocks, axis=0)
            u, nmk = _leading_score(X, cfg.ped_min_markers, cfg.ped_min_valid_frac)
            if u is None:
                continue
            u = u - u.mean()
            nrm = np.sqrt((u ** 2).sum())
            if nrm < 1e-9:
                continue
            r = (u / nrm) @ self.Y[fi]                  # (P+1,)
            n = len(u)
            T = n * r ** 2
            fam_rows.setdefault(fi, []).append((role, chrom, w, T))
            raw_tracks.setdefault((self.fams[fi].id, role), {}).setdefault(chrom, []).append((w, u / nrm))
            e = acc.setdefault((role, chrom, w), {"T": np.zeros(cfg.n_perm + 1), "k": 0, "r": [], "m": 0})
            e["T"] += T
            e["k"] += 1
            e["r"].append(abs(r[0]))
            e["m"] += nmk
        info = {"skipped": self.skipped,
                "informative": {f"{self.fams[fi].id}:{role}": n for (fi, role), n in self.n_informative.items()}}
        for key, per_chrom in raw_tracks.items():
            self.tracks[key] = {c: (np.array([w for w, _ in sorted(v, key=lambda t: t[0])]),
                                    np.vstack([u for _, u in sorted(v, key=lambda t: t[0])]))
                                for c, v in per_chrom.items()}
        for f in self.fams:
            self.family_meta[f.id] = dict(names=[self.unit.names[i] for i in f.offspring], sex=f.sex.copy())
        if not acc:
            return pd.DataFrame(), {}, info
        keys = sorted(acc, key=lambda k: (k[1], k[2], k[0]))
        T = np.vstack([acc[k]["T"] for k in keys])
        p, padj = maxz_adjust(T)
        long = pd.DataFrame({
            "role": [k[0] for k in keys], "chrom": [k[1] for k in keys], "widx": [k[2] for k in keys],
            "T": T[:, 0], "p": p, "padj": padj,
            "k": [acc[k]["k"] for k in keys], "absr": [float(np.mean(acc[k]["r"])) for k in keys],
            "n_markers": [acc[k]["m"] for k in keys]})
        wide = None
        for role in ("sire", "dam"):
            sub = long[long.role == role].drop(columns="role").rename(columns={
                c: f"ped_{role}_{c}" for c in ("T", "p", "padj", "k", "absr", "n_markers")})
            wide = sub if wide is None else wide.merge(sub, on=["chrom", "widx"], how="outer")
        # per-family tables (each family is tested on its own, then combined above)
        fam_tables = {}
        for fi, rows in fam_rows.items():
            Tm = np.vstack([r[3] for r in rows])
            pu, pa = maxz_adjust(Tm)
            df = pd.DataFrame({"role": [r[0] for r in rows], "chrom": [r[1] for r in rows],
                               "widx": [r[2] for r in rows], "T": Tm[:, 0], "p": pu, "padj": pa})
            fam_tables[self.fams[fi].id] = df
        return wide, fam_tables, info
