"""Data-dependent QC collected during the single VCF pass, and the warnings derived from it."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps

from .findings import Findings, short_list


class QCAccumulator:
    name = "qc"

    def __init__(self, unit, cfg):
        self.unit, self.cfg = unit, cfg
        n = len(unit.names)
        self.dp_sum = np.zeros(n)
        self.dp_cnt = np.zeros(n)
        self.called = np.zeros(n)
        self.het = np.zeros(n)
        self.n = 0
        self.have_dp = False
        self.fams = unit.families
        self.inf = [np.zeros(len(f.offspring)) for f in self.fams]
        self.err = [np.zeros(len(f.offspring)) for f in self.fams]

    def update(self, chunk):
        if chunk.dp is not None:
            self.have_dp = True
            v = ~np.isnan(chunk.dp)
            self.dp_sum += np.where(v, chunk.dp, 0).sum(0)
            self.dp_cnt += v.sum(0)
        G = chunk.dos_role if chunk.dos_role is not None else chunk.dos
        c = G >= 0
        self.called += c.sum(0)
        self.het += (G == 1).sum(0)
        self.n += len(chunk.pos)
        for k, f in enumerate(self.fams):
            s, d = G[:, f.sire], G[:, f.dam]
            O = G[:, f.offspring]
            # Only sites where BOTH parents are homozygous for the SAME allele: the offspring must carry it too. Sites with
            # opposite homozygous parents are skipped on purpose - there a low-depth heterozygous offspring is miscalled
            # homozygous, which would look like a Mendelian error although nothing is wrong with the parentage.
            mask = ((s == 0) | (s == 2)) & (s == d)
            if not mask.any():
                continue
            exp = s[mask][:, None]
            Om = O[mask]
            ok = Om >= 0
            self.inf[k] += ok.sum(0)
            self.err[k] += (ok & (Om != exp)).sum(0)

    def sample_table(self):
        u = self.unit
        with np.errstate(invalid="ignore", divide="ignore"):
            depth = np.where(self.dp_cnt > 0, self.dp_sum / self.dp_cnt, np.nan)
            par, exc = set(u.parent_idx), set(u.excluded_idx)
            ped = {i for f in u.families for i in f.offspring}
            return pd.DataFrame({
                "sample": u.names, "sex": ["M" if x == 1 else "F" for x in u.sex],
                "role": ["parent" if i in par else ("offspring, transmission only" if i in exc and i in ped else
                         "excluded from contrasts" if i in exc else "analysed") for i in range(len(u.names))],
                "mean_depth": depth, "missing_rate": 1 - self.called / max(self.n, 1),
                "het_rate": self.het / np.maximum(self.called, 1)})


def qc_findings(qc: QCAccumulator, tables, cfg, F: Findings) -> pd.DataFrame:
    u = qc.unit
    st = qc.sample_table().set_index("sample")
    an = [u.names[i] for i in u.contrast_idx]
    par = [u.names[i] for i in u.parent_idx]
    ped_off = [u.names[i] for f in u.families for i in f.offspring]
    others = [s for s in u.names if s not in set(par) | set(ped_off)]
    if not qc.have_dp:
        F.info("no_dp", "The VCF has no DP field: coverage, sex-depth and heterozygosity-depth checks skipped.")
    else:
        d = st["mean_depth"]
        # family offspring are judged by the (tolerant) transmission-test thresholds, everyone else by the strict ones
        bad = [s for s in ped_off if d[s] < cfg.ped_cov_bad] + [s for s in others if d[s] < cfg.cov_bad]
        low = [s for s in ped_off if cfg.ped_cov_bad <= d[s] < cfg.ped_cov_warn] + [
            s for s in others if cfg.cov_bad <= d[s] < cfg.cov_warn]
        if bad:
            F.warn("unusable_coverage", f"{len(bad)} samples have mean VCF depth below the usable minimum "
                   f"({cfg.ped_cov_bad:g}x for family offspring, {cfg.cov_bad:g}x for other animals): "
                   + short_list((f"{s} ({d[s]:.1f}x)" for s in bad), 14), bad, cap="B",
                   fix="Remove these samples from the VCF or re-sequence them.",
                   rows=[dict(sample=s, column="mean VCF depth", current=f"{d[s]:.1f}x", needed="usable depth") for s in bad])
        if low:
            F.warn("low_coverage", f"{len(low)} samples have low mean VCF depth: "
                   + short_list((f"{s} ({d[s]:.1f}x)" for s in low), 14), low)
        lowp = [s for s in par if cfg.cov_bad <= d[s] < cfg.parent_cov_warn]
        if lowp:
            F.warn("low_parent_depth", f"Parents below {cfg.parent_cov_warn:g}x in the VCF: "
                   + short_list(f"{s} ({d[s]:.1f}x)" for s in lowp), lowp,
                   why="Parent genotypes decide which markers are informative, so parents should be well covered.")
        dd = d[an + par].dropna()
        if len(dd) >= 3 and dd.min() > 0 and dd.max() / dd.min() > cfg.cov_spread_fold and not u.families:
            F.warn("depth_spread", f"VCF depth ranges {dd.min():.1f}x-{dd.max():.1f}x across samples "
                   f"(>{cfg.cov_spread_fold:g}-fold).")
        # sex vs depth
        dm = d[[s for s in an if st.at[s, "sex"] == "M"]].dropna()
        df_ = d[[s for s in an if st.at[s, "sex"] == "F"]].dropna()
        if len(dm) >= 3 and len(df_) >= 3 and df_.median() > 0:
            ratio = dm.median() / df_.median()
            p = sps.mannwhitneyu(dm, df_).pvalue
            if ratio < cfg.sex_depth_ratio_lo or ratio > cfg.sex_depth_ratio_hi or (
                    p < cfg.sex_depth_p and (ratio < 0.8 or ratio > 1.25)):
                F.warn("sex_depth_difference", f"Median depth differs between sexes (M/F = {ratio:.2f}, "
                       f"Mann-Whitney p={p:.3g}). Heterozygosity and coverage contrasts are unreliable; "
                       "confidence capped at tier B.", cap="B")
        both = st.loc[an].dropna(subset=["mean_depth"])
        if len(both) >= 8:
            rho, p = sps.spearmanr(both["mean_depth"], both["het_rate"])
            if np.isfinite(rho) and abs(rho) >= cfg.het_depth_rho and p < 0.01:
                F.warn("het_tracks_depth", f"Per-sample heterozygosity tracks depth (Spearman rho={rho:.2f}); "
                       "low-depth heterozygotes are undercalled, which distorts male-female heterozygosity contrasts.")
    miss = [s for s in dict.fromkeys(an + par + ped_off) if st.at[s, "missing_rate"] > cfg.sample_missing_warn]
    if miss:
        F.warn("high_missingness", f"{len(miss)} samples are missing > {cfg.sample_missing_warn:.0%} of genotypes "
               "after QC: " + short_list(f"{s} ({st.at[s, 'missing_rate']:.0%})" for s in miss), miss)
    # table coverage vs VCF depth
    al = tables.alignment
    if al is not None and qc.have_dp and "Coverage" in al:
        a = al.drop_duplicates("SampleID").set_index("SampleID")
        common = [s for s in st.index if s in a.index and np.isfinite(st.at[s, "mean_depth"])
                  and pd.notna(a.at[s, "Coverage"]) and a.at[s, "Coverage"] > 0]
        off = [s for s in common if not 0.4 <= st.at[s, "mean_depth"] / a.at[s, "Coverage"] <= 2.5]
        if off:
            F.warn("coverage_mismatch", "VCF depth and alignment-table coverage disagree by > 2.5-fold (wrong "
                   "sample ID mapping or a different alignment?): " + short_list(
                       f"{s} (VCF {st.at[s, 'mean_depth']:.1f}x vs table {a.at[s, 'Coverage']:.1f}x)" for s in off), off)
    # Mendelian consistency of declared trios
    for k, f in enumerate(qc.fams):
        inf, err = qc.inf[k], qc.err[k]
        ok = inf >= cfg.trio_min_sites
        if not ok.any():
            F.info("trio_uninformative", f"{f.id}: too few sites with both parents homozygous for the same allele "
                   f"(< {cfg.trio_min_sites}) to check parentage.")
            continue
        rate = np.where(ok, err / np.maximum(inf, 1), np.nan)
        badk = np.flatnonzero(rate > cfg.trio_error_warn)
        if len(badk):
            F.warn("mendelian_errors", f"{f.id}: {len(badk)} offspring have Mendelian-error rates > "
                   f"{cfg.trio_error_warn:.0%} against the listed parents (wrong parents or sample swap): "
                   + short_list(f"{u.names[f.offspring[i]]} ({rate[i]:.1%})" for i in badk),
                   [u.names[f.offspring[i]] for i in badk], cap="B")
    return st.reset_index()
