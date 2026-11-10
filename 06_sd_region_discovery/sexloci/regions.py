"""Combine per-window evidence into candidate sex-linked regions and grade them.

Evidence classes
  * transmission  - parent-to-offspring inheritance of a parental haplotype predicts offspring sex
  * contrast      - SNP/window differences between the sexes (association, Fst, delta-heterozygosity, diagnostic SNPs)
  * depth         - male:female read-depth ratio (supporting only)
Association and Fst/dHet are the same information, so they never count as separate confirmation.
Confidence comes from *independent units* (families that share no parent, and the pool of non-family
animals) each reaching family-wise significance in the region:
  tier A  >= 2 independent units significant (replication)
  tier B  the combined test is significant (permutation FWER) but replication is not shown
  tier C  population-genetic outlier only
A tier can be capped downward by data-quality findings (e.g. sex confounded with batch).
"""
from __future__ import annotations

import re
from typing import Optional

import numpy as np
import pandas as pd

from .windows import robust_z


def merge_windows(*tables, window_size: int) -> pd.DataFrame:
    tabs = [t for t in tables if t is not None and len(t)]
    if not tabs:
        return pd.DataFrame(columns=["chrom", "widx", "start", "end"])
    df = tabs[0]
    for t in tabs[1:]:
        df = df.merge(t, on=["chrom", "widx"], how="outer")
    df["start"] = df["widx"] * window_size + 1
    df["end"] = (df["widx"] + 1) * window_size
    return df.sort_values(["chrom", "widx"]).reset_index(drop=True)


def _peak(w, padj_col, stat_col, alpha, frac):
    if padj_col not in w:
        return pd.Series(False, index=w.index), pd.Series(False, index=w.index)
    sig = (w[padj_col] <= alpha).fillna(False)
    mx = w[stat_col].where(sig).groupby(w["chrom"]).transform("max")
    return sig, sig & (w[stat_col] >= frac * mx)


def flag_evidence(win: pd.DataFrame, cfg) -> pd.DataFrame:
    """Per-window flags. `sig_*` = significant anywhere on the linked chromosome; `ev_*` = on the plateau
    (>= peak_frac of the best statistic on that chromosome), the best estimate of the locus position."""
    w = win.copy()
    g = lambda c: w[c] if c in w else pd.Series(np.nan, index=w.index)
    w["sig_assoc"], w["ev_assoc"] = _peak(w, "assoc_padj", "assoc_max_chi2", cfg.alpha, cfg.peak_frac)
    ped_sig = pd.Series(False, index=w.index)
    ped_ev = pd.Series(False, index=w.index)
    for role in ("sire", "dam"):
        s_, e_ = _peak(w, f"ped_{role}_padj", f"ped_{role}_T", cfg.alpha, cfg.peak_frac)
        w[f"sig_ped_{role}"], w[f"ev_ped_{role}"] = s_, e_
        ped_sig |= s_
        ped_ev |= e_
    w["sig_ped"], w["ev_ped"] = ped_sig, ped_ev
    enough = (w["pg_n_markers"] >= cfg.min_evidence_markers) if "pg_n_markers" in w else pd.Series(True, index=w.index)
    if "fst" in w:
        w["fst_z"] = robust_z(w["fst"])
        w["dhet_z"] = robust_z(w["dhet"])
        w["ev_fst"] = ((w["fst_z"] >= cfg.z_thr) & (w["fst"] >= cfg.min_fst) & enough).fillna(False)
        w["ev_dhet"] = ((w["dhet_z"].abs() >= cfg.z_thr) & (w["dhet"].abs() >= cfg.min_dhet) & enough).fillna(False)
        diag = w[["n_diag_xy", "n_diag_zw"]].max(axis=1)
        thr = max(cfg.min_diag, np.nanquantile(diag.dropna(), 0.995) if diag.notna().any() else 0)
        w["ev_diag"] = ((diag >= thr) & enough).fillna(False)
    else:
        for c in ("ev_fst", "ev_dhet", "ev_diag"):
            w[c] = False
    if "cov_log2_mf" in w:
        w["cov_z"] = robust_z(w["cov_log2_mf"])
        w["ev_cov"] = ((w["cov_z"].abs() >= cfg.z_thr) & (w["cov_log2_mf"].abs() >= cfg.min_cov_log2) & enough).fillna(False)
    else:
        w["ev_cov"] = False
    w["n_popgen_ev"] = w[["ev_fst", "ev_dhet", "ev_diag", "ev_cov"]].sum(axis=1)
    w["candidate"] = w["ev_assoc"] | w["ev_ped"] | (w["n_popgen_ev"] >= 2)
    return w


def pedigree_chromosome_summary(w: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for role in ("sire", "dam"):
        if f"ped_{role}_T" not in w:
            continue
        for chrom, sub in w.groupby("chrom", sort=False):
            sub = sub.dropna(subset=[f"ped_{role}_T"])
            if sub.empty:
                continue
            best = sub.loc[sub[f"ped_{role}_T"].idxmax()]
            rows.append(dict(parent=role, chrom=chrom, n_windows=len(sub),
                             n_sig_windows=int(sub[f"sig_ped_{role}"].sum()), max_T=float(best[f"ped_{role}_T"]),
                             best_window_start=int(best["start"]), best_padj=float(best[f"ped_{role}_padj"]),
                             mean_absr=float(sub[f"ped_{role}_absr"].mean())))
    return pd.DataFrame(rows)


def _unit_key(unit, label):
    return "pool" if not label.startswith("fam:") else f"family-group {unit.component.get(label, label)}"


def _support(chrom, wlo, whi, unit, strata_tables, fam_tables, alpha):
    """Which independent units reach FWER significance inside [wlo, whi] on `chrom`."""
    sig, labels = {}, set()
    def rng(t):
        return t[(t.chrom == chrom) & (t.widx >= wlo) & (t.widx <= whi)]
    for label, t in list((strata_tables or {}).items()) + list((fam_tables or {}).items()):
        k = _unit_key(unit, label)
        hit = bool((rng(t).padj <= alpha).any())
        sig[k] = sig.get(k, False) or hit
        if hit:
            labels.add(label)
    return sig, labels


def decide_model(sire_sig, dam_sig, dhet, dhet_evidence, n_xy, n_zw, diag_evidence, cfg):
    """XY / ZW. The pedigree decides when it can: whichever parent's inherited haplotype predicts sex is the heterogametic one.
    The sign of dHet and the diagnostic-SNP counts are only used when there is no pedigree evidence, and then only when they were
    outliers in the region. With a single sire the sign of dHet is NOT diagnostic: alleles private to the sire's X make daughters
    heterozygous (negative dHet) in a perfectly ordinary XY region."""
    ped = []
    if sire_sig and not dam_sig:
        ped.append("pedigree:XY")
    elif dam_sig and not sire_sig:
        ped.append("pedigree:ZW")
    elif sire_sig and dam_sig:
        ped.append("pedigree:both")
    pop = []
    if dhet_evidence and np.isfinite(dhet) and abs(dhet) >= cfg.min_dhet:
        pop.append("dhet:XY" if dhet > 0 else "dhet:ZW")
    if diag_evidence and n_xy + n_zw >= cfg.min_diag and n_xy != n_zw:
        pop.append("diag:XY" if n_xy > n_zw else "diag:ZW")
    votes = ped + pop
    use = [v for v in ped if not v.endswith("both")] or pop
    systems = {v.split(":")[1] for v in use}
    model = "undetermined" if not systems else (next(iter(systems)) if len(systems) == 1 else "conflicting")
    return model, votes


def call_regions(w, hits, cfg, unit, strata_tables, fam_tables, tier_cap=None, genes=None,
                 alpha_unit=None) -> pd.DataFrame:
    alpha_unit = cfg.alpha if alpha_unit is None else alpha_unit
    cols = ["chrom", "start", "end", "size_bp", "n_windows", "tier", "tier_raw", "model", "model_votes",
            "evidence", "units_supporting", "units_share_parent", "units_total", "supporting_units", "assoc_max_chi2",
            "assoc_min_padj", "ped_sire_min_padj", "ped_dam_min_padj", "ped_sire_absr", "ped_dam_absr",
            "mean_dhet", "mean_fst", "n_diag_xy", "n_diag_zw", "mean_cov_log2_mf", "core_start", "core_end", "linked_start", "linked_end",
            "n_sig_markers", "top_marker_pos", "sig_marker_span", "flag", "genes"]
    cand = w[w["candidate"]]
    if cand.empty:
        return pd.DataFrame(columns=cols)
    rows = []
    for chrom, sub in cand.groupby("chrom", sort=False):
        sub = sub.sort_values("widx")
        grp = (sub["widx"].diff().fillna(0) > cfg.merge_gap + 1).cumsum()
        for _, blk in sub.groupby(grp):
            lo, hi = int(blk["start"].min()), int(blk["end"].max())
            span = w[(w.chrom == chrom) & (w.start >= lo) & (w.end <= hi)]
            rows.append(_summarise(chrom, lo, hi, blk, span, hits, cfg, unit, strata_tables, fam_tables,
                                   tier_cap, genes, w, alpha_unit))
    out = pd.DataFrame(rows, columns=cols)
    order = {"A": 0, "B": 1, "C": 2}
    out = out.sort_values(by=["tier", "chrom", "start"], key=lambda s: s.map(order) if s.name == "tier" else s)
    return out.reset_index(drop=True)


def _nmin(s):
    return float(s.min()) if s.notna().any() else np.nan


def _nmax(s):
    return float(s.max()) if s.notna().any() else np.nan


def _nmean(s):
    return float(s.mean()) if s.notna().any() else np.nan


def _summarise(chrom, lo, hi, blk, span, hits, cfg, unit, strata_tables, fam_tables, tier_cap, genes, w,
               alpha_unit):
    get = lambda c: span[c] if c in span else pd.Series(dtype=float)
    ev = [c[3:] for c in ("ev_assoc", "ev_ped", "ev_fst", "ev_dhet", "ev_diag", "ev_cov") if c in blk and blk[c].any()]
    combined_sig = bool(blk["ev_assoc"].any() or blk["ev_ped"].any())
    wlo, whi = (lo - 1) // cfg.window_size, (hi - 1) // cfg.window_size
    sup, sup_labels = _support(chrom, wlo, whi, unit, strata_tables, fam_tables, alpha_unit)
    shared_by = [p_ for p_, labs in getattr(unit, "shared_parents", {}).items() if len(set(labs) & sup_labels) >= 2]
    n_sup = sum(sup.values())
    if combined_sig and n_sup >= cfg.min_support_units:
        raw = "A"
    elif combined_sig:
        raw = "B"
    else:
        raw = "C"
    tier = max(raw, tier_cap) if tier_cap else raw

    sire_p, dam_p = _nmin(get("ped_sire_padj")), _nmin(get("ped_dam_padj"))
    dhet = _nmean(get("dhet"))
    nxy = int(get("n_diag_xy").sum()) if "n_diag_xy" in span else 0
    nzw = int(get("n_diag_zw").sum()) if "n_diag_zw" in span else 0
    model, votes = decide_model(sire_p <= cfg.alpha, dam_p <= cfg.alpha, dhet,
                                bool(span["ev_dhet"].any()) if "ev_dhet" in span else False, nxy, nzw,
                                bool(span["ev_diag"].any()) if "ev_diag" in span else False, cfg)
    # differentiated core: windows where male-female heterozygosity is an outlier (the heterogametic signature). The
    # transmission plateau is wider because linkage fades slowly with few recombinants.
    core = span[span["ev_dhet"]] if "ev_dhet" in span else span.iloc[0:0]
    core_start = int(core["start"].min()) if len(core) else np.nan
    core_end = int(core["end"].max()) if len(core) else np.nan

    h = hits[(hits.chrom == chrom) & (hits.pos >= lo) & (hits.pos <= hi)] if len(hits) else hits
    sigm = h[h["padj"] <= cfg.alpha] if len(h) else h
    flag = ""
    if hi - lo + 1 >= cfg.block_flag_mb * 1e6:
        flag = ("extended_block: large contiguous sex-associated region; consistent with a non-recombining "
                "segment (e.g. inversion) but verify with long reads / LD / PCA - not assumed")
    if len(sigm) and (sigm.pos.max() - sigm.pos.min()) < cfg.narrow_spike_bp and len(sigm) <= cfg.narrow_spike_max_markers \
            and (hi - lo + 1) <= 2 * cfg.window_size:
        neigh = w[(w.chrom == chrom) & (((w.start == lo - cfg.window_size)) | ((w.end == hi + cfg.window_size)))]
        if not (neigh["sig_assoc"] | neigh["sig_ped"]).any():
            flag = (flag + " | " if flag else "") + (
                f"narrow_spike: only {len(sigm)} significant SNPs within {(sigm.pos.max() - sigm.pos.min()) / 1e3:.1f} kb and "
                "no shoulders in the neighbouring windows. A linked locus in a family raises neighbouring windows too, so "
                "this looks like Y/W-linked sequence mapping here, a mapping artifact or a small structural variant - "
                "inspect the reads in each sex before treating it as a locus")
    gl = ""
    if genes is not None and len(genes):
        gg = genes[(genes.chrom == chrom) & (genes.end >= lo) & (genes.start <= hi)]
        gl = ";".join(gg["name"].head(60)) + (";..." if len(gg) > 60 else "")
    # linked interval: extend through adjacent significant windows
    cw = w[w.chrom == chrom].sort_values("widx")
    sigw = (cw["sig_assoc"] | cw["sig_ped"]).to_numpy()
    wi = cw["widx"].to_numpy()
    a_ = int(np.flatnonzero((cw["start"] >= lo).to_numpy())[0])
    b_ = int(np.flatnonzero((cw["end"] <= hi).to_numpy())[-1])
    while a_ > 0 and sigw[a_ - 1] and wi[a_] - wi[a_ - 1] == 1:
        a_ -= 1
    while b_ < len(wi) - 1 and sigw[b_ + 1] and wi[b_ + 1] - wi[b_] == 1:
        b_ += 1
    return dict(
        chrom=chrom, start=lo, end=hi, size_bp=hi - lo + 1, n_windows=int(len(blk)), tier=tier, tier_raw=raw,
        model=model, model_votes=",".join(votes), evidence=",".join(ev), units_supporting=int(n_sup), units_share_parent=",".join(shared_by),
        units_total=len(sup), supporting_units=";".join(k for k, v in sup.items() if v),
        assoc_max_chi2=_nmax(get("assoc_max_chi2")), assoc_min_padj=_nmin(get("assoc_padj")),
        ped_sire_min_padj=sire_p, ped_dam_min_padj=dam_p, ped_sire_absr=_nmean(get("ped_sire_absr")),
        ped_dam_absr=_nmean(get("ped_dam_absr")), mean_dhet=dhet, mean_fst=_nmean(get("fst")), n_diag_xy=nxy,
        n_diag_zw=nzw, mean_cov_log2_mf=_nmean(get("cov_log2_mf")), core_start=core_start, core_end=core_end,
        linked_start=int(cw["start"].iloc[a_]),
        linked_end=int(cw["end"].iloc[b_]), n_sig_markers=int(len(sigm)),
        top_marker_pos=int(h.iloc[0]["pos"]) if len(h) else np.nan,
        sig_marker_span=f"{int(sigm.pos.min())}-{int(sigm.pos.max())}" if len(sigm) else "", flag=flag, genes=gl)


def read_gff_genes(path: str) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            rows.append((f[0], int(f[3]), int(f[4]),
                         re.sub(r"^gene-", "", attr.get("Name") or attr.get("gene") or attr.get("ID", "."))))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "name"])


def read_gff_regions(path: str) -> dict:
    """seqid -> linkage-group label from the chromosome 'region' records of an NCBI GFF3 (chromosome= / Name=)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "region":
                continue
            attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if attr.get("genome", "chromosome") != "chromosome":
                continue
            lab = attr.get("chromosome") or attr.get("Name")
            if not lab or lab.lower() in ("mt", "unknown", "un"):
                continue
            out[f[0]] = f"LG{lab}" if lab.isdigit() else ("LG" + lab[2:] if lab.upper().startswith("LG") else lab)
    return out


def match_annotation(vcf_contigs, cfg, gff_labels):
    """Labels for plots and the VCF-contig -> annotation-contig aliases.

    labels: the config's table (e.g. UMD2a accessions -> LG) plus the annotation's own (GT3a accessions -> LG).
    aliases: a VCF contig whose linkage group is found in the annotation under a different name (YH was called against GT3 with
    UMD2a chromosome names) is mapped to that annotation chromosome, so genes can be listed."""
    labels = {**gff_labels, **(cfg.contig_labels or {})}
    by_lg = {}
    for seq, lg in gff_labels.items():
        by_lg.setdefault(lg, seq)
    aliases = {}
    if cfg.annotation_alias_by_linkage_group:
        for c in vcf_contigs:
            lg = labels.get(c)
            if lg and c not in gff_labels and lg in by_lg:
                aliases[c] = by_lg[lg]
    return labels, aliases


def alias_genes(genes: pd.DataFrame, aliases: dict) -> pd.DataFrame:
    """Add copies of the gene rows under the VCF's own chromosome names."""
    if genes is None or genes.empty or not aliases:
        return genes
    extra = [genes[genes.chrom == seq].assign(chrom=vcf) for vcf, seq in aliases.items()]
    return pd.concat([genes] + extra, ignore_index=True)
