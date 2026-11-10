"""Table-based checks that need no genotype data: alignment records, coverage, batch confounding."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps

from .findings import Findings, compress_ids, short_list
from .samples import Unit, design_checks
from .stats import perm_chi2_pvalue
from .windows import robust_z

QC_METRICS = ["UnmappedReads", "DiscordantReads", "InversionReads", "DuplicationReads", "ClippedReads",
              "ChimericReads"]


def table_checks(unit: Unit, tables, cfg, F: Findings) -> None:
    design_checks(unit, cfg, F)
    names = unit.names
    parents = {names[i] for i in unit.parent_idx}
    al = tables.alignment
    if al is None:
        F.info("no_alignment_table", "No alignment table supplied: coverage/QC pre-checks skipped "
               "(coverage is still estimated from VCF depth).")
    else:
        a = al[al.SampleID.isin(names)].drop_duplicates("SampleID").set_index("SampleID")
        miss = [s for s in names if s not in a.index]
        if miss:
            F.info("no_alignment_record", f"{len(miss)}/{len(names)} samples have no alignment record "
                   "(coverage checked from VCF depth only): " + short_list(miss), miss)
        if "GenomeVersion" in a:
            bad = a.index[a["GenomeVersion"].astype(str) != tables.genome_version]
            if len(bad):
                found = a.loc[bad, "GenomeVersion"].astype(str).value_counts().to_dict()
                F.block("genome_version_mismatch", f"{len(bad)} of {len(a)} samples with alignment records were aligned "
                        f"to a different genome than '{tables.genome_version}' (found: {found}): "
                        + compress_ids(bad) + ".", list(bad),
                        why="Coordinates, coverage and mapping quality are only comparable when every sample was "
                            "aligned to the same reference genome.",
                        fix=f"Realign those samples to {tables.genome_version}, or run with --genome-version <the genome "
                            "the VCF was really called against> (per file: manifest.tsv).",
                        rows=[dict(sample=s, column="GenomeVersion", current=str(a.at[s, "GenomeVersion"]),
                                   needed=tables.genome_version) for s in bad])
        if "Coverage" in a:
            cov = pd.to_numeric(a["Coverage"], errors="coerce").dropna()
            fam_off = {names[i] for f in unit.families for i in f.offspring}
            # family offspring feed the transmission test, which tolerates low depth; everyone else needs strict depth
            low_ped = [s for s in cov.index if s in fam_off and cov[s] < cfg.ped_cov_warn]
            contrast_only = [s for s in cov.index if s in fam_off and cfg.ped_cov_warn <= cov[s] < cfg.cov_warn]
            low = [s for s in cov.index if s not in fam_off and s not in parents and cov[s] < cfg.cov_warn]
            if low_ped:
                F.warn("low_coverage_pedigree", f"{len(low_ped)} family offspring have alignment coverage < "
                       f"{cfg.ped_cov_warn:g}x: " + short_list((f"{s} ({cov[s]:.1f}x)" for s in low_ped), 14), low_ped,
                       why="The transmission test tolerates low depth, but below about "
                           f"{cfg.ped_cov_warn:g}x most offspring genotypes are missing and the inherited haplotype "
                           "cannot be read reliably.",
                       fix="Remove these offspring from the VCF or re-sequence them.",
                       rows=[dict(sample=s, column="Coverage", current=f"{cov[s]:.1f}x", needed=f">= {cfg.ped_cov_warn:g}x")
                             for s in low_ped])
            if contrast_only:
                F.info("low_coverage_contrast_only", f"{len(contrast_only)} family offspring have alignment coverage "
                       f"< {cfg.cov_warn:g}x ({short_list((f'{s} ({cov[s]:.1f}x)' for s in contrast_only), 14)}). They are "
                       "used by the transmission test (depth >= "
                       f"{cfg.ped_min_dp}) but are expected to be left out of the sex contrasts (association, Fst, "
                       "heterozygosity), which need confident genotypes.", contrast_only)
            if low:
                kept = [sps.poisson.sf(cfg.min_dp - 1, cov[s]) for s in low]
                F.warn("low_coverage_table", f"{len(low)} samples have alignment coverage < {cfg.cov_warn:g}x "
                       f"(< {cfg.cov_bad:g}x is unusable for genotyping): " + short_list(
                           (f"{s} ({cov[s]:.1f}x)" for s in low), 14), low,
                       why=f"Genotypes with depth < {cfg.min_dp} are discarded, so at these depths only about "
                           f"{min(kept):.0%}-{max(kept):.0%} of a sample's genotypes survive (Poisson estimate); samples "
                           f"keeping fewer than {1 - cfg.max_sample_missing_contrast:.0%} are left out of the sex contrasts.",
                       fix="Remove the lowest-depth samples from the VCF, re-sequence them, or lower min_dp in the config "
                           "(more genotypes kept, more genotyping error).",
                       rows=[dict(sample=s, column="Coverage", current=f"{cov[s]:.1f}x", needed=f">= {cfg.cov_warn:g}x")
                             for s in low])
            lowp = [s for s in cov.index if s in parents and cov[s] < cfg.parent_cov_warn]
            if lowp:
                F.warn("low_parent_coverage", f"Parents below {cfg.parent_cov_warn:g}x (miscalled parent "
                       "genotypes add noise to the transmission test): " + short_list(
                           f"{s} ({cov[s]:.1f}x)" for s in lowp), lowp)
            if len(cov) >= 3 and cov.min() > 0 and cov.max() / cov.min() > cfg.cov_spread_fold:
                (F.info if unit.families else F.warn)("coverage_spread", f"Coverage ranges {cov.min():.1f}x-{cov.max():.1f}x "
                       f"(>{cfg.cov_spread_fold:g}-fold); heterozygosity-based signals depend on depth.")
        for col in ("minimap2_version", "gatk_version", "pysam_version"):
            if col in a and a[col].dropna().astype(str).nunique() > 1:
                F.info("tool_versions", f"Samples were processed with different {col} values: "
                       + ", ".join(sorted(a[col].dropna().astype(str).unique())))
        if len(a) >= 8:
            for m in QC_METRICS:
                if m not in a:
                    continue
                x = pd.to_numeric(a[m], errors="coerce")
                z = pd.Series(robust_z(x.to_numpy()), index=x.index)
                med = x.median()
                out = x.index[(z >= cfg.qc_outlier_z) & (x > 1.5 * med) & (x - med > 0.01)]
                if len(out):
                    F.warn("alignment_qc_outlier", f"{m} is an outlier for: " + short_list(
                        f"{s} ({x[s]:.3f} vs median {med:.3f})" for s in out), list(out))
    eco_ref = cfg.genome_ecogroup.get(tables.genome_version)
    if eco_ref and unit.ecogroup and unit.ecogroup != eco_ref:
        F.warn("reference_divergence", f"Reads were mapped to a {eco_ref} reference ({tables.genome_version}) but the "
               f"samples are {unit.ecogroup}; divergent reads map poorly, which biases heterozygosity and coverage.")
    _confounding(unit, tables, cfg, F)


def _confounding(unit: Unit, tables, cfg, F: Findings) -> None:
    idx = unit.contrast_idx
    names = [unit.names[i] for i in idx]
    sex = unit.sex[idx]
    variables = {}
    if tables.reads is not None:
        rd = tables.reads[tables.reads.SampleID.isin(names)]
        for col in ("ProjectID", "Instrument"):
            if col in rd:
                variables[col] = rd.groupby("SampleID")[col].agg(lambda x: "+".join(sorted(set(map(str, x.dropna())))))
    for col in ("LabReared", "BroodID"):
        if col in unit.rows:
            variables[col] = unit.rows[col]
    testable = []
    for var, ser in variables.items():
        v = ser.reindex(names)
        ok = v.notna() & (v.astype(str).str.strip() != "")
        v2 = v[ok].astype(str)
        counts = v2.value_counts()
        keep = v2.isin(counts[counts >= 2].index)
        if ok.sum() >= 8 and keep.sum() >= 8 and v2[keep].nunique() >= 2:
            testable.append(var)
    thr = cfg.confound_p / max(1, len(testable))          # Bonferroni across the variables examined
    for var in testable:
        v = variables[var].reindex(names)
        ok = v.notna() & (v.astype(str).str.strip() != "")
        if ok.sum() < 8:
            continue
        v = v[ok].astype(str)
        s = sex[ok.to_numpy()]
        counts = v.value_counts()
        keep = v.isin(counts[counts >= 2].index)
        if keep.sum() < 8 or v[keep].nunique() < 2:
            continue
        p = perm_chi2_pvalue(s[keep.to_numpy()], v[keep].to_numpy(), seed=cfg.seed)
        if p < thr:
            tab = pd.crosstab(v[keep], s[keep.to_numpy()])
            tab.columns = ["F" if c == 0 else "M" for c in tab.columns]
            desc = "; ".join(f"{k}: {int(r.get('M', 0))}M/{int(r.get('F', 0))}F" for k, r in tab.iterrows())
            F.warn("sex_confounded", f"Sex is associated with {var} (permutation p={p:.4f}): {desc[:300]}. Sex "
                   "differences may reflect this variable, not genetics; confidence capped at tier B.", cap="B")
