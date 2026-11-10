"""Preflight (table checks, BLOCKs) -> one pass over the VCF -> QC warnings -> region calling -> outputs."""
from __future__ import annotations

import json
import os
import time

import numpy as np
import pandas as pd

import dataclasses

from . import __version__
from .association import AssociationAccumulator
from .config import Config, config_dict
from .findings import Findings, short_list
from .pedigree import PedigreeAccumulator
from .plots import chromosome_summary_plot, genome_plot, haplotype_genome_plot, region_plot
from .popgen import PopgenAccumulator
from .preflight import table_checks
from .qc import QCAccumulator, qc_findings
from .regions import (alias_genes, call_regions, flag_evidence, match_annotation, merge_windows,
                      pedigree_chromosome_summary, read_gff_genes, read_gff_regions)
from .report import findings_frame, fix_list_frame, write_preflight_report, write_report
from .samples import build_unit
from .tables import Tables
from .vcf import VariantFilter, iter_chunks, mask_by_sample, read_header_info, read_samples, scan_sites


class PreflightBlocked(Exception):
    def __init__(self, findings):
        super().__init__("; ".join(f.message for f in findings if f.severity == "BLOCK"))
        self.findings = findings


def header_filter_findings(vcf: str, cfg: Config, F: Findings) -> None:
    """Header-only heuristics: can bad variants be filtered out of this VCF at all?"""
    if not cfg.apply_variant_filter:
        F.info("variant_filter_off", "Variant filtering is switched off (apply_variant_filter: false).")
        return
    hdr = read_header_info(vcf)
    vf = VariantFilter(cfg, hdr["infos"])
    usable = [h for h in vf.applicable() if h[0] != "QUAL"]
    if not usable and not (hdr["filters"] - {"PASS"}):
        F.warn("no_variant_filtering_possible", "The VCF declares none of the INFO fields used for hard filters "
               f"({', '.join(sorted({h[0] for h in vf.hard if h[0] != 'QUAL'}))}) and no FILTER labels, so bad "
               "variants cannot be removed here. Annotate/filter the VCF (e.g. GATK VariantFiltration) first.")
    elif vf.not_applicable():
        F.info("filters_not_applicable", "Hard filters skipped because the INFO field is absent: "
               + ", ".join(vf.label(h) for h in vf.not_applicable()))


def preflight(vcf: str, tables: Tables, cfg: Config, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    F = Findings()
    vcf_samples = read_samples(vcf)
    unit = build_unit(vcf_samples, tables, cfg, F)
    if unit is not None:
        table_checks(unit, tables, cfg, F)
    header_filter_findings(vcf, cfg, F)
    F.acknowledge(cfg.acknowledged_codes)
    findings_frame(F).to_csv(os.path.join(outdir, "findings.tsv"), sep="\t", index=False)
    fix = fix_list_frame(F)
    fix_path = os.path.join(outdir, "fix_list.tsv")
    if len(fix):
        fix.to_csv(fix_path, sep="\t", index=False)
    elif os.path.exists(fix_path):
        os.remove(fix_path)                      # a rerun that no longer needs fixes must not leave a stale list
    write_preflight_report(os.path.join(outdir, "preflight_report.md"), F, unit, tables)
    return unit, F, vcf_samples


def _prescan_missing(vcf, samples, cfg, vf) -> np.ndarray:
    """Strict genotype-missing rate of every sample over the first `prescan_variants` sites that pass the variant filter."""
    n = len(samples)
    called, seen = np.zeros(n), 0
    gen = iter_chunks(vcf, list(range(n)), samples, cfg, vf.clone())
    try:
        for chunk in gen:
            called += (chunk.dos >= 0).sum(axis=0)
            seen += len(chunk.pos)
            if seen >= cfg.prescan_variants:
                break
    finally:
        gen.close()
    return 1 - called / max(seen, 1)


def _rewrite_preflight(outdir, F, unit, tables):
    findings_frame(F).to_csv(os.path.join(outdir, "findings.tsv"), sep="\t", index=False)
    fix = fix_list_frame(F)
    if len(fix):
        fix.to_csv(os.path.join(outdir, "fix_list.tsv"), sep="\t", index=False)
    write_preflight_report(os.path.join(outdir, "preflight_report.md"), F, unit, tables)


def run_pipeline(vcf: str, tables: Tables, outdir: str, cfg: Config, quiet: bool = False) -> dict:
    t0 = time.time()
    log = (lambda *a: None) if quiet else (lambda *a: print("[sexloci]", *a, flush=True))
    unit, F, vcf_samples = preflight(vcf, tables, cfg, outdir)
    if F.blocked or unit is None:
        raise PreflightBlocked(F)
    # ---- variant filter (optionally a fast sites-only scan for the INFO/DP quantile bounds)
    hdr = read_header_info(vcf)
    vf = VariantFilter(cfg, hdr["infos"])
    if cfg.apply_variant_filter and cfg.variant_depth_quantiles:
        if "DP" in hdr["infos"]:
            scan_vf = VariantFilter(dataclasses.replace(cfg, apply_variant_filter=False), hdr["infos"])
            log("scanning INFO/DP for depth bounds ...")
            scan_sites(vcf, cfg, scan_vf)
            vf.set_depth_bounds(scan_vf.dp_values, cfg.variant_depth_quantiles)
        else:
            F.info("no_info_dp", "INFO/DP is absent: the depth-quantile site filter was skipped.")

    # ---- genotype thresholds per role: the transmission test is tolerant of low depth, everything else is strict
    n_s = len(unit.names)
    min_dp_v, min_gq_v = np.full(n_s, float(cfg.min_dp)), np.full(n_s, float(cfg.min_gq))
    for f in unit.families:
        min_dp_v[f.offspring], min_gq_v[f.offspring] = cfg.ped_min_dp, cfg.ped_min_gq
    min_dp_v[unit.parent_idx], min_gq_v[unit.parent_idx] = cfg.ped_parent_min_dp, cfg.ped_parent_min_gq

    # ---- which samples can take part in the sex contrasts? Measured on the first sites of the VCF.
    miss = _prescan_missing(vcf, vcf_samples, cfg, vf)
    drop = [int(i) for i in unit.contrast_idx if miss[i] > cfg.max_sample_missing_contrast]
    if drop:
        ped_off = {int(i) for f in unit.families for i in f.offspring}
        pool_drop = [i for i in drop if i not in ped_off]
        names_d = [unit.names[i] for i in drop]
        (F.warn if pool_drop else F.info)(
            "samples_excluded_from_contrasts",
            f"{len(drop)} samples keep fewer than {1 - cfg.max_sample_missing_contrast:.0%} of their genotypes at depth >= "
            f"{cfg.min_dp}" + (f" / GQ >= {cfg.min_gq}" if cfg.min_gq else "") + " and are left out of the sex contrasts "
            f"(association, Fst, heterozygosity): " + short_list((f"{unit.names[i]} ({miss[i]:.0%} missing)" for i in drop), 14)
            + (f". {len(drop) - len(pool_drop)} of them are family offspring and still count for the transmission test, "
               f"which only needs depth >= {cfg.ped_min_dp}." if len(drop) > len(pool_drop) else "."), names_d,
            why="The sex contrasts compare called genotypes between the sexes; a sample that is mostly missing adds "
                "noise and, through the call-rate filter, removes sites for everyone.",
            fix="Nothing is needed for family offspring. For other animals, deeper sequencing or a lower min_dp in the "
                "config would bring them back.",
            rows=[dict(sample=unit.names[i], column="genotype retention", current=f"{1 - miss[i]:.0%}",
                       needed=f">= {1 - cfg.max_sample_missing_contrast:.0%}") for i in pool_drop])
        unit.excluded_idx = np.array(drop, dtype=int)
        unit.contrast_idx = np.array([i for i in unit.contrast_idx if i not in set(drop)], dtype=int)
    sx = unit.sex[unit.contrast_idx]
    n_m, n_f = int((sx == 1).sum()), int((sx == 0).sum())
    contrasts_ok = min(n_m, n_f) >= cfg.min_per_sex_block
    if drop and contrasts_ok and min(n_m, n_f) < cfg.min_per_sex_warn and not any(f.code == "low_per_sex" for f in F):
        F.warn("low_per_sex", f"After leaving out low-retention samples, {n_m} males / {n_f} females remain for the sex "
               f"contrasts: power is low and permutation p-values are coarse (>= {cfg.min_per_sex_good} per sex recommended).")
    if not contrasts_ok:
        if not unit.families:
            F.block("too_few_usable_samples", f"Only {n_m} males / {n_f} females keep enough genotypes for the sex contrasts "
                    f"(at least {cfg.min_per_sex_block} of each are needed) and there is no pedigree to fall back on.",
                    why="With too few usable animals of a sex no between-sex comparison is possible.",
                    fix="Add deeper-sequenced animals, or lower min_dp in the config.")
            _rewrite_preflight(outdir, F, unit, tables)
            raise PreflightBlocked(F)
        F.info("contrasts_off", f"Only {n_m} males / {n_f} females keep enough genotypes: the sex contrasts (association, "
               "Fst, heterozygosity) are skipped and the result rests on the transmission test, which needs only depth >= "
               f"{cfg.ped_min_dp} in offspring. Replication across independent families is still assessed.")
    log(f"{unit.species}: {n_m} M / {n_f} F in the sex contrasts, {len(unit.parent_idx)} parents, "
        f"{len(unit.families)} families with a transmission test; {sum(f.severity == 'WARN' for f in F)} warnings")

    assoc = pop = None
    if contrasts_ok:
        try:
            assoc = AssociationAccumulator(unit, cfg)
        except ValueError as e:
            F.warn("association_off", str(e))
        if assoc is not None and assoc.skipped:
            F.info("strata_skipped", "Strata too small / single-sex for the association test: " + ", ".join(assoc.skipped))
        if min(n_m, n_f) >= cfg.min_per_sex_called:
            pop = PopgenAccumulator(unit, cfg)
    ped = PedigreeAccumulator(unit, cfg) if unit.families else None
    if assoc is None and pop is None and ped is None:
        raise ValueError("No analysis can run: no usable sex contrast and no usable pedigree")
    qc = QCAccumulator(unit, cfg)
    cmods = [m for m in (assoc, pop) if m is not None]       # strict genotypes, sites need min_call_rate
    rmods = [m for m in (ped, qc) if m is not None]          # role-specific genotypes

    n_seen = n_var = n_cr_seen = n_cr_dropped = 0
    all_idx = list(range(len(vcf_samples)))
    for chunk in iter_chunks(vcf, all_idx, vcf_samples, cfg, vf):
        if len(chunk.pos) == 0:
            continue
        n_seen += len(chunk.pos)
        chunk.dos_role = mask_by_sample(chunk.dos_raw, chunk.dp, chunk.gq, min_dp_v, min_gq_v)
        for m in rmods:
            m.update(chunk)
        if cmods:
            ok = (chunk.dos[:, unit.contrast_idx] >= 0).mean(axis=1) >= cfg.min_call_rate
            n_cr_seen += len(ok)
            n_cr_dropped += int((~ok).sum())
            sub = chunk if ok.all() else chunk.subset(ok)
            if len(sub.pos):
                n_var += len(sub.pos)
                for m in cmods:
                    m.update(sub)
    log(f"processed {n_seen:,} sites ({n_var:,} used by the sex contrasts)")
    if n_seen == 0:
        raise ValueError("No variants passed filters (check --contigs / QC thresholds)")
    if cmods and n_var == 0 and ped is None:
        raise ValueError("No site passed the call-rate filter for the sex contrasts (see min_call_rate / min_dp)")

    vs = vf.summary()
    if cfg.apply_variant_filter and vs["seen"]:
        top = ", ".join(f"{k}: {v:,}" for k, v in sorted(vs["reasons"].items(), key=lambda kv: -kv[1])[:6])
        F.info("variant_filter", f"Variant filter kept {vs['kept']:,} of {vs['seen']:,} SNPs "
               f"({vs['removed_fraction']:.1%} removed)" + (f"; reasons: {top}" if top else ""))
        if vs["removed_fraction"] > cfg.max_filtered_fraction_warn:
            F.warn("heavy_variant_filtering", f"{vs['removed_fraction']:.0%} of sites were removed by the variant "
                   "filter; the thresholds may not suit this VCF (run `sexloci describe` to inspect INFO fields).")
    cr_frac = n_cr_dropped / n_cr_seen if n_cr_seen else 0.0
    if n_cr_seen:
        F.info("call_rate_filter", f"Call-rate filter (sex contrasts only) removed {n_cr_dropped:,} of {n_cr_seen:,} sites "
               f"({cr_frac:.1%}) because fewer than {cfg.min_call_rate:.0%} of the {len(unit.contrast_idx)} contrast samples "
               f"had a genotype (depth >= {cfg.min_dp}); {n_var:,} sites were analysed. The transmission test is not "
               "subject to this filter.")
        if cr_frac > cfg.max_call_rate_drop_warn:
            F.warn("heavy_call_rate_filtering", f"{cr_frac:.0%} of sites were removed by the call-rate filter.",
                   why="Too many samples lack a genotype at most sites; sex-linked sites are removed with the rest.",
                   fix="Lower min_call_rate or min_dp in the config, or remove the most incomplete samples (sample_qc.tsv).")
    if n_seen < 1000:
        F.warn("few_variants", f"Only {n_seen:,} variants remain after filtering; mapping resolution will be poor.")
    sample_qc = qc_findings(qc, tables, cfg, F)
    F.acknowledge(cfg.acknowledged_codes)
    a_win, a_hits, a_info, strata_tables = assoc.finalize() if assoc else (pd.DataFrame(), pd.DataFrame(), {}, {})
    p_win = pop.finalize() if pop else pd.DataFrame()
    ped_win, fam_tables, ped_info = (ped.finalize() if ped else (pd.DataFrame(), {}, {}))
    # two primary combined tests (association, transmission) -> Bonferroni so the dataset-level FWER stays at alpha
    cfg_p = dataclasses.replace(cfg, alpha=cfg.alpha / (2 if (ped is not None and assoc is not None) else 1))
    win = flag_evidence(merge_windows(a_win, p_win, ped_win, window_size=cfg.window_size), cfg_p)
    gff = cfg.gff or tables.gff
    genes, gff_labels = None, {}
    if gff and os.path.exists(gff):
        genes, gff_labels = read_gff_genes(gff), read_gff_regions(gff)
    labels, aliases = match_annotation(sorted(win.chrom.unique()) if len(win) else [], cfg, gff_labels)
    cfg = dataclasses.replace(cfg, contig_labels=labels)                # plots and reports use the merged labels
    if aliases:
        genes = alias_genes(genes, aliases)
        F.info("annotation_matched_by_linkage_group", f"{len(aliases)} chromosomes of this VCF ({', '.join(list(aliases)[:3])}"
               f"{', ...' if len(aliases) > 3 else ''}) are not named in the annotation; genes are listed by matching linkage groups "
               "(e.g. GT3 sequence called under UMD2a chromosome names). This assumes the VCF coordinates are those of the annotated "
               "assembly.")
    elif genes is None:
        F.info("no_annotation", "No gene annotation file was available: regions do not list genes.")
    regions = call_regions(win, a_hits, cfg_p, unit, strata_tables, fam_tables, F.tier_cap(), genes,
                           alpha_unit=cfg.alpha)

    j = lambda name: os.path.join(outdir, name)
    import glob                              # results of an earlier run in this folder must not be mistaken for this run's
    for pattern in ("region_*.png", "region_*.pdf", "genome_haplotypes.*", "genome_by_chromosome.*", "family_transmission_stats.tsv",
                    "pedigree_chromosome_summary.tsv", "stratum_window_stats.tsv"):
        for old in glob.glob(os.path.join(outdir, pattern)):
            os.remove(old)
    win.to_csv(j("window_stats.tsv"), sep="\t", index=False)
    a_hits.to_csv(j("sex_associated_markers.tsv"), sep="\t", index=False)
    regions.to_csv(j("sex_linked_regions.tsv"), sep="\t", index=False)
    sample_qc.to_csv(j("sample_qc.tsv"), sep="\t", index=False)
    pd.DataFrame(sorted(list(vs["reasons"].items()) + [(f"LOW_CALL_RATE<{cfg.min_call_rate:g}", n_cr_dropped)],
                        key=lambda kv: -kv[1]), columns=["reason", "sites"]).assign(
        total_seen=vs["seen"], kept=vs["kept"]).to_csv(j("variant_filter_summary.tsv"), sep="\t", index=False)
    findings_frame(F).to_csv(j("findings.tsv"), sep="\t", index=False)
    long = [t.assign(stratum=k) for k, t in strata_tables.items()]
    if long:
        pd.concat(long).to_csv(j("stratum_window_stats.tsv"), sep="\t", index=False)
    if fam_tables:
        pd.concat([t.assign(family=k) for k, t in fam_tables.items()]).to_csv(
            j("family_transmission_stats.tsv"), sep="\t", index=False)
    if ped is not None and len(ped_win):
        pedigree_chromosome_summary(win).to_csv(j("pedigree_chromosome_summary.tsv"), sep="\t", index=False)
    thr = a_info.get("threshold_chi2", np.nan)
    genome_plot(win, regions, thr, cfg, j("genome_overview.png"))
    if ped is not None and ped.tracks:
        haplotype_genome_plot(ped.tracks, ped.family_meta, win, regions, cfg, j("genome_haplotypes.png"))
        chromosome_summary_plot(fam_tables, ped.family_meta, cfg, j("genome_by_chromosome.png"))
    for i, r in regions[regions.tier.isin(["A", "B"])].head(8).iterrows():
        region_plot(win, r, thr, cfg, j(f"region_{i + 1}_{r.chrom}.png"))

    out = dict(version=__version__, runtime_sec=round(time.time() - t0, 1), vcf=vcf,
               vcf_size_bytes=os.path.getsize(vcf), vcf_modified=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(vcf))),
               genome_version=tables.genome_version, n_variants_seen=n_seen, n_variants_used=n_var, variant_filter=vs, call_rate_filter=dict(seen=n_cr_seen, dropped=n_cr_dropped),
               samples=dict(n_male=n_m, n_female=n_f, n_parents=int(len(unit.parent_idx))),
               assoc=a_info, pedigree=ped_info, modules=[m.name for m in (assoc, pop, ped) if m is not None],
               n_regions=int(len(regions)), tier_cap=F.tier_cap(), config=config_dict(cfg),
               findings=[dict(severity=f.severity, code=f.code, message=f.message) for f in F.sorted()])
    with open(j("run_summary.json"), "w") as fh:
        json.dump(out, fh, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    write_report(j("report.md"), out, regions, F, unit, cfg)
    out.update(regions=regions, windows=win, findings_obj=F, unit=unit,
               haplotype_tracks=(ped.tracks if ped is not None else {}))
    return out
