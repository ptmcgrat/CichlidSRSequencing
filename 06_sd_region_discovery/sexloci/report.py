"""Markdown / TSV reporting of findings and called regions."""
from __future__ import annotations

import pandas as pd

from .findings import BLOCK, Findings


def findings_frame(F: Findings) -> pd.DataFrame:
    cols = ["severity", "code", "message", "why", "fix", "n_samples", "samples", "tier_cap"]
    return pd.DataFrame([dict(severity=f.severity, code=f.code, message=f.message, why=f.why, fix=f.fix,
                              n_samples=len(f.samples), samples=";".join(f.samples), tier_cap=f.cap or "")
                         for f in F.sorted()], columns=cols)


def fix_list_frame(F: Findings) -> pd.DataFrame:
    """One row per affected sample: what to change in the sample sheet."""
    rows = [dict(problem=f.code, **r) for f in F.sorted() for r in f.rows]
    return pd.DataFrame(rows, columns=["problem", "sample", "column", "current", "needed"])


def _findings_md(F: Findings):
    if not len(F):
        return ["No warnings.", ""]
    L = []
    for sev, title in (("BLOCK", "Blocking problems"), ("WARN", "Warnings"), ("INFO", "Notes")):
        items = [f for f in F.sorted() if f.severity == sev]
        if not items:
            continue
        L += [f"### {title}", ""]
        for f in items:
            L.append(f"- **{f.code}** - {f.message}")
            if f.why:
                L.append(f"  - *Why:* {f.why}")
            if f.fix:
                L.append(f"  - *How to fix:* {f.fix}")
            if f.rows:
                L.append(f"  - *Affected samples:* {len(f.rows)} (listed in `fix_list.tsv`)")
        L.append("")
    return L


def write_preflight_report(path, F: Findings, unit=None, tables=None):
    L = ["# sexloci preflight report", ""]
    if F.blocked:
        L += ["**The run was stopped: fix the blocking problems below and rerun.**", ""]
    else:
        L += ["No blocking problems.", ""]
    L += _findings_md(F)
    with open(path, "w") as fh:
        fh.write("\n".join(L))


def write_report(path, out, regions, F: Findings, unit, cfg):
    sm = out["samples"]
    L = ["# sexloci report", "",
         f"Species: **{unit.species}**  |  genome: {out['genome_version']}  |  {sm['n_male']} males, "
         f"{sm['n_female']} females analysed (+{sm['n_parents']} parents)  |  {out['n_variants_used']:,} variants  |  "
         f"modules: {', '.join(out['modules'])}", ""]
    cap = F.tier_cap()
    if cap:
        L += [f"> Confidence is capped at **tier {cap}** by the data-quality warnings below.", ""]
    L += ["## Data-quality findings", ""] + _findings_md(F)
    L += ["## Candidate sex-linked regions", ""]
    if regions.empty:
        L += ["No region passed the evidence thresholds. That is expected if sex is not controlled by a detectable "
              "genetic locus here, if sex determination is polygenic or environmental, or if the cohort is too small. "
              "Inspect `window_stats.tsv` and `genome_overview.png`.", ""]
    else:
        L += ["Tier A: >= 2 independent units (families sharing no parent / the pool of non-family animals) each "
              "significant (replication). B: the combined test is significant but replication is not shown. "
              "C: population-genetic outlier only (a lead, not a result).", "",
              "| tier | region (transmission plateau) | differentiated core (male-female heterozygosity) | model | independent units | evidence | top SNP |",
              "|---|---|---|---|---|---|---|"]
        for _, r in regions.iterrows():
            top = "" if pd.isna(r.top_marker_pos) else f"{int(r.top_marker_pos):,}"
            note = f" (raw {r.tier_raw})" if r.tier != r.tier_raw else ""
            core = "" if pd.isna(r.core_start) else f"{int(r.core_start):,}-{int(r.core_end):,} ({(r.core_end - r.core_start + 1) / 1e6:.2f} Mb)"
            L.append(f"| {r.tier}{note} | {r.chrom}:{r.start:,}-{r.end:,} ({r.size_bp / 1e6:.2f} Mb) | {core} | {r.model} | "
                     f"{r.units_supporting}/{r.units_total} | {r.evidence} | {top} |")
            if r.flag:
                L.append(f"|  | {r.flag} | | | | | |")
            if isinstance(r.units_share_parent, str) and r.units_share_parent:
                L.append(f"|  | note: the supporting broods share {r.units_share_parent}, so this replication is across different "
                         "partners of the same parent, not across unrelated parents | | | | | |")
        L.append("")
    L += ["## Reading the results", "",
          "- `XY`: the sire-inherited haplotype predicts sex and/or males are more heterozygous; `ZW` the reverse.",
          "- In small pedigrees linkage spans much of a chromosome: `start-end` is the peak plateau, "
          "`linked_start-linked_end` the wider significant interval, `sig_marker_span` the span of individually "
          "significant SNPs.",
          "- Association and Fst/heterozygosity are the same information and never count as separate confirmation; "
          "only independent units do.",
          "- No inversion is assumed. `extended_block` only reflects size.",
          "- Y/W-specific sequence absent from the reference cannot appear as SNPs; see the depth panel.", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(L))
