"""`describe` (inspect a VCF's filter annotations) and `filter-compare` (raw vs already-filtered VCF)."""
from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

from .vcf import VariantFilter, _contig_filter, _open, read_header_info, scan_sites


def describe_vcf(path, cfg) -> dict:
    hdr = read_header_info(path)
    vf = VariantFilter(cfg, hdr["infos"])
    labels = scan_sites(path, cfg, vf, record_dp=True)
    s = vf.summary()
    n = sum(labels.values())
    n_pass = labels.get("PASS", 0) + labels.get(".", 0)
    non_pass = {k: v for k, v in labels.items() if k not in ("PASS", ".")}
    if non_pass:
        kind = ("FILTER column is annotated: " + f"{n - n_pass:,} of {n:,} sites carry a non-PASS label "
                f"({', '.join(sorted(non_pass)[:5])}); the pipeline drops them")
    elif labels.get("PASS", 0) == n and n:
        kind = "every site is PASS: the file looks already filtered (bad sites removed) or never annotated"
    else:
        kind = "FILTER column is '.' everywhere: the file looks RAW (unfiltered)"
    dp = np.asarray(vf.dp_values, float)
    return dict(path=path, n_samples=len(hdr["samples"]), n_snps=n, kind=kind, filter_labels=dict(labels),
                info_fields_present=sorted(hdr["infos"] & ({h[0] for h in vf.hard} | {"DP"})),
                hard_filters_applicable=[vf.label(h) for h in vf.applicable()],
                hard_filters_skipped=[vf.label(h) for h in vf.not_applicable()],
                would_remove=s["removed"], would_remove_fraction=round(s["removed_fraction"], 4),
                per_reason=s["reasons"],
                info_dp_quantiles={q: float(np.quantile(dp, q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95, 0.99)}
                if len(dp) else {})


def _sites(path, cfg, keep, vf=None):
    """Yield (chrom, pos, key, passed_by_label, our_keep) for biallelic SNP records."""
    with _open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if not keep(f[0]):
                continue
            ref, alt = f[3], f[4]
            if "," in alt or alt in (".", "*") or len(ref) != 1 or len(alt) != 1:
                continue
            label_pass = f[6] in ("PASS", ".")
            ours = vf.accept(f[6], f[5], f[7]) if vf is not None else None
            yield f[0], int(f[1]), hash((f[0], f[1], ref, alt)), label_pass, ours


def compare_filters(raw, filtered, cfg, window=None) -> dict:
    """How well does our site filter applied to the RAW VCF reproduce the sites kept in the FILTERED VCF?
    Also reports where the filtered VCF lost sites (per contig / window) - e.g. sex-linked regions that a
    population-genomic filter (InbreedingCoeff, depth) would delete."""
    keep = _contig_filter(cfg)
    window = window or cfg.window_size
    hdr = read_header_info(raw)
    vf = VariantFilter(cfg, hdr["infos"])
    kept_f = np.array([k for _, _, k, ok, _ in _sites(filtered, cfg, keep) if ok], dtype=np.int64)
    rk, ours, ck, cp = [], [], [], []
    for chrom, pos, k, _, o in _sites(raw, cfg, keep, vf):
        rk.append(k); ours.append(bool(o)); ck.append(chrom); cp.append(pos // window)
    rk, ours = np.array(rk, dtype=np.int64), np.array(ours)
    in_f = np.isin(rk, kept_f)
    both = int((ours & in_f).sum())
    res = dict(raw_sites=int(len(rk)), filtered_file_kept_sites=int(len(kept_f)),
               filtered_sites_found_in_raw=int(in_f.sum()), our_filter_kept=int(ours.sum()),
               recall_of_filtered_file=both / max(int(in_f.sum()), 1),
               precision_vs_filtered_file=both / max(int(ours.sum()), 1),
               ours_only=int((ours & ~in_f).sum()), filtered_file_only=int((~ours & in_f).sum()))
    df = pd.DataFrame(dict(chrom=ck, window=cp, raw=1, kept_filtered_file=in_f.astype(int), kept_ours=ours.astype(int)))
    by_contig = df.groupby("chrom")[["raw", "kept_filtered_file", "kept_ours"]].sum()
    by_contig["retained_filtered_file"] = by_contig.kept_filtered_file / by_contig.raw
    by_contig["retained_ours"] = by_contig.kept_ours / by_contig.raw
    by_win = df.groupby(["chrom", "window"])[["raw", "kept_filtered_file", "kept_ours"]].sum().reset_index()
    by_win["retained_filtered_file"] = by_win.kept_filtered_file / by_win.raw
    by_win["retained_ours"] = by_win.kept_ours / by_win.raw
    med = by_win.loc[by_win.raw >= 20, "retained_filtered_file"].median()
    by_win["depleted_in_filtered_file"] = (by_win.raw >= 20) & (by_win.retained_filtered_file < 0.5 * med)
    res["median_window_retention_filtered_file"] = float(med) if np.isfinite(med) else None
    res["n_depleted_windows"] = int(by_win.depleted_in_filtered_file.sum())
    return dict(summary=res, by_contig=by_contig.reset_index(), by_window=by_win)


def audit_filter_by_window(path, cfg, window=None) -> pd.DataFrame:
    """Where along the genome does the site filter bite? Per window: sites seen, sites removed and, for every reason, the
    number of sites it flags (a site can have several reasons). The INFO/DP quantile bound is not included."""
    keep = _contig_filter(cfg)
    window = window or cfg.window_size
    hdr = read_header_info(path)
    vf = VariantFilter(cfg, hdr["infos"])
    rows = {}
    with _open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if not keep(f[0]):
                continue
            ref, alt = f[3], f[4]
            if "," in alt or alt in (".", "*") or (cfg.snps_only and (len(ref) != 1 or len(alt) != 1)):
                continue
            reasons, _ = vf.reasons(f[6], f[5], f[7]) if vf.enabled else ([], None)
            d = rows.setdefault((f[0], (int(f[1]) - 1) // window), Counter())
            d["sites"] += 1
            if reasons:
                d["removed"] += 1
            for r in reasons:
                d[r] += 1
    df = pd.DataFrame([dict(chrom=c, widx=i, **v) for (c, i), v in rows.items()]).fillna(0)
    df["start"] = df.widx * window + 1
    df["end"] = (df.widx + 1) * window
    df["retained"] = 1 - df.removed / df.sites
    return df.sort_values(["chrom", "widx"]).reset_index(drop=True)
