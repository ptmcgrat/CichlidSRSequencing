#!/usr/bin/env python3
"""chr7ins_combine.py - one chr7-insertion genotype per sample from the two reference copies of the site.

`sexloci sv-genotype` genotypes chr7ins_c1 and chr7ins_c2 separately. A carrier has one copy of the gsdf unit per chromosome and GT3a has two, so the
carrier's reads at this locus split between the two reference copies. This adds the informative pairs of the two sites and calls the genotype once from
the total, with exactly the model sexloci uses (P(insertion read) = 0.02 for 0/0, 0.5 for 0/1, 0.98 for 1/1; GQ is the phred-scaled gap between the
best and second-best genotype; fewer than 5 informative pairs gives ./.).

    python3 chr7ins_combine.py --genotypes ~/Temp/SexDeterminationAnalysis/sv_chr7ins_alleles/sv_genotypes.tsv

Writes chr7ins_combined.tsv with SampleID, chr7ins_GT, chr7ins_ref, chr7ins_alt, chr7ins_GQ, chr7ins_reads, chr7ins_status and the per-copy counts. The
columns follow the sexloci naming, so chr7ins_concordance.py --genotypes chr7ins_combined.tsv --name chr7ins compares them with Munby's calls.
"""
import argparse
import math
import os

import pandas as pd


def call(n_ref, n_alt, min_pairs=5, false_alt=0.02, false_ref=0.02):
    """The sexloci model, from counts. Returns (GT, GQ)."""
    if n_ref + n_alt < min_pairs:
        return "./.", 0.0
    lls = {}
    for g, p in (("0/0", false_alt), ("0/1", 0.5), ("1/1", 1.0 - false_ref)):
        p = min(max(p, 1e-12), 1 - 1e-12)
        lls[g] = n_alt * math.log(p) + n_ref * math.log(1 - p)
    best = max(lls, key=lls.get)
    s = sorted(lls.values(), reverse=True)
    return best, (s[0] - s[1]) * (10.0 / math.log(10))


def combine(df, sites, name="chr7ins", min_pairs=5):
    num = lambda c: pd.to_numeric(df[c], errors="coerce") if c in df.columns else pd.Series(0.0, index=df.index)
    ok = pd.concat([(df[f"{s}_status"] == "OK") if f"{s}_status" in df.columns else pd.Series(False, index=df.index) for s in sites], axis=1)
    n_ref = sum(num(f"{s}_ref").fillna(0) for s in sites)
    n_alt = sum(num(f"{s}_alt").fillna(0) for s in sites)
    out = df[[c for c in ("SampleID", "Species", "Ecogroup") if c in df.columns]].copy()
    gts = [call(r, a, min_pairs) if o else ("", float("nan")) for r, a, o in zip(n_ref, n_alt, ok.any(axis=1))]
    out[f"{name}_GT"] = [g for g, _ in gts]
    out[f"{name}_ref"] = n_ref
    out[f"{name}_alt"] = n_alt
    out[f"{name}_GQ"] = [round(q, 1) if q == q else q for _, q in gts]
    out[f"{name}_reads"] = sum(num(f"{s}_reads").fillna(0) for s in sites)
    out[f"{name}_status"] = ["OK" if o else "no_data" for o in ok.any(axis=1)]
    for s in sites:
        out[f"{s}_ref"], out[f"{s}_alt"] = num(f"{s}_ref"), num(f"{s}_alt")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", default="~/Temp/SexDeterminationAnalysis/sv_chr7ins_alleles/sv_genotypes.tsv")
    ap.add_argument("--sites", nargs="+", default=["chr7ins_c1", "chr7ins_c2"])
    ap.add_argument("--name", default="chr7ins")
    ap.add_argument("--min-pairs", type=int, default=5)
    ap.add_argument("--min-gq", type=float, default=10.0)
    ap.add_argument("--out", default="chr7ins_combined.tsv")
    a = ap.parse_args()
    df = pd.read_csv(os.path.expanduser(a.genotypes), sep="\t", dtype=str)
    out = combine(df, a.sites, a.name, a.min_pairs)
    out.to_csv(a.out, sep="\t", index=False)
    gt = out[f"{a.name}_GT"]
    conf = out[gt.isin(["0/0", "0/1", "1/1"]) & (pd.to_numeric(out[f"{a.name}_GQ"], errors="coerce") >= a.min_gq)]
    print(f"{len(out)} samples; {int((gt != '').sum())} with data; {len(conf)} confident calls (GQ >= {a.min_gq:g}); ./. for {int((gt == './.').sum())}")
    if "Ecogroup" in conf.columns:
        print(pd.crosstab(conf["Ecogroup"].replace("", "(none)"), conf[f"{a.name}_GT"]).to_string())
    else:
        print(conf[f"{a.name}_GT"].value_counts().to_string())
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
