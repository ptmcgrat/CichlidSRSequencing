#!/usr/bin/env python3
"""chr7ins_concordance.py - does the candidate junction track Munby's chr7-ins calls?

Compares the read counts that `sexloci sv-genotype` produced at the candidate site (columns chr7ins_cand_ref / _alt / _GT in sv_genotypes.tsv) with
Munby's genotypes for the same fish (chr7ins_validation_calls.tsv). It works from the counts, not the GT calls, because at about 17x depth the
site often has too few informative pairs for a confident call.

    python3 chr7ins_concordance.py
    python3 chr7ins_concordance.py --genotypes ~/Temp/SexDeterminationAnalysis/sv_chr7ins_candidate/sv_genotypes.tsv --truth chr7ins_validation_calls.tsv

What it reports
  * a per-fish table (Munby's call, our pairs, our call), also written to chr7ins_concordance.tsv
  * per group (homozygote / heterozygote / non-carrier): how many fish show insertion-allele pairs and how many
  * the probability that a random carrier has more insertion pairs than a random non-carrier (AUC; 0.5 = no relationship, 1 = perfect), with a permutation p
  * whether homozygotes show more insertion pairs than heterozygotes, and whether our insertion fraction rises with Munby's
  * the same split inside Lake Kingiri, where carriers and non-carriers share a population
  * a plain-language reading. The thresholds in it are rules of thumb, not a statistical test.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

DEFAULT_GENOTYPES = "~/Temp/SexDeterminationAnalysis/sv_chr7ins_candidate/sv_genotypes.tsv"


def auc(pos, neg):
    """P(a random positive scores higher than a random negative); ties count half."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return (gt + 0.5 * eq) / (len(pos) * len(neg))


def perm_p(pos, neg, n=10000, seed=1):
    """One-sided permutation p for AUC >= the observed value."""
    rng = np.random.default_rng(seed)
    obs = auc(pos, neg)
    allv = np.concatenate([pos, neg])
    k = len(pos)
    hit = 0
    for _ in range(n):
        rng.shuffle(allv)
        hit += auc(allv[:k], allv[k:]) >= obs - 1e-12
    return (hit + 1) / (n + 1)


def rank(x):
    return pd.Series(x).rank().to_numpy()


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 4 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(rank(x), rank(y))[0, 1])


def build_table(gen: pd.DataFrame, truth: pd.DataFrame, name: str) -> pd.DataFrame:
    for c in (f"{name}_ref", f"{name}_alt", f"{name}_GT"):
        if c not in gen.columns:
            sys.exit(f"column {c} not found in the genotypes file; the site name differs? use --name (columns present: "
                     f"{[x for x in gen.columns if 'chr7' in x or 'cand' in x][:6]})")
    g = gen[["SampleID", f"{name}_GT", f"{name}_ref", f"{name}_alt"] + [c for c in (f"{name}_GQ", "Coverage", "Status") if c in gen.columns]].copy()
    g = g.rename(columns={f"{name}_GT": "our_GT", f"{name}_ref": "our_ref", f"{name}_alt": "our_alt", f"{name}_GQ": "our_GQ"})
    for c in ("our_ref", "our_alt"):
        g[c] = pd.to_numeric(g[c], errors="coerce")
    t = truth.merge(g, left_on="BioSample", right_on="SampleID", how="left")
    t["our_pairs"] = t.our_ref + t.our_alt
    t["our_alt_fraction"] = np.where(t.our_pairs > 0, t.our_alt / t.our_pairs, np.nan)
    return t


def summarise(t: pd.DataFrame) -> str:
    out = []
    have = t[t.our_alt.notna()]
    missing = t[t.our_alt.isna()]
    out.append(f"{len(have)} of {len(t)} fish have a result" + (f"; no result yet for {len(missing)}: {', '.join(missing.BioSample.head(8))}"
                                                                    f"{' ...' if len(missing) > 8 else ''}" if len(missing) else ""))
    out.append("")
    out.append("by Munby's call (pairs are informative read pairs at the candidate junction):")
    out.append(f"   {'group':12s} {'n':>3s}  {'with >=1 alt pair':>18s}  {'with >=3':>9s}  {'median alt':>10s}  {'median ref':>10s}  our calls (0/0, 0/1, 1/1, ./.)")
    for grp in ("hom", "het", "non-carrier", "unknown"):
        d = have[have.group == grp]
        if not len(d):
            continue
        calls = d.our_GT.value_counts()
        out.append(f"   {grp:12s} {len(d):3d}  {int((d.our_alt >= 1).sum()):18d}  {int((d.our_alt >= 3).sum()):9d}  {d.our_alt.median():10.1f}  {d.our_ref.median():10.1f}  "
                   f"{calls.get('0/0', 0)}, {calls.get('0/1', 0)}, {calls.get('1/1', 0)}, {calls.get('./.', 0)}")
    car = have[have.group.isin(["hom", "het"])]
    non = have[have.group == "non-carrier"]
    out.append("")
    verdict = []
    if len(car) and len(non):
        a = auc(car.our_alt, non.our_alt)
        p = perm_p(car.our_alt.to_numpy(float), non.our_alt.to_numpy(float))
        out.append(f"carriers vs non-carriers, insertion pairs: AUC {a:.2f} (permutation p = {p:.4f}); "
                   f"non-carriers with >=2 insertion pairs: {int((non.our_alt >= 2).sum())} of {len(non)}; carriers with >=2: {int((car.our_alt >= 2).sum())} of {len(car)}")
        verdict.append(("auc", a, int((non.our_alt >= 2).sum()), len(non)))
        k = have[have.location == "Lake_Kingiri"]
        kc, kn = k[k.group.isin(["hom", "het"])], k[k.group == "non-carrier"]
        if len(kc) and len(kn):
            out.append(f"   inside Lake Kingiri only: AUC {auc(kc.our_alt, kn.our_alt):.2f} ({len(kc)} carriers, {len(kn)} non-carriers); "
                       f"non-carriers with >=1 insertion pair: {int((kn.our_alt >= 1).sum())} of {len(kn)}")
    hom, het = have[have.group == "hom"], have[have.group == "het"]
    if len(hom) and len(het):
        out.append(f"homozygotes vs heterozygotes: median insertion pairs {hom.our_alt.median():.1f} (n={len(hom)}) vs {het.our_alt.median():.1f} (n={len(het)}); "
                   f"mean insertion fraction {hom.our_alt_fraction.mean():.2f} vs {het.our_alt_fraction.mean():.2f}")
        verdict.append(("dose", hom.our_alt.median(), het.our_alt.median()))
    enough = car[(car.our_pairs >= 3) & car.Munby_alt_fraction.astype(str).ne("")]
    if len(enough) >= 5:
        r = spearman(pd.to_numeric(enough.Munby_alt_fraction), enough.our_alt_fraction)
        out.append(f"carriers with >=3 informative pairs (n={len(enough)}): Spearman correlation between Munby's insertion fraction and ours = {r:.2f}")
    out.append("")
    out.append("reading:")
    thr = 2                                                  # a fish counts as showing the junction with >= 2 insertion pairs (about the noise floor at this depth)
    if len(car) and len(non):
        sens = float((car.our_alt >= thr).mean())
        fpr = float((non.our_alt >= thr).mean())
        k = have[have.location == "Lake_Kingiri"]
        kn = k[k.group == "non-carrier"]
        kn_pos = int((kn.our_alt >= thr).sum())
        other_non = non[non.location != "Lake_Kingiri"]
        other_pos = int((other_non.our_alt >= thr).sum())
        dose_ok = len(hom) and len(het) and hom.our_alt_fraction.mean() > het.our_alt_fraction.mean() and hom.our_alt.median() > het.our_alt.median()
        out.append(f"   {sens:.0%} of Munby's carriers show the junction (>= {thr} insertion pairs); {fpr:.0%} of his non-carriers do "
                   f"({kn_pos} of {len(kn)} in Kingiri, {other_pos} of {len(other_non)} elsewhere).")
        if sens >= 0.7 and fpr <= 0.1 and kn_pos <= 1:
            out.append("   The candidate junction tracks Munby's chr7 insertion calls: carriers show it, non-carriers, including those from the carriers' own "
                       "population, do not." + (" Homozygotes also carry more of it than heterozygotes, as they should." if dose_ok else
                                               " Homozygotes do not clearly exceed heterozygotes, so the dose is still unsettled."))
        elif len(kn) >= 3 and kn_pos >= 2 and other_pos <= max(1, round(0.1 * len(other_non))):
            out.append("   The junction is present in Munby's non-carriers from Lake Kingiri as well, but not elsewhere: it looks like a variant of that population, "
                       "not his chr7 insertion (or a different, nearby event). That would argue against the candidate.")
        elif sens < 0.3:
            out.append("   Few of Munby's carriers show the junction: at this depth it does not detect his carriers, so either the candidate is not his insertion "
                       "or the signal is too thin to test.")
        else:
            out.append("   Partial agreement. Look at which fish disagree in the table before concluding; low-depth fish and the read split between the two "
                       "reference copies can blur a real signal. More depth (or the second reference copy as a site) would sharpen it.")
    out.append("   With ~17x data and half the reads going to the other reference copy, expect about 5-7 informative insertion pairs per heterozygote and "
               "10-14 per homozygote, and 0 in non-carriers.")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", default=DEFAULT_GENOTYPES, help="sv_genotypes.tsv written by sexloci sv-genotype")
    ap.add_argument("--truth", default="chr7ins_validation_calls.tsv", help="Munby's calls for the fish (this file ships with the script)")
    ap.add_argument("--name", default="chr7ins_cand", help="site name used with --indel-names (default chr7ins_cand)")
    ap.add_argument("--out", default="chr7ins_concordance.tsv")
    a = ap.parse_args()
    gen = pd.read_csv(os.path.expanduser(a.genotypes), sep="\t", dtype=str)
    truth = pd.read_csv(os.path.expanduser(a.truth), sep="\t", dtype=str).fillna("")
    t = build_table(gen, truth, a.name)
    order = {"hom": 0, "het": 1, "non-carrier": 2, "unknown": 3}
    t = t.sort_values(by=["group", "location", "our_alt"], key=lambda c: c.map(order) if c.name == "group" else c, ascending=[True, True, False])
    show = t[["BioSample", "group", "location", "Munby_ref", "Munby_alt", "our_ref", "our_alt", "our_GT", "our_GQ"] if "our_GQ" in t.columns else
             ["BioSample", "group", "location", "Munby_ref", "Munby_alt", "our_ref", "our_alt", "our_GT"]]
    print(show.to_string(index=False, na_rep="-"))
    print()
    print(summarise(t))
    t.to_csv(a.out, sep="\t", index=False)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
