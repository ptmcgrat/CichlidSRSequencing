#!/usr/bin/env python3
"""gsdf_cn.py - genotype the gsdf-region duplication (copy number) from short reads aligned to the M. zebra GT3a reference.

The reference carries TWO tandem copies of the 22 kb unit (NC_135173.1:50,751,000-50,794,727). Reads from every genomic copy map to
both reference copies, so a sample's depth over the unit relative to its own flanks is proportional to its total copy number:

    ratio = (total copies in the diploid genome) / 4          1 -> 0.25   2 -> 0.50   3 -> 0.75   4 -> 1.00 (like the reference)

In practice there is a small offset, so the script calibrates: by default on the MC family (daughters = 2 copies, sons = 3), or on
samples of known state you provide (--anchors). 2 = no extra copy, 3 = one chromosome carries the duplication (e.g. a male whose Y does),
4 = both chromosomes do. Depth alone cannot say WHICH chromosome carries it.

Input: the depth table made with samtools (columns sample, bin_start, mean_depth; 5 kb bins over NC_135173.1:50,700,000-50,830,000):

    printf "sample\\tbin_start\\tmean_depth\\n" > gsdf_depth.tsv
    for b in *_gsdf.bam; do s=${b%_gsdf.bam}
      samtools depth -a -r NC_135173.1:50700000-50830000 "$b" | awk -v s="$s" 'BEGIN{OFS="\\t"} {bin=int($2/5000)*5000; sum[bin]+=$3; n[bin]++} END{for(k in sum) print s,k,sum[k]/n[k]}' >> gsdf_depth.tsv
    done

Usage:  python3 gsdf_cn.py gsdf_depth.tsv [--sex sex.tsv] [--anchors anchors.tsv] [--out gsdf_cn]
  sex.tsv / anchors.tsv: tab-separated with a header; sex.tsv has columns sample,sex; anchors.tsv has sample,copies (a few samples of
  known copy number from at least two different classes, e.g. 2 and 3).
Writes <out>.tsv (ratio, estimated copies, call, flag) and <out>.png.
"""
import argparse
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

UNIT = (50_755_000, 50_785_000)          # bins lying fully inside the two copies
FLANK = ((0, 50_735_000), (50_800_000, 50_815_000))   # outside the block, before the repeat-rich bump at 50.825 Mb


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("depth")
    ap.add_argument("--sex")
    ap.add_argument("--anchors")
    ap.add_argument("--theory", action="store_true", help="no calibration: ratio = copies / 4")
    ap.add_argument("--min-flank-depth", type=float, default=3.0)
    ap.add_argument("--out", default="gsdf_cn")
    a = ap.parse_args()

    d = pd.read_csv(a.depth, sep="\t")
    W = d.pivot_table(index="sample", columns="bin_start", values="mean_depth")
    unit = [b for b in W.columns if UNIT[0] <= b <= UNIT[1]]
    flank = [b for b in W.columns if any(lo <= b <= hi for lo, hi in FLANK)]
    if len(unit) < 5 or len(flank) < 6:
        sys.exit(f"the depth table needs 5 kb bins covering {UNIT[0]:,}-{UNIT[1]:,} and the flanks (found {len(unit)} unit and {len(flank)} flank bins)")
    R = pd.DataFrame({"flank_depth": W[flank].mean(axis=1), "unit_depth": W[unit].mean(axis=1)})
    R["ratio"] = R.unit_depth / R.flank_depth

    # calibration: ratio = r2 + (copies - 2) * step
    if a.theory:
        r2, step, how = 0.50, 0.25, "theory (copies/4)"
    elif a.anchors:
        A = pd.read_csv(a.anchors, sep="\t").set_index("sample").copies
        A = A[A.index.isin(R.index)]
        if A.nunique() < 2:
            sys.exit("--anchors needs samples of at least two different copy numbers that are present in the depth table")
        step, icpt = np.polyfit(A.values.astype(float), R.loc[A.index, "ratio"].values, 1)
        r2, how = icpt + 2 * step, f"fitted on {len(A)} anchor samples"
    else:
        r2, step, how = 0.530, 0.268, "MC family (daughters = 2 copies, sons = 3), bwa/GT3a, 4-8x"
    R["est_copies"] = 2 + (R.ratio - r2) / step
    R["call"] = np.clip(np.rint(R.est_copies), 1, 6).astype(int)
    R["flag"] = ""
    R.loc[(R.est_copies - R.call).abs() > 0.35, "flag"] = "uncertain"
    R.loc[R.flank_depth < a.min_flank_depth, "flag"] = "low depth"
    if a.sex:
        S = pd.read_csv(a.sex, sep="\t").set_index("sample").sex
        R["sex"] = S.reindex(R.index)
    R = R.sort_values("ratio")
    R.round(3).to_csv(a.out + ".tsv", sep="\t")
    print(f"calibration: {how}  (ratio at 2 copies {r2:.3f}, per extra copy {step:.3f})")
    print(f"{len(R)} samples; unit bins {len(unit)}, flank bins {len(flank)}")
    summ = R.groupby(["call"] + (["sex"] if a.sex else [])).size().rename("n").reset_index()
    print("calls (total copies in the diploid genome):\n" + summ.to_string(index=False))
    print(f"flagged: {(R.flag == 'uncertain').sum()} uncertain, {(R.flag == 'low depth').sum()} low depth")

    fig, ax = plt.subplots(figsize=(11, 4.4))
    col = {"M": "#4393c3", "F": "#d6604d"}
    x = np.arange(len(R))
    ax.scatter(x, R.ratio, s=28, c=[col.get(s, "0.4") for s in R.sex] if a.sex else "0.3", edgecolor="k", lw=0.3)
    for k in range(1, 6):
        ax.axhline(r2 + (k - 2) * step, c="0.8", lw=0.8, ls="--")
        ax.text(len(R) + 0.3, r2 + (k - 2) * step, f"{k} copies", va="center", fontsize=8, color="0.45")
    ax.set_xlabel("samples, sorted by ratio" + ("   (blue: male, red: female)" if a.sex else ""))
    ax.set_ylabel("depth over unit / depth in flanks")
    ax.set_title("gsdf-region copy number from read depth", fontsize=10)
    ax.set_xlim(-1, len(R) + 4)
    fig.tight_layout()
    fig.savefig(a.out + ".png", dpi=150)
    print(f"wrote {a.out}.tsv and {a.out}.png")


if __name__ == "__main__":
    main()
