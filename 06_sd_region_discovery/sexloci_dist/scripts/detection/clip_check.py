#!/usr/bin/env python3
"""Count reads that are soft-clipped (>= 20 bp) right at TE1 / TE2 in the region BAMs that `sexloci sv-genotype` keeps.

An insertion that is not in the reference leaves a pile of soft-clipped reads at its position, whatever the inserted sequence is. So this
separates three situations that the genotype call alone cannot:
  * clean reference at the site: almost no clipped reads;
  * the insertion sexloci looks for: clipped reads at the site AND insertion pairs in the genotype (TE_alt > 0);
  * a DIFFERENT insertion or rearrangement at the same site: clipped reads at the site but no insertion pairs.

usage:  python3 clip_check.py ASW_615_1 ASW_557_1 AG_552_1 ASG_236_1
        (a TE-negative female and a TE-positive male make good controls; BAM folder: --bamdir)
needs samtools on the PATH; the BAMs need no index.
"""
import argparse, os, re, subprocess, sys

SITES = {"TE1": ("NC_135176.1", 24_921_178), "TE2": ("NC_135176.1", 24_921_394)}
OPS = re.compile(r"(\d+)([MIDNSHP=X])")


def clip_boundaries(pos, cigar, min_clip=20):
    """Reference positions where a soft clip of >= min_clip bases begins: (position, 'left'|'right')."""
    ops = [(int(n), o) for n, o in OPS.findall(cigar)]
    if not ops:
        return []
    out = []
    if ops[0][1] == "S" and ops[0][0] >= min_clip:
        out.append((pos, "left"))                                   # alignment starts here; the clipped part lies before it
    if ops[-1][1] == "S" and ops[-1][0] >= min_clip:
        out.append((pos + sum(n for n, o in ops if o in "MDN=X"), "right"))   # first reference base after the alignment
    return out


def count_clips(sam_lines, sites=SITES, window=25, min_clip=20):
    res = {k: dict(left=0, right=0) for k in sites}
    n_reads = {k: 0 for k in sites}
    for line in sam_lines:
        f = line.split("\t")
        if len(f) < 11 or f[0][0] == "@" or f[5] == "*":
            continue
        chrom, pos, flag = f[2], int(f[3]), int(f[1])
        if flag & (0x4 | 0x100 | 0x400 | 0x800):
            continue
        for k, (c, p) in sites.items():
            if chrom != c:
                continue
            if abs(pos - p) <= 300:
                n_reads[k] += 1
            for b, side in clip_boundaries(pos, f[5], min_clip):
                if abs(b - p) <= window:
                    res[k][side] += 1
    return res, n_reads


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples", nargs="+")
    ap.add_argument("--bamdir", default=os.path.expanduser("~/Temp/SexDeterminationAnalysis/sv_genotypes/region_bams"))
    ap.add_argument("--window", type=int, default=25, help="a clip counts if it begins within this many bp of the site (default 25)")
    a = ap.parse_args()
    print(f"{'sample':14s} {'site':4s} {'clipped reads at the site (left/right)':40s} reads within 300 bp")
    for s in a.samples:
        bam = os.path.join(a.bamdir, f"{s}_sv.bam")
        if not os.path.exists(bam):
            print(f"{s:14s} no BAM at {bam}")
            continue
        p = subprocess.run(["samtools", "view", bam], capture_output=True, text=True)
        if p.returncode:
            sys.exit(p.stderr)
        res, n = count_clips(p.stdout.splitlines(), window=a.window)
        for k in SITES:
            r = res[k]
            print(f"{s:14s} {k:4s} {r['left'] + r['right']:3d} clipped ({r['left']}/{r['right']}){'':24s} {n[k]}")


if __name__ == "__main__":
    main()
