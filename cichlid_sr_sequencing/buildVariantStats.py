"""Summarise the per-sample genotypes into per-variant statistics.

One pass over every sample VCF produces, for each of the 60,396 variants: how
many samples were called, the allele frequency, the genotype counts, and the
correlation between the variant's allele dosage and chr10 inversion dosage.

That correlation is the test the candidate table cannot make on its own. Names
like Y_INS_24921394 record which haplotype the variant sat on in the assembly
comparison; the population genotypes say whether it actually segregates with the
inverted haplotype. A Y-labelled insertion should correlate positively with
Inversion10, an X-labelled one negatively. Disagreement means the assembly-based
label and the population data tell different stories -- which is worth seeing,
especially for transposons, where a mobile element can be present on both
haplotypes or have moved since.

Output is per-variant, so the browser needs no genotype matrix: a few numbers per
variant instead of one per sample per variant, which is the difference between a
2 MB payload and a 200 MB one.

    python3 buildVariantStats.py --vcf-dir chr10_XY
"""

import argparse
import glob
import gzip
import json
import math
import os
import sys
from collections import Counter, defaultdict

import pandas as pd

from helper_modules.file_manager import FileManager as FM
from helper_modules import pipeline_checks as pc
from helper_modules.pipeline_checks import PipelineError, log, warn

GT_MAP = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1,
          "1/1": 2, "1|1": 2, "0/2": 1, "0|2": 1, "1/2": 2, "1|2": 2,
          "2/2": 2, "2|2": 2}


def parse_args():
    p = argparse.ArgumentParser(description="Per-variant genotype statistics.")
    p.add_argument("--vcf-dir", default="chr10_XY",
                   help="Folder of per-sample VCFs, relative to the Nikesh "
                        "directory (or an absolute path).")
    p.add_argument("--variants", default="candidateQTNs_chr10_XY.tsv")
    p.add_argument("--source-dir", default="QTG_Candidates")
    p.add_argument("--contig", default="NC_135176.1")
    p.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    p.add_argument("--min-dp", type=int, default=6)
    p.add_argument("--out", default=None)
    p.add_argument("--genotypes", action="store_true", default=True,
                   help="Also write a packed genotype matrix for the browser's "
                        "per-variant sample distribution (default on).")
    p.add_argument("--no-genotypes", dest="genotypes", action="store_false")
    p.add_argument("--no-upload", action="store_true")
    return p.parse_args()


def load_samples(fm_obj):
    """Sample metadata keyed by SampleID, including inversion state."""
    fm_obj.readSampleDatabase()
    sdt = fm_obj.sample_dt
    meta = {}
    for _, r in sdt.iterrows():
        meta[r.SampleID] = {
            "eco": r.Ecogroup if pd.notna(r.Ecogroup) else "",
            "sub": r.Subgroup if pd.notna(r.get("Subgroup")) else "",
            "cat": r.Category if pd.notna(r.get("Category")) else "",
            "sex": r.Sex if pd.notna(r.Sex) else "",
            "inv": int(r.Inversion10) if pd.notna(r.get("Inversion10")) else -1,
            "lab": int(r.LabReared) if pd.notna(r.get("LabReared")) else -1,
        }
    return meta


def verified_samples(vcf_dir):
    """Only load VCFs whose manifest says the run verified.

    Filenames are reused across runs, so a file from an aborted attempt looks
    exactly like a good one. The manifest is the record of whether it was
    checked.
    """
    ok = {}
    for path in sorted(glob.glob(os.path.join(vcf_dir, "*.manifest.json"))):
        try:
            m = json.load(open(path))
        except Exception:
            continue
        if m.get("status") == "ok" and m.get("sample_id"):
            ok[m["sample_id"]] = m
    out = []
    for vcf in sorted(glob.glob(os.path.join(vcf_dir, "*.vcf.gz"))):
        base = os.path.basename(vcf)
        sid = base.rsplit("_", 1)[0] if "_" in base else base.replace(".vcf.gz", "")
        # Filenames are <SampleID>_<run-name>.vcf.gz; match against manifests.
        match = next((s for s in ok if base.startswith(s + "_")), None)
        if match:
            out.append((match, vcf))
    return out


def main():
    args = parse_args()
    fm_obj = FM(genome_version=args.genome_version)

    vcf_dir = args.vcf_dir
    if not os.path.isabs(vcf_dir):
        vcf_dir = fm_obj.localNikeshDir + vcf_dir.strip("/") + "/"
    if not os.path.isdir(vcf_dir):
        try:
            fm_obj.downloadData(vcf_dir.rstrip("/"))
        except Exception as e:
            raise PipelineError(f"could not fetch {vcf_dir}: {e}")
    pairs = verified_samples(vcf_dir)
    if not pairs:
        raise PipelineError(f"no verified sample VCFs in {vcf_dir}")
    log(f"{len(pairs)} verified sample VCF(s)")

    meta_all = load_samples(fm_obj)
    meta = meta_all
    missing_meta = [s for s, _ in pairs if s not in meta]
    if missing_meta:
        warn(f"{len(missing_meta)} sample(s) have no database row: "
             f"{missing_meta[:5]}")
    inv = [meta.get(s, {}).get("inv", -1) for s, _ in pairs]
    n_typed = sum(1 for v in inv if v >= 0)
    log(f"{n_typed} sample(s) have an Inversion10 call "
        f"({Counter(v for v in inv if v >= 0)})")

    var = args.variants
    if not os.path.isabs(var):
        rel = var if "/" in var else args.source_dir.strip("/") + "/" + var
        var = fm_obj.localNikeshDir + rel
        if not os.path.exists(var):
            fm_obj.downloadData(var)
    dt = pd.read_csv(var, sep="\t")
    dt = dt[dt.Chromosome == args.contig]
    key_to_i = {int(p): i for i, p in enumerate(dt.Position)}
    names = dt.Name.tolist()
    positions = dt.Position.tolist()
    classes = dt.Notes.str.extract(r"CLASS=([^;]+)")[0].fillna("unknown").tolist()
    NV = len(dt)
    log(f"{NV:,} variants in the candidate table")

    # Accumulators, per variant
    n_called = [0] * NV
    n_het = [0] * NV
    n_hom = [0] * NV
    # Running sums for the dosage correlation against inversion state
    sx = [0.0] * NV      # sum of inversion dosage
    sy = [0.0] * NV      # sum of variant dosage
    sxx = [0.0] * NV
    syy = [0.0] * NV
    sxy = [0.0] * NV
    n_corr = [0] * NV
    # Carriers and called counts split by inversion state. This is the statistic
    # the scorecard actually uses -- a Y-restricted variant is present in
    # heterozygotes (inv=1) and absent from homozygotes (inv=2) -- and it is NOT
    # a linear function of inversion dosage, so the Pearson r above cannot
    # express it.
    carr = [[0, 0, 0] for _ in range(NV)]
    call_by_inv = [[0, 0, 0] for _ in range(NV)]

    # Packed genotypes, two bits each: 0=0/0, 1=0/1, 2=1/1, 3=no call.
    # 60,396 variants x 211 samples is 12.7M calls -- 3.2 MB packed, against
    # 64 MB as one byte per field. Small enough to ship inside the browser.
    NS = len(pairs)
    import array
    # 0xFF sets every 2-bit field to 3 (no call) in one go. Looping over 12.7M
    # fields to do the same thing takes minutes.
    gt_packed = array.array("B", b"\xff" * ((NV * NS + 3) // 4))

    def set_gt(vi, si, code):
        k = vi * NS + si
        byte, shift = k >> 2, (k & 3) * 2
        gt_packed[byte] = (gt_packed[byte] & ~(3 << shift)) | (code << shift)

    for si, (sid, path) in enumerate(pairs, 1):
        inv_s = meta.get(sid, {}).get("inv", -1)
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                i = key_to_i.get(int(f[1]))
                if i is None:
                    continue
                fmt, val = f[8].split(":"), f[9].split(":")
                d = dict(zip(fmt, val))
                g = GT_MAP.get(d.get("GT", "./."))
                if g is None:
                    continue
                dp = d.get("DP", "")
                if dp.isdigit() and int(dp) < args.min_dp:
                    continue
                n_called[i] += 1
                if args.genotypes:
                    set_gt(i, si - 1, g)
                if g == 1:
                    n_het[i] += 1
                elif g == 2:
                    n_hom[i] += 1
                if inv_s >= 0:
                    call_by_inv[i][inv_s] += 1
                    if g > 0:
                        carr[i][inv_s] += 1
                    n_corr[i] += 1
                    sx[i] += inv_s; sy[i] += g
                    sxx[i] += inv_s * inv_s; syy[i] += g * g
                    sxy[i] += inv_s * g
        if si % 25 == 0 or si == len(pairs):
            log(f"  {si}/{len(pairs)} samples read")

    rows = []
    for i in range(NV):
        n = n_corr[i]
        r = 0.0
        if n > 2:
            num = n * sxy[i] - sx[i] * sy[i]
            den = math.sqrt(max(0.0, n * sxx[i] - sx[i] ** 2)) * \
                  math.sqrt(max(0.0, n * syy[i] - sy[i] ** 2))
            r = num / den if den > 0 else 0.0
        called = n_called[i]
        af = ((n_het[i] + 2 * n_hom[i]) / (2 * called)) if called else None
        hap = str(names[i]).split("_")[0]

        # phi between "inversion heterozygote" and "carries the alt allele",
        # over inv=1 and inv=2 samples only. phi = 1 means the variant is
        # present in every heterozygote and absent from every homozygote --
        # the Y-restricted pattern. Deliberately separate from inv_r: the two
        # answer different questions and a variant can score high on one and
        # low on the other.
        a = carr[i][1]                      # inv=1 carriers
        b = call_by_inv[i][1] - a           # inv=1 non-carriers
        c = carr[i][2]                      # inv=2 carriers
        d = call_by_inv[i][2] - c           # inv=2 non-carriers
        den = math.sqrt(float((a + b) * (c + d) * (a + c) * (b + d)))
        phi_het = ((a * d - b * c) / den) if den > 0 else 0.0

        f1 = a / (a + b) if (a + b) else None
        f2 = c / (c + d) if (c + d) else None
        f0 = (carr[i][0] / call_by_inv[i][0]) if call_by_inv[i][0] else None

        rows.append({
            "name": names[i], "pos": positions[i], "hap": hap,
            "cls": classes[i],
            "n_called": called, "n_het": n_het[i], "n_hom": n_hom[i],
            "af": round(af, 4) if af is not None else "",
            "inv_r": round(r, 4), "n_inv_typed": n,
            "phi_het": round(phi_het, 4),
            "freq_inv0": round(f0, 4) if f0 is not None else "",
            "freq_inv1": round(f1, 4) if f1 is not None else "",
            "freq_inv2": round(f2, 4) if f2 is not None else "",
            "n_inv0": call_by_inv[i][0], "n_inv1": call_by_inv[i][1],
            "n_inv2": call_by_inv[i][2],
        })

    df = pd.DataFrame(rows)
    log(f"call rate: median {df.n_called.median():.0f}/{len(pairs)} samples")
    log(f"{(df.inv_r.abs() >= 0.5).sum():,} variants track inversion DOSAGE "
        f"at |r| >= 0.5")
    log(f"{(df.phi_het >= 0.7).sum():,} variants show the het-restricted pattern "
        f"at phi >= 0.7 (present in inv=1, absent from inv=2)")
    log(f"{(df.phi_het >= 0.9).sum():,} at phi >= 0.9")

    for cls in ("TE_insertion", "large_indel", "SNV", "small_indel"):
        sub = df[df.cls == cls]
        if not len(sub):
            continue
        log(f"  {cls:14s} n={len(sub):6,}  "
            f"|r|>=0.5: {(sub.inv_r.abs() >= 0.5).sum():5,}  "
            f"phi>=0.7: {(sub.phi_het >= 0.7).sum():5,}")

    # Where the haplotype label and the population data disagree, by class.
    for cls in ("TE_insertion",):
        sub = df[(df.cls == cls) & (df.phi_het >= 0.7)]
        if len(sub):
            log(f"  {cls} with phi>=0.7 by label: "
                f"{dict(Counter(sub.hap))}")

    if args.genotypes:
        gt_path = (fm_obj.localNikeshDir + "WebServer/"
                   + f"{args.contig}_genotypes.bin")
        os.makedirs(os.path.dirname(gt_path), exist_ok=True)
        with open(gt_path, "wb") as fh:
            gt_packed.tofile(fh)
        meta = {
            "n_variants": NV, "n_samples": NS,
            "positions": [int(p) for p in positions],
            "samples": [{"id": sid,
                         "eco": meta_all.get(sid, {}).get("eco", ""),
                         "sub": meta_all.get(sid, {}).get("sub", ""),
                         "cat": meta_all.get(sid, {}).get("cat", ""),
                         "sex": meta_all.get(sid, {}).get("sex", ""),
                         "inv": meta_all.get(sid, {}).get("inv", -1),
                         "lab": meta_all.get(sid, {}).get("lab", -1)}
                        for sid, _ in pairs],
        }
        with open(gt_path.replace(".bin", "_samples.json"), "w") as fh:
            json.dump(meta, fh)
        log(f"wrote {gt_path} ({os.path.getsize(gt_path)/1e6:.1f} MB packed, "
            f"{NV:,} x {NS})")

    out = args.out or (fm_obj.localNikeshDir + "WebServer/"
                       + f"{args.contig}_variant_stats.tsv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    log(f"wrote {out} ({os.path.getsize(out)/1e6:.1f} MB)")

    if not args.no_upload:
        try:
            fm_obj.uploadData(out)
            if args.genotypes:
                fm_obj.uploadData(gt_path)
                fm_obj.uploadData(gt_path.replace(".bin", "_samples.json"))
            log("uploaded")
        except Exception as e:
            warn(f"upload failed: {e}")


if __name__ == "__main__":
    main()