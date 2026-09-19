"""Predict protein consequences for X and Y haplotype variants.

Runs `bcftools csq` in haplotype-aware mode. That matters here: X and Y have
diverged enough that two substitutions can fall in the same codon, and their
combined effect is not the sum of their separate effects. csq reconstructs each
haplotype's protein and compares, which annotating variants one at a time cannot
do.

The candidate table is already in the right shape for it. Each variant carries a
phased genotype -- X_ variants are 1|0, Y_ are 0|1, XY_ are 1|2 with two
comma-separated ALT alleles -- so a single pseudo-sample whose two haplotypes
are the X and the Y reproduces exactly the comparison of interest.

Two conversions are needed first:

  1. csq accepts only Ensembl-style GFF3: ID=gene:<id>, ID=transcript:<id>,
     Parent=transcript:<id>, biotype=<x>. NCBI writes ID=gene-LOC…,
     Parent=rna-…, gene_biotype=<x>, so the file is rewritten here. Feeding csq
     the NCBI file fails with "Parent=transcript: not present".
  2. The TSV becomes a phased VCF with one sample.

    python3 buildConsequences.py
    python3 buildConsequences.py --variants candidateQTNs_chr10_XY.tsv
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

import pandas as pd

from helper_modules.file_manager import FileManager as FM
from helper_modules import pipeline_checks as pc
from helper_modules.pipeline_checks import PipelineError, log, warn


# Transcript-level feature types in NCBI GFF3 that should become
# "ID=transcript:" lines for csq.
TRANSCRIPT_TYPES = {
    "mRNA", "transcript", "lnc_RNA", "lncRNA", "tRNA", "rRNA", "ncRNA",
    "snRNA", "snoRNA", "primary_transcript", "miRNA", "guide_RNA",
    "antisense_RNA", "RNase_P_RNA", "SRP_RNA", "telomerase_RNA", "Y_RNA",
    "vault_RNA", "scRNA", "misc_RNA", "pseudogenic_transcript",
}
CHILD_TYPES = {"CDS", "exon", "three_prime_UTR", "five_prime_UTR"}


def parse_args():
    p = argparse.ArgumentParser(description="Predict X/Y protein consequences.")
    p.add_argument("--gff", default=None, help="Per-contig NCBI GFF3")
    p.add_argument("--variants", default="candidateQTNs_chr10_XY.tsv")
    p.add_argument("--source-dir", default="QTG_Candidates")
    p.add_argument("--contig", default="NC_135176.1")
    p.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    p.add_argument("--workdir", default=None)
    p.add_argument("--out", default=None,
                   help="Output TSV of per-variant consequences. Default: "
                        "<NikeshDir>/WebServer/<contig>_consequences.tsv")
    p.add_argument("--ncsq", type=int, default=64,
                   help="Max consequences recorded per site. bcftools defaults "
                        "to 16 and truncates in gene-dense regions where a "
                        "variant hits many transcripts (default 64).")
    p.add_argument("--no-upload", action="store_true")
    return p.parse_args()


def open_maybe_gz(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def attr(s, key):
    m = re.search(rf"(?:^|;){key}=([^;]*)", s)
    return m.group(1) if m else None


def ncbi_gff_to_ensembl(src, out_path, contig):
    """Rewrite an NCBI GFF3 into the dialect bcftools csq parses.

    csq looks for CDS, exon and UTR lines, finds their parent via
    "Parent=transcript:", the gene via the transcript's "Parent=gene:", and the
    biotype from "biotype=". NCBI uses different ID conventions and puts the
    biotype only on gene lines, so transcripts inherit it here.
    """
    gene_biotype, tx_gene = {}, {}
    with open_maybe_gz(src) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != contig:
                continue
            if f[2] in ("gene", "pseudogene"):
                gid = attr(f[8], "ID")
                if gid:
                    gene_biotype[gid] = attr(f[8], "gene_biotype") or "protein_coding"
            elif f[2] in TRANSCRIPT_TYPES:
                tid, parent = attr(f[8], "ID"), attr(f[8], "Parent")
                if tid and parent:
                    tx_gene[tid] = parent

    kept = Counter()
    with open_maybe_gz(src) as fh, open(out_path, "w") as out:
        out.write("##gff-version 3\n")
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != contig:
                continue
            kind, a = f[2], f[8]

            if kind in ("gene", "pseudogene"):
                gid = attr(a, "ID")
                if not gid:
                    continue
                name = attr(a, "Name") or gid
                bt = gene_biotype.get(gid, "protein_coding")
                f[8] = f"ID=gene:{gid};biotype={bt};Name={name}"
            elif kind in TRANSCRIPT_TYPES:
                tid, parent = attr(a, "ID"), attr(a, "Parent")
                if not tid or not parent:
                    continue
                bt = gene_biotype.get(parent, "protein_coding")
                f[8] = f"ID=transcript:{tid};Parent=gene:{parent};biotype={bt}"
                f[2] = "transcript"
            elif kind in CHILD_TYPES:
                parent = attr(a, "Parent")
                if not parent:
                    continue
                f[8] = f"Parent=transcript:{parent}"
            else:
                continue

            kept[kind] += 1
            out.write("\t".join(f) + "\n")

    return kept


def tsv_to_phased_vcf(variants_tsv, contig, out_vcf, ref_fai):
    """Write the candidate table as a phased VCF with one pseudo-sample.

    Sample "XY" carries the phased genotype already present in the table:
    haplotype 1 is the X, haplotype 2 is the Y.
    """
    dt = pd.read_csv(variants_tsv, sep="\t")
    dt = dt[dt.Chromosome == contig].copy()
    dt["GTph"] = dt.Notes.str.extract(r"GT=([^;]+)")[0]
    missing = dt.GTph.isna().sum()
    if missing:
        warn(f"{missing} variants have no GT= in Notes; they are dropped")
        dt = dt[dt.GTph.notna()]

    contigs = []
    with open(ref_fai) as fh:
        for line in fh:
            name, length = line.split("\t")[:2]
            if name == contig:
                contigs.append(f"##contig=<ID={name},length={length}>")

    dt = dt.sort_values("Position")
    with open(out_vcf, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        for c in contigs:
            out.write(c + "\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                   "FILTER", "INFO", "FORMAT", "XY"]) + "\n")
        for _, r in dt.iterrows():
            out.write("\t".join([
                r.Chromosome, str(r.Position), str(r.Name),
                str(r.Reference).upper(), str(r.Alt).upper(),
                ".", "PASS", ".", "GT", str(r.GTph),
            ]) + "\n")
    return len(dt)


def run_csq(vcf, ref_fa, gff, out_vcf, ncsq=64):
    cmd = ["bcftools", "csq", "-f", ref_fa, "-g", gff,
           "-p", "a", "-l", "--ncsq", str(ncsq),
           "-O", "v", "-o", out_vcf, vcf]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise PipelineError(
            f"bcftools csq failed (exit {r.returncode}):\n{r.stderr[-3000:]}")
    for line in r.stderr.splitlines():
        if line.strip():
            log("csq: " + line.strip())
    return out_vcf


# Consequences that alter the protein. Everything else (synonymous, intronic,
# UTR, intergenic) leaves the coding sequence intact.
ALTERING = {
    "missense", "stop_gained", "stop_lost", "start_lost", "frameshift",
    "inframe_insertion", "inframe_deletion", "inframe_altering",
    "splice_acceptor", "splice_donor", "coding_sequence",
}


def parse_csq(out_vcf):
    """Pull BCSQ annotations out and flatten to one row per variant."""
    rows = []
    with open(out_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            m = re.search(r"BCSQ=([^;\t]+)", f[7])
            if not m:
                continue
            for entry in m.group(1).split(","):
                parts = entry.lstrip("*@").split("|")
                csq = parts[0] if parts else ""
                gene = parts[1] if len(parts) > 1 else ""
                tx = parts[2] if len(parts) > 2 else ""
                aa = parts[5] if len(parts) > 5 else ""
                rows.append({
                    "name": f[2], "pos": int(f[1]),
                    "consequence": csq, "gene": gene, "transcript": tx,
                    "aa_change": aa,
                    "altering": any(k in csq for k in ALTERING),
                })
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    pc.require_tools(("bcftools",))
    fm_obj = FM(genome_version=args.genome_version)

    gff = args.gff
    if gff is None:
        gdir = getattr(fm_obj, "localGenomeDir", None)
        gff = os.path.join(gdir, f"{args.contig}.gff3.gz")
    if not os.path.exists(gff):
        fm_obj.downloadData(gff)
    pc.require_file(gff, "annotation GFF3")

    var = args.variants
    if not os.path.isabs(var):
        rel = var if "/" in var else args.source_dir.strip("/") + "/" + var
        var = fm_obj.localNikeshDir + rel
        if not os.path.exists(var):
            fm_obj.downloadData(var)
    pc.require_file(var, "candidate table")

    ref = fm_obj.localGenomeFile
    if not os.path.exists(ref):
        fm_obj.downloadData(ref)
    if not os.path.exists(ref + ".fai"):
        subprocess.run(["samtools", "faidx", ref], check=True)

    work = args.workdir or (fm_obj.localNikeshDir + "WebServer/csq_work/")
    os.makedirs(work, exist_ok=True)

    ens_gff = os.path.join(work, f"{args.contig}.ensembl.gff3")
    kept = ncbi_gff_to_ensembl(gff, ens_gff, args.contig)
    log(f"converted GFF3 for csq: {dict(kept.most_common(6))}")
    if not kept.get("CDS"):
        raise PipelineError("no CDS features survived conversion; "
                            "protein consequences cannot be predicted")

    in_vcf = os.path.join(work, "xy_haplotypes.vcf")
    n = tsv_to_phased_vcf(var, args.contig, in_vcf, ref + ".fai")
    log(f"wrote {n:,} phased variants as sample 'XY' (hap1=X, hap2=Y)")

    csq_vcf = os.path.join(work, "xy_haplotypes.csq.vcf")
    run_csq(in_vcf, ref, ens_gff, csq_vcf, ncsq=args.ncsq)

    df = parse_csq(csq_vcf)
    if df.empty:
        warn("csq produced no BCSQ annotations -- check that the converted GFF3 "
             "contains transcripts overlapping the variants")
    else:
        log(f"{len(df):,} consequence records over "
            f"{df.name.nunique():,} variants")
        log("top consequences: "
            f"{dict(Counter(df.consequence).most_common(8))}")
        alt = df[df.altering]
        log(f"{alt.name.nunique():,} variants are protein-altering, "
            f"in {alt.gene.nunique():,} genes")

    out = args.out or (fm_obj.localNikeshDir + "WebServer/"
                       + f"{args.contig}_consequences.tsv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    log(f"wrote {out}")

    if not args.no_upload:
        try:
            fm_obj.uploadData(out)
            log("uploaded")
        except Exception as e:
            warn(f"upload failed: {e}")


if __name__ == "__main__":
    main()