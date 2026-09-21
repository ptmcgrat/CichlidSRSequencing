"""Write a VCF of X/Y-distinguishing SNVs suitable for allele-specific expression.

Keeps a candidate only if it is:

  1. an SNV -- indels make allelic read assignment unreliable in RNA-seq;
  2. biallelic -- ASE counters (GATK ASEReadCounter, phASER) count REF vs ALT
     and do not handle a third allele. XY_ sites, where X and Y each differ from
     the reference, are dropped by default for this reason and counted in the log;
  3. inside a predicted transcript (any exon of any isoform), since only
     transcribed sequence yields RNA-seq reads;
  4. heterozygous in the majority of Yellow Head fish carrying one copy of the
     inversion -- the XY individuals whose expression you would actually
     measure. A marker that is not het in those fish cannot separate the two
     alleles in them, whatever the assembly comparison says.

Each record carries which allele sits on the X and which on the Y, so ASE results
can be reported as X- versus Y-expressed rather than REF versus ALT. For an X_
variant the X carries ALT; for a Y_ variant the Y does -- the mapping flips, so
without this tag the direction of any imbalance is ambiguous.

    python3 buildASEsites.py
    python3 buildASEsites.py --min-het-frac 0.8 --min-called 20
"""

import argparse
import gzip
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

import pandas as pd

from helper_modules.file_manager import FileManager as FM
from helper_modules import pipeline_checks as pc
from helper_modules.pipeline_checks import PipelineError, log, warn

YH_CATS = ("YH_Brood1", "YH_Brood2", "OtherYHs")


def parse_args():
    p = argparse.ArgumentParser(description="ASE marker VCF from X/Y SNVs.")
    p.add_argument("--contig", default="NC_135176.1")
    p.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    p.add_argument("--variants", default="candidateQTNs_chr10_XY.tsv")
    p.add_argument("--source-dir", default="QTG_Candidates")
    p.add_argument("--min-het-frac", type=float, default=0.5,
                   help="Minimum fraction of called inv=1 YH samples that must be "
                        "heterozygous. Default 0.5 means a strict majority.")
    p.add_argument("--min-called", type=int, default=10,
                   help="Minimum inv=1 YH samples with a call. A majority of three "
                        "is not evidence (default 10).")
    p.add_argument("--require-hom-inv2", type=float, default=None,
                   help="Optionally also require this fraction of inv=2 YH samples "
                        "to be homozygous for the X allele, confirming the site "
                        "really separates X from Y. Off by default.")
    p.add_argument("--include-multiallelic", action="store_true",
                   help="Keep XY_ sites with two ALT alleles. Most ASE tools "
                        "cannot use them.")
    p.add_argument("--out", default=None)
    p.add_argument("--no-upload", action="store_true")
    return p.parse_args()


def open_maybe_gz(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def attr(s, key):
    m = re.search(rf"(?:^|;){key}=([^;]*)", s)
    return m.group(1) if m else None


def exon_index(gff, contig):
    """Exons with their transcript and gene, sorted for lookup."""
    gene_name, tx_gene = {}, {}
    exons = []
    with open_maybe_gz(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != contig:
                continue
            kind, a = f[2], f[8]
            if kind in ("gene", "pseudogene"):
                gid = attr(a, "ID")
                if gid:
                    gene_name[gid] = attr(a, "Name") or gid
            elif kind in ("mRNA", "transcript", "lnc_RNA", "lncRNA"):
                tid, parent = attr(a, "ID"), attr(a, "Parent")
                if tid and parent:
                    tx_gene[tid] = parent
            elif kind == "exon":
                parent = attr(a, "Parent")
                if parent:
                    exons.append((int(f[3]), int(f[4]), parent))
    exons.sort()
    return exons, gene_name, tx_gene


def main():
    args = parse_args()
    fm_obj = FM(genome_version=args.genome_version)
    web = fm_obj.localNikeshDir + "WebServer/"

    def fetch(path, what):
        if not os.path.exists(path):
            try:
                fm_obj.downloadData(path)
            except Exception as e:
                raise PipelineError(f"could not fetch {what} at {path}: {e}")
        pc.require_file(path, what)
        return path

    var = args.variants
    if not os.path.isabs(var):
        rel = var if "/" in var else args.source_dir.strip("/") + "/" + var
        var = fm_obj.localNikeshDir + rel
    fetch(var, "candidate table")
    gff = fetch(os.path.join(fm_obj.localGenomeDir, f"{args.contig}.gff3.gz"),
                "annotation")
    gt_bin = fetch(web + f"{args.contig}_genotypes.bin", "packed genotypes")
    gt_json = fetch(web + f"{args.contig}_genotypes_samples.json",
                    "genotype sample list")

    # --- genotypes ---------------------------------------------------------
    meta = json.load(open(gt_json))
    NV, NS = meta["n_variants"], meta["n_samples"]
    samples = meta["samples"]
    raw = open(gt_bin, "rb").read()
    if len(raw) != (NV * NS + 3) // 4:
        raise PipelineError("genotype matrix size does not match its sample list")

    def gt(vi, si):
        k = vi * NS + si
        return (raw[k >> 2] >> ((k & 3) * 2)) & 3

    yh1 = [i for i, s in enumerate(samples)
           if s.get("cat") in YH_CATS and s.get("inv") == 1]
    yh2 = [i for i, s in enumerate(samples)
           if s.get("cat") in YH_CATS and s.get("inv") == 2]
    by_cat = Counter(samples[i]["cat"] for i in yh1)
    log(f"YH samples with one inversion copy (inv=1): {len(yh1)} {dict(by_cat)}")
    log(f"YH samples with two copies (inv=2): {len(yh2)}")
    if not yh1:
        raise PipelineError(
            "no YH samples with Inversion10 = 1. Check that the sample list "
            "carries OtherYHs and the brood categories -- re-run "
            "buildVariantStats.py if it predates that change.")

    # --- candidates, in the same order the matrix was built ----------------
    dt = pd.read_csv(var, sep="\t")
    dt = dt[dt.Chromosome == args.contig].reset_index(drop=True)
    if len(dt) != NV or dt.Position.tolist() != meta["positions"]:
        raise PipelineError(
            "candidate table order does not match the genotype matrix. The "
            "matrix was built from a different version of the table -- rerun "
            "buildVariantStats.py against this one.")
    dt["cls"] = dt.Notes.str.extract(r"CLASS=([^;]+)")[0].fillna("")
    dt["GTph"] = dt.Notes.str.extract(r"GT=([^;]+)")[0].fillna("")
    dt["hap"] = dt.Name.astype(str).str.split("_").str[0]

    exons, gene_name, tx_gene = exon_index(gff, args.contig)
    starts = [e[0] for e in exons]
    log(f"{len(exons):,} exon records indexed")

    from bisect import bisect_right

    def transcripts_at(pos):
        """Transcripts whose exons cover this position."""
        hits = set()
        j = bisect_right(starts, pos) - 1
        # exons are sorted by start; walk back over any that could still cover pos
        while j >= 0 and pos - exons[j][0] < 250000:
            a, b, tid = exons[j]
            if a <= pos <= b:
                hits.add(tid)
            j -= 1
        return hits

    tally = Counter()
    kept = []
    for i, r in dt.iterrows():
        if r.cls != "SNV":
            tally["not an SNV"] += 1
            continue
        alts = str(r.Alt).split(",")
        if len(alts) > 1 and not args.include_multiallelic:
            tally["multiallelic (XY_)"] += 1
            continue
        txs = transcripts_at(int(r.Position))
        if not txs:
            tally["outside any transcript"] += 1
            continue

        het = called = 0
        for si in yh1:
            g = gt(i, si)
            if g == 3:
                continue
            called += 1
            if g == 1:
                het += 1
        if called < args.min_called:
            tally[f"fewer than {args.min_called} inv=1 YH called"] += 1
            continue
        frac = het / called
        if frac <= args.min_het_frac:
            tally[f"het in <= {args.min_het_frac:.0%} of inv=1 YH"] += 1
            continue

        # Which allele is on which haplotype. hap1 of the phased GT is the X,
        # hap2 the Y, matching the candidate table's encoding.
        ph = r.GTph.split("|") if "|" in r.GTph else ["", ""]
        alle = [str(r.Reference).upper()] + [a.upper() for a in alts]
        try:
            x_al = alle[int(ph[0])]
            y_al = alle[int(ph[1])]
        except (ValueError, IndexError):
            x_al = y_al = "."

        hom2 = called2 = 0
        if yh2:
            x_idx = int(ph[0]) if ph[0].isdigit() else None
            for si in yh2:
                g = gt(i, si)
                if g == 3:
                    continue
                called2 += 1
                # homozygous for the X allele: 0/0 if X is REF, 1/1 if X is ALT
                if (x_idx == 0 and g == 0) or (x_idx and x_idx > 0 and g == 2):
                    hom2 += 1
        if args.require_hom_inv2 is not None:
            if not called2 or hom2 / called2 < args.require_hom_inv2:
                tally["inv=2 YH not homozygous for the X allele"] += 1
                continue

        genes = sorted({gene_name.get(tx_gene.get(t, ""), "") for t in txs} - {""})
        kept.append({
            "pos": int(r.Position), "id": r.Name,
            "ref": str(r.Reference).upper(), "alt": ",".join(a.upper() for a in alts),
            "hap": r.hap, "gt": r.GTph, "x": x_al, "y": y_al,
            "genes": genes, "txs": sorted(txs),
            "het": het, "called": called, "frac": frac,
            "hom2": hom2, "called2": called2,
        })
        tally["KEPT"] += 1

    log("filter breakdown:")
    for k, v in sorted(tally.items(), key=lambda kv: (kv[0] != "KEPT", -kv[1])):
        log(f"    {k:42s} {v:7,}")

    if not kept:
        raise PipelineError("no sites passed every filter")

    genes_hit = Counter(g for k in kept for g in k["genes"])
    log(f"{len(kept):,} ASE marker SNVs across {len(genes_hit):,} genes")
    log(f"markers per gene: median "
        f"{sorted(genes_hit.values())[len(genes_hit)//2]}, "
        f"max {max(genes_hit.values())}")
    log(f"X allele is REF at {sum(1 for k in kept if k['x']==k['ref']):,} sites, "
        f"ALT at {sum(1 for k in kept if k['x']!=k['ref']):,}")

    # --- write the VCF -----------------------------------------------------
    out = args.out or (web + f"{args.contig}_ASE_markers.vcf")
    fai = fm_obj.localGenomeFile + ".fai"
    contig_line = f"##contig=<ID={args.contig}>"
    if os.path.exists(fai):
        for line in open(fai):
            n, L = line.split("\t")[:2]
            if n == args.contig:
                contig_line = f"##contig=<ID={n},length={L}>"
    with open(out, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##source=buildASEsites.py min_het_frac={args.min_het_frac} "
                 f"min_called={args.min_called}\n")
        fh.write(contig_line + "\n")
        for tag, num, typ, desc in [
            ("GENE", ".", "String", "Gene(s) whose transcripts contain this site"),
            ("TX", ".", "String", "Transcript(s) whose exons contain this site"),
            ("HAP", "1", "String", "Haplotype label from the assembly comparison (X or Y)"),
            ("XALLELE", "1", "String", "Allele carried on the X haplotype"),
            ("YALLELE", "1", "String", "Allele carried on the Y haplotype"),
            ("YH_HET", "1", "Integer", "inv=1 Yellow Head samples heterozygous here"),
            ("YH_CALLED", "1", "Integer", "inv=1 Yellow Head samples with a call"),
            ("YH_HETFRAC", "1", "Float", "YH_HET / YH_CALLED"),
            ("YH2_HOMX", "1", "Integer", "inv=2 Yellow Head samples homozygous for the X allele"),
            ("YH2_CALLED", "1", "Integer", "inv=2 Yellow Head samples with a call"),
        ]:
            fh.write(f'##INFO=<ID={tag},Number={num},Type={typ},'
                     f'Description="{desc}">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for k in sorted(kept, key=lambda k: k["pos"]):
            info = ";".join([
                f"GENE={','.join(k['genes']) or '.'}",
                f"TX={','.join(k['txs'])}",
                f"HAP={k['hap']}",
                f"XALLELE={k['x']}", f"YALLELE={k['y']}",
                f"YH_HET={k['het']}", f"YH_CALLED={k['called']}",
                f"YH_HETFRAC={k['frac']:.3f}",
                f"YH2_HOMX={k['hom2']}", f"YH2_CALLED={k['called2']}",
            ])
            fh.write("\t".join([args.contig, str(k["pos"]), k["id"], k["ref"],
                                k["alt"], ".", "PASS", info]) + "\n")

    subprocess.run(["bgzip", "-f", out], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", out + ".gz"], check=True)
    log(f"wrote {out}.gz ({os.path.getsize(out + '.gz')/1e3:.0f} KB)")

    if not args.no_upload:
        try:
            fm_obj.uploadData(out + ".gz")
            fm_obj.uploadData(out + ".gz.tbi")
            log("uploaded")
        except Exception as e:
            warn(f"upload failed: {e}")


if __name__ == "__main__":
    main()