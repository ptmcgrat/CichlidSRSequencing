"""Aggregate candidate variants by gene, for the gene-centric browser.

Produces one record per gene in the region of interest:

  - the gene's own body, and its DISPLAY INTERVAL, which runs from the previous
    gene's end to the next gene's start. Intergenic variants therefore appear
    under both flanking genes, which is intended: a variant between two genes is
    potentially relevant to either, and no variant falls through a gap.
  - variant counts split by class (SNV / small indel / TE insertion / large
    indel) and by haplotype (X / Y / XY), both within the gene body and across
    the whole display interval.
  - ASE markers: variants landing inside a predicted mRNA. Every variant in this
    set distinguishes the X haplotype from the Y -- X_ variants are 1|0, Y_ are
    0|1, XY_ are 1|2, so all three carry different alleles on the two
    haplotypes. Any of them inside a transcript can in principle assign an
    RNA-seq read to a parental allele. SNVs are counted separately because
    indels make allelic read assignment considerably harder.

Two further annotations are left as placeholders, since each needs an input this
script does not have: protein consequences (a `bcftools csq` pass over the CDS
features) and TE-haplotype correlation (per-sample genotypes).

    python3 buildGeneTable.py --gff <contig>.gff3.gz --variants candidates.tsv \\
        --region 11816878-29898132 --out gene_table.json
"""

import argparse
import gzip
import json
import os
import re
import sys
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict

import pandas as pd

from helper_modules.file_manager import FileManager as FM
from helper_modules import pipeline_checks as pc
from helper_modules.pipeline_checks import PipelineError, log, warn


def parse_args():
    p = argparse.ArgumentParser(description="Aggregate candidate variants by gene.")
    p.add_argument("--gff", default=None,
                   help="Per-contig GFF3. Default: <GenomeDir>/<contig>.gff3.gz "
                        "from prepareAnnotation.py, fetched from cloud storage "
                        "if not already local.")
    p.add_argument("--variants", default="candidateQTNs_chr10_XY.tsv",
                   help="Candidate TSV. A bare filename is taken as living in "
                        "--source-dir and downloaded; an absolute path is used "
                        "as-is.")
    p.add_argument("--source-dir", default="QTG_Candidates",
                   help="Cloud folder holding the candidate TSVs.")
    p.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    p.add_argument("--no-download", action="store_true",
                   help="Use local copies only.")
    p.add_argument("--no-upload", action="store_true",
                   help="Do not push the result back to cloud storage.")
    p.add_argument("--consequences", default=None,
                   help="Per-variant consequence TSV from buildConsequences.py. "
                        "Default: <NikeshDir>/WebServer/<contig>_consequences.tsv")
    p.add_argument("--variant-stats", default=None,
                   help="Per-variant statistics TSV from buildVariantStats.py. "
                        "Default: <NikeshDir>/WebServer/<contig>_variant_stats.tsv")
    p.add_argument("--phi-cut", type=float, default=0.7,
                   help="phi_het at or above which a variant counts as "
                        "Y-restricted (default 0.7).")
    p.add_argument("--contig", default="NC_135176.1")
    p.add_argument("--region", default=None,
                   help="start-end to restrict the OUTPUT to. Genes just outside "
                        "are still read, so the boundary genes get correct "
                        "display intervals. Default: the variants' own span.")
    p.add_argument("--biotypes", nargs="+",
                   default=["protein_coding", "lncRNA"],
                   help="Gene biotypes to include (default protein_coding and "
                        "lncRNA; tRNA/snRNA/snoRNA are excluded as they are short, "
                        "numerous, and rarely the object of interest here).")
    p.add_argument("--all-biotypes", action="store_true",
                   help="Include every biotype regardless of --biotypes.")
    p.add_argument("--out", default=None,
                   help="Output path. Default: <NikeshDir>/WebServer/"
                        "<contig>_gene_table.json")
    return p.parse_args()


def open_maybe_gz(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def attr(s, key):
    m = re.search(rf"(?:^|;){key}=([^;]*)", s)
    return m.group(1) if m else None


def read_annotation(path, contig):
    """Collect genes with their transcripts' exon spans.

    Exons are gathered per gene as a merged set of intervals across every
    transcript, so a variant in any isoform's exon counts as ASE-testable.
    """
    genes = {}
    rna_to_gene = {}
    tx_strand = {}
    exons = defaultdict(list)
    cds = defaultdict(list)

    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != contig:
                continue
            kind, start, end, strand, a = f[2], int(f[3]), int(f[4]), f[6], f[8]

            if kind == "gene" or kind == "pseudogene":
                gid = attr(a, "ID")
                if not gid:
                    continue
                genes[gid] = {
                    "id": gid,
                    "name": attr(a, "Name") or attr(a, "gene") or gid,
                    "biotype": attr(a, "gene_biotype") or "",
                    "start": start, "end": end, "strand": strand,
                    "exons": [], "transcripts": [],
                }
            elif kind in ("mRNA", "transcript", "lncRNA", "tRNA", "ncRNA",
                          "snoRNA", "snRNA", "rRNA", "primary_transcript"):
                rid, parent = attr(a, "ID"), attr(a, "Parent")
                if rid and parent:
                    rna_to_gene[rid] = parent
                    tx_strand[rid] = strand
            elif kind == "exon":
                parent = attr(a, "Parent")
                if parent:
                    exons[parent].append((start, end))
            elif kind == "CDS":
                parent = attr(a, "Parent")
                if parent:
                    cds[parent].append((start, end))

    # Per-transcript structure, kept separate from the merged per-gene exon set.
    # The browser draws transcripts individually -- isoforms differ, and which
    # exon a variant lands in depends on which transcript you are looking at.
    for rid in set(list(exons) + list(cds)):
        gid = rna_to_gene.get(rid)
        if not gid or gid not in genes:
            continue
        ex = sorted(exons.get(rid, []))
        cd = sorted(cds.get(rid, []))
        # Protein length from total coding bases: CDS/3 minus the stop codon.
        cds_bp = sum(b - a + 1 for a, b in cd)
        aa = max(0, cds_bp // 3 - 1) if cds_bp else 0
        genes[gid]["transcripts"].append({
            "id": rid, "strand": tx_strand.get(rid, genes[gid]["strand"]),
            "exons": [[a, b] for a, b in ex],
            "cds": [[a, b] for a, b in cd],
            "aa": aa,
            "start": min([a for a, _ in ex] or [genes[gid]["start"]]),
            "end": max([b for _, b in ex] or [genes[gid]["end"]]),
        })
        genes[gid]["exons"].extend(ex)

    for g in genes.values():
        merged = []
        for a_, b_ in sorted(g["exons"]):
            if merged and a_ <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], b_)
            else:
                merged.append([a_, b_])
        g["exons"] = merged
        g["transcripts"].sort(key=lambda t: (-t["aa"], t["id"]))
        # The longest-coding isoform stands in for the gene's protein.
        g["aa"] = g["transcripts"][0]["aa"] if g["transcripts"] else 0
        g["gene_bp"] = g["end"] - g["start"] + 1
        g["mrna_bp"] = sum(b - a + 1 for a, b in g["exons"])
        g["n_tx"] = len(g["transcripts"])
    return list(genes.values())


def read_variants(path, contig):
    dt = pd.read_csv(path, sep="\t")
    dt = dt[dt.Chromosome == contig].copy()
    dt["cls"] = dt.Notes.str.extract(r"CLASS=([^;]+)")[0].fillna("unknown")
    dt["hap"] = dt.Name.astype(str).str.split("_").str[0]
    dt["reflen"] = dt.Reference.astype(str).str.len()
    dt["altlen"] = dt.Alt.astype(str).str.len()
    dt["big"] = (dt.reflen > 10) | (dt.altlen > 10)
    dt = dt.sort_values("Position").reset_index(drop=True)
    return dt


def assign(genes, dt):
    """Build per-gene records with variant counts over body and display interval."""
    genes = sorted(genes, key=lambda g: (g["start"], g["end"]))
    positions = dt.Position.to_numpy()
    cls = dt.cls.to_numpy()
    hap = dt.hap.to_numpy()
    names = dt.Name.to_numpy()

    out = []
    for i, g in enumerate(genes):
        # Display interval: previous gene's end to next gene's start. Overlapping
        # neighbours would otherwise invert the interval, so clamp to the body.
        prev_end = genes[i - 1]["end"] if i > 0 else g["start"] - 10000
        next_start = genes[i + 1]["start"] if i + 1 < len(genes) else g["end"] + 10000
        disp_start = min(prev_end + 1, g["start"])
        disp_end = max(next_start - 1, g["end"])

        lo = bisect_left(positions, disp_start)
        hi = bisect_right(positions, disp_end)
        idx = range(lo, hi)

        body_lo = bisect_left(positions, g["start"])
        body_hi = bisect_right(positions, g["end"])

        exon_hits, exon_snv = [], 0
        for j in range(body_lo, body_hi):
            p = positions[j]
            if any(a <= p <= b for a, b in g["exons"]):
                exon_hits.append(names[j])
                if cls[j] == "SNV":
                    exon_snv += 1

        rec = {
            "id": g["id"], "name": g["name"], "biotype": g["biotype"],
            "strand": g["strand"], "start": g["start"], "end": g["end"],
            "disp_start": disp_start, "disp_end": disp_end,
            "n_exons": len(g["exons"]),
            "exons": g["exons"],
            "transcripts": g.get("transcripts", []),
            "aa": g.get("aa", 0),
            "n_tx": g.get("n_tx", 0),
            "gene_bp": g.get("gene_bp", g["end"] - g["start"] + 1),
            "mrna_bp": g.get("mrna_bp", 0),
            "n_body": int(body_hi - body_lo),
            "n_interval": int(hi - lo),
            "by_class": dict(Counter(cls[lo:hi])),
            "by_hap": dict(Counter(hap[lo:hi])),
            "body_by_class": dict(Counter(cls[body_lo:body_hi])),
            # ASE: a variant inside a predicted transcript can assign an RNA-seq
            # read to the X or the Y allele. SNVs counted separately because
            # indels make allelic assignment harder in short-read data.
            "ase_markers": len(exon_hits),
            "ase_snv_markers": exon_snv,
            "ase_examples": exon_hits[:5],
            # Filled in later, from bcftools csq and from per-sample genotypes.
            "protein_altering": None,
            "te_haplotype_correlated": None,
        }
        out.append(rec)
    return out


def fetch_inputs(fm_obj, args):
    """Resolve the GFF and candidate table, downloading whichever are missing.

    Both live in cloud storage already -- the GFF subset put there by
    prepareAnnotation.py, the candidate table beside the other TSVs -- so there
    is no reason to make anyone locate them by hand.
    """
    gff = args.gff
    if gff is None:
        gdir = getattr(fm_obj, "localGenomeDir", None)
        if not gdir:
            raise PipelineError("FileManager has no localGenomeDir; pass --gff")
        gff = os.path.join(gdir, f"{args.contig}.gff3.gz")
    if not os.path.exists(gff) and not args.no_download:
        try:
            fm_obj.downloadData(gff)
        except FileNotFoundError:
            raise PipelineError(
                f"no annotation at {gff}. Run prepareAnnotation.py --contig "
                f"{args.contig} first, or pass --gff.")
    pc.require_file(gff, "annotation GFF3")
    log(f"annotation: {gff}")

    var = args.variants
    if os.path.isabs(var):
        pc.require_file(var, "candidate table")
    else:
        rel = var if "/" in var else args.source_dir.strip("/") + "/" + var
        var = fm_obj.localNikeshDir + rel
        if not os.path.exists(var) and not args.no_download:
            try:
                fm_obj.downloadData(var)
            except FileNotFoundError:
                raise PipelineError(f"could not download {var}")
        pc.require_file(var, "candidate table")
    log(f"variants: {var}")
    return gff, var


def optional_table(fm_obj, explicit, default_name, what):
    """Load an optional annotation TSV, fetching it from cloud storage if needed."""
    path = explicit or (fm_obj.localNikeshDir + "WebServer/" + default_name)
    if not os.path.exists(path):
        try:
            fm_obj.downloadData(path)
        except Exception:
            pass
    if not os.path.exists(path):
        warn(f"no {what} at {path}; those fields stay empty")
        return None
    log(f"{what}: {path}")
    return pd.read_csv(path, sep="\t")


def enrich(recs, dt, csq, stats, phi_cut):
    """Attach consequence and genotype annotations to each gene record.

    Both are joined on variant, then aggregated over the gene's DISPLAY interval
    -- so an intergenic variant contributes to both flanking genes, consistent
    with how the intervals are drawn.
    """
    pos_of = dict(zip(dt.Name.astype(str), dt.Position))
    cls_of = dict(zip(dt.Name.astype(str),
                      dt.Notes.str.extract(r"CLASS=([^;]+)")[0].fillna("unknown")))

    altering_pos, altering_by_gene = set(), Counter()
    if csq is not None and len(csq):
        alt = csq[csq.altering == True] if "altering" in csq.columns else csq.iloc[0:0]
        for _, r in alt.iterrows():
            p = pos_of.get(str(r["name"]))
            if p is not None:
                altering_pos.add(int(p))
            if isinstance(r.get("gene"), str) and r["gene"]:
                altering_by_gene[r["gene"]] += 1

    phi_pos, te_phi_pos, best_phi = {}, set(), {}
    if stats is not None and len(stats):
        for _, r in stats.iterrows():
            p = int(r["pos"])
            phi = float(r.get("phi_het") or 0)
            best_phi[p] = max(best_phi.get(p, -1.0), phi)
            if phi >= phi_cut:
                phi_pos[p] = phi
                if str(r.get("cls")) == "TE_insertion":
                    te_phi_pos.add(p)

    for rec in recs:
        lo, hi = rec["disp_start"], rec["disp_end"]
        blo, bhi = rec["start"], rec["end"]
        in_iv = [p for p in phi_pos if lo <= p <= hi]
        rec["n_y_restricted"] = len(in_iv)
        rec["n_te_y_restricted"] = sum(1 for p in in_iv if p in te_phi_pos)
        body_phis = [v for p, v in best_phi.items() if blo <= p <= bhi]
        iv_phis = [v for p, v in best_phi.items() if lo <= p <= hi]
        rec["max_phi_body"] = round(max(body_phis), 3) if body_phis else None
        rec["max_phi_interval"] = round(max(iv_phis), 3) if iv_phis else None
        # Consequences are counted on the gene BODY: a coding change only makes
        # sense inside the transcript it alters, unlike the interval-wide counts.
        rec["protein_altering"] = (altering_by_gene.get(rec["name"], 0)
                                   if altering_by_gene else None)
        rec["te_haplotype_correlated"] = rec["n_te_y_restricted"]
    return recs


def main():
    args = parse_args()

    fm_obj = FM(genome_version=args.genome_version)
    gff, variants = fetch_inputs(fm_obj, args)

    out_path = args.out
    if out_path is None:
        web = fm_obj.localNikeshDir + "WebServer/"
        os.makedirs(web, exist_ok=True)
        out_path = web + f"{args.contig}_gene_table.json"

    genes = read_annotation(gff, args.contig)
    print(f"{len(genes)} gene records on {args.contig}")

    dt = read_variants(variants, args.contig)
    print(f"{len(dt):,} variants, {dt.Position.min():,}-{dt.Position.max():,}")

    if args.region:
        lo, hi = (int(x) for x in args.region.split("-"))
    else:
        lo, hi = int(dt.Position.min()), int(dt.Position.max())
    print(f"region: {lo:,}-{hi:,}")

    if not args.all_biotypes:
        before = len(genes)
        keep = set(args.biotypes)
        # Filter AFTER interval assignment would shift neighbours, so filter here
        # and accept that excluded biotypes do not break up intervals.
        genes = [g for g in genes if g["biotype"] in keep]
        print(f"biotype filter {sorted(keep)}: {len(genes)} of {before} genes")

    recs = assign(genes, dt)
    in_region = [r for r in recs if r["end"] >= lo and r["start"] <= hi]
    print(f"{len(in_region)} genes in region")

    csq = optional_table(fm_obj, args.consequences,
                         f"{args.contig}_consequences.tsv", "consequences")
    stats = optional_table(fm_obj, args.variant_stats,
                           f"{args.contig}_variant_stats.tsv", "variant statistics")
    in_region = enrich(in_region, dt, csq, stats, args.phi_cut)

    if csq is not None:
        n = sum(1 for r in in_region if (r["protein_altering"] or 0) > 0)
        print(f"{n} genes have >=1 protein-altering variant")
    if stats is not None:
        n = sum(1 for r in in_region if r["n_y_restricted"] > 0)
        t = sum(1 for r in in_region if r["n_te_y_restricted"] > 0)
        print(f"{n} genes have >=1 Y-restricted variant "
              f"(phi >= {args.phi_cut}); {t} have a Y-restricted TE insertion")
        three = sum(1 for r in in_region
                    if (r["protein_altering"] or 0) > 0
                    and r["ase_snv_markers"] > 0
                    and r["n_te_y_restricted"] > 0)
        print(f"{three} genes have all three: coding change, an ASE marker, "
              f"and a Y-restricted TE")

    with_ase = [r for r in in_region if r["ase_markers"]]
    print(f"{len(with_ase)} genes have >=1 ASE marker inside a transcript")
    print(f"{sum(1 for r in in_region if r['ase_snv_markers'])} have an exonic SNV marker")
    tot = sum(r["ase_markers"] for r in in_region)
    print(f"{tot:,} exonic variants total")

    sizes = sorted(r["n_interval"] for r in in_region)
    if sizes:
        print(f"variants per display interval: min {sizes[0]}, "
              f"median {sizes[len(sizes)//2]}, max {sizes[-1]}")

    with open(out_path, "w") as fh:
        json.dump({"contig": args.contig, "region": [lo, hi],
                   "genes": in_region}, fh)
    print(f"wrote {out_path} ({os.path.getsize(out_path)/1e6:.2f} MB)")

    if not args.no_upload:
        try:
            fm_obj.uploadData(out_path)
            print("uploaded to cloud storage")
        except Exception as e:
            warn(f"upload failed: {e}")


if __name__ == "__main__":
    main()