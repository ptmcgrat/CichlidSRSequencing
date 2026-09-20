"""Build the gene-centric chr10 X/Y browser: one self-contained HTML file.

Combines everything the earlier steps produced:

    prepareAnnotation.py  -> per-chromosome GFF3          (gene models)
    buildGeneTable.py     -> gene_table.json              (per-gene aggregates)
    buildConsequences.py  -> consequences.tsv             (protein effects)
    buildVariantStats.py  -> variant_stats.tsv            (genotype statistics)

The browser has two levels. A landing table of every gene in the region, sortable
and filterable on the three annotations, and a per-gene view showing that gene's
display interval -- previous gene's end to next gene's start -- with the gene
models, a variant track, and the variants in that interval.

Per-variant statistics rather than a genotype matrix keep this to a few MB. 211
samples x 60,396 variants would be 200 MB of genotypes; the summary statistics
that answer the questions being asked are closer to 4 MB.

    python3 buildXYBrowser.py
"""

import argparse
import gzip
import json
import os
import re
import sys
from collections import Counter, defaultdict

import pandas as pd

from helper_modules.file_manager import FileManager as FM
from helper_modules import pipeline_checks as pc
from helper_modules.pipeline_checks import PipelineError, log, warn

TEMPLATE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "helper_modules", "xy_browser_template.html")

# Ranked worst-first: a variant hitting several transcripts is summarised by its
# most severe consequence, which is what the gene table counts.
SEVERITY = ["stop_gained", "frameshift", "stop_lost", "start_lost",
            "splice_acceptor", "splice_donor", "inframe_deletion",
            "inframe_insertion", "inframe_altering", "missense",
            "splice_region", "synonymous", "non_coding", "intron",
            "5_prime_utr", "3_prime_utr", "intergenic"]
SEV_RANK = {c: i for i, c in enumerate(SEVERITY)}


def parse_args():
    p = argparse.ArgumentParser(description="Build the chr10 X/Y gene browser.")
    p.add_argument("--contig", default="NC_135176.1")
    p.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    p.add_argument("--variants", default="candidateQTNs_chr10_XY.tsv")
    p.add_argument("--source-dir", default="QTG_Candidates")
    p.add_argument("--out-name", default="xy_browser.html")
    p.add_argument("--phi-cut", type=float, default=0.7)
    p.add_argument("--no-upload", action="store_true")
    return p.parse_args()


def need(fm_obj, path, what):
    if not os.path.exists(path):
        try:
            fm_obj.downloadData(path)
        except Exception as e:
            raise PipelineError(f"could not fetch {what} at {path}: {e}")
    pc.require_file(path, what)
    log(f"{what}: {path}")
    return path


def worst_consequence(csq):
    """One consequence per variant: the most severe across all transcripts."""
    best = {}
    if csq is None or not len(csq):
        return best
    for name, c in zip(csq["name"].astype(str), csq["consequence"].astype(str)):
        # csq can emit compound terms such as "missense&splice_region"
        terms = re.split(r"[&+]", c)
        rank = min((SEV_RANK.get(t, 99) for t in terms), default=99)
        if name not in best or rank < best[name][0]:
            best[name] = (rank, c)
    return {k: v[1] for k, v in best.items()}


def main():
    args = parse_args()
    fm_obj = FM(genome_version=args.genome_version)
    web = fm_obj.localNikeshDir + "WebServer/"
    os.makedirs(web, exist_ok=True)

    gene_json = need(fm_obj, web + f"{args.contig}_gene_table.json", "gene table")
    csq_tsv = need(fm_obj, web + f"{args.contig}_consequences.tsv", "consequences")
    stats_tsv = need(fm_obj, web + f"{args.contig}_variant_stats.tsv",
                     "variant statistics")

    var = args.variants
    if not os.path.isabs(var):
        rel = var if "/" in var else args.source_dir.strip("/") + "/" + var
        var = fm_obj.localNikeshDir + rel
    need(fm_obj, var, "candidate table")

    genes = json.load(open(gene_json))
    log(f"{len(genes['genes'])} genes")
    if genes["genes"] and "exons" not in genes["genes"][0]:
        warn("gene table has no exon coordinates -- re-run buildGeneTable.py so "
             "the detail view can draw gene models")

    dt = pd.read_csv(var, sep="\t")
    dt = dt[dt.Chromosome == args.contig].copy()
    dt["cls"] = dt.Notes.str.extract(r"CLASS=([^;]+)")[0].fillna("unknown")
    dt["hap"] = dt.Name.astype(str).str.split("_").str[0]
    dt = dt.sort_values("Position").reset_index(drop=True)

    csq = pd.read_csv(csq_tsv, sep="\t")
    worst = worst_consequence(csq)
    altering_names = set()
    if "altering" in csq.columns:
        altering_names = set(csq.loc[csq.altering == True, "name"].astype(str))
    log(f"{len(worst):,} variants have a consequence call, "
        f"{len(altering_names):,} protein-altering")

    st = pd.read_csv(stats_tsv, sep="\t").set_index("name")
    log(f"{len(st):,} variants have genotype statistics")

    # Which variants sit inside a transcript -- the ASE marker test.
    exonic = set()
    for g in genes["genes"]:
        for a, b in g.get("exons", []):
            lo = dt.Position.searchsorted(a, "left")
            hi = dt.Position.searchsorted(b, "right")
            exonic.update(dt.Name.iloc[lo:hi].astype(str))
    log(f"{len(exonic):,} variants fall inside a predicted transcript")

    cls_levels = sorted(dt.cls.unique())
    csq_levels = sorted({c for c in worst.values()})
    cls_idx = {c: i for i, c in enumerate(cls_levels)}
    csq_idx = {c: i for i, c in enumerate(csq_levels)}

    # Parallel arrays: markedly smaller than an array of objects once repeated
    # keys are gone, and the browser indexes them directly.
    V = {"pos": [], "name": [], "cls": [], "hap": [], "phi": [], "r": [],
         "af": [], "nc": [], "csq": [], "exon": [], "alt": []}
    for row in dt.itertuples():
        n = str(row.Name)
        s = st.loc[n] if n in st.index else None
        V["pos"].append(int(row.Position))
        V["name"].append(n)
        V["cls"].append(cls_idx[row.cls])
        V["hap"].append({"X": 0, "Y": 1, "XY": 2}.get(row.hap, 3))
        V["phi"].append(round(float(s["phi_het"]), 3) if s is not None else None)
        V["r"].append(round(float(s["inv_r"]), 3) if s is not None else None)
        V["af"].append(round(float(s["af"]), 3)
                       if s is not None and str(s["af"]) != "nan" else None)
        V["nc"].append(int(s["n_called"]) if s is not None else 0)
        V["csq"].append(csq_idx.get(worst.get(n), -1))
        V["exon"].append(1 if n in exonic else 0)
        V["alt"].append(1 if n in altering_names else 0)

    payload = {
        "contig": args.contig,
        "region": genes["region"],
        "phi_cut": args.phi_cut,
        "cls_levels": cls_levels,
        "csq_levels": csq_levels,
        "genes": genes["genes"],
        "variants": V,
        "n_samples": int(st["n_called"].max()) if len(st) else 0,
    }

    pc.require_file(TEMPLATE, "browser template")
    html = open(TEMPLATE).read().replace(
        "__DATA__", json.dumps(payload, separators=(",", ":")))

    out = web + args.out_name
    with open(out, "w") as fh:
        fh.write(html)
    log(f"wrote {out} ({os.path.getsize(out)/1e6:.1f} MB)")

    if not args.no_upload:
        try:
            fm_obj.uploadData(out)
            log(f"uploaded to WebServer/{args.out_name}")
        except Exception as e:
            warn(f"upload failed: {e}")


if __name__ == "__main__":
    main()