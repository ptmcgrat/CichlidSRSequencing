"""Extract one chromosome from the genome GFF3 and cache it in cloud storage.

The full annotation is ~463 MB. The browser build needs gene, transcript, exon,
CDS and UTR features for a single chromosome -- roughly 3% of it. Downloading
the whole file on every build is wasteful, so this filters once and uploads the
subset; later runs fetch ~17 MB instead.

    python3 prepareAnnotation.py                       # chr10, default
    python3 prepareAnnotation.py --contig NC_135177.1
    python3 prepareAnnotation.py --force               # rebuild even if cached

The whole chromosome is kept rather than just the region of interest. Gene
intervals in the browser run from the previous gene's end to the next gene's
start, so the genes flanking the region are needed to define the intervals at
its edges -- and the region of interest has already moved once.

Output: <GenomeDir>/<contig>.gff3.gz, bgzipped and tabix-indexed, in the same
feature order and with the same attributes as the source, so `bcftools csq` can
consume it directly for protein-consequence prediction.
"""

import argparse
import gzip
import os
import subprocess
import sys
from collections import Counter

from helper_modules.file_manager import FileManager as FM
from helper_modules import pipeline_checks as pc
from helper_modules.pipeline_checks import PipelineError, log, warn


def parse_args():
    p = argparse.ArgumentParser(description="Cache a per-chromosome GFF3 subset.")
    p.add_argument("--contig", default="NC_135176.1")
    p.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    p.add_argument("--gff", default=None,
                   help="Override the source GFF path. Defaults to FileManager's "
                        "GFF attribute in the genome directory.")
    p.add_argument("--force", action="store_true",
                   help="Rebuild the subset even if a cached copy exists.")
    p.add_argument("--no-upload", action="store_true")
    return p.parse_args()


def open_maybe_gz(path):
    """Open text, detecting gzip by magic bytes rather than by filename."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def find_source_gff(fm_obj, override):
    """Locate the whole-genome GFF.

    FileManager may expose it under any of several attribute names depending on
    when it was added, so try the likely ones before falling back to a scan of
    the genome directory.
    """
    if override:
        return override
    for attr in ("localGFFFile", "localGff3File", "localGFF3File",
                 "localAnnotationFile", "localGTFFile"):
        path = getattr(fm_obj, attr, None)
        if path and str(path).endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
            log(f"source GFF from FileManager.{attr}")
            return path
    gdir = getattr(fm_obj, "localGenomeDir", None)
    if gdir and os.path.isdir(gdir):
        for fn in sorted(os.listdir(gdir)):
            if fn.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
                log(f"source GFF found by scanning {gdir}: {fn}")
                return os.path.join(gdir, fn)
    raise PipelineError(
        "could not locate the genome GFF. Pass --gff with its path, or tell me "
        "which FileManager attribute holds it.")


def filter_contig(src, contig, out_gz):
    """Write every feature on `contig`, preserving header and attributes.

    Comment lines are kept except ##sequence-region records for other
    chromosomes, which would otherwise declare sequences the file no longer
    describes. Attributes are copied verbatim: bcftools csq needs the ID and
    Parent fields to rebuild the gene/transcript/CDS hierarchy.
    """
    raw = out_gz[:-3] if out_gz.endswith(".gz") else out_gz
    kinds = Counter()
    n_kept = n_total = 0
    with open_maybe_gz(src) as fh, open(raw, "w") as out:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) > 1 and parts[1] != contig:
                        continue
                out.write(line)
                continue
            n_total += 1
            if line.split("\t", 1)[0] != contig:
                continue
            n_kept += 1
            kinds[line.split("\t")[2]] += 1
            out.write(line)

    subprocess.run(["bgzip", "-f", raw], check=True)
    # csq and the browser both read by region, so index it.
    try:
        subprocess.run(["tabix", "-f", "-p", "gff", out_gz], check=True)
    except subprocess.CalledProcessError:
        warn("tabix failed -- the file may not be coordinate sorted. "
             "`sort -k1,1 -k4,4n` before bgzip would fix it.")
    return n_total, n_kept, kinds


def main():
    args = parse_args()
    pc.require_tools(("bgzip", "tabix"))

    fm_obj = FM(genome_version=args.genome_version)
    gdir = getattr(fm_obj, "localGenomeDir", None)
    if not gdir:
        raise PipelineError("FileManager has no localGenomeDir attribute")
    os.makedirs(gdir, exist_ok=True)

    out_gz = os.path.join(gdir, f"{args.contig}.gff3.gz")

    if os.path.exists(out_gz) and not args.force:
        log(f"cached subset already present: {out_gz}")
        log("pass --force to rebuild")
        return
    if not args.force:
        try:
            fm_obj.downloadData(out_gz)
            if os.path.exists(out_gz):
                log(f"downloaded cached subset from cloud storage: {out_gz}")
                return
        except FileNotFoundError:
            log("no cached subset in cloud storage; building it")

    src = find_source_gff(fm_obj, args.gff)
    if not os.path.exists(src):
        log(f"downloading the full annotation ({src}) -- this is the slow part")
        fm_obj.downloadData(src)
    pc.require_file(src, "genome GFF", min_bytes=1000)
    log(f"source: {src} ({os.path.getsize(src)/1e6:.0f} MB)")

    n_total, n_kept, kinds = filter_contig(src, args.contig, out_gz)
    if n_kept == 0:
        raise PipelineError(
            f"no features on {args.contig}. Check the contig name matches the "
            f"GFF's first column.")

    log(f"{n_kept:,} of {n_total:,} features kept "
        f"({n_kept/max(1,n_total)*100:.1f}%)")
    log(f"feature types: {dict(kinds.most_common(8))}")
    log(f"wrote {out_gz} ({os.path.getsize(out_gz)/1e6:.1f} MB)")

    for want in ("gene", "mRNA", "exon", "CDS"):
        if want not in kinds:
            warn(f"no '{want}' features found -- "
                 + ("protein-consequence prediction needs CDS"
                    if want == "CDS" else
                    "ASE marker detection needs exon/mRNA features"))

    if not args.no_upload:
        fm_obj.uploadData(out_gz)
        for ext in (".tbi", ".csi"):
            if os.path.exists(out_gz + ext):
                fm_obj.uploadData(out_gz + ext)
        log("uploaded; later builds will fetch this instead of the full GFF")


if __name__ == "__main__":
    main()