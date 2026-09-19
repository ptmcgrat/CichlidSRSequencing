"""Extract one chromosome from a sample's whole-genome BAM.

Produces a coordinate-sorted, indexed slice containing every read that a region
query on that chromosome would return from the full BAM -- mapped reads, and
unmapped reads whose mate is mapped there, which aligners store at the mate's
coordinate. That is precisely the set IndelReadClassifier and bcftools mpileup
read, so genotyping from the slice gives identical results to genotyping from
the full BAM.

This replaces the signal-type split (unmapped / discordant / clipped / chimeric)
for genotyping purposes. That split partitioned reads by why they were
interesting, which meant every consumer had to reassemble the pieces -- and the
reassembly was getting it wrong in two ways: the discordant file held only
inter-chromosomal pairs, blind to a transposon inserting near another copy of
itself, and the unmapped reads carrying inserted sequence were in a file nobody
read. A region slice has neither problem.

Deliberately dropped: read pairs where BOTH mates are unmapped. They have no
coordinate, so no region query can return them, and they cannot be localised to
a variant site anyway.

    python3 -m unit_scripts.sliceBam <genome_version> <SampleID> <contig> <out_bam>
"""

import argparse
import os
import shutil
import subprocess
import sys

sys.path.append("..")

import pysam

from helper_modules.file_manager import FileManager as FM
from helper_modules import pipeline_checks as pc
from helper_modules.pipeline_checks import Manifest, PipelineError, log, warn


def parse_args():
    p = argparse.ArgumentParser(description="Slice one chromosome out of a sample BAM.")
    p.add_argument("genome_version", type=str)
    p.add_argument("SampleID", type=str)
    p.add_argument("contig", type=str, help="Chromosome to extract, e.g. NC_135176.1")
    p.add_argument("OUT_BAM", type=str, help="Output slice path (.bam)")
    p.add_argument("--keep-bams", action="store_true",
                   help="Keep the downloaded whole-genome BAM afterwards.")
    p.add_argument("--min-reads", type=int, default=1000,
                   help="Fail if the slice holds fewer reads than this. Guards "
                        "against a silently truncated download (default 1000).")
    return p.parse_args()


def slice_contig(in_bam, contig, out_bam):
    """Copy every read on `contig` into a new sorted, indexed BAM.

    Uses pysam's region iterator rather than shelling out to `samtools view`, so
    the selection is exactly the one the genotyping code will later perform.
    Unmapped reads sitting at a mapped mate's coordinate are included by the
    region query itself -- no special handling needed.
    """
    src = pysam.AlignmentFile(str(in_bam), "rb")
    if contig not in src.references:
        src.close()
        raise PipelineError(
            f"{contig} is not in the BAM header. Present: "
            f"{list(src.references)[:5]}...")

    n_total = n_unmapped = 0
    with pysam.AlignmentFile(str(out_bam), "wb", template=src) as dst:
        for r in src.fetch(contig):
            dst.write(r)
            n_total += 1
            if r.is_unmapped:
                n_unmapped += 1
    src.close()

    pysam.index(str(out_bam))
    return n_total, n_unmapped


def verify_slice(full_bam, slice_bam, contig, n_spot=3):
    """Confirm the slice returns the same reads as the full BAM in sample windows.

    Cheap insurance: if the slice were built wrongly, every downstream genotype
    would be quietly based on less evidence. Compares read-name sets over a few
    windows rather than trusting the copy.
    """
    a = pysam.AlignmentFile(str(full_bam), "rb")
    b = pysam.AlignmentFile(str(slice_bam), "rb")
    length = dict(zip(a.references, a.lengths))[contig]
    mismatches = []
    for i in range(1, n_spot + 1):
        start = int(length * i / (n_spot + 1))
        end = start + 50000
        ka = {(r.query_name, r.is_read1, r.reference_start)
              for r in a.fetch(contig, start, end)}
        kb = {(r.query_name, r.is_read1, r.reference_start)
              for r in b.fetch(contig, start, end)}
        if ka != kb:
            mismatches.append((start, len(ka), len(kb)))
    a.close()
    b.close()
    return mismatches


def main():
    args = parse_args()
    man = Manifest(sample_id=args.SampleID)
    fm_obj = None

    try:
        man.tool_versions = pc.require_tools(("samtools",))
        fm_obj = FM(genome_version=args.genome_version)
        fm_obj.createSampleFiles(args.SampleID, reads=False)

        fm_obj.downloadData(fm_obj.localSampleBamDir)
        pc.require_file(fm_obj.localBamFile, "sample BAM", min_bytes=10000)
        pc.require_index(fm_obj.localBamFile, "sample BAM")

        size_before = os.path.getsize(fm_obj.localBamFile)
        os.makedirs(os.path.dirname(args.OUT_BAM), exist_ok=True)

        log(f"{args.SampleID}: slicing {args.contig}")
        n_total, n_unmapped = slice_contig(
            fm_obj.localBamFile, args.contig, args.OUT_BAM)
        man.observed_out = n_total
        log(f"{args.SampleID}: {n_total:,} reads ({n_unmapped:,} unmapped) "
            f"on {args.contig}")

        if n_total < args.min_reads:
            raise PipelineError(
                f"slice holds only {n_total} reads, below --min-reads "
                f"{args.min_reads}. The source BAM may be truncated or the "
                f"contig name wrong.")

        bad = verify_slice(fm_obj.localBamFile, args.OUT_BAM, args.contig)
        if bad:
            raise PipelineError(
                f"slice does not match the full BAM in {len(bad)} sampled "
                f"window(s): {bad}. Not uploading it.")

        size_after = os.path.getsize(args.OUT_BAM)
        man.mean_depth = round(size_after / max(1, size_before) * 100, 2)
        log(f"{args.SampleID}: {size_before/1e9:.1f} GB -> "
            f"{size_after/1e9:.2f} GB ({man.mean_depth}% of the original)")

        fm_obj.uploadData(args.OUT_BAM)
        for ext in (".bai", ".csi"):
            if os.path.exists(args.OUT_BAM + ext):
                fm_obj.uploadData(args.OUT_BAM + ext)
        man.status = "ok"

    except Exception as e:
        man.status = "failed"
        man.errors.append(f"{type(e).__name__}: {e}")
        log(f"{args.SampleID}: FAILED -- {e}", "ERROR")
        man.write(args.OUT_BAM + ".manifest.json")
        if fm_obj is not None and not args.keep_bams:
            shutil.rmtree(getattr(fm_obj, "localSampleBamDir", "") or "/nonexistent",
                          ignore_errors=True)
        sys.exit(1)

    man_path = man.write(args.OUT_BAM + ".manifest.json")
    try:
        fm_obj.uploadData(man_path)
    except Exception as e:
        warn(f"{args.SampleID}: manifest upload failed: {e}")

    if not args.keep_bams:
        # The whole-genome BAM is the thing we are trying to stop downloading;
        # keeping it after slicing would defeat the point.
        shutil.rmtree(fm_obj.localSampleBamDir, ignore_errors=True)
        log(f"{args.SampleID}: removed the whole-genome BAM")


if __name__ == "__main__":
    main()