"""Genotype candidate QTNs for a single sample.

Two callers contribute to the output: bcftools force-calling for variants whose
ref and alt are both short, and IndelReadClassifier for everything longer. Each
stage now asserts that the number of records it produced equals the number of
sites it was asked about, and the script exits non-zero with a diagnostic if not.

Previously all three stages could produce a valid, indexed, empty VCF and the
run would report success.
"""

import argparse
import os
import shutil
import sys
import subprocess

sys.path.append("..")

import pysam
import pandas as pd

from helper_modules.file_manager import FileManager as FM
from helper_modules.smallVariantGenotyper import genotype_at_sites, BcftoolsError
from helper_modules.InsertionGenotyper import (
    IndelReadClassifier, ClassifierCall, concat_sample_vcfs,
)
from helper_modules import pipeline_checks as pc
from helper_modules.pipeline_checks import Manifest, PipelineError, log, warn


def parse_args():
    p = argparse.ArgumentParser(
        description="Genotype candidate QTNs (small + large variants) for one sample.")
    p.add_argument("SV_VCF", type=str, help="Normalized, bgzipped, indexed VCF of small variants")
    p.add_argument("LV_CSV", type=str, help="CSV of large variants")
    p.add_argument("OUT_VCF", type=str, help="Output bgzipped VCF path")
    p.add_argument("genome_version", type=str, help="Genome version key for FileManager")
    p.add_argument("SampleID", type=str, help="Sample to genotype")
    p.add_argument("--max-missing", type=int, default=0,
                   help="Tolerate up to this many UNEXPLAINED missing sites -- ones "
                        "that had reads passing --min-MQ yet still produced no record. "
                        "Sites with no usable reads are always tolerated and recorded "
                        "separately: a low-coverage sample missing sites it had no data "
                        "for is not a pipeline failure. Default 0.")
    p.add_argument("--allow-partial", action="store_true",
                   help="Write output even when a stage produced fewer records than "
                        "expected. Off by default: a partial VCF that looks complete "
                        "is worse than no VCF.")
    p.add_argument("--keep-bams", action="store_true",
                   help="Keep the downloaded BAM and temp files after the sample "
                        "finishes. Off by default: BAMs are ~10 GB each and 370 of "
                        "them will not fit anywhere.")
    p.add_argument("--min-mapped-reads", type=int, default=1000,
                   help="Minimum reads in the candidate region before genotyping is "
                        "considered meaningful (default 1000)")
    return p.parse_args()


def write_classifier_vcf_v2(calls, sample_name, output_vcf, template_vcf):
    """Write classifier calls with GT:AD:DP:GQ instead of GT:AD.

    IndelReadClassifier already computes quality, n_equal and n_uninformative;
    the original writer discarded all of it. Without GQ there is no way to tell
    "12 ref pairs, 0 alt, 2 unresolvable" from "12 ref, 0 alt, 40 unresolvable" --
    the two look identical downstream but mean very different things.

    Alleles are uppercased here. The master TSV stores them lowercase, so mixing
    these records with the uppercase bcftools records in one file would make
    `bcftools norm -d` and any downstream merge treat g and G as distinct alleles.
    """
    if not output_vcf.endswith(".vcf.gz"):
        raise ValueError("output_vcf must end in .vcf.gz")

    header = subprocess.run(["bcftools", "view", "-h", template_vcf],
                            capture_output=True, text=True, check=True)
    contigs = [l for l in header.stdout.splitlines() if l.startswith("##contig=")]
    if not contigs:
        warn(f"template {template_vcf} contributed no ##contig lines; "
             f"bcftools concat may reject the result")

    calls = sorted(calls, key=lambda c: (c.chrom, c.pos))
    raw_path = output_vcf[:-3]
    with open(raw_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in contigs:
            f.write(c + "\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,'
                'Description="Allelic depths (ref,alt), read pairs">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,'
                'Description="Informative read pairs (ref+alt)">\n')
        f.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,'
                'Description="Phred-scaled genotype quality">\n')
        f.write('##INFO=<ID=SOURCE,Number=1,Type=String,Description="Caller of origin">\n')
        f.write('##INFO=<ID=NEQUAL,Number=1,Type=Integer,'
                'Description="Pairs aligning equally well to both alleles">\n')
        f.write('##INFO=<ID=NUNINF,Number=1,Type=Integer,'
                'Description="Pairs carrying no breakpoint information">\n')
        f.write("#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                 "FILTER", "INFO", "FORMAT", sample_name]) + "\n")
        for c in calls:
            dp = c.ad_ref + c.ad_alt
            gq = int(round(getattr(c, "gq", 0) or 0))
            f.write("\t".join([
                c.chrom, str(c.pos), ".", c.ref.upper(), c.alt.upper(), ".", "PASS",
                f"SOURCE=IndelReadClassifier;NEQUAL={getattr(c, 'n_equal', 0)};"
                f"NUNINF={getattr(c, 'n_uninformative', 0)}",
                "GT:AD:DP:GQ",
                f"{c.gt}:{c.ad_ref},{c.ad_alt}:{dp}:{gq}",
            ]) + "\n")

    subprocess.run(["bgzip", "-f", raw_path], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", output_vcf], check=True)


class RichCall(ClassifierCall):
    """ClassifierCall plus the evidence fields the original tuple dropped."""
    def __new__(cls, gq=0, n_equal=0, n_uninformative=0, **kwargs):
        self = super().__new__(cls, **kwargs)
        self.gq = gq
        self.n_equal = n_equal
        self.n_uninformative = n_uninformative
        return self


def cleanup_sample(fm_obj, sample_id, keep=False):
    """Remove this sample's downloaded BAM and scratch directory.

    Runs whether the sample succeeded or failed. A run that fails on 300 samples
    while keeping every 10 GB BAM fills the disk just as effectively as one that
    succeeds. Everything worth keeping -- the output VCF, the manifest, the error
    log -- lives outside these directories.
    """
    if keep:
        log(f"{sample_id}: keeping BAM and temp files (--keep-bams)")
        return
    freed = 0
    for d in (getattr(fm_obj, "localSampleBamDir", None),
              getattr(fm_obj, "localSampleTempDir", None)):
        if not d or not os.path.isdir(d):
            continue
        try:
            for root, _, files in os.walk(d):
                for fn in files:
                    try:
                        freed += os.path.getsize(os.path.join(root, fn))
                    except OSError:
                        pass
            shutil.rmtree(d, ignore_errors=True)
        except Exception as e:
            warn(f"{sample_id}: could not remove {d}: {e}")
    if freed:
        log(f"{sample_id}: freed {freed / 1e9:.1f} GB")


def main():
    args = parse_args()
    man = Manifest(sample_id=args.SampleID)
    fm_obj = None

    try:
        # ---------------- preflight ----------------
        man.tool_versions = pc.require_tools()
        # Detected once and passed explicitly, so every sample in a run uses the
        # same indel model. Letting it default to True would make the run fail
        # outright on a build without the flag; letting it silently default to
        # False would change the calls without saying so.
        use_indels_2 = pc.mpileup_supports("--indels-2.0")
        if not use_indels_2:
            man.add_warning("this bcftools build does not advertise --indels-2.0; "
                            "falling back to the default indel model, which is less "
                            "accurate for small indels in repeats")

        fm_obj = FM(genome_version=args.genome_version)
        fm_obj.createSampleFiles(args.SampleID, reads=False)

        pc.require_file(fm_obj.localGenomeFile, "reference FASTA", min_bytes=1000)
        pc.require_index(fm_obj.localGenomeFile, "reference FASTA")
        pc.require_file(args.SV_VCF, "small-variant sites VCF")
        pc.require_index(args.SV_VCF, "small-variant sites VCF")
        pc.require_file(args.LV_CSV, "large-variant CSV")

        fm_obj.downloadData(fm_obj.localSampleBamDir)

        pc.require_file(fm_obj.localBamFile, "sample BAM", min_bytes=10000)
        pc.require_index(fm_obj.localBamFile, "sample BAM")
        if not os.path.exists(fm_obj.localDiscordantBamFile):
            man.add_warning("discordant BAM missing; large-indel sensitivity will drop")

        overlap = pc.contig_overlap(fm_obj.localBamFile, fm_obj.localGenomeFile)
        if overlap["shared"] == 0:
            raise PipelineError(
                f"BAM and reference share no contig names. "
                f"BAM has {overlap['bam_only']}, reference has {overlap['ref_only']}")

        expected_sv = pc.count_vcf_records(args.SV_VCF)
        lv_dt = pd.read_csv(args.LV_CSV)
        expected_lv = len(lv_dt)
        man.expected_sv, man.expected_lv = expected_sv, expected_lv
        log(f"{args.SampleID}: expecting {expected_sv} small + {expected_lv} large "
            f"= {expected_sv + expected_lv} records")

        if expected_sv == 0:
            raise PipelineError(f"sites VCF {args.SV_VCF} contains no records")
        if expected_lv == 0:
            man.add_warning(f"{args.LV_CSV} contains no large variants")

        refObj = pysam.FastaFile(fm_obj.localGenomeFile)

        # ---------------- stage 1: small variants ----------------
        sv_temp_vcf = fm_obj.localSampleTempDir + args.SampleID + ".sv.vcf.gz"
        try:
            genotype_at_sites(fm_obj.localBamFile, args.SV_VCF,
                              fm_obj.localGenomeFile, sv_temp_vcf,
                              indels_2=use_indels_2)
        except BcftoolsError as e:
            raise PipelineError(f"small-variant genotyping raised: {e}")

        observed_sv = pc.count_vcf_records(sv_temp_vcf)
        man.observed_sv = observed_sv
        log(f"{args.SampleID}: small-variant stage produced {observed_sv}/{expected_sv}")

        if observed_sv != expected_sv:
            # This is the failure that used to pass silently. Diagnose it here,
            # while the BAM is still local, rather than making someone re-derive it.
            man.diagnostic = pc.diagnose_empty_genotyping(
                fm_obj.localBamFile, args.SV_VCF, fm_obj.localGenomeFile)
            ref_check = pc.check_reference_alleles(args.SV_VCF, fm_obj.localGenomeFile)
            man.diagnostic["reference_alleles"] = ref_check
            # Ask bcftools directly which allele casing it will accept, using this
            # sample's own BAM, rather than leaving it as a manual experiment.
            man.diagnostic["allele_casing"] = pc.probe_allele_casing(
                fm_obj.localBamFile, args.SV_VCF, fm_obj.localGenomeFile,
                fm_obj.localSampleTempDir + "casing_probe/")
            man.missing_sites, man.rerepresented_sites = pc.missing_sites(
                args.SV_VCF, sv_temp_vcf, bam_file=fm_obj.localBamFile)
            msg = (f"small-variant stage produced {observed_sv} of {expected_sv} "
                   f"expected records.\n"
                   f"  unconstrained pileup at the same sites: "
                   f"{man.diagnostic.get('unconstrained_records', 'n/a')} records\n"
                   f"  interpretation: {man.diagnostic.get('interpretation')}\n"
                   f"  sites VCF REF vs genome: {ref_check['n_exact']} exact, "
                   f"{ref_check['n_case_only']} differ only in case, "
                   f"{ref_check['n_mismatch']} genuinely mismatched\n"
                   f"  targets file: {sv_temp_vcf}.targets.tsv.gz")
            casing = man.diagnostic.get("allele_casing", {})
            if casing.get("conclusion"):
                msg += f"\n  allele casing probe: {casing['conclusion']}"
                for mode in ("upper", "lower", "genome"):
                    if mode in casing:
                        msg += (f"\n    {mode:7s} -> {casing[mode]['records']}"
                                f"/{casing[mode]['of']} records")
            if man.rerepresented_sites:
                msg += (f"\n  {len(man.rerepresented_sites)} site(s) returned with "
                        f"different alleles than requested (present in the output, "
                        f"joinable on chrom+pos):")
                for e in man.rerepresented_sites[:5]:
                    msg += (f"\n    {e['chrom']}:{e['pos']} asked {e['ref']}>{e['alt']} "
                            f"[{e['type']}], got {','.join(e['output_alleles'])}")
            if man.missing_sites:
                by_type = {}
                for e in man.missing_sites:
                    by_type[e["type"]] = by_type.get(e["type"], 0) + 1
                msg += f"\n  {len(man.missing_sites)} site(s) ABSENT, by type: {by_type}"
                for e in man.missing_sites[:8]:
                    msg += (f"\n    {e['chrom']}:{e['pos']} {e['ref']}>{e['alt']} "
                            f"[{e['type']}] reads={e.get('depth_raw', '?')} "
                            f"usable(MQ>=20)={e.get('depth_mq20', '?')} "
                            f"MQ0={e.get('frac_mq0', '?')}")
                nodata = [e for e in man.missing_sites if not e.get("depth_raw")]
                unmappable = [e for e in man.missing_sites
                              if e.get("depth_mq20") == 0 and e.get("depth_raw", 0) > 0]
                if nodata:
                    msg += (f"\n  {len(nodata)} with no reads at all -- coverage gap "
                            f"in this sample.")
                if unmappable:
                    msg += (f"\n  {len(unmappable)} with reads but none passing "
                            f"--min-MQ 20 -- multi-mapping repeat sequence, not a "
                            f"caller failure.")
                if man.unexplained_missing:
                    msg += (f"\n  {len(man.unexplained_missing)} had usable reads and "
                            f"still produced no record -- these are the ones worth "
                            f"looking at.")

            # A site with no usable reads was never callable in this sample; that is
            # a coverage fact, not a pipeline fault, and counting it against a fixed
            # budget just penalises low-coverage samples for being low-coverage.
            man.unexplained_missing = [e for e in man.missing_sites
                                       if e.get("depth_mq20", 0) > 0]
            n_unexplained = len(man.unexplained_missing)

            if n_unexplained <= args.max_missing and observed_sv > 0:
                man.add_warning(
                    f"tolerating {len(man.missing_sites)} missing site(s), "
                    f"{n_unexplained} of them unexplained "
                    f"(--max-missing {args.max_missing})\n" + msg)
            elif not args.allow_partial:
                raise PipelineError(msg)
            else:
                man.add_warning(msg)

        # ---------------- stage 2: large variants ----------------
        lv_temp_vcf = fm_obj.localSampleTempDir + args.SampleID + ".lv.vcf.gz"
        calls = []
        for i, row in lv_dt.iterrows():
            classifier = IndelReadClassifier.from_vcf_record(
                ref_genome=refObj, chrom=row.Chromosome, position=row.Position,
                ref=row.Reference.upper(), alt=row.Alt.upper(), flanking=250,
            )
            reads = classifier.fetch_reads_near_insertion(
                bam=fm_obj.localBamFile,
                discordant_bam=(fm_obj.localDiscordantBamFile
                                if os.path.exists(fm_obj.localDiscordantBamFile) else None),
            )
            pair_results = classifier.classify_read_pairs(reads)
            call = classifier.call_genotype(pair_results)
            calls.append(RichCall(
                chrom=row.Chromosome, pos=int(row.Position),
                ref=str(row.Reference), alt=str(row.Alt),
                gt=call.genotype, ad_ref=call.n_ref, ad_alt=call.n_alt,
                gq=call.quality, n_equal=call.n_equal,
                n_uninformative=call.n_uninformative,
            ))

        if len(calls) != expected_lv:
            raise PipelineError(
                f"built {len(calls)} classifier calls from {expected_lv} CSV rows")

        write_classifier_vcf_v2(calls, args.SampleID, lv_temp_vcf,
                                template_vcf=sv_temp_vcf)

        observed_lv = pc.count_vcf_records(lv_temp_vcf)
        man.observed_lv = observed_lv
        if observed_lv != expected_lv:
            raise PipelineError(
                f"large-variant VCF has {observed_lv} records, expected {expected_lv}. "
                f"Duplication or loss happened inside write_classifier_vcf_v2.")

        # ---------------- stage 3: concat ----------------
        sv_names = pc.vcf_sample_names(sv_temp_vcf)
        lv_names = pc.vcf_sample_names(lv_temp_vcf)
        if sv_names != lv_names:
            raise PipelineError(
                f"sample columns disagree before concat: small={sv_names}, "
                f"large={lv_names}. bcftools concat will not combine them correctly.")

        concat_sample_vcfs(sv_temp_vcf, lv_temp_vcf, args.OUT_VCF)

        observed_out = pc.count_vcf_records(args.OUT_VCF)
        man.observed_out = observed_out
        expected_out = observed_sv + observed_lv

        keys = pc.vcf_keys(args.OUT_VCF)
        dups = pc.duplicate_keys(keys)
        man.duplicates = len(dups)
        man.format_fields = sorted(pc.format_fields(args.OUT_VCF))
        man.sources = pc.info_sources(args.OUT_VCF)

        if dups:
            raise PipelineError(
                f"output contains {len(dups)} duplicate (chrom,pos,ref,alt) records, "
                f"e.g. {dups[:3]}. Inputs held {observed_sv} + {observed_lv} unique "
                f"records, so concat is emitting each record more than once.")

        if observed_out != expected_out:
            raise PipelineError(
                f"output has {observed_out} records but inputs held "
                f"{observed_sv} + {observed_lv} = {expected_out}")

        if len(man.sources) < 2 and expected_lv and observed_sv:
            man.add_warning(f"output records all came from one caller: {man.sources}")

        # ---------------- summary stats ----------------
        gt_counts, depths = {}, []
        with pc._open_maybe_gz(args.OUT_VCF) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                fmt, val = f[8].split(":"), f[9].split(":")
                d = dict(zip(fmt, val))
                gt_counts[d.get("GT", "?")] = gt_counts.get(d.get("GT", "?"), 0) + 1
                if "AD" in d and "," in d["AD"]:
                    try:
                        depths.append(sum(int(x) for x in d["AD"].split(",")))
                    except ValueError:
                        pass
        man.genotype_counts = gt_counts
        man.mean_depth = round(sum(depths) / len(depths), 2) if depths else 0.0

        nocall = gt_counts.get("./.", 0)
        if observed_out and nocall / observed_out > 0.5:
            man.add_warning(f"{nocall}/{observed_out} sites are no-calls; "
                            f"mean depth {man.mean_depth}")

        fm_obj.uploadData(args.OUT_VCF)
        man.status = "ok"
        log(f"{args.SampleID}: OK -- {observed_out} records, "
            f"mean depth {man.mean_depth}, genotypes {gt_counts}")

    except Exception as e:
        man.status = "failed"
        man.errors.append(f"{type(e).__name__}: {e}")
        log(f"{args.SampleID}: FAILED -- {e}", "ERROR")
        man.write(args.OUT_VCF + ".manifest.json")
        if fm_obj is not None:
            cleanup_sample(fm_obj, args.SampleID, keep=args.keep_bams)
        sys.exit(1)

    man.write(args.OUT_VCF + ".manifest.json")
    cleanup_sample(fm_obj, args.SampleID, keep=args.keep_bams)


if __name__ == "__main__":
    main()