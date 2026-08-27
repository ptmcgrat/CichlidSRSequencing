"""Driver: genotype candidate QTNs across all aligned Deep and Shallow Benthic samples.

Structure:
  1. Preflight   -- validate the reference, the candidate table, and the sample
                    database before launching anything.
  2. Smoke test  -- run ONE sample to completion and check its manifest. If the
                    pipeline is broken, this catches it in a few minutes instead
                    of after several hundred parallel jobs have each produced a
                    valid-looking, mostly-empty VCF.
  3. Fan out     -- run the rest, reading each manifest back rather than trusting
                    exit codes.
  4. Report      -- write a run-level summary and keep the logs of anything that
                    failed or came back short.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter
from types import SimpleNamespace

import pandas as pd

from helper_modules.file_manager import FileManager as FM
from helper_modules.smallVariantGenotyper import normalize_sites
from helper_modules import pipeline_checks as pc
from helper_modules.pipeline_checks import PipelineError, log, warn


def parse_args():
    p = argparse.ArgumentParser(description="Genotype candidate QTNs across samples.")
    p.add_argument("--candidates", default=None,
                   help="Master candidate table (TSV). Defaults to "
                        "<localNikeshDir>/QTG_Candidates/candidateQTNs_all.tsv, "
                        "downloaded via FileManager. Pass a path to use a local file "
                        "instead (no download attempted).")
    p.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    p.add_argument("--ecogroups", nargs="+", default=["Deep_Benthic", "Shallow_Benthic"])
    p.add_argument("--threshold", type=int, default=10,
                   help="Ref/alt length at or below which a variant goes to bcftools")
    p.add_argument("--num-parallel", type=int, default=48)
    p.add_argument("--keep-bams", action="store_true",
                   help="Keep each sample's downloaded BAM after it finishes. "
                        "Off by default -- BAMs average 10 GB and 370 of them is "
                        "roughly 4.8 TB.")
    p.add_argument("--max-missing", type=int, default=0,
                   help="Per sample, tolerate up to this many requested small-variant "
                        "sites being absent from the output. They are still recorded "
                        "in each manifest. Default 0.")
    p.add_argument("--skip-smoke-test", action="store_true",
                   help="Skip the single-sample validation run. Not recommended.")
    p.add_argument("--preflight-only", action="store_true",
                   help="Run the checks and exit without genotyping anything.")
    p.add_argument("--force", action="store_true",
                   help="Continue past preflight warnings that would otherwise stop the run.")
    return p.parse_args()


def load_candidates(fm_obj, args):
    """Resolve and read the master candidate table.

    With no --candidates argument the canonical copy is pulled from cloud storage,
    so every run starts from the same table rather than whatever happens to be in
    the working directory. An explicit path is used as-is and never downloaded,
    which is what you want when testing a modified table.
    """
    if args.candidates is None:
        path = fm_obj.localNikeshDir + "QTG_Candidates/candidateQTNs_all.tsv"
        try:
            fm_obj.downloadData(path)
        except FileNotFoundError as e:
            raise PipelineError(
                f"could not download the candidate table: {e}. Pass --candidates "
                f"to point at a local file instead.")
        log(f"candidate table downloaded to {path}")
    else:
        path = args.candidates
        pc.require_file(path, "candidate table")
        log(f"using local candidate table {path}")

    args.candidates = path
    dt = pd.read_csv(path, sep="\t")
    if dt.empty:
        raise PipelineError(f"candidate table {path} has no rows")
    return dt


def ensure_reference(fm_obj):
    """Download the reference and make sure its .fai is present.

    FileManager.downloadData copies a single file, so asking for the FASTA does
    not bring the index with it. Try the cloud copy first, then build one locally
    -- faidx on an existing FASTA is cheap and deterministic.
    """
    fm_obj.downloadData(fm_obj.localGenomeFile)
    fai = fm_obj.localGenomeFile + ".fai"
    if os.path.exists(fai):
        return
    try:
        fm_obj.downloadData(fai)
        log("downloaded reference .fai")
        return
    except FileNotFoundError:
        log("no .fai in cloud storage; building one locally with samtools faidx")
    r = subprocess.run(["samtools", "faidx", fm_obj.localGenomeFile],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise PipelineError(f"samtools faidx failed: {r.stderr}")


def print_sv(outfile, dt):
    """Write the small-variant sites VCF.

    Alleles are uppercased to match VCF convention. Note that the reference FASTA
    may be soft-masked, in which case these uppercase alleles will not string-match
    the lowercase genome under `bcftools call -C alleles`. Preflight checks for
    exactly that.
    """
    with open(outfile, "w") as fp:
        print("##fileformat=VCFv4.2", file=fp)
        print("##source=MSAUniqueVariantCaller", file=fp)
        print("##sample=Y_reg", file=fp)
        print("##contig=<ID=NC_135176.1>", file=fp)
        print('##INFO=<ID=TYPE,Number=1,Type=String,'
              'Description="SUBSTITUTE, INSERTION, or DELETION">', file=fp)
        print('##INFO=<ID=ALN_START,Number=1,Type=Integer,'
              'Description="1-based alignment column where the variant starts">', file=fp)
        print("\t".join(["#CHROM", "POS", "ID", "REF", "ALT",
                         "QUAL", "FILTER", "INFO"]), file=fp)
        for i, row in dt.iterrows():
            name = row.Name if isinstance(row.Name, str) and row.Name.strip() else "."
            print("\t".join([row.Chromosome, str(row.Position), name,
                             row.Reference.upper(), row.Alt.upper(),
                             str(row.Q), "PASS", row.Info]), file=fp)


def preflight(args, fm_obj, dt):
    """Everything that can be checked before a single BAM is touched."""
    problems, cautions = [], []

    log("=== preflight ===")
    versions = pc.require_tools()
    if not pc.mpileup_supports("--indels-2.0"):
        cautions.append("bcftools build does not advertise --indels-2.0")

    pc.require_file(fm_obj.localGenomeFile, "reference FASTA", min_bytes=1000)
    pc.require_index(fm_obj.localGenomeFile, "reference FASTA")

    # --- candidate table ---
    log(f"candidate table: {len(dt)} rows")
    if dt.Name.astype(str).str.strip().isin([".", "", "nan"]).all():
        cautions.append("every candidate has Name='.', so records can only be joined "
                        "back on chrom+pos+ref+alt after normalization")
    dup_pos = dt.duplicated(subset=["Chromosome", "Position", "Reference", "Alt"]).sum()
    if dup_pos:
        problems.append(f"{dup_pos} duplicate rows in {args.candidates}")

    for col in ["Chromosome", "Position", "Reference", "Alt", "Q", "Info"]:
        if col not in dt.columns:
            problems.append(f"candidate table is missing column {col}")
    if dt.Reference.isna().any() or dt.Alt.isna().any():
        problems.append("candidate table has null Reference or Alt values")

    return problems, cautions, versions


def preflight_sites(sv_norm_vcf_file, genome_file, n_input_sv):
    """Checks that need the normalized sites VCF to exist."""
    problems, cautions = [], []

    n_norm = pc.count_vcf_records(sv_norm_vcf_file)
    log(f"normalized sites VCF: {n_norm} records (from {n_input_sv} input rows)")
    if n_norm != n_input_sv:
        cautions.append(f"bcftools norm changed the record count: "
                        f"{n_input_sv} in, {n_norm} out")

    dups = pc.duplicate_keys(pc.vcf_keys(sv_norm_vcf_file))
    if dups:
        problems.append(f"normalized sites VCF has {len(dups)} duplicate records, "
                        f"e.g. {dups[:3]}")

    # The check that would have caught this run's failure before it started.
    ref_check = pc.check_reference_alleles(sv_norm_vcf_file, genome_file)
    log(f"sites REF vs genome: {ref_check['n_exact']} exact, "
        f"{ref_check['n_case_only']} case-only, {ref_check['n_mismatch']} mismatched")

    if ref_check["n_mismatch"]:
        problems.append(
            f"{ref_check['n_mismatch']} sites have REF alleles that disagree with the "
            f"reference genome. Examples: {ref_check['examples'][:3]}. "
            f"These will be dropped silently by `bcftools call -C alleles`.")

    if ref_check["n_case_only"]:
        problems.append(
            f"{ref_check['n_case_only']} sites still disagree with the reference in case "
            f"after case-matching. That should not happen -- inspect "
            f"{sv_norm_vcf_file}.")

    mask = pc.softmask_summary(sv_norm_vcf_file, genome_file)
    log(f"candidate sites on masked sequence: {mask['masked']}/"
        f"{mask['masked'] + mask['unmasked']}")

    return problems, cautions


def preflight_samples(fm_obj, ecogroups):
    """Sample-database sanity. These are cautions, not blockers -- but they change
    how the resulting genotypes can be interpreted, so they get reported loudly."""
    cautions = []
    sample_dt, alignment_dt = fm_obj.sample_dt, fm_obj.alignment_dt

    dupes = sample_dt[sample_dt.SampleID.duplicated(keep=False)]
    if len(dupes):
        cautions.append(
            f"{dupes.SampleID.nunique()} SampleID(s) appear more than once in the sample "
            f"database with conflicting metadata: "
            f"{sorted(dupes.SampleID.unique())[:5]}. Deduplicate before relying on "
            f"Subgroup/Category groupings.")

    in_scope = sample_dt[sample_dt.Ecogroup.isin(ecogroups)]
    aligned = alignment_dt[alignment_dt.SampleID.isin(in_scope.SampleID)]
    aligned_meta = in_scope[in_scope.SampleID.isin(aligned.SampleID)]

    log(f"samples in scope: {len(in_scope)} in database, {len(aligned_meta)} with alignments")

    for col in ["Subgroup", "Category"]:
        if col not in sample_dt.columns:
            cautions.append(f"sample database has no {col} column")
            continue
        missing = aligned_meta[col].isna().sum()
        if missing:
            cautions.append(f"{missing}/{len(aligned_meta)} aligned samples have no {col}")
        log(f"{col}: {aligned_meta[col].nunique()} distinct values")

    sexes = aligned_meta.Sex.fillna("(blank)").value_counts().to_dict()
    usable = sum(v for k, v in sexes.items() if k in ("M", "F"))
    log(f"sex: {sexes}")
    if usable < len(aligned_meta) * 0.8:
        cautions.append(
            f"only {usable}/{len(aligned_meta)} aligned samples have a definite M/F sex. "
            f"Any sex-association statistic will be computed on that subset, not the "
            f"full cohort.")

    # Sex recorded in the Sex column vs sex implied by a Category label.
    if "Category" in aligned_meta.columns:
        for cat, expected in [("YH_Males", "M"), ("YH_Females", "F"),
                              ("MC_males", "M"), ("MC_females", "F"),
                              ("CV_males", "M"), ("CV_females", "F")]:
            sub = aligned_meta[aligned_meta.Category == cat]
            bad = sub[sub.Sex.isin(["M", "F"]) & (sub.Sex != expected)]
            if len(bad):
                cautions.append(
                    f"{len(bad)} sample(s) with Category={cat} have Sex="
                    f"{sorted(bad.Sex.unique())}: {sorted(bad.SampleID)[:5]}")

    missing_bam = set(aligned.SampleID) - set(alignment_dt.SampleID)
    if missing_bam:
        cautions.append(f"{len(missing_bam)} samples lack alignment records")

    if "Coverage" in alignment_dt.columns:
        cov = alignment_dt[alignment_dt.SampleID.isin(aligned_meta.SampleID)].Coverage
        low = (cov < 5).sum()
        log(f"coverage: median {cov.median():.1f}x, {low} samples below 5x")
        if low:
            cautions.append(f"{low} samples have coverage below 5x and will produce "
                            f"mostly no-calls")

    return sorted(aligned_meta.SampleID.unique()), cautions



def preflight_disk(fm_obj, sample_ids, num_parallel, keep_bams):
    """Will the concurrent BAM downloads fit on disk?

    Each worker downloads its sample's BAM directory before genotyping, so peak
    usage is roughly num_parallel BAMs at once -- and they are not uniform: the
    largest in this cohort is six times the median. Sizing on the median is how
    you fill a filesystem at 3am.
    """
    cautions = []
    adt = fm_obj.alignment_dt
    if "BamSize" not in adt.columns:
        cautions.append("alignment database has no BamSize column; cannot estimate "
                        "disk usage")
        return cautions

    sizes = pd.to_numeric(
        adt[adt.SampleID.isin(sample_ids)].BamSize, errors="coerce").dropna()
    if sizes.empty:
        return cautions

    target = fm_obj.localBamfilesDir
    while target and not os.path.isdir(target):
        target = os.path.dirname(target.rstrip("/"))
    free = shutil.disk_usage(target or "/").free

    # 1.35x covers the discordant BAM and indexes downloaded alongside the main one.
    p90 = float(sizes.quantile(0.9)) * 1.35
    peak = p90 * num_parallel
    total = float(sizes.sum()) * 1.35

    log(f"BAM sizes: median {sizes.median()/1e9:.1f} GB, "
        f"max {sizes.max()/1e9:.1f} GB, {len(sizes)} samples")
    log(f"disk free at {target}: {free/1e9:.0f} GB")
    log(f"estimated peak usage with {num_parallel} concurrent: {peak/1e9:.0f} GB")

    if keep_bams:
        log(f"--keep-bams is set: total retained would be {total/1e12:.1f} TB")
        if total > free:
            cautions.append(
                f"--keep-bams needs about {total/1e12:.1f} TB but only "
                f"{free/1e9:.0f} GB is free. The run will fill the disk.")
    elif peak > free * 0.8:
        safe = max(1, int(free * 0.6 / p90))
        cautions.append(
            f"peak BAM usage (~{peak/1e9:.0f} GB with {num_parallel} concurrent) is "
            f"close to or above the {free/1e9:.0f} GB free. Consider "
            f"--num-parallel {safe}.")
    return cautions


def run_one(sampleID, command, error_file):
    """Run a single sample synchronously and return (returncode, manifest_or_None)."""
    with open(error_file, "w") as fp:
        proc = subprocess.run(command, stderr=fp, stdout=subprocess.DEVNULL)
    return proc.returncode


def check_manifest(out_vcf, sampleID, max_missing=0):
    """Read back what the worker actually produced. Exit codes are not enough --
    the original failure mode was a clean exit over an incomplete file.

    Sites the worker was told to tolerate are subtracted from the expected count,
    otherwise the driver rejects exactly the output it authorised. What must still
    hold is that every requested site is accounted for: written, or recorded as
    missing. A record that is neither is unexplained loss.
    """
    path = out_vcf + ".manifest.json"
    if not os.path.exists(path):
        return None, f"{sampleID}: no manifest written"
    m = json.load(open(path))
    if m["status"] != "ok":
        return m, f"{sampleID}: status={m['status']} errors={m['errors']}"

    expected = m["expected_sv"] + m["expected_lv"]
    n_missing = len(m.get("missing_sites", []))

    if n_missing > max_missing:
        return m, (f"{sampleID}: {n_missing} sites missing, above --max-missing "
                   f"{max_missing}")
    if m["observed_out"] + n_missing != expected:
        unexplained = expected - n_missing - m["observed_out"]
        return m, (f"{sampleID}: {m['observed_out']} records + {n_missing} recorded "
                   f"missing != {expected} expected ({unexplained} unaccounted for)")
    if m["duplicates"]:
        return m, f"{sampleID}: {m['duplicates']} duplicate records"
    return m, None


def main():
    args = parse_args()

    fm_obj = FM(genome_version=args.genome_version)
    fm_obj.readSampleDatabase()
    fm_obj.readAlignmentDatabase()

    out_dir = fm_obj.localNikeshDir + "QTG_Candidates/"
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(fm_obj.localErrorsDir, exist_ok=True)
    ensure_reference(fm_obj)

    dt = load_candidates(fm_obj, args)

    problems, cautions, versions = preflight(args, fm_obj, dt)

    # Build the two site files.
    sv_vcf_file = out_dir + "candidateQTNs_sv.vcf"
    sv_norm_vcf_file = out_dir + "candidateQTNs_sv.norm.vcf.gz"
    lv_csv_file = out_dir + "candidateQTNs_lv.csv"

    is_small = (dt.Alt.str.len() <= args.threshold) & (dt.Reference.str.len() <= args.threshold)
    sv_dt, lv_dt = dt[is_small], dt[~is_small]
    if len(sv_dt) + len(lv_dt) != len(dt):
        problems.append("small/large partition does not cover the candidate table")
    log(f"partition at threshold {args.threshold}: {len(sv_dt)} small, {len(lv_dt)} large")

    print_sv(sv_vcf_file, sv_dt)
    normalize_sites(sv_vcf_file, fm_obj.localGenomeFile, sv_norm_vcf_file)
    lv_dt.to_csv(lv_csv_file, index=False)

    # Match REF allele case to the reference before anything reads this file.
    # The genome is soft-masked; a sites VCF written with .upper() disagrees with
    # it at every repeat-masked position. Harmless if bcftools is case-insensitive,
    # necessary if it is not -- so just do it rather than testing for it. The
    # reference itself is never modified.
    sv_cased_vcf_file = out_dir + "candidateQTNs_sv.norm.cased.vcf.gz"
    case_stats = pc.case_match_sites(sv_norm_vcf_file, fm_obj.localGenomeFile,
                                     sv_cased_vcf_file)
    log(f"allele case vs reference: {case_stats['already_matching']} already matched, "
        f"{case_stats['case_corrected']} corrected, "
        f"{case_stats['true_mismatch']} genuinely mismatched")
    sv_norm_vcf_file = sv_cased_vcf_file

    p2, c2 = preflight_sites(sv_norm_vcf_file, fm_obj.localGenomeFile, len(sv_dt))
    problems += p2
    cautions += c2

    aligned_samples, c3 = preflight_samples(fm_obj, args.ecogroups)
    cautions += c3
    cautions += preflight_disk(fm_obj, aligned_samples, args.num_parallel,
                               args.keep_bams)

    log("=== preflight summary ===")
    for c in cautions:
        warn(c)
    for p in problems:
        log(p, "ERROR")

    if problems and not args.force:
        log(f"{len(problems)} blocking problem(s). Fix them, or re-run with --force "
            f"to proceed anyway.", "ERROR")
        sys.exit(1)
    if args.preflight_only:
        log("preflight complete (--preflight-only)")
        sys.exit(0)
    if not aligned_samples:
        log("no samples to genotype", "ERROR")
        sys.exit(1)

    def cmd_for(sampleID):
        out_vcf = out_dir + sampleID + "_candidate_QTNs.vcf.gz"
        cmd = ["python", "-m", "unit_scripts.genotypeCandidates",
               sv_norm_vcf_file, lv_csv_file, out_vcf,
               args.genome_version, sampleID,
               "--max-missing", str(args.max_missing)]
        if args.keep_bams:
            cmd.append("--keep-bams")
        return out_vcf, cmd

    # ---------------- smoke test ----------------
    if not args.skip_smoke_test:
        probe = aligned_samples[0]
        log(f"=== smoke test on {probe} ===")
        out_vcf, command = cmd_for(probe)
        err = fm_obj.localErrorsDir + "QTGFinder_" + probe + "_errors.txt"
        rc = run_one(probe, command, err)
        man, problem = check_manifest(out_vcf, probe, args.max_missing)
        if rc != 0 or problem:
            log(f"smoke test failed: {problem or f'exit code {rc}'}", "ERROR")
            log(f"worker stderr: {err}", "ERROR")
            if man and man.get("status") != "ok" and man.get("diagnostic"):
                log("diagnostic: " + json.dumps(man["diagnostic"], indent=2), "ERROR")
            log("Not launching the remaining samples.", "ERROR")
            sys.exit(1)
        n_missing = len(man.get("missing_sites", []))
        log(f"smoke test passed: {man['observed_out']} records "
            f"({man['observed_sv']} small + {man['observed_lv']} large), "
            f"mean depth {man['mean_depth']}"
            + (f", {n_missing} site(s) uncallable and recorded" if n_missing else ""))
        aligned_samples = aligned_samples[1:]

    # ---------------- fan out ----------------
    commands = []
    for sampleID in aligned_samples:
        out_vcf, command = cmd_for(sampleID)
        commands.append(SimpleNamespace(
            sampleID=sampleID, command=command, out_vcf=out_vcf,
            error_file=fm_obj.localErrorsDir + "QTGFinder_" + sampleID + "_errors.txt",
            process=None, error_fp=None))

    pending = list(commands)
    running, results, failures = [], {}, []
    total = len(pending)
    log(f"=== genotyping {total} samples, {args.num_parallel} at a time ===")

    def launch(data):
        data.error_fp = open(data.error_file, "w")
        data.process = subprocess.Popen(data.command, stderr=data.error_fp,
                                        stdout=subprocess.DEVNULL)
        running.append(data)

    while pending and len(running) < args.num_parallel:
        launch(pending.pop(0))

    done = 0
    while running:
        time.sleep(1)
        for data in [x for x in running if x.process.poll() is not None]:
            data.error_fp.close()
            running.remove(data)
            done += 1

            man, problem = check_manifest(data.out_vcf, data.sampleID,
                                          args.max_missing)
            results[data.sampleID] = man
            if data.process.returncode != 0 or problem:
                failures.append(problem or f"{data.sampleID}: exit "
                                           f"{data.process.returncode}")
                log(f"[{done}/{total}] {data.sampleID} FAILED -- "
                    f"{problem or data.process.returncode} (log: {data.error_file})",
                    "ERROR")
            else:
                # Only discard the log when the output has been verified.
                if os.path.exists(data.error_file):
                    os.remove(data.error_file)
                n_missing = len(man.get("missing_sites", []))
                log(f"[{done}/{total}] {data.sampleID} ok "
                    f"({man['observed_out']} records, depth {man['mean_depth']}"
                    + (f", {n_missing} missing" if n_missing else "") + ")")

            if pending:
                launch(pending.pop(0))

    # ---------------- report ----------------
    ok = {k: v for k, v in results.items() if v and v.get("status") == "ok"}
    log("=== run summary ===")
    log(f"{len(ok)}/{total} samples produced verified output")

    if ok:
        counts = Counter(v["observed_out"] for v in ok.values())
        log(f"record counts across samples: {dict(counts)}")
        if len(counts) > 1:
            warn("samples disagree on record count -- the matrix will be ragged")
        depths = [v["mean_depth"] for v in ok.values()]
        log(f"mean depth: min {min(depths):.1f}, median "
            f"{sorted(depths)[len(depths) // 2]:.1f}, max {max(depths):.1f}")
        nocall = [(k, v["genotype_counts"].get("./.", 0) / max(1, v["observed_out"]))
                  for k, v in ok.items()]
        bad = [k for k, f in nocall if f > 0.5]
        if bad:
            warn(f"{len(bad)} samples are more than half no-calls: {bad[:10]}")

        # Which sites go missing, and in how many samples. A site missing
        # everywhere is a property of the site; a site missing in a few samples
        # is a property of those samples' coverage.
        site_counts, rerep_counts = Counter(), Counter()
        for v in ok.values():
            for e in v.get("missing_sites", []):
                site_counts[(e["chrom"], e["pos"], e["type"])] += 1
            for e in v.get("rerepresented_sites", []):
                rerep_counts[(e["chrom"], e["pos"], e["type"])] += 1

        if rerep_counts:
            log(f"{len(rerep_counts)} site(s) returned with re-anchored alleles "
                f"(present in output; join on chrom+pos, not on alleles):")
            for (c, p, t), n in rerep_counts.most_common(10):
                log(f"    {c}:{p} [{t}] in {n}/{len(ok)} samples")

        if site_counts:
            log(f"{len(site_counts)} distinct site(s) ABSENT in at least one sample:")
            for (c, p, t), n in site_counts.most_common(15):
                log(f"    {c}:{p} [{t}] absent in {n}/{len(ok)} samples")
            universal = [k for k, n in site_counts.items() if n == len(ok)]
            if universal:
                warn(f"{len(universal)} site(s) absent in EVERY sample -- uncallable by "
                     f"this pipeline rather than sample-specific dropouts")

    if failures:
        log(f"{len(failures)} failures:", "ERROR")
        for f in failures:
            log("  " + f, "ERROR")
        log(f"worker logs kept in {fm_obj.localErrorsDir}", "ERROR")

    summary_path = out_dir + "run_summary.json"
    with open(summary_path, "w") as fh:
        json.dump({
            "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "tool_versions": versions,
            "expected_records": len(sv_dt) + len(lv_dt),
            "n_requested": total,
            "n_ok": len(ok),
            "failures": failures,
            "preflight_problems": problems,
            "preflight_cautions": cautions,
            "samples": results,
        }, fh, indent=2)
    log(f"summary written to {summary_path}")

    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()