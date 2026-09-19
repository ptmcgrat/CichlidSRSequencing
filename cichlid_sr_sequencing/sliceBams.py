"""Build per-chromosome BAM slices for every aligned sample, once.

Genotyping chr10 currently downloads a ~10 GB whole-genome BAM to read ~40 Mb of
it. Slicing that chromosome out once, into cloud storage, makes every subsequent
run over it ~30x cheaper to fetch -- which matters less for one run than for the
many re-runs that follow a changed variant set, a changed threshold, or a second
pass over the same region.

    python3 sliceBams.py --contig NC_135176.1
    python3 sliceBams.py --contig NC_135176.1 --resume --shard 1/2

Slices land in <localNikeshDir>/BamSlices/<contig>/<SampleID>.<contig>.bam and
are uploaded with their indexes. Each sample writes a manifest, and --resume
reads those back rather than trusting exit codes, exactly as the genotyping
drivers do.

Structured deliberately like genotypeY_QTNs.py: same preflight, sync, resume,
shard, manifest and reporting machinery, because the failure modes are the same
and the operational habits should be too.
"""

import argparse
import glob
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
from helper_modules import pipeline_checks as pc
from helper_modules.pipeline_checks import PipelineError, log, warn


def parse_args():
    p = argparse.ArgumentParser(description="Slice one chromosome from every sample BAM.")
    p.add_argument("--contig", default="NC_135176.1",
                   help="Chromosome to extract (default NC_135176.1, the sex "
                        "chromosome).")
    p.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    p.add_argument("--ecogroups", nargs="+", default=None,
                   help="Restrict to these ecogroups. Default is every aligned "
                        "sample, since a slice is useful to any later analysis.")
    p.add_argument("--num-parallel", type=int, default=24,
                   help="Concurrent samples. Slicing is I/O-bound rather than "
                        "CPU-bound, so this is usually limited by network and "
                        "disk, not cores (default 24).")
    p.add_argument("--out-subdir", default="BamSlices",
                   help="Folder under the Nikesh directory for the slices.")
    p.add_argument("--resume", action="store_true",
                   help="Skip samples that already have a verified slice.")
    p.add_argument("--no-sync", action="store_true",
                   help="Skip the initial listing of existing slices in the cloud.")
    p.add_argument("--shard", default=None, metavar="i/n",
                   help="Process only shard i of n, so several servers can share "
                        "the work without duplicating it.")
    p.add_argument("--keep-bams", action="store_true",
                   help="Keep whole-genome BAMs after slicing. Off by default.")
    p.add_argument("--preflight-only", action="store_true")
    return p.parse_args()


def preflight_disk(fm_obj, sample_ids, num_parallel):
    """Peak disk is driven by the whole-genome BAMs being downloaded, not the
    slices -- the slices are the thing that makes this unnecessary next time."""
    cautions = []
    adt = fm_obj.alignment_dt
    if "BamSize" not in adt.columns:
        return cautions
    sizes = pd.to_numeric(adt[adt.SampleID.isin(sample_ids)].BamSize,
                          errors="coerce").dropna()
    if sizes.empty:
        return cautions
    target = fm_obj.localBamfilesDir
    while target and not os.path.isdir(target):
        target = os.path.dirname(target.rstrip("/"))
    free = shutil.disk_usage(target or "/").free
    p90 = float(sizes.quantile(0.9)) * 1.35
    peak = p90 * num_parallel
    log(f"BAM sizes: median {sizes.median()/1e9:.1f} GB, max {sizes.max()/1e9:.1f} GB")
    log(f"disk free: {free/1e9:.0f} GB | estimated peak with {num_parallel} "
        f"concurrent: {peak/1e9:.0f} GB")
    log(f"total to transfer once: {sizes.sum()/1e12:.1f} TB")
    if peak > free * 0.8:
        safe = max(1, int(free * 0.6 / p90))
        cautions.append(f"peak usage is close to free space; consider "
                        f"--num-parallel {safe}")
    return cautions


def check_manifest(out_bam, sampleID):
    path = out_bam + ".manifest.json"
    if not os.path.exists(path) or not os.path.exists(out_bam):
        return None, f"{sampleID}: no slice or manifest"
    try:
        m = json.load(open(path))
    except Exception as e:
        return None, f"{sampleID}: unreadable manifest ({e})"
    if m.get("status") != "ok":
        return m, f"{sampleID}: status={m.get('status')} {m.get('errors')}"
    if not m.get("observed_out"):
        return m, f"{sampleID}: slice recorded zero reads"
    return m, None


def apply_shard(samples, spec):
    if not spec:
        return samples
    try:
        i, n = (int(x) for x in spec.split("/"))
    except ValueError:
        raise PipelineError(f"--shard must look like i/n, got {spec!r}")
    picked = [s for k, s in enumerate(sorted(samples)) if k % n == i - 1]
    log(f"--shard {i}/{n}: {len(picked)} of {len(samples)} samples here")
    return picked


def main():
    args = parse_args()
    pc.require_tools(("samtools",))

    fm_obj = FM(genome_version=args.genome_version)
    fm_obj.readSampleDatabase()
    fm_obj.readAlignmentDatabase()

    out_dir = (fm_obj.localNikeshDir + args.out_subdir.strip("/") + "/"
               + args.contig + "/")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(fm_obj.localErrorsDir, exist_ok=True)

    if not args.no_sync:
        log("checking cloud storage for existing slices")
        try:
            fm_obj.downloadData(out_dir.rstrip("/"))
            n = len(glob.glob(out_dir + "*.bam"))
            log(f"{n} slice(s) already present")
        except FileNotFoundError:
            log("no slices in cloud storage yet")
        except Exception as e:
            warn(f"sync failed ({e}); --resume will see local files only")

    sample_dt, alignment_dt = fm_obj.sample_dt, fm_obj.alignment_dt
    in_scope = sample_dt
    if args.ecogroups:
        in_scope = sample_dt[sample_dt.Ecogroup.isin(args.ecogroups)]
    samples = sorted(set(in_scope.SampleID) & set(alignment_dt.SampleID))
    log(f"{len(samples)} aligned sample(s) in scope for {args.contig}")

    cautions = preflight_disk(fm_obj, samples, args.num_parallel)
    for c in cautions:
        warn(c)
    if args.preflight_only:
        log("preflight complete (--preflight-only)")
        sys.exit(0)
    if not samples:
        log("nothing to do", "ERROR")
        sys.exit(1)

    def out_for(sid):
        return out_dir + f"{sid}.{args.contig}.bam"

    if args.resume:
        done = [s for s in samples if check_manifest(out_for(s), s)[1] is None]
        if done:
            samples = [s for s in samples if s not in done]
            log(f"--resume: skipping {len(done)} verified slice(s); "
                f"{len(samples)} to build")
        if not samples:
            log("nothing left to do")
            sys.exit(0)

    samples = apply_shard(samples, args.shard)
    if not samples:
        log("no samples in this shard")
        sys.exit(0)

    commands = []
    for sid in samples:
        cmd = ["python", "-m", "unit_scripts.sliceBam",
               args.genome_version, sid, args.contig, out_for(sid)]
        if args.keep_bams:
            cmd.append("--keep-bams")
        commands.append(SimpleNamespace(
            sampleID=sid, command=cmd, out_bam=out_for(sid),
            error_file=fm_obj.localErrorsDir + f"Slice_{sid}_errors.txt",
            process=None, error_fp=None))

    pending, running, results, failures = list(commands), [], {}, []
    total = len(pending)
    log(f"=== slicing {total} sample(s), {args.num_parallel} at a time ===")

    def launch(d):
        d.error_fp = open(d.error_file, "w")
        d.process = subprocess.Popen(d.command, stderr=d.error_fp,
                                     stdout=subprocess.DEVNULL)
        running.append(d)

    while pending and len(running) < args.num_parallel:
        launch(pending.pop(0))

    done_n = 0
    while running:
        time.sleep(1)
        for d in [x for x in running if x.process.poll() is not None]:
            d.error_fp.close()
            running.remove(d)
            done_n += 1
            man, problem = check_manifest(d.out_bam, d.sampleID)
            results[d.sampleID] = man
            if d.process.returncode != 0 or problem:
                failures.append(problem or f"{d.sampleID}: exit "
                                           f"{d.process.returncode}")
                log(f"[{done_n}/{total}] {d.sampleID} FAILED -- "
                    f"{problem or d.process.returncode} "
                    f"(log: {d.error_file})", "ERROR")
            else:
                if os.path.exists(d.error_file):
                    os.remove(d.error_file)
                sz = os.path.getsize(d.out_bam) / 1e9
                log(f"[{done_n}/{total}] {d.sampleID} ok "
                    f"({man['observed_out']:,} reads, {sz:.2f} GB)")
            if pending:
                launch(pending.pop(0))

    ok = {k: v for k, v in results.items() if v and v.get("status") == "ok"}
    log("=== summary ===")
    log(f"{len(ok)}/{total} slice(s) built")
    if ok:
        sizes = [os.path.getsize(out_for(k)) for k in ok if os.path.exists(out_for(k))]
        if sizes:
            log(f"slice size: min {min(sizes)/1e9:.2f} GB, "
                f"median {sorted(sizes)[len(sizes)//2]/1e9:.2f} GB, "
                f"max {max(sizes)/1e9:.2f} GB")
        pcts = [v.get("mean_depth", 0) for v in ok.values()]
        if pcts:
            log(f"slices average {sum(pcts)/len(pcts):.1f}% of the whole-genome BAM")
    if failures:
        log(f"{len(failures)} failure(s):", "ERROR")
        for f in failures:
            log("  " + f, "ERROR")

    summary = out_dir + "slice_summary.json"
    with open(summary, "w") as fh:
        json.dump({"finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
                   "contig": args.contig, "n_requested": total,
                   "n_ok": len(ok), "failures": failures,
                   "samples": results}, fh, indent=2)
    try:
        fm_obj.uploadData(summary)
    except Exception as e:
        warn(f"summary upload failed: {e}")
    log(f"summary written to {summary}")

    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()