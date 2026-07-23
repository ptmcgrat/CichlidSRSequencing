#!/usr/bin/env python3
"""
find_diagnostic_indels.py

Find ~200 bp indel length-polymorphism markers between two cichlid species, in
STANDARD-REFERENCE coordinates, for PCR-based genotyping of recombinants in an
F1 cross.

Strategy
--------
1. Align species A and species B *each* to the standard reference (minimap2),
   call indels (paftools.js), and left-align/normalize them (bcftools norm).
   Normalization is load-bearing: the "same" indel called from two independent
   alignments must be represented identically or set operations won't see it as
   shared.
2. Intersect the two call sets (bcftools isec). An indel PRIVATE to one species
   is a candidate marker (that species carries it; the other is reference-like
   there). A SHARED indel is NOT a marker -- both species differ from the
   reference identically, so there is no A-vs-B length difference to score.
3. Keep private indels in the PCR-resolvable size window (default 150-300 bp),
   restricted to the region of interest.
4. Guard: the size you see on the gel is the A-vs-B difference, which equals the
   A-vs-ref size ONLY when the other species truly matches the reference at that
   locus. Drop any candidate where the other species has an indel nearby, since
   its predicted band size would be wrong.

Outputs (in --outdir)
---------------------
  diagnostic_indels.vcf   sites-only VCF, reference coords, INFO: CARRIER/GELDIFF/VARTYPE
  diagnostic_indels.tsv   eyeball-friendly marker table
  ambiguous_dropped.tsv   candidates dropped by the proximity guard (rescue by hand if wanted)
  align/                  per-species normalized, indexed call sets
  isec/                   bcftools isec output (0000=private A, 0001=private B, 0002/3=shared)
  run.log                 full stderr from every tool invocation

Requires on PATH: minimap2, paftools.js (ships with minimap2; needs the k8 engine),
                  samtools, bcftools, sort

Example
-------
  ./find_diagnostic_indels.py \
      --ref cichlid_ref.fa \
      --species-a speciesA.fa --species-b speciesB.fa \
      --region chr7:12000000-13500000 \
      --outdir markers_chr7/
"""

import argparse
import gzip
import logging
import shutil
import subprocess
import sys
from pathlib import Path

REQUIRED_TOOLS = ["minimap2", "paftools.js", "samtools", "bcftools", "sort"]


# --------------------------------------------------------------------------- #
# subprocess helpers
# --------------------------------------------------------------------------- #
def check_tools(tools):
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        sys.exit("ERROR: required tools not found on PATH: " + ", ".join(missing))


def run(cmd, log_fh):
    """Run a single command; stderr -> log file; raise on non-zero exit."""
    cmd = [str(c) for c in cmd]
    logging.info("RUN: %s", " ".join(cmd))
    log_fh.write("\n# RUN: " + " ".join(cmd) + "\n")
    log_fh.flush()
    proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=log_fh)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed (exit {proc.returncode}): {' '.join(cmd)} -- see run.log")


def run_pipeline(stages, log_fh):
    """Run a list of commands as a shell pipeline via chained Popen.

    Each stage's stdout feeds the next stage's stdin. The final stage is
    expected to write its own output file (e.g. `bcftools norm -o ...`), so its
    stdout is discarded. Every stage's stderr goes to the log file, which also
    avoids the classic PIPE-buffer deadlock on chatty tools like minimap2.
    """
    stages = [[str(c) for c in s] for s in stages]
    logging.info("PIPE: %s", " | ".join(" ".join(s) for s in stages))
    log_fh.write("\n# PIPE: " + " | ".join(" ".join(s) for s in stages) + "\n")
    log_fh.flush()

    procs = []
    prev = None
    for i, cmd in enumerate(stages):
        last = i == len(stages) - 1
        proc = subprocess.Popen(
            cmd,
            stdin=(prev.stdout if prev is not None else None),
            stdout=(subprocess.DEVNULL if last else subprocess.PIPE),
            stderr=log_fh,
        )
        if prev is not None:
            prev.stdout.close()  # let upstream receive SIGPIPE if downstream dies
        procs.append(proc)
        prev = proc

    rets = [p.wait() for p in procs]
    for cmd, ret in zip(stages, rets):
        if ret != 0:
            raise RuntimeError(f"pipeline stage failed (exit {ret}): {' '.join(cmd)} -- see run.log")


# --------------------------------------------------------------------------- #
# pipeline stages
# --------------------------------------------------------------------------- #
def align_and_call(label, species_fa, ref, align_dir, args, log_fh):
    """minimap2 -> sort -> paftools call (VCF) -> bcftools norm; index. Returns .vcf.gz path."""
    out_gz = align_dir / f"{label}.norm.vcf.gz"
    stages = [
        # -c: emit CIGAR in PAF; --cs: base-level diff string that paftools parses
        ["minimap2", "-cx", args.preset, "--cs", "-t", args.threads, ref, species_fa],
        # paftools requires the PAF sorted by target name (col 6), target start (col 8)
        ["sort", "-k6,6", "-k8,8n"],
        # -f enables VCF output; -L/-l lowered so a small region's alignment isn't skipped
        ["paftools.js", "call", "-s", label,
         "-L", args.min_aln_var, "-l", args.min_aln_cov, "-f", ref, "-"],
        # left-align + normalize so shared indels match across the two call sets
        ["bcftools", "norm", "-f", ref, "-Oz", "-o", out_gz, "-"],
    ]
    run_pipeline(stages, log_fh)
    run(["bcftools", "index", "-f", out_gz], log_fh)
    return out_gz


def run_isec(a_gz, b_gz, isec_dir, log_fh):
    """Partition into private/shared. -c some tolerates near-identical representations."""
    run(["bcftools", "isec", "-c", "some", "-p", isec_dir, a_gz, b_gz], log_fh)
    # 0000 = records private to first input (A); 0001 = private to second (B)
    return isec_dir / "0000.vcf", isec_dir / "0001.vcf"


# --------------------------------------------------------------------------- #
# VCF parsing / filtering (done in Python for clarity)
# --------------------------------------------------------------------------- #
def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_indel_spans(vcf_path, chrom=None):
    """Return {chrom: sorted[(start, end)]} for every indel (for the proximity guard)."""
    spans = {}
    with _open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            c, pos, ref, alt = f[0], int(f[1]), f[3], f[4]
            if chrom is not None and c != chrom:
                continue
            if len(ref) == len(alt):  # SNV/MNV, not an indel
                continue
            spans.setdefault(c, []).append((pos, pos + len(ref) - 1))
    for c in spans:
        spans[c].sort()
    return spans


def overlaps_any(pos, end, intervals, window):
    """True if [pos-window, end+window] overlaps any (start, end) in the sorted list."""
    lo, hi = pos - window, end + window
    for s, e in intervals:
        if s > hi:
            break
        if e >= lo:
            return True
    return False


def collect_candidates(private_vcf, carrier, other_spans, region, args):
    """Return (kept, dropped) records. Each record: (chrom,pos,ref,alt,info,carrier,size,vartype)."""
    chrom_f, start_f, end_f = region
    kept, dropped = [], []
    with _open(private_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            c, pos, ref, alt, info = f[0], int(f[1]), f[3], f[4], f[7]
            if len(ref) == len(alt):
                continue
            size = abs(len(ref) - len(alt))
            if not (args.min_size <= size <= args.max_size):
                continue
            end = pos + len(ref) - 1
            if chrom_f is not None and (c != chrom_f or end < start_f or pos > end_f):
                continue
            vartype = "DEL" if len(ref) > len(alt) else "INS"
            rec = (c, pos, ref, alt, info, carrier, size, vartype)
            if overlaps_any(pos, end, other_spans.get(c, []), args.guard_window):
                dropped.append(rec)  # other species has a nearby indel -> band size unreliable
            else:
                kept.append(rec)
    return kept, dropped


# --------------------------------------------------------------------------- #
# output
# --------------------------------------------------------------------------- #
def write_vcf(records, ref, out_vcf):
    contigs = []
    with open(str(ref) + ".fai") as fh:
        for line in fh:
            name, length = line.split("\t")[:2]
            contigs.append((name, length))
    with open(out_vcf, "w") as v:
        v.write("##fileformat=VCFv4.2\n")
        for name, length in contigs:
            v.write(f"##contig=<ID={name},length={length}>\n")
        v.write('##INFO=<ID=CARRIER,Number=1,Type=String,Description="Species carrying the indel relative to reference">\n')
        v.write('##INFO=<ID=GELDIFF,Number=1,Type=Integer,Description="Predicted A-vs-B amplicon size difference (bp); equals indel size because the other species is reference-like here">\n')
        v.write('##INFO=<ID=VARTYPE,Number=1,Type=String,Description="DEL or INS relative to reference">\n')
        v.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, pos, ref_a, alt, info, carrier, size, vartype in sorted(records, key=lambda r: (r[0], r[1])):
            prefix = "" if info in (".", "") else info + ";"
            newinfo = f"{prefix}CARRIER={carrier};GELDIFF={size};VARTYPE={vartype}"
            v.write(f"{c}\t{pos}\t.\t{ref_a}\t{alt}\t.\t.\t{newinfo}\n")


def write_tsv(records, out_tsv):
    with open(out_tsv, "w") as t:
        t.write("chrom\tpos\tcarrier\tvartype\tref_len\talt_len\tgel_diff_bp\n")
        for c, pos, ref_a, alt, info, carrier, size, vartype in sorted(records, key=lambda r: (r[0], r[1])):
            t.write(f"{c}\t{pos}\t{carrier}\t{vartype}\t{len(ref_a)}\t{len(alt)}\t{size}\n")


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def parse_region(s):
    if not s:
        return (None, None, None)
    chrom, rest = s.split(":")
    a, b = rest.replace(",", "").split("-")
    return (chrom, int(a), int(b))


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", required=True, help="standard reference FASTA")
    p.add_argument("--species-a", required=True, help="species A assembly FASTA")
    p.add_argument("--species-b", required=True, help="species B assembly FASTA")
    p.add_argument("--region", default=None, help="restrict to region, e.g. chr7:12000000-13500000")
    p.add_argument("--outdir", default="diagnostic_indels_out", help="output directory")
    p.add_argument("--label-a", default="A", help="name for species A in outputs")
    p.add_argument("--label-b", default="B", help="name for species B in outputs")
    p.add_argument("--preset", default="asm5", help="minimap2 asm preset (asm5/asm10/asm20 by divergence)")
    p.add_argument("--threads", type=int, default=4, help="minimap2 threads")
    p.add_argument("--min-size", type=int, default=150, help="min indel size (bp) to keep")
    p.add_argument("--max-size", type=int, default=300, help="max indel size (bp) to keep")
    p.add_argument("--guard-window", type=int, default=50,
                   help="drop a candidate if the other species has an indel within this many bp")
    p.add_argument("--min-aln-var", type=int, default=1000, help="paftools -L: min alignment length to call variants")
    p.add_argument("--min-aln-cov", type=int, default=1000, help="paftools -l: min alignment length for coverage")
    return p.parse_args()


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    align_dir = outdir / "align"
    isec_dir = outdir / "isec"
    for d in (outdir, align_dir, isec_dir):
        d.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s", datefmt="%H:%M:%S")
    check_tools(REQUIRED_TOOLS)

    ref = Path(args.ref)
    region = parse_region(args.region)

    with open(outdir / "run.log", "w", buffering=1) as log_fh:
        if not Path(str(ref) + ".fai").exists():
            run(["samtools", "faidx", ref], log_fh)

        logging.info("Aligning + calling: %s", args.label_a)
        a_gz = align_and_call(args.label_a, Path(args.species_a), ref, align_dir, args, log_fh)
        logging.info("Aligning + calling: %s", args.label_b)
        b_gz = align_and_call(args.label_b, Path(args.species_b), ref, align_dir, args, log_fh)

        logging.info("Intersecting call sets")
        priv_a, priv_b = run_isec(a_gz, b_gz, isec_dir, log_fh)

        # For the guard, load the OTHER species' indels (restricted to region chrom if given)
        chrom_f = region[0]
        logging.info("Loading indel spans for proximity guard")
        b_spans = load_indel_spans(b_gz, chrom=chrom_f)
        a_spans = load_indel_spans(a_gz, chrom=chrom_f)

        logging.info("Filtering private candidates (%d-%d bp)", args.min_size, args.max_size)
        kept_a, drop_a = collect_candidates(priv_a, args.label_a, b_spans, region, args)
        kept_b, drop_b = collect_candidates(priv_b, args.label_b, a_spans, region, args)

    kept = kept_a + kept_b
    dropped = drop_a + drop_b

    write_vcf(kept, ref, outdir / "diagnostic_indels.vcf")
    write_tsv(kept, outdir / "diagnostic_indels.tsv")
    write_tsv(dropped, outdir / "ambiguous_dropped.tsv")

    logging.info("Done.")
    logging.info("  %s-private markers kept: %d", args.label_a, len(kept_a))
    logging.info("  %s-private markers kept: %d", args.label_b, len(kept_b))
    logging.info("  dropped by proximity guard: %d (see ambiguous_dropped.tsv)", len(dropped))
    logging.info("  -> %s", outdir / "diagnostic_indels.tsv")


if __name__ == "__main__":
    main()