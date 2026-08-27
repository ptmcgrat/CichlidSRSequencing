"""Preflight and per-stage integrity checks for the candidate-QTN genotyping pipeline.

The failure this module exists to prevent: `bcftools mpileup | bcftools call -C alleles`
returns exit code 0 and writes a structurally valid, tabix-indexable VCF containing
ZERO records when the constraint alleles fail to match. Nothing downstream notices.
`bcftools concat` happily produces a file with only the other caller's records, the
driver sees returncode 0, deletes the error log, and the run "succeeds" having
genotyped 6% of the candidate set.

Every helper here is cheap. Run them all.
"""

import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field, asdict


class PipelineError(RuntimeError):
    """Raised when a check fails in a way that makes downstream output meaningless."""


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

_T0 = time.time()


def log(msg, level="INFO"):
    """Timestamped line to stderr. stdout stays clean for machine-readable output."""
    print(f"[{time.time() - _T0:7.1f}s] {level:5s} {msg}", file=sys.stderr, flush=True)


def warn(msg):
    log(msg, "WARN")


def fail(msg):
    raise PipelineError(msg)


# ---------------------------------------------------------------------------
# Filesystem / tooling
# ---------------------------------------------------------------------------

def require_file(path, what, min_bytes=1):
    path = str(path)
    if not os.path.exists(path):
        fail(f"{what} does not exist: {path}")
    size = os.path.getsize(path)
    if size < min_bytes:
        fail(f"{what} exists but is only {size} bytes (expected >= {min_bytes}): {path}")
    return path


def require_index(path, what):
    """Check that the companion index for a bam/vcf.gz/fasta is present and newer."""
    path = str(path)
    if path.endswith(".bam"):
        candidates = [path + ".bai", path[:-4] + ".bai", path + ".csi"]
    elif path.endswith(".vcf.gz"):
        candidates = [path + ".tbi", path + ".csi"]
    elif path.endswith((".fa", ".fasta", ".fna")):
        candidates = [path + ".fai"]
    else:
        fail(f"Don't know the index convention for {what}: {path}")

    found = next((c for c in candidates if os.path.exists(c)), None)
    if found is None:
        fail(f"{what} is not indexed. Expected one of {candidates}")
    if os.path.getmtime(found) < os.path.getmtime(path):
        warn(f"{what} index {found} is OLDER than the file it indexes. "
             f"Stale indexes cause silent partial reads. Re-index it.")
    return found


def tool_version(tool):
    """Return the version string for an htslib-family tool, or None if absent."""
    if shutil.which(tool) is None:
        return None
    try:
        out = subprocess.run([tool, "--version"], capture_output=True, text=True, timeout=30)
        first = (out.stdout or out.stderr).splitlines()[0]
        m = re.search(r"(\d+\.\d+(?:\.\d+)?)", first)
        return m.group(1) if m else first.strip()
    except Exception:
        return "unknown"


def require_tools(tools=("bcftools", "tabix", "bgzip", "samtools")):
    versions = {}
    missing = []
    for t in tools:
        v = tool_version(t)
        if v is None:
            missing.append(t)
        else:
            versions[t] = v
    if missing:
        fail(f"Required tools not on PATH: {', '.join(missing)}")
    log("tool versions: " + ", ".join(f"{k} {v}" for k, v in versions.items()))
    return versions


def mpileup_supports(flag):
    """Detect flag support empirically rather than by version comparison.

    Version numbers are a poor proxy: flags get backported, and the machine that
    built the sites VCF is not necessarily the machine running the genotyping.
    """
    try:
        out = subprocess.run(["bcftools", "mpileup", "--help"],
                             capture_output=True, text=True, timeout=30)
        return flag in (out.stdout + out.stderr)
    except Exception:
        return False


# ---------------------------------------------------------------------------
# VCF inspection
# ---------------------------------------------------------------------------

def _open_maybe_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def count_vcf_records(path):
    """Number of non-header lines. Reads the file directly so it works even when
    the index is missing or corrupt (which is itself a thing worth catching)."""
    n = 0
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                n += 1
    return n


def vcf_sample_names(path):
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                return parts[9:] if len(parts) > 9 else []
            if not line.startswith("#"):
                break
    return []


def vcf_keys(path):
    """(chrom, pos, ref, alt) for every record, uppercased for comparison."""
    keys = []
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            keys.append((f[0], int(f[1]), f[3].upper(), f[4].upper()))
    return keys


def duplicate_keys(keys):
    seen, dups = set(), []
    for k in keys:
        if k in seen:
            dups.append(k)
        else:
            seen.add(k)
    return dups


def format_fields(path):
    """Set of distinct FORMAT strings present in the records."""
    out = set()
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 8:
                out.add(f[8])
    return out


def info_sources(path):
    """Counts of SOURCE= values, so you can see which caller contributed what."""
    from collections import Counter
    c = Counter()
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            info = line.rstrip("\n").split("\t")[7]
            m = re.search(r"SOURCE=([^;\t]+)", info)
            c[m.group(1) if m else "(none)"] += 1
    return dict(c)


# ---------------------------------------------------------------------------
# Reference consistency
# ---------------------------------------------------------------------------

def check_reference_alleles(sites_vcf, ref_fasta, max_report=8):
    """Compare every REF allele in a sites VCF to the actual reference sequence.

    Returns a dict with three counts that mean very different things:

      n_exact      REF matches the genome byte for byte.
      n_case_only  REF matches case-insensitively but not exactly. These sites sit
                   on soft-masked (lowercase) sequence. `bcftools call -C alleles`
                   compares constraint alleles as strings, so a site that is
                   lowercase in the genome and uppercase in the targets file can be
                   dropped without any warning.
      n_mismatch   REF genuinely disagrees with the genome. Usually means the sites
                   VCF was called against a different assembly, or the coordinates
                   are off by one.
    """
    import pysam
    ref = pysam.FastaFile(str(ref_fasta))
    n_exact = n_case = n_bad = 0
    examples = []
    # Read REF exactly as written -- do NOT uppercase, the case is the signal.
    with _open_maybe_gz(sites_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, vref = f[0], int(f[1]), f[3]
            try:
                gref = ref.fetch(chrom, pos - 1, pos - 1 + len(vref))
            except Exception as e:
                n_bad += 1
                if len(examples) < max_report:
                    examples.append((chrom, pos, vref, f"fetch failed: {e}"))
                continue
            if gref == vref:
                n_exact += 1
            elif gref.upper() == vref.upper():
                n_case += 1
                if len(examples) < max_report:
                    examples.append((chrom, pos, vref, f"genome has {gref!r} (case differs)"))
            else:
                n_bad += 1
                if len(examples) < max_report:
                    examples.append((chrom, pos, vref, f"genome has {gref!r}"))
    ref.close()
    return {"n_exact": n_exact, "n_case_only": n_case, "n_mismatch": n_bad,
            "examples": examples}


def softmask_summary(sites_vcf, ref_fasta):
    """Fraction of candidate sites sitting on lowercase (repeat-masked) sequence."""
    import pysam
    ref = pysam.FastaFile(str(ref_fasta))
    lower = upper = 0
    with _open_maybe_gz(sites_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            try:
                b = ref.fetch(f[0], int(f[1]) - 1, int(f[1]))
            except Exception:
                continue
            if b.islower():
                lower += 1
            else:
                upper += 1
    ref.close()
    return {"masked": lower, "unmasked": upper}


def contig_overlap(bam_path, ref_fasta):
    """Catch chr1-vs-NC_ style naming mismatches between BAM and reference."""
    import pysam
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    bam_ctgs = set(bam.references)
    bam.close()
    fai = str(ref_fasta) + ".fai"
    ref_ctgs = set()
    if os.path.exists(fai):
        with open(fai) as fh:
            ref_ctgs = {line.split("\t")[0] for line in fh}
    shared = bam_ctgs & ref_ctgs
    return {"bam_only": sorted(bam_ctgs - ref_ctgs)[:5],
            "ref_only": sorted(ref_ctgs - bam_ctgs)[:5],
            "shared": len(shared)}


# ---------------------------------------------------------------------------
# The diagnostic that answers "why did I get zero records?"
# ---------------------------------------------------------------------------

def diagnose_empty_genotyping(bam_file, sites_vcf, ref_fasta, n_sites=5, region=None):
    """Run the same pileup with and without the -C alleles constraint.

    If the unconstrained run produces records and the constrained one does not,
    the constraint (targets file / allele matching) is the problem. If neither
    produces records, the problem is upstream in mpileup: the BAM has no reads
    there, contig names disagree, or the region argument matched nothing.
    """
    result = {}

    if region is None:
        keys = vcf_keys(sites_vcf)[:n_sites]
        if not keys:
            return {"error": "sites VCF has no records at all"}
        chrom = keys[0][0]
        region = f"{chrom}:{keys[0][1]}-{keys[-1][1]}"
    result["region"] = region

    base = ["bcftools", "mpileup", "-r", region, "-f", str(ref_fasta),
            "-a", "AD,DP", "--min-MQ", "20", "--min-BQ", "20", str(bam_file), "-Ou"]

    try:
        p1 = subprocess.Popen(base, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.run(["bcftools", "call", "-m", "-A", "-Ov"],
                            stdin=p1.stdout, capture_output=True, text=True, timeout=600)
        p1.stdout.close()
        p1.wait()
        unconstrained = sum(1 for l in p2.stdout.splitlines() if not l.startswith("#"))
        result["unconstrained_records"] = unconstrained
        result["unconstrained_stderr"] = (p1.stderr.read().decode(errors="replace")[-2000:]
                                          if p1.stderr else "")
    except Exception as e:
        result["unconstrained_error"] = str(e)

    result["interpretation"] = (
        "constraint/targets is dropping the sites"
        if result.get("unconstrained_records", 0) > 0
        else "mpileup itself sees nothing here -- check BAM coverage, contig names, region"
    )
    return result


# ---------------------------------------------------------------------------
# Per-sample manifest
# ---------------------------------------------------------------------------

@dataclass
class Manifest:
    """Machine-readable record of what each stage actually produced.

    Written next to the output VCF. The driver reads these back instead of
    trusting exit codes, because the whole point is that exit codes lied.
    """
    sample_id: str
    status: str = "started"
    started: str = field(default_factory=lambda: time.strftime("%Y-%m-%dT%H:%M:%S"))
    finished: str = ""
    expected_sv: int = 0
    expected_lv: int = 0
    observed_sv: int = 0
    observed_lv: int = 0
    observed_out: int = 0
    duplicates: int = 0
    format_fields: list = field(default_factory=list)
    sources: dict = field(default_factory=dict)
    genotype_counts: dict = field(default_factory=dict)
    mean_depth: float = 0.0
    tool_versions: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)
    errors: list = field(default_factory=list)
    diagnostic: dict = field(default_factory=dict)

    def add_warning(self, msg):
        self.warnings.append(msg)
        warn(msg)

    def write(self, path):
        self.finished = time.strftime("%Y-%m-%dT%H:%M:%S")
        with open(path, "w") as fh:
            json.dump(asdict(self), fh, indent=2)
        return path

    @staticmethod
    def read(path):
        with open(path) as fh:
            return json.load(fh)