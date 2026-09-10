"""
Genotype small variants (SNVs and small indels) at known sites using bcftools.

Wraps `bcftools mpileup | bcftools call` for force-calling alleles from
a sites VCF, plus normalization and filtering helpers. Requires
bcftools, tabix, and bgzip on PATH.
"""

import subprocess
from typing import Sequence


class BcftoolsError(RuntimeError):
    """Raised when a bcftools subprocess returns non-zero."""


def _run(cmd: Sequence) -> None:
    """Run a command, letting stderr stream to the terminal, raising on failure."""
    result = subprocess.run([str(x) for x in cmd])
    if result.returncode != 0:
        raise BcftoolsError(f"Command failed: {' '.join(str(x) for x in cmd)}")


def _build_targets_tsv(sites_vcf, output_tsv) -> None:
    """Build the CHROM\\tPOS\\tREF,ALT TSV expected by bcftools call -C alleles -T.

    bcftools call with -C alleles does NOT accept a VCF for -T even
    though most other bcftools commands' -T flags do. It requires a
    bgzipped, tabix-indexed TSV with the literal 3-column layout
    CHROM<TAB>POS<TAB>REF,ALT. This helper produces that file from a
    sites VCF.
    """
    query_cmd = ["bcftools", "query",
                 "-f", "%CHROM\\t%POS\\t%REF,%ALT\\n",
                 str(sites_vcf)]

    with open(str(output_tsv), "wb") as out_fh:
        query = subprocess.Popen(query_cmd, stdout=subprocess.PIPE)
        bgzip = subprocess.Popen(["bgzip", "-c"], stdin=query.stdout, stdout=out_fh)
        query.stdout.close()
        bgzip.wait()
        query.wait()

    if query.returncode != 0:
        raise BcftoolsError(f"bcftools query failed (exit {query.returncode})")
    if bgzip.returncode != 0:
        raise BcftoolsError(f"bgzip failed (exit {bgzip.returncode})")

    _run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", str(output_tsv)])


def normalize_sites(sites_vcf, reference_fasta, output_vcf, index: bool = True) -> None:
    """Left-normalize a sites VCF against the reference.

    Accepts either an uncompressed (.vcf) or bgzipped (.vcf.gz) input —
    bcftools auto-detects the format. Output is always bgzipped and (by
    default) tabix-indexed, since the downstream ``genotype_at_sites``
    call needs an indexed VCF for its ``-R`` / ``-T`` arguments.

    bcftools `-C alleles` silently fails to match sites whose alleles
    are not left-normalized, so this step is critical before genotyping.
    """
    _run([
        "bcftools", "norm",
        "-f", reference_fasta,
        "-o", output_vcf,
        "-O", "z",
        sites_vcf,
    ])
    if index:
        _run(["tabix", "-p", "vcf", output_vcf])


def genotype_at_sites(
    bam_file,
    sites_vcf,
    reference_fasta,
    output_vcf,
    *,
    ploidy: int = 2,
    min_mq: int = 20,
    min_bq: int = 20,
    annotations: Sequence[str] = ("AD", "DP", "SP"),
    indels_2: bool = True,
    index: bool = True,
) -> None:
    """Force-call the alleles from sites_vcf at the positions they describe.

    Equivalent to::

        bcftools mpileup -R sites.vcf.gz -f ref.fa -a AD,DP,SP \\
            --indels-2.0 --min-MQ 20 --min-BQ 20 sample.bam \\
        | bcftools call -m -A -C alleles -T sites.vcf.gz \\
            --ploidy 2 -o output.vcf.gz -O z

    Parameters
    ----------
    bam_file : path
        Sample BAM. Must be coord-sorted and indexed (.bai next to it).
    sites_vcf : path
        Sites VCF: bgzipped, tabix-indexed, left-normalized against the
        same reference passed here. Run ``normalize_sites`` first if
        you're not sure.
    reference_fasta : path
        Reference FASTA. Must be faidx-indexed.
    output_vcf : path
        Output bgzipped VCF.
    ploidy : int
        Ploidy for variant calling (default 2).
    min_mq, min_bq : int
        Minimum mapping and base quality thresholds (defaults 20).
    annotations : sequence of str
        FORMAT fields to add. AD is required for the depth-based filter
        applied by ``filter_genotypes``.
    indels_2 : bool
        Use the newer indel-aware model (--indels-2.0). Materially
        better in repetitive regions for small indels.
    index : bool
        Tabix-index the output (default True).
    """
    # bcftools call -C alleles needs CHROM\tPOS\tREF,ALT (NOT a VCF).
    # Build it as a sibling of the output so the user can inspect or
    # reuse it.
    targets_tsv = str(output_vcf) + ".targets.tsv.gz"
    _build_targets_tsv(sites_vcf, targets_tsv)

    mpileup_cmd = [
        "bcftools", "mpileup",
        "-R", sites_vcf,
        "-f", reference_fasta,
        "-a", ",".join(annotations),
        "--min-MQ", min_mq,
        "--min-BQ", min_bq,
        bam_file,
    ]
    if indels_2:
        mpileup_cmd.append("--indels-2.0")

    call_cmd = [
        "bcftools", "call",
        "-m",
        "-A",
        "-C", "alleles",
        "-T", targets_tsv,
        "--ploidy", ploidy,
        "-o", output_vcf,
        "-O", "z",
    ]

    # str() everything for subprocess
    mpileup_cmd = [str(x) for x in mpileup_cmd]
    call_cmd = [str(x) for x in call_cmd]

    mpileup = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    call = subprocess.Popen(call_cmd, stdin=mpileup.stdout,
                            stderr=subprocess.PIPE)
    # Close our handle to mpileup's stdout so mpileup gets SIGPIPE if
    # call exits early; the actual pipe between the two processes is
    # still open.
    mpileup.stdout.close()

    # Drain stderr concurrently so neither process blocks on a full
    # stderr pipe buffer.
    import threading
    mpileup_err: list[bytes] = []
    call_err: list[bytes] = []
    t_m = threading.Thread(target=lambda: mpileup_err.append(mpileup.stderr.read()))
    t_c = threading.Thread(target=lambda: call_err.append(call.stderr.read()))
    t_m.start(); t_c.start()

    call.wait()
    mpileup.wait()
    t_m.join(); t_c.join()

    m_err = mpileup_err[0].decode(errors="replace") if mpileup_err else ""
    c_err = call_err[0].decode(errors="replace") if call_err else ""

    # Check `call` first: when mpileup exits with code 13 it's almost
    # always because call closed the pipe early after its own error.
    if call.returncode != 0:
        raise BcftoolsError(
            f"bcftools call failed (exit {call.returncode})\n"
            f"call stderr:\n{c_err}\n"
            f"mpileup stderr:\n{m_err}"
        )
    if mpileup.returncode != 0:
        raise BcftoolsError(
            f"bcftools mpileup failed (exit {mpileup.returncode})\n"
            f"mpileup stderr:\n{m_err}"
        )

    if index:
        _run(["tabix", "-p", "vcf", output_vcf])


def filter_genotypes(
    input_vcf,
    output_vcf,
    *,
    min_alt_depth: int = 3,
    min_total_depth: int = 10,
    index: bool = True,
) -> None:
    """Keep calls supported by enough reads.

    Filters by alt-allele depth (FORMAT/AD index 1) and total depth.
    QUAL is intentionally not filtered on -- for small indels it's
    noisy and AD/DP are the practical thresholds.
    """
    expr = (f"FORMAT/AD[0:1] >= {min_alt_depth} "
            f"&& FORMAT/DP >= {min_total_depth}")
    _run([
        "bcftools", "view",
        "-i", expr,
        "-o", output_vcf,
        "-O", "z",
        input_vcf,
    ])
    if index:
        _run(["tabix", "-p", "vcf", output_vcf])


if __name__ == "__main__":
    # Input can be uncompressed (e.g. straight from MSAUniqueVariantCaller.write_vcf);
    # normalize_sites bgzips and indexes the output for downstream steps.
    sites      = "sites.vcf"
    sites_norm = "sites.norm.vcf.gz"
    bam        = "sample.bam"
    ref        = "reference.fasta"
    raw        = "genotypes.vcf.gz"
    filtered   = "genotypes.filtered.vcf.gz"

    normalize_sites(sites, ref, sites_norm)
    genotype_at_sites(bam, sites_norm, ref, raw)
    filter_genotypes(raw, filtered)