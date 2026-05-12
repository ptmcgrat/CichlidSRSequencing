"""Classify short Illumina reads (and read pairs) as supporting the reference
or the alternative allele for an indel — either an insertion (e.g. transposon)
or a deletion. Both variant types are handled by the same class via the
standard VCF "anchor base" convention.

Setup
-----
Both alleles include the same flanking buffer sequence on each side; the only
difference between them is the variant content (the inserted bases for an
insertion, or the deleted bases for a deletion):

    insertion:
        ref_allele = left_buffer + anchor + right_buffer
        alt_allele = left_buffer + anchor + INSERTED_SEQ + right_buffer
    deletion:
        ref_allele = left_buffer + anchor + DELETED_SEQ  + right_buffer
        alt_allele = left_buffer + anchor + right_buffer

`left_buffer` and `right_buffer` must be identical between the two alleles.
The buffer is typically 250 bp on each side.

For insertions, the ref allele has one breakpoint (where the inserted
sequence would go) and the alt allele has two (the start and end of the
inserted sequence). For deletions, this flips: ref has two breakpoints
(start and end of the deleted region) and alt has one (where the two flanks
meet). The classifier treats both cases uniformly — it just needs to know
the breakpoint positions in each allele.

A read is informative only if its alignment to the chosen allele crosses
one of those breakpoints with enough sequence on both sides. Reads aligned
inside the variant region without crossing a junction are flagged
uninformative.

Usage
-----
The most common entry point is the VCF-style classmethod, which builds the
local alleles directly from a pysam reference-genome object and a VCF row::

    import pysam
    refObj = pysam.FastaFile("genome.fa")

    # Insertion: ref is the anchor base, alt is anchor + inserted sequence
    classifier = IndelReadClassifier.from_vcf_record(
        ref_genome=refObj, chrom="NC_135176.1", position=24852331,
        ref="g", alt="gcacattttatatttattgtctttttttgctct",
        flanking=250,
    )

    # Deletion: ref is anchor + deleted sequence, alt is the anchor base
    classifier = IndelReadClassifier.from_vcf_record(
        ref_genome=refObj, chrom="NC_135176.1", position=24849727,
        ref="aaaaagacatcactgaattcctcacctgtgctgcc", alt="a",
        flanking=250,
    )

    # Pull reads from your aligned BAM (and optionally a discordant/unmapped BAM)
    reads = classifier.fetch_reads_near_insertion(
        bam="sample.bam",
        discordant_bam="sample.disc.bam",   # optional
    )

    pair_results = classifier.classify_read_pairs(reads)
    call = classifier.call_genotype(pair_results)
    print(call.genotype, call.quality)
"""

from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Iterable
import os
import subprocess
from typing import Iterable, NamedTuple

import parasail


# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

class Support(str, Enum):
    REF = "ref"
    ALT = "alt"
    EQUAL = "equal"
    UNINFORMATIVE = "uninformative"


@dataclass
class ReadResult:
    read_id: str
    support: Support
    score_ref: int
    score_alt: int
    aln_ref_span: tuple[int, int] | None
    aln_alt_span: tuple[int, int] | None
    spans_breakpoint: bool
    note: str = ""


@dataclass
class PairResult:
    pair_id: str                          # base read name (no /1, /2 suffix)
    support: Support
    mate1: ReadResult | None              # None if mate1 wasn't found
    mate2: ReadResult | None              # None if mate2 wasn't found
    note: str = ""


@dataclass
class GenotypeCall:
    genotype: str            # "0/0", "0/1", "1/1", or "./."
    n_ref: int               # ref-supporting pair count
    n_alt: int               # alt-supporting pair count
    n_equal: int             # pairs that aligned equally well to both
    n_uninformative: int     # pairs that didn't carry useful breakpoint info
    log_likelihoods: dict    # {"0/0": float, "0/1": float, "1/1": float}
    quality: float           # phred-scaled gap between best and second-best (GQ)
    note: str = ""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _aligned_span_on_target(result) -> tuple[int, int]:
    """Return (start, end_exclusive) of the matched alignment region on the
    target. parasail's `sw_trace` CIGAR includes leading/trailing 'D' op
    blocks for the unaligned target prefix/suffix; we skip those to find
    the actually-matched region."""
    end = result.end_ref + 1

    cigar_str = result.cigar.decode.decode("ascii")
    ops: list[tuple[int, str]] = []
    num = ""
    for ch in cigar_str:
        if ch.isdigit():
            num += ch
        else:
            ops.append((int(num), ch))
            num = ""

    i = 0
    while i < len(ops) and ops[i][1] == "D":
        i += 1
    j = len(ops)
    while j > i and ops[j - 1][1] == "D":
        j -= 1

    ref_consumed = sum(length for length, op in ops[i:j]
                       if op in ("M", "=", "X", "D"))
    return end - ref_consumed, end


def _crosses_with_overhang(span: tuple[int, int],
                           breakpoints: list[int],
                           overhang: int) -> bool:
    """True iff the alignment span straddles any breakpoint with at least
    `overhang` bp on either side."""
    start, end = span
    return any(start + overhang <= bp <= end - overhang for bp in breakpoints)


def _split_pair_id(read_id: str) -> tuple[str, str | None]:
    """Strip a trailing '/1' or '/2' from a read id; returns (base, mate_tag)."""
    if read_id.endswith("/1"):
        return read_id[:-2], "1"
    if read_id.endswith("/2"):
        return read_id[:-2], "2"
    return read_id, None


def _common_prefix_len(a: str, b: str) -> int:
    """Length of the longest common prefix of two strings."""
    n = 0
    for ca, cb in zip(a, b):
        if ca != cb:
            break
        n += 1
    return n


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_DNAFULL_MATCH_SCORE = 5  # parasail.dnafull match score for ACGT


# ---------------------------------------------------------------------------
# Pair-combination rules
# ---------------------------------------------------------------------------

def _combine_pair(s1: Support, s2: Support) -> tuple[Support, str]:
    """Combine two mate Supports into one pair-level Support, plus a short
    explanatory note. The rules:

    REF + (REF | EQUAL | UNINFORMATIVE) -> REF
    ALT + (ALT | EQUAL | UNINFORMATIVE) -> ALT
    EQUAL + EQUAL -> EQUAL
    REF + ALT -> UNINFORMATIVE  (mates disagree; conservative call)
    anything else with UNINFORMATIVE -> UNINFORMATIVE
    """
    pair = frozenset({s1, s2})
    if Support.REF in pair and Support.ALT in pair:
        return Support.UNINFORMATIVE, "mates disagree (one REF, one ALT)"
    if Support.REF in pair:
        return Support.REF, ""
    if Support.ALT in pair:
        return Support.ALT, ""
    if pair == {Support.EQUAL}:
        return Support.EQUAL, ""
    return Support.UNINFORMATIVE, ""


# ---------------------------------------------------------------------------
# Classifier
# ---------------------------------------------------------------------------

class IndelReadClassifier:
    """Classify reads (and read pairs) against a (ref, alt) pair where alt
    differs from ref by either an insertion (alt longer than ref) or a
    deletion (alt shorter than ref). Both variant types use the same logic:
    the only thing that differs is the breakpoint geometry, which is computed
    automatically from the input.

    Parameters
    ----------
    ref_allele, alt_allele : str
        Local reference and alternative allele sequences. Both must include
        the same flanking buffer on each side; only the variant region (the
        inserted or deleted bases) differs between them.
    contig : str
        Name of the chromosome / contig in the BAM.
    insertion_pos : int
        0-based start position of the variant region on `contig`. For an
        insertion this is the position right after the VCF anchor base. For
        a deletion this is the first deleted base.
    variant_end_pos : int, optional
        0-based half-open end position of the variant region. For deletions,
        this is one past the last deleted base; for insertions and SNVs,
        defaults to `insertion_pos` (single-point variant). Used to define
        the BAM-fetch window so reads spanning either breakpoint of a long
        deletion are captured.

    Breakpoint geometry (one of these two pairs is required)
    --------------------------------------------------------
    buffer_size, transposon_length : int, int
        For the simple "buffer + insertion + buffer" layout (no VCF anchor
        base). Most users should prefer :meth:`from_vcf_record` instead.
    ref_breakpoints, alt_breakpoints : list of int
        Explicit breakpoint positions within ref_allele and alt_allele.
        :meth:`from_vcf_record` populates these for you.

    Scoring knobs (all keyword-only)
    -------------
    gap_open, gap_extend : int  (default 5, 2)
    score_margin : int  (default 20)
    min_score_frac : float  (default 0.6)
    breakpoint_overhang : int  (default 15)
    """

    def __init__(self,
                 ref_allele: str,
                 alt_allele: str,
                 contig: str,
                 insertion_pos: int,
                 *,
                 variant_end_pos: int | None = None,
                 buffer_size: int | None = None,
                 transposon_length: int | None = None,
                 ref_breakpoints: list[int] | None = None,
                 alt_breakpoints: list[int] | None = None,
                 gap_open: int = 5,
                 gap_extend: int = 2,
                 score_margin: int = 20,
                 min_score_frac: float = 0.6,
                 breakpoint_overhang: int = 15) -> None:
        # Two ways to specify breakpoint geometry:
        #   (a) buffer_size + transposon_length  -> simple buffer+TE+buffer
        #   (b) ref_breakpoints + alt_breakpoints -> explicit (used by from_vcf_record)
        if ref_breakpoints is None or alt_breakpoints is None:
            if buffer_size is None or transposon_length is None:
                raise ValueError(
                    "must provide either (buffer_size, transposon_length) or "
                    "(ref_breakpoints, alt_breakpoints)"
                )
            expected_ref_len = 2 * buffer_size
            expected_alt_len = 2 * buffer_size + transposon_length
            if len(ref_allele) != expected_ref_len:
                raise ValueError(
                    f"ref_allele has length {len(ref_allele)}, but expected "
                    f"{expected_ref_len} (= 2 * buffer_size). Make sure the "
                    f"reference allele includes the flanking buffer on both sides."
                )
            if len(alt_allele) != expected_alt_len:
                raise ValueError(
                    f"alt_allele has length {len(alt_allele)}, but expected "
                    f"{expected_alt_len} (= 2 * buffer_size + transposon_length). "
                    f"Make sure the alt allele includes the buffer on both sides."
                )
            ref_breakpoints = [buffer_size]
            alt_breakpoints = [buffer_size, buffer_size + transposon_length]

        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        self.buffer_size = buffer_size
        self.transposon_length = (transposon_length
                                  if transposon_length is not None
                                  else len(alt_allele) - len(ref_allele))
        self.contig = contig
        # Variant region in BAM coordinates (0-based half-open). For
        # insertions and SNVs, end == start (zero-width point). For deletions,
        # end > start.
        self.variant_start = insertion_pos
        self.variant_end = (variant_end_pos
                            if variant_end_pos is not None
                            else insertion_pos)
        # Backwards-compatible alias
        self.insertion_pos = insertion_pos

        self.ref_breakpoints = ref_breakpoints
        self.alt_breakpoints = alt_breakpoints

        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.score_margin = score_margin
        self.min_score_frac = min_score_frac
        self.breakpoint_overhang = breakpoint_overhang

    # -- Alternate constructor: build from a VCF-style record ---------------

    @classmethod
    def from_vcf_record(cls,
                        ref_genome,
                        chrom: str,
                        position: int,
                        ref: str,
                        alt: str,
                        flanking: int = 250,
                        **kwargs) -> "IndelReadClassifier":
        """Build a classifier directly from VCF-style fields and a pysam
        reference-genome object.

        Works for both insertions and deletions, provided they're written in
        the standard VCF "anchor base" convention:

            insertion:  ref="g", alt="gcacatt..."     (ref is just the anchor)
            deletion:   ref="aaaag...", alt="a"       (alt is just the anchor)

        The anchor is the shared prefix of ref and alt — typically a single
        base in the reference at `position` — and the variant content is
        everything after the anchor. Records where ref and alt share NO
        prefix (e.g. alt="-" placeholders) will raise an error; preprocess
        those into the standard form first.

        Parameters
        ----------
        ref_genome : pysam.FastaFile (or anything that supports ``[chrom][a:b]``).
        chrom : str
            Chromosome / contig name (must match the BAM and the FASTA).
        position : int
            1-based position of the first base of `ref` (standard VCF convention).
        ref : str
            VCF REF field. For an insertion this is the anchor base alone; for
            a deletion it's the anchor + the deleted bases.
        alt : str
            VCF ALT field. For an insertion it's the anchor + the inserted
            bases; for a deletion it's the anchor base alone.
        flanking : int, default 250
            Number of flanking bases to include on each side of the variant.
        **kwargs : forwarded to the main constructor.

        Notes
        -----
        Coordinate conversion: VCF is 1-based, pysam is 0-based half-open. So
        ``position - 1`` is the 0-based position of the first base of `ref`.

        For deletions, the BAM-fetch window will cover both breakpoints (the
        start AND end of the deleted region in the reference) plus margin,
        so reads spanning either junction get pulled in.
        """
        if not ref or not alt:
            raise ValueError(
                f"ref and alt must be non-empty strings (got ref={ref!r}, "
                f"alt={alt!r}). For deletions, use VCF anchor convention: "
                f"ref='<anchor><deleted seq>', alt='<anchor>' instead of "
                f"alt='-' or empty."
            )

        # Determine the anchor (common prefix). Standard VCF gives length-1
        # anchors but we accept any non-zero length.
        anchor_len = _common_prefix_len(ref, alt)
        if anchor_len == 0:
            raise ValueError(
                f"ref={ref!r} and alt={alt!r} have no common prefix. "
                f"This is not a valid VCF-anchored indel. For SNVs use a "
                f"different tool; for indels write them as e.g. ref='g', "
                f"alt='gATC' (insertion) or ref='gATC', alt='g' (deletion)."
            )

        # Convert VCF 1-based to pysam 0-based half-open
        ref_start_0 = position - 1                 # first base of ref
        ref_end_0   = ref_start_0 + len(ref)       # one past last base of ref

        # Pull flanking sequence on each side
        left_buffer  = ref_genome[chrom][ref_start_0 - flanking : ref_start_0]
        right_buffer = ref_genome[chrom][ref_end_0 : ref_end_0 + flanking]
        if len(left_buffer) != flanking or len(right_buffer) != flanking:
            raise ValueError(
                f"could not fetch full {flanking}bp flanks at {chrom}:{position} "
                f"(got {len(left_buffer)} / {len(right_buffer)} bp). "
                f"The variant may be too close to the contig edge."
            )

        # Canonicalize to upper case so case mismatches between text file and
        # FASTA don't trip up the aligner
        left_buffer = str(left_buffer).upper()
        right_buffer = str(right_buffer).upper()
        ref = str(ref).upper()
        alt = str(alt).upper()

        ref_allele = left_buffer + ref + right_buffer
        alt_allele = left_buffer + alt + right_buffer

        # Breakpoint geometry. With the anchor convention, the variant region
        # in each allele runs from `anchor_len` to `len(content)`. For an
        # insertion (len(ref) == anchor_len) the ref breakpoint pair collapses
        # to one position; for a deletion (len(alt) == anchor_len) the alt
        # breakpoint pair collapses to one position. Sorted+deduped naturally:
        ref_bps = sorted(set([
            len(left_buffer) + anchor_len,
            len(left_buffer) + len(ref),
        ]))
        alt_bps = sorted(set([
            len(left_buffer) + anchor_len,
            len(left_buffer) + len(alt),
        ]))

        # Genomic coordinates of the variant region (used for the fetch window).
        # For insertions, start == end (zero-width point right after the anchor).
        # For deletions, start..end brackets the deleted bases.
        variant_start = ref_start_0 + anchor_len
        variant_end = ref_end_0

        return cls(
            ref_allele=ref_allele,
            alt_allele=alt_allele,
            contig=chrom,
            insertion_pos=variant_start,
            variant_end_pos=variant_end,
            ref_breakpoints=ref_bps,
            alt_breakpoints=alt_bps,
            **kwargs,
        )

    # -- Internal alignment helper ------------------------------------------

    def _best_alignment(self, read: str, target: str):
        """Align read in both orientations, return the better of the two."""
        fwd = parasail.sw_trace_striped_16(
            read, target, self.gap_open, self.gap_extend, parasail.dnafull
        )
        rev = parasail.sw_trace_striped_16(
            revcomp(read), target, self.gap_open, self.gap_extend, parasail.dnafull
        )
        return fwd if fwd.score >= rev.score else rev

    # -- Per-read classification --------------------------------------------

    def classify_read(self, read_id: str, read_seq: str) -> ReadResult:
        """Classify a single read against the configured alleles."""
        if not read_seq:
            return ReadResult(read_id, Support.UNINFORMATIVE, 0, 0, None, None,
                              spans_breakpoint=False, note="empty sequence")

        aln_ref = self._best_alignment(read_seq, self.ref_allele)
        aln_alt = self._best_alignment(read_seq, self.alt_allele)
        score_ref, score_alt = aln_ref.score, aln_alt.score
        span_ref = _aligned_span_on_target(aln_ref)
        span_alt = _aligned_span_on_target(aln_alt)

        min_score = self.min_score_frac * len(read_seq) * _DNAFULL_MATCH_SCORE
        if max(score_ref, score_alt) < min_score:
            return ReadResult(read_id, Support.UNINFORMATIVE,
                              score_ref, score_alt, span_ref, span_alt,
                              spans_breakpoint=False,
                              note="below min_score on both alleles")

        diff = score_ref - score_alt

        if abs(diff) < self.score_margin:
            spans_either = (
                _crosses_with_overhang(span_ref, self.ref_breakpoints,
                                       self.breakpoint_overhang)
                or _crosses_with_overhang(span_alt, self.alt_breakpoints,
                                          self.breakpoint_overhang)
            )
            return ReadResult(read_id, Support.EQUAL,
                              score_ref, score_alt, span_ref, span_alt,
                              spans_breakpoint=spans_either,
                              note=f"score difference {diff} within margin "
                                   f"{self.score_margin}")

        if diff > 0:
            if _crosses_with_overhang(span_ref, self.ref_breakpoints,
                                      self.breakpoint_overhang):
                return ReadResult(read_id, Support.REF,
                                  score_ref, score_alt, span_ref, span_alt,
                                  spans_breakpoint=True)
            return ReadResult(read_id, Support.UNINFORMATIVE,
                              score_ref, score_alt, span_ref, span_alt,
                              spans_breakpoint=False,
                              note="ref-preferred but no breakpoint crossed "
                                   "(read falls in flanking buffer only)")
        else:
            if _crosses_with_overhang(span_alt, self.alt_breakpoints,
                                      self.breakpoint_overhang):
                return ReadResult(read_id, Support.ALT,
                                  score_ref, score_alt, span_ref, span_alt,
                                  spans_breakpoint=True)
            return ReadResult(read_id, Support.UNINFORMATIVE,
                              score_ref, score_alt, span_ref, span_alt,
                              spans_breakpoint=False,
                              note="alt-preferred but inside transposon body "
                                   "(no junction spanned)")

    def classify_reads(self,
                       reads: Iterable[tuple[str, str]]
                       ) -> list[ReadResult]:
        """Classify many reads. `reads` is an iterable of (read_id, sequence).
        Returns a list of ReadResult in the same order."""
        return [self.classify_read(rid, seq) for rid, seq in reads]

    def call_genotype(self,
                      pair_results: Iterable[PairResult],
                      *,
                      false_alt_rate: float = 0.02,
                      false_ref_rate: float = 0.02,
                      min_informative_pairs: int = 5,
                      ) -> GenotypeCall:
        """Predict a diploid genotype (0/0, 0/1, 1/1, or ./.) from pair-level
        classifications using a binomial likelihood model.

        Likelihood model
        ----------------
        For each candidate genotype, the probability that a randomly chosen
        informative pair appears ALT-supporting is:

            P(alt | 0/0) = false_alt_rate
            P(alt | 0/1) = 0.5
            P(alt | 1/1) = 1 - false_ref_rate

        Given n_alt alt-supporting and n_ref ref-supporting pairs (the
        "informative" pairs), each genotype's likelihood is binomial:

            L(G) = p_alt(G)^n_alt  *  (1 - p_alt(G))^n_ref

        We take the log, pick the genotype with highest log-likelihood, and
        report the gap to the second-best as a phred-scaled quality (GQ).

        Parameters
        ----------
        pair_results : iterable of PairResult
            Output of :meth:`classify_read_pairs`.
        false_alt_rate : float, default 0.02
            Probability that a true REF-genotype pair gets classified ALT.
            This is NOT the per-base sequencing error rate — it's the
            classifier's empirical false-positive rate for ALT calls,
            absorbing mapping artifacts, paralog contamination, and the
            handful of pairs that score borderline against both alleles.
            For real-world tuning, run the classifier on a known REF/REF
            sample and use the observed alt fraction.
        false_ref_rate : float, default 0.02
            Symmetric: probability that a true ALT-genotype pair gets called
            REF (e.g., the breakpoint-spanning evidence got missed).
        min_informative_pairs : int, default 5
            Minimum REF+ALT pair count required to make a call. Below this,
            returns "./." regardless of what the likelihoods say.

        EQUAL and UNINFORMATIVE pairs do not enter the likelihood — they
        confirm there's coverage at the locus but don't distinguish genotypes.
        Their counts are still reported in the output for context.
        """
        import math

        # Tally the four classes
        counts = {Support.REF: 0, Support.ALT: 0,
                  Support.EQUAL: 0, Support.UNINFORMATIVE: 0}
        for pr in pair_results:
            counts[pr.support] += 1
        n_ref = counts[Support.REF]
        n_alt = counts[Support.ALT]
        n_equal = counts[Support.EQUAL]
        n_uninf = counts[Support.UNINFORMATIVE]
        n_inf = n_ref + n_alt

        # Below the call threshold → no-call with empty likelihoods
        if n_inf < min_informative_pairs:
            return GenotypeCall(
                genotype="./.",
                n_ref=n_ref, n_alt=n_alt, n_equal=n_equal, n_uninformative=n_uninf,
                log_likelihoods={"0/0": float("nan"),
                                 "0/1": float("nan"),
                                 "1/1": float("nan")},
                quality=0.0,
                note=f"only {n_inf} informative pairs; need >= "
                     f"{min_informative_pairs}",
            )

        # Compute log-likelihood for each genotype.
        # Using log probabilities directly avoids combinatorial-coefficient
        # tracking — it cancels when we take differences for GQ.
        p_alt_by_g = {
            "0/0": false_alt_rate,
            "0/1": 0.5,
            "1/1": 1.0 - false_ref_rate,
        }
        log_lls: dict[str, float] = {}
        for g, p in p_alt_by_g.items():
            # Guard against log(0) when p == 0 or 1, which only happens if
            # someone sets the error rates to 0 exactly.
            p = min(max(p, 1e-12), 1.0 - 1e-12)
            log_lls[g] = n_alt * math.log(p) + n_ref * math.log(1 - p)

        # Pick the maximum-likelihood genotype
        best_g = max(log_lls, key=log_lls.get)
        # Quality = phred-scaled gap to second best:
        #   log_diff in nats -> bits -> phred via ln(10)/10 conversion
        sorted_lls = sorted(log_lls.values(), reverse=True)
        ln_gap = sorted_lls[0] - sorted_lls[1]
        gq = ln_gap * (10.0 / math.log(10))

        return GenotypeCall(
            genotype=best_g,
            n_ref=n_ref, n_alt=n_alt, n_equal=n_equal, n_uninformative=n_uninf,
            log_likelihoods=log_lls,
            quality=float(gq),
        )

    # -- Per-pair classification --------------------------------------------

    def classify_read_pairs(self,
                            reads: Iterable[tuple[str, str]]
                            ) -> list[PairResult]:
        """Group reads by pair (using /1, /2 suffixes on read IDs) and return
        one PairResult per unique pair.

        Combination rules (mate1, mate2 -> pair):
            REF   + REF/EQUAL/UNINFORMATIVE  -> REF
            ALT   + ALT/EQUAL/UNINFORMATIVE  -> ALT
            EQUAL + EQUAL                    -> EQUAL
            REF   + ALT                      -> UNINFORMATIVE (mates disagree)
            anything else                    -> UNINFORMATIVE

        Singleton reads (mate not found) are classified using just that read.
        """
        # Bucket by base name, preserving first-seen order for stable output
        buckets: dict[str, dict[str, tuple[str, str]]] = {}
        order: list[str] = []
        for rid, seq in reads:
            base, mate = _split_pair_id(rid)
            if base not in buckets:
                buckets[base] = {}
                order.append(base)
            # If duplicate mate tag (or no tag), the second one wins; warn-via-note
            # could be added but for typical BAM input this won't happen.
            buckets[base][mate or "1"] = (rid, seq)

        results: list[PairResult] = []
        for base in order:
            mates = buckets[base]
            r1 = self.classify_read(*mates["1"]) if "1" in mates else None
            r2 = self.classify_read(*mates["2"]) if "2" in mates else None

            if r1 is not None and r2 is not None:
                support, note = _combine_pair(r1.support, r2.support)
            elif r1 is not None:
                support, note = r1.support, "singleton (mate2 missing)"
            elif r2 is not None:
                support, note = r2.support, "singleton (mate1 missing)"
            else:
                support, note = Support.UNINFORMATIVE, "no mates"

            results.append(PairResult(base, support, r1, r2, note))

        return results

    # -- BAM-fetch helper ---------------------------------------------------

    def fetch_reads_near_insertion(self,
                                   bam,
                                   discordant_bam=None,
                                   *,
                                   window: int | None = None,
                                   min_mapq: int = 0
                                   ) -> list[tuple[str, str]]:
        """Pull reads from a BAM around the configured insertion site.

        Parameters
        ----------
        bam : str | Path | pysam.AlignmentFile
            The main aligned BAM (must be coordinate-sorted and indexed).
        discordant_bam : str | Path | pysam.AlignmentFile, optional
            A pre-extracted BAM containing discordant and/or unmapped reads
            (e.g. from `samtools view -f 4` or samblaster's discordants
            output). When provided, unmapped mates are pulled from here much
            faster than scanning the main BAM. If omitted, unmapped mates are
            recovered from the main BAM via pysam's `mate()` (slower).
        window : int, optional
            Half-window around `insertion_pos` to fetch reads from. Defaults
            to `buffer_size + 200` (which catches every read that could span
            a breakpoint plus margin for typical insert sizes).
        min_mapq : int, default 0
            Minimum mapping quality. Reads below this are skipped (but unmapped
            mates retrieved via the discordant BAM are kept regardless of mapq).

        Filters applied (hard-coded based on user preference):
          - Skip duplicates.
          - Skip secondary and supplementary alignments.

        Returns
        -------
        list of (read_id, sequence) tuples, with read_id in 'name/1' or 'name/2'
        form, ready to pass to :meth:`classify_reads` or
        :meth:`classify_read_pairs`. Deduplicated by read_id.
        """
        import pysam   # imported lazily so the module loads without pysam installed

        if window is None:
            # Default fetch window: enough to catch reads anchored anywhere in
            # the local ref allele plus margin for typical insert sizes.
            # The "left buffer" is everything in ref_allele before the first
            # ref breakpoint.
            left_buffer_size = self.ref_breakpoints[0]
            window = left_buffer_size + 200

        win_start = max(0, self.variant_start - window)
        win_end = self.variant_end + window

        # Open BAMs if paths given, remember to close what we opened
        main_bam, main_opened = self._open_bam(bam, pysam)
        disc_bam, disc_opened = (None, False)
        if discordant_bam is not None:
            disc_bam, disc_opened = self._open_bam(discordant_bam, pysam)

        try:
            if not main_bam.has_index():
                raise ValueError(
                    f"main BAM has no index. Run `samtools index <file>` or "
                    f"`pysam.index(...)` first."
                )

            collected: dict[str, str] = {}  # full_read_id -> sequence; dedupes

            # 1. Anchored reads in the window
            for r in main_bam.fetch(self.contig, win_start, win_end):
                if not self._keep_alignment(r, min_mapq):
                    continue
                if r.query_sequence is None:
                    continue   # no sequence stored
                key = self._read_key(r)
                collected[key] = r.query_sequence

            # 2. Unmapped mates whose anchored mate is in the window
            if disc_bam is not None:
                if not disc_bam.has_index():
                    # Allow unindexed discordant BAM — just iterate it
                    iterator = disc_bam.fetch(until_eof=True)
                else:
                    iterator = disc_bam.fetch(until_eof=True)
                for r in iterator:
                    if r.is_duplicate or r.is_secondary or r.is_supplementary:
                        continue
                    # Keep if read OR mate position is in window on contig.
                    # For unmapped reads, BAM convention puts r.reference_start
                    # at the mate's position so reads sort together.
                    own_in_window = (
                        not r.is_unmapped
                        and r.reference_name == self.contig
                        and win_start <= r.reference_start < win_end
                    )
                    mate_in_window = (
                        not r.mate_is_unmapped
                        and r.next_reference_name == self.contig
                        and win_start <= r.next_reference_start < win_end
                    )
                    if not (own_in_window or mate_in_window):
                        continue
                    if r.query_sequence is None:
                        continue
                    if not r.is_unmapped and r.mapping_quality < min_mapq:
                        # Mapped read in the discordant file: still respect mapq
                        continue
                    key = self._read_key(r)
                    if key not in collected:
                        collected[key] = r.query_sequence
            else:
                # Fallback: recover unmapped mates from the main BAM
                # (slower but works without a separate file)
                for r in main_bam.fetch(self.contig, win_start, win_end):
                    if not self._keep_alignment(r, min_mapq):
                        continue
                    if not r.mate_is_unmapped:
                        continue
                    try:
                        mate = main_bam.mate(r)
                    except ValueError:
                        continue   # mate not found
                    if mate.query_sequence is None:
                        continue
                    if mate.is_duplicate or mate.is_secondary or mate.is_supplementary:
                        continue
                    key = self._read_key(mate)
                    if key not in collected:
                        collected[key] = mate.query_sequence

            return list(collected.items())

        finally:
            if main_opened:
                main_bam.close()
            if disc_opened and disc_bam is not None:
                disc_bam.close()

    # -- BAM-fetch helpers (private) ----------------------------------------

    @staticmethod
    def _open_bam(bam_or_path, pysam):
        """Return (AlignmentFile, opened_by_us). Accepts an open file or path."""
        if isinstance(bam_or_path, (str, Path)):
            return pysam.AlignmentFile(str(bam_or_path), "rb"), True
        return bam_or_path, False

    @staticmethod
    def _keep_alignment(r, min_mapq: int) -> bool:
        """Apply the standard quality filters to a primary alignment."""
        if r.is_unmapped:
            return False
        if r.is_duplicate or r.is_secondary or r.is_supplementary:
            return False
        if r.mapping_quality < min_mapq:
            return False
        return True

    @staticmethod
    def _read_key(r) -> str:
        """Return read name with /1 or /2 suffix to disambiguate mates."""
        if r.is_paired:
            mate = "1" if r.is_read1 else "2"
        else:
            mate = "1"
        return f"{r.query_name}/{mate}"


class ClassifierCall(NamedTuple):
    """A per-sample genotype call for one long-indel site."""
    chrom: str
    pos: int          # 1-based, VCF-style, of the anchor base
    ref: str
    alt: str
    gt: str           # "0/0", "0/1", "1/1", or "./." for no-call
    ad_ref: int      # reads supporting the reference allele
    ad_alt: int      # reads supporting the alt allele
 
 
def write_classifier_vcf(
    calls: Iterable[ClassifierCall],
    sample_name: str,
    output_vcf: str,
    template_vcf: str,
) -> None:
    """Write classifier calls to a bgzipped, tabix-indexed VCF.
 
    Parameters
    ----------
    calls : iterable of ClassifierCall
        One per long-indel site for this sample.
    sample_name : str
        Sample column name. MUST match the sample column in the
        short-variant VCF for the same sample, or `bcftools concat`
        will not combine them.
    output_vcf : str
        Output path, must end in ``.vcf.gz``.
    template_vcf : str
        Path to any VCF for the same reference (e.g. the sites VCF or
        the short-variant VCF for this sample). Its ``##contig`` lines
        are copied to the output header so `bcftools concat` and
        `bcftools merge` recognise the same coordinate system.
    """
    if not output_vcf.endswith(".vcf.gz"):
        raise ValueError("output_vcf must end in .vcf.gz")
 
    # Pull ##contig lines from the template so chromosome ordering and
    # lengths match the short-variant VCF.
    header_out = subprocess.run(
        ["bcftools", "view", "-h", template_vcf],
        capture_output=True, text=True, check=True,
    )
    contig_lines = [ln for ln in header_out.stdout.splitlines()
                    if ln.startswith("##contig=")]
 
    # Sort calls by genomic position so the output is sorted.
    calls = sorted(calls, key=lambda c: (c.chrom, c.pos))
 
    raw_path = output_vcf[:-3]  # strip .gz; bgzip will recreate it
    with open(raw_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for cl in contig_lines:
            f.write(cl + "\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,'
                'Description="Allelic depths (ref,alt)">\n')
        f.write('##INFO=<ID=SOURCE,Number=1,Type=String,'
                'Description="Caller of origin">\n')
        f.write("#" + "\t".join([
            "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", sample_name,
        ]) + "\n")
        for c in calls:
            f.write("\t".join([
                c.chrom, str(c.pos), ".",
                c.ref, c.alt, ".", "PASS",
                "SOURCE=IndelReadClassifier",
                "GT:AD",
                f"{c.gt}:{c.ad_ref},{c.ad_alt}",
            ]) + "\n")
 
    # bgzip + tabix
    subprocess.run(["bgzip", "-f", raw_path], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", output_vcf], check=True)
 
 
def concat_sample_vcfs(
    short_vcf: str,
    long_vcf: str,
    output_vcf: str,
) -> None:
    """Combine short-variant + long-indel VCFs for one sample.
 
    Both inputs must use the same sample name in their genotype column
    and the same chromosome naming. The output is sorted, bgzipped,
    and tabix-indexed.
    """
    tmp = output_vcf + ".unsorted.vcf.gz"
    # -a allows records from the two files to overlap (rare but
    # possible if a short and long indel anchor at the same position).
    subprocess.run([
        "bcftools", "concat", "-a", "-Oz", "-o", tmp, short_vcf, long_vcf,
    ], check=True)
    subprocess.run([
        "bcftools", "sort", "-Oz", "-o", output_vcf, tmp,
    ], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", output_vcf], check=True)
    os.remove(tmp)
 