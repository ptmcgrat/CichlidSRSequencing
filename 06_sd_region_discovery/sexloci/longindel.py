"""Classify short Illumina reads (and read pairs) as supporting the reference or the alternative allele of a long indel (e.g. a
transposon insertion), and call a diploid genotype from them.

This is the lab's IndelReadClassifier, kept as it was written (same scoring, same breakpoint-overhang rule, same pair-combination rules,
same binomial genotype model). Two things changed so that it can run inside sexloci:

  * the Smith-Waterman calls sit behind a small `backend` object. `ParasailBackend` is the original code (parasail's striped SIMD
    alignment); `NumpyBackend` is a slower pure-numpy equivalent used when parasail is not installed and by the tests;
  * reads can come from `samtools view` text (`reads_from_samtools_lines`), so no pysam is needed. The pysam-based
    `fetch_reads_near_insertion` is still there for notebook use.

Setup, as before: both alleles carry the same flanking buffer on each side and differ only by the variant content,

    insertion:  ref = left + anchor + right        alt = left + anchor + INSERTED + right
    deletion:   ref = left + anchor + DELETED + right    alt = left + anchor + right

A read is informative only if its alignment to the allele it prefers crosses a breakpoint with enough sequence on both sides.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np


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
    aln_ref_span: Optional[Tuple[int, int]]
    aln_alt_span: Optional[Tuple[int, int]]
    spans_breakpoint: bool
    note: str = ""


@dataclass
class PairResult:
    pair_id: str                          # base read name (no /1, /2 suffix)
    support: Support
    mate1: Optional[ReadResult]           # None if mate1 wasn't found
    mate2: Optional[ReadResult]           # None if mate2 wasn't found
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


def call_from_counts(n_ref: int, n_alt: int, *, false_alt_rate: float = 0.02, false_ref_rate: float = 0.02,
                     min_informative_pairs: int = 5) -> Tuple[str, float, dict]:
    """The diploid genotype model from counts of informative read pairs: a binomial with P(alt | 0/0) = false_alt_rate, P(alt | 0/1) = 0.5 and
    P(alt | 1/1) = 1 - false_ref_rate. Returns (genotype, GQ, log-likelihoods): GQ is the phred-scaled gap between the best and second-best
    genotype, and fewer than `min_informative_pairs` REF + ALT pairs gives ('./.', 0.0, NaNs). IndelReadClassifier.call_genotype uses this, and
    so do the combined calls over several sites (svgeno.add_group_columns), which keeps the two exactly alike."""
    import math
    n_inf = n_ref + n_alt
    if n_inf < min_informative_pairs:
        return "./.", 0.0, {"0/0": float("nan"), "0/1": float("nan"), "1/1": float("nan")}
    p_alt_by_g = {"0/0": false_alt_rate, "0/1": 0.5, "1/1": 1.0 - false_ref_rate}
    log_lls = {}
    for g, p in p_alt_by_g.items():
        p = min(max(p, 1e-12), 1.0 - 1e-12)
        log_lls[g] = n_alt * math.log(p) + n_ref * math.log(1 - p)
    best_g = max(log_lls, key=log_lls.get)
    sorted_lls = sorted(log_lls.values(), reverse=True)
    gq = (sorted_lls[0] - sorted_lls[1]) * (10.0 / math.log(10))
    return best_g, float(gq), log_lls


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _aligned_span_on_target(result) -> Tuple[int, int]:
    """Return (start, end_exclusive) of the matched alignment region on the target. parasail's `sw_trace` CIGAR includes leading/trailing
    'D' op blocks for the unaligned target prefix/suffix; we skip those to find the actually-matched region."""
    end = result.end_ref + 1

    cigar_str = result.cigar.decode.decode("ascii")
    ops: List[Tuple[int, str]] = []
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

    ref_consumed = sum(length for length, op in ops[i:j] if op in ("M", "=", "X", "D"))
    return end - ref_consumed, end


def _crosses_with_overhang(span: Tuple[int, int], breakpoints: List[int], overhang: int) -> bool:
    """True iff the alignment span straddles any breakpoint with at least `overhang` bp on either side."""
    start, end = span
    return any(start + overhang <= bp <= end - overhang for bp in breakpoints)


def _split_pair_id(read_id: str) -> Tuple[str, Optional[str]]:
    """Strip a trailing '/1' or '/2' from a read id; returns (base, mate_tag)."""
    if read_id.endswith("/1"):
        return read_id[:-2], "1"
    if read_id.endswith("/2"):
        return read_id[:-2], "2"
    return read_id, None


def _common_prefix_len(a: str, b: str) -> int:
    n = 0
    for ca, cb in zip(a, b):
        if ca != cb:
            break
        n += 1
    return n


_DNAFULL_MATCH_SCORE = 5  # parasail.dnafull match score for ACGT


def _combine_pair(s1: Support, s2: Support) -> Tuple[Support, str]:
    """Combine two mate Supports into one pair-level Support:
    REF + (REF | EQUAL | UNINFORMATIVE) -> REF;  ALT + (ALT | EQUAL | UNINFORMATIVE) -> ALT;  EQUAL + EQUAL -> EQUAL;
    REF + ALT -> UNINFORMATIVE (mates disagree); anything else with UNINFORMATIVE -> UNINFORMATIVE."""
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
# Alignment backends
# ---------------------------------------------------------------------------

class ParasailBackend:
    """The original alignment calls (parasail striped 16-bit Smith-Waterman, DNAFULL scores)."""
    name = "parasail"

    def __init__(self, gap_open: int, gap_extend: int):
        import parasail                                   # raises ImportError if not installed
        self._p, self.go, self.ge = parasail, gap_open, gap_extend

    def score(self, read: str, target: str) -> Tuple[int, bool]:
        """Best score over both orientations, without traceback. Returns (score, used_revcomp)."""
        p = self._p
        fwd = p.sw_striped_16(read, target, self.go, self.ge, p.dnafull)
        rev = p.sw_striped_16(revcomp(read), target, self.go, self.ge, p.dnafull)
        return (fwd.score, False) if fwd.score >= rev.score else (rev.score, True)

    def span(self, read: str, target: str, use_revcomp: bool) -> Tuple[int, int]:
        """(start, end_exclusive) of the best local alignment on the target, one known orientation."""
        p = self._p
        query = revcomp(read) if use_revcomp else read
        res = p.sw_trace_striped_16(query, target, self.go, self.ge, p.dnafull)
        return _aligned_span_on_target(res)


class NumpyBackend:
    """Affine-gap Smith-Waterman (match +5, mismatch -4, like DNAFULL for ACGT) in numpy: rows of the read are processed one at a time,
    vectorised over the target. Slower than parasail on long targets but needs nothing beyond numpy."""
    name = "numpy"
    MATCH, MISMATCH = 5, -4
    NEG = -(10 ** 6)

    def __init__(self, gap_open: int, gap_extend: int):
        self.go, self.ge = gap_open, gap_extend

    @staticmethod
    def _arr(s: str) -> np.ndarray:
        return np.frombuffer(s.upper().encode("ascii"), np.uint8)

    def _sub_row(self, base: int, t: np.ndarray) -> np.ndarray:
        s = np.where(t == base, self.MATCH, self.MISMATCH).astype(np.int32)
        if base == 78:                                         # N never matches
            s[:] = self.MISMATCH
        s[t == 78] = self.MISMATCH
        return s

    def _dp(self, read: str, target: str, keep: bool):
        r, t = self._arr(read), self._arr(target)
        m, n = len(r), len(t)
        go, ge, NEG = self.go, self.ge, self.NEG
        idx = np.arange(n, dtype=np.int64)
        H_prev = np.zeros(n + 1, np.int32)
        F_prev = np.full(n + 1, NEG, np.int32)
        best, bi, bj = 0, 0, 0
        if keep:
            Hm = np.zeros((m + 1, n + 1), np.int32)
            Em = np.full((m + 1, n + 1), NEG, np.int32)
            Fm = np.full((m + 1, n + 1), NEG, np.int32)
        for i in range(1, m + 1):
            s = self._sub_row(int(r[i - 1]), t)
            F_row = np.maximum(F_prev - ge, H_prev - go)                 # gap in the target (consumes read bases)
            Hpre = np.maximum(np.maximum(H_prev[:-1] + s, F_row[1:]), 0)  # columns 1..n
            acc = np.maximum.accumulate(Hpre.astype(np.int64) + idx * ge)
            E = np.full(n, NEG, np.int64)                                # gap in the read (consumes target bases)
            if n > 1:
                E[1:] = acc[:-1] - go - (idx[1:] - 1) * ge
            Hrow = np.maximum(Hpre, E).astype(np.int32)
            H_cur = np.concatenate(([0], Hrow)).astype(np.int32)
            j = int(Hrow.argmax())
            if Hrow[j] > best:
                best, bi, bj = int(Hrow[j]), i, j + 1
            if keep:
                Hm[i], Fm[i] = H_cur, F_row
                Em[i, 1:] = E.astype(np.int32)
            H_prev, F_prev = H_cur, F_row.astype(np.int32)
        if keep:
            return best, bi, bj, Hm, Em, Fm, r, t
        return best, bi, bj

    def score(self, read: str, target: str) -> Tuple[int, bool]:
        fwd = self._dp(read, target, False)[0]
        rev = self._dp(revcomp(read), target, False)[0]
        return (fwd, False) if fwd >= rev else (rev, True)

    def span(self, read: str, target: str, use_revcomp: bool) -> Tuple[int, int]:
        q = revcomp(read) if use_revcomp else read
        best, i, j, H, E, F, r, t = self._dp(q, target, True)
        if best <= 0:
            return 0, 0
        end = j
        go, ge = self.go, self.ge
        state = "H"
        while i > 0 and j > 0:
            if state == "H":
                h = int(H[i, j])
                if h == 0:
                    break
                s = self.MATCH if (r[i - 1] == t[j - 1] and r[i - 1] != 78) else self.MISMATCH
                if h == int(H[i - 1, j - 1]) + s:
                    i, j = i - 1, j - 1
                elif h == int(E[i, j]):
                    state = "E"
                elif h == int(F[i, j]):
                    state = "F"
                else:
                    break
            elif state == "E":                                   # gap in the read: moves along the target
                opened = int(E[i, j]) == int(H[i, j - 1]) - go
                j -= 1
                if opened:
                    state = "H"
            else:                                                # gap in the target: moves along the read
                opened = int(F[i, j]) == int(H[i - 1, j]) - go
                i -= 1
                if opened:
                    state = "H"
        return j, end


def make_backend(name: str, gap_open: int, gap_extend: int):
    """name: 'auto' (parasail if installed, else numpy), 'parasail' or 'numpy'."""
    if name in ("auto", "parasail"):
        try:
            return ParasailBackend(gap_open, gap_extend)
        except ImportError:
            if name == "parasail":
                raise ImportError("the parasail backend was requested but parasail is not installed (pip install parasail)") from None
    return NumpyBackend(gap_open, gap_extend)


# ---------------------------------------------------------------------------
# Reads from `samtools view` text
# ---------------------------------------------------------------------------

def reads_from_samtools_lines(lines: Iterable[str], min_mapq: int = 0, include_unmapped: bool = True) -> List[Tuple[str, str]]:
    """(read_id, sequence) from SAM text lines, the way fetch_reads_near_insertion collects them: duplicates, secondary and supplementary
    alignments are skipped (pass `-F 3328` to samtools to do that for you), unmapped reads are kept (they carry the inserted sequence and
    sit at their mate's position, so a region query returns them), ids are name/1 or name/2, deduplicated."""
    out, seen = [], set()
    for line in lines:
        if not line or line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        flag, mapq, seq = int(f[1]), int(f[4]), f[9]
        if flag & (0x400 | 0x100 | 0x800) or seq == "*":
            continue
        if flag & 0x4:
            if not include_unmapped:
                continue
        elif mapq < min_mapq:
            continue
        mate = ("1" if flag & 0x40 else "2") if flag & 0x1 else "1"
        key = f"{f[0]}/{mate}"
        if key not in seen:
            seen.add(key)
            out.append((key, seq))
    return out


# ---------------------------------------------------------------------------
# Classifier
# ---------------------------------------------------------------------------

class IndelReadClassifier:
    """Classify reads (and read pairs) against a (ref, alt) pair where alt differs from ref by an insertion or a deletion. See the module
    docstring; parameters are as in the lab's original class, plus `backend` ('auto', 'parasail' or 'numpy', or a backend object).

    Scoring knobs: gap_open=5, gap_extend=2, score_margin=20, min_score_frac=0.6, breakpoint_overhang=15."""

    def __init__(self, ref_allele: str, alt_allele: str, contig: str, insertion_pos: int, *,
                 variant_end_pos: Optional[int] = None, buffer_size: Optional[int] = None, transposon_length: Optional[int] = None,
                 ref_breakpoints: Optional[List[int]] = None, alt_breakpoints: Optional[List[int]] = None,
                 gap_open: int = 5, gap_extend: int = 2, score_margin: int = 20, min_score_frac: float = 0.6,
                 breakpoint_overhang: int = 15, backend="auto") -> None:
        if ref_breakpoints is None or alt_breakpoints is None:
            if buffer_size is None or transposon_length is None:
                raise ValueError("must provide either (buffer_size, transposon_length) or (ref_breakpoints, alt_breakpoints)")
            expected_ref_len = 2 * buffer_size
            expected_alt_len = 2 * buffer_size + transposon_length
            if len(ref_allele) != expected_ref_len:
                raise ValueError(f"ref_allele has length {len(ref_allele)}, but expected {expected_ref_len} (= 2 * buffer_size). "
                                 f"Make sure the reference allele includes the flanking buffer on both sides.")
            if len(alt_allele) != expected_alt_len:
                raise ValueError(f"alt_allele has length {len(alt_allele)}, but expected {expected_alt_len} "
                                 f"(= 2 * buffer_size + transposon_length). Make sure the alt allele includes the buffer on both sides.")
            ref_breakpoints = [buffer_size]
            alt_breakpoints = [buffer_size, buffer_size + transposon_length]

        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        self.buffer_size = buffer_size
        self.transposon_length = transposon_length if transposon_length is not None else len(alt_allele) - len(ref_allele)
        self.contig = contig
        self.variant_start = insertion_pos
        self.variant_end = variant_end_pos if variant_end_pos is not None else insertion_pos
        self.insertion_pos = insertion_pos
        self.ref_breakpoints = ref_breakpoints
        self.alt_breakpoints = alt_breakpoints
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.score_margin = score_margin
        self.min_score_frac = min_score_frac
        self.breakpoint_overhang = breakpoint_overhang
        self.backend = make_backend(backend, gap_open, gap_extend) if isinstance(backend, str) else backend

    # -- Alternate constructor: build from a VCF-style record ---------------

    @classmethod
    def from_vcf_record(cls, ref_genome, chrom: str, position: int, ref: str, alt: str, flanking: int = 250, **kwargs) -> "IndelReadClassifier":
        """Build a classifier from VCF-style fields (1-based position, anchor-base convention) and a reference object that supports
        ``ref_genome[chrom][a:b]`` (0-based half-open), e.g. a pysam.FastaFile or sexloci.svgeno.RefWindows.

            insertion:  ref="g", alt="gcacatt..."     deletion:  ref="aaaag...", alt="a"
        """
        if not ref or not alt:
            raise ValueError(f"ref and alt must be non-empty strings (got ref={ref!r}, alt={alt!r}). For deletions, use VCF anchor "
                             f"convention: ref='<anchor><deleted seq>', alt='<anchor>' instead of alt='-' or empty.")
        anchor_len = _common_prefix_len(ref.upper(), alt.upper())
        if anchor_len == 0:
            raise ValueError(f"ref={ref!r} and alt={alt!r} have no common prefix. This is not a valid VCF-anchored indel. "
                             f"Write indels as ref='g', alt='gATC' (insertion) or ref='gATC', alt='g' (deletion).")
        ref_start_0 = position - 1
        ref_end_0 = ref_start_0 + len(ref)
        left_buffer = ref_genome[chrom][ref_start_0 - flanking: ref_start_0]
        right_buffer = ref_genome[chrom][ref_end_0: ref_end_0 + flanking]
        if len(left_buffer) != flanking or len(right_buffer) != flanking:
            raise ValueError(f"could not fetch full {flanking}bp flanks at {chrom}:{position} (got {len(left_buffer)} / {len(right_buffer)} bp). "
                             f"The variant may be too close to the contig edge.")
        left_buffer = str(left_buffer).upper()
        right_buffer = str(right_buffer).upper()
        ref = str(ref).upper()
        alt = str(alt).upper()
        ref_allele = left_buffer + ref + right_buffer
        alt_allele = left_buffer + alt + right_buffer
        ref_bps = sorted(set([len(left_buffer) + anchor_len, len(left_buffer) + len(ref)]))
        alt_bps = sorted(set([len(left_buffer) + anchor_len, len(left_buffer) + len(alt)]))
        return cls(ref_allele=ref_allele, alt_allele=alt_allele, contig=chrom, insertion_pos=ref_start_0 + anchor_len,
                   variant_end_pos=ref_end_0, ref_breakpoints=ref_bps, alt_breakpoints=alt_bps, **kwargs)

    # -- Per-read classification --------------------------------------------

    def classify_read(self, read_id: str, read_seq: str) -> ReadResult:
        """Classify a single read against the configured alleles."""
        if not read_seq:
            return ReadResult(read_id, Support.UNINFORMATIVE, 0, 0, None, None, spans_breakpoint=False, note="empty sequence")

        score_ref, rev_ref = self.backend.score(read_seq, self.ref_allele)
        score_alt, rev_alt = self.backend.score(read_seq, self.alt_allele)

        min_score = self.min_score_frac * len(read_seq) * _DNAFULL_MATCH_SCORE
        if max(score_ref, score_alt) < min_score:
            return ReadResult(read_id, Support.UNINFORMATIVE, score_ref, score_alt, None, None,
                              spans_breakpoint=False, note="below min_score on both alleles")

        diff = score_ref - score_alt

        if abs(diff) < self.score_margin:
            span_ref = self.backend.span(read_seq, self.ref_allele, rev_ref)
            span_alt = self.backend.span(read_seq, self.alt_allele, rev_alt)
        elif diff > 0:
            span_ref = self.backend.span(read_seq, self.ref_allele, rev_ref)
            span_alt = None
        else:
            span_ref = None
            span_alt = self.backend.span(read_seq, self.alt_allele, rev_alt)

        if abs(diff) < self.score_margin:
            spans_either = (_crosses_with_overhang(span_ref, self.ref_breakpoints, self.breakpoint_overhang)
                            or _crosses_with_overhang(span_alt, self.alt_breakpoints, self.breakpoint_overhang))
            return ReadResult(read_id, Support.EQUAL, score_ref, score_alt, span_ref, span_alt, spans_breakpoint=spans_either,
                              note=f"score difference {diff} within margin {self.score_margin}")

        if diff > 0:
            if _crosses_with_overhang(span_ref, self.ref_breakpoints, self.breakpoint_overhang):
                return ReadResult(read_id, Support.REF, score_ref, score_alt, span_ref, span_alt, spans_breakpoint=True)
            return ReadResult(read_id, Support.UNINFORMATIVE, score_ref, score_alt, span_ref, span_alt, spans_breakpoint=False,
                              note="ref-preferred but no breakpoint crossed (read falls in flanking buffer only)")
        if _crosses_with_overhang(span_alt, self.alt_breakpoints, self.breakpoint_overhang):
            return ReadResult(read_id, Support.ALT, score_ref, score_alt, span_ref, span_alt, spans_breakpoint=True)
        return ReadResult(read_id, Support.UNINFORMATIVE, score_ref, score_alt, span_ref, span_alt, spans_breakpoint=False,
                          note="alt-preferred but inside transposon body (no junction spanned)")

    def classify_reads(self, reads: Iterable[Tuple[str, str]]) -> List[ReadResult]:
        return [self.classify_read(rid, seq) for rid, seq in reads]

    def call_genotype(self, pair_results: Iterable[PairResult], *, false_alt_rate: float = 0.02, false_ref_rate: float = 0.02,
                      min_informative_pairs: int = 5) -> GenotypeCall:
        """Diploid genotype from pair-level classifications with a binomial model: P(alt | 0/0) = false_alt_rate, P(alt | 0/1) = 0.5,
        P(alt | 1/1) = 1 - false_ref_rate; EQUAL and UNINFORMATIVE pairs do not enter the likelihood. Quality = phred-scaled gap between
        the best and second-best genotype. Fewer than `min_informative_pairs` REF+ALT pairs gives './.'."""
        counts = {Support.REF: 0, Support.ALT: 0, Support.EQUAL: 0, Support.UNINFORMATIVE: 0}
        for pr in pair_results:
            counts[pr.support] += 1
        n_ref, n_alt = counts[Support.REF], counts[Support.ALT]
        n_equal, n_uninf = counts[Support.EQUAL], counts[Support.UNINFORMATIVE]
        gt, gq, log_lls = call_from_counts(n_ref, n_alt, false_alt_rate=false_alt_rate, false_ref_rate=false_ref_rate,
                                           min_informative_pairs=min_informative_pairs)
        note = f"only {n_ref + n_alt} informative pairs; need >= {min_informative_pairs}" if gt == "./." else ""
        return GenotypeCall(genotype=gt, n_ref=n_ref, n_alt=n_alt, n_equal=n_equal, n_uninformative=n_uninf, log_likelihoods=log_lls, quality=gq, note=note)

    # -- Per-pair classification --------------------------------------------

    def classify_read_pairs(self, reads: Iterable[Tuple[str, str]]) -> List[PairResult]:
        """Group reads by pair (/1, /2 suffixes) and return one PairResult per pair (rules in `_combine_pair`). Singleton reads are
        classified using just that read."""
        buckets: dict = {}
        order: List[str] = []
        for rid, seq in reads:
            base, mate = _split_pair_id(rid)
            if base not in buckets:
                buckets[base] = {}
                order.append(base)
            buckets[base][mate or "1"] = (rid, seq)
        results: List[PairResult] = []
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

    # -- Where to fetch reads -----------------------------------------------

    def fetch_window(self, window: Optional[int] = None) -> Tuple[int, int]:
        """(start0, end0) of the genomic window the lab's fetch uses: the left buffer plus 200 bp on each side of the variant."""
        if window is None:
            window = self.ref_breakpoints[0] + 200
        return max(0, self.variant_start - window), self.variant_end + window

    # -- pysam fetch (unchanged, for notebook use; sexloci streams with samtools instead) ---------------

    def fetch_reads_near_insertion(self, bam, discordant_bam=None, *, window: Optional[int] = None, min_mapq: int = 0,
                                   disc_region_fetch: bool = False, include_unmapped: bool = False) -> List[Tuple[str, str]]:
        """Pull reads from an indexed BAM around the variant with pysam. See the lab's original docstring: duplicates, secondary and
        supplementary alignments are skipped; with include_unmapped the unmapped mates stored at the anchored mate's position are kept."""
        import pysam   # imported lazily so the module loads without pysam installed
        win_start, win_end = self.fetch_window(window)
        main_bam, main_opened = self._open_bam(bam, pysam)
        disc_bam, disc_opened = (None, False)
        if discordant_bam is not None:
            disc_bam, disc_opened = self._open_bam(discordant_bam, pysam)
        try:
            if not main_bam.has_index():
                raise ValueError("main BAM has no index. Run `samtools index <file>` or `pysam.index(...)` first.")
            collected: dict = {}
            for r in main_bam.fetch(self.contig, win_start, win_end):
                if r.is_duplicate or r.is_secondary or r.is_supplementary:
                    continue
                if r.is_unmapped:
                    if not include_unmapped:
                        continue
                elif r.mapping_quality < min_mapq:
                    continue
                if r.query_sequence is None:
                    continue
                collected[self._read_key(r)] = r.query_sequence
            if disc_bam is not None:
                iterator = None
                if disc_region_fetch:
                    try:
                        if disc_bam.has_index():
                            iterator = disc_bam.fetch(self.contig, win_start, win_end)
                    except (ValueError, AttributeError):
                        iterator = None
                if iterator is None:
                    iterator = disc_bam.fetch(until_eof=True)
                for r in iterator:
                    if r.is_duplicate or r.is_secondary or r.is_supplementary:
                        continue
                    own_in_window = (not r.is_unmapped and r.reference_name == self.contig and win_start <= r.reference_start < win_end)
                    mate_in_window = (not r.mate_is_unmapped and r.next_reference_name == self.contig
                                      and win_start <= r.next_reference_start < win_end)
                    if not (own_in_window or mate_in_window) or r.query_sequence is None:
                        continue
                    if not r.is_unmapped and r.mapping_quality < min_mapq:
                        continue
                    key = self._read_key(r)
                    if key not in collected:
                        collected[key] = r.query_sequence
            return list(collected.items())
        finally:
            if main_opened:
                main_bam.close()
            if disc_opened and disc_bam is not None:
                disc_bam.close()

    @staticmethod
    def _open_bam(bam_or_path, pysam):
        if isinstance(bam_or_path, (str, Path)):
            return pysam.AlignmentFile(str(bam_or_path), "rb"), True
        return bam_or_path, False

    @staticmethod
    def _read_key(r) -> str:
        mate = ("1" if r.is_read1 else "2") if r.is_paired else "1"
        return f"{r.query_name}/{mate}"
