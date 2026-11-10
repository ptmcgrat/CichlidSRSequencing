#!/usr/bin/env python3
"""find_insertion.py - locate an insertion that is absent from the reference (e.g. the chr7 insertion near gsdf) from short reads, and
assemble the two ends of the inserted sequence.

How it works
  An insertion that the reference lacks leaves soft-clipped reads on BOTH sides of its position: reads aligned up to the junction with the
  inserted sequence clipped off (right clips), and reads that start at the second junction with the inserted sequence clipped off the front
  (left clips). If the insertion is flanked by a target-site duplication (TSD) the two boundaries overlap by its length.
  1. For each carrier and control, only the region (default: the 130 kb gsdf window) is streamed from the BAM, exactly as sv-genotype does.
  2. Clip boundaries are counted per sample. Candidates are right/left clip pairs that are supported in the carriers and (nearly) absent
     from every control.
  3. For each candidate the script counts reads that span the site without clipping, so every carrier gets an insertion fraction:
     a homozygous carrier should be near 1, a heterozygous one near 0.5. A variant private to the carriers' lineage (a different kind of
     difference from the reference) does not show that hom/het pattern, so candidates are flagged when the pattern fits.
  4. The reference may hold the region twice (the gsdf unit is duplicated in GT3a, but a carrier with two gsdf copies has one copy per
     chromosome), so a carrier's reads from that one locus map ambiguously to both reference copies and split their evidence. Candidates
     whose inserted ends share sequence are therefore merged and their reads pooled ("seen at N loci").
  5. For the best candidates the clipped tails (insertion start) and clipped heads (insertion end) are collapsed to a consensus and
     extended with unmapped mates by a greedy k-mer walk. The result is written as FASTA. Short reads reach only a few hundred bases
     into the insertion from each end; a 5 kb element needs long reads (or the carrier assembly) to be recovered in full, but the ends,
     the TSD and any terminal inverted repeat are usually enough to BLAST against a transposon library (Dfam / Repbase).

run, in the environment where the `sexloci` command works:
    conda activate ncbi
    python3 find_insertion.py                              # the two Kingiri carriers against 8 benthic females
    python3 find_insertion.py --region NC_135173.1:50500000-51000000        # a wider window if nothing is found
    python3 find_insertion.py --carriers SAMEA117806260:hom SAMEA117806261:het --controls AG_552_1 ASG_160_1 YH_006_f
Outputs (in --outdir): candidates.tsv, report.txt, <rank>_insertion_ends.fasta; region BAMs are cached, rerun with --refetch to replace them.
Needs samtools and rclone on the PATH.
"""
import argparse
import collections
import os
import re
import sys
from bisect import bisect_left, bisect_right
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_CIG = re.compile(r"(\d+)([MIDNSHP=X])")
FLANK = 30                                  # aligned bases kept next to each clip, to tell the insertion from the reference
DEFAULT_REGION = "NC_135173.1:50700000-50830000"
DEFAULT_CARRIERS = ["SAMEA117806260:hom", "SAMEA117806261:het"]
# deep-benthic females with GT3a alignments at 11-27x: chr7-ins has only been seen in A. calliptera and mbuna
DEFAULT_CONTROLS = ["AG_552_1", "AG_558_1", "AG_584_1", "ASG_160_1", "ASG_195_1", "YH_006_f", "YH_010_f", "YH_013_f"]


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


# ------------------------------------------------------------------------------------------------ reads
@dataclass
class Rd:
    name: str
    flag: int
    chrom: str
    pos: int                 # 1-based first aligned base
    mapq: int
    ops: List[Tuple[int, str]]
    mchrom: str
    mpos: int
    tlen: int
    seq: str

    @property
    def aligned_len(self) -> int:
        return sum(n for n, o in self.ops if o in "MDN=X")

    @property
    def end(self) -> int:    # last aligned base
        return self.pos + self.aligned_len - 1


def parse_sam(lines):
    for line in lines:
        if not line or line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        ops = [(int(n), o) for n, o in _CIG.findall(f[5])] if f[5] != "*" else []
        yield Rd(f[0], int(f[1]), f[2], int(f[3]), int(f[4]), ops, f[6], int(f[7]), int(f[8]), f[9])


def clip_events(r: Rd, min_clip: int, min_mapq: int, min_aligned: int):
    """[('R', last aligned base, clipped tail), ('L', first aligned base, clipped head)] for soft clips of at least min_clip bases."""
    if r.flag & (0x4 | 0x100 | 0x400 | 0x800) or r.mapq < min_mapq or not r.ops or r.aligned_len < min_aligned:
        return []
    out = []
    n0, o0 = r.ops[0]
    if o0 == "S" and n0 >= min_clip:
        out.append(("L", r.pos, r.seq[:n0], r.seq[n0:n0 + FLANK]))                    # clipped head, first aligned bases
    n1, o1 = r.ops[-1]
    if o1 == "S" and n1 >= min_clip:
        out.append(("R", r.end, r.seq[-n1:], r.seq[max(len(r.seq) - n1 - FLANK, 0):len(r.seq) - n1]))   # clipped tail, last aligned bases
    return out


@dataclass
class SampleData:
    name: str
    R: Dict[int, List[str]] = field(default_factory=lambda: collections.defaultdict(list))
    L: Dict[int, List[str]] = field(default_factory=lambda: collections.defaultdict(list))
    Rflank: Dict[int, List[str]] = field(default_factory=lambda: collections.defaultdict(list))
    Lflank: Dict[int, List[str]] = field(default_factory=lambda: collections.defaultdict(list))
    spans: List[Tuple[int, int]] = field(default_factory=list)          # (first, last aligned base) of every mapped read
    anomalous: List[int] = field(default_factory=list)                  # positions of reads whose mate is unmapped / elsewhere / far
    unmapped: List[str] = field(default_factory=list)                   # sequences of unmapped reads placed in the region
    n_reads: int = 0

    def finish(self):
        self.spans.sort()
        self.anomalous.sort()
        self._starts = [s for s, _ in self.spans]


def load_sample(name: str, sam_lines, min_clip=20, min_mapq=20, min_aligned=30) -> SampleData:
    d = SampleData(name)
    for r in parse_sam(sam_lines):
        d.n_reads += 1
        if r.flag & 0x4:
            d.unmapped.append(r.seq)
            continue
        if r.flag & (0x100 | 0x400 | 0x800):
            continue
        if r.aligned_len and r.mapq >= min_mapq:
            d.spans.append((r.pos, r.end))
        if (r.flag & 0x1) and ((r.flag & 0x8) or (r.mchrom not in ("=", r.chrom, "*")) or abs(r.tlen) > 3000):
            d.anomalous.append(r.pos)
        for side, pos, seq, flank in clip_events(r, min_clip, min_mapq, min_aligned):
            (d.R if side == "R" else d.L)[pos].append(seq)
            (d.Rflank if side == "R" else d.Lflank)[pos].append(flank)
    d.finish()
    return d


# ------------------------------------------------------------------------------------------------ candidates
def near(counts: Dict[int, List[str]], pos: int, jitter: int) -> List[str]:
    out = []
    for p in range(pos - jitter, pos + jitter + 1):
        out.extend(counts.get(p, ()))
    return out


def cluster_positions(weights: Dict[int, int], jitter: int, min_count: int):
    """Group positions that are within `jitter` bp; return [(representative position (the most clipped), total count)]."""
    out, cur = [], []
    for p in sorted(weights):
        if cur and p - cur[-1] > jitter:
            out.append(cur)
            cur = []
        cur.append(p)
    if cur:
        out.append(cur)
    res = []
    for g in out:
        tot = sum(weights[p] for p in g)
        if tot >= min_count:
            res.append((max(g, key=lambda p: weights[p]), tot))
    return res


def spanning(d: SampleData, p_l: int, p_r: int, margin: int = 10) -> int:
    """Reads aligned across the whole site (with `margin` bases either side) with no clipping at it."""
    lo, hi = min(p_l, p_r) - margin, max(p_l, p_r) + margin
    i = bisect_right(d._starts, lo)
    return sum(1 for s, e in d.spans[:i] if e >= hi)


def anomalous_near(d: SampleData, p: int, w: int = 1000) -> int:
    return bisect_right(d.anomalous, p + w) - bisect_left(d.anomalous, p - w)


def evaluate(members, carriers: Dict[str, SampleData], expected: Dict[str, str], controls: List[SampleData], *, jitter=3, min_side=3,
             max_control=1, allowed_dirty_controls=1) -> dict:
    """One candidate = one or more loci (right-clip boundary, left-clip boundary); several loci are the same insertion seen in different reference
    copies. Evidence is pooled across the loci."""
    sup, fit, per = 0, True, {}
    totR = totL = 0
    for name, d in carriers.items():
        nR = nL = sp = 0
        for p_r, p_l in members:
            nR += len(near(d.R, p_r, jitter)) if p_r is not None else 0
            nL += len(near(d.L, p_l, jitter)) if p_l is not None else 0
            pr = p_r if p_r is not None else p_l - 1
            pl = p_l if p_l is not None else p_r + 1
            sp += spanning(d, pl, pr)
        totR += nR
        totL += nL
        per[name] = (nR, nL, sp)
    two = totR >= min_side and totL >= min_side
    for name, (nR, nL, sp) in per.items():
        ins = (nR + nL) / (2 if two else 1)
        frac = ins / (ins + sp) if ins + sp else 0.0
        per[name] = (nR, nL, sp, frac)
        sup += nR + nL
        exp = expected.get(name, "")
        if exp == "hom" and frac < 0.75:
            fit = False
        if exp == "het" and not 0.2 <= frac <= 0.8:
            fit = False
    dirty, cmax, anom_c = 0, 0, 0
    for c in controls:
        k = sum((len(near(c.R, p_r, jitter)) if p_r is not None else 0) + (len(near(c.L, p_l, jitter)) if p_l is not None else 0) for p_r, p_l in members)
        cmax = max(cmax, k)
        dirty += k > max_control
        anom_c += sum(anomalous_near(c, ((p_r if p_r is not None else p_l - 1) + (p_l if p_l is not None else p_r + 1)) // 2) for p_r, p_l in members)
    mid = lambda p_r, p_l: ((p_r if p_r is not None else p_l - 1) + (p_l if p_l is not None else p_r + 1)) // 2
    anom_by = {name: sum(anomalous_near(d, mid(p_r, p_l)) for p_r, p_l in members) for name, d in carriers.items()}
    anom_ca = sum(anom_by.values()) / max(len(carriers), 1)
    p_r0, p_l0 = members[0]
    tsd = (p_r0 - p_l0 + 1) if (p_r0 is not None and p_l0 is not None) else ""
    row = dict(right_clip_pos=p_r0 if p_r0 is not None else "", left_clip_pos=p_l0 if p_l0 is not None else "", two_sided=two, tsd_len=tsd,
               n_loci=len(members), loci=" + ".join(f"R{a if a is not None else '-'}/L{b if b is not None else '-'}" for a, b in members),
               carrier_reads=sup, hom_het_pattern_fits=fit, controls_with_clips=dirty, max_clips_in_a_control=cmax, controls_clean=dirty <= allowed_dirty_controls,
               anomalous_mates_carriers=round(anom_ca, 1), anomalous_mates_controls=round(anom_c / max(len(controls), 1), 1), low_complexity=False, _members=list(members))
    for name, (nR, nL, sp, frac) in per.items():
        row[f"{name}_right_clips"], row[f"{name}_left_clips"], row[f"{name}_spanning"], row[f"{name}_ins_fraction"] = nR, nL, sp, round(frac, 2)
        row[f"{name}_anomalous_mates"] = anom_by[name]
    return row


def rank_key(r):
    return (not r["low_complexity"], r["controls_clean"], r["two_sided"], r["hom_het_pattern_fits"], r["carrier_reads"])


def find_candidates(carriers: Dict[str, SampleData], expected: Dict[str, str], controls: List[SampleData], *, jitter=3, min_side=3,
                    max_tsd=50, max_gap=10, max_control=1, allowed_dirty_controls=1, allow_one_sided=True):
    wR, wL = collections.Counter(), collections.Counter()
    for d in carriers.values():
        for p, v in d.R.items():
            wR[p] += len(v)
        for p, v in d.L.items():
            wL[p] += len(v)
    Rc = cluster_positions(wR, jitter, min_side)
    Lc = cluster_positions(wL, jitter, min_side)
    cands, used_L = [], set()
    for p_r, n_r in Rc:
        partners = [(p_l, n_l) for p_l, n_l in Lc if p_r - max_tsd <= p_l <= p_r + max_gap]
        if partners:
            p_l, n_l = max(partners, key=lambda x: x[1])
            used_L.add(p_l)
            cands.append((p_r, p_l))
        elif allow_one_sided and n_r >= 2 * min_side:
            cands.append((p_r, None))
    if allow_one_sided:
        for p_l, n_l in Lc:
            if p_l not in used_L and n_l >= 2 * min_side:
                cands.append((None, p_l))
    rows = [evaluate([c], carriers, expected, controls, jitter=jitter, min_side=min_side, max_control=max_control,
                     allowed_dirty_controls=allowed_dirty_controls) for c in cands]
    rows.sort(key=rank_key, reverse=True)
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    return rows


def anomaly_peaks(carriers: Dict[str, SampleData], controls: List[SampleData], expected: Dict[str, str], bin_size=300, step=150, min_carrier=6,
                  min_excess=6, top=10):
    """Region-wide scan that does not need clean clipped reads. Near a transposon insertion, reads have mates that are unmapped or map to other
    chromosomes (the element's interior maps to its other copies). Bins where the carriers have many such reads and no control does mark the site;
    a homozygote should have about twice the heterozygote's count. Returns merged peaks, best first."""
    pos = [p for d in carriers.values() for p in (d.spans[0][0], d.spans[-1][1]) if d.spans]
    if not pos or not controls:
        return []
    lo, hi = min(pos), max(pos)
    bins = []
    for a in range(lo, hi, step):
        b = a + bin_size
        cc = {n: bisect_right(d.anomalous, b) - bisect_left(d.anomalous, a) for n, d in carriers.items()}
        ct = [bisect_right(c.anomalous, b) - bisect_left(c.anomalous, a) for c in controls]
        mean_c = sum(cc.values()) / len(cc)
        excess = mean_c - max(ct)
        if mean_c >= min_carrier and excess >= min_excess:
            bins.append((a, b, cc, max(ct), sum(ct) / len(ct), excess))
    peaks, cur = [], []
    for bn in bins:
        if cur and bn[0] > cur[-1][1]:
            peaks.append(cur)
            cur = []
        cur.append(bn)
    if cur:
        peaks.append(cur)
    out = []
    for g in peaks:
        best = max(g, key=lambda x: x[5])
        hom = [n for n in carriers if expected.get(n) == "hom"]
        het = [n for n in carriers if expected.get(n) == "het"]
        dose = ""
        if hom and het and sum(best[2][n] for n in het):
            dose = round(sum(best[2][n] for n in hom) / sum(best[2][n] for n in het), 2)
        row = dict(start=g[0][0], end=g[-1][1], best_bin_start=best[0], best_bin_end=best[1], controls_max=best[3], controls_mean=round(best[4], 1),
                   excess=round(best[5], 1), hom_over_het=dose)
        for n, v in best[2].items():
            row[f"{n}_anomalous_mates"] = v
        out.append(row)
    out.sort(key=lambda r: r["excess"], reverse=True)
    return out[:top]


def _kset(seq: str, k=20):
    return {min(seq[i:i + k], rc(seq[i:i + k])) for i in range(len(seq) - k + 1)}


def merge_candidates(rows, carriers, expected, controls, args, top=20, min_shared=2, ref_seq=None):
    """Merge candidates (among the best `top`) whose assembled inserted ends share sequence: the same insertion at several loci."""
    head = rows[:top]
    asm = {r["rank"]: assemble(r["_members"], carriers, args.jitter, ref_seq=ref_seq) for r in head}
    ks = {rk: _kset(a["left_ext"]) | _kset(a["right_ext"]) for rk, a in asm.items()}
    parent = {r["rank"]: r["rank"] for r in head}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, r in enumerate(head):
        for r2 in head[i + 1:]:
            if len(ks[r["rank"]] & ks[r2["rank"]]) >= min_shared:
                parent[find(r2["rank"])] = find(r["rank"])
    groups = collections.defaultdict(list)
    for r in head:
        groups[find(r["rank"])].append(r)
    out = []
    for g in groups.values():
        if len(g) == 1:
            r = g[0]
            r["_asm"] = asm[r["rank"]]
        else:
            members = [m for r in sorted(g, key=lambda x: x["rank"]) for m in r["_members"]]
            r = evaluate(members, carriers, expected, controls, jitter=args.jitter, min_side=args.min_side, max_control=args.max_control,
                         allowed_dirty_controls=args.dirty_controls)
            r["_asm"] = assemble(members, carriers, args.jitter, ref_seq=ref_seq)
        out.append(r)
    out.extend(rows[top:])
    return out


# ------------------------------------------------------------------------------------------------ the inserted sequence
def consensus_prefix(seqs: List[str], min_depth=2, min_frac=0.7, max_len=400) -> str:
    out = []
    for i in range(max_len):
        col = [s[i] for s in seqs if len(s) > i]
        if len(col) < min_depth:
            break
        b, c = collections.Counter(col).most_common(1)[0]
        if c / len(col) < min_frac:
            break
        out.append(b)
    return "".join(out)


def consensus_suffix(seqs: List[str], **kw) -> str:
    return consensus_prefix([s[::-1] for s in seqs], **kw)[::-1]


def kmer_table(seqs, k: int) -> Dict[str, int]:
    t = collections.Counter()
    for s in seqs:
        for x in (s, rc(s)):
            for i in range(len(x) - k + 1):
                t[x[i:i + k]] += 1
    return t


def extend_right(seed: str, table: Dict[str, int], k: int, max_ext=3000, min_count=2) -> str:
    cur, seen = seed, set()
    while len(cur) - len(seed) < max_ext:
        pre = cur[-(k - 1):]
        if pre in seen:
            break
        seen.add(pre)
        opts = sorted(((table.get(pre + b, 0), b) for b in "ACGT"), reverse=True)
        (c1, b1), (c2, _) = opts[0], opts[1]
        if c1 < min_count or c1 < 2 * c2:
            break
        cur += b1
    return cur


def join_overlap(a: str, b: str, min_ov=25) -> Optional[str]:
    for ov in range(min(len(a), len(b)), min_ov - 1, -1):
        if a[-ov:] == b[:ov]:
            return a + b[ov:]
    return None


def identity(a: str, b: str) -> float:
    n = min(len(a), len(b))
    return sum(x == y for x, y in zip(a[:n], b[:n])) / n if n else 0.0


def low_complexity(seq: str, frac=0.8, min_len=8) -> bool:
    """Mostly one base, or a short repeat unit (poly-A runs and microsatellites are the usual false junctions)."""
    s = (seq or "").upper()
    if len(s) < min_len:
        return False
    if max(s.count(b) for b in "ACGT") / len(s) >= frac:
        return True
    dis = {s[i:i + 2] for i in range(len(s) - 1)}
    return max(sum(1 for i in range(len(s) - 1) if s[i:i + 2] == d) for d in dis) * 2 / len(s) >= frac


def inverted_repeat(left: str, right: str, max_off=6, window=34, min_len=12):
    """Longest start-of-insertion / end-of-insertion match that reads as an inverted repeat: left[:n] against revcomp(right)[:n], allowing the
    two to be offset by a few bases (clip boundaries move by a base or two). Returns (identity, length, offset) or None."""
    r = rc(right)
    best = None
    for d in range(-max_off, max_off + 1):
        a = left[max(0, -d):max(0, -d) + window]
        b = r[max(0, d):max(0, d) + window]
        for n in range(min(len(a), len(b)), min_len - 1, -1):
            idn = sum(x == y for x, y in zip(a[:n], b[:n])) / n
            if idn >= 0.9 and a[n - 3:n] == b[n - 3:n]:
                if best is None or (n, idn) > (best[1], best[0]):
                    best = (idn, n, d)
                break
    return best


def describe_ends(left: str, right: str, tsd_len, tsd_seq: str) -> List[str]:
    """Plain-language hints about the kind of element; they suggest, they do not identify."""
    hints = []
    if left and right:
        ir = inverted_repeat(left, right)
        if ir:
            hints.append(f"terminal inverted repeat: about {ir[1]} bases of the start match the reverse complement of the end at {ir[0]:.0%} "
                         f"(offset {ir[2]:+d}); a DNA transposon is likely")
        if re.search(r"A{8,}$", right) or re.search(r"^T{8,}", left):
            hints.append("a poly-A / poly-T run at an end (a LINE or SINE-like retrotransposon)")
    if tsd_len not in ("", None) and tsd_len > 0:
        seq = (tsd_seq or "").upper()
        family = ""
        if tsd_len == 2 and seq == "TA":
            family = " (TA: typical of Tc1/mariner)"
        elif tsd_len == 4 and seq == "TTAA":
            family = " (TTAA: typical of piggyBac)"
        elif tsd_len == 8:
            family = " (8 bp: typical of hAT)"
        hints.append(f"a {tsd_len} bp target-site duplication" + (f" {tsd_seq}" if tsd_seq else " (its sequence needs the reference FASTA)") + family)
    elif tsd_len == 0:
        hints.append("no target-site duplication (the two junctions abut)")
    return hints


def trim_into_reference(ext: str, flank: str, seed_len: int, end: str, full=12, partial=6):
    """Cut an extension where it runs into reference. `flank` is the reference next to the junction as seen in the clipped reads: for the 5' end
    (end='right') it is the sequence that follows the insertion, whose first bases appear after the insertion's last base; for the 3' end
    (end='left') it is the sequence that precedes the insertion, whose last bases appear before the insertion's first base. A run of `full` flank
    bases is looked for anywhere in the extension past the seed; a shorter overlap (>= `partial`) only counts right at the extension's tip."""
    if len(flank) < partial:
        return ext, False
    if end == "right":
        i = ext.find(flank[:full], max(seed_len // 2, 1)) if len(flank) >= full else -1
        if i > 0:
            return ext[:i], True
        for m in range(min(full - 1, len(flank)), partial - 1, -1):
            if ext.endswith(flank[:m]) and len(ext) - m > seed_len // 2:
                return ext[:-m], True
    else:
        j = ext.find(flank[-full:], 0, max(len(ext) - seed_len // 2, 1) + full) if len(flank) >= full else -1
        if j >= 0 and j + full < len(ext):
            return ext[j + full:], True
        for m in range(min(full - 1, len(flank)), partial - 1, -1):
            if ext.startswith(flank[-m:]) and len(ext) - m > seed_len // 2:
                return ext[m:], True
    return ext, False


def trim_with_reference(ext: str, window: str, seed_len: int, end: str, k=12) -> Tuple[str, bool]:
    """Cut an extension where it enters the reference. `window` is the reference around the junction(s): reference sequence cannot be part of the
    inserted element, so the first k-mer past the seed that is found in it marks where the extension left the insertion."""
    if len(window) < k or len(ext) < k:
        return ext, False
    ws = {window[i:i + k] for i in range(len(window) - k + 1)}
    if end == "right":                                  # walking away from the 5' junction: look for the first reference k-mer
        for i in range(max(seed_len // 2, 1), len(ext) - k + 1):
            if ext[i:i + k] in ws:
                return ext[:i], True
    else:                                               # 3' end: the extension grew leftwards; look for the last reference k-mer
        for j in range(len(ext) - max(seed_len // 2, 1), k - 1, -1):
            if ext[j - k:j] in ws:
                return ext[j:], True
    return ext, False


def assemble(members, carriers: Dict[str, SampleData], jitter=3, k=25, ref_seq=None):
    """Consensus of the inserted sequence at its two ends, pooled over the loci in `members` [(right-clip pos or None, left-clip pos or None)]."""
    tails = [s for d in carriers.values() for p_r, _ in members if p_r is not None for s in near(d.R, p_r, jitter)]
    heads = [s for d in carriers.values() for _, p_l in members if p_l is not None for s in near(d.L, p_l, jitter)]
    left = consensus_prefix(tails)
    right = consensus_suffix(heads)
    pool = [s for d in carriers.values() for s in d.unmapped] + tails + heads
    table = kmer_table(pool, k)
    left_ext = extend_right(left, table, k) if len(left) >= k else left
    right_ext = rc(extend_right(rc(right), table, k)) if len(right) >= k else right
    # junction reads with a short anchor are left unmapped, so their flank bases enter the pool: cut each extension where it reaches reference
    lflank = consensus_suffix([s for d in carriers.values() for p_r, _ in members if p_r is not None for s in near(d.Rflank, p_r, jitter)], min_depth=3)
    rflank = consensus_prefix([s for d in carriers.values() for _, p_l in members if p_l is not None for s in near(d.Lflank, p_l, jitter)], min_depth=3)
    trimmed = []
    left_ext, cut = trim_into_reference(left_ext, rflank, len(left), end="right")
    if cut:
        trimmed.append("5' end")
    right_ext, cut = trim_into_reference(right_ext, lflank, len(right), end="left")
    if cut:
        trimmed.append("3' end")
    if ref_seq:                                          # the reference itself, when the reads alone were too few to show the flank
        rw = "".join(ref_seq((p_r if p_r is not None else p_l) - 60, (p_r if p_r is not None else p_l) + 120) for p_r, p_l in members)
        lw = "".join(ref_seq((p_l if p_l is not None else p_r) - 120, (p_l if p_l is not None else p_r) + 60) for p_r, p_l in members)
        left_ext, cut = trim_with_reference(left_ext, rw, len(left), "right")
        if cut and "5' end" not in trimmed:
            trimmed.append("5' end")
        right_ext, cut = trim_with_reference(right_ext, lw, len(right), "left")
        if cut and "3' end" not in trimmed:
            trimmed.append("3' end")
        both = [(p_r, p_l) for p_r, p_l in members if p_r is not None and p_l is not None]
        tsd_len = (both[0][0] - both[0][1] + 1) if both else 0          # estimate a missing boundary from the duplication length of a two-sided locus
        for p_r, p_l in members:                         # short overshoots at the very tip: each locus's own reference flanks
            if p_r is None:
                p_r = p_l + tsd_len - 1
            if p_l is None:
                p_l = p_r - tsd_len + 1
            left_ext, cut = trim_into_reference(left_ext, ref_seq(p_l, p_l + 29), len(left), end="right")
            if cut and "5' end" not in trimmed:
                trimmed.append("5' end")
            right_ext, cut = trim_into_reference(right_ext, ref_seq(p_r - 29, p_r), len(right), end="left")
            if cut and "3' end" not in trimmed:
                trimmed.append("3' end")
    joined = join_overlap(left_ext, right_ext) if left_ext and right_ext else None
    return dict(left=left, right=right, left_ext=left_ext, right_ext=right_ext, joined=joined, n_tails=len(tails), n_heads=len(heads),
                n_pool=len(pool), trimmed=trimmed, left_flank=lflank, right_flank=rflank)


# ------------------------------------------------------------------------------------------------ streaming
def fetch_region(sample, base_url, outdir, region, run, refetch=False, retries=2):
    """Region BAM for one sample (cached in outdir/bams). Returns (path or None, message)."""
    from sexloci.gsdf import _classify_failure
    import time
    os.makedirs(os.path.join(outdir, "bams"), exist_ok=True)
    bam = os.path.join(outdir, "bams", f"{sample}.region.bam")
    if os.path.exists(bam) and os.path.getsize(bam) > 0 and not refetch:
        return bam, "cached"
    url = f"{base_url}/{sample}/{sample}.all.bam"
    err = ""
    for attempt in range(retries + 1):
        rc_, _, err = run(["samtools", "view", "-b", "-o", bam, url, region])
        if rc_ == 0:
            return bam, "fetched"
        kind = _classify_failure(err)
        if kind != "error" or attempt == retries:
            return None, kind
        time.sleep(2 * (attempt + 1))
    return None, "error"


def read_bam_text(path, run):
    rc_, out, err = run(["samtools", "view", path])
    if rc_ != 0:
        raise RuntimeError(f"samtools view {path}: {err.strip()[:150]}")
    return out.splitlines()


# ------------------------------------------------------------------------------------------------ reporting
def write_outputs(outdir, rows, carriers, controls, ref_lookup, log, show=5):
    import csv
    os.makedirs(outdir, exist_ok=True)
    clean = [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]
    if clean:
        cols = ["rank"]
        for r in clean:
            cols += [c for c in r if c not in cols]
        with open(os.path.join(outdir, "candidates.tsv"), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(cols)
            for r in clean:
                w.writerow([r.get(c, "") for c in cols])
    lines = [f"carriers: {', '.join(carriers)}", f"controls: {', '.join(c.name for c in controls)}", ""]
    if not rows:
        lines.append("No candidate junctions passed the filters. Try a wider --region, --min-mapq 0, fewer controls, or a lower --min-side.")
    for r in rows[:show]:
        pr, pl = r["right_clip_pos"], r["left_clip_pos"]
        flags = []
        if r["low_complexity"]:
            flags.append("LOW COMPLEXITY (repeat/poly-A: probably not an insertion)")
        if not r["controls_clean"]:
            flags.append(f"clipped reads in {r['controls_with_clips']}/{len(controls)} controls")
        lines.append(f"#{r['rank']}  {'seen at %d loci: ' % r['n_loci'] if r['n_loci'] > 1 else ''}{r['loci']}  "
                     f"{'two-sided' if r['two_sided'] else 'ONE-SIDED'}, carrier reads {r['carrier_reads']}, "
                     f"controls with clips {r['controls_with_clips']}/{len(controls)} (max {r['max_clips_in_a_control']}), "
                     f"hom/het pattern {'fits' if r['hom_het_pattern_fits'] else 'does NOT fit'}")
        if flags:
            lines.append("      warning: " + "; ".join(flags))
        for name in carriers:
            lines.append(f"      {name}: right clips {r[name + '_right_clips']}, left clips {r[name + '_left_clips']}, "
                         f"spanning reads {r[name + '_spanning']}, insertion fraction {r[name + '_ins_fraction']}, anomalous mates {r[name + '_anomalous_mates']}")
        lines.append(f"      anomalous mates within 1 kb: carriers {r['anomalous_mates_carriers']}, controls {r['anomalous_mates_controls']}")
        a = r.get("_asm")
        if a:
            tsd = r["tsd_len"]
            tsd_seq = a["right_flank"][:tsd] if (tsd not in ("", None) and tsd > 0 and len(a["right_flank"]) >= tsd) else ""
            fl = ref_lookup(pl, pr) if ref_lookup and pl != "" and pr != "" else None
            if fl:
                tsd_seq = fl["tsd"] or tsd_seq
                lines.append(f"      reference: ...{fl['left']} | insertion | {fl['right']}...   TSD {fl['tsd'] or '-'}")
            lines.append(f"      insertion start ({len(a['left'])} bp from {a['n_tails']} clipped tails; {len(a['left_ext'])} bp after extension): {a['left_ext'][:90]}"
                         f"{'...' if len(a['left_ext']) > 90 else ''}")
            lines.append(f"      insertion end   ({len(a['right'])} bp from {a['n_heads']} clipped heads; {len(a['right_ext'])} bp after extension): "
                         f"{'...' if len(a['right_ext']) > 90 else ''}{a['right_ext'][-90:]}")
            if a["joined"]:
                lines.append(f"      the two ends overlap: full insertion recovered, {len(a['joined'])} bp")
            for h in describe_ends(a["left_ext"], a["right_ext"], tsd, tsd_seq):
                lines.append(f"      hint: {h}")
            with open(os.path.join(outdir, f"{r['rank']}_insertion_ends.fasta"), "w") as fh:
                loc = r["loci"].replace(" ", "")
                if a["left_ext"]:
                    fh.write(f">insertion_5prime_end loci={loc} length={len(a['left_ext'])}\n{a['left_ext']}\n")
                if a["right_ext"]:
                    fh.write(f">insertion_3prime_end loci={loc} length={len(a['right_ext'])}\n{a['right_ext']}\n")
                if a["joined"]:
                    fh.write(f">insertion_full loci={loc} length={len(a['joined'])}\n{a['joined']}\n")
        lines.append("")
    text = "\n".join(lines)
    with open(os.path.join(outdir, "report.txt"), "w") as fh:
        fh.write(text + "\n")
    log(text)


def analyse(carrier_lines: Dict[str, List[str]], expected: Dict[str, str], control_lines: Dict[str, List[str]], args, ref_lookup=None, log=print, ref_seq=None):
    carriers = {n: load_sample(n, l, args.min_clip, args.min_mapq, args.min_aligned) for n, l in carrier_lines.items()}
    controls = [load_sample(n, l, args.min_clip, args.min_mapq, args.min_aligned) for n, l in control_lines.items()]
    for d in list(carriers.values()) + controls:
        log(f"[scan] {d.name}: {d.n_reads:,} reads, {sum(map(len, d.R.values())):,} right-clipped, {sum(map(len, d.L.values())):,} left-clipped, "
            f"{len(d.unmapped):,} unmapped")
    rows = find_candidates(carriers, expected, controls, jitter=args.jitter, min_side=args.min_side, max_tsd=args.max_tsd, max_gap=args.max_gap,
                           max_control=args.max_control, allowed_dirty_controls=args.dirty_controls)
    rows = merge_candidates(rows, carriers, expected, controls, args, top=args.merge_top, ref_seq=ref_seq)
    for r in rows:                                                       # flag repeats and poly-A using the assembled ends of the best candidates
        a = r.get("_asm")
        if a:
            tsd = r["tsd_len"]
            tsd_seq = a["right_flank"][:tsd] if (tsd not in ("", None) and tsd > 0) else ""
            r["low_complexity"] = bool(low_complexity(tsd_seq) or low_complexity(a["left_ext"]) or low_complexity(a["right_ext"])
                                       or (not a["left_ext"] and not a["right_ext"]))
    rows.sort(key=rank_key, reverse=True)
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    for r in rows[:args.top]:
        r.setdefault("_asm", assemble(r["_members"], carriers, args.jitter, ref_seq=ref_seq))
    write_outputs(args.outdir, rows, list(carriers), controls, ref_lookup, log, show=args.top)
    peaks = anomaly_peaks(carriers, controls, expected)
    if peaks:
        import csv
        with open(os.path.join(args.outdir, "anomalous_mate_peaks.tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(peaks[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(peaks)
    log("regions where the carriers have many reads with unmapped / far / other-chromosome mates and no control does "
        "(hom_over_het is about 2 for a homozygote : heterozygote pair at a real insertion):")
    for r in peaks[:6]:
        log(f"   {r['start']}-{r['end']}  " + "  ".join(f"{k.replace('_anomalous_mates', '')} {v}" for k, v in r.items() if k.endswith("_anomalous_mates"))
            + f"  | controls max {r['controls_max']} mean {r['controls_mean']}  | hom/het {r['hom_over_het']}")
    if not peaks:
        log("   none")
    return rows, {r["rank"]: r.get("_asm") for r in rows[:args.top]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--carriers", nargs="+", default=DEFAULT_CARRIERS, help="SampleID:hom|het (default: the two Kingiri carriers)")
    ap.add_argument("--controls", nargs="+", default=DEFAULT_CONTROLS, help="SampleIDs of fish without the insertion that have GT3a alignments")
    ap.add_argument("--region", default=DEFAULT_REGION)
    ap.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    ap.add_argument("--outdir", default=os.path.expanduser("~/Temp/SexDeterminationAnalysis/insertion_scan"))
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--refetch", action="store_true")
    ap.add_argument("--reference-fasta", help="GT3a FASTA, to print the target-site duplication and flanks (default: found in the genome folder)")
    ap.add_argument("--min-clip", type=int, default=20, help="soft clip length that counts (default 20)")
    ap.add_argument("--min-mapq", type=int, default=20, help="mapping quality for clipped and spanning reads; use 0 in a duplicated region (default 20)")
    ap.add_argument("--min-aligned", type=int, default=30)
    ap.add_argument("--jitter", type=int, default=3, help="boundaries within this many bp are the same junction (default 3)")
    ap.add_argument("--min-side", type=int, default=3, help="clipped reads needed on each side, summed over carriers (default 3)")
    ap.add_argument("--max-tsd", type=int, default=50)
    ap.add_argument("--max-gap", type=int, default=10, help="bases by which the left-clip boundary may lie beyond the right-clip boundary (default 10)")
    ap.add_argument("--merge-top", type=int, default=20, help="candidates compared for shared inserted sequence (the same insertion at several loci)")
    ap.add_argument("--max-control", type=int, default=1, help="clipped reads a control may have at a junction (default 1)")
    ap.add_argument("--dirty-controls", type=int, default=1, help="controls allowed to break that rule (default 1)")
    ap.add_argument("--top", type=int, default=3, help="candidates to assemble (default 3)")
    a = ap.parse_args()
    expected = {}
    for c in a.carriers:
        name, _, gt = c.partition(":")
        expected[name] = gt
    from sexloci.config import Config
    from sexloci.gsdf import HttpServer, default_runner
    from sexloci.labdata import LabData
    lab = LabData.from_config(Config())
    log = lambda m: print(m, flush=True)
    chrom = a.region.split(":")[0]
    ref_lookup = ref_seq = None
    try:
        from sexloci.svgeno import find_reference_fasta, load_reference_windows
        lo, hi = (int(x) for x in a.region.split(":")[1].replace(",", "").split("-"))
        fasta = find_reference_fasta(lab.genome_dir(a.genome_version), a.reference_fasta)
        win = load_reference_windows(fasta, {chrom: (max(lo - 1, 0), hi)})
        ref_lookup = lambda pl, pr: dict(left=win[chrom][pr - 40:pr], right=win[chrom][pl - 1:pl + 39], tsd=win[chrom][pl - 1:pr] if pr >= pl else "")
        ref_seq = lambda a, b: win[chrom][max(a, 1) - 1:b]
    except Exception as e:
        log(f"[scan] reference FASTA not used ({str(e)[:90]})")
    samples = list(expected) + [c for c in a.controls if c not in expected]
    os.makedirs(a.outdir, exist_ok=True)
    paths = {}
    with HttpServer(lab.cloud_master_dir + f"Bamfiles/{a.genome_version}/", None) as base:
        with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as ex:
            futs = {ex.submit(fetch_region, s, base, a.outdir, a.region, default_runner, a.refetch): s for s in samples}
            for f in as_completed(futs):
                s = futs[f]
                path, msg = f.result()
                paths[s] = path
                log(f"[fetch] {s}: {msg}")
    missing = [s for s in expected if paths.get(s) is None]
    if missing:
        sys.exit(f"cannot continue without the carrier BAM(s): {', '.join(missing)}")
    ctrl = {s: read_bam_text(paths[s], default_runner) for s in a.controls if paths.get(s)}
    if len(ctrl) < 3:
        sys.exit("fewer than 3 controls could be read; give more with --controls")
    car = {s: read_bam_text(paths[s], default_runner) for s in expected}
    analyse(car, expected, ctrl, a, ref_lookup, log, ref_seq)
    log(f"\nwrote candidates.tsv, report.txt and the FASTA files in {a.outdir}")


if __name__ == "__main__":
    main()
