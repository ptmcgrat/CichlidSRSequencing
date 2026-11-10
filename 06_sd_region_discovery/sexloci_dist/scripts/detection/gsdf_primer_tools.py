#!/usr/bin/env python3
"""gsdf_primer_tools.py - primers to genotype the gsdf tandem duplication, and an in-silico PCR check on any assembly.

Needs only numpy. Two subcommands.

1) design: pick primers for three assays from the sequence of the locus (GT3a coordinates, NC_135173.1)

     python3 gsdf_primer_tools.py design --fasta GT3a_genome.fna --out gsdf     # from the reference sequence
     python3 gsdf_primer_tools.py design --bam MC_050_gsdf.bam --out gsdf       # or from the reads of one fish (consensus; also
                                                                                # marks positions where the reads disagree)
   Assays (the unit is the 21,728 bp segment that is repeated head-to-tail; J = where copy 1 ends and copy 2 begins):
     J  junction PCR: a forward primer at the END of the unit and a reverse primer at the START of the unit. They face each other only
        across J, so they amplify only chromosomes that carry the tandem duplication. On a single-copy chromosome they point away from
        each other and give nothing. Presence/absence of the duplication (cannot tell one duplicated chromosome from two).
     FL single-copy control in the left flank: amplifies every sample (DNA/PCR control; also the normaliser for the qPCR below).
     U  copy-number qPCR / ddPCR amplicon inside the unit (both copies identical there): amplifies every copy.
        U / FL = (total copies)/2 : 1.0 for 2 copies, 1.5 for 3, 2.0 for 4. This is the assay that tells -/-, dup/- and dup/dup apart.

2) check: where would the primers amplify in an assembly? (in-silico PCR, up to 3 mismatches, exact 3' end)

     python3 gsdf_primer_tools.py check gsdf_primers.tsv assembly1.fa assembly2.fa.gz ...

   Predicts the junction product in assemblies that carry the duplication, none in those that do not, and warns where a primer site
   differs from the assembly (e.g. other species).
"""
import argparse
import gzip
import math
import os
import re
import struct
import sys

import numpy as np

CHROM, LO, HI = "NC_135173.1", 50_700_000, 50_830_000
J_DEFAULT = 50_773_050          # junction between copy 1 and copy 2 (GT3a), from the paired-read fragment-coverage minimum
PERIOD = 21_728
COPY1_START = J_DEFAULT - PERIOD
COPY2_END = COPY1_START + 2 * PERIOD - 1
FL_WIN = (50_736_000, 50_748_000)     # unique left flank, outside the unit
U_WIN = (50_757_000, 50_770_000)      # inside the unit, middle of copy 1 (identical in copy 2)
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def rc(s):
    return s.translate(COMP)[::-1]


# ------------------------------------------------------------------------------------------------ sequence input
def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_region_fasta(path, chrom, start, end):
    parts, on, pos = [], False, 0
    with open_text(path) as fh:
        for line in fh:
            if line[0] == ">":
                if on:
                    break
                on = line[1:].split()[0] == chrom
                continue
            if on:
                line = line.rstrip("\n")
                n = len(line)
                lo, hi = max(start - 1 - pos, 0), min(end - pos, n)
                if hi > lo:
                    parts.append(line[lo:hi])
                pos += n
                if pos >= end:
                    break
    if not parts:
        sys.exit(f"{chrom}:{start}-{end} not found in {path}")
    s = "".join(parts).upper()
    return s, np.zeros(len(s), bool)


def consensus_from_bam(path, lo, hi, min_depth=3):
    """Consensus of the reads of one fish over lo..hi (N where depth < min_depth) and a mask of positions where the reads disagree."""
    data = gzip.open(path, "rb").read()
    o = 4
    (lt,) = struct.unpack_from("<i", data, o)
    o += 4 + lt
    (nr,) = struct.unpack_from("<i", data, o)
    o += 4
    for _ in range(nr):
        (ln,) = struct.unpack_from("<i", data, o)
        o += 4 + ln + 4
    cig_ops, bases = "MIDNSHP=X", "=ACMGRSVTWYHKDBN"
    n = hi - lo + 1
    cnt = np.zeros((n, 4), int)
    idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    while o < len(data):
        (bs,) = struct.unpack_from("<i", data, o)
        _ref, pos, lrn, _mq, _bin, ncig, flag, lseq, _nr, _np, _tl = struct.unpack_from("<iiBBHHHiiii", data, o + 4)
        p = o + 4 + 32 + lrn
        if not flag & 4:
            cig = [(struct.unpack_from("<I", data, p + 4 * i)[0] >> 4, cig_ops[struct.unpack_from("<I", data, p + 4 * i)[0] & 15]) for i in range(ncig)]
            q0 = p + 4 * ncig
            seq = "".join(bases[(data[q0 + i // 2] >> (4 if i % 2 == 0 else 0)) & 15] for i in range(lseq))
            q, r = 0, pos + 1
            for L, op in cig:
                if op in "M=X":
                    for i in range(L):
                        b, x = idx.get(seq[q + i]), r + i - lo
                        if b is not None and 0 <= x < n:
                            cnt[x, b] += 1
                    q += L
                    r += L
                elif op in "IS":
                    q += L
                elif op in "DN":
                    r += L
        o += 4 + bs
    depth = cnt.sum(1)
    order = np.sort(cnt, 1)
    cons = "".join("ACGT"[m] if d >= min_depth else "N" for m, d in zip(cnt.argmax(1), depth))
    variable = (order[:, 2] >= 2) & (order[:, 2] >= 0.2 * np.maximum(depth, 1))
    return cons, variable


# ------------------------------------------------------------------------------------------------ primer properties
NN = {"AA": (-7.9, -22.2), "TT": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3), "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
      "GT": (-8.4, -22.4), "AC": (-8.4, -22.4), "CT": (-7.8, -21.0), "AG": (-7.8, -21.0), "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
      "CG": (-10.6, -27.2), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9), "CC": (-8.0, -19.9)}


def tm_nn(s, na=0.05, ct=2.5e-7):
    """SantaLucia 1998 nearest-neighbour Tm (50 mM Na+, 250 nM primer). A guide: confirm with your oligo vendor's calculator."""
    dh = ds = 0.0
    for a, b in zip(s, s[1:]):
        h, e = NN[a + b]
        dh += h
        ds += e
    for end in (s[0], s[-1]):
        dh, ds = (dh + 0.1, ds - 2.8) if end in "GC" else (dh + 2.3, ds + 4.1)
    ds += 0.368 * (len(s) - 1) * math.log(na)
    return dh * 1000 / (ds + 1.987 * math.log(ct / 4)) - 273.15


def gc(s):
    return 100 * sum(c in "GC" for c in s) / len(s)


def has_run(s, n=4):
    return any(len(set(s[i:i + n])) == 1 for i in range(len(s) - n + 1))


def has_dinuc_repeat(s, reps=4):
    return any(s[i:i + 2] * reps == s[i:i + 2 * reps] for i in range(len(s) - 2 * reps + 1))


def self_complementary_3p(s, k=4):
    return rc(s[-k:]) in s                      # the 3' end could prime on the primer itself


def hairpin(s, stem=5):
    for i in range(len(s) - stem):
        if rc(s[i:i + stem]) in s[i + stem + 3:]:
            return True
    return False


def pair_3p_dimer(a, b, k=4):
    return rc(a[-k:]) in b or rc(b[-k:]) in a


def good_primer(s, vmask=None, lo=19, hi=25, tm_lo=57.0, tm_hi=63.0):
    if "N" in s or not (lo <= len(s) <= hi):
        return False
    if not (40 <= gc(s) <= 60) or not (tm_lo <= tm_nn(s) <= tm_hi):
        return False
    last5 = sum(c in "GC" for c in s[-5:])
    if last5 < 1 or last5 > 3 or s[-1] not in "GC":
        return False
    if has_run(s) or has_dinuc_repeat(s) or self_complementary_3p(s) or hairpin(s):
        return False
    if vmask is not None and vmask.any():
        return False
    return True


# ------------------------------------------------------------------------------------------------ in-silico PCR
def find_sites(primer, seq, max_mm=3, seed=9, exact_3p=5):
    """Binding sites of `primer` in `seq`: (start0, end0, strand). '+' = primer equals the top strand and extends rightwards."""
    L, sites = len(primer), []
    for strand, p in (("+", primer), ("-", rc(primer))):
        s = p[-seed:] if strand == "+" else p[:seed]            # the primer's 3' end is the right end of p on '+', the left end on '-'
        i = seq.find(s)
        while i != -1:
            start = i - (L - seed) if strand == "+" else i
            if start >= 0 and start + L <= len(seq):
                w = seq[start:start + L]
                mm = [a != b for a, b in zip(w, p)]
                tail = mm[-exact_3p:] if strand == "+" else mm[:exact_3p]
                if sum(mm) <= max_mm and not any(tail):
                    sites.append((start, start + L, strand, sum(mm)))
            i = seq.find(s, i + 1)
    return sites


def amplicons(primers, seq, max_size=5000, **kw):
    """All products between any two primers facing each other. primers: {name: seq}. Returns (name_left, name_right, start0, end0, mm_left, mm_right)."""
    plus, minus = [], []
    for name, p in primers.items():
        for st, en, strand, mm in find_sites(p, seq, **kw):
            (plus if strand == "+" else minus).append((name, st, en, mm))
    out = []
    for n1, s1, e1, m1 in plus:
        for n2, s2, e2, m2 in minus:
            if s2 >= s1 + 10 and e2 - s1 <= max_size:
                out.append((n1, n2, s1, e2, m1, m2))
    return sorted(out, key=lambda x: (x[2], x[3]))


# ------------------------------------------------------------------------------------------------ design
def candidates(seq, vmask, a, b, lo=19, hi=25):
    """Primers (top-strand and bottom-strand) fully inside seq[a:b] that pass the filters. a, b are 0-based offsets."""
    out = []
    for L in range(lo, hi + 1):
        for i in range(a, b - L + 1):
            s = seq[i:i + L]
            if good_primer(s, vmask[i:i + L]):
                out.append(("F", i, i + L, s))
            r = rc(s)
            if good_primer(r, vmask[i:i + L]):
                out.append(("R", i, i + L, r))
    return out


def score(s):
    return abs(tm_nn(s) - 60) + abs(gc(s) - 50) / 10


def best_pair(seq, vmask, f_win, r_win, size_rng, offset, name, single_seq=None, n_keep=3, expect_sites=1, target=None):
    """Enumerate forward x reverse candidates, keep pairs whose Tm differ by <= 2 C, have no 3' dimer, a product in the size range, and the
    expected in-silico PCR behaviour; rank by score. expect_sites: how many binding sites each primer should have in `seq`."""
    fs = [c for c in candidates(seq, vmask, *f_win) if c[0] == "F"]
    rs = [c for c in candidates(seq, vmask, *r_win) if c[0] == "R"]
    fs.sort(key=lambda c: score(c[3]))
    rs.sort(key=lambda c: score(c[3]))
    res = []
    for f in fs[:80]:
        for r in rs[:80]:
            size = r[2] - f[1]
            if not (size_rng[0] <= size <= size_rng[1]) or abs(tm_nn(f[3]) - tm_nn(r[3])) > 2 or pair_3p_dimer(f[3], r[3]):
                continue
            res.append((score(f[3]) + score(r[3]) + abs(tm_nn(f[3]) - tm_nn(r[3])) + (abs(size - target) / 50 if target else 0), f, r, size))
    res.sort(key=lambda x: x[0])
    kept = []
    for sc, f, r, size in res:
        P = {"F": f[3], "R": r[3]}
        a = amplicons(P, seq)
        ok = a and all(x[0] == "F" and x[1] == "R" for x in a)
        row = dict(assay=name, F=f[3], R=r[3], F_start=offset + f[1], F_end=offset + f[2] - 1, R_start=offset + r[1], R_end=offset + r[2] - 1,
                   F_tm=tm_nn(f[3]), R_tm=tm_nn(r[3]), F_gc=gc(f[3]), R_gc=gc(r[3]), size=size, score=sc,
                   products_in_duplicated_allele=len(a),
                   products_in_single_copy_allele=len(amplicons(P, single_seq)) if single_seq is not None else "",
                   sites_F=len([x for x in find_sites(f[3], seq)]), sites_R=len([x for x in find_sites(r[3], seq)]))
        if not ok:
            continue
        if single_seq is not None and name == "J" and (row["products_in_duplicated_allele"] != 1 or row["products_in_single_copy_allele"] != 0
                                                       or row["sites_F"] > 2 or row["sites_R"] > 2):
            continue
        if name != "J" and (row["sites_F"] != expect_sites or row["sites_R"] != expect_sites):
            continue
        if any(abs(row["F_start"] - k["F_start"]) < 15 and abs(row["R_start"] - k["R_start"]) < 15 for k in kept):
            continue
        kept.append(row)
        if len(kept) == n_keep:
            break
    return kept


def design(a):
    if a.bam:
        seq, vmask = consensus_from_bam(a.bam, LO, HI)
        src = f"reads of {os.path.basename(a.bam)} (consensus; N = depth < 3; positions where the reads disagree are avoided)"
    else:
        seq, vmask = read_region_fasta(a.fasta, CHROM, LO, HI)
        src = f"reference sequence {os.path.basename(a.fasta)}"
    J = a.junction
    off = LO - 1                                    # coordinate = off + 0-based index + 0  (1-based genomic coordinate = idx + LO)
    j0 = J - LO                                     # 0-based index of J
    c2_end0 = COPY2_END - LO
    # single-copy allele: delete copy 2 (J .. end of copy 2)
    single = seq[:j0] + seq[c2_end0 + 1:]
    nN = seq[j0 - 800:j0 + 800].count("N")
    print(f"sequence source: {src}\nunit {COPY1_START:,}-{J - 1:,} (copy 1) and {J:,}-{COPY2_END:,} (copy 2), junction J = {J:,}; N bases within 800 bp of J: {nN}")
    jf = (j0 - 800, j0 - 180)                       # forward primer: end of copy 1, at least 180 bp before J
    jr = (j0 + 180, j0 + 800)                       # reverse primer: start of copy 2
    rows = []
    rows += best_pair(seq, vmask, jf, jr, (380, 900), LO, "J", single_seq=single, expect_sites=2, target=550)
    fl = (FL_WIN[0] - LO, FL_WIN[1] - LO)
    rows += best_pair(seq, vmask, fl, fl, (90, 160), LO, "FL", expect_sites=1, target=120)
    u1 = (U_WIN[0] - LO, U_WIN[1] - LO)
    rows += best_pair(seq, vmask, u1, u1, (90, 140), LO, "U", expect_sites=2, target=110)
    if not rows:
        sys.exit("no primer pairs passed the filters; check the sequence input")
    names = {"J": "gsdf_J", "FL": "gsdf_FL", "U": "gsdf_U"}
    out_rows, per = [], {}
    for r in rows:
        per[r["assay"]] = per.get(r["assay"], 0) + 1
        tag = f"{names[r['assay']]}_{per[r['assay']]}"
        r["pair"] = tag
        out_rows.append(r)
    with open(a.out + "_primers.tsv", "w") as fh:
        fh.write("pair\tassay\tprimer\tsequence_5to3\tstart\tend\tstrand\tlength\tTm_C\tGC_pct\tproduct_bp\tproducts_in_duplicated_allele\tproducts_in_single_copy_allele\n")
        for r in out_rows:
            for role, strand in (("F", "+"), ("R", "-")):
                s = r[role]
                fh.write("\t".join(map(str, [r["pair"], r["assay"], r["pair"] + "_" + role, s, r[role + "_start"], r[role + "_end"], strand, len(s),
                                              round(r[role + "_tm"], 1), round(r[role + "_gc"]), r["size"], r["products_in_duplicated_allele"],
                                              r["products_in_single_copy_allele"]])) + "\n")
    with open(a.out + "_region.fa", "w") as fh:
        fh.write(f">{CHROM}:{LO}-{HI} {src}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")
    print(f"\n{'pair':10s} {'assay':5s} {'forward 5->3':28s} {'reverse 5->3':28s} {'product':>7s}  Tm F/R   products dup/single")
    for r in out_rows:
        print(f"{r['pair']:10s} {r['assay']:5s} {r['F']:28s} {r['R']:28s} {r['size']:7d}  {r['F_tm']:.1f}/{r['R_tm']:.1f}   {r['products_in_duplicated_allele']}/{r['products_in_single_copy_allele']}")
    print(f"\nwrote {a.out}_primers.tsv and {a.out}_region.fa")


# ------------------------------------------------------------------------------------------------ check
def read_primers(path):
    prim = {}
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(head, line.rstrip("\n").split("\t")))
            prim[f["primer"]] = f["sequence_5to3"].upper()
    return prim


def iter_fasta(path):
    name, parts = None, []
    with open_text(path) as fh:
        for line in fh:
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts).upper()


def check(a):
    prim = read_primers(a.primers)
    pairs = {}
    for n in prim:
        pairs.setdefault(n.rsplit("_", 1)[0], {})[n.rsplit("_", 1)[1]] = prim[n]
    for path in a.genomes:
        print(f"\n== {os.path.basename(path)}")
        per_pair = {p: [] for p in pairs}
        nsites = {n: 0 for n in prim}
        for contig, seq in iter_fasta(path):
            for pname, pr in pairs.items():
                for n1, n2, s, e, m1, m2 in amplicons({pname + "_F": pr["F"], pname + "_R": pr["R"]}, seq, max_size=a.max_size):
                    per_pair[pname].append((contig, s + 1, e, e - s, m1, m2, n1.split("_")[-1] + "+" + n2.split("_")[-1]))
        for pname, lst in per_pair.items():
            if not lst:
                print(f"  {pname}: no product")
            for contig, s, e, size, m1, m2, kind in lst:
                print(f"  {pname}: {size} bp  {contig}:{s:,}-{e:,}  ({kind}; mismatches {m1}/{m2})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("design")
    g = d.add_mutually_exclusive_group(required=True)
    g.add_argument("--fasta", help="reference genome FASTA (GT3a, chromosome names NC_...)")
    g.add_argument("--bam", help="region BAM of one fish (reads aligned to GT3a over NC_135173.1:50.70-50.83 Mb)")
    d.add_argument("--junction", type=int, default=J_DEFAULT, help="coordinate where copy 2 starts (default 50,773,050)")
    d.add_argument("--out", default="gsdf")
    c = sub.add_parser("check")
    c.add_argument("primers")
    c.add_argument("genomes", nargs="+")
    c.add_argument("--max-size", type=int, default=5000)
    a = ap.parse_args()
    (design if a.cmd == "design" else check)(a)


if __name__ == "__main__":
    main()
