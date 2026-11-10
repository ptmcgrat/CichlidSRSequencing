#!/usr/bin/env python3
"""gsdf_dup_check.py - is the tandem duplication around gsdf present in a genome assembly?

Needs only numpy. Two steps.

1) Build the probe once from the M. zebra GT3a reference (left flank, the 22 kb duplicated unit, right flank):

   python3 gsdf_dup_check.py probe GT3a_genome.fna --out gsdf_probe.fa

   The defaults are the coordinates found on NC_135173.1 (unit = copy 1, 50,751,000-50,772,999; copy 2 follows 21,728 bp later).

2) Check any number of assemblies (plain or .gz FASTA):

   python3 gsdf_dup_check.py check gsdf_probe.fa genomeA.fa genomeB.fa.gz --out results

   Writes results.tsv (one row per genome) and results_copies.tsv (every copy found, with its contig and coordinates).
   Run it on the reference itself as a positive control: it must report DUPLICATED with 2 copies.

How it works: every 25-mer of the probe is looked up in the assembly. Matches that fall on one diagonal (same offset between probe and
contig) form one copy of the unit; flank matches show where the locus is. The number of copies of the unit between the two flanks is
the answer. The distance between the flanks gives an independent estimate (each extra copy adds one repeat period).

Calls: DUPLICATED (2 copies), MULTI (3+), SINGLE (1 copy), each with a '?' when a flank is missing on that contig (the copy number may be
underestimated if the contig ends near the locus); MIXED (haplotypes differ, e.g. a phased assembly with one duplicated and one single
haplotype); UNRESOLVED (a contig breaks inside the locus, so the copy number cannot be read); ALTERED (a copy is partly deleted between
intact flanks); FRAGMENTED (unit only partly present); NOT FOUND.
IMPORTANT: a nearly identical 22 kb tandem duplicate is exactly what assemblers collapse, so SINGLE in an assembly is not proof of absence;
confirm with read depth over the unit (males carry about 1.5x the female depth if the Y carries the duplication).
"""
import argparse
import gzip
import os
import sys

import numpy as np

K = 25
LUT = np.full(256, 4, np.uint8)
for _i, _c in enumerate("ACGT"):
    LUT[ord(_c)] = _i
    LUT[ord(_c.lower())] = _i


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def iter_fasta(path):
    """yield (full header line without '>', sequence)"""
    head, parts = None, []
    with open_text(path) as fh:
        for line in fh:
            if line[0] == ">":
                if head is not None:
                    yield head, "".join(parts)
                head, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
    if head is not None:
        yield head, "".join(parts)


def codes(digits, k):
    """canonical k-mer code, strand flag (+1: forward is canonical) and validity for every window of a digit array"""
    a = digits.astype(np.int64)
    bad = a > 3
    a[bad] = 0
    fw = np.correlate(a, 4 ** np.arange(k - 1, -1, -1, dtype=np.int64), "valid")
    rc = np.correlate(3 - a, 4 ** np.arange(k, dtype=np.int64), "valid")
    cs = np.concatenate(([0], np.cumsum(bad)))
    valid = ((cs[k:] - cs[:-k]) == 0) & (fw != rc)
    return np.minimum(fw, rc), np.where(fw <= rc, 1, -1).astype(np.int8), valid


# ------------------------------------------------------------------------------------------------ probe
def read_region(path, chrom, start, end):
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
        sys.exit(f"{chrom}:{start}-{end} not found in {path} (chromosome names: grep '>' {path} | head)")
    return "".join(parts).upper()


def rng(text):
    a, b = text.replace(",", "").split("-")
    return int(a), int(b)


def cmd_probe(a):
    (l0, l1), (u0, u1), (r0, r1) = rng(a.left), rng(a.unit), rng(a.right)
    with open(a.out, "w") as fh:
        for name, (s, e) in (("L", (l0, l1)), ("U", (u0, u1)), ("R", (r0, r1))):
            seq = read_region(a.reference, a.chrom, s, e)
            extra = f" period={a.period} ref_span={r0 - l1}" if name == "U" else ""
            fh.write(f">{name} chrom={a.chrom} start={s} end={e}{extra}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
            print(f"{name}: {a.chrom}:{s:,}-{e:,} ({len(seq):,} bp, {seq.count('N')} N)")
    print(f"wrote {a.out}; reference span between the flanks (2 copies): {r0 - l1:,} bp, repeat period {a.period:,} bp")


class Probe:
    pass


def load_probe(path, k):
    P = Probe()
    P.names, P.length, P.meta = [], [], {}
    cs, ids, qs, ss = [], [], [], []
    for head, seq in iter_fasta(path):
        name = head.split()[0]
        kv = dict(x.split("=", 1) for x in head.split()[1:] if "=" in x)
        digits = LUT[np.frombuffer(seq.upper().encode(), np.uint8)]
        can, strand, valid = codes(digits, k)
        j = np.flatnonzero(valid)
        cs.append(can[j]); ids.append(np.full(len(j), len(P.names))); qs.append(j); ss.append(strand[j])
        P.names.append(name); P.length.append(len(seq)); P.meta[name] = kv
    code, pid, q, strand = map(np.concatenate, (cs, ids, qs, ss))
    u, inv, cnt = np.unique(code, return_inverse=True, return_counts=True)
    keep = cnt[inv] == 1                                    # k-mers that occur once within the probe
    code, pid, q, strand = code[keep], pid[keep], q[keep], strand[keep]
    o = np.argsort(code)
    P.code, P.pid, P.q, P.strand = code[o], pid[o], q[o], strand[o]
    P.nk = np.bincount(P.pid, minlength=len(P.names))
    h = (P.code.astype(np.uint64) * np.uint64(11400714819323198485)) >> np.uint64(40)
    P.bloom = np.zeros(1 << 24, bool)
    P.bloom[h.astype(np.int64)] = True
    U = P.meta.get("U", {})
    P.period = int(U.get("period", 21728))
    P.ref_span = int(U.get("ref_span", 0))
    return P


# ------------------------------------------------------------------------------------------------ scanning
def scan_contig(digits, P, k, chunk=8_000_000):
    out = []
    n = len(digits)
    for s in range(0, max(n - k + 1, 0), chunk):
        can, strand, valid = codes(digits[s:min(s + chunk + k - 1, n)], k)
        h = ((can.astype(np.uint64) * np.uint64(11400714819323198485)) >> np.uint64(40)).astype(np.int64)
        j = np.flatnonzero(valid & P.bloom[h])
        if not len(j):
            continue
        ii = np.searchsorted(P.code, can[j])
        ii[ii >= len(P.code)] = len(P.code) - 1
        ok = P.code[ii] == can[j]
        j, ii = j[ok], ii[ok]
        out.append((P.pid[ii], strand[j] == P.strand[ii], P.q[ii], s + j))
    if not out:
        return None
    return tuple(np.concatenate(x) for x in zip(*out))


def copies_from_matches(contig, m, P, k, min_frac, gap=3000):
    pid, same, q, t = m
    items = []
    for p in np.unique(pid):
        for sm in (True, False):
            sel = (pid == p) & (same == sm)
            if not sel.any():
                continue
            qq, tt = q[sel], t[sel]
            d = tt - qq if sm else tt + qq
            o = np.argsort(d)
            d, qq, tt = d[o], qq[o], tt[o]
            cuts = np.flatnonzero(np.diff(d) > gap) + 1
            for a, b in zip(np.r_[0, cuts], np.r_[cuts, len(d)]):
                nq = len(np.unique(qq[a:b]))
                frac = nq / P.nk[p]
                if frac >= min_frac[P.names[p]]:
                    items.append(dict(contig=contig, probe=P.names[p], strand="+" if sm else "-", start=int(tt[a:b].min()) + 1,
                                      end=int(tt[a:b].max()) + k, frac=frac, cover=(qq[a:b].max() - qq[a:b].min() + k) / P.length[p]))
    return items


def scan_genome(path, P, k, min_frac):
    items = []
    for head, seq in iter_fasta(path):
        if len(seq) < k:
            continue
        m = scan_contig(LUT[np.frombuffer(seq.encode(), np.uint8)], P, k)
        if m is not None:
            items += copies_from_matches(head.split()[0], m, P, k, min_frac)
    return items


# ------------------------------------------------------------------------------------------------ interpretation
def loci(items, join=80_000):
    by = {}
    for it in items:
        by.setdefault(it["contig"], []).append(it)
    groups = []
    for c, lst in by.items():
        lst.sort(key=lambda x: x["start"])
        cur = [lst[0]]
        for it in lst[1:]:
            if it["start"] - max(x["end"] for x in cur) <= join:
                cur.append(it)
            else:
                groups.append(cur)
                cur = [it]
        groups.append(cur)
    return [g for g in groups if any(x["probe"] == "U" for x in g)]


def evaluate(g, P):
    Us = sorted([x for x in g if x["probe"] == "U"], key=lambda x: x["start"])
    # the same stretch matched on both strands (inverted repeat) is not a second copy
    uniq = []
    for u in Us:
        if all(min(u["end"], v["end"]) - max(u["start"], v["start"]) < 0.5 * (u["end"] - u["start"]) for v in uniq):
            uniq.append(u)
    Us = uniq
    Ls = [x for x in g if x["probe"] == "L"]
    Rs = [x for x in g if x["probe"] == "R"]
    L = max(Ls, key=lambda x: x["frac"]) if Ls else None
    R = max(Rs, key=lambda x: x["frac"]) if Rs else None
    notes, n = [], len(Us)
    flanks_ok, span_n = False, None
    if L and R:
        fwd = L["start"] < R["start"]
        first, last = (L, R) if fwd else (R, L)
        inside = first["end"] <= Us[0]["start"] + 3000 and last["start"] >= Us[-1]["end"] - 3000
        flanks_ok = inside
        if inside and P.ref_span:
            span = last["start"] - first["end"]
            span_n = 2 + (span - P.ref_span) / P.period
            if abs(span_n - n) > 0.4:
                notes.append(f"span between flanks suggests {span_n:.1f} copies")
        if not inside:
            notes.append("flanks do not surround the unit as expected")
    else:
        notes.append("left flank missing" if not L else "right flank missing")
        if not L and not R:
            notes[-1] = "both flanks missing"
    offs = [Us[i + 1]["start"] - Us[i]["start"] for i in range(n - 1)]
    if offs and any(not (0.6 * P.period < o < 1.5 * P.period) for o in offs):
        notes.append("copies are not one repeat period apart: " + ",".join(f"{o:,}" for o in offs))
    n_part = sum(u["cover"] < 0.75 for u in Us)
    partial = n_part == n                                  # no complete copy at all
    if n_part and not partial:
        notes.append(f"{n_part} of {n} copies only partly present")
    if partial:
        notes.append("unit only partly present")
    arrangement = " ".join(f'{x["probe"]}{x["strand"]}' for x in sorted(g, key=lambda x: x["start"]))
    ident = [u["frac"] ** (1 / K) for u in Us]
    return dict(n=n, n_part=n_part, flanks_ok=flanks_ok, span_n=span_n, notes=notes, arrangement=arrangement, partial=partial,
                contig=g[0]["contig"], start=min(x["start"] for x in g), end=max(x["end"] for x in g), ident=ident)


def call_from(evs):
    if not evs:
        return "NOT FOUND"
    if all(e["partial"] for e in evs):
        return "FRAGMENTED (unit only partly present)"
    altered = [e for e in evs if e["n_part"] and e["flanks_ok"]]
    if altered:
        return "ALTERED (a copy of the unit is partly missing between intact flanks)"
    broken = [e for e in evs if e["n_part"] and not e["flanks_ok"]]
    if broken:
        return ("UNRESOLVED (the assembly breaks inside the locus: "
                + " | ".join(f'{e["contig"]}: {e["n"] - e["n_part"]} full + {e["n_part"]} partial' for e in evs) + ")")

    def one(e):
        q = "" if e["flanks_ok"] else "?"
        return ("SINGLE" if e["n"] == 1 else "DUPLICATED" if e["n"] == 2 else "MULTI") + q
    sizes = {e["n"] for e in evs}
    if len(sizes) > 1:
        return "MIXED (" + " | ".join(f'{e["n"]} cop{"y" if e["n"] == 1 else "ies"}' for e in evs) + ")"
    return one(evs[0])


def cmd_check(a):
    P = load_probe(a.probe, a.k)
    min_frac = {"L": a.min_flank_frac, "U": a.min_unit_frac, "R": a.min_flank_frac}
    rows, copy_rows = [], []
    for path in a.genomes:
        name = os.path.basename(path)
        print(f"scanning {name} ...", flush=True)
        items = scan_genome(path, P, a.k, min_frac)
        evs = [evaluate(g, P) for g in loci(items)]
        call = call_from(evs)
        row = dict(genome=name, call=call, copies="|".join(str(e["n"]) for e in evs) or "0",
                   copies_by_flank_span="|".join("-" if e["span_n"] is None else f'{e["span_n"]:.1f}' for e in evs) or "-",
                   locus="|".join(f'{e["contig"]}:{e["start"]:,}-{e["end"]:,}' for e in evs) or "-",
                   arrangement="|".join(e["arrangement"] for e in evs) or "-",
                   unit_identity="|".join("/".join(f"{x:.3f}" for x in e["ident"]) for e in evs) or "-",
                   notes="; ".join(n for e in evs for n in e["notes"]) or "-")
        rows.append(row)
        for it in items:
            copy_rows.append(dict(genome=name, contig=it["contig"], probe=it["probe"], strand=it["strand"], start=it["start"], end=it["end"],
                                  fraction_of_probe_matched=round(it["frac"], 3), est_identity=round(it["frac"] ** (1 / a.k), 3)))
        print(f"  {call}   copies={row['copies']}  {row['locus']}  {row['arrangement']}  {('(' + row['notes'] + ')') if row['notes'] != '-' else ''}", flush=True)
    import csv
    with open(a.out + ".tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(a.out + "_copies.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome", "contig", "probe", "strand", "start", "end", "fraction_of_probe_matched", "est_identity"], delimiter="\t")
        w.writeheader()
        w.writerows(copy_rows)
    print(f"\nwrote {a.out}.tsv and {a.out}_copies.tsv")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("probe", help="build the probe from the reference genome")
    p.add_argument("reference")
    p.add_argument("--chrom", default="NC_135173.1")
    p.add_argument("--left", default="50735000-50749000", help="left flank, ends before the unit (default 14 kb)")
    p.add_argument("--unit", default="50752500-50771500", help="the duplicated unit: copy 1, trimmed by 1.5 kb at each end")
    p.add_argument("--right", default="50797000-50811000", help="right flank, starts after copy 2")
    p.add_argument("--period", type=int, default=21728, help="distance between the two copies in the reference")
    p.add_argument("--out", default="gsdf_probe.fa")
    c = sub.add_parser("check", help="check assemblies")
    c.add_argument("probe")
    c.add_argument("genomes", nargs="+")
    c.add_argument("--out", default="gsdf_dup_results")
    c.add_argument("--k", type=int, default=K)
    c.add_argument("--min-unit-frac", type=float, default=0.2, help="fraction of the unit's k-mers needed to call a copy")
    c.add_argument("--min-flank-frac", type=float, default=0.1)
    a = ap.parse_args()
    (cmd_probe if a.cmd == "probe" else cmd_check)(a)


if __name__ == "__main__":
    main()
