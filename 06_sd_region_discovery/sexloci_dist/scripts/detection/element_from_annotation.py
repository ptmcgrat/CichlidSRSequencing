#!/usr/bin/env python3
"""element_from_annotation.py - reconstruct the full-length transposon from intact copies already present in the reference.

A young polymorphic insertion should be nearly identical to intact copies elsewhere in the genome. RepeatMasker marks the two ends of such a copy as
two library families (here rnd-1_family-458 and rnd-1_family-764). Where one family of each sits a few kb from the other, with the pair arranged the
same way in most copies, the genome holds a full-length copy. This script

  1. finds those pairs and reports the arrangements and the length distribution of the spans between them (how long is the element?),
  2. lists the repeat families that lie between the two ends, i.e. the element's body, with their library consensus lengths,
  3. extracts the best candidate copies from the genome FASTA, oriented so the first family comes first, and ranks them by how closely their two
     ends match the junction sequences recovered from the reads (the 69 bp start and 51 bp end by default),
  4. writes the candidates and the best template to FASTA.

    python3 element_from_annotation.py            # fetches the annotation (as mates_vs_repeatmasker.py does) and scans GT3a
    python3 element_from_annotation.py --families rnd-1_family-458 rnd-1_family-764 --max-span 12000

Run it from the folder that holds mates_vs_repeatmasker.py (it reuses that script's annotation download and name translation).
"""
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mates_vs_repeatmasker as mr  # noqa: E402

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
FLIP = {"+": "C", "C": "+"}
DEFAULT_FAMILIES = ["rnd-1_family-458", "rnd-1_family-764"]
DEFAULT_FASTA = "~/Temp/CichlidSequencingData/Genomes/Mzebra_GT3_NCBI/GCF_041146795.1_Mzebra_GT3a_genomic.fna"
# junction sequences recovered from the carriers' reads (candidate #1): the end next to the left flank (family-458 side) and the end next to the right flank
DEFAULT_START = "CTCAGAAAAATAAAGGTGCCAACTGGAACCAAAAGTGGTTCTTTAGACTGATGCCATAGAAGAACCTTT"
DEFAULT_END = "AACCATGAAGACATCTGAAGAACCATTTAAGAACCATAATTTTTCTGAGTG"


def rc(s):
    return s.translate(COMP)[::-1]


# ------------------------------------------------------------------------------------------------ pairs of end hits
def end_hits(ivs, fam_a, fam_b):
    """{chrom: ([A hits], [B hits])}; a hit is (start, end, strand). Strand is the annotation's: '+' or 'C'."""
    out = collections.defaultdict(lambda: ([], []))
    for chrom, lst in ivs.items():
        for s, e, name, _cls, ex in lst:
            strand = ex[0] if ex else "+"
            if name == fam_a:
                out[chrom][0].append((s, e, strand))
            elif name == fam_b:
                out[chrom][1].append((s, e, strand))
    return out


def find_pairs(hits, min_gap, max_gap):
    """Pair every A hit with the nearest B hit on the same chromosome, each B used once. Spans are oriented so family A comes first."""
    pairs = []
    for chrom, (A, B) in hits.items():
        cand = []
        for ai, a in enumerate(A):
            for bi, b in enumerate(B):
                gap = b[0] - a[1] if b[0] > a[1] else a[0] - b[1] if a[0] > b[1] else None
                if gap is not None and min_gap <= gap <= max_gap:
                    cand.append((gap, ai, bi))
        used_a, used_b = set(), set()
        for gap, ai, bi in sorted(cand):
            if ai in used_a or bi in used_b:
                continue
            used_a.add(ai)
            used_b.add(bi)
            a, b = A[ai], B[bi]
            if a[1] < b[0]:                                   # A first as the genome reads: forward
                pairs.append(dict(chrom=chrom, start=a[0], end=b[1], flip=False, a=a, b=b, gap=gap, patt=(a[2], b[2])))
            else:                                             # B first: reverse-complement so that A comes first
                pairs.append(dict(chrom=chrom, start=b[0], end=a[1], flip=True, a=a, b=b, gap=gap, patt=(FLIP[a[2]], FLIP[b[2]])))
    return pairs


def histogram(values, width=500, top=None):
    c = collections.Counter((v // width) * width for v in values)
    keys = sorted(c)
    return [(k, c[k]) for k in keys]


# ------------------------------------------------------------------------------------------------ genome
def extract(fasta, wanted):
    """wanted: {chrom: [(key, start, end)]} (1-based inclusive). One streaming pass; returns {key: sequence}."""
    out = {}
    cur, buf = None, []
    def flush():
        if cur in wanted:
            seq = "".join(buf).upper()
            for key, s, e in wanted[cur]:
                out[key] = seq[s - 1:e]
    with open(os.path.expanduser(fasta)) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                cur, buf = line[1:].split()[0], []
            elif cur in wanted:
                buf.append(line.strip())
    flush()
    return out


def best_ungapped(seed, window):
    """Best ungapped identity of `seed` against `window` (either strand of the seed); returns (identity, offset, strand)."""
    best = (0.0, 0, "+")
    n = len(seed)
    for strand, q in (("+", seed), ("-", rc(seed))):
        for i in range(0, len(window) - n + 1):
            m = sum(1 for a, b in zip(q, window[i:i + n]) if a == b)
            if m / n > best[0]:
                best = (m / n, i, strand)
    return best


# ------------------------------------------------------------------------------------------------ report
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", nargs=2, default=DEFAULT_FAMILIES, metavar=("FAMILY_A", "FAMILY_B"),
                    help="the two end families; A is the one on the left-of-site side (family-458 for candidate #1)")
    ap.add_argument("--repeatmasker", help="a local RepeatMasker .out instead of fetching the lab's annotation")
    ap.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    ap.add_argument("--annotation-remote")
    ap.add_argument("--annotation-dir")
    ap.add_argument("--fasta", default=DEFAULT_FASTA, help="the genome FASTA (names must match the annotation after the name map)")
    ap.add_argument("--min-gap", type=int, default=300, help="smallest distance between the two end hits (default 300)")
    ap.add_argument("--max-span", type=int, default=12000, help="largest distance between the two end hits (default 12000)")
    ap.add_argument("--start-seq", default=DEFAULT_START, help="junction sequence at the A end (default: the 69 bp from the carriers)")
    ap.add_argument("--end-seq", default=DEFAULT_END, help="junction sequence at the B end (default: the 51 bp from the carriers)")
    ap.add_argument("--per-class", type=int, default=60, help="copies extracted and scored from each length class (default 60)")
    ap.add_argument("--keep", type=int, default=12, help="copies written to the candidates FASTA (default 12)")
    ap.add_argument("--outdir", default=os.path.expanduser("~/Temp/SexDeterminationAnalysis/element_template"))
    a = ap.parse_args()
    fam_a, fam_b = a.families
    os.makedirs(a.outdir, exist_ok=True)

    # annotation (same download / name translation as mates_vs_repeatmasker.py)
    if a.repeatmasker:
        annot, nm_path = a.repeatmasker, None
    else:
        files = mr.fetch_annotation(a.annotation_remote or mr.REMOTE_TEMPLATE.format(genome=a.genome_version),
                                    a.annotation_dir or mr.LOCAL_DIR_TEMPLATE.format(genome=a.genome_version))
        annot, nm_path = files.get("Mzebra_GT3a.fasta.out"), files.get("name_map.tsv")
        if not annot:
            sys.exit("the .out file is needed (the .gff has no consensus coordinates or classes)")
    ivs = mr.read_annotation(annot)
    if nm_path:
        mapping, msg = mr.load_name_map(nm_path, ivs.keys())
        print("[annotation]", msg)
        if mapping:
            moved = collections.defaultdict(list)
            for c, lst in ivs.items():
                moved[mapping.get(c, c)] += lst
            ivs = moved
    idx = mr.Index(ivs)
    cons = mr.consensus_lengths(ivs)
    hits = end_hits(ivs, fam_a, fam_b)
    nA, nB = sum(len(v[0]) for v in hits.values()), sum(len(v[1]) for v in hits.values())
    print(f"{fam_a}: {nA:,} hits (consensus {cons.get(fam_a, 0):,} bp); {fam_b}: {nB:,} hits (consensus {cons.get(fam_b, 0):,} bp)\n")
    if not nA or not nB:
        sys.exit("one of the families has no hits; check the names")

    pairs = find_pairs(hits, a.min_gap, a.max_span)
    print(f"{len(pairs):,} {fam_a} / {fam_b} pairs within {a.min_gap:,}-{a.max_span:,} bp of each other")
    patt = collections.Counter(p["patt"] for p in pairs)
    print("arrangement of the pair, oriented so that family A is first (strand of A, strand of B; '+' = same strand as the library consensus):")
    for pt, n in patt.most_common():
        sp = [p["end"] - p["start"] + 1 for p in pairs if p["patt"] == pt]
        sp.sort()
        print(f"   A {pt[0]:>1s}  B {pt[1]:>1s}: {n:5d} pairs; span median {sp[len(sp) // 2]:,} bp (quartiles {sp[len(sp) // 4]:,} - {sp[3 * len(sp) // 4]:,})")
    if not pairs:
        sys.exit("no pairs: try a larger --max-span")
    dom = patt.most_common(1)[0][0]
    sel = [p for p in pairs if p["patt"] == dom]
    spans = [p["end"] - p["start"] + 1 for p in sel]
    print(f"\nlength of the span between the outer ends of the pair, arrangement A {dom[0]} B {dom[1]} ({len(sel)} copies; 500 bp bins):")
    h = histogram(spans)
    mx = max(c for _, c in h)
    for k, c in h:
        print(f"   {k:6,d}-{k + 499:6,d} {c:4d} {'#' * max(1, round(40 * c / mx))}")
    mode = max(h, key=lambda kc: kc[1])[0] + 250
    print(f"   most common length about {mode:,} bp")

    # the body: families lying between the two ends in the dominant arrangement
    body = collections.Counter()
    for p in sel:
        lo, hi = (p["a"][1], p["b"][0]) if not p["flip"] else (p["b"][1], p["a"][0])
        seen = set()
        for h_ in idx.overlapping(p["chrom"], lo + 1, hi - 1):
            key = (h_[2], h_[3])
            if key not in seen and h_[2] not in (fam_a, fam_b):
                seen.add(key)
                body[key] += 1
    print(f"\nrepeat families found between the two ends (copies of {len(sel)} that contain them):")
    for (name, cls), n in body.most_common(12):
        print(f"   {name:24s} {cls:22s} {n:4d} copies   consensus {cons.get(name, 0):,} bp")

    # length classes relative to the modal length: shorter / the mode / longer (e.g. a nested insertion) / much longer
    def klass(span):
        d = span - mode
        return "shorter" if d < -600 else "modal" if d <= 600 else "longer (+0.6 to +3 kb)" if d <= 3000 else "much longer (> +3 kb)"
    classes = collections.defaultdict(list)
    for p in sel:
        classes[klass(p["end"] - p["start"] + 1)].append(p)
    order = ["shorter", "modal", "longer (+0.6 to +3 kb)", "much longer (> +3 kb)"]
    cand = []
    for c in order:
        lst = sorted(classes.get(c, []), key=lambda p: abs((p["end"] - p["start"] + 1) - mode))[:a.per_class]
        for p in lst:
            p["klass"] = c
        cand += lst
    wanted = collections.defaultdict(list)
    for i, p in enumerate(cand):
        wanted[p["chrom"]].append((i, p["start"], p["end"]))
    fasta = os.path.expanduser(a.fasta)
    if not os.path.exists(fasta):
        sys.exit(f"genome FASTA {fasta} not found (use --fasta)")
    print(f"\nreading {len(cand)} copies from {os.path.basename(fasta)} and scoring both ends against the junction sequences ...", flush=True)
    seqs = extract(fasta, wanted)
    scored = []
    zone = 700
    for i, p in enumerate(cand):
        s_ = seqs.get(i)
        if not s_:
            continue
        if p["flip"]:
            s_ = rc(s_)
        sa = best_ungapped(a.start_seq, s_[:zone])
        sb = best_ungapped(a.end_seq, s_[-zone:])
        scored.append((min(sa[0], sb[0]), sa, sb, p, s_))
    scored.sort(key=lambda x: -x[0])

    print("\nby length class (copies with both junction sequences at >= 95% identity, and the families found between the ends):")
    for c in order:
        rows = [x for x in scored if x[3]["klass"] == c]
        if not rows:
            continue
        lens = sorted(len(x[4]) for x in rows)
        good = sum(1 for x in rows if x[0] >= 0.95)
        fams = collections.Counter()
        for x in rows:
            p = x[3]
            lo, hi = (p["a"][1], p["b"][0]) if not p["flip"] else (p["b"][1], p["a"][0])
            seen = set()
            for h_ in idx.overlapping(p["chrom"], lo + 1, hi - 1):
                if h_[2] not in (fam_a, fam_b) and h_[2] not in seen:
                    seen.add(h_[2])
                    fams[(h_[2], h_[3])] += 1
        print(f"   {c:26s} {len(classes[c]):4d} copies (scored {len(rows):3d}); lengths {lens[0]:,}-{lens[-1]:,} (median {lens[len(lens) // 2]:,}); "
              f"both ends >= 95%: {good} of {len(rows)}")
        for (name, cls), n in fams.most_common(4):
            print(f"        {name:22s} {cls:20s} in {n:3d} of {len(rows)} scored copies  (consensus {cons.get(name, 0):,} bp)")

    print(f"\n{'rank':>4s}  {'copy':36s} {'length':>7s}  {'start 69 bp':>11s} {'end 51 bp':>10s}  class")
    for r, (score, sa, sb, p, s_) in enumerate(scored[:a.keep], 1):
        print(f"{r:4d}  {p['chrom']}:{p['start']}-{p['end']:<13d} {len(s_):7,d}  {sa[0]:10.0%}{'' if sa[2] == '+' else '*'} {sb[0]:9.0%}{'' if sb[2] == '+' else '*'}  {p['klass']}")
    print("   (* = the junction sequence matches the reverse complement; both seeds are expected on the + strand when the copy is oriented A first)")
    fa = os.path.join(a.outdir, "element_candidates.fasta")
    with open(fa, "w") as fh:
        for r, (score, sa, sb, p, s_) in enumerate(scored[:a.keep], 1):
            fh.write(f">candidate{r} {p['chrom']}:{p['start']}-{p['end']} orientation={'rc' if p['flip'] else 'fwd'} length={len(s_)} class={p['klass'].split(' ')[0]} "
                     f"start_match={sa[0]:.2f} end_match={sb[0]:.2f}\n{s_}\n")
    # one best copy per class, so the length classes can be compared
    with open(os.path.join(a.outdir, "element_best_per_class.fasta"), "w") as fh:
        for c in order:
            rows = [x for x in scored if x[3]["klass"] == c]
            if rows:
                x = rows[0]
                fh.write(f">best_{c.split(' ')[0]} {x[3]['chrom']}:{x[3]['start']}-{x[3]['end']} length={len(x[4])} start_match={x[1][0]:.2f} end_match={x[2][0]:.2f}\n{x[4]}\n")
    if scored:
        # template = the most typical intact copy: among modal-length copies whose ends match the junction sequences, the one with the highest mean
        # 16-mer Jaccard similarity to the others (ties at 100% on the junction sequences say nothing about which copy is representative)
        pool = [x for x in scored if x[3]["klass"] == "modal" and x[0] >= 0.95]
        basis = "modal-length copies with both ends >= 95%"
        if len(pool) < 3:
            pool = [x for x in scored if x[3]["klass"] == "modal"] or scored
            basis = "modal-length copies" if pool is not scored else "all scored copies"
        pool = pool[:60]
        kset = [{x[4][i:i + 16] for i in range(len(x[4]) - 15)} for x in pool]
        def mean_j(i):
            return sum(len(kset[i] & kset[j]) / max(len(kset[i] | kset[j]), 1) for j in range(len(pool)) if j != i) / max(len(pool) - 1, 1)
        k_best = max(range(len(pool)), key=mean_j)
        best = pool[k_best]
        with open(os.path.join(a.outdir, "element_template.fasta"), "w") as fh:
            fh.write(f">element_template {best[3]['chrom']}:{best[3]['start']}-{best[3]['end']} length={len(best[4])} start_match={best[1][0]:.2f} "
                     f"end_match={best[2][0]:.2f} medoid_of={len(pool)}\n{best[4]}\n")
        with open(os.path.join(a.outdir, "element_modal_copies.fasta"), "w") as fh:        # every candidate for template_check_with_mates.py --templates
            for x in (pool if len(pool) >= 3 else scored[:a.keep]):
                fh.write(f">{x[3]['chrom']}:{x[3]['start']}-{x[3]['end']} orientation={'rc' if x[3]['flip'] else 'fwd'} length={len(x[4])} "
                         f"start_match={x[1][0]:.2f} end_match={x[2][0]:.2f}\n{x[4]}\n")
        print(f"\ntemplate: {best[3]['chrom']}:{best[3]['start']}-{best[3]['end']} ({len(best[4]):,} bp), the most typical of {len(pool)} {basis} "
              f"(mean 16-mer Jaccard to the others {mean_j(k_best):.2f}; junction sequences matched at {best[1][0]:.0%} and {best[2][0]:.0%})")
        print(f"wrote element_candidates.fasta, element_best_per_class.fasta, element_modal_copies.fasta ({len(pool)} copies) and element_template.fasta in {a.outdir}")
        if best[0] < 0.85:
            print("note: even the best copy matches the junction sequences below 85%: the carriers' element may belong to a subfamily that GT3a "
                  "lacks, or the junction sequences carry errors (the 69 bp start came from only 5 reads).")


if __name__ == "__main__":
    main()
