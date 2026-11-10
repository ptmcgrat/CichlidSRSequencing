#!/usr/bin/env python3
"""template_check_with_mates.py - test and polish the element template with the carriers' own reads.

Read pairs with one mate anchored in the flank at the insertion site have the other mate inside the element, within about one fragment length of
a junction. Those mates (and the clipped junction sequences) are reads of THIS insertion, so pooled over all carriers they cover the first and last
~500 bp of the element many times over. This script

  1. fetches, for each carrier, the reads at the two reference copies of the site, the mates of the anomalous ones (by position, from the full BAM),
     and the soft-clipped junction sequences,
  2. aligns every mate to the template (either strand),
  3. reports how well they fit (identity), where on the template they fall for reads left vs right of the site, and how the clipped junction
     sequences line up with the template's two ends,
  4. calls positions where the carriers consistently differ from the template, and writes a polished template with those edits applied.

Only the covered ends (about 500 bp each) can be checked or polished; the interior stays as in the template.

    python3 template_check_with_mates.py            # the 16 carriers, both reference copies, the template from element_from_annotation.py
    python3 template_check_with_mates.py --carriers SAMEA117806259 SAMEA117806260 SAMEA117806270

Run it in the environment where `sexloci` works (it reuses its Dropbox streaming); needs samtools and rclone.
"""
import argparse
import collections
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
DEFAULT_TEMPLATE = "~/Temp/SexDeterminationAnalysis/element_template/element_template.fasta"
DEFAULT_CARRIERS = ["SAMEA117806259", "SAMEA117806260", "SAMEA117806270", "SAMEA117806261", "SAMEA117806262", "SAMEA117806263", "SAMEA117806264",
                    "SAMEA117806265", "SAMEA117806266", "SAMEA117806267", "SAMEA117806271", "SAMEA117806273", "SAMEA117806076", "SAMEA117806077",
                    "SAMEA2661383", "SAMEA2661384"]
DEFAULT_REGIONS = ["NC_135173.1:50764650-50765450", "NC_135173.1:50786400-50787200"]
DEFAULT_SITES = [50765052, 50786786]


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta1(path):
    seq = []
    with open(os.path.expanduser(path)) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


# ------------------------------------------------------------------------------------------------ placing reads on the template
class Template:
    def __init__(self, seq, k=12):
        self.seq, self.k, self.L = seq, k, len(seq)
        self.index = collections.defaultdict(list)
        for i in range(self.L - k + 1):
            self.index[seq[i:i + k]].append(i)

    def place(self, read, min_overlap=80, min_identity=0.9, min_votes=2, band=10):
        """Best gapped local placement of `read` (either strand). The diagonal comes from k-mer votes; a banded Smith-Waterman around it handles
        small indels. Returns dict(strand, start, overlap, identity, bases=[(template pos, base)], indels=[(kind, template pos, length, seq)]) or None."""
        best = None
        for strand, q in (("+", read), ("-", rc(read))):
            votes = collections.Counter()
            for i in range(len(q) - self.k + 1):
                for j in self.index.get(q[i:i + self.k], ()):
                    votes[j - i] += 1
            if not votes:
                continue
            diag, v = max(votes.items(), key=lambda kv: kv[1])
            if v < min_votes:
                continue
            aln = banded_sw(q, self.seq, diag, band)
            if aln is None:
                continue
            cols = aln["cols"]
            aligned_q = sum(1 for c in cols if c[0] in "=X" or c[0] == "I")
            matches = sum(1 for c in cols if c[0] == "=")
            span = len(cols)
            ident = matches / span if span else 0
            bases = [(c[1], c[2]) for c in cols if c[0] in "=X"]
            if len(bases) < min_overlap or ident < min_identity:
                continue
            cand = dict(strand=strand, start=bases[0][0], overlap=len(bases), identity=ident, bases=bases, indels=aln["indels"], matches=matches,
                        score=aln["score"])
            if best is None or cand["score"] > best["score"]:
                best = cand
        return best


def banded_sw(q, ref, diag, band=10, match=2, mism=-3, gap_open=5, gap_ext=2):
    """Local alignment of q against ref restricted to |j - i - diag| <= band (0-based indices), with affine gap costs (a gap of n bases costs
    gap_open + (n-1) * gap_ext, so one contiguous indel is preferred to several scattered ones). Returns dict(score, cols, indels) or None.
    cols: ('=' | 'X', template pos, query base) for aligned columns, ('I', template pos before, query base) for query bases absent from the template,
    ('D', template pos, '-') for template bases absent from the read."""
    n, m, W = len(q), len(ref), 2 * band + 1
    NEG = -10 ** 9
    M = [[NEG] * W for _ in range(n + 1)]
    E = [[NEG] * W for _ in range(n + 1)]          # a template base with no query base (deletion in the read)
    F = [[NEG] * W for _ in range(n + 1)]          # a query base with no template base (insertion in the read)
    PM = [[0] * W for _ in range(n + 1)]
    PE = [[0] * W for _ in range(n + 1)]
    PF = [[0] * W for _ in range(n + 1)]
    best, bi, bk = 0, 0, 0
    for i in range(1, n + 1):
        base_j = i + diag - band
        qi = q[i - 1]
        Mp, Ep, Fp = M[i - 1], E[i - 1], F[i - 1]
        Mi, Ei, Fi = M[i], E[i], F[i]
        for k in range(W):
            j = base_j + k
            if j < 1 or j > m:
                continue
            # match / mismatch: from the diagonal predecessor (same band offset), or a fresh start
            prev, src = 0, 0
            if not (i == 1 or j == 1):
                for val, sc_ in ((Mp[k], 1), (Ep[k], 2), (Fp[k], 3)):
                    if val > prev:
                        prev, src = val, sc_
            Mi[k] = prev + (match if qi == ref[j - 1] else mism)
            PM[i][k] = src
            # insertion in the read: query base i with no template base; comes from (i-1, j) = band offset k+1 in the previous row
            if i > 1 and k + 1 < W:
                a_ = Mp[k + 1] - gap_open if Mp[k + 1] > NEG // 2 else NEG
                b_ = Fp[k + 1] - gap_ext if Fp[k + 1] > NEG // 2 else NEG
                if a_ >= b_:
                    Fi[k], PF[i][k] = a_, 1
                else:
                    Fi[k], PF[i][k] = b_, 2
            # deletion in the read: template base j with no query base; comes from (i, j-1) = band offset k-1 in this row
            if k >= 1:
                a_ = Mi[k - 1] - gap_open if Mi[k - 1] > NEG // 2 else NEG
                b_ = Ei[k - 1] - gap_ext if Ei[k - 1] > NEG // 2 else NEG
                if a_ >= b_:
                    Ei[k], PE[i][k] = a_, 1
                else:
                    Ei[k], PE[i][k] = b_, 2
            h = max(Mi[k], Ei[k], Fi[k], 0)
            if h > best:
                best, bi, bk = h, i, k
    if best <= 0:
        return None
    cols = []
    i, k = bi, bk
    h = max(M[i][k], E[i][k], F[i][k])
    state = "M" if M[i][k] == h else "E" if E[i][k] == h else "F"
    while i > 0:
        j = i + diag - band + k
        if state == "M":
            cols.append(("=" if q[i - 1] == ref[j - 1] else "X", j - 1, q[i - 1]))
            src = PM[i][k]
            i -= 1
            if src == 0:
                break
            state = {1: "M", 2: "E", 3: "F"}[src]
        elif state == "F":
            cols.append(("I", j - 1, q[i - 1]))
            src = PF[i][k]
            i -= 1
            k += 1
            state = "M" if src == 1 else "F"
        else:
            cols.append(("D", j - 1, "-"))
            src = PE[i][k]
            k -= 1
            state = "M" if src == 1 else "E"
    cols.reverse()
    indels, run = [], []
    for c in cols:                                               # collect indel events; same-kind columns within 2 template bases are one event
        if c[0] in "ID":
            if run and run[0][0] == c[0] and 0 <= (c[1] - run[-1][1]) <= 2:
                run.append(c)
                continue
            if run:
                indels.append((run[0][0], run[0][1], len(run), "".join(x[2] for x in run)))
            run = [c]
        elif run:
            indels.append((run[0][0], run[0][1], len(run), "".join(x[2] for x in run)))
            run = []
    if run:
        indels.append((run[0][0], run[0][1], len(run), "".join(x[2] for x in run)))
    return dict(score=best, cols=cols, indels=indels)


def read_fasta_multi(path):
    out, name, seq = [], None, []
    with open(os.path.expanduser(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(seq)))
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip().upper())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def element_ends(tseq, o_s, o_e, max_back=6, max_fwd=4):
    """Where the element really begins and ends. The clipped junction sequences can sit a base or two inside the element (a flank base that matches
    the element's first base is absorbed into the alignment), so the clip offsets (o_s from the start, o_e from the end) are widened until the two ends
    pair up as an inverted repeat. The smallest widening that comes within one base of the longest exact repeat is taken, because a palindromic flank
    such as TA pairs with itself and would otherwise drag the boundary outwards. Returns (start index, end index exclusive, exact run, longest run)."""
    comp = str.maketrans("ACGT", "TGCA")
    L = len(tseq)
    rows = []
    for d in range(0, min(max_back, o_s) + 1):
        for f in range(0, max_fwd + 1):
            st, en = o_s - d, L - o_e + f
            if en > L or en - st < 100:
                continue
            e = tseq[st:en]
            run = 0
            for i in range(min(40, len(e) // 2)):
                if e[i] != e[-1 - i].translate(comp):
                    break
                run += 1
            rows.append((run, d + f, st, en))
    if not rows:
        return o_s, L - o_e, 0, 0
    longest = max(r[0] for r in rows)
    ok = [r for r in rows if r[0] >= longest - 1]
    run, _, st, en = min(ok, key=lambda r: (r[1], -r[0]))
    return st, en, run, longest


def place_end(seq, template, end, max_off=8):
    """Compare a clipped junction sequence with the template's start ('start') or end ('end'), allowing the boundary to sit a few bases inside the
    element (a duplicated base or two can be assigned to the flank). Returns (identity, offset)."""
    best = (0.0, 0)
    for off in range(0, max_off + 1):
        if end == "start":
            ref = template.seq[off:off + len(seq)]
            q = seq[:len(ref)]
        else:
            ref = template.seq[max(template.L - off - len(seq), 0):template.L - off]
            q = seq[len(seq) - len(ref):]
        if len(ref) < 20:
            continue
        ident = sum(a == b for a, b in zip(q, ref)) / len(ref)
        if ident > best[0]:
            best = (ident, off)
    return best


# ------------------------------------------------------------------------------------------------ fetching reads
def fetch_sample(sample, base, regions, sites, run, far=3000, min_clip=20, near=6):
    """Reads of one carrier at the site: [{kind: 'mate'|'unmapped'|'tail'|'head', seq, side, site}]"""
    url = f"{base}/{sample}/{sample}.all.bam"
    reads, anchors, notes = [], [], []
    for region, site in zip(regions, sites):
        rc_, out, err = run(["samtools", "view", url, region])
        if rc_ != 0:
            notes.append(f"{sample}: {region} failed ({err.strip()[:80]})")
            continue
        for line in out.splitlines():
            f = line.split("\t")
            if len(f) < 11:
                continue
            flag, rname, pos, cigar, rnext, pnext, tlen, seq = int(f[1]), f[2], int(f[3]), f[5], f[6], int(f[7]), int(f[8]), f[9]
            if flag & (0x100 | 0x400 | 0x800):
                continue
            ops = [(int(n), o) for n, o in CIG.findall(cigar)] if cigar != "*" else []
            alen = sum(n for n, o in ops if o in "MDN=X") if ops else len(seq)
            side = "left" if pos + alen // 2 < site else "right"      # side of the site the ALIGNED part lies on (not the clipped part)
            if flag & 0x4:
                reads.append(dict(kind="unmapped", seq=seq, side=side, site=site))
                continue
            if flag & 1 and ((rnext not in ("=", "*") and rnext != rname) or (rnext in ("=", rname) and abs(tlen) > far)) and not flag & 0x8:
                anchors.append((f[0], rname if rnext == "=" else rnext, pnext, side, site))
            if ops:
                if ops[-1][1] == "S" and ops[-1][0] >= min_clip and abs(pos + alen - 1 - site) <= near + 4:
                    reads.append(dict(kind="tail", seq=seq[-ops[-1][0]:], side="left", site=site))     # insertion start after the left junction
                if ops[0][1] == "S" and ops[0][0] >= min_clip and abs(pos - site) <= near + 4:
                    reads.append(dict(kind="head", seq=seq[:ops[0][0]], side="right", site=site))      # insertion end before the right junction
    done = set()
    for qn, mc, mp, side, site in anchors:
        if (qn, mc, mp) in done:
            continue
        done.add((qn, mc, mp))
        rc_, out, err = run(["samtools", "view", url, f"{mc}:{mp}-{mp}"])
        if rc_ != 0:
            continue
        for line in out.splitlines():
            f = line.split("\t")
            if f[0] == qn and len(f) >= 11:
                reads.append(dict(kind="mate", seq=f[9], side=side, site=site))
                break
    return reads, notes


def rank_templates(all_reads, cands, args):
    """Score every candidate copy by the summed alignment score of the carriers' mates against it (a read that does not place scores 0)."""
    reads = [r["seq"] for rs in all_reads.values() for r in rs if r["kind"] in ("mate", "unmapped")]
    rows = []
    for name, seq in cands:
        t = Template(seq)
        placed, score, aligned, errs = 0, 0, 0, 0.0
        for q in reads:
            p = t.place(q, args.min_overlap, args.min_identity)
            if p:
                placed += 1
                score += p["score"]
                aligned += p["overlap"]
                errs += p["overlap"] * (1 - p["identity"])
        rows.append(dict(name=name, length=len(seq), placed=placed, score=score, errors_per_kb=1000 * errs / max(aligned, 1)))
    rows.sort(key=lambda r: (-r["score"], r["errors_per_kb"]))
    return rows


# ------------------------------------------------------------------------------------------------ analysis
def analyse(all_reads, template, args, log=print):
    L = template.L
    placed, unplaced = [], []
    for sample, reads in all_reads.items():
        for r in reads:
            if r["kind"] in ("mate", "unmapped"):
                p = template.place(r["seq"], args.min_overlap, args.min_identity)
                if p:
                    placed.append((sample, r["side"], p))
                else:
                    unplaced.append((sample, r))
    n_mates = len(placed) + len(unplaced)
    log(f"template {L:,} bp; {n_mates} mates / unmapped reads from {len(all_reads)} carriers; "
        f"{len(placed)} ({len(placed) / max(n_mates, 1):.0%}) align to the template at >= {args.min_identity:.0%} identity over >= {args.min_overlap} bp "
        f"(gapped; {sum(1 for _, _, p in placed if p['indels'])} of them cross an indel)")
    out = [""]
    if unplaced:
        c = collections.Counter((r["side"], r["kind"]) for _, r in unplaced)
        out.append(f"{len(unplaced)} did not place ({', '.join(f'{k[0]} {k[1]}: {v}' for k, v in sorted(c.items()))}); written to unplaced_reads.fasta")
    for side in ("left", "right"):
        pl = [p for _, s_, p in placed if s_ == side]
        if not pl:
            out.append(f"reads on the {side} of the site: none placed")
            continue
        mid = [p["start"] + p["overlap"] // 2 for _, s_, p in placed if s_ == side]
        first = sum(1 for m in mid if m < L / 2)
        ids = sorted(p["identity"] for p in pl)
        out.append(f"reads on the {side} of the site: {len(pl)} placed; {first} in the first half of the template, {len(pl) - first} in the second; "
                   f"identity median {ids[len(ids) // 2]:.1%} (10th percentile {ids[len(ids) // 10]:.1%}); strand: "
                   f"{sum(1 for p in pl if p['strand'] == '+')} forward, {sum(1 for p in pl if p['strand'] == '-')} reverse")
    out.append("   (with the template oriented A first, reads left of the site belong to its first half and reads right of the site to its second half)")
    # pileup and depth
    cov = collections.Counter()
    pile = collections.defaultdict(collections.Counter)
    who = collections.defaultdict(lambda: collections.defaultdict(set))
    for sample, side, p in placed:
        for j, b in p["bases"]:
            pile[j][b] += 1
            who[j][b].add(sample)
            cov[j // 100] += 1
    out.append("\ndepth along the template (reads pooled over carriers; mean per base in 100 bp bins):")
    nb = (L + 99) // 100
    mx = max(cov.values()) / 100 if cov else 1
    for i in range(nb):
        d = cov.get(i, 0) / 100
        if d > 0 or i < 8 or i > nb - 9:
            out.append(f"   {i * 100:5d}-{min(i * 100 + 99, L - 1):5d} {d:6.1f} {'#' * max(0, round(40 * d / mx)) if d else ''}")
    # clipped junction sequences and the insertion they imply
    offs = {}
    bounds = None
    for kind, end, label in (("tail", "start", "insertion start (clipped after the left junction) vs template start"),
                             ("head", "end", "insertion end (clipped before the right junction) vs template end")):
        rows = [(sample, r) for sample, rs in all_reads.items() for r in rs if r["kind"] == kind]
        if rows:
            res = [place_end(r["seq"], template, end) for _, r in rows]
            c = collections.Counter(o for _, o in res)
            ids = sorted(i for i, _ in res)
            offs[end] = c.most_common(1)[0][0]
            out.append(f"\n{label}: {len(rows)} clipped sequences; identity median {ids[len(ids) // 2]:.1%}; "
                       f"boundary offset inside the template (most common first): {', '.join(f'{o} ({n})' for o, n in c.most_common(3))}")
    out.append(f"\ntemplate begins {template.seq[:14]}... and ends ...{template.seq[-14:]}")
    if "start" in offs and "end" in offs:
        o_s, o_e = offs["start"], offs["end"]
        ins = template.seq[o_s:L - o_e]
        out.append(f"the clipped sequences place the insertion at template positions {o_s + 1}-{L - o_e} ({len(ins):,} bp); it begins {ins[:20]}... and ends ...{ins[-20:]}")
        st, en, run, longest = element_ends(template.seq, o_s, o_e)
        bounds = (st, en) if run >= 6 else (o_s, L - o_e)
        if run >= 6:
            elem = template.seq[st:en]
            out.append(f"   widening the boundaries until the ends pair up gives an exact inverted repeat of {run} bp: the element is template positions {st + 1}-{en} "
                       f"({len(elem):,} bp), beginning {elem[:16]}... and ending ...{elem[-16:]}"
                       + (f"; the pairing can run {longest - run} base(s) further into the flanks, which a palindromic flank such as TA allows" if longest > run else ""))
            out.append(f"   the template's flanking bases there: 5' side '{template.seq[:st]}', 3' side '{template.seq[en:]}'"
                       + (f" (the same two bases on both sides would be a {st}-bp flank duplication)" if st and template.seq[en:en + st] == template.seq[:st] else ""))
        else:
            out.append(f"   no exact inverted repeat of 6 bp or more between the two ends at any of the offsets tried (best {run} bp)")
    # differences from the template: substitutions need support from at least two carriers
    variants = []
    for j in sorted(pile):
        c = pile[j]
        tot = sum(c.values())
        base, n = c.most_common(1)[0]
        if tot >= args.min_depth and base != template.seq[j] and n / tot >= args.min_alt_fraction and len(who[j][base]) >= 2:
            variants.append((j, template.seq[j], base, n, tot))
    depth_ok = sum(1 for j in pile if sum(pile[j].values()) >= args.min_depth)
    out.append(f"\nsubstitutions where the carriers consistently differ from the template (depth >= {args.min_depth}, >= {args.min_alt_fraction:.0%} of reads, "
               f">= 2 carriers): {len(variants)} of {depth_ok} positions checked ({len(variants) / max(depth_ok, 1):.2%})")
    for j, t, b, n, tot in variants[:40]:
        out.append(f"   position {j + 1}: template {t} -> carriers {b}  ({n} of {tot} reads)")
    # indels
    ind = collections.defaultdict(lambda: [0, set()])
    for sample, side, p in placed:
        for kind, pos, ln, sq in p["indels"]:
            ind[(kind, pos, ln, sq)][0] += 1
            ind[(kind, pos, ln, sq)][1].add(sample)
    ind_calls = [(k, v) for k, v in ind.items() if v[0] >= 3 and len(v[1]) >= 2]
    ind_calls.sort(key=lambda kv: kv[0][1])
    out.append(f"\nindels relative to the template seen in >= 3 reads from >= 2 carriers: {len(ind_calls)} (the well-supported ones are applied only to the final estimate)")
    for (kind, pos, ln, sq), (n, smp) in ind_calls[:30]:
        cover = max([sum(pile[q_].values()) for q_ in range(pos - 3, pos + ln + 3) if q_ in pile] or [0])      # depth of the neighbouring bases
        out.append(f"   {'deletion of' if kind == 'D' else 'insertion of'} {ln} bp ({sq}) at template position {pos + 1}: {n} reads from {len(smp)} carriers (depth nearby ~{cover})")
    polished = list(template.seq)
    for j, t, b, n, tot in variants:
        polished[j] = b
    # indels well enough supported to apply: >= min depth reads, >= min alt fraction of the local depth, >= 3 carriers
    applied = []
    for (kind, pos, ln, sq), (n, smp) in ind_calls:
        near = max([sum(pile[q_].values()) for q_ in range(pos - 3, pos + ln + 3) if q_ in pile] or [0])
        if n >= args.min_depth and near and n / near >= args.min_alt_fraction and len(smp) >= 3:
            applied.append((kind, pos, ln, sq, n, near, len(smp)))
    final = list(polished)
    for kind, pos, ln, sq, n, near, ns in sorted(applied, key=lambda x: -x[1]):          # from the right so earlier coordinates stay valid
        if kind == "D":
            del final[pos:pos + ln]
        else:
            final[pos + 1:pos + 1] = list(sq)
    final = "".join(final)
    if bounds:
        final = final[bounds[0]:len(final) - (L - bounds[1])]
    out.append(f"\nindels applied to the final estimate (>= {args.min_depth} reads, >= {args.min_alt_fraction:.0%} of the local depth, >= 3 carriers): {len(applied)}")
    for kind, pos, ln, sq, n, near, ns in applied:
        out.append(f"   {'deletion of' if kind == 'D' else 'insertion of'} {ln} bp ({sq if kind == 'I' else template.seq[pos:pos + ln]}) at template position {pos + 1}: "
                   f"{n} of ~{near} reads, {ns} carriers")
    out.append(f"final estimate: {len(final):,} bp (template trimmed to the element boundaries, {len(variants)} substitutions and {len(applied)} indels applied where the reads reach; "
               f"the interior between the covered ends is the template's)")
    return out, variants, "".join(polished), pile, unplaced, final


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--template", default=DEFAULT_TEMPLATE)
    ap.add_argument("--templates", help="a multi-FASTA of candidate copies (element_modal_copies.fasta from element_from_annotation.py): they are ranked by how "
                    "well the carriers' reads fit each one, and the best is used for the rest of the analysis")
    ap.add_argument("--carriers", nargs="+", default=DEFAULT_CARRIERS)
    ap.add_argument("--regions", nargs="+", default=DEFAULT_REGIONS, help="window around the site in each reference copy")
    ap.add_argument("--sites", nargs="+", type=int, default=DEFAULT_SITES, help="the junction position in each copy (same order as --regions)")
    ap.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--min-overlap", type=int, default=80)
    ap.add_argument("--min-identity", type=float, default=0.9)
    ap.add_argument("--min-depth", type=int, default=8, help="depth needed to call a difference from the template (default 8)")
    ap.add_argument("--min-alt-fraction", type=float, default=0.7)
    ap.add_argument("--outdir", default=os.path.expanduser("~/Temp/SexDeterminationAnalysis/element_template/mate_check"))
    a = ap.parse_args()
    from sexloci.config import Config
    from sexloci.gsdf import HttpServer, default_runner
    from sexloci.labdata import LabData
    lab = LabData.from_config(Config())
    all_reads = {}
    with HttpServer(lab.cloud_master_dir + f"Bamfiles/{a.genome_version}/", None) as base:
        with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as ex:
            futs = {ex.submit(fetch_sample, s, base, a.regions, a.sites, default_runner): s for s in a.carriers}
            for f in as_completed(futs):
                s = futs[f]
                reads, notes = f.result()
                all_reads[s] = reads
                k = collections.Counter(r["kind"] for r in reads)
                print(f"[fetch] {s}: {k['mate']} mates, {k['unmapped']} unmapped, {k['tail']} tails, {k['head']} heads" + (f"  ({'; '.join(notes)})" if notes else ""), flush=True)
    os.makedirs(a.outdir, exist_ok=True)
    label = "template"
    if a.templates:
        cands = read_fasta_multi(a.templates)
        print(f"\nranking {len(cands)} candidate copies by how well the carriers' mates fit them ...", flush=True)
        ranking = rank_templates(all_reads, cands, a)
        best_score = ranking[0]["score"]
        print(f"{'rank':>4s}  {'copy':44s} {'length':>7s} {'placed':>7s} {'differences per kb':>19s} {'score vs best':>14s}")
        for i, r in enumerate(ranking[:12], 1):
            print(f"{i:4d}  {r['name'].split()[0]:44s} {r['length']:7,d} {r['placed']:7d} {r['errors_per_kb']:19.1f} {r['score'] / best_score:14.1%}")
        with open(os.path.join(a.outdir, "template_ranking.tsv"), "w") as fh:
            fh.write("rank\tcopy\tlength\tplaced\tdifferences_per_kb\tscore_vs_best\n")
            for i, r in enumerate(ranking, 1):
                fh.write(f"{i}\t{r['name']}\t{r['length']}\t{r['placed']}\t{r['errors_per_kb']:.2f}\t{r['score'] / best_score:.4f}\n")
        best_name = ranking[0]["name"]
        tpl = Template(dict(cands)[best_name])
        label = best_name.split()[0]
        print(f"\nusing the closest copy, {label}, for the analysis below\n")
    else:
        tpl = Template(read_fasta1(a.template))
    lines, variants, polished, pile, unplaced, final = analyse(all_reads, tpl, a)
    print("\n".join(lines))
    with open(os.path.join(a.outdir, "unplaced_reads.fasta"), "w") as fh:
        for i, (sample, r) in enumerate(unplaced):
            fh.write(f">{sample}_{r['kind']}_{r['side']}_{i}\n{r['seq']}\n")
    with open(os.path.join(a.outdir, "template_differences.tsv"), "w") as fh:
        fh.write("position\ttemplate\tcarriers\treads_supporting\tdepth\n")
        for j, t, b, n, tot in variants:
            fh.write(f"{j + 1}\t{t}\t{b}\t{n}\t{tot}\n")
    with open(os.path.join(a.outdir, "element_template_polished.fasta"), "w") as fh:
        fh.write(f">element_template_polished source={label} length={len(polished)} substitutions={len(variants)} (applied only where carriers' reads cover the ends; indels not applied)\n{polished}\n")
    with open(os.path.join(a.outdir, "element_final_estimate.fasta"), "w") as fh:
        fh.write(f">element_final_estimate source={label} length={len(final)} (the closest GT3a copy, trimmed to the element, with the carriers' substitutions and "
                 f"well-supported indels applied where their reads reach; the interior is that copy's)\n{final}\n")
    print(f"\nwrote template_differences.tsv, element_template_polished.fasta and element_final_estimate.fasta in {a.outdir}")


if __name__ == "__main__":
    main()
