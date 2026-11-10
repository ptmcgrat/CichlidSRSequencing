#!/usr/bin/env python3
"""build_chr7ins_alleles.py - turn the final element estimate into genotyping alleles for the chr7 (gsdf-upstream) insertion.

Why junction alleles and not the whole element: short reads can only decide between "insertion" and "no insertion" at the two junctions, where a read
crosses from flank into element (or back). So the insertion allele only has to carry the first and last few hundred bases of the element, which the
carriers' own reads have verified; the unverified interior would add nothing a read could use. The reference in GT3a holds the gsdf unit twice, and a
carrier has one copy per chromosome, so its reads split between the two reference copies: the same allele is therefore built at both copies and the two
sets of counts are summed afterwards (chr7ins_combine.py).

How the allele is written. The insertion creates a 2 bp duplication of the target site (TG here), so against the reference the insertion is
    REF = the base that ends the first copy of the 2 bp duplication (the G of TG)
    ALT = that base + the element's first --left-bp bases + the element's last --right-bp bases + the duplicated 2 bp (TG)
and the reference then carries on after the anchor. The two element pieces are joined in the middle; no read can cross that seam, only the two
junctions matter. Output is the VCF layout `sexloci sv-genotype --indel-vcf` reads.

    python3 build_chr7ins_alleles.py --final ~/Temp/SexDeterminationAnalysis/element_template/mate_check/element_final_estimate.fasta
    sexloci sv-genotype --indel-vcf chr7ins_alleles.vcf --no-gsdf --outdir ~/Temp/SexDeterminationAnalysis/sv_chr7ins_alleles --jobs 4
    python3 chr7ins_combine.py --genotypes ~/Temp/SexDeterminationAnalysis/sv_chr7ins_alleles/sv_genotypes.tsv
"""
import argparse
import hashlib
import os
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
DEFAULT_FASTA = "~/Temp/CichlidSequencingData/Genomes/Mzebra_GT3_NCBI/GCF_041146795.1_Mzebra_GT3a_genomic.fna"
DEFAULT_FINAL = "~/Temp/SexDeterminationAnalysis/element_template/mate_check/element_final_estimate.fasta"


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta1(path):
    seq = []
    with open(os.path.expanduser(path)) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def read_chrom(path, chrom):
    parts, cur = [], None
    with open(os.path.expanduser(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur == chrom:
                    break
                cur = line[1:].split()[0]
            elif cur == chrom:
                parts.append(line.strip())
    if not parts:
        sys.exit(f"{chrom} not found in {path}")
    return "".join(parts).upper()


def inverted_repeat_run(e):
    n = 0
    for i in range(min(40, len(e) // 2)):
        if e[i] != e[-1 - i].translate(COMP):
            break
        n += 1
    return n


def find_copy2(ref, anchor, offset, search, flank=40, min_identity=0.85):
    """Anchor (1-based) of the same site in the second copy of the unit: best ungapped match of the flanks either side of the anchor."""
    best = None
    L, R = ref[anchor - flank:anchor], ref[anchor:anchor + flank]
    for shift in range(offset - search, offset + search + 1):
        a2 = anchor + shift
        cl, cr = ref[a2 - flank:a2], ref[a2:a2 + flank]
        if len(cl) < flank or len(cr) < flank:
            continue
        ident = (sum(x == y for x, y in zip(L, cl)) + sum(x == y for x, y in zip(R, cr))) / (2 * flank)
        if best is None or ident > best[0]:
            best = (ident, a2)
    return best if best and best[0] >= min_identity else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--final", default=DEFAULT_FINAL, help="element_final_estimate.fasta from template_check_with_mates.py (oriented left end first)")
    ap.add_argument("--fasta", default=DEFAULT_FASTA, help="the GT3a genome FASTA")
    ap.add_argument("--chrom", default="NC_135173.1")
    ap.add_argument("--anchor", type=int, default=50765050, help="1-based position of the base that ends the first copy of the 2 bp duplication, in GT3a copy 1 "
                    "(default 50,765,050: the G of TG; predicted 50,765,051 from Munby's coordinates, observed junction 50,765,049-52)")
    ap.add_argument("--trim-start", type=int, default=0, help="bases to remove from the start of --final (use only if it carries flank bases)")
    ap.add_argument("--trim-end", type=int, default=0, help="bases to remove from the end of --final. element_template_polished.fasta (without the indel edits) "
                    "begins at the element start but carries one flank base at its end: use --trim-end 1 with that file")
    ap.add_argument("--tsd-len", type=int, default=2)
    ap.add_argument("--left-bp", type=int, default=300, help="bases of the element's start kept in the allele (default 300; the reads verify about 300)")
    ap.add_argument("--right-bp", type=int, default=300, help="bases of the element's end kept in the allele (default 300; the reads verify about 400)")
    ap.add_argument("--copy-offset", type=int, default=21728, help="distance to the second reference copy of the gsdf unit (default 21,728)")
    ap.add_argument("--search", type=int, default=400, help="how far either side of the expected offset to look for the second copy")
    ap.add_argument("--no-copy2", action="store_true")
    ap.add_argument("--out", default="chr7ins_alleles.vcf")
    ap.add_argument("--fasta-out", default="chr7ins_alleles.fasta")
    a = ap.parse_args()

    elem = read_fasta1(a.final)
    if a.trim_start or a.trim_end:
        elem = elem[a.trim_start:len(elem) - a.trim_end]
        print(f"trimmed {a.trim_start} base(s) from the start and {a.trim_end} from the end of {os.path.basename(a.final)}")
    print(f"element: {len(elem):,} bp; begins {elem[:16]}... ends ...{elem[-16:]}")
    run = inverted_repeat_run(elem)
    print(f"   exact terminal inverted repeat: {run} bp" + ("" if run >= 8 else "  <-- expected ~13: the element must begin CACTCAG and end CTGAGTG exactly; "
                                                                       "with element_template_polished.fasta add --trim-end 1"))
    if run < 8:
        sys.exit("stopping: the element's ends do not pair up, so the junction alleles would be misplaced (see the note above)")
    if a.left_bp + a.right_bp > len(elem):
        sys.exit("--left-bp + --right-bp exceed the element length")
    ref = read_chrom(a.fasta, a.chrom)
    sites = []
    c1 = a.anchor
    sites.append(("chr7ins_c1", c1, 0.0))
    if not a.no_copy2:
        f2 = find_copy2(ref, c1, a.copy_offset, a.search)
        if f2:
            sites.append(("chr7ins_c2", f2[1], f2[0]))
            print(f"second copy: anchor {f2[1]:,} ({f2[1] - c1:+,} from copy 1; the flanks match copy 1 at {f2[0]:.0%} over 80 bp)")
        else:
            print("second copy: no good match for the flanks near the expected offset; only copy 1 is written (use --search to widen, or --no-copy2)")
    lines, fasta = [], []
    left, right = elem[:a.left_bp], elem[-a.right_bp:]
    for name, anc, ident in sites:
        base = ref[anc - 1]
        tsd = ref[anc - a.tsd_len:anc]
        alt = base + left + right + tsd
        ctx = f"{ref[anc - 12:anc]}|{ref[anc:anc + 12]}"
        note = "" if name == "chr7ins_c1" else "" if tsd == ref[c1 - a.tsd_len:c1] else f"  (this copy's duplicated bases are {tsd}, not {ref[c1 - a.tsd_len:c1]})"
        print(f"{name}: {a.chrom}:{anc:,} anchor {base}; reference context {ctx}; duplicated bases {tsd}; ALT {len(alt)} bp "
              f"(md5 {hashlib.md5(alt.encode()).hexdigest()[:12]}){note}")
        lines.append(f"{a.chrom}\t{anc}\t{name}\t{base}\t{alt}\t100\t*\tTYPE=INSERTION")
        fasta.append(f">{name}_ref {a.chrom}:{anc - 59}-{anc + 60}\n{ref[anc - 60:anc + 60]}\n>{name}_alt (reference to the anchor, the allele, then the reference after it)\n"
                     f"{ref[anc - 60:anc - 1]}{alt}{ref[anc:anc + 60]}\n")
    with open(a.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(a.fasta_out, "w") as fh:
        fh.write("".join(fasta))
    print(f"\nwrote {a.out} and {a.fasta_out}")
    try:                                                              # confirm that sexloci reads it and that REF matches the reference
        from sexloci.svgeno import build_classifiers, load_indel_sites, load_reference_windows
        ss = load_indel_sites(a.out)
        win = load_reference_windows(os.path.expanduser(a.fasta), {a.chrom: (min(s.pos for s in ss) - 700, max(s.pos for s in ss) + 700)})
        build_classifiers(ss, win)
        print("sexloci loads the sites and the REF alleles match the reference.")
    except ImportError:
        print("(sexloci is not importable here, so the sites were not test-loaded)")


if __name__ == "__main__":
    main()
