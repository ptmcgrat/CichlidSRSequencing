"""Known sex-determination loci of Lake Malawi cichlids, for marking on the genome plots.

Positions are GT3a (Mzebra_GT3_NCBI) coordinates. Each locus is stored with its linkage group as well as its GT3a accession, so the
marker is also placed when the VCF uses the UMD2a chromosome names for GT3 sequence (the YH-style VCFs): a chromosome matches if it is
named like the accession or if its plot label is the same linkage group. As elsewhere in sexloci, a UMD2a-named VCF is assumed to carry
GT3a coordinates; a VCF in true UMD2a coordinates would put the markers in the wrong place.
"""
from dataclasses import dataclass
from typing import Callable, Iterable, List, Tuple


@dataclass(frozen=True)
class KnownLocus:
    name: str          # label on the plot
    lg: str            # linkage group
    chrom: str         # GT3a accession
    start: int         # GT3a, 1-based inclusive
    end: int
    kind: str          # 'gene' or 'insertion'
    color: str
    note: str = ""

    @property
    def mid(self) -> float:
        return (self.start + self.end) / 2


KNOWN_LOCI: List[KnownLocus] = [
    KnownLocus("gsdf", "LG7", "NC_135173.1", 50_754_902, 50_776_300, "gene", "#7b2cbf",
               "gsdf (gonadal soma derived factor); a 21.7 kb unit containing it is tandemly duplicated on the Y (unit 50,751,322-50,773,049 and copy 2 to 50,794,777)"),
    KnownLocus("TE1", "LG10", "NC_135176.1", 24_921_178, 24_921_178, "insertion", "#00a86b",
               "candidate transposon insertion (+1,139 bp), Deep Benthics"),
    KnownLocus("TE2", "LG10", "NC_135176.1", 24_921_394, 24_921_394, "insertion", "#00a86b",
               "candidate transposon insertion (+9,868 bp), 216 bp from TE1"),
]


def place_known_loci(chroms: Iterable[str], label_of: Callable[[str], str],
                     loci: List[KnownLocus] = None) -> List[Tuple[KnownLocus, str]]:
    """(locus, the plotted chromosome it belongs to). A chromosome named like the locus's accession wins; otherwise the chromosome whose
    plot label is the locus's linkage group. Loci whose chromosome is not in `chroms` are left out."""
    chroms = list(chroms)
    out = []
    for loc in (KNOWN_LOCI if loci is None else loci):
        if loc.chrom in chroms:
            out.append((loc, loc.chrom))
            continue
        same_lg = [c for c in chroms if label_of(c) == loc.lg]
        if len(same_lg) == 1:
            out.append((loc, same_lg[0]))
    return out
