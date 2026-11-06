# LG5 (orange blotch / ZW) inversion: how to genotype it

Status: design and the data-preparation step (`python -m sd_genotypes.lg5inv prepare`); no marker has been defined or tested on real reads. Items marked **(unverified)** rest on a table or an
alignment coordinate, not on sequence or reads.

## What is known (module 01 run 20261105, `inversion_breakpoints`; Behrens et al. 2025, Mol Ecol)

* Breakpoints in GT3 (NC_135171.1): left 19,657,805 (19,657,762-19,657,848), right 26,726,680 (26,726,676-26,726,685); 7,068,875 bp; both OB assemblies agree (confidence high).
* **GT3 carries the standard (BB) arrangement.** The chromosome-scale contigs of LatrZW (CM105746.1) and MezeOBm (CM105983.1) are inverted between the breakpoints (minus strand); fAulStu2, fDipLim2 and
  fLetChi2 are collinear with GT3 (interior and flanks on the same strand). `reference_is_standard` in `inversion_breakpoints.yaml` is still `null` for LG5.
* **The standard interior is not in the OB assemblies** (only the two main contigs align to it with megabases): GT3 and the three collinear assemblies are the standard-arrangement sequences.
* **Left junction of the inverted haplotype (unverified, from strand runs and PAF coordinates):** GT3 left flank (to 19,657,848) -> ~5.3 kb that does not align -> a plus-strand copy of GT3
  19,660,106-19,670,124 (10,062 / 10,118 bp, both assemblies) -> the start of the inverted segment (GT3 26,726,685). The in-place copy of that 10 kb window is at the other end of the inverted segment.
  The right junction has no gap.
* **Mappability (pilot slice SRR9672204, one B+ M. lombardoi female, ~10x):** the left breakpoint is at the edge of a repeat. From ~19,657,500 to ~19,670,000 the fraction of reads with MAPQ >= 30 is 0.00-0.06
  (a stretch around 19,661,000-19,662,500 reaches 0.34-0.78); outside it 0.94-1.00. The 86 bp uncertainty of the left breakpoint equals the annotated repeat feature 19,657,763-19,657,848. **The right
  breakpoint is clean**: MAPQ >= 30 share 0.96-1.00 over 21 kb and no annotated repeat between 26,722,675 and 26,727,366. Consequence: window-over-flank depth of the 10 kb duplication is NOT a usable marker;
  junction evidence is more promising at the right breakpoint.
* OB colour is not the inversion: M. callainos Luwino and Lupingu carry the OB pigment allele without the inversion (Behrens et al.). Inverted haplotypes of L. trewavasae and M. zebra differ in length
  and TE content; only one repeat expansion is shared between LatrZW and MezeOBm, so most TE insertions in the inverted haplotypes may be lineage-private.

## Approach (in this order)

1. **Right-breakpoint junction evidence.** In an inverted chromosome, reads in unique sequence at >= 26,726,676 are soft-clipped on the left and reads at <= 26,726,685 on the right; a standard chromosome has
   spanning reads and no clips. Count spanning reads and clips (and the junction k-mers) per sample; the pools calibrate bias (BB pools 0, OB female pools ~0.5 inverted).
2. **Interior markers, defined with reads.** Candidate sites: both inverted haplotypes share an allele, GT3 and the three collinear assemblies share the other (needs a `--cs` alignment of LG5; the existing PAFs
   carry `cg` but not `cs`). Keep only unique sequence outside the TE-expansion zones. Each candidate must be confirmed in BB reads (reference allele) and OB reads (inverted allele at the expected frequency).
   Alignment-free alternative: lineage-shared arrangement-specific 31-mers (needs a k-mer counter; not installed).
3. **Insertions.** Insertions shared by both OB assemblies (present/empty site k-mers), if enough are shared.
4. **Cohort checks** (later, from the joint VCF of module 05): diagnostic SNVs, local PCA of the region, and the sex pattern of the calls (female-biased in OB populations, nearly absent in males).

## Data (module `lg5inv`)

Five pools from PRJNA808791 (Behrens et al. 2025, Table 1), all females: M. callainos Luwino white (SRR18083720, hard negative), M. zebra Nkhata Bay OB (SRR18083724) and BB (SRR18083730),
L. trewavasae Maison Reef OB (SRR18083717; BioSample SAMN26105339 = pool of 17 OB females, the ENA title "OB male pool" is a typo) and BB (SRR18083718). Per pool the BAM is cut to NC_135171.1:19.40-26.98 Mb
and to all unmapped reads before it is deleted, and the pool is uploaded as `lg5inv/pool_<run>/` as soon as it is done. Pools give frequencies, not genotypes.

## Open

* The 5.3 kb sequence and the exact junction coordinates (needs the sequences of LatrZW CM105746.1 and MezeOBm CM105983.1, on Dropbox under `genomes/`).
* Whether the junction and marker approaches transfer between genera (only one OB pool per genus here).
* Size of the unmapped-read BAMs (unmeasured).
