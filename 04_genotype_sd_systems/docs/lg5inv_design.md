# LG5 (orange blotch / ZW) inversion: how to genotype it

Status: data preparation done for five pools (`python -m sd_genotypes.lg5inv prepare`, 2026-11-07/08); junction counts calibrated in two genera; interior marker discovery written and tested on synthetic data only (`sd_genotypes.lg5inv.markers`), not yet run on the real alignments. Items marked **(unverified)** rest on a table or an
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

## Results so far (pools on D05; junction counts from `junction_count.py`, an alignment-free count of the junction sequences in the region slice and unmapped reads)

The inverted chromosome has two junction points at the right breakpoint: junction a = RC(GT3 19,657,763+) | GT3 26,726,677+ (reads clipped left at 26,726,677) and junction b = GT3 ..26,726,685 | GT3 19,657,763+ forward (reads clipped
right at 26,726,685); both were found in LatrZW (CM105746.1:23,076,449 and 14,558,115) and MezeOBm (CM105983.1:9,701,365 and 1,954,626); the standard junction is in neither assembly. Inverted fraction
p = ((A+B)/2) / ((A+B)/2 + C), C = reads across the standard sequence.

| pool | A | B | C | p |
|---|---|---|---|---|
| SRR18083724 M. zebra Nkhata Bay OB females | 12 | 12 | 9 | 0.57 |
| SRR18083730 M. zebra Nkhata Bay BB females | 0 | 0 | 31 | 0.00 |
| SRR18083717 L. trewavasae Maison OB females | 8 | 11 | 16 | 0.37 |
| SRR18083718 L. trewavasae Maison BB females | 0 | 0 | 29 | 0.00 |
| SRR18083720 M. callainos Luwino white females | 1 | 0 | 29 | 0.017 |

Reading: the junction is present in OB pools of both genera (expected ~0.5 if every OB female is a ZW heterozygote; standard errors ~0.1), absent from BB pools; one read in Luwino cannot separate a rare inverted
chromosome (1 in 40 would give about that) from background. The 12 B+ carrier females of the B chromosome analysis and 6 male pools showed no junction reads at the clipping level. Expected per single fish at ~10x: std/std ~8 C reads
and 0 A/B; std/inv ~4 A, ~4 B, ~4 C; inv/inv ~8 A, ~8 B (extrapolated from pool yields; false-positive rate per fish unmeasured).

**SRR18083730 has no base qualities** (every quality is Q0 in its slice): the SRA copy was rejected because its qualities differed from ENA's in all 1,000 reads checked, and the fallback ENA FASTQ has none. Sequences, mapping and
junction counts are unaffected; nothing may filter on base quality for the pools (MarkDuplicates and its duplicate fraction, 4.75% against 2.3-3.6%, may be affected). The fetch guard does not flag a quality-less fallback yet.

LatrZW aligned to GT3 (`minimap2 -cx asm10 --cs`, inverted segment with 150 kb flank): one chain of minus-strand blocks from GT3 26.73 to 19.66 Mb between plus-strand flanks, with a ~332 kb insertion in the assembly at GT3 ~26.075 Mb
(contiguous in GT3), probably the repeat expansion Behrens et al. describe.

## Marker discovery (`sd_genotypes.lg5inv.markers`, tests/test_lg5inv_markers.py)

`candidates` (variants seen identically in both OB assemblies, left-aligned, in sequence covered exactly once by each, with >= 12 matching bases either side and alignment divergence <= 5%), `count` (SNV allele counts in
the pools with `samtools mpileup`, MAPQ >= 20, no base-quality filter) and `select` (>= 15x in every pool, <= 1 alt read in every negative pool, alt fraction 0.2-0.8 in every positive pool). Negatives: SRR18083730,
SRR18083718, SRR18083720; positives: SRR18083724, SRR18083717. Indels are listed but not counted from reads yet. Not yet verified: the orientation of alt alleles of minus-strand alignments (an inversion allele
showing up in the OB pools would confirm it).

## Held-out and distribution checks in the pools (markers selected on some pools, evaluated on others)

Aggregate alt fraction over the selected sites (selection on pools of one genus, evaluation on the other):

| selected on | markers | trewavasae OB | trewavasae BB | zebra OB | zebra BB | Luwino |
|---|---|---|---|---|---|---|
| zebra pair | 1,741 | 0.526 (held out) | 0.068 (held out) | 0.509 | 0.005 | 0.070 (held out) |
| trewavasae pair | 3,180 | 0.502 | 0.002 | 0.574 (held out) | 0.174 (held out) | 0.196 (held out) |
| both BB + both OB pools | 1,421 | 0.504 | 0.001 | 0.508 | 0.005 | 0.061 (held out) |
| saved set (also clean in Luwino) | 954 | 0.502 | 0.001 | 0.507 | 0.004 | 0.013 (selected) |

The inverted alleles transfer between genera (held-out OB pools 0.53 and 0.57). The held-out BB-type pools are NOT clean (0.06-0.20): many alt alleles also occur on standard chromosomes of other species. Of 10,052 SNV
candidates: 1,597 under 15x in some pool, 2,645 fixed in both BB pools (GT3-private), 4,787 with alt reads in a negative pool, 8 absent from both OB pools, 61 with an OB fraction outside 0.2-0.8, 954 pass.

Per-site distribution of alt reads (observed against a binomial with the same mean): in the OB pool it is a narrow binomial (a real haplotype is uniform). In the L. trewavasae BB pool (zebra-selected sites) 82.5% of sites have 0 alt
reads, ~10% have 6 or more and almost none lie between: standing variation at ~10% of sites and nothing else. Luwino (1,421 sites selected without it) has 40.6% zero, 27.7% one, 14.0% two, 10.3% three to five, 7.5% six or
more: the same ~7.5% standing-variation subset plus a bulk of roughly 1 alt read per site (~3%) that the L. trewavasae BB pool does not have. A uniform ~3% component is what one inverted chromosome in 40 (2.5%) would give,
and it fits the single junction read; cross-contamination is the alternative. Unresolved. Consequence for genotyping a diploid fish: use the MEDIAN per-site alt fraction (robust to the standing-variation sites), not the aggregate.

## Per-sample measurement (`lg5inv.measure`, `lg5inv.call`, config/lg5panel.yaml)

Same shape as the other panel modules: only samtools, a local indexed BAM aligned to Mzebra_GT3_NCBI, a definition hash. Stored per sample: ref and alt read counts at the markers (MAPQ >= 20, no base-quality filter) and the counts of
reads containing the junction sequences (A, B inverted; C standard). Call: n inverted chromosomes (0/1/2) from the median per-site alt fraction (cuts 0.25 and 0.75, PROVISIONAL), with flags for low depth, few markers, near a cut,
aggregate or junction disagreeing with the class, uninformative junction. The definition (`sd_genotypes/lg5inv/data/lg5_markers.tsv`, `lg5_junctions.json`) is built from a `markers select` output and the reference by
`python -m sd_genotypes.lg5inv.definition` and committed. It runs standalone (`python -m sd_genotypes.lg5inv.measure`, which registers itself with the panel at run time); adding it to `panel/measure.py` MODULES is a one-line change
that was deliberately not made. Untested on real BAMs: how `samtools mpileup` without `-f` prints bases (the code handles both letters and `.`/`,`) and `-a` with `-l` and `-r`.
