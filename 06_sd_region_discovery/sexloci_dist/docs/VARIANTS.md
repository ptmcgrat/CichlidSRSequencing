# The four variants genotyped by sexloci: detection, genetic change, and genotyping against GT3a

Handover document for sexloci 0.1.26 (status 8 November 2026). It covers the four loci that `sexloci sv-genotype` genotypes:

1. the **gsdf tandem duplication** (LG7),
2. **TE1** and 3. **TE2**, two transposon insertions on LG10, 216 bp apart,
4. **chr7ins**, a transposon insertion directly upstream of gsdf (LG7).

For each: what is known about the genetic change, how it was detected in this project, and how to genotype it against the GT3a reference.

## 0. Conventions, and how to read the claims

* **GT3a** is the reference `Mzebra_GT3_NCBI` (GCF_041146795.1, *Metriaclima zebra*). All GT3a coordinates are 1-based. LG7 is `NC_135173.1`, LG10 is `NC_135176.1`.
  **AC** is the *Astatotilapia calliptera* reference fAstCal1.2, in which Munby et al. report their coordinates.
* Insertions use the VCF anchor convention: `POS` is the anchor base, `REF` is that base, and `ALT` is the anchor base followed by the inserted sequence.
  sexloci compares the anchor base with the FASTA and stops on a mismatch.
* Evidence labels follow the lab's summary page. **established**: replicated in a pedigree or by an independent method. **strong**: clear statistical
  support, small sample or no replication. **suggestive**: fits, but could be chance or an artifact. **open**: not yet tested.
* **NOT RECORDED** marks something this repository, and the sessions that built it, do not document. Each has an action in section 7.
* Statements about the literature come from the papers listed in Appendix B. Statements about lab data come from the analyses run while building
  sexloci (20 October to 8 November 2026). Where a statement is an inference, it says so.

## 1. At a glance

| | gsdf duplication | TE1 | TE2 | chr7ins |
|---|---|---|---|---|
| Chromosome | LG7, NC_135173.1 | LG10, NC_135176.1 | LG10, NC_135176.1 | LG7, NC_135173.1 |
| GT3a position | unit 50,751,322-50,794,777; the two copies join near 50,773,050 | 24,921,178 | 24,921,394 | 50,765,050 (copy 1) and 50,786,787 (copy 2) |
| Genetic change | head-to-tail duplication of a 21.7 kb unit containing gsdf; **GT3a already carries two copies** | insertion, +1,139 bp | insertion, +9,868 bp | DNA transposon (CACTA-type), ~3.8 kb, 2 bp target-site duplication |
| Detected by | MC pedigree (LG7 XY) and read depth, then self-alignment of GT3a; the variant was already known from Munby et al. | supplied by the lab as a VCF-style record; **calling method NOT RECORDED** | as TE1 | reported by Munby et al. in AC coordinates; located in GT3a by a carrier-versus-control junction-read scan |
| Genotyped by | depth over the unit relative to its flanks, plus a junction check | read pairs crossing the junctions | read pairs crossing the junctions | read pairs crossing the junctions at both copies, summed into one call |
| Sex link | Y allele in MC (*Mchenga conophoros*) and in Masoko *A. calliptera* | male-linked in YH, saulosi and gertrudae | male-linked in YH and saulosi | Y allele in *A. calliptera* and some mbuna (Munby) |
| Standing | established as male-limited in MC; cause strong | association strong, cause open | association strong, cause open | genotyping validated against Munby's calls; role as in Munby |

Two of the four sit at the same gene: gsdf is raised by the duplication and, independently, by chr7ins (Munby et al. describe both as Y alleles that increase gsdf expression).
chr7ins lies inside the duplicated unit, which is why GT3a has two copies of its site.

## 2. Two genotyping methods shared by all four

All four are genotyped by `sexloci sv-genotype`, which streams from each sample's cloud BAM only the regions it needs (one `samtools view` request per
sample: the 130 kb gsdf region plus a window of flank + 200 bp on each side of every insertion; overlapping windows merge). Nothing is downloaded whole
and nothing is written to the cloud. Only samples that have an alignment to GT3a in the AlignmentDatabase can be genotyped.

### 2.1 Depth method (the gsdf duplication only)

Reads from every genomic copy of the unit map to both reference copies, so a sample's depth over the unit, relative to its own flanks, is proportional
to its total copy number in the diploid genome:

    ratio = copies / 4        (2 copies 0.50, 3 copies 0.75, 4 copies 1.00 = the reference's own state)

Constants are in `gsdf.py`. Region read: `NC_135173.1:50,700,000-50,830,000`. Depth is `samtools depth -a` (duplicates excluded by default), averaged
in 5 kb bins. Unit depth is the mean of bins starting in 50,755,000-50,785,000. Flank depth is the mean of bins in 50,700,000-50,735,000 and
50,800,000-50,815,000, which lie outside the block and before a repeat-rich bump at 50.825 Mb. Then

    EstCopies  = 2 + (ratio - 0.530) / 0.268          (0.530 and 0.268 calibrated on the MC family)
    CopyNumber = round(EstCopies), clipped to 1-6
    -/- = 2 copies,  dup/- = 3,  dup/dup = 4,  CNn = anything else

A call is flagged `low depth` (flank depth under 3x), `uncertain` (EstCopies more than 0.35 from an integer) or `outside calibrated range`
(EstCopies below 0.6 or above 5.4). Alternatives to the default calibration: `--theory` (ratio = copies/4) and `--anchors FILE` (fit on samples of known
copy number, which must include two different states).

A second, independent measure is `JunctionJ`: uniquely mapped coverage (MAPQ of at least 30) across the junction between the two copies,
`50,772,650-50,773,150`, divided by unique flank coverage. Reads there exist only on chromosomes that carry the tandem duplication.

### 2.2 Junction-pair method (TE1, TE2, chr7ins)

Short reads can decide between "insertion" and "no insertion" only where a read crosses from flank into element (or back). The method, in
`longindel.py` (the lab's `IndelReadClassifier`) and `svgeno.py`:

1. For each site build two local alleles from the GT3a FASTA: the reference allele (250 bp of flank each side, `--flanking`) and the insertion allele
   (the same flanks around the inserted sequence).
2. Align every read pair to both alleles (Smith-Waterman; parasail if installed, else a slower numpy implementation with the same scores).
   Unmapped mates stored at their mate's position are kept, because they carry inserted sequence. Duplicates and secondary and supplementary alignments are skipped.
3. A pair is **informative** only if one read aligns clearly better to one allele (score margin of at least 20) and crosses a junction with at least
   15 bp on each side. Pairs that score alike on both alleles, do not cross a junction, or whose mates disagree are counted as equal or uninformative.
4. The counts of REF and ALT pairs go into a binomial model, with P(ALT pair) = 0.02 for 0/0, 0.5 for 0/1 and 0.98 for 1/1. The genotype is the most likely
   one, and **GQ** is the phred-scaled gap between the best and second-best.
   The model is one function, `longindel.call_from_counts(n_ref, n_alt)`.
5. Fewer than 5 informative pairs (`--min-informative-pairs`) gives `./.`. The Sample Database receives only calls with GQ of at least 10 (`--min-gq`).

Consequences to keep in mind. Only junction pairs inform, so no-calls are common at low coverage: about 10x is needed for 5 pairs at a single junction in a
homozygous-reference fish, and a homozygous-insertion fish has two junctions, so it gets there sooner. At about 17x with half the reads going to a second
reference copy, chr7ins heterozygotes gave 5 to 7 informative insertion pairs and homozygotes 10 to 14. Repetitive junction sequence can make pairs score equal
on both alleles (TE2 begins with a 25 bp poly-T run and TE1 with a GAATG microsatellite), which shows up as no-calls at good coverage. The method gives
neither the length of the insertion nor phase: it cannot say whether TE1 and TE2 lie on the same chromosome.

### 2.3 Running it

    sexloci sv-genotype --limit 10                         # trial run
    sexloci sv-genotype --jobs 4                           # the whole database; resumable; all built-in sites (gsdf, TE1, TE2, chr7ins)
    sexloci sv-genotype --builtin-sites chr7_gsdf_ins --no-gsdf      # one set only
    sexloci sv-genotype --jobs 4 --write-sample-database   # DRY RUN: shows the plan; add --yes to write

* **Output** `sv_genotypes.tsv`: the gsdf call (`Genotype`, `CopyNumber`, `EstCopies`, `DepthRatio`, `FlankDepth`, `JunctionJ`, `Flag`) and, per site, `<name>_GT`,
  `_ref`, `_alt`, `_GQ`, `_reads`, `_status`. chr7ins also has a combined `chr7ins_*` set (section 5.5). The Sample Database columns are `gsdf number`, `TE1`, `TE2` and `chr7ins`.
* **Caches** `gsdf_measurements.tsv` and `te_measurements.tsv` are keyed by sample and **site name**, not by the allele sequence. A later run fetches only what is
  missing. **If an allele is ever changed, delete that site's rows from `te_measurements.tsv`**, or the old results are reused silently. `--refresh` re-measures only the selected samples.
* **Pinned alleles.** Each built-in allele has a checksum that the log prints and a test pins (Appendix A), so an accidental edit fails the test suite.

## 3. The gsdf tandem duplication (LG7)

### 3.1 What is known

* **The gene and its annotation.** gsdf (gonadal soma-derived factor) is on LG7 (`NC_135173.1`; GeneID 143412435). The GT3a annotation has one gsdf model, 50,754,902-50,776,300 on the minus strand.
  That model is 21.4 kb, almost exactly one repeat period, and straddles the junction between the copies, so it looks like a merged model of both copies.
  Distances to "gsdf" computed from its coordinates are not meaningful.
* **The variant.** A head-to-tail tandem duplication of a 21.7 kb unit that contains the whole gene. Munby et al. (2021) identify it as one of two Y alleles on chromosome 7
  that raise gsdf expression in Lake Masoko *A. calliptera* (the other is chr7ins), and report it in much of the Lake Malawi haplochromine radiation. In their Table S5 it is absent from
  *Diplotaxodon* and *Rhamphochromis*.
* **GT3a already has two copies.** Self-alignment of the reference gives a unit at `NC_135173.1:50,751,322-50,794,777`, period 21,728 bp, about 99.4% identical, with the copies
  joined near 50,773,050, no spacer of note between them and nothing inverted. GT3a is therefore in the "four copies" state, and samples are measured against it.
* **Y-linkage in the MC pedigree (*Mchenga conophoros*; established).** All 24 sons carry one extra copy and none of 24 daughters do (p = 3e-9). The XY region on LG7 is 50.25-53.0 Mb.
  The male-minus-female depth difference spans about 50,750-50,795 kb, so the extra copy on the Y covers the whole unit, not just the gene.
  In the database, 27 of 28 *M. conophoros* males carry it (26 with one extra copy, 1 with two) against 0 of 43 females; the one male without it is unexplained.
* **Elsewhere in the lake.** Across the database 111 of 423 sexed males carry the duplication against 4 of 165 females. YH (63 animals) and saulosi (24 animals) all have two copies, so gsdf is not involved in their LG10 systems.
  Two gertrudae males (AG_563_1, AG_619_1) carry it, both also TE1-positive; that is too few to test.
* **Not known.** Whether the extra copy on the Y has the same boundaries as the reference's duplication (this needs long reads or an assembly of a carrier). The AC allele's breakpoints have not been compared
  with GT3a's. Munby's coverage panel suggests an elevated block about 21-22 kb wide, which agrees in size only.

### 3.2 How it was detected in this project

1. **Mapping.** sexloci's pedigree and association analysis of the MC family placed an XY signal on LG7 at 50.25-53.0 Mb (tier A; both broods significant). SNP calls inside the block are unreliable, because reads from two copies pile onto one,
   and the strongest clean SNPs lie just outside it (50,743,244, 7.8 kb before copy 1, and 50,807,442, 12.7 kb after copy 2). That makes a structural variant the more likely sex-linked variant (an inference).
2. **Depth.** Read depth over the region, relative to flanks, was 0.530 in daughters (SD 0.027) and 0.799 in sons (SD 0.063): a step of 0.268 against the 0.25 expected for one extra copy of four.
3. **Read pairs.** Pairs spanning the two copies had spacing of 21,577-21,650 bp, which gave the tandem arrangement and a period near 21.7 kb.
4. **Self-alignment of GT3a.** `scripts/detection/selfdot.py` draws a k-mer dot plot of `NC_135173.1:50,700,000-50,830,000` (500 bp windows) and reports the period and both copies. It found one line parallel to the main
   diagonal and nothing inverted. Exact boundaries come from `minimap2 -X -c` of the region against itself, which gave the 50,751,322-50,794,777 span.
5. **Junction reads.** Uniquely mapped reads across the junction (section 2.1) exist only on duplication carriers, and confirm the head-to-tail arrangement.

### 3.3 Genotyping against GT3a

Use the depth method of section 2.1. In practice: run `sexloci sv-genotype` and read `Genotype` / `CopyNumber` (and the `gsdf number` column), ignoring flagged calls unless you have checked them.
Depth gives **total copy number, not which chromosome carries the extra copy**; calling a duplicated Y needs sex or pedigree. Use `JunctionJ` as a cross-check: it should be clearly above zero in 3-copy
and 4-copy samples and near zero in 2-copy ones.

### 3.4 How well it works

* MC family (51 samples): all 46 unflagged calls match the expected state (females 2 copies, males 3). The two wrong calls (MC_049 called 2, MC_060 called 4) are both flagged.
* These numbers use the same samples as the calibration, so they flatter the method. Calibrating on only three samples still gave 50 of 51 right. With no calibration (`--theory`), females come out at 2 and males at 3-4.
* Copy numbers agree with Munby et al. Table S5 for 98.9% of 629 shared fish.

### 3.5 Limits

* Needs depth: about 5x in the flanks at least (three samples below 3x were flagged). Telling 3 from 4 copies is harder than 2 from 3; wait for about 8x before trusting a 4.
* The default calibration (0.530, 0.268) is for the lab's pipeline (bwa to GT3a, 4-8x). For another pipeline calibrate with `--anchors`.
* Partial deletions inside the unit could confuse a depth-only call. The labels assume each chromosome carries one or two copies.
* Junction-PCR and qPCR assays were designed (`gsdf_primer_tools.py`, `data/gsdf_primers_MC050.tsv`) but **have not been tested at the bench**.
* *A. calliptera* fish have no GT3a alignments in the database, so they cannot be genotyped here.

## 4. TE1 and TE2 (LG10)

### 4.1 What is known

* **The records.** TE1 is an insertion of 1,139 bp at `NC_135176.1:24,921,178` (anchor base `t`). TE2 is an insertion of 9,868 bp at `NC_135176.1:24,921,394` (anchor base `a`).
  They are 216 bp apart in GT3a. The full inserted sequences are in `known_sites.py` (Appendix A gives checksums).
* **The LG10 sex-linked block.** In YH (the lab's code for *Aulonocara* "Nkhata Yellow Head", which is *Aulonocara* sp. 'chitande type north' Nkhata Bay in Kumar et al., the species behind their chr10 inversion work; their YH_GT1 assemblies) males are heterozygous and females homozygous across 11.6-30.0 Mb (18.4 Mb, 8,679 diagnostic SNPs).
  That matches the roughly 18 Mb chr10 inversion that Kumar et al. report in the same species, whose cross had 47 of 48 offspring fitting XY by inversion genotype. TE1 and TE2 at 24.9 Mb lie inside this block.
* **Which haplotype is the Y.** In Kumar et al., females were homozygous for the inversion and males heterozygous, and the GT3a orientation is the ancestral one. So the Y carries the non-inverted (GT3a-orientation) haplotype, as Blumer et al.
  also describe for chromosome 10. TE1 and TE2 are absent from YH females, so in GT3a coordinates the insertions are on that non-inverted Y haplotype (an inference from those two sets of results).
* **Three *Aulonocara* groups, one locus** (all called against GT3a; the male-specific LG10 regions differ):

  | Group | LG10 region | TE1 | TE2 | Standing |
  |---|---|---|---|---|
  | YH (two families, 52 pedigree animals) | block 11.6-30.0 Mb | 21 of 24 males, 0 of 17 females (p = 7.5e-9) | 22 of 26 males, 0 of 20 females (p = 1.9e-9) | XY on LG10 established (tier A); insertions travel with the Y haplotype (strong) |
  | gertrudae (21 animals, no pedigree) | narrow region 24.0-26.5 Mb | 11 of 11 males, 1 of 7 females (p = 3.8e-4) | absent from all 14 sexed animals | tier B, the ceiling for unrelated animals; strong |
  | saulosi (24 animals, no pedigree) | block 13.5-25.0 Mb, interval at 24.8 Mb or beyond | 8 of 8 males, 0 of 9 females (p = 4.1e-5) | 8 of 8 males, 0 of 6 females (p = 3.3e-4) | tier C by window statistics; strong at animal level |

  In YH the two insertions agree on carrier status in 41 of 42 animals called for both, so association cannot separate them. In saulosi the block shares no male-specific alleles with YH's (48% agreement against 51% expected by chance), and its edges differ.
* **Two saulosi recombinants** (ASG_588_1, a male, and ASG_209_1, a female) place the sex-determining interval at 24.8 Mb or beyond; the interval contains TE1 and TE2. The recorded sex follows the insertion genotypes, not the block haplotype (strong, but two animals).
* **Across the lake.** TE1 occurs in 26 species (21 of them Deep Benthics) and TE2 in 12. TE2 was never seen without TE1: 28 animals in 12 species have TE1 without TE2, and the 13 TE2 carriers without a confident TE1 call all have at least one TE1 insertion pair.
* **Counter-examples** that keep this short of a causal claim: YH_004_m, YH_009_m and YH_014_m (19-29x) carry neither insertion, and YH_069, the one pedigree male that looks XX at the diagnostic SNPs, lacks TE2 (consistent with Kumar et al.'s single XX male). A gertrudae female, AG_567_1, carries TE1 (2 reference and 5 insertion pairs, GQ 13) on an X-like haplotype.
  In *A. stuartgranti* the carriers include all three confidently called females, the opposite of the saulosi pattern. Whether the insertions are the sex determiner or passengers on the Y haplotype is **open**.

### 4.2 How they were detected: NOT RECORDED

TE1 and TE2 were **supplied to sexloci by the lab** as two VCF-style records, together with the lab's earlier genotyping script (now `longindel.py`). The repository does not record how they were first called: the assembly, the
caller and its parameters are not documented here. What the records themselves show:

* Format: tab-separated, no header, lower-case bases, `QUAL` 100, `FILTER` `*`, a quoted `INFO` with `TYPE=INSERTION` and `ALN_START`, and a trailing column `Special`. That is typical of an alignment-based call of a long-read assembly against the reference (an inference from the format).
* `ALN_START` is 197,394 for TE1 and 198,749 for TE2. The difference is 1,355 bp, which equals TE1's insert (1,139 bp) plus the 216 bp between the sites in GT3a. **If `ALN_START` is a coordinate on the assembly contig, both insertions lie on one assembled haplotype, in the order TE1 then TE2**
  (an inference from this arithmetic, not a recorded fact; check it against the source assembly). If so, the assembly shows the two in cis, which the short-read genotypes cannot (section 2.2).
* Sequence properties computed from the stored alleles (not from an annotation): TE1 is 51% GC and begins with a GAATGAATGAATGAATGAATGC microsatellite; TE2 is 41.5% GC, begins 10 bp after the anchor with a 25 bp poly-T run and contains a 32 bp poly-A run.
  Neither has an exact terminal inverted repeat or an exact direct terminal repeat, so neither looks like the CACTA-type element of chr7ins. **Their transposon family has never been assigned.**

### 4.3 Genotyping against GT3a

Use the junction-pair method (section 2.2) with the built-in set `deep_benthic_TE`. TE1 and TE2 are separate sites, each with its own local alleles, and their fetch windows overlap and are merged.
Because they are 216 bp apart, a fragment can in principle span both, which is the basis for phasing them (section 7). `clip_check.py` helps interpret odd calls: soft-clipped reads at the site with no insertion pairs point to a **different** insertion or rearrangement at the same position, not the one sexloci looks for.

### 4.4 Limits

* No-calls are common (section 2.2). The calls used in the summaries above are GQ of at least 10, and only samples with a GT3a alignment are genotyped.
* The Y haplotype in YH comes from only two sires. The saulosi block is not known to be an inversion (short reads cannot tell it from a diverged XY region); HiFi sequencing is the planned test.
* Genotypes cannot say whether TE1 and TE2 are on the same chromosome.

## 5. chr7ins, the insertion upstream of gsdf (LG7)

### 5.1 What is known

* **Munby et al. (2021)** found it in Lake Masoko *A. calliptera* as a third Y-like signal on chromosome 7: an insertion at `chr7:18,086,980` in AC, about 2.5 kb upstream of gsdf, exclusive to males, independent of the duplication, reported as about 5 kb.
  Two subregions of 638 bp and 510 bp occur at 19 and 18 places in the AC reference, and RepeatModeler assigns them to the two ends of one unknown repeat family. Their assembled contigs matched either end of a **3,947 bp** unknown transposable element at 94% and 97% identity, so the "5 kb" is probably a rough estimate and nobody assembled the whole insertion.
  They genotyped it by mapping reads to the assembled insertion with 1 kb flanks, calling 0/1 when both alleles each had at least 10% of reads.
* **Carriers in their Table S5.** chr7-ins: 108 carriers (82 het, 26 hom), in *A. calliptera* (74) and mbuna (34, at Nkhata Bay and Chilumba). At Lake Itamba all 32 sampled males carry it (22 hom, 10 het) with no gsdf duplication, and no females were sampled.
* **What we reconstructed (section 5.4).** A CACTA-type DNA transposon, 13 bp exact terminal inverted repeat (`CACTCAGAAAAAT ... ATTTTTCTGAGTG`), 2 bp target-site duplication `TG`, estimated length 3,841 bp.
* **Not known.** The true length and interior (3,841 bp, Munby's 3,947 bp and the "5 kb" disagree, and the interior of a carrier's element has not been sequenced); whether it is causal; its effect on gsdf expression in the lab's fish.

### 5.2 Locating it in GT3a

* **Coordinate lift.** gsdf is on opposite strands in AC and GT3a, so the lift reverses orientation. Matching five gsdf exons gives, at this locus, `GT3a = 68,852,031 - AC` (for example AC 18,090,125 maps to 50,761,906, an exact GT3a exon end).
  The insertion is 2,564 bp upstream of the first coding exon in AC and so maps to **NC_135173.1:50,765,051** in copy 1. (An earlier estimate from the AC duplication boundary was 172 bp off; do not use it.)
* **Two copies.** Copy 2 is copy 1 plus 21,728 bp: predicted 50,786,779, observed junction at 50,786,786.
* **Junction scan.** The clip boundaries found by `find_insertion.py` (section 5.3) were 50,765,052 (right clips) and 50,765,049 (left clips), within 1 and 2 bp of the lift.

### 5.3 How it was detected in GT3a

`scripts/detection/find_insertion.py` looks for an insertion that the reference lacks. An insertion leaves soft-clipped reads on both sides of its position (right clips at the first junction, left clips at the second), with boundaries overlapping by the target-site duplication.
It streams the 130 kb gsdf window from each carrier and control, counts clip boundaries per sample, keeps right/left clip pairs that are supported in carriers and (nearly) absent from controls, merges the same insertion at several loci, and assembles the two ends.

* **Carriers and controls.** Two Kingiri carriers from Munby's calls (SAMEA117806260, homozygous; SAMEA117806261, heterozygous) against eight non-carriers with GT3a alignments (AG_552_1, AG_558_1, AG_584_1, ASG_160_1, ASG_195_1, YH_006_f, YH_010_f, YH_013_f).
* **Result.** Candidate #1 was the only plausible insertion: the same insertion found at both reference copies and merged; 0 of 8 controls with meaningful clipping; 15 supporting carrier reads, 3-6 per side per carrier; anomalous mates 34 per carrier against 6 per control.
  Candidate #2 was rejected (most reads span it; controls had more anomalous mates than carriers). Candidate #3 was rejected: its "duplication" was 40 bp and identical on both sides, a short tandem repeat, and 2 of 8 controls showed clipping.
  An earlier candidate (paired breakpoints at 50,744,293 and 50,744,301) was also withdrawn: it was about four independent fragments with local mates and a G/A-rich stretch.
* **Pilot genotyping.** The junction sequences of candidate #1 (a 124 bp stub allele) were run through `sv-genotype` on 47 fish with Munby's calls: AUC for carriers against non-carriers 0.97 (permutation p = 0.0001), 15 of 16 carriers with at least 2 insertion pairs, 0 of 24 non-carriers, and 702 fish outside *A. calliptera* all `0/0` or `./.`.

### 5.4 Identifying and reconstructing the element

1. **Family.** `mates_vs_repeatmasker.py` takes the carriers' reads whose mates map far away or to other chromosomes and looks up the GT3a RepeatMasker annotation at the mates' positions. Mates of the left junction fall in `rnd-1_family-458` (863 bp consensus), mates of the right junction in `rnd-1_family-764` (1,240 bp consensus). The element was classified DNA/CMC-EnSpm (CACTA-type) from this.
   The split is clean in both carriers and absent in the control AG_552_1. It is consistent with Munby's two end regions, which occur at 19 and 18 places in AC.
2. **Intact copies.** `element_from_annotation.py` pairs family-458 and family-764 hits a few kb apart in GT3a. The body between the ends in every modal copy is `rnd-1_family-940` (Unknown, 1,632 bp) and `rnd-5_family-1892` (Unknown, 787 bp); longer copies (5-6.5 kb) also contain a nested Tc1 element of about 1.6 kb.
   Modal-length intact copies are about 3.9 kb (median of the class 3,896 bp), which fits Munby's 3,947 bp and not "5 kb". 25 modal copies were kept (`element_modal_copies.fasta`).
3. **Check against the carriers.** `template_check_with_mates.py` aligns the carriers' mates (default: Munby's 16 carriers, both reference copies) to a template, gets median 99% identity in the covered ends, and applies the 9 substitutions and 2 indels (1 bp deletion at template position 235, 1 bp insertion at 333) that the carriers consistently show, giving a **3,841 bp final estimate**.
   All 25 modal GT3a copies differ from the carriers by 14-16 bases per kb, so no GT3a copy is the source. Only about 500 bp at each end can be checked: **the interior is that of the template, not of the carriers.**

### 5.5 Genotyping against GT3a

Short reads can only use the junctions, so the allele carries only what a read can use. `build_chr7ins_alleles.py` writes, at each of the two reference copies, an insertion allele made of the anchor base, the first 300 bp of the element, the last 300 bp, and the duplicated `TG` (603 bp; the anchor is the G of the `TG` duplication).
Both sites share the same allele (checksum in Appendix A):

| Site | Position | Anchor |
|---|---|---|
| `chr7ins_c1` | `NC_135173.1:50,765,050` | G |
| `chr7ins_c2` | `NC_135173.1:50,786,787` | G |

A carrier has one gsdf unit per chromosome and GT3a has two, so its reads split between the two copies. sexloci therefore genotypes both and **sums the informative pairs into one combined call** (`SITE_GROUPS` in `known_sites.py`): the columns `chr7ins_GT`, `_ref`, `_alt`, `_GQ`, `_reads`, `_status`,
and the Sample Database column `chr7ins`. The combined call needs both sites; renaming them with `--indel-names` switches it off. In the lab's pipeline the combination equals `scripts/detection/chr7ins_combine.py`.

### 5.6 Validation (47 fish with Munby's calls)

* All 16 of Munby's carriers show the junction and none of his 24 non-carriers does (4 of those non-carriers are from the carriers' own population, Lake Kingiri); AUC 1.00 (permutation p = 0.0001).
* Genotypes: 3 of 3 homozygotes `1/1`, 12 of 13 heterozygotes `0/1`, 24 of 24 non-carriers `0/0`, and 7 of 7 fish without a Munby call `0/0`. Homozygotes carry a median of 11 insertion pairs against 6 for heterozygotes (mean insertion fraction 0.97 against 0.44).
* Outside *A. calliptera*, all 874 confident calls among 1,043 aligned fish (benthics, *Diplotaxodon*, *Rhamphochromis*, utaka, hybrids) are `0/0`.
* The one discordance, SAMEA117806076 (Munby het, 78:19; ours 63:6, `0/0`), has three gsdf copies. See 5.7.

### 5.7 Limits

* **Copy number.** The model expects half of the informative pairs from the insertion in a heterozygote. A fish with three or four gsdf copies has extra reference-like copies of this locus that dilute the insertion allele (about a third at three copies). SAMEA117806076's insertion fraction (9%; Munby's own ratio 20%) is below even a third, so the dilution probably contributes but does not fully explain it.
  Treat a `0/0` with a few insertion pairs in a fish with more than two gsdf copies as uncertain. A copy-number-aware model is not implemented.
* **Only fish aligned to GT3a are genotyped.** Munby's Itamba and Masoko carriers and most mbuna carriers have no GT3a alignment, so the lake-wide frequency, and its presence in mbuna, have not been tested here.
* **Length and interior** are unresolved (5.1). Flank-to-flank PCR across the insertion in a homozygote would settle 3.9 kb against 5 kb (69% of Itamba males are homozygous); long reads from one homozygote would give the actual element. Neither has been done.

## 6. Related locus not included: chr19-ins

Munby et al. report a third Y allele, a roughly 700 bp insertion on chromosome 19 that matches an LTR/Unknown element (97% identity, 99% coverage), at `chr19:21,572,413` in AC, 1.7 kb upstream of id3. In their Table S5 there are 89 carriers (86 het, 3 hom), all *A. calliptera*, 65 of them in Lake Masoko.
It is **not** built into sexloci: its GT3a position has not been located. The next step is `grep -i -w -E "id3|e2f2" genomic.gff` on the GT3a annotation to map the nearby genes, then the same junction-scan approach as section 5.3.

## 7. Reproducing the detection, and open items

The scripts that did the detection are in `scripts/detection/` (see its `README.md` for the order, inputs and outputs). They are standalone, were written for the lab laptop (paths under `~/Temp/SexDeterminationAnalysis`, the lab's `ptm_dropbox:` remote, `samtools`, `rclone`),
and are **not covered by the test suite** apart from a check that each file parses.

| Variant | Detection and validation steps | Scripts |
|---|---|---|
| gsdf duplication | self dot plot, then depth and junction genotyping; assembly check; primers | `selfdot.py`, `gsdf_cn.py`, `gsdf_dup_check.py`, `gsdf_primer_tools.py` |
| TE1, TE2 | genotyping only (the calls came from the lab); soft-clip check | `clip_check.py`, `cov_check.py` |
| chr7ins | junction scan, family lookup, element reconstruction, allele build, validation | `find_insertion.py`, `mates_vs_repeatmasker.py`, `element_from_annotation.py`, `template_check_with_mates.py`, `build_chr7ins_alleles.py`, `chr7ins_combine.py`, `chr7ins_concordance.py`, `alignment_check.py` |

Open items, most useful first:

1. **TE1/TE2 provenance.** Record the assembly, caller and parameters that produced the two records, and check in the assembly that `ALN_START` really places them in cis (section 4.2).
2. **TE1/TE2 family.** Assign it (RepeatMasker annotation at the inserted sequences, as was done for chr7ins) and look for a target-site duplication in the source assembly.
3. **Phase TE1 and TE2** with read pairs that bridge both sites (216 bp apart).
4. **chr7ins reach.** Align Munby's Itamba, Masoko and mbuna carriers to GT3a and rerun `sv-genotype` to test the lake-wide distribution; add copy-number-aware calling; settle the length by PCR or long reads.
5. **gsdf.** Test the junction PCR and qPCR assays at the bench; compare the Y-linked allele's boundaries with the reference's using long reads.
6. **chr19-ins** (section 6), and the saulosi HiFi test of the LG10 block (section 4.4).

## Appendix A. Allele fingerprints (pinned by tests)

| Site | Anchor | Inserted | Allele fingerprint |
|---|---|---|---|
| TE1 | `NC_135176.1:24,921,178` t | 1,139 bp | ALT 1,140 bp, md5 1957b81b6617 |
| TE2 | `NC_135176.1:24,921,394` a | 9,868 bp | ALT 9,869 bp, md5 d2b75a3b30e2 |
| chr7ins_c1 | `NC_135173.1:50,765,050` G | 602 bp (junctions only) | ALT 603 bp, md5 8b7b5cf1d8c7 |
| chr7ins_c2 | `NC_135173.1:50,786,787` G | 602 bp (junctions only) | ALT 603 bp, md5 8b7b5cf1d8c7 |

The chr7ins allele begins `GCACTCAGAAAAATAAAGGT` and ends `CCATAATTTTTCTGAGTGTG`. `sexloci write-builtin-vcf out.vcf --set all` writes every site as a VCF.

## Appendix B. Sources

* Munby et al. (2021), bioRxiv, "Differential use of multiple genetic sex determination systems in divergent ecomorphs of an African crater lake cichlid", doi 10.1101/2021.08.05.455235 (abstract and Table S5 used here).
* Kumar, Cooper, Kocher, Streelman and McGrath (2025), eLife 14:RP104923, "Large inversions in Lake Malawi cichlids are associated with habitat preference, lineage, and sex determination".
* Blumer et al. (2025), Science 388:eadr9961, "Introgression dynamics of sex-linked chromosomal inversions shape the Malawi cichlid radiation".
* The lab's analyses while building sexloci 0.1.0-0.1.25, summarised in the "Sex determination in YH, MC, gertrudae and saulosi" page (status 1 November 2026).
