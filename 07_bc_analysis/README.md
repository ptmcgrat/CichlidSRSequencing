# 07 Backcross genotyping, recombination, TD, QTL mapping, fine mapping

Specification: [`SPEC.md`](../SPEC.md), module 07. Package: `bc_analysis` (the number is only for ordering folders).

```
python -m bc_analysis bc_genotyping provenance|flagstat|original|v2-sites|v2-counts|v2-hmm|validate|compare
python -m bc_analysis genotype sites|counts|hmm|compare|simulate        # the new genotyper
```

Tests: `python tools/run_tests.py 07` (or `python -m bc_analysis selftest`). Configuration: `config/`. Environment: `../envs/`.

Current state: the NEW genotyper (`python -m bc_analysis genotype`, below) replaces the original `GenotypeBCs.py` approach for the BC fish; the earlier `bc_genotyping` steps (original, v2) stay until the new one has run on real data; the analyses after the genotyping are not written.

## The new genotyper: `python -m bc_analysis genotype` (increment 1, 2026-11-12)

**What it replaces.** `GenotypeBCs.py` (the original, run on another reference): markers by the YH sire and the brood's MC mothers, every read one observation of a two-state HMM with a switch probability of 1e-22, Viterbi, interpolation to a 100 kb grid.
What was found in it: the one uncommented parent filter (`GT[4]="AA" && GT[1]="RR" & ...`) selects MC-5G11G, not the YH sire, and mixes `&` with `&&`, so the parents VCF of the original run cannot be reproduced from it; the pileup has `-Q 0` and no `-q`; the three MC mothers have no depth requirement;
`AD[:,:,0:2]` assumes the order of the alleles; the shuffle is unseeded; the final inner join cuts the grid at the shortest fish; the F1 sires are never used; the alleles are read as REF = MC, which was true only on the old reference (on Mzebra_GT3 the YH allele is the reference allele at about half the sites,
and the unchanged script then calls at chance, accuracy 0.50 in the simulation). **Kept, as decided by the PI:** sites are the alleles alternatively fixed in the YH father and the MC mothers of THAT brood (a heterozygous site in an MC mother who is not the dam of the brood does not matter), the resolution is 100 kb, and short interior
recombinant segments are treated as artifacts (repeats, calling error), not as gene conversion: as a setting that can be switched on (below).

**Steps.** `genotype sites --vcf <bc_family VCF>`: per dam, the sites where `YH_008_m` is homozygous for one allele and `MC-1R6R-f` (the F1 sires' mother) and the dam are homozygous for the other; depths relative to each parent's median (the old `> 12` is above the 9 to 13x effective depth of the new alignments), a purity rule, polarity kept per site.
The parents VCF comes from module 05: `python -m variant_calling cohort build bc_family` (the seven parents, bcftools, per sample). `genotype counts [--run]`: for each BC fish whose dam is configured (the four "BC1 hybrid" fish of another cross are left out), the counts of the MC and the YH allele at its dam's sites, read BY ALLELE NAME from a
`bcftools mpileup -T` (mapping quality 20, base quality 20, duplicates excluded, indels skipped); the BAM is downloaded, counted and deleted; a fish is counted again only when the definition key or its BAM changed. `genotype hmm`: reads per site capped, summed per 100 kb bin, an HMM over the bins with a switch hazard of `rho` (1 cM per Mb),
the error rate and the reference bias fitted per fish by EM and applied by each site's polarity; calls, posteriors, breakpoints with intervals and a flag for terminal calls under 0.5 Mb, per-fish QC, and `end_support.tsv` (the reads in the first and last megabase of every chromosome and whether the end calls are confident, so that you can see which ends are callable at all).
`genotype compare --b <folder>`: what a second run changed. `genotype simulate`: the comparison below.

**The first real run: no interior limit, terminal minimum 0.2 Mb** (the PI's decision: measure how large the problem is before limiting anything). `hmm` then lists every interior segment shorter than 2 Mb in `short_segments.tsv` and prints how many there are per 100 fish-chromosomes; that is what a minimum would remove.
The sensitivity run is `genotype hmm --interior-min-mb 2 --tag _int2` and then `genotype compare --b hmm_int2`. An interior minimum is a penalty of about 70 nats, not a ban. A terminal minimum STRETCHES a strongly supported short terminal segment to the minimum and drops a weakly supported one.

**Why a minimum for interior segments but not for the ends.** A flip-flop needs two switches, so a penalty on short INTERIOR segments removes it; one switch near an end is an ordinary crossover and must stay cheap. A switch penalty of 1e-22 treats both alike. The simulation (`genotype simulate`; 480 chromosomes per coverage; markers sparse near the ends, mixed polarity, subtelomeric artifacts; 40% of the crossovers near an end)
gives, as recall of the true crossovers by distance to the nearest chromosome end (the number of crossovers in brackets):

| | under 0.5 Mb | 0.5 to 1 Mb | 1 to 2 Mb | 2 Mb and beyond |
|---|---|---|---|---|
| 0.3x, the original's model with the alleles named | 0.04 (23) | 0.25 (36) | 0.65 (54) | 1.0 (187) |
| 0.3x, bins, no limits | 0.74 | 0.97 | 1.0 | 1.0 |
| 0.3x, bins, first-run setting (terminal 0.2 Mb) | 0.70 | 0.97 | 1.0 | 1.0 |
| 0.3x, bins, interior limit 2 Mb | 0.74 | 0.97 | 1.0 | 1.0 |
| 0.5x, the original's model with the alleles named | 0.00 (24) | 0.56 (32) | 0.97 (34) | 1.0 (219) |
| 0.5x, bins, no limits | 0.67 | 1.0 | 1.0 | 1.0 |
| 0.5x, bins, first-run setting (terminal 0.2 Mb) | 0.58 | 1.0 | 1.0 | 1.0 |
| 0.5x, bins, interior limit 2 Mb | 0.67 | 1.0 | 1.0 | 1.0 |

Flip-flops per 100 chromosomes: the original 0; bins with no interior limit 20 to 22 (the first-run setting 20 to 22); with the interior limit 0. Spurious terminal calls per 100 chromosomes: the original 0; bins 4 to 10 (the terminal minimum changes that by less than 1). Bin accuracy 0.998 to 0.999 for all bin variants. At an even marker density (no sparse ends, no end artifacts) the original loses only
about a third of the crossovers in the last 0.5 Mb at 0.5x (0.67 against 0.93). **These numbers rest on assumptions about the real data** (marker density near the ends, how many repeat artifacts there are and how strong, how many crossovers lie near an end); the real data will say, and `short_segments.tsv` and `end_support.tsv` are made to say it.

**Tests: 57 in `tests/test_genotyper.py`** (the module in all 58): the sites on a planted parents VCF (one site of every kind), the counts (parsing by allele name, a fake tool chain, resumability), the HMM against ENUMERATION OF ALL PATHS for eight parameter sets, the limits in bins, the likelihood of the four read classes, the simulation, and the commands end to end (a planted 400 kb artifact appears in `short_segments.tsv`
without the limit and vanishes with it; `compare` reports 4 bins, 2 breakpoints and 1 segment). **Mutation checks: 41 mutations, all caught in the end;** six survived at first and each exposed a missing test (the reference-bias term for sites where the YH allele is the reference, the fraction rule of the purity check, the absolute depth floor, a failed fish counting as current, and two simulation variants that did not use the settings in their names).
**Bugs found while building, and fixed:** the terminal limit did nothing when the interior limit was off (the age dimension was capped at the interior limit and the left end borrowed the interior hazard: the combination of the first run); an off-by-one in both limits; the LG label looked at the wrong characters; `choose_contigs` gives no end for whole chromosomes. `bc_analysis/cli.py` had no `import glob` (used by `publish_bc`).

**Not exercised:** a real parents VCF, a real BAM, `bcftools mpileup -T` with a plain tab-delimited targets file (the command is built and its output parsed from canned text, not run), the run time and memory of `counts` and `hmm` on 267 fish, and the figures. **Not decided:** whether the QTL step takes hard calls at 100 kb (as before, `genotypes_calls.csv` keeps the layout of `OutputGenotypes.csv`) or the posteriors (R/qtl2 or own code, SPEC open decision);
the thresholds of the QC flags (`flag_max_eps`, `flag_max_switches_per_chromosome`) are guesses until a real run.

### Increment 2 (2026-11-12, later): backcross broods only, and what the first real sites look like

**Backcross broods only.** `genotype counts` took every fish whose dam is one of the three configured mothers; on the real sample table that included the 8 fish of `MCYH-F1-23-02-21` (dam MC-3R6R-f), an F1 brood: heterozygous everywhere, nothing to map. A fish is now a backcross fish only if its BroodID starts with a configured prefix (`sites.brood_prefixes: [MCYH-BC1-]`, the five broods of the original script), the plan shows the sire next to brood and dam, and it lists the fish of a configured dam that it leaves out.
The rule rests on the brood name, not on the sire's name, because the Parents column spells the sires in more than one way (the real plan shows which). 267 fish with an alignment remain, in five broods (48, 20, 54, 84, 61); one more backcross fish (MCYHBC1-926-1) has no alignment row yet, so 268 in all. The left-out brood has the sire YH_005_m, a pure YH male: its fish are YH x MC F1s.

**The parents VCF** (`bc_family`, built 2026-11-12 08:41 on mzebra from the seven parents' BAMs; 10,895,752 records, 9,539,050 biallelic SNPs) and `genotype sites` (135 s): 5,913,423 SNPs with YH_008_m homozygous; 476,598 of them with MC-1R6R-f homozygous for the other allele; **296,787 / 295,918 / 381,771 informative sites for the dams MC-1G2G-f / MC-3R6R-f / MC-5G11G-f**
(about 32 / 32 / 41 per 100 kb bin of a 0.93 Gb genome). **Polarity is mixed: the YH allele is the reference allele at 49%, 56% and 57% of the sites**, so the original script's reading of REF as MC would be wrong at about half the sites on this reference.
Density: interior bins have a median of 22 / 22 / 30 sites; 11.8 / 11.8 / 10.1% of them have under 5 sites and 7.3 / 7.1 / 6.7% none. The first and last megabase of every chromosome have a MEAN density equal to the interior's (104 / 96 / 97%) but are twice as often empty (none in 14.8 / 15.0 / 14.5% of the end bins); sites begin and end 0.04 and 0.05 Mb from the chromosome ends (median; worst chromosome ends 0.78 and 0.32 Mb).
**This corrects an assumption of the simulation:** scenario S2 assumed sparse markers towards the ends; the real ends are not uniformly sparse, they are patchy (about one end bin in seven is blind). They are closer to S1, so the gain at the ends over the original's model will be smaller than the S2 table above (at S1 and 0.5x: under 0.5 Mb, 67% against 77%); the larger effect on this reference is the polarity.
**Pedigree check** at the 295,918 sites of MC-3R6R-f (genotypes of MCYHF1-002 and -003 in the cohort VCF): both heterozygous at 52.5%; heterozygous in at least 66% (-002) and 74% (-003) of the sites (lower bounds: only the six most frequent genotype pairs were looked at); the other calls are homozygous, twice as often for the reference allele as for the alternative. That fits depth of 4 to 5x effective and a reference bias in the genotype calls, and no sire looks like anything but an F1;
it has not been shown to be the depth (a table of the heterozygous fraction by depth would).

### Increment 3 (2026-11-12, evening): the first real runs, and a mask for loci that recur

**Counts** (`genotype counts --run`, 263 fish after a 4-fish trial; 57.9 min with 4 jobs, no failures): 70,000 to 150,000 sites with reads per fish (a few fish with 0.8 to 2.4 million reads took 5 to 10 min each); the share of YH reads is 17 to 32%, mean about 24%, where a backcross fish should show 25%.

**The first HMM run** (no interior limit, terminal minimum 0.2 Mb; 1.0 min with 16 jobs): 5,212 breakpoints, a mean map of 88.7 cM per chromosome, 680 interior segments under 2 Mb in 245 of 267 fish (11.6 per 100 fish-chromosomes; the simulation, with artifacts as common as it assumed, gave 20 to 22). One fish flagged (MCYHBC1-809-1, 36 switches). **With an interior minimum of 2 Mb** (`--interior-min-mb 2 --tag _int2`): 4,341 breakpoints, 73.9 cM, 47 short segments left,
3,759 of 2,496,984 bins (0.15%) called differently. The fitted reference bias is about 0.5 (no reference bias to speak of).

**The artifacts are systematic, not random.** 445 of the 680 short segments lie on LG7, LG3 and LG9; the same 100 kb bins recur: LG7 7.4 to 8.2 Mb in 136 to 147 fish (a locus that shows the same signal in every fish appears as a short segment only in the fish whose neighbouring stretch is in the other state, about half of them), LG9 31.9 Mb in 89; 172 bins are hit in 5 or more fish. The segments are called with a median posterior of 0.94, and the fish with 4 or more of them do not have less coverage (115k reads against 104k):
this is a strong signal at fixed loci, not noise. At LG7 7.4 to 8.3 Mb the parents have 0.6 to 0.76 times their usual depth (YH_008_m 6.2, MC-1R6R-f 5.8, MC-3R6R-f 6.1, MC-1G2G-f 8.5, MC-5G11G-f 8.5; MCYHF1-002 2.7, -003 3.3) and 17,221 SNPs in 0.9 Mb (19 per kb, about twice the genome-wide density): not a collapsed duplication (that would be above 1.5 times), more like a hypervariable or poorly mappable region. What it is has not been established; LG3 and LG7 are not covered by module 01, LG9 has the 18.1 Mb inversion and an open suspected misassembly.
**Map length per linkage group** (2 Mb run, cM = breakpoints per fish x 100): 13 linkage groups lie between 48 and 58 cM (LG5 50, LG13 48, LG21 50, LG19 51, LG18 51, LG22 51, LG8 52, LG17 52, LG4 54, LG6 54, LG10 54, LG15 55, LG12 58), which is one crossover per chromosome (an obligate chiasma gives 50 cM in a male meiosis) and the mean of the original (48.8); nine have an excess: LG1 146, LG11 125, LG7 123, LG2 115, LG14 104, LG9 104, LG16 82, LG3 78, LG20 72. **The whole excess of the mean map (73.9 against about 52) sits on these nine.**
Breakpoints recur too: 26 bins hold a breakpoint in 13 or more fish (5%): LG14 40.7 Mb in 118 fish, LG1 42.0 in 81, LG7 7.4 in 54, LG16 0.3 and 0.2 in 51 and 49, LG9 32.7 in 43, LG11 8.3, 18.4 and 30.4 in 42, 40 and 39, LG7 9.4 in 38. A crossover hotspot of 100 kb carrying the crossovers of 44% of the meioses is not credible. **A retraction:** the earlier remark that 9.7% of the breakpoints lie within 2 Mb of a chromosome end, as uniform crossovers would give, was premature; the near-end breakpoints are concentrated (LG16 0.2 to 0.3 Mb in about 100 fish; LG14 40.7 and LG1 42.0 may be chromosome ends too: the chromosome lengths were not looked at).
So on this data the new model finds about the same map as the original on the 13 unaffected chromosomes; its gains are the polarity and the reproducibility, not a large number of extra crossovers.

**`genotype mask`** (then `genotype hmm --mask --tag _mask`): the bins that hold a short interior segment in at least 5% of the fish or a breakpoint in at least 10% (`mask:` in `config/genotyper.yaml`; `--min-segment-fraction`, `--min-breakpoint-fraction`, `--pad`), plus one bin either side, are written to `<genome folder>/mask/mask_bins.tsv` and their counts are left out of the HMM step (no recount; the EM sees no reads there either). The report lists the intervals with the most fish. Best made from a run WITHOUT the interior limit. `genotype compare --a hmm_int2 --b hmm_mask` shows what the mask changed;
whether the interior minimum is still needed after masking is what the comparison of `hmm --mask` with `hmm --mask --interior-min-mb 2` will say. The thresholds are guesses; a mask that removes much of a chromosome is a sign to look at that chromosome, not to mask more.

**Tests: 64 in `tests/test_genotyper.py`** (65 for the module); 16 more mutation checks on the mask code, all caught in the end (one survived: my docstring example of floating point noise was wrong, `0.1 * 30` is exactly 3.0; the real cases are `0.07 * 100`, `0.14 * 50`).

