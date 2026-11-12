# 06 Sex-linked regions from cohort VCFs (sexloci)

Specification: [`SPEC.md`](../SPEC.md), module 06 (decision 58). Packages: `sexloci` (vendored, unchanged) and `sd_regions` (the module's own command line and the older scripts).

```
python -m sd_regions sexloci list-vcfs | run-folder | run | check | describe | filter-audit | filter-compare | simulate     # sexloci 0.1.26, every option as in its README
python -m sd_regions sexloci run-folder --check-only                  # download and preflight every VCF in Outputs/PedigreeVCFs/, no analysis
python -m sd_regions sexloci run --vcf <local vcf> --n-perm 1000
python -m sd_regions parental_sex_region callset|depth|scan|summarize [--species YH,MC]      # the older scripts, until parity (SPEC module 06, Work 4)
python -m sd_regions sex_association fisher [--table T --out D]
```

Environment: `../envs/sex_mapping.yml` (`conda env create -f envs/sex_mapping.yml`, then `source tools/env.sh`). Tests: `python tools/run_tests.py 06`.

## sexloci, vendored

`sexloci/` is the package, `sexloci_dist/` the rest of sexloci 0.1.26 (README, METHODS and VARIANTS documents, CHANGELOG, configuration, examples, detection scripts and its own 105 pytest tests). **Both are copies, unchanged**:
`SEXLOCI_MANIFEST.tsv` lists the size and md5 of every file, and `tests/test_vendored_sexloci.py` fails when one changes or when a file is added. Change sexloci in its own project; to take a new release run
`python 06_sd_region_discovery/vendor_sexloci.py <sexloci-x_y_z.zip>` (it replaces both folders and rewrites the manifest), then the tests, and read `sexloci_dist/CHANGELOG.md`.

Read first: `sexloci_dist/README.md` (the rules, the outputs, the confidence tiers) and `sexloci_dist/docs/METHODS.md`. One VCF is one species is one analysis; module 05 builds the VCFs and checks a cohort against the same rules first.

Its genotyping commands (`gsdf-genotype`, `sv-genotype`, `gsdf-summary`, `write-builtin-vcf`, and the guarded sheet writer) come along unchanged but belong to module 04, which holds their parity-tested equivalents: use
`python -m sd_genotypes` for them (SPEC open decision 4 is who writes the genotype columns back).

## region_block: is a block of sex-linked signal real?

```
python -m sd_regions region_block --vcf <vcf> --qc <sexloci output>/sample_qc.tsv --region NC_135176.1:10000000-30000000 --core 13500000-25000000 --outdir <dir>
```

For a block that sexloci calls (or any region you want to examine): can it be made by low coverage, by inbreeding, or by family structure that lines up with sex? The VCF is read **once, through sexloci's own reader**
(same QUAL and MQ filters, the same `--min-dp`, the same call-rate rule; one difference: the filter of INFO/DP above its 99th percentile is not applied, which would need a second pass and removes about 1% of the sites).
Sample and sex come from sexloci's `sample_qc.tsv` (`--qc`) or any table with the columns `sample` and `sex` (`--sex-file`). `summary.md` states the numbers in words and how to read them; the files are listed in
`sd_regions/region_block/cli.py`.

| question | what is computed | files |
|---|---|---|
| Are the females (or males) homozygous over the block, and is that inbreeding? | heterozygosity per sample per 250 kb window; runs of homozygosity (4+ adjacent windows with at most 3% heterozygous calls: a threshold, not a model); F_ROH inside and **outside** the region; F from heterozygosity | `het_heatmap.png`, `samples.tsv`, `roh_segments.tsv`, `window_het.tsv` |
| Are the fish related, and does relatedness line up with sex? | genomic relationship matrix from the genome outside the region's chromosome; same-sex minus mixed-sex mean, permutation p; groups of related fish and their sexes | `kinship.png`, `kinship.tsv`, `kinship_pairs.tsv` |
| Does the block split the fish in a way the rest of the genome does not? | PCA of the region's SNPs against the genome-wide PCA from the relationship matrix; AUC of each for sex | `pca.png`, `pca_region.tsv`, `pca_genome.tsv` |
| Is there a block of linkage disequilibrium, in which sex? | r^2 of up to 300 thinned SNPs in the region and in the rest of the chromosome, for all fish, males, females; r^2 by distance | `ld.png`, `ld_decay.tsv` |
| How robust is the male-diagnostic count? | sexloci's rule (up to 10% exceptions) and a strict rule (none), at depth 3, 5 and 8; random subsets of fewer fish (default 9+7 and 8+8); the heterozygous share of every fish at the diagnostic sites | `diagnostic.png`, `diagnostic.tsv`, `diagnostic_subsample.tsv` |
| Is the region the only place where the sexes differ? | male and female heterozygosity of every chromosome | `chromosome_sex_het.tsv` |

`--subregion START-END` (repeatable) repeats the PCA and the diagnostic counts for parts of the region (to find where a block starts or ends). Why the diagnostic count is not compared between runs: with 12 fish of a sex the rule allows one
exception per sex, with 9 or 7 none, so the same data give very different counts (the first saulosi run, 9 males and 7 females at depth 5, found 16 in LG10's block; the run of 12 and 12 at depth 3 found 605, while the window heterozygosity
profiles of the two runs correlate at r = 0.97).

Limits, so they are not forgotten: the relationship coefficients are biased by about -1/(n-1) off the diagonal (allele frequencies come from the same few fish: compare pairs, do not read pedigree values); the PCA p-values pick the best of
three components; runs of homozygosity are called by a threshold on window heterozygosity (resolution one window), not by a hidden Markov model, so cross-check a surprising result with `bcftools roh`; allele balance of the heterozygous calls
(a mapping-artifact check) is not included. It tests alternative explanations; it does not decide that there is an inversion (breakpoints, an independent pedigree, a male assembly). Run time on a whole-genome VCF is **not measured yet**
(one pass; the `allel` backend is recommended, as for sexloci).

## Tests

`tests/test_region_block.py` (46 tests): every statistic on data with a known answer (windows, runs of homozygosity, relatedness and its sex test, PCA, LD, the diagnostic rule and its dependence on the number of fish), the VCF reader (depth mask, call-rate rule, exact region boundaries, which sites feed the relatedness), and the command end to end on a planted VCF (a Y haplotype, homozygous females, one inbred female) with its outputs and summary. The mutation checks and what they found are in the increment notes below.
`tests/test_vendored_sexloci.py`: the manifest; the import and version; the command-line delegation; sexloci's simulator (an XY locus planted on chr2 is recovered at tier A by all 3 families; a wrong genome version is
refused by the preflight); and, **when pytest is installed**, sexloci's own 105 tests (skipped otherwise: they were not run when this module was built, because the sandbox had no pytest).
Not exercised, as in sexloci's own README: the scikit-allel backend, the live rclone and the live Google sheets.

## Current state

Copied in and tested as above. Not done yet (SPEC module 06, Work): parity on the positive controls (the chr10 Y inversion in the YH pedigree, gsdf dosage in *M. conophoros*, the LG5 locus in the orange-blotch
species), a configuration for bcftools VCFs, retiring the two older scripts, replacing sexloci's own sheet access (`LabData`) by the shared file manager, and a location for the results on Dropbox.

## Increment: region_block (2026-11-12)

Built after the saulosi run at `min_dp` 3 (24 fish, 12 and 12): its LG10 block (13.5 to 25.0 Mb, seven tier C windows, 605 male-diagnostic SNPs) raised three objections, the sex-linked haplotype being an artifact of low coverage, of inbreeding or of family
structure, and nobody had computed runs of homozygosity or relatedness. `python -m sd_regions region_block` (above) answers them from the same VCF. What the two window tables of the saulosi runs (24 fish at `min_dp` 3, 16 at `min_dp` 5) already showed,
and why the tool measures what it measures: the heterozygosity and Fst profiles of LG10 are the same in the 16 high-depth fish alone (window dHet r = 0.97), so the 8 shallow fish do not create the block; the diagnostic count is **not** comparable (16 against 605);
female heterozygosity on LG10 is very low (26 of the 35 windows of the genome below 0.12 are on LG10) while male heterozygosity is only moderately raised, which is what runs of homozygosity in females would look like, hence the per-sample runs.

**Tests: 46, all passing; module 06 in all 54.** Mutation checks (19 mutations of the code, each must fail a test): 3 survived at first and each exposed a missing test, now added: a run of homozygosity continuing across a window that holds no SNP, a depth exactly equal to
the threshold (min_dp 5 keeps depth 5), and the exact start and end of the region (the region is half-open). Looking for the last one also found that the relatedness matrix of the sites *inside* the region, which the reader computed, was never used: removed.
Another finding of the planted-data runs: the first LD figure showed a blank panel and no panel where the SNPs do not vary within a sex (males all heterozygous, females all alike), which is the Y-haplotype signature itself; the SNPs are now chosen by
variation within the group, and an empty group shows as a labelled panel and as `n_sites = 0` in `ld_decay.tsv`.

**Not exercised:** the `allel` backend (not installed in the sandbox; its code path returns the same `Chunk`, with the raw dosage and the depth, and the tool reads only those), a real VCF, the run time of a whole-genome VCF, the figures on real data.
**Not done:** allele balance of the heterozygous calls (a check for mapping artifacts); a hidden Markov model for runs of homozygosity (`bcftools roh` as the cross-check); comparing two runs' window tables (done by hand above); sub-block scans of
the whole chromosome (local PCA in sliding windows). Open question of the module unchanged: SPEC open decision 4 (who writes results back).
