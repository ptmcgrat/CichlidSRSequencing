# Methods

One VCF = one species = one analysis unit. Genotypes are reduced to alt-allele dosage (0/1/2); genotypes
below `min_dp`/`min_gq` become missing; markers with call rate < `min_call_rate` are dropped.

## 0. Preflight (before any genotype is read)
Blocking: VCF sample missing from the sample table; any sample whose Sex is not exactly M/F; more than one
Species value (or a blank one); alignment records made against a different genome than `--genome-version`;
fewer than 3 samples of either sex. Warnings and notes: see the README table.
After the genotype pass the QC accumulator adds depth, missingness, sex-depth, heterozygosity-depth,
table-vs-VCF coverage and Mendelian-error warnings.

## 0b. Variant filtering
Applied while reading, before any statistic: FILTER labels other than PASS/'.', hard INFO thresholds (QUAL<30, QD<2,
FS>40, MQ<50, SOR>3, MQRankSum<-12.5, ReadPosRankSum<-8; skipped if the field is absent) and INFO/DP above the 99th
percentile (found with a fast sites-only scan). AF, InbreedingCoeff/ExcessHet and a lower depth cutoff are *not* used
because they remove sex-linked sites (all-male heterozygosity; hemizygous low-depth X/Y or Z/W sequence).
`filter-compare` quantifies the loss in an existing filtered VCF.

## 0c. Genotype thresholds by role
The transmission test and the sex contrasts do not need the same genotype quality. In a pedigree the question is which of
a parent's two haplotypes each offspring received; with rare recombination that is constant over long blocks, so each
window pools hundreds of informative markers (parent heterozygous, other parent homozygous) into one consensus. Low depth
produces one-sided errors (a heterozygous offspring looks homozygous when every read shows the same allele, p = 0.5^(depth-1))
that the consensus absorbs. The transmission test therefore uses offspring genotypes at depth >= 2 with no GQ filter and
parents at depth >= 6 / GQ >= 20, and is not subject to the call-rate filter. The sex contrasts compare called genotypes
(heterozygosity in particular is under-called at low depth), so they keep depth >= 5; samples that keep < 50% of their
genotypes (measured on the first 20,000 sites) are excluded from them, which prevents a few low-depth animals from removing
sites for everyone. Parentage is checked only at sites where both parents are homozygous for the same allele, because at
sites with opposite homozygous parents a low-depth heterozygous offspring is miscalled and would look like a Mendelian error.
In simulation (3 families x 16 offspring, depth-dependent heterozygote under-calling) the locus is recovered as tier A with all
three families supporting it with offspring at 3x and 2x (XY and ZW); the earlier single-threshold pipeline failed at 5x and 3x.

## 0d. Inheritance tracks (whole-genome haplotype view)
For every family, parent and window the transmission test already computes a consensus score per offspring (first principal
component of the informative-marker matrix). The sign of that score is arbitrary in each window, so for display each window is
flipped to agree with the previous one and the state is the sign of the aligned score. The result is the inherited parental
haplotype of every offspring along every chromosome; orientation is arbitrary per chromosome but constant along it. This is
display only: the statistics do not use the alignment.

## 1. Strata = units of replication
Offspring of the same sire x dam form a stratum; all other (non-parent) animals form one `pool` stratum.
Parents are used by the transmission test only. Families that share a parent are merged into one
*independent family group* (union-find over parents) because they are not independent replicates.

## 2. Association with sex (`association.py`)
Per stratum and SNP: `chi2_s = n_s * max(r_dosage^2, r_het^2)`. Squared correlations because the allele
linked to sex differs between families. The combined statistic is the sum over strata. Sex labels are
permuted within strata. The genome-wide maximum per permutation gives FWER-adjusted p for windows and SNPs;
per-stratum windows are corrected by max-Z (below) to decide whether each stratum supports a region.

## 3. Pedigree transmission (`pedigree.py`)
Informative markers = parent heterozygous, other parent homozygous; the transmitted allele is read directly
from the offspring. Per window, the first principal component of the marker x offspring matrix (sign
invariant) is the consensus "which parental haplotype was inherited", tested against offspring sex:
`T = sum_families n * r^2`, within-family permutation, Westfall-Young max-T over all windows and both roles.
Sire role significant only -> XY-like; dam only -> ZW-like. Each family is also tested on its own.

## 4. Population-genetic contrasts (`popgen.py`)
Windowed Weir-Cockerham Fst, male-minus-female heterozygosity, sex-diagnostic SNP counts and the male:female
depth ratio from variant-site DP (a coarse indicator). Outliers use robust z-scores against the genome.

## 5. Calling and grading (`regions.py`)
Windows are candidates if the combined association or transmission test is significant on the *plateau*
(>= 90% of the chromosome's best statistic) or >= 2 population-genetic outliers coincide. The two primary
tests are Bonferroni-corrected (alpha/2) so the dataset-level FWER stays near alpha.
Independent units supporting a region = family groups and the pool that each reach FWER significance there.
* **A** combined test significant and >= 2 independent units significant (replication)
* **B** combined test significant, replication not shown
* **C** population-genetic outlier only
Association and Fst/heterozygosity are the same information and are never counted as confirmation of each
other. Data-quality findings can cap the tier (sex confounded with batch, sex-depth bias, Mendelian errors,
unusable coverage).

## Validation on simulated data (`tests/`)
XY and ZW loci are recovered from three replicated families (tier A, correct system), from unrelated animals
(tier B, by design a single pool cannot show replication) and from a single family (tier B). In 70 null datasets
(no sex locus) 5 produced a tier B call (7%, mostly transmission-only), none a tier A; with offspring at 3x, 3 of 30 null
datasets produced a tier A/B call. Both are consistent with the 5% target but too few runs to show it is exactly 5%. The simulations are idealised (no population
structure, no linkage disequilibrium beyond the block, simple error model); thresholds should be revisited on
real cohorts such as the *Aulonocara* chromosome 10 pedigree.

## Limits
Y/W sequence missing from the reference is invisible to SNP calling. Permutation p-values resolve to
1/(n_perm+1). Tiny cohorts are underpowered. Polygenic or environmental sex determination appears as weak
regions or none. The `allel` VCF backend, the live rclone/Dropbox download and the live Google-sheets read were not exercised
in the development environment (the lab-data layer was tested over a fake rclone and fake Google modules).
