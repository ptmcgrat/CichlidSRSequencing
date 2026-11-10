# 05 VCFs for cohorts of 4 to 100 fish

Specification: [`SPEC.md`](../SPEC.md), module 05 (decisions 58 to 68). Package: `variant_calling` (the number is only for ordering folders).

```
python -m variant_calling cohort list
python -m variant_calling cohort members <name>
python -m variant_calling cohort check <name> [--json report.json]
python -m variant_calling cohort check <name> --sample-table SampleDatabase.tsv --alignment-table AlignmentDatabase.tsv     # offline, from exports
python -m variant_calling callset build <name> [--dry-run]                                                                   # the GATK route of the BC work
```

Tests: `python tools/run_tests.py 05` (or `python -m variant_calling selftest`). Configuration: `config/callsets.yaml` (the cohorts), `config/cohort_rules.yaml` (thresholds). Environment: `../envs/analysis.yml`.

## What is built (increment 1): the cohort file and the cohort check

A cohort is an entry of `config/callsets.yaml` with a `kind` (`family`, `sex_set` or `bc_parents`), a `genome`, a `caller` (`bcftools` by default, after Blumer et al. 2025) and `groups`
(`per_sample` by default: no shared Hardy-Weinberg prior for inbred lab and farm stock, decision 63). Selectors are in the docstring of `variant_calling/cohort.py`: explicit `samples`; `broods` (the children are the samples with that BroodID, and their parents are the two samples their `Parents` column names);
`species` (the species name of the sheet, case and spacing ignored, so a shorter name is another species) with optional `ecogroups` / `subgroups`; `exclude_samples`. **The sheet's Category is not used to define a cohort**
(a definition that has `categories` is BLOCKED with a message). Fish whose Ecogroup is Error or Exclude are skipped from species-selected cohorts, as in the alignment pipeline.

Cohorts shipped in `config/callsets.yaml`: `yh_broods`, `mc_broods` (family), `bc_family` (bc_parents), `gertrudae_orange_collar` (*Aulonocara gertrudae orange collar Tanzania*) and
`saulosi_green_face_peacock` (*Aulonocara saulosi green face peacock*), both sex_set.

`cohort check` reads only the SampleDatabase and the AlignmentDatabase (the ledger): no BAM, seconds, read-only. A call costs hours, so the check comes first. Exit code 0 ok, 1 BLOCKED, 2 usage.

| rule | family | sex_set | bc_parents |
|---|---|---|---|
| every member in the sheet (BLOCK), one row each (BLOCK) | yes | yes | yes |
| every member aligned to the genome: a row in the genome's ledger, `GenomeVersion` matching (BLOCK) | yes | yes | yes |
| Sex exactly M or F for every member, parents included (BLOCK) | yes | yes | no: unknown sex is a WARN |
| one species (BLOCK); the declared `species` matches (BLOCK) | yes | yes | no: a cross has two species and F1s |
| at least 3 analysed males and 3 females, parents not counted (BLOCK); fewer than 6 (WARN) | yes | yes | no |
| offspring list their parents, and the parents are in the cohort (BLOCK) | yes | no | the pedigree must be complete inside the cohort (BLOCK) |
| one parent unknown, a malformed `Parents` field, swapped parent sexes, a brood under 12 offspring (WARN) | yes | no | no |
| coverage: offspring < 2x, parents < 10x, other fish < 5x (< 3x unusable) (WARN) | yes | yes | below `bc_parents_min_coverage` BLOCK (provisional, an open decision), below 10x WARN |
| the sheet lists relatives (siblings, parent and offspring) in the pool (WARN) | no | yes | no |
| sex perfectly confounded with BroodID, Project or LabReared (WARN) | no | yes | no |
| a member marked Ecogroup Error or Exclude, included by name (WARN) | yes | yes | yes |
| a selector that matches nothing: a brood, a species, an ecogroup or a subgroup (WARN) | yes | yes | yes |
| a BroodID's children name different parent pairs, or some leave Parents blank (WARN) | yes | no | no |
| lab-reared fish with no BroodID or Parents: the sheet cannot rule out relatives (INFO) | no | yes | no |

`family` and `sex_set` follow sexloci 0.1.26's preflight, which module 06 applies again to the VCF. **Where 05 is stricter on purpose:** a member with no alignment BLOCKS here (sexloci only notes it, because
05 calls from the BAM), and a `family` cohort BLOCKS when a listed parent is missing or when nobody lists parents (sexloci warns, because it can still use the offspring). `bc_parents` exists because the BC
cross is what sexloci refuses (more than one species in a VCF); its VCF is for module 07.

Checked on the live sheets (yellowhead, 2026-11-10, commit d799aa84): `yh_broods` (52 samples: 2 broods of 24 children and their 4 parents), `mc_broods` (51: two broods sharing a sire, 3 parents) and `bc_family` (7:
5 founders and 2 F1s of three species) pass with no warning. The two species sets were checked on an older sheet export with a stand-in ledger: *A. gertrudae orange collar Tanzania* 21 fish (12 M, 9 F) and
*A. saulosi green face peacock* 24 fish (12 M, 12 F), all lab-reared, Deep_Benthic, with no BroodID or Parents, pass the sheet-side rules; whether they are aligned is for the live check to say (`no_alignment`).
On the older export, 3 of 11 broods (the BC broods) have children with a blank `Parents` field, which the new brood rule reports.

## Decisions made while building (SPEC.md is not edited for increments)

* 2026-11-10: **Category is not used to define cohorts.** The two `*_plus_adults` cohorts, which added adults by Category, were removed: they added nobody on the live sheet, and no label of the sheet identifies such adults.
* 2026-11-10: new cohorts are named by the species name of the sheet; the first two are *Aulonocara gertrudae orange collar Tanzania* and *Aulonocara saulosi green face peacock*.
* 2026-11-10: a BroodID is one family of two parents (the children are found by BroodID, the parents by the Parents column); a brood whose children disagree is a WARN, not a BLOCK.
* The old GATK builder (`callset.py`) still has its own resolver that can read `categories`; no cohort uses it.

## Not built yet

The bcftools builder (`cohort build`), the run-of-homozygosity track, the genotype-level heterozygote filter, the accessibility mask, kinship and F_ROH in `qc`, and the brood benchmark (SPEC decisions 59 to 68).
The old `callset build` (GATK) still works and still reads the same file.
