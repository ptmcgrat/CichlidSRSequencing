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
(`per_sample` by default: no shared Hardy-Weinberg prior for inbred lab and farm stock, decision 63). Selectors are in the docstring of `variant_calling/cohort.py`: explicit `samples`, `broods` (with
their parents), `species` / `ecogroups` / `subgroups` / `categories`, `exclude_samples`. Fish whose Ecogroup is Error or Exclude are skipped from criteria-selected cohorts, as in the alignment pipeline.

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

`family` and `sex_set` follow sexloci 0.1.26's preflight, which module 06 applies again to the VCF. **Where 05 is stricter on purpose:** a member with no alignment BLOCKS here (sexloci only notes it, because
05 calls from the BAM), and a `family` cohort BLOCKS when a listed parent is missing or when nobody lists parents (sexloci warns, because it can still use the offspring). `bc_parents` exists because the BC
cross is what sexloci refuses (more than one species in a VCF); its VCF is for module 07.

Checked on the project's older sample sheet export with a stand-in 20x ledger (2026-11-10): `yh_broods` (52 samples), `mc_broods` (51) and `bc_family` (7) pass; `yh_broods_plus_adults` is BLOCKED by one adult,
`YH_002_m_`, whose Sex is `U` although its name says male: a sheet fix.

## Not built yet

The bcftools builder (`cohort build`), the run-of-homozygosity track, the genotype-level heterozygote filter, the accessibility mask, kinship and F_ROH in `qc`, and the brood benchmark (SPEC decisions 59 to 68).
The old `callset build` (GATK) still works and still reads the same file.
