# 05 VCFs for cohorts of 4 to 100 fish

Specification: [`SPEC.md`](../SPEC.md), module 05 (decisions 58 to 68). Package: `variant_calling` (the number is only for ordering folders).

```
python -m variant_calling cohort list
python -m variant_calling cohort members <name>
python -m variant_calling cohort check <name> [--json report.json]
python -m variant_calling cohort check <name> --sample-table SampleDatabase.tsv --alignment-table AlignmentDatabase.tsv     # offline, from exports
python -m variant_calling cohort build <name>                          # the PLAN: BAMs and GB, contigs, the exact bcftools commands, outputs, preflight problems; writes nothing
python -m variant_calling cohort build <name> --run                    # stage the BAMs, call, merge, upload and verify
python -m variant_calling callset build <name> [--dry-run]            # the GATK route of the BC work
```

Tests: `python tools/run_tests.py 05` (or `python -m variant_calling selftest`; 110 tests, about 20 s). Configuration: `config/callsets.yaml` (the cohorts), `config/cohort_rules.yaml` (check thresholds),
`config/cohort_build.yaml` (calling parameters). Environment: `../envs/variant_calling.yml` for `cohort build` (bcftools 1.24 and samtools 1.24); `list`, `members` and `check` need only the sheets' libraries (`align` works).

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

## `cohort build` (increment 2)

`cohort build <name>` runs the cohort check first (BLOCK stops it), then prints the **plan** and runs the **preflight**: the tools are on the PATH, this bcftools has `call --group-samples` and `mpileup --ignore-RG`,
the reference and its `.fai` exist, every member has a `BamSize` in the ledger, and the disk holds the BAMs times 1.2. `--run` then does, resumably:

1. **stage** the members' `<s>.all.bam` and `.all.bai` from `Bamfiles/<genome>/<s>/` into `<master>/Cohorts/<name>/staged/` (4 at a time; the size is compared with the ledger and `samtools quickcheck` is run);
2. **call** each chromosome-scale contig (the `NC_` contigs of at least 5 Mb), 8 at a time: `bcftools mpileup -q 13 -Q 20 -p -m 2 -F 0.05 -a FORMAT/AD,FORMAT/DP --ignore-RG | bcftools call -m -v -P 0.001 -G - -f GQ`
   (Blumer et al. 2025; no shared Hardy-Weinberg prior, decision 63). A failed contig fails the build with its stderr; the finished contigs are kept and a rerun calls only the rest;
3. **merge**: `concat`, `reheader` (the columns, one per BAM in cohort order, are renamed to the SampleIDs after the count is checked), `norm`, `+fill-tags` (AN, AC, NS, F_MISSING, ExcHet: annotated, nothing is filtered), the biallelic-SNP subset;
4. **manifest** (`<name>.manifest.json`): the members and their BAM sizes and coverage, the check report, every parameter, the bcftools and samtools versions, the commit, per-contig counts and times;
5. **upload and verify** the VCF, the SNP VCF, their indexes and the manifest to `Outputs/PedigreeVCFs/` (the folder `sexloci run-folder` reads), then delete the staged BAMs. An existing output, local or on Dropbox, is not replaced without `--force`.

Exit codes of `build`: 0 plan or built, 1 BLOCKED, 3 the plan has preflight problems, 4 build error. `--workers`, `--keep-staged`, `--keep-work`, `--no-upload`, `--contigs`.

**A trial first.** `--max-samples N` takes the N smallest BAMs and `--region CONTIG:START-END` one region: a cheap test of the plumbing, written under `<master>/Cohorts/<name>_trial/trial_out/` and **never uploaded**, so it cannot reach
the folder that module 06 analyses. Suggested first runs, on yellowhead, in `variant_calling`: the plan of a cohort (the contig names are in it); then
`python -m variant_calling cohort build saulosi_green_face_peacock --max-samples 4 --region <a contig>:1-3000000 --run`; then look at the trial VCF (`bcftools query -l`: the SampleIDs; `bcftools stats`; `bcftools view -h | tail`);
then a real build of a small cohort (`bc_family`, 7 samples), then the rest.

**What the tests establish and what they do not.** The tests run the builder against FAKE bcftools, samtools and rclone (`tests/fake_tools.py`): they show that it issues the intended command lines (asserted exactly) and that staging, resume
after a failed contig, renaming, the manifest, the trial's isolation and the refusals work. They say nothing about how the real tools behave, because none was installed where this was built. **Unverified until the first real run:** that
`-G -` gives each sample its own group (`bcftools call` usage and manual say), that `--ignore-RG` gives one column per BAM in list order (the builder stops if the column count is wrong), the `-a FORMAT/AD,FORMAT/DP` and `-f GQ` tag names of the
pinned version, that htslib finds the `.all.bai` index (a `.all.bam.bai` link is made as well), and the time and the size of a build (not measured).

## Not built yet

The run-of-homozygosity track, the genotype-level heterozygote filter and the per-sex heterozygosity annotation, the accessibility mask, kinship and F_ROH in `qc`, hardy-weinberg groups other than `per_sample` (refused for now), and the
brood benchmark (SPEC decisions 59 to 68). The old `callset build` (GATK) still works and still reads the same file.

## Decisions made while building (continued)

* 2026-11-10: the first builder stages whole BAMs (rclone) and calls from local disk; streaming by region from Dropbox is a later option. Contigs are the `NC_` contigs of at least 5 Mb (the mitochondrion and unplaced pieces are left out).
* 2026-11-10: a trial is local only by construction, and an existing output is never replaced without `--force`.
* 2026-11-10: `bcftools` lives in its own environment (`envs/variant_calling.yml`), pinned to the release of the samtools of `align`, so the alignment toolset is not touched.
