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

## Effective coverage (`effcov`, 2026-11-11)

```
python -m variant_calling effcov measure (--cohort NAME | --samples-file F | -s ID ... | --all [--exclude-ecogroup Mbuna])      # the plan
python -m variant_calling effcov measure ... --run                                                                          # download, measure, delete, append
python -m variant_calling effcov report [--tsv F] [selection]
```

The AlignmentDatabase's `Coverage` is raw yield. What the caller sees is much less and differs between libraries, so `effcov` measures it per sample from the sample's `all.bam` (downloaded, measured, deleted) and appends one row per
sample to `<local master>/Outputs/EffectiveCoverage/effective_coverage.<genome>.tsv`. **It does not write to any Google sheet** (reads are read-only); writing the numbers to the AlignmentDatabase is a later step that reads this file.

The definition (`effcov.py` docstring; parameters in `config/effcov.yaml`; the caller's own parameters come from `config/cohort_build.yaml`, so the measurement cannot drift from the calling):
`EffectiveCoverage` is the mean FORMAT/DP of `bcftools mpileup` run exactly as `cohort build` runs it (one function builds both commands) over one 1 Mb window on each of the 22 chromosomes (about 2% of the genome), at the positions
the caller emits, so it is the depth at covered sites, which is what a genotype is made from. Also: `WindowRawDepth` (INFO/DP: reads passing the read-level filters), `BaseSurvivalFraction` (the base-level filters: quality, BAQ, mate overlap),
`WindowCoveredFraction`, and over the whole BAM `MapqPassFraction` and `DuplicateFraction`. Every row carries a `DefinitionKey`; a sample measured under another key is measured again. A rerun skips samples already measured, a failed sample is
recorded with its reason and does not stop the others, `--redo` measures again, and the default of `measure` is the plan.

Why this definition (measured 2026-11-11 on three BAMs with the same pileup over six windows; "caller" is the depth in the windows, "table" the AlignmentDatabase's Coverage):

| sample | table | duplicates | MAPQ >= 13 | overlap (insert sizes) | base survival | caller | caller / table |
|---|---|---|---|---|---|---|---|
| ASG_212_1 (saulosi) | 9.9 | 29.1% | 72.6% | 34% | 68% | 4.1 | 41% |
| YH_005_m (YH parent) | 19.9 | 2.7% | 80.9% | 0.8% | 95% | 16.9 | 85% |
| SAMEA117805283 (Blumer) | 14.6 | 8.2% | 81.6% | 1.2% | 92% | 11.9 | 81% |

The saulosi library keeps 41% of its nominal coverage (short inserts, heavy PCR), the other two 81 to 85%, so the usable depth cannot be derived from the table. A formula from counts (non-duplicate coverage x MAPQ pass x (1 - overlap)) under-predicted the
caller by 12 to 34% (the windows are covered sites, the genome-wide mean is diluted by positions that are uncovered after filtering), and overlap is only part of the base loss (the two other libraries lose 5 to 8% of bases with about 1% overlap), so
the pileup itself is the measurement. **(unverified)** that 22 windows of 1 Mb represent the genome: the 24 saulosi samples are the check (effective coverage against sexloci's per-sample mean depth of 3.8 to 7.9).

Tests (37, with fake tools that describe a BAM by a few numbers, so every expected value is exact): arithmetic and parsing, the definition key, the output file, one sample (including that the pileup command equals the build's), a whole run with a
failure and a resume, the preflight, and the command line. They found a real defect while being written: `15.0e9` in a YAML file is a string, not a number (the config loader now coerces and rejects).
**The first real run is the check of the real tools** (`bcftools mpileup` over 22 comma-separated regions with the build's options, `samtools flagstat` and `view -c`): try one sample first, ASG_212_1, which has known numbers.

Agreed next: drop `UsableCoverage`, `NRuns`, `Instrument` and `ReadLength` from the sheet **only after batches 3 and 4 have finished and the new code is deployed**
(a running pipeline refuses to write a row when the header lacks a column it writes); the same pass over each BAM then runs the SD panel (gsdf, the three TE sites, LG5, B chromosome) and fills the SampleDatabase, fill-empty-cells-only.

### Writing to the AlignmentDatabase (`effcov write`, 2026-11-11)

```
python -m variant_calling effcov write                      # the PLAN: what would happen to every cell; writes nothing
python -m variant_calling effcov write --apply --limit 3    # a first trial on the live sheet
python -m variant_calling effcov write --apply              # the rest
```

It reads the output file, takes the newest successful row of every sample measured with the CURRENT definition (older rows are never written), and fills `EffectiveCoverage`, `MapqPassFraction` and `BaseSurvivalFraction` (rounded to 4 decimals)
in the worksheet of the genome. The writer is `common/cichlid_common/sheet_writer.py` (also for the SampleDatabase later); its rules: only the named cells; the default is a plan; a cell that already holds a different value is
left alone (`conflict`) unless `--overwrite`; a sample in no row or in more than one row is skipped and named; missing columns are added at the END of the header (what `align_pipeline setup-database` does, and running pipelines do not mind: they lay a row out
by the header they read at that moment, so a column they do not know stays blank in their rows); it re-reads and re-plans immediately before writing, writes `RAW` values in batched calls (gspread_retry repeats them after a network error or HTTP 429),
reads the sheet again and checks that every cell sits in the row of the right sample, and keeps a journal (`<output folder>/sheet_writes/`, planned and verified lines with the old values). Sheets has no conditional write, so a row inserted or deleted by
a person in the few milliseconds between the final read and the write could put a value one row off: the check catches it and names the cells (docstring of the writer). For these new columns that cannot destroy anything; for the SampleDatabase it is why
the first writes there will be `--limit` trials.
**The old columns are not touched.** `UsableCoverage`, `NRuns`, `Instrument` and `ReadLength` stay in the sheet until batches 3 and 4 have finished and the new code is deployed (a running pipeline refuses to write a row when the header lacks a column it writes).
52 + 168 tests with a worksheet test double (`cichlid_common/fake_sheet.py`); the guards were mutation-tested (a verification that always passes, a plan that treats a different value as writable, a row off by one and writing rows of every definition each make tests fail).

### Selecting samples: `--lab-reared`

`--lab-reared` selects the samples of the AlignmentDatabase that the SampleDatabase marks LabReared (1, true or yes, as the cohort rules read it; samples marked but without a row in the AlignmentDatabase are counted and skipped). It combines with
`--exclude-ecogroup Mbuna`. In the snapshot of the sample table (613 LabReared samples, the live sheet has more): 267 backcross hybrids, 91 *Mchenga conophoros*, 65 YH, 52 *Copadichromis virginalis*, 24 saulosi, 24 MCxCV F1, 21 gertrudae and others; 22 are Mbuna.

## Not built yet

The run-of-homozygosity track, the genotype-level heterozygote filter and the per-sex heterozygosity annotation, the accessibility mask, kinship and F_ROH in `qc`, hardy-weinberg groups other than `per_sample` (refused for now), and the
brood benchmark (SPEC decisions 59 to 68). The old `callset build` (GATK) still works and still reads the same file.

## Decisions made while building (continued)

* 2026-11-10: the first builder stages whole BAMs (rclone) and calls from local disk; streaming by region from Dropbox is a later option. Contigs are the `NC_` contigs of at least 5 Mb (the mitochondrion and unplaced pieces are left out).
* 2026-11-10: a trial is local only by construction, and an existing output is never replaced without `--force`.
* 2026-11-10: `bcftools` lives in its own environment (`envs/variant_calling.yml`), pinned to the release of the samtools of `align`, so the alignment toolset is not touched.
