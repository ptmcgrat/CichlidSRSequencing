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

It reads the output file, takes the newest successful row of every sample measured with the CURRENT definition (older rows are never written), and fills `EffectiveCoverage`, `MapqPassFraction`, `BaseSurvivalFraction` and `DuplicationRate` (rounded to 4 decimals)
in the worksheet of the genome, and, since 2026-11-11 (later the same day), **`DuplicationRate`**. The writer is `common/cichlid_common/sheet_writer.py` (also for the SampleDatabase later); its rules: only the named cells; the default is a plan; a cell that already holds a different value is
left alone (`conflict`) unless `--overwrite`; a sample in no row or in more than one row is skipped and named; missing columns are added at the END of the header (what `align_pipeline setup-database` does, and running pipelines do not mind: they lay a row out
by the header they read at that moment, so a column they do not know stays blank in their rows); it re-reads and re-plans immediately before writing, writes `RAW` values in batched calls (gspread_retry repeats them after a network error or HTTP 429),
reads the sheet again and checks that every cell sits in the row of the right sample, and keeps a journal (`<output folder>/sheet_writes/`, planned and verified lines with the old values). Sheets has no conditional write, so a row inserted or deleted by
a person in the few milliseconds between the final read and the write could put a value one row off: the check catches it and names the cells (docstring of the writer). For these new columns that cannot destroy anything; for the SampleDatabase it is why
the first writes there will be `--limit` trials.
**The old columns are not touched.** `UsableCoverage`, `NRuns`, `Instrument` and `ReadLength` stay in the sheet until batches 3 and 4 have finished and the new code is deployed (a running pipeline refuses to write a row when the header lacks a column it writes).
52 + 168 tests with a worksheet test double (`cichlid_common/fake_sheet.py`); the guards were mutation-tested (a verification that always passes, a plan that treats a different value as writable, a row off by one and writing rows of every definition each make tests fail).

**`DuplicationRate`** is the column the alignment pipeline already fills for the samples it aligns (`DuplicationRate` <- `DuplicateFraction`: Picard's PERCENT_DUPLICATION, or the flag-based value when there are no Picard metrics). Picard
examines MAPPED reads only (unmapped reads are in neither the numerator nor the denominator), and so does the pipeline's own flag-based value (`reads duplicated / mapped`). `effcov` therefore writes `primary duplicates / primary MAPPED reads` (the new
TSV columns `PrimaryMappedReads` and `DuplicationRate`, from the `primary mapped` line of `samtools flagstat`), NOT its older `DuplicateFraction`, which divides by all primary reads and is lower by the unmapped share (29.1% against about 29.8% for a library with
2.2% unmapped). The sheet column is filled in EMPTY cells only: the old-pipeline rows get a value, the pipeline's own rows keep Picard's, and where the two differ the plan prints `conflict` with both values so the difference can be judged
(`--overwrite` would replace them, and it applies to all four columns). **Rows measured before this column existed cannot be corrected afterwards** (the mapped count was not stored): they are not written for `DuplicationRate` and a plain `measure` of
their samples measures them again (the plan says so); the output file, whose header is older, is migrated automatically before the first new row is appended (a copy is kept as `<file>.bak-<time>`).
`DuplicationRate` is not `DuplicationReads`, the legacy column of the old sheet (reads that suggest a structural duplication).

`effcov report --by COLUMN [--summary-only]` groups the samples by any column of the SampleDatabase or the AlignmentDatabase (Species, RunIDs, BroodID, Ecogroup ...) and prints, per group, the number of samples and the medians of the duplicate rate (with its range),
the MAPQ pass, the bases kept and effective / table coverage, and how many samples are below the minimum depth: the view for judging which libraries are poor. A sample sequenced in several runs has a RunIDs like `run1,,run2` and forms its own group: its
duplicate rate is measured on the merged BAM and cannot be split between the runs.

### Selecting samples: `--lab-reared`

`--lab-reared` selects the samples of the AlignmentDatabase that the SampleDatabase marks LabReared (1, true or yes, as the cohort rules read it; samples marked but without a row in the AlignmentDatabase are counted and skipped). The live sheets reach the code through gspread, which turns the text `1` into the number 1, `TRUE` into a boolean and a blank into NaN: every cell is read through one helper that handles all of them (a first version read only text and found no
LabReared sample at all on 2026-11-11, and the plan said only "no samples are selected"); an empty selection now prints what the column holds. It combines with
`--exclude-ecogroup Mbuna`. In the snapshot of the sample table (613 LabReared samples, the live sheet has more): 267 backcross hybrids, 91 *Mchenga conophoros*, 65 YH, 52 *Copadichromis virginalis*, 24 saulosi, 24 MCxCV F1, 21 gertrudae and others; 22 are Mbuna.

## Not built yet

The run-of-homozygosity track, the genotype-level heterozygote filter and the per-sex heterozygosity annotation, the accessibility mask, kinship and F_ROH in `qc`, hardy-weinberg groups other than `per_sample` (refused for now), and the
brood benchmark (SPEC decisions 59 to 68). The old `callset build` (GATK) still works and still reads the same file.

## Decisions made while building (continued)

* 2026-11-10: the first builder stages whole BAMs (rclone) and calls from local disk; streaming by region from Dropbox is a later option. Contigs are the `NC_` contigs of at least 5 Mb (the mitochondrion and unplaced pieces are left out).
* 2026-11-10: a trial is local only by construction, and an existing output is never replaced without `--force`.
* 2026-11-10: `bcftools` lives in its own environment (`envs/variant_calling.yml`), pinned to the release of the samtools of `align`, so the alignment toolset is not touched.


### Increment (2026-11-13): builds on another genome keep their own names

`cohort build <name> --genome G` already used G for the BAMs, the reference and the ledger, but the output names (`Outputs/PedigreeVCFs/<name>.vcf.gz`) and the work folder (`Cohorts/<name>/`) did not carry the genome: a build on a second assembly would have stopped at the first genome's VCF, and with `--force` overwritten it, and its resumable `per_contig/` folder would have been taken for its own. Now anything built on a genome other than Mzebra_GT3_NCBI gets the genome in its names (`<name>.<genome>.vcf.gz`, `Cohorts/<name>.<genome>/`, also for trials), by `cichlid_common.common_light.genome_tag`; the names on Mzebra_GT3_NCBI are unchanged, so nothing that exists is renamed or lost.
`contigs: chromosomes` (the NC_ contigs of at least `min_contig_bp`) is now an error that says what to do when it finds none (a reference whose chromosomes are not called NC_), instead of an empty selection. Tests: 190 in this module (4 new: coexistence of the two genomes, a failed build on one is not resumed by a build on the other, trials, the error).

### Increment (2026-11-17): `panel`, module 04's panel over existing BAMs and its calls into the SampleDatabase

The alignment pipeline (03) makes the panel's measurements for the samples it aligns, in the same pass. The 573 mbuna (5.8 TB on Dropbox, aligned to Mzebra_GT3_NCBI before the panel existed) need the same measurements made later, one BAM at a time, and the panel's calls written into the SampleDatabase (the handoff of the module 04 chat, `PILOT_HANDOFF.md`, decisions 76 and 77). It is built here, next to `effcov` (the same plan / `--run` / resumable-log / injected-fakes pattern), and calls the other modules, copying none of their logic:
module 03's `align_pipeline/panel.py` (the commands, the preflight that stops a run whose definition keys are not the frozen ones, the rule that says whether a JSON belongs to a BAM), module 04's `sd_genotypes panel measure` as a subprocess (`--panel-python`: the interpreter with its dependencies, e.g. the `sdpanel` environment), `build.stage_bam` for the download, and the guarded sheet writer for the cells.

```
python -m variant_calling panel measure (--samples-file F | -s S ... | --ecogroup Mbuna) [--run] [--limit N] [--jobs N] [--redo] [--no-upload] [--keep-bams] --panel-python /path/to/sdpanel/bin/python
python -m variant_calling panel write   [selection] [--json-dir D] [--tab COPY_OF_THE_TAB] [--apply] [--limit N]
python -m variant_calling panel table   [selection] [--json-dir D] [--out F]
```
`config/panelpass.yaml` holds the frozen keys (`expect_keys`), `jobs` (2: parallel downloads are untested, the first batch prints each BAM's and the whole run's MB/s), the disk factor, the places (under the local master, so on utaka under /Data and never /home) and `sheet.write_columns`.

**measure**, per sample: a free-disk check, download of `<sample>.all.bam` and `.all.bai` from `Bamfiles/<genome>/<sample>/` (size checked against the ledger, `samtools quickcheck`), `panel measure --modules all` into a `.part` file, the check that the JSON names THIS BAM (name and size), the JSON copied to Dropbox next to the BAM and its size read back there (the durable record: Temp is not a record), the BAM deleted. One row per sample goes to `Outputs/PanelJSON/<genome>/panel_pass.<genome>.tsv` as it finishes. A sample is measured again when its last row is not `ok`, the definition keys changed, the BAM's size changed (a realigned sample keeps its name and not its size), its JSON is gone, or with `--redo`.
A failed sample is recorded with its reason and never stops the others; running the same command again retries only those. A non-zero exit of the panel is a failure even if it left a JSON. `--limit N` takes the first N samples still to do, so the full set runs in resumable batches. Without `--run` it prints the plan (what is done, what is to do, the keys, the GB, the places) and measures nothing.

**write**, from the JSONs only (never a BAM): the calls of `sd_genotypes.panel.calls`, filtered to `gsdf number`, `chr7ins`, `TE1`, `TE2` and their flags columns; **Inversion5 and B score are never written** (module 04's hold until about 50 mbuna have been looked at), they stay in the JSONs. Rules of the handoff, all in the code and the tests: empty cells only; a call with a blocking flag writes its flags cell and not its value; the flags cell starts with the call tag (`pilot-v1`); a different value already in a cell is reported and never overwritten, and there is
NO overwrite option; when a value cell conflicts, its flags cell is left alone too (a panel tag next to a hand-entered, different value would claim a call that is not in the sheet); a measurement made under another definition (`DefinitionMismatch`) is refused and listed. New columns are added at the end of the header by the guarded writer. A plan unless `--apply`; `--tab` names a copy of the SampleDatabase tab for trials; a journal of every change is kept in `<json dir>/sheet_writes/`.
**table** writes every call of every JSON (all six columns, the flags as `sheet_row` makes them, the numbers of each call's `details`: LG5 median site fractions, B scores, read counts) into one TSV and prints, per column, the values of the confident calls and the commonest flags: what module 04 asks for after about 50 mbuna.

**Tests**: 56 in `tests/test_panelpass.py`, with injected fakes for Dropbox (a folder), rclone, samtools and the panel subprocess; and 5 in `tests/test_panelpass_integration.py` against module 04's REAL call layer, with panel JSONs made under the real (frozen) definitions by module 04's own fake-BAM helpers: they check what the other tests' stand-in only describes (the interface, the calls of 2, 3 and 4 copies, the dilution flag, no-calls, a non-confident gsdf, the refusal of a measurement of another definition, the table). They are SKIPPED when module 04's call layer is not in the checkout. No test reads a real sheet.
30 mutation checks on the rules above (identity, upload check, BAM deletion, free disk, held columns, flags cells, conflicts, resume rules, `--limit`, plan-only defaults, the column order, the table's tolerance of an unreadable sheet, the import of the call layer), all caught in the end; the ones that survived at first each exposed a missing test (a panel that writes its JSON and then exits non-zero; the import by full name; the table without a readable sheet).
Bugs the tests found while building: the log line divided by a download time of zero (a BAM kept from an earlier run), which would have stopped the whole pass after the first such sample; `table` stopped when the SampleDatabase could not be read; the order of new columns depended on which sample came first.
**The shared writer** (`common/cichlid_common/sheet_writer.py`) gained an OPTIONAL `column_order` argument for the order of NEW columns (the panel gives value column, its flags column, next value column ...); without it nothing changes (the order of first appearance), so `effcov write` is untouched.
**Module 04 must be committed first.** The git HEAD of 2026-11-17 has no `sd_genotypes/panel/calls.py`, none of the `status`/`version`/`sheet` keys of the four panel configs, no `docs/PILOT_HANDOFF.md` and no `tests/test_panel_calls.py`: they were delivered by the module 04 chat and exist only in the copy on utaka (`m04_utaka.zip`: a superset of HEAD, 4 configs changed by additive keys, one doc changed, three files new, nothing of HEAD lost; its 122 tests pass). `panel write` and `panel table` need it.
**Not tested:** the real rclone, samtools and panel subprocess on a real BAM (`--panel-python` through `env PYTHONPATH=...`), the real Google sheets (the live tab path was only tested with a stand-in file manager), and whether parallel downloads speed anything up.

### Increment (2026-11-18): `panel write` puts `./.` in a genotype cell that could not be called, and the flags columns are optional

Both are settings of `config/panelpass.yaml` (`sheet:`); module 04's call layer is untouched.

**`nocall`** (on in the shipped config): a genotype cell (`chr7ins`, `TE1`, `TE2`) whose call has the flag `no_call` or `low_gq` (too little data to call, or a call below the GQ cut) and no `site_error` gets `./.` instead of staying empty. A `site_error` is an error of the site, not a lack of data: that cell stays empty and the flag is in `panel table`. The numbers (`gsdf number`, `Inversion5`, `B score`) are never `./.`: with no confident call they stay empty. An information-only flag (for example `dilution_possible`) does not change anything: the value is written. The writer still never overwrites: a hand-entered value next to a `./.` is a reported conflict. Remove the `nocall` block to go back to empty cells. The 36 withheld cells of the first mbuna write (gsdf 10, chr7ins 8, TE1 7, TE2 11) were empty, so a second `panel write` fills in the genotype ones that are `no_call` or `low_gq`.

**`flags_in_sheet`** (`true` in the shipped config, which is the old behaviour): `false` writes only the values, with no `<column> flags` cells. The flags and the numbers behind every call stay in the JSONs and in `panel table`; the call tag (`pilot-v1`) is in the file name of each journal under `<json dir>/sheet_writes/` (`SampleDatabase.<time>.pilot-v1.tsv`), which also holds the old and new value of every cell written. The columns already in the sheet are not removed. What the sheet then no longer says: why a number is empty (the table does), a `dilution_possible` next to a `chr7ins` `0/0`, and which settings made a cell (the journal does).

**Tests**: 11 more in `tests/test_panelpass.py` (the rule on hand-made calls: no_call, low_gq, `no_call` together with `site_error`, a number with a nocall-type flag, a confident value that carries `low_gq`, a withheld call with an unrelated flag; the config checks; the flags option; the journal name; `panel table` keeps every flag) and 3 in `tests/test_panelpass_integration.py` against the REAL call layer (M1, M3, M4 of its scenarios). 11 deliberate breakages of the new rules, all caught by the tests. 262 + 3 tests in the module.
**Update (2026-11-18, later)**: module 04 lifted the hold (decision 78) and the shipped config now writes `Inversion5` and `B score` (with their flags columns). See the next increment.

### Increment (2026-11-18, later): `Inversion5` and `B score` written; `scope`; `./.` for Inversion5; a read-back script

- `sheet.write_columns` is now `[gsdf number, chr7ins, TE1, TE2, Inversion5, B score]`.
- **`scope`** (new): `scope: {Inversion5: [Mbuna], B score: [Mbuna]}` writes a column (value AND flags cell) only for samples whose `Ecogroup` is listed, because module 04 validated LG5 and the B blocks on the mbuna only (decision 78: do not apply them elsewhere). The Ecogroup of each sample is read from the same tab the cells are written to (the tab must have an `Ecogroup` column, otherwise nothing is written and the command stops with exit 3). A sample whose Ecogroup is empty or unknown is outside every scope. The plan says how many samples were left out of each column. gsdf number, chr7ins, TE1 and TE2 are written for every ecogroup. Remove the block (or add an ecogroup to a list) to change it.
- **`nocall`** now also covers `Inversion5` (an inversion genotype: 0, 1 or 2 inverted chromosomes) with the flags `low_depth` and `few_markers` added: a fish that cannot be genotyped gets `./.`. `exclude_flags` gains `unvalidated_value`: an Inversion5 of 2 is held back on purpose (decision 79) and stays empty, with its reason in the flags cell. `B score` and `gsdf number` are numbers and stay empty. A sample that was never measured (no BAM, no JSON) is not touched: `./.` means measured and not callable, not "not looked at".
- **`scripts/panel_readback.py`**: recomputes every call from the stored JSONs and compares it with an export of the live tab; the expected cell is derived from the call and the config (value, `./.`, empty, or out of scope), independently of the writer. `PANEL=<json folder> SHEET=<export> [CONFIG=...] python 05_variant_calling/scripts/panel_readback.py` after `source tools/env.sh`; it ends with PASS or FAIL. Module 04's own read-back script treats a `./.` as an unexpected value and knows nothing of the scope: use this one once `nocall` and `scope` are on.
- Tests: the writer tests run on an explicit four-column config without nocall or scope (they are about the writer); 7 new tests for the scope and the extended nocall on hand-made calls, 4 against the REAL call layer (the six columns, a non-mbuna fish, module 04's guard simulated, a tab without an Ecogroup column) and 1 for the read-back script; 25 deliberate breakages of the rules, all caught.
