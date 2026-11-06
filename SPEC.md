# CichlidSRSequencing: specification

Status: draft v1 (2026-11-06). Owner: McGrath lab. This file lives at the repository root and is the source of truth for what each folder is for, what it consumes, what it produces, and when it is done.
Where a statement is not yet verified it is marked **(unverified)**; where it is a decision still to be made it is marked **(open)**.

---

## 1. Purpose

Take public and in-house short-read data for Lake Malawi cichlids from a read archive to a set of results about sex determination (SD): where the SD loci are, what they are made of (inversions, transposon insertions, a gene duplication, a B chromosome), how they behave in backcross families, and which variants are candidate causal changes. Everything that goes into the paper must be reproducible from versioned code and recorded inputs.

The work is ten modules. Each is one folder. They communicate only through **data contracts** (section 4), never by importing each other's internals, so any module can be re-run, replaced or tested alone.

## 2. Principles that apply to every module

1. **One folder, one job.** A module has a README, a command-line entry point, a config folder, tests, and an environment file. Numbered folders are for humans; the importable package inside has no number.
2. **Raw measurements are stored; calls are made later.** Anything computed from a BAM or VCF stores raw numbers first (JSON/TSV) and makes calls in a separate, offline step, so a changed threshold or calibration never needs the big files again.
3. **Immutable, recorded runs.** Every step writes a run folder (inputs, parameters, tool versions, git commit, outputs, checksums) to Dropbox with a `LATEST` pointer, as the current `Run` helper does. A figure or table must name the run it came from.
4. **Definitions are hashed.** A measurement carries a hash of what was measured (region, windows, thresholds that define it). A caller refuses to mix measurements made under different definitions.
5. **Fail loudly, never silently.** Non-zero exit codes for failures, status files saying which stage failed and why, no database row for work that did not finish, no `sys.exit()` with status 0 after an error.
6. **Safe on shared resources.** Google Sheets are written cell-by-cell for the changed rows, never by rewriting a whole worksheet; runs are claimed before they are processed; Dropbox uploads are verified; disk is reserved before work starts.
7. **Tests without the real tools.** Logic is tested with fakes (fake `samtools`, in-memory sheet, planted signal) so that the whole suite runs in seconds on a laptop; the real tools are exercised by small, documented pilots.
8. **Real data decides.** Assumptions are written down together with the measurement that supports them (for example: ENA's `read_count` means *reads* for some runs and *pairs* for others, so no code assumes one meaning).
9. **Nothing is kept that need not be.** Reads and large intermediates are deleted once the compact result exists; the result records what was deleted.

## 3. Repository layout

```
CichlidSRSequencing/
  SPEC.md                         this file
  README.md                       how to get started; machine inventory; where data lives
  common/                         shared library (package: cichlid_common)
  01_structural_landmarks/        centromeres and inversion breakpoints               (package: landmarks)
  02_download_reads/              ENA / SRA download, download database               (package: download_reads)
  03_align_pipeline/              align, mark duplicates, split, slice, QC, upload    (package: align_pipeline)
  04_genotype_sd_systems/         gsdf duplication, transposon insertions, inversions, B chromosome   (package: sd_genotypes)
  05_variant_calling/             VCFs for families and species                       (package: variant_calling)
  06_sd_region_discovery/         SD regions in many species                          (package: sd_regions)
  07_bc_analysis/                 backcross genotyping, recombination, TD, QTL        (package: bc_analysis)
  08_ase_variants/                variants for allele-specific expression             (package: ase_variants)
  09_candidate_qtn_webserver/     candidate causal variants, allele frequencies, web server   (package: candidate_qtn)
  10_paper_figures/               figures and data tables                             (package: paper_figures)
  envs/                           one conda environment file per module
  docs/                           design notes, decision log, method write-ups
  legacy/                         retired scripts (kept for reference, not run)
```

Inside each module folder:

```
NN_name/
  README.md        purpose, how to run, current state
  <package>/       importable code; cli.py is the only entry point (python -m <package> <command>)
  config/          yaml: every threshold, window, path and accession that defines the method
  tests/           unit tests with fakes; one small integration test per command
  env.yml          the conda environment this module needs
```

Command-line convention: `python -m <package> <command> [options]`; every command that writes results accepts `--out-dir` (plain files, no Dropbox, no Project) and, without it, writes an immutable run folder.

## 4. Dependencies and data contracts

```
 01 landmarks ──────────────┐
                            ▼
 02 download ──▶ 03 align ──▶ 04 genotype SD systems ──────────────┐
                    │  │                                           │
                    │  └──────────▶ 05 variant calling ──▶ 06 SD regions
                    │                    │   │                     │
                    │                    │   └──▶ 08 ASE variants  │
                    ▼                    ▼                         │
                 07 BC analysis ◀── (03 and 05)                    │
                    │                                              │
                    └────────▶ 09 candidate QTN + web server ◀─────┘  (uses 03, 04, 07)
 every module ───────────────▶ 10 paper figures and tables
```

| Producer | Artifact (the contract) | Consumers |
|---|---|---|
| 01 | `landmarks/<genome>/centromeres.tsv`, `inversions.tsv`, `inversion_breakpoint_windows.bed`, `known_misassemblies.tsv` (schemas in 01) | 03 (windows to slice), 04, 07, 09, 10 |
| 02 | uBAM + adapter metrics per run on Dropbox (`Reads/<project>/<run>.unmapped_marked_adapters.bam`); rows in `DNAReadsDatabase`; `DownloadDatabase` tabs | 03 |
| 03 | `Bamfiles/<genome>/<sample>/<sample>.all.bam` + index, slices, `qc.json`; one row per sample in `AlignmentDatabase` | 04, 05, 07 |
| 04 | `<sample>.panel.json` (raw measurements, with definition hash); call tables (`gsdf_calls.tsv`, `te_calls.tsv`, `inversion_calls.tsv`, `bchrom_calls.tsv`); calls written back to `SampleDatabase` columns | 05 (cohort definitions), 06, 07, 09, 10 |
| 05 | per-cohort VCF (bgzip + tabix) with a cohort manifest | 06, 07, 08, 09 |
| 06 | SD candidate regions per species (BED + statistics TSV + plots) | 09, 10 |
| 07 | BC genotype matrices, breakpoints, recombination maps, TD tables, QTL results | 09, 10 |
| 08 | per-gene, per-family ASE variant lists | 10 (and the RNA-seq analysis) |
| 09 | candidate-QTN tables, allele and genotype frequency tables, the web server's data store | 10 |

**Shared stores.** Dropbox `Apps/CichlidSequencingData/` (reads, BAMs, genomes, analysis run folders). Google Sheets: `SampleDatabase` (tabs `SampleDatabase`, `DNAReadsDatabase`), `AlignmentDatabase` (one tab per genome), `DownloadDatabase` (one tab per Run Info file). The reference is `Mzebra_GT3_NCBI` unless a module says otherwise.

## 5. Shared library: `common/` (package `cichlid_common`)

| Part | Content | Today's code |
|---|---|---|
| File manager | local/cloud paths per machine, rclone transfers with verification, Google authentication | `helper_modules/file_manager.py` |
| Sheets | cell-level, batched, journalled writes; claims; backends for Google and for tests | `download_manifest.py` (manifest, `GspreadBackend`, `MemoryBackend`) |
| Runs | immutable run folders, manifests (git commit, parameters, versions, checksums), `LATEST` | `YH_SD_Paper/yhsd/runner.py` |
| Genome registry | all genomes and their accessions, known misassemblies, LG-to-accession map | `YH_SD_Paper/config/genomes.yaml`, `genomes_registry.py` |
| Sample table | typed access to `SampleDatabase` (sex, species, ecogroup, brood, parents, genotype columns) | `sample_table`, `Project.sample_table()` |
| Tool runners | subprocess with injectable fakes, resumable downloads, md5, disk ledger | `ena_fetch.py`, `disk_ledger.py` |

**Constraint.** The alignment run on the server uses an older `file_manager.py` that can overwrite another machine's sheet changes. `common/` therefore re-exports the existing file manager at first, and replaces it only after that run has finished and every machine has been updated.

---

## 6. Modules

### 01_structural_landmarks: centromeres and inversion locations (LG2, 5, 9, 10, 11, 13, 20)

**Goal.** For each target chromosome, the centromere and every inversion between the reference and the other assemblies: breakpoints with uncertainty intervals, size, orientation, pericentric or paracentric, and which assemblies support it.

**Inputs.** Genome assemblies (registry), repeat annotation, published centromere table (Chromosome Research 2026, Table 3) for validation, optional read evidence from 04.

**Outputs.** The `landmarks/<genome>/` tables of section 4; to-scale chromosome figure data; a record of which regions of the reference are known misassemblies (currently LG2 4.0 to 14.4 Mb; LG9 9.7 to 12.3 Mb is suspected, not declared).

**Interface.** `landmarks centromeres arrays|scan|lift`, `landmarks inversions align|call`, `landmarks report`.

**Current state.** `centromeres` (satellite-array method, density scan, lift of published centromeres) and `inversion_breakpoints` (assembly alignment, interval-aware consensus, tandem segments) exist for LG2, 5, 9, 10, 11, 13 in `YH_SD_Paper`; run 20261105 gave the breakpoint table and centromeres. The code has been tested on synthetic data and the run on the cluster; the LG9 four-assembly consensus, the LG11 double inversion and the corrected fLetChi2 configuration still need the re-run. **LG20 is not covered yet.**

**Work.**
1. Add LG20 (assemblies carrying it, expected orientation, centromere). **(open: what is known about LG20 and which species or sex system motivated it.)**
2. Decide the LG9 suspected misassembly with the published-assembly check; declare or dismiss it.
3. Resolve the LG13 centromere ambiguity (two comparable satellite arrays).
4. A second assembly for LG2, 10 and 13 (today: YH only) so the calls are not single-assembly.
5. Write the schemas of the output tables and the per-inversion breakpoint windows used by 03 (slices) and 04 (read evidence).
6. Confirm GT3 and GT3_NCBI coordinate equivalence; confirm the orientation convention per chromosome (`reference_is_standard`, unset for LG5).

**Done when.** Every target chromosome has a centromere and an inversion table row with an interval for each breakpoint; the six existing chromosomes reproduce run 20261105 within the stated uncertainty; the centromeres agree with the published lift or the disagreement is documented; the tables validate against their schemas; tests pass offline.

### 02_download_reads: ENA / SRA download and the download database

**Goal.** Get every needed run into Dropbox as a uBAM with adapters marked, recording what happened, without anyone babysitting it.

**Inputs.** Run Info files (tabs of `DownloadDatabase`), `DNAReadsDatabase`, ENA's portal API for the file lists, NCBI's AWS copy of SRA.

**Outputs.** uBAM + metrics on Dropbox; rows in `DNAReadsDatabase`; the status of every run in `DownloadDatabase` (status, when, source, machine, attempts, failure kind, integrity checks, size, speed).

**Interface.** `download_reads run --manifest TAB ...`, `download_reads manifest import|status|reconcile|release`, `download_reads bench`, `download_reads add-pending`.

**Current state (v9c, `lab_repo_patch_v9c_ena`, 131 tests with fakes).** Two sources (ENA with Aspera or https, and SRA from AWS with `aria2c`, `vdb-validate`, `fasterq-dump`), tried in turn per run; every SRA run is checked (integrity, read count against ENA, first reads against ENA's FASTQ); Aspera probe and stall detection; two pools (downloads, conversions) with an adaptive download limit; claims, cell-level sheet writes and a journal; `--no-db` with `--add-pending`; `--max-runs`. Measured: ENA Aspera and https are very slow from three networks; the AWS copy is fast (11 MB/s with 8 connections on D07; 40 MB/s single connection on utaka and mzebra); a run's SRA file is 70% of its FASTQ; one run took 470 s to download, 38 s to convert to FASTQ and 12 min to become a uBAM. One real run completed end to end on D07. Pilots on utaka and mzebra and a scaling test are outstanding. **(unverified: parallel-download scaling; behaviour of the sheet with several machines at once.)**

**Work.**
1. Finish the pilots; set `--dl-workers`, `--workers`, `--sra-connections` per machine.
2. Run the campaign (328 pending runs, about 2.5 TB to transfer), then `--add-pending` once from one machine and `reconcile`.
3. Decide how to use the 0.09% count excess (SRA has more spots than ENA's `read_count` for ERR14734229).
4. Move from `helper_modules/` into this folder with thin wrappers under the old script names for one release (section 9).
5. Retire `downloadShortReadData.py` and `grabENA.py` to `legacy/`.

**Done when.** The sheet shows no `pending` runs apart from deliberate skips; every uBAM on Dropbox has a database row and a `DownloadDatabase` row; a second invocation does nothing; a failure leaves no phantom rows and no leftover files; tests pass offline.

### 03_align_pipeline: align, mark duplicates, split, slice, QC, upload

**Goal.** Turn uBAMs into aligned, duplicate-marked BAMs for a registered genome, with the splits and slices later modules need, and the per-sample measurements of 04, in one pass while the BAM is on local disk.

**Inputs.** uBAMs (02), the reference and its minimap2 index, `slice_regions` (from 01), the panel definitions (from 04).

**Outputs.** `Bamfiles/<genome>/<sample>/` (BAM, index, splits, slices, `qc.json`, `<sample>.panel.json`); one row per sample in `AlignmentDatabase`; no uBAM kept after use.

**Interface.** `align_pipeline run <genome> [-e ecogroup] [-u subgroups] ...`, `align_pipeline setup-database`, `align_pipeline genotype-only <samples>` (download a BAM, run 04's measurement, keep only the JSON, delete the BAM).

**Current state (`lab_repo_patch_v6`, 76 tests).** Per-sample pipeline without batch barriers (download, align, mark duplicates, split, stats, slice, QC, upload, database, cleanup), continuous space-aware downloads, a separate upload pool, drain-first scheduling, atomic database append with read-back, status lines every 10 minutes, `--setup-database`, `--include-excluded`; whole LG10 sliced; HaplotypeCaller removed. The old `alignFastQReads.py` is untouched. **(unverified: real-data validation.)** Validation on one sample with `validate_against_old.py` (seven BAMs must be identical) and the small-batch run have not been done; a BC run on utaka is still in its HaplotypeCaller stage and the old file manager is still running.

**Work.**
1. Validate on real data as above; fix whatever it shows.
2. Add the per-sample panel stage that calls `sd_genotypes panel measure` as a subprocess on the local BAM (no logic copied; 04 owns the definitions).
3. Add the genotype-only mode for the re-download campaign (all existing BAMs are downloaded once more for LG5 and B).
4. Split and slice definitions read from 01's tables, not hard-coded.
5. Two-machine procedure for the mbuna samples (about 477, 7 TB of uBAMs), documented and rehearsed.
6. Retire the old script after the mbuna run.

**Done when.** One test sample produces BAMs identical to the old pipeline; a batch of mbuna runs without intervention on two machines; every sample has an `AlignmentDatabase` row and a panel JSON; a killed run resumes without redoing finished samples; tests pass offline.

### 04_genotype_sd_systems: gsdf duplication, transposon insertions, inversions, B chromosome

**Goal.** Per-sample genotypes of the SD systems, using user-provided definitions, normally produced by 03 while the BAM is local, and re-runnable on samples that missed it.

**Inputs.** A local indexed BAM (or a stored slice); definitions: the gsdf region and calibration, the transposon insertion definitions for gsdf and chr10 (**provided by the user**), inversion breakpoint windows (01), B-chromosome core blocks (its own discovery step).

**Outputs.** `<sample>.panel.json`; call tables; columns in `SampleDatabase` (for example `gsdf number`, `TE1`, `TE2`, `Inversion*`, B status) with a source and flag column for each.

**Interface.** `sd_genotypes panel measure --bam B --sample S --out-json J [--modules ...]`; `sd_genotypes call <module> --measurements DIR`; `sd_genotypes validate <module>`; `sd_genotypes bchrom prepare|call|compare`; `sd_genotypes sheet-update`.

**Current state.**
- *gsdf duplication* (`panel` and `gsdf`, tests with fakes): raw 5 kb bins and junction reads stored, copy number called with the calibration of the MC pedigree (0.530 and 0.268), flags for low depth, non-integer, range, junction disagreement; validation against the existing 910 `gsdf number` values and sex-by-species tables. **(unverified: real BAMs; the calibrated range of 2 to 4 copies is my reading.)**
- *B chromosome* (`bchrom`): downloads Clark et al.'s public reads, aligns them with the lab recipe, finds blocks and a core set in reference coordinates, evaluates each carrier leave-one-out, separates tier A (normal male-pool coverage) from tier B (low baseline). Not yet run on real data. The per-sample genotype step (a score over the core blocks, a threshold from the data) is not written.
- *Inversions*: Blumer et al. (Table S7) give chr2, 9, 10, 11, 13 for 1,101 samples; all 440 mbuna among them are 0/0 on all five. These are **validation labels**, not an input. The LG5 (orange blotch / ZW) inversion is polymorphic in the mbuna and has no labels; PRJNA808791 (Behrens et al.) has BB and OB pools of known composition; `bchrom prepare` keeps slices around both LG5 breakpoints for them. A breakpoint-read genotyper is designed but not written.
- *Transposon insertions*: previously produced by a separate tool (`sv-genotype`, columns `TE1`, `TE2`); to be re-implemented here from user-provided definitions. **(open: the definitions, the evidence to use, and which samples.)**

**Work.** Finish the B-chromosome run and the genotype score; write the inversion genotyper (breakpoint-read evidence with windows, validated on Blumer labels, BC families and the pools); implement the transposon genotyper; validate gsdf on about 40 samples with known values; one `panel` JSON schema for all modules; the sheet write-back; the retrospective mode (BAM or slice in, JSON out).

**Done when.** gsdf: unflagged calls agree with the existing column for at least 98% of the samples and females are 2 and males 3 in the MC family; B chromosome: carriers score above the threshold leave-one-out, controls and mbuna males below it; LG5: BB pools about 0 and OB pools about 0.5 inverted; inversions: at least 98% agreement with Blumer on labelled samples, with hets reported separately; transposons: agreement with the old `TE1`/`TE2` columns where they exist; every measurement carries its definition hash.

### 05_variant_calling: VCFs for families and species

**Goal.** Reproducible VCFs for the cohorts that the SD analyses need: broods and families with parents, and groups of related species where the sex system is unknown.

**Inputs.** BAMs and slices (03), cohort definitions from the sample table (species, brood, parents, sex, ecogroup), a regions list (whole genome or candidate regions), the reference.

**Outputs.** One bgzipped, tabix-indexed VCF per cohort with a manifest (samples, caller, version, filters, regions, sites and their attrition), plus per-sample coverage tracks needed downstream.

**Interface.** `variant_calling cohort build <name>`, `variant_calling cohort list`, `variant_calling qc <name>`.

**Current state.** `callset build` in `YH_SD_Paper` (callsets in `config/callsets.yaml`: `bc_family`, `yh_broods`, `mc_broods`) from the BC genotyping work. HaplotypeCaller was removed from the alignment pipeline on purpose. **(open: the caller, joint versus per-family calling, region-restricted versus whole genome, the filter set, and how unknown species are grouped.)**

**Work.** Define the cohort file format; choose and benchmark a caller on one brood against the existing VCFs; cohort builder that needs no manual sample lists; QC (missingness, depth, Ti/Tv, Mendelian errors in families); the manifest.

**Done when.** The BC and brood cohorts rebuilt from the sheet reproduce the existing VCFs' calls at the sites they share (stated concordance); every cohort VCF has a manifest; Mendelian-error rates are reported.

### 06_sd_region_discovery: SD regions in many species

**Goal.** For each species or cohort, the regions associated with sex, with the statistics and the evidence that identifies the kind of system (XY or ZW, inversion, dosage).

**Inputs.** VCFs (05), the sample table, 04's calls, 01's landmarks.

**Outputs.** Candidate regions per cohort (BED and statistics: association, Fst, heterozygosity by sex, male/female coverage ratio, parent-offspring inheritance) with plots, and a species-by-region summary.

**Interface.** `sd_regions scan <cohort>`, `sd_regions summarize`.

**Current state.** `parental_sex_region` (inheritance, Fst, heterozygosity, depth ratio per family) and `sex_association fisher` exist and were tested on synthetic data; the Fisher tests have not run on real data. **(open: species list and the minimum sample sizes; multiple-testing and permutation scheme.)**

**Work.** One statistics core for all methods; permutation thresholds; sex-specific coverage; classification of the system; a result table that 09 and 10 read.

**Done when.** It recovers the known systems (the chr10 Y inversion in the YH lineage, the gsdf dosage pattern in *M. conophoros*, the LG5 locus in the orange-blotch species) as positive controls and reports nothing in species known to have none.

### 07_bc_analysis: backcross genotyping, recombination, transmission distortion, QTL mapping, fine mapping

**Goal.** Genotype the BC individuals from their reads, then analyse them: recombination, transmission distortion, QTL mapping, and fine mapping of the intervals.

**Inputs.** BAMs (03), cohort VCF and sites (05), pedigree (sample table), 01's landmarks, traits.

**Outputs.** Genotype and posterior matrices, breakpoints, per-fish quality, recombination maps (including suppression in inversions), TD tests, QTL scans, fine-mapping tables.

**Interface.** `bc_analysis genotype <step>`, `bc_analysis recombination`, `bc_analysis td`, `bc_analysis qtl`, `bc_analysis finemap`.

**Current state.** `bc_genotyping` (about 266 fish, `OutputGenotypes.csv`, 9,336 markers; steps `provenance`, `flagstat`, `original`, `v2-sites`, `v2-counts`, `v2-hmm`, `validate`, `compare`) is the part to **fix**; the analyses after it are not written.

**Work.** Fix the genotyping (the v2 site counts and the HMM against the original calls); recombination analysis including inversions as suppressed blocks; TD; QTL mapping (**open: R/qtl2 or own code**); fine mapping.

**Done when.** The genotyper agrees with the original at the shared markers, recovers the expected F1 phase and the known inversion state of each fish, and shows suppressed recombination across the chr10 inversion; each analysis reproduces from a single command.

### 08_ase_variants: variants for allele-specific expression

**Goal.** For a gene in a family, a list of variants that can be used to tell the two alleles apart in RNA-seq, with the information needed to trust them.

**Inputs.** VCFs (05), annotation, mappability, parent and F1 genotypes.

**Outputs.** Per-gene, per-family tables: position, alleles, parental origin, transcript context, coverage, mappability, reference-bias risk, rank.

**Interface.** `ase_variants suggest --gene G --family F`, `ase_variants batch --genes FILE`.

**Current state.** Not started. **(open: the annotation source, the minimum coverage, and how reference bias is judged.)**

**Done when.** For ten test genes it returns the same variants as a hand-made list, flags variants in repeats or paralogs, and says when a gene has none.

### 09_candidate_qtn_webserver: candidate causal variants and the web server

**Goal.** From the mapped intervals, find the genetic differences between parental assemblies, annotate them, check how their alleles and genotypes vary in other species with and without the SD system, and present it.

**Inputs.** Assemblies and alignments between them, QTL and fine-mapping results (07), 04's genotypes, the cohort VCFs (05), SD regions (06), annotation.

**Outputs.** Candidate-QTN tables (SNPs, indels, structural variants, transposon insertions with their effect), allele and genotype frequency tables across species, correlation with SD, and a read-only web server over them.

**Interface.** `candidate_qtn diff --interval I`, `candidate_qtn frequencies`, `candidate_qtn build-site`, `candidate_qtn serve`.

**Current state.** Not started. **(open: hosting and who may see it; a static site over precomputed files, or an API over a database; whether a genome browser (igv.js) is embedded.)**

**Done when.** For one known interval the page lists every difference between the assemblies, shows its frequency by species, and links to the evidence; rebuilding the site from the run folders is a single command.

### 10_paper_figures: figures and data tables

**Goal.** Every figure and supplementary table of the paper regenerated by one command, each naming the run folder it was built from.

**Inputs.** Run folders of 01 to 09.

**Outputs.** SVG and PDF figures, TSV and XLSX tables, and a `figure_sources.tsv` mapping each to its inputs.

**Interface.** `paper_figures build <figure>`, `paper_figures build --all`, `paper_figures sources`.

**Current state.** The to-scale inversion figure (six chromosomes) and the introduction figure (panel B and four icons) exist and are drawn from real coordinates. Figures 2 to 4 (the transposon insertions) come after 04. **(open: the figure list.)**

**Done when.** `build --all` on a clean machine reproduces every figure and table from Dropbox run folders, and a changed upstream run is detected as a stale figure.

---

## 7. Machines and environments

| Machine | Role | Notes |
|---|---|---|
| utaka | alignment and heavy work (96 cores) | runs the old file manager; downloads only from a separate clone |
| mzebra | downloads and conversions (96 cores, 503 GB RAM, 3.5 TB free) | login `bio-mcgrath-dropbox`; environment `enadl` works |
| D07 (workstation) | downloads, development (32 cores, 62 GB RAM) | Aspera blocked; SRA about 11 MB/s |
| yellowhead | validation of the pipeline | student machine |
| laptop | development, tests, publishing | tests run offline |

Each module has an `env.yml`. A known conflict: `aria2` and `sra-tools` cannot be solved in one conda environment on macOS; on Linux they can.

## 8. Order of work

1. **Now, in parallel:** finish 02 (pilots, campaign); validate 03 on real data; finish 04's definitions before any BAM is downloaded again.
2. **Before the mbuna alignment run:** 04 (gsdf validated, B-chromosome blocks and genotype, LG5 genotyper, transposon definitions) so every BAM is measured once; 01 finished for LG20 and the open items.
3. **Then:** 05, then 06 and 08 (independent of each other), then 07.
4. **Last:** 09 and 10, which read everything.

Critical path to the mbuna run: 02 complete, 03 validated, the 04 panel defined.

## 9. Migration from today's layout

| Today | Becomes |
|---|---|
| `YH_SD_Paper/yhsd/analyses/centromeres`, `inversion_breakpoints` | 01 |
| `helper_modules/ena_fetch.py`, `sra_fetch.py`, `ena_scheduler.py`, `grab_ena_v2.py`, `download_manifest.py`; `downloadShortReadData_v2.py`, `downloadManifest.py`, `unit_scripts/ena_bench.py`, `patch_file_manager.py` | 02 (and `common/` for the sheet backend) |
| `alignFastQReads_v2.py`, `seq_pipeline.py`, `bam_split.py`, `bam_slices.py`, `bam_stats.py`, `alignment_db.py`, `disk_ledger.py`, `download_manager.py` | 03 (`disk_ledger` also in `common/`) |
| `YH_SD_Paper/yhsd/analyses/panel`, `gsdf`, `bchrom`, `sample_table` | 04 |
| `YH_SD_Paper` callset builder (`callset build`) | 05 |
| `parental_sex_region`, `sex_association` | 06 |
| `bc_genotyping` | 07 |
| `YH_SD_Paper/figures` | 10 |
| `YH_SD_Paper/yhsd/runner.py`, `common.py`, `genomes_registry.py`, `config/genomes.yaml` | `common/` |
| `downloadShortReadData.py`, `grabENA.py`, `alignFastQReads.py` | `legacy/` once retired |

**Rules for the move.** It happens in one change, after the pilots and while no download or alignment is running on any machine; the downloader launches its workers by file path, and the alignment run imports the old modules. Thin wrapper scripts under today's names stay for one release. The test suites are the safety net: they import by name and fail loudly if anything was missed. The package name `yhsd` is retired in favour of the module packages above.

## 10. Open decisions

1. LG20: what is known, and what it adds to the paper. (01)
2. The transposon insertion definitions and the evidence to use. (04)
3. Variant caller, joint or per-family calling, filters, region scope. (05)
4. Species list and sample-size rules for SD discovery. (06)
5. QTL software. (07)
6. Annotation and mappability sources for ASE. (08)
7. Web server hosting, audience and data store. (09)
8. The figure and table list. (10)
9. When to run the migration, and who updates each machine.
10. Whether `SPEC.md` is edited directly in the repo, or changes go through a decision log in `docs/`.
