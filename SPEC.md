# CichlidSRSequencing: specification

Status: v3 (2026-11-07). Owner: McGrath lab. This file lives at the repository root and is the source of truth for what each part of the repository is for, what it consumes, what it produces, and when it is done.
Where a statement is not yet verified it is marked **(unverified)**; where a decision is still to be made it is marked **(open)**.

## 0. What changed since v1

* **Module 02 (downloads) is finished for the two Run Info files.** 1,832 of 1,833 runs are in the database; 334 of them came from NCBI's copy of SRA because ENA was unusable (section 6, module 02).
* **Module 03 (alignment) is set against v6.** v6 (per-sample flow, no GVCFs, `samtools stats`, atomic row append) is written and has offline tests but has never run on real data. The changes agreed on 2026-11-07 are: a leaner ledger with DuplicationRate and Species added, a rerun rule that takes RunIDs into account, updates in place, a journal, claims, supervised processes (section 6, module 03). The ledger columns still have to be reconciled with the ~20 statistics columns v6 already writes.
* **The SampleDatabase has a section of its own** (section 7): the Blumer update, the species policy, the new `SRASpecies` column, and the open question of 186 samples that have reads but no sample row.
* **New tools are in the repository**, still in the original layout (section 3): the downloader with an ENA and an SRA source, the DownloadDatabase tool, `verify_ubam`, `audit_alignments`, `register_alignments`, `sampledb`, a retry of every Google call.
* **The reorganised layout of v1 exists as a tested tree but has not been applied** (section 10).
* **Principles 10 to 14 are new, and 6 is refined:** a whole-sheet rewrite is allowed only when a single writer is guaranteed (principle 11).
* **A decision log** with the measurements behind the decisions (section 11) and a text for the methods (appendix A).

---

## 1. Purpose

Take public and in-house short-read data for Lake Malawi cichlids from a read archive to a set of results about sex determination (SD): where the SD loci are, what they are made of (inversions, transposon insertions, a gene duplication, a B chromosome), how they behave in backcross families, and which variants are candidate causal changes. Everything that goes into the paper must be reproducible from versioned code and recorded inputs.

The work is ten modules. Each is one folder in the target layout. They communicate only through **data contracts** (section 4), never by importing each other's internals, so any module can be re-run, replaced or tested alone. A module may import the package of an upstream module it declares a dependency on (04 uses `download_reads.sources`; 06 and 07 use `variant_calling.callset`; all use `cichlid_common`), never its internals.

## 2. Principles that apply to every module

1. **One folder, one job.** A module has a README, a command-line entry point, a config folder, tests, and an environment file. Numbered folders are for humans; the importable package inside has no number.
2. **Raw measurements are stored; calls are made later.** Anything computed from a BAM or VCF stores raw numbers first (JSON/TSV) and makes calls in a separate, offline step, so a changed threshold or calibration never needs the big files again.
3. **Immutable, recorded runs.** Every step writes a run folder (inputs, parameters, tool versions, git commit, outputs, checksums) to Dropbox with a `LATEST` pointer, as the current `Run` helper does. A figure or table must name the run it came from.
4. **Definitions are hashed.** A measurement carries a hash of what was measured (region, windows, thresholds that define it). A caller refuses to mix measurements made under different definitions.
5. **Fail loudly, never silently.** Non-zero exit codes for failures, status files saying which stage failed and why, no database row for work that did not finish, no `sys.exit()` with status 0 after an error.
6. **Safe on shared resources.** Google Sheets are written cell-by-cell (or by appending one row) for the changed rows; a whole-worksheet rewrite is allowed only under principle 11; runs are claimed before they are processed; Dropbox uploads are verified; disk is reserved before work starts.
7. **Tests without the real tools.** Logic is tested with fakes (fake `samtools`, in-memory sheet, planted signal) so that the whole suite runs in seconds on a laptop; the real tools are exercised by small, documented pilots.
8. **Real data decides.** Assumptions are written down together with the measurement that supports them (for example: ENA's `read_count` means *reads* for some runs and *pairs* for others, so no code assumes one meaning).
9. **Nothing is kept that need not be.** Reads and large intermediates are deleted once the compact result exists; the result records what was deleted.
10. **Report before write.** Every tool that changes a sheet writes a report of what it would change first (`--dry-run` or the default mode), saves a backup of what it overwrites, writes once, and reads the result back cell by cell.
11. **Whole-sheet rewrites only with a single writer.** `FileManager` rewrites whole worksheets from a DataFrame (formulas become values, concurrent writers overwrite each other). Bulk tools that do so say so, and are run only when nothing else writes that sheet. Anything several machines may do at once uses cell-level writes and claims.
12. **A convenience layer must never be able to break the program that imports it.** The Google retry is installed at import by `file_manager.py`; its installation cannot raise.
13. **A hung process is a bug.** Every subprocess is supervised: a timeout or a stall limit, a process group that is killed with its parent, a status someone can read. (On utaka an `rclone` child of a dead downloader sat for 88 days.)
14. **An archive is checked, not trusted.** Reads from a second archive are compared with the first (read counts, the first reads' sequence and quality) and the result is recorded per run.

## 3. Repository layout

### Target layout (module folders)

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
  envs/  docs/  legacy/
```

Inside each module: `README.md`, the package (`cli.py` is the only entry point: `python -m <package> <command>`), `config/`, `tests/`, `env.yml`. Every command that writes results accepts `--out-dir` (plain files, no Dropbox).

### Where the code is today (branch `YH_SD_Paper`, original layout)

| Piece | Location in `cichlid_sr_sequencing/` |
|---|---|
| Downloader (ENA and SRA sources, `--manifest`, `--add-pending [--bulk]`, `--max-runs`) | `downloadShortReadData_v2.py` with `helper_modules/ena_fetch.py`, `sra_fetch.py`, `ena_scheduler.py`, `grab_ena_v2.py` |
| DownloadDatabase tool (import, status, reconcile, release) | `downloadManifest.py`, `helper_modules/download_manifest.py` |
| Retry of every Google call | `helper_modules/gspread_retry.py`, installed by `helper_modules/file_manager.py` at import |
| File manager patcher (download-database methods; `--retry-only`) | `unit_scripts/patch_file_manager.py` |
| uBAM recheck against ENA | `unit_scripts/verify_ubam.py` |
| Alignment audit (read-only) and registration of uploaded-but-unregistered alignments | `unit_scripts/audit_alignments.py`, `unit_scripts/register_alignments.py` |
| SampleDatabase checks and the Blumer update | `unit_scripts/sampledb.py`, input tables in `data/blumer/` |
| Network benchmark | `unit_scripts/ena_bench.py` |
| Tests (offline) | `tests/` |
| Documentation of the tools | `docs/DOWNLOAD_AND_DATABASE_TOOLS.md` |
| Alignment pipeline, old: `alignFastQReads.py`, `helper_modules/alignment_worker.py`, with the download queue of patches v3/v4 (`helper_modules/download_queue.py`) | `./`, `helper_modules/` |
| Alignment pipeline, v6 (offline tests, never run on real data): `alignFastQReads_v2.py` with `seq_pipeline.py`, `bam_split.py`, `bam_slices.py`, `bam_stats.py`, `alignment_db.py`, `disk_ledger.py`, `download_manager.py`, `slice_regions/`, `unit_scripts/validate_against_old.py`; described in `PIPELINE_V2_README.md` | `./`, `helper_modules/` |
| Analyses of the paper (centromeres, inversion breakpoints, panel, gsdf, bchrom, sample table, BC genotyping, callsets, figures) | `YH_SD_Paper/` (package `yhsd`) |

The utaka checkout runs a different branch, `mcgrath_dev` (last commit 2026-09-21), with the old alignment pipeline and an older `file_manager.py` (410 lines, no download-database methods, no retry).

## 4. Dependencies, data contracts and shared stores

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
| 01 | `landmarks/<genome>/centromeres.tsv`, `inversions.tsv`, `inversion_breakpoint_windows.bed`, `known_misassemblies.tsv` | 03, 04, 07, 09, 10 |
| 02 | uBAM + adapter metrics per run on Dropbox (`Reads/<project>/<run>.unmapped_marked_adapters.bam`); rows in `DNAReadsDatabase` and `SampleDatabase`; `DownloadDatabase` tabs | 03 |
| 03 | `Bamfiles/<genome>/<sample>/<sample>.all.bam` + index, splits, slices, `qc.json`; one row per sample in `AlignmentDatabase` | 04, 05, 07 |
| 04 | `<sample>.panel.json` (raw measurements with definition hash); call tables; calls written back to `SampleDatabase` | 05, 06, 07, 09, 10 |
| 05 | per-cohort VCF (bgzip + tabix) with a manifest | 06, 07, 08, 09 |
| 06 | SD candidate regions per species | 09, 10 |
| 07 | BC genotype matrices, breakpoints, recombination maps, TD tables, QTL results | 09, 10 |
| 08 | per-gene, per-family ASE variant lists | 10 |
| 09 | candidate-QTN tables, frequency tables, the web server's data store | 10 |

### Shared stores

**Dropbox** `Apps/CichlidSequencingData/` (`ptm_dropbox:/COS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/`). A local path maps to its Dropbox path by replacing the master directory (`/Data/mcgrath-lab/Temp/CichlidSequencingData/` on utaka, `$HOME/Temp/CichlidSequencingData/` on mzebra); code never builds a Dropbox path by hand. Folders: `Reads/<project>/`, `Reads/ReadDownloadFiles/` (Run Info CSVs), `Bamfiles/<genome>/<sample>/`, `Genomes/<genome>/`, `Temp/` (never uploaded), and the analysis project folder `YH_SD_Paper/<genome>/<analysis>/<step>/<run_id>/` with `LATEST` pointers. The full description is in `docs/data_layout.md` (reorganised tree).

**Google sheets** (opened with the service account; ids are attributes of `FileManager`):

| Spreadsheet (attribute) | Tabs | Content |
|---|---|---|
| SampleDatabase (`s_ID`) | `SampleDatabase`, `DNAReadsDatabase`, genome tabs | samples and their metadata; one row per sequencing run |
| AlignmentDatabase (`a_ID`) | one per genome version | one row per aligned sample (section 6, module 03) |
| DownloadDatabase (`dl_ID`) | one per Run Info file (`Blumer`, `Malinsky`) | one row per run: its metadata and what happened to it |

`SampleDatabase` tab, 22 columns on 2026-11-07: `SampleID`, `Sex`, `Species`, `DoB`, `BroodID`, `Parents`, `Ecogroup`, `Subgroup`, `LabReared`, `Category`, `gsdf number`, `TE1`, `TE2`, `BlumerSequenceID`, `BlumerSpecimenID`, `Inversion2`, `Inversion9`, `Inversion10`, `Inversion11`, `Inversion13`, `InversionSource`, `SRASpecies`, followed by empty trailing columns.
`DNAReadsDatabase`: `SampleID`, `ProjectID`, `RunID`, `ReadLength` (the Run Info's AvgSpotLen: both mates), `TotalBases`, `Instrument`, `LibraryID`, `LibraryLayout`, `LibrarySource`, `Platform`, `FileLocations` (relative to `Reads/`), `FileSize`.
`DownloadDatabase` tab: the 14 Run Info columns, then `DownloadStatus` (pending, claimed, downloading, converting, uploaded, in_database, failed, skipped), `InDatabase`, `DownloadedOn`, `Source` (ENA-aspera, ENA-https, SRA-aws, legacy), `Machine`, `ClaimedBy`, `ClaimedAt`, `UpdatedAt`, `Attempts`, `FailureKind` (permanent, transient, throttle, source_unavailable), `LastError`, `Verified`, `QualityCheck`, `ReadCountCheck`, `ReadsSource`, `ReadsInUBAM`, `UBamBytes`, `UBamPath`, `DownloadSeconds`, `MBPerSec`, `ConvertSeconds`, `UploadSeconds`, `PipelineVersion`, `Notes`, `Skip`.

The reference is `Mzebra_GT3_NCBI` unless a module says otherwise.

## 5. Shared library: `common/` (package `cichlid_common`)

| Part | Content | Today's code |
|---|---|---|
| File manager | local/cloud paths per machine, rclone transfers with verification, Google authentication | `helper_modules/file_manager.py` |
| Google retry | retry of connection resets, timeouts and HTTP 429 / 5xx in every gspread call; never repeats a request that may already have reached Google unless it is idempotent | `helper_modules/gspread_retry.py` |
| Sheets | cell-level, batched, journalled writes; claims; backends for Google and for tests | `helper_modules/download_manifest.py` |
| Runs | immutable run folders, manifests, `LATEST` | `YH_SD_Paper/yhsd/runner.py` |
| Genome registry | all genomes, their accessions, known misassemblies, annotations | `YH_SD_Paper/config/genomes.yaml`, `genomes_registry.py` |
| Sample table | typed access to `SampleDatabase` | `sample_table`, `unit_scripts/sampledb.py` |
| Data manifest | one index of every run folder and file on Dropbox: the newest run per step, steps expected but missing, runs that need attention, what exists only locally | reorganised tree only (`common/cichlid_common/data_manifest.py`) |

**File manager facts that cost time and are now designed around.** `returnCloudFiles` is an uncached `rclone lsf` of a whole folder and `returnFileSize` an `rclone size`: one call of each per run, minutes for a few hundred runs; bulk code uses one `rclone lsjson` per folder. `setSampleDatabase` rewrites both worksheets completely and does not enlarge the grid; `setAlignmentDatabase` retried its write three times but `readAlignmentDatabase` never retried.

---

## 6. Modules

### 01_structural_landmarks: centromeres and inversion locations (LG2, 5, 9, 10, 11, 13, 20)

**Goal.** For each target chromosome, the centromere and every inversion between the reference and the other assemblies: breakpoints with uncertainty intervals, size, orientation, pericentric or paracentric, and which assemblies support it.

**Inputs.** Genome assemblies (registry), repeat annotation, published centromere table (Chromosome Research 2026, Table 3) for validation, optional read evidence from 04.

**Outputs.** The `landmarks/<genome>/` tables of section 4; to-scale chromosome figure data; a record of known misassemblies of the reference (currently LG2 4.0 to 14.4 Mb; LG9 9.7 to 12.3 Mb is suspected, not declared).

**Interface.** `landmarks centromeres arrays|scan|lift`, `landmarks inversions align|call`, `landmarks report`.

**Current state.** `centromeres` and `inversion_breakpoints` exist for LG2, 5, 9, 10, 11, 13 in `YH_SD_Paper`; run 20261105 gave the breakpoint table and centromeres. The code is tested on synthetic data and the run on the cluster; the LG9 four-assembly consensus, the LG11 double inversion and the corrected fLetChi2 configuration still need the re-run. **LG20 is not covered yet.** The whole-genome alignments (PAF) of each assembly to the reference are kept only in the `tmp/` folder of the machine that ran the step (utaka); only the blocks on the target chromosomes are uploaded (`targets.paf`, `contig_map.tsv`).

**Work.** (1) Add LG20 **(open: what is known about LG20 and which species or sex system motivated it)**. (2) Decide the LG9 suspected misassembly. (3) Resolve the LG13 centromere ambiguity. (4) A second assembly for LG2, 10 and 13. (5) Write the schemas of the output tables and the per-inversion breakpoint windows used by 03 and 04. (6) Confirm GT3 and GT3_NCBI coordinate equivalence and the orientation convention per chromosome (`reference_is_standard`, unset for LG5).

**Done when.** Every target chromosome has a centromere and an inversion table row with an interval for each breakpoint; the six existing chromosomes reproduce run 20261105 within the stated uncertainty; the centromeres agree with the published lift or the disagreement is documented; the tables validate against their schemas; tests pass offline.

### 02_download_reads: ENA / SRA download and the download database

**Goal.** Get every needed run into Dropbox as a uBAM with adapters marked, recording what happened, without anyone babysitting it.

**Inputs.** Run Info files (tabs of `DownloadDatabase`), `DNAReadsDatabase`, ENA's portal API for the file lists, NCBI's copy of SRA on AWS.

**Outputs.** uBAM + metrics on Dropbox; rows in `DNAReadsDatabase` and `SampleDatabase`; the status of every run in `DownloadDatabase`.

**Interface.** The downloader: `--manifest TAB ...` (or Run Info files, `--ena`, `--samples-file`), `--source ena|sra|auto`, `--max-runs N`, `--no-db`, `--dry-run`, `--add-pending FILE [--bulk]`. The sheet tool: `downloadManifest.py import|status|reconcile|release|tabs`. Also `unit_scripts/verify_ubam.py --skipped [--write]`.

**State: complete for the two Run Info files (2026-11-07).**

| | |
|---|---|
| Runs in `BlumerPaperData.csv` and `MalinskyPaperData.csv` | 1,160 and 673 (1,833) |
| In `DNAReadsDatabase` | 1,159 and 673 (1,832) |
| Not in the database | 1: ERR14734399, `Skip` = TRUE with a note (below) |
| Obtained by the earlier ENA scripts (`Source` legacy) | 857 and 641 (1,498) |
| Obtained in this campaign from SRA (`Source` SRA-aws) | 302 and 32 (334), among them 8 runs re-downloaded from SRA after their older uBAMs (7 on Dropbox without a database row, 1 pilot) were deleted |
| `DNAReadsDatabase` / `SampleDatabase` rows after the campaign | 2,763 / 2,154 |

**What was measured and decided** (details in section 11):
* ENA moved about 1 Mbit/s over Aspera and about 0.3 MB/s over https from three unrelated networks, at the same time as NCBI's copy ran at 17 to 40 MB/s per connection. Not a Georgia Tech firewall problem (unverified with IT; no request was made). The downloader therefore has two sources and `--source auto` (the Aspera probe orders them).
* On mzebra one SRA download ran at about 104 MB/s, three at once at about 150 MB/s in total; D07 reached about 11 MB/s. An SRA file is about 70% of a run's FASTQ size; `fasterq-dump` took 38 s on the pilot run; conversion to a uBAM took 19 to 23 minutes on mzebra (12 to 16 parallel workers; about 12 minutes on D07) and the upload to Dropbox 5 to 12 minutes at 20 to 25 MB/s per stream, about 4 streams at a time.
* Every SRA run is checked: `vdb-validate`; the read count against ENA's `read_count` (which means reads for some runs and pairs for others; both readings are tried; SRA has 0.01% to 1.71% MORE spots than ENA in every run measured, never fewer, accepted up to 5% and recorded in `ReadCountCheck`); and the first 1,000 reads (sequence and quality) against ENA's FASTQ. Of 327 runs at the time of the check, 323 were `identical-to-ENA` and 4 were accepted as `skipped (ENA unreachable)`.
* ERR14734399 was rejected: all of the first 1,000 sequences differ from ENA's (the cause was not investigated). Its ENA copy was not fetched (ENA's transfer speed); the run is excluded.
* The database rows are added with `--no-db` first and `--add-pending --bulk` at the end: one read and one write of the two sheets instead of two reads and two whole-sheet writes per run, which hit Google's 60 reads per minute (HTTP 429) after 26 rows, and one `rclone lsjson` per project folder instead of two rclone calls per run.

**Open.** (1) The 4 unverified runs: `verify_ubam.py --skipped` rechecks them from the finished uBAMs. (2) ERR14734399: leave excluded, or fetch from ENA when its transfer is usable. (3) The 4 failures of older tests on mzebra and the laptop are missing modules (`psutil`, `pysam`) of the alignment code, not of this module.

**Done when (met):** no `pending` runs apart from deliberate skips; every uBAM on Dropbox has a database row and a `DownloadDatabase` row; a second invocation does nothing; a failure leaves no phantom row; tests pass offline.

### 03_align_pipeline: align, mark duplicates, split, slice, QC, upload

**Goal.** Turn uBAMs into aligned, duplicate-marked BAMs for a registered genome, with the splits and slices later modules need and the per-sample measurements of 04, one sample at a time, with a ledger (the AlignmentDatabase) that says what exists, how it was made and from which reads.

**Inputs.** uBAMs (02), the reference and its minimap2 index, `slice_regions/<genome>.tsv` (01; today all of chromosome 10, `NC_135176.1`), the panel definitions (04), the SampleDatabase (selection and metadata).

**Outputs.** `Bamfiles/<genome>/<sample>/` (BAM, index, splits, slices, QC, later the panel JSON); one row per sample in the AlignmentDatabase; no uBAM kept after use.

**Interface.** Today: `alignFastQReads.py <genome> -e|-u|-p|-c|-s ...` (old) and `alignFastQReads_v2.py <genome> [--samples-in-flight N --align-threads N --uploads N --list --no-upload --setup-database ...]` (v6). Target: `align_pipeline run|setup-database|genotype-only`.

#### What exists

* **The old pipeline** (`alignFastQReads.py`, `alignment_worker.py`; utaka and yellowhead, the latter with the download queue of patches v3/v4): every stage is a barrier over all samples, alignment is one sample at a time with `minimap2 -t 96`, HaplotypeCaller GVCFs are made for every sample, and the upload and the database rows come last, one row at a time through `FileManager.addAlignmentRow`.
* **v6** (`alignFastQReads_v2.py`; replaces the old script, offline tests, **never run on real data**). Per `PIPELINE_V2_README.md`: each sample moves on its own (download, align, mark duplicates, split / stats / slice, qc, verified upload, database row, cleanup); several samples in flight with a shared core and RAM budget; space-aware continuous downloads with a disk ledger; a separate upload pool; stage checkpoints so a rerun resumes; `bash -o pipefail` in the alignment pipe; **no GVCFs**; one pass over the BAM for the six splits; `samtools stats` and `idxstats` instead of seven `gatk CountReads`; uBAMs deleted once the duplicate-marked BAM is verified; **one row appended atomically** (`AppendingDatabaseWriter`: append, read the SampleID column back, retry up to four times, warn on a duplicate) so two machines finishing at the same moment lose no row; `--setup-database` adds its columns to the header once; a status line every 10 minutes; samples with Ecogroup `Error` or `Exclude` are skipped unless `--include-excluded`. The recipe (SamToFastq, `minimap2 -x sr`, MergeBamAlignment, MarkDuplicates, the six split definitions) is unchanged on purpose.
* **Tools** (`unit_scripts/`): `audit_alignments.py` (read-only: the ledger against Dropbox against the local machine; classes ok, needs_row, partial_upload, needs_upload, local_incomplete, row_only, not_started); `register_alignments.py` (adds rows for uploaded alignments: counts with GATK as the old loop did, in parallel, one read and one write of the sheet, backup, read-back; `--check-existing` compares every computed column with existing rows; `--counter pysam|samtools` is faster; `--allow-unlisted` for samples missing from the SampleDatabase); `validate_against_old.py` (the seven BAMs of a sample from v6 against the old pipeline's); `gspread_retry.py`.

#### What happened, and why

| Observation | Source |
|---|---|
| An alignment run on utaka died in `addAlignmentRow` (`Connection reset by peer`): the old script uploads every sample first, then counts reads and adds the rows one by one; `readAlignmentDatabase` never retries; 216 samples were left uploaded and complete but without a row | traceback and `audit_alignments`, 2026-11-07 |
| A batch of about 45 samples on yellowhead (old script): download 291 min (77 uBAMs), alignment 1,971 min (one sample at a time; most 30 to 50 min, four 80 to 144 min), duplicate marking 76 min, splitting 56 min, then GVCFs, then upload and database | its log |
| Seven `gatk CountReads` per sample: three samples took 233 s in parallel, 178 would take hours one by one | `register_alignments` |
| Samples that are not in the SampleDatabase are skipped silently by the old loop ("Weird Error ... Skipping") | `alignment_worker.py` |
| A sample is done if it has a row; new runs of the same sample do not trigger a rerun (v6 also leaves a sample that has a row alone) | both selection logics |
| An `rclone` child of a dead downloader hung for 88 days | utaka |
| Audit of 2026-11-07 (genome `Mzebra_GT3_NCBI`): 1,073 ok; 216 uploaded without a row (178 Hybrid with BAMs on utaka, 38 not in the SampleDatabase and without local BAMs); 3 partial folders (`PM30_backups`, `SAMEA117806170`, `SAMEA2661230`); 1,078 never started (477 Mbuna, 45 AC, 23 Utaka, 20 Shallow_Benthic, 1 Exclude, 186 without a sample row, 326 with a blank Ecogroup before the Blumer update) | `audit_alignments` |

#### Agreed in discussion on 2026-11-07, set against v6

1. **No GVCF stage.** v6 has none; the old script still makes them. Downstream code builds them on demand: `callset.py` and the BC code call `ensure_gvcf` (a cloud GVCF if there is one, HaplotypeCaller from the BAM otherwise; `paper.yaml`: "only used if a GVCF is missing"). **(open: whether the 45 yellowhead samples keep the GVCF stage or the run is stopped and finished by hand: upload, then `register_alignments`.)**
2. **The ledger keeps** SampleID, GenomeVersion, RunIDs, Coverage, TotalReads, the version columns (minimap2, gatk, pysam) and BamSize; **adds** DuplicationRate and Species; **drops** the six split-read fractions. The discussion was held without reference to v6's own new columns, which are many (below), so the reconciliation is **open**.

| Column | Agreed or proposed | Relation to v6 |
|---|---|---|
| DuplicationRate | agreed | v6 `DuplicateFraction` (Picard's definition, pooled over libraries); not the old `DuplicationReads`, which counted the `.duplication.bam` split (structural-duplication evidence) |
| Species | agreed | new: a copy of `SampleDatabase.Species` at alignment time; the audit flags a mismatch, a refresh command updates it |
| MappedFraction, QC, PipelineVersion | proposed, welcomed | v6 has `MappedFraction`, `QCWarnings` and `PipelineVersion` (the text `v2`; to become the git commit) |
| Usable coverage | proposed | v6 measures `CoverageNonDuplicate` (and `CoverageMapped`) instead of a derived formula |
| NRuns, InputFingerprint, AlignedOn, Machine, Instrument, ReadLength | proposed | new (the fingerprint is a hash of the sorted (RunID, uBAM size) pairs) |
| about 20 other v6 statistics (`PrimaryReads`, `UnmappedFraction`, `UnmappedNoPositionFraction`, `UnmappedOrphanFraction`, `DuplicateFractionFromBam`, `OpticalDuplicateFraction`, `MeanReadLength`, `MappedBases`, `ProperlyPairedFraction`, `MQ0Fraction`, `InsertSizeMean`, `InsertSizeSD`, `ErrorRate`, `UBamReads`) | **(open)** | proposal: keep them in the per-sample `qc.json` beside the BAM and keep the sheet lean; `QCWarnings` stays in the sheet |

   Dropping the six fractions means changing v6's required legacy columns (`LEGACY_KEYS`) and the stage that counts the splits. Whether the split BAMs themselves are still produced is **(open)**; the old structural-variant genotyping reads them. Panel and genotype results do not belong in the ledger: a changed panel definition must never trigger a realignment.
3. **The rerun rule (new in either pipeline).** A sample needs (re)alignment if it has no row; or its run set now differs from the row's (runs added or removed; the old RunIDs separator is a double comma and is parsed as a set); or its InputFingerprint differs (a uBAM re-created under the same RunID, for example re-downloaded from SRA instead of ENA); `--force` overrides. A run that is in the database but whose uBAM is not yet on Dropbox makes the sample wait and be reported, never be aligned from a subset. A rerun replaces: upload the new set, verify, update the row in place (v6 only appends). A change of tool versions or reference is recorded but triggers nothing unless asked. The panel JSON carries the alignment fingerprint so stale genotype results are detectable.
4. **Ledger writes.** v6's append with read-back is the base. Still to add: update in place for reruns; a local journal replayed when Google is unreachable (v6 raises after four attempts and leaves an uploaded sample without a row, which the audit and `register_alignments` then repair); claims so two machines do not align the same sample (nothing in v6 prevents it). **(open: ledger only, or also a queue with status and claims, as the DownloadDatabase has.)**
5. **Supervision and reported skips.** Every subprocess has a timeout or stall limit and dies with its parent; selection lists the samples it leaves out and why (no SampleDatabase row, blank Ecogroup, missing uBAM). **(unverified: how far v6 already does this; `--list` previews the selection.)**

**Work.** (1) Validate v6 on real data: the BAMs of the yellowhead run are the reference (`validate_against_old.py`: the seven BAMs must be identical); keep them until then. (2) Settle the ledger columns, then change `LEGACY_KEYS`, the header setup and the writer. (3) The rerun rule with tests. (4) Update in place, the journal, claims. (5) The panel stage (04) after duplicate marking, so each BAM is read once. (6) Check supervision and the reporting of skipped samples. (7) Retire the old script after the mbuna run.

**Done when.** One sample gives BAMs identical to the old pipeline's; a batch of mbuna samples runs on two machines without intervention; a killed run resumes without redoing finished samples; the ledger row of every sample states its runs, its fingerprint and its pipeline version; a sample with a new run is selected for rerun and one without is not; a Google outage delays the ledger but loses nothing; tests pass offline.

### 04_genotype_sd_systems: gsdf duplication, transposon insertions, inversions, B chromosome

**Goal.** Per-sample genotypes of the SD systems, using user-provided definitions, normally produced by 03 while the BAM is local, and re-runnable on samples that missed it.

**Inputs.** A local indexed BAM (or a stored slice); definitions: the gsdf region and calibration, the transposon insertion definitions for gsdf and chr10 (**provided by the user**), inversion breakpoint windows (01), B-chromosome core blocks (its own discovery step).

**Outputs.** `<sample>.panel.json`; call tables; columns in `SampleDatabase` (for example `gsdf number`, `TE1`, `TE2`, `Inversion*`, B status) with a source and flag column for each.

**Interface.** `sd_genotypes panel measure --bam B --sample S --out-json J [--modules ...]`; `sd_genotypes call <module> --measurements DIR`; `sd_genotypes validate <module>`; `sd_genotypes bchrom prepare|call|compare`; `sd_genotypes sheet-update`.

**Current state.**
- *gsdf duplication* (`panel` and `gsdf`, tests with fakes): raw 5 kb bins and junction reads stored, copy number called with the calibration of the MC pedigree (0.530 and 0.268), flags for low depth, non-integer, range, junction disagreement; validation against the existing 910 `gsdf number` values and sex-by-species tables. **(unverified: real BAMs; the calibrated range of 2 to 4 copies is an interpretation.)**
- *B chromosome* (`bchrom`): downloads Clark et al.'s public reads (through the SRA source of module 02, `download.source: auto`; all 21 runs exist on NCBI's AWS copy, 280 GB of `.sra`; reorganised tree), aligns them with the lab recipe, finds blocks and a core set in reference coordinates, evaluates each carrier leave-one-out, separates tier A (normal male-pool coverage) from tier B (low baseline). Not yet run on real data; the work continues in a separate session on D07. The per-sample genotype step (a score over the core blocks, a threshold from the data) is not written.
- *Inversions*: Blumer et al. (Table S7) give chr2, 9, 10, 11, 13 for 1,375 samples, 1,364 of them now in the SampleDatabase (section 7); all 440 mbuna among the first 1,101 are 0/0 on all five. These are **validation labels**, not an input. The LG5 (orange blotch / ZW) inversion is polymorphic in the mbuna and has no labels; PRJNA808791 (Behrens et al.) has BB and OB pools of known composition; `bchrom prepare` keeps slices around both LG5 breakpoints for them. A breakpoint-read genotyper is designed but not written.
- *Transposon insertions*: previously produced by a separate tool (`sv-genotype`, columns `TE1`, `TE2`); to be re-implemented here from user-provided definitions. **(open: the definitions, the evidence to use, and which samples.)**

**Work.** Finish the B-chromosome run and the genotype score; write the inversion genotyper (breakpoint-read evidence with windows, validated on Blumer labels, BC families and the pools); implement the transposon genotyper; validate gsdf on about 40 samples with known values; one `panel` JSON schema for all modules; the sheet write-back; the retrospective mode (BAM or slice in, JSON out).

**Done when.** gsdf: unflagged calls agree with the existing column for at least 98% of the samples and females are 2 and males 3 in the MC family; B chromosome: carriers score above the threshold leave-one-out, controls and mbuna males below it; LG5: BB pools about 0 and OB pools about 0.5 inverted; inversions: at least 98% agreement with Blumer on labelled samples, with hets reported separately; transposons: agreement with the old `TE1`/`TE2` columns where they exist; every measurement carries its definition hash.

### 05_variant_calling: VCFs for families and species

**Goal.** Reproducible VCFs for the cohorts the SD analyses need: broods and families with parents, and groups of related species where the sex system is unknown.

**Inputs.** BAMs and slices (03), cohort definitions from the sample table (species, brood, parents, sex, ecogroup), a regions list, the reference.

**Outputs.** One bgzipped, tabix-indexed VCF per cohort with a manifest (samples, caller, version, filters, regions, sites and their attrition).

**Interface.** `variant_calling cohort build <name>`, `cohort list`, `qc <name>`.

**Current state.** `callset build` in `YH_SD_Paper` (callsets `bc_family`, `yh_broods`, `mc_broods`): per-sample GVCFs (cloud or HaplotypeCaller from the BAM through `ensure_gvcf`), `GenomicsDBImport`, `GenotypeGVCFs`. HaplotypeCaller is no longer part of alignment (module 03). **(open: the caller, joint versus per-family calling, region-restricted versus whole genome, the filter set, and how unknown species are grouped.)**

**Work.** The cohort file format; a caller benchmarked on one brood against the existing VCFs; a cohort builder that needs no manual sample lists; QC (missingness, depth, Ti/Tv, Mendelian errors); the manifest.

**Done when.** The BC and brood cohorts rebuilt from the sheet reproduce the existing VCFs' calls at the shared sites (stated concordance); every cohort VCF has a manifest; Mendelian-error rates are reported.

### 06_sd_region_discovery: SD regions in many species

**Goal.** For each species or cohort, the regions associated with sex, with the statistics and the evidence that identifies the kind of system (XY or ZW, inversion, dosage).

**Inputs.** VCFs (05), the sample table, 04's calls, 01's landmarks.

**Outputs.** Candidate regions per cohort (BED and statistics: association, Fst, heterozygosity by sex, male/female coverage ratio, parent-offspring inheritance) with plots, and a species-by-region summary.

**Interface.** `sd_regions scan <cohort>`, `sd_regions summarize`.

**Current state.** `parental_sex_region` and `sex_association fisher` exist and were tested on synthetic data; not run on real data. **(open: species list, minimum sample sizes, multiple-testing and permutation scheme.)** The SampleDatabase now carries Sex for the Blumer samples (section 7).

**Done when.** It recovers the known systems (the chr10 Y inversion in the YH lineage, the gsdf dosage pattern in *M. conophoros*, the LG5 locus in the orange-blotch species) as positive controls and reports nothing in species known to have none.

### 07_bc_analysis: backcross genotyping, recombination, transmission distortion, QTL mapping, fine mapping

**Goal.** Genotype the BC individuals from their reads, then analyse them: recombination, transmission distortion, QTL mapping, and fine mapping.

**Inputs.** BAMs (03), cohort VCF and sites (05), pedigree (sample table), 01's landmarks, traits.

**Outputs.** Genotype and posterior matrices, breakpoints, per-fish quality, recombination maps (including suppression in inversions), TD tests, QTL scans, fine-mapping tables.

**Interface.** `bc_analysis genotype <step>`, `recombination`, `td`, `qtl`, `finemap`.

**Current state.** `bc_genotyping` (about 266 fish, 9,336 markers; steps `provenance`, `flagstat`, `original`, `v2-sites`, `v2-counts`, `v2-hmm`, `validate`, `compare`) is the part to **fix**; the analyses after it are not written. **(open: QTL software.)**

**Done when.** The genotyper agrees with the original at the shared markers, recovers the expected F1 phase and the known inversion state of each fish, and shows suppressed recombination across the chr10 inversion; each analysis reproduces from one command.

### 08_ase_variants: variants for allele-specific expression

**Goal.** For a gene in a family, variants that tell the two alleles apart in RNA-seq, with the information needed to trust them. **Inputs.** VCFs (05), annotation, mappability, parent and F1 genotypes. **Outputs.** Per-gene, per-family tables (position, alleles, parental origin, transcript context, coverage, mappability, reference-bias risk, rank). **Interface.** `ase_variants suggest --gene G --family F`, `batch --genes FILE`. **Current state.** Not started **(open: annotation source, minimum coverage, how reference bias is judged)**. **Done when.** For ten test genes it returns the same variants as a hand-made list, flags variants in repeats or paralogs, and says when a gene has none.

### 09_candidate_qtn_webserver: candidate causal variants and the web server

**Goal.** From the mapped intervals, the genetic differences between parental assemblies, annotated, with their allele and genotype frequencies in other species with and without the SD system, presented in a web server. **Inputs.** Assemblies and their alignments, QTL and fine-mapping results (07), 04's genotypes, cohort VCFs (05), SD regions (06), annotation. **Outputs.** Candidate-QTN tables, frequency tables, correlation with SD, a read-only server. **Interface.** `candidate_qtn diff --interval I`, `frequencies`, `build-site`, `serve`. **Current state.** Not started **(open: hosting and audience; a static site or an API over a database; an embedded genome browser)**. **Done when.** For one known interval the page lists every difference between the assemblies, shows its frequency by species and links to the evidence; rebuilding the site from the run folders is one command.

### 10_paper_figures: figures and data tables

**Goal.** Every figure and supplementary table regenerated by one command, each naming the run folder it was built from. **Inputs.** Run folders of 01 to 09. **Outputs.** SVG and PDF figures, TSV and XLSX tables, `figure_sources.tsv`. **Interface.** `paper_figures build <figure>`, `build --all`, `sources`. **Current state.** The to-scale inversion figure (six chromosomes) and the introduction figure exist; figures 2 to 4 (transposon insertions) come after 04 **(open: the figure list)**. **Done when.** `build --all` on a clean machine reproduces every figure and table from Dropbox run folders, and a changed upstream run is detected as a stale figure.

---

## 7. The SampleDatabase

**Rule for every tool that edits it:** a report first, a backup of the tab, one write of the tab only, a cell-by-cell read-back. Blanks are filled; a value that differs is kept and listed, because the sheet carries manual curation (`Error`, `Exclude`, a sample labelled Utaka although Blumer's clade is Shallow). The one deliberate overwrite is `Species` for samples in Blumer's table.

**The Blumer update** (`unit_scripts/sampledb.py blumer`; tables in `data/blumer/`; written 2026-11-07, verified cell by cell). For samples whose SampleID is a Blumer BioSample (1,364 of Blumer's 1,375 matched):
* `Species` = Blumer's `genus species` (388 changed, 129 distinct old-to-new pairs). The SRA name is not lost: `SRASpecies` records the organism of the sample's runs in the DownloadDatabase, for every sample that has runs there (1,492 samples). The sheet's earlier names were SRA names (the downloader writes the organism into `Species` of every new sample); 29% of the 1,101 samples checked first differ from Blumer's, for example *Cynotilapia afra* (SRA) against *Cynotilapia zebroides* (Blumer) ×66 and *Maylandia zebra* against *Maylandia callainos* ×26. Informal taxa are spelled differently (`Tropheops sp. 'yellow gular'` against `Tropheops yellow-gular`); samples outside Blumer's table keep their SRA names.
* `Sex`, `Ecogroup` filled if blank. Blumer's clades map one to one onto Ecogroups: Mbuna, Shallow → Shallow_Benthic, Deep → Deep_Benthic, Utaka, Diplo → Diplotaxodon, Rhampho → Rhamphochromis, AstCal → AC. On the 1,101 samples checked first, sex agreed for all and Ecogroup differed for 4 (`Error` ×2, `Exclude`, one `Utaka`), kept.
* `LabReared` 0 (every Blumer sample is wild-caught; all 1,101 checked had 0).
* `BlumerSequenceID` (= S1 `sequence_id` = S7 `sample_id`), `BlumerSpecimenID`, `Inversion2/9/10/11/13` from Table S7 (codes 0, 1, 2; Inversion10 already agreed for all 1,101), and `InversionSource` with the value the existing Blumer rows carry (`Blumer2025_TableS7`).

**State after the update (to be confirmed by `sampledb.py check`):** about 71 samples with a blank Ecogroup and LabReared and about 246 with a blank Sex remain; they are non-Blumer samples (the Malinsky data and others) for which no metadata table is in the repository; 176 samples already had a blank Sex in the earlier export. `Subgroup` and `Category`, which the lab curates and some analyses group by, are blank for the new samples. 11 samples of Blumer's table are not in the sheet (10 BioSamples and one row without a BioSample).

**186 samples have reads but no sample row (open).** 184 carry the project label `LakeMalawiDurbin` (150 of them with the SRA organism `Astatotilapia` and no species; the rest riverine haplochromines such as *Haplochromis tweddlei*, *Serranochromis robustus*, *Astatotilapia bloyeti*), 2 are in `PRJEB21343` but not in Blumer's table. Without a sample row a sample cannot be selected for alignment. **(open: create rows for them with the SRA species and a new Ecogroup value such as Outgroup, create rows for some, or leave them out and document it.)** `sampledb.py check` writes the list (`sampledb_orphans.tsv`).

**Checks** (`sampledb.py check`, read-only): duplicate SampleIDs and RunIDs; orphans in both directions; DownloadDatabase runs missing from `DNAReadsDatabase` or whose SampleID is neither the BioSample nor the library name; codings of Sex and LabReared; blank key fields; taxa written in more than one way; comparison with Blumer.

## 8. Machines and environments

| Machine | Role | Notes |
|---|---|---|
| utaka | alignment and heavy work (96 cores) | branch `mcgrath_dev`, Python 3.7 env `CichlidSRSequencing`; the old file manager (no download-database methods, no retry); the audit of 2026-11-07 ran here; the 178 Hybrid rows are still to register |
| mzebra | downloads, conversions, sheet maintenance (96 cores, 503 GB RAM) | login `bio-mcgrath-dropbox`; environment `enadl`; the download campaign ran here and the checks of this repository run here |
| D07 (workstation) | the B-chromosome work, development (32 cores, 62 GB RAM) | Aspera blocked; SRA about 11 MB/s with 8 connections (2.7 MB/s with one); needs minimap2, mosdepth, scipy, matplotlib in its environment |
| yellowhead | validation of the pipeline (student machine) | on 2026-11-07 a batch of about 45 samples was running with the old script and the v3/v4 download queue, in the GVCF stage; keep its local BAMs until v6 is validated |
| laptop | development, tests, publishing | tests run offline; the `syri_env` environment has an installed package called `tests` that shadows the repository's, so run `python3 -m unittest discover -s tests`, not `tests.<module>`; `psutil` and `pysam` are missing there and on mzebra, so three older tests of the old alignment code cannot import |

Each module has an `env.yml`. A known conflict: `aria2` and `sra-tools` cannot be solved in one conda environment on macOS; on Linux they can.

## 9. Order of work

1. **Done:** module 02 for both Run Info files; the Blumer update of the SampleDatabase; the tools of section 3.
2. **Now:** the post-update `sampledb.py check`; the decision about the 186 samples; the recheck of the 4 unverified runs; on utaka the registration of the 178 Hybrid rows (`register_alignments.py --check-existing`, then a three-sample real write, then the rest); the outcome of the yellowhead batch.
3. **Next (module 03):** validate v6 against the yellowhead BAMs; the ledger columns and upsert writer; the rerun rule; the panel stage. This gates the mbuna alignment run (about 477 samples).
4. **In parallel (module 04):** the B chromosome on D07; then the LG5 genotyper, the transposon genotyper, gsdf validation.
5. **Then:** 05 and 06, 08, 07; 09 and 10 read everything.

Critical path to the mbuna run: 03 validated and redesigned, the panel defined (so every BAM is measured once), the Ecogroup of every sample set.

## 10. Migration from today's layout

| Today | Becomes |
|---|---|
| `YH_SD_Paper/yhsd/analyses/centromeres`, `inversion_breakpoints` | 01 |
| `helper_modules/ena_fetch.py`, `sra_fetch.py`, `ena_scheduler.py`, `grab_ena_v2.py`, `download_manifest.py`; `downloadShortReadData_v2.py`, `downloadManifest.py`, `unit_scripts/ena_bench.py`, `patch_file_manager.py`, `verify_ubam.py` | 02 (and `common/` for the sheet backend) |
| `helper_modules/gspread_retry.py` | `common/` |
| `unit_scripts/audit_alignments.py`, `register_alignments.py` | 03 |
| `unit_scripts/sampledb.py`, `data/blumer/` | 04 (`sample_table`) |
| `alignFastQReads_v2.py`, `seq_pipeline.py`, `bam_split.py`, `bam_slices.py`, `bam_stats.py`, `alignment_db.py`, `disk_ledger.py`, `download_manager.py` | 03 (`disk_ledger` also in `common/`) |
| `YH_SD_Paper/yhsd/analyses/panel`, `gsdf`, `bchrom`, `sample_table` | 04 |
| `YH_SD_Paper` callset builder (`callset build`) | 05 |
| `parental_sex_region`, `sex_association` | 06 |
| `bc_genotyping` | 07 |
| `YH_SD_Paper/figures` | 10 |
| `YH_SD_Paper/yhsd/runner.py`, `common.py`, `genomes_registry.py`, `config/genomes.yaml` | `common/` |
| `downloadShortReadData.py`, `grabENA.py`, `alignFastQReads.py`, `alignment_worker.py`, `download_queue.py` | `legacy/` once retired |

**Rules for the move.** It happens in one change, after the pilots and while no download or alignment is running on any machine; the downloader launches its workers by file path, and the alignment run imports the old modules. Thin wrapper scripts under today's names stay for one release. The test suites are the safety net: they import by name and fail loudly if anything was missed. The package name `yhsd` is retired in favour of the module packages above.

**Status (2026-11-07).** The reorganised tree (module folders, a shared `cichlid_common` with the data manifest, one `cli.py` per module, tests per module, a migration script) is built and passes its tests, but it is **not applied**: it is meant to live on a branch (`reorg_spec`) until no download or alignment job runs anywhere. Module 03 is not moved with the others, because the alignment runs on utaka and yellowhead import the old modules and the old file manager. The original layout keeps working through the table in section 3; this document describes both.

## 11. Decisions and measurements

1. **Reads come from NCBI's copy on AWS when ENA is slow (2026-11-06).** ENA: about 1 Mbit/s over Aspera and 0.3 MB/s over https from utaka, D07 and a home connection, four files, three protocols; NCBI's copy 17 to 40 MB/s per connection from the same places. Not a firewall problem; no request to IT.
2. **Every SRA run is checked against ENA** (integrity, counts, first reads). ENA's `read_count` has no single meaning (ERR13389718: 109,512,218 = 2 × 54,756,109 spots; ERR14734229: 40,613,306 against 40,650,256 spots). SRA has 0.01% to 1.71% more spots than ENA, never fewer; accepted up to 5%, rejected if more than 0.1% fewer. A hypothesis, not verified: ENA's FASTQ leaves out some reads that SRA keeps.
3. **Binned qualities are not a sign of simplified-quality SRA files.** ERR13389718 (HiSeq 2500) has 6 distinct quality symbols at ENA and at NCBI. Only one or two symbols is treated as simplified.
4. **The DownloadDatabase is a spreadsheet of its own**, one tab per Run Info file, written cell by cell with claims and a journal; an unreachable sheet never stops a download.
5. **Database rows are added at the end, in bulk** (2026-11-07): `--no-db` then `--add-pending --bulk`. The first run-by-run attempt stopped at Google's quota after 26 rows; the bulk mode makes identical sheets (tested) with 2 reads and 1 write.
6. **A bad SRA run exists and the check caught it:** ERR14734399; excluded.
7. **Speeds:** mzebra 104 MB/s for one SRA download and about 150 MB/s for three; D07 about 11 MB/s; utaka about 40 MB/s (one connection). Conversion 19 to 23 minutes per run; Dropbox upload 20 to 25 MB/s per stream with about 4 streams in parallel.
8. **No GVCF in alignment** (2026-11-07): v6 has none; downstream code builds them on demand (`ensure_gvcf`).
9. **The alignment ledger** keeps SampleID, GenomeVersion, RunIDs, Coverage, TotalReads, versions and BamSize; adds DuplicationRate and Species; drops the six split-read fractions (agreed in discussion). Reconciliation with the ~20 statistics columns that v6 already writes is open (section 6, module 03). `DuplicationReads` (the `.duplication.bam` split: structural-duplication evidence) and DuplicationRate (PCR and optical duplicates) mean different things.
10. **Rerun when the runs change**, by run set and by fingerprint (neither pipeline does it today).
11. **Blumer is the correct species; the SRA name is kept in `SRASpecies`** (taken from the DownloadDatabase for every sample that has runs there). Blanks are filled, differing values kept and listed, the one overwrite is `Species`.
12. **Blumer's tables live in the repository** (`data/blumer/`) so the update runs from a checkout with no files to upload.
13. **A convenience layer must not be able to break the program** (the Google retry cannot raise at import) and **a hung process is a bug** (the 88-day `rclone`).
14. **The B chromosome is found in the project reference's own coordinates** (blocks and a core set found there, each carrier evaluated leave-one-out), following Clark et al.'s coverage-ratio idea, instead of lifting his blocks from the older assembly; its reads come through the SRA source.
15. **Module 03 is not moved with the others** while an alignment run uses the old file manager.

## 12. Open decisions

1. LG20: what is known, and what it adds to the paper. (01)
2. The transposon insertion definitions and the evidence to use. (04)
3. Variant caller, joint or per-family calling, filters, region scope. (05)
4. Species list and sample-size rules for SD discovery. (06)
5. QTL software. (07)
6. Annotation and mappability sources for ASE. (08)
7. Web server hosting, audience and data store. (09)
8. The figure and table list. (10)
9. When to apply the reorganised layout, and who updates each machine.
10. Whether `SPEC.md` is edited directly in the repo, or changes go through a decision log in `docs/`.
11. GVCFs for the 45 yellowhead samples: let the old run finish its GVCF stage, or stop it and finish by hand (upload, then `register_alignments`). (03)
12. The ledger columns: reconcile the agreed set with v6's ~20 statistics columns (sheet or `qc.json`). (03)
13. Whether the split BAMs stay. (03)
14. Whether runs can be excluded individually, or only whole samples. (03)
15. The ledger as a result sheet only, or also a queue with claims. (03)
16. How a stale `Species` in the ledger is refreshed. (03)
17. The 186 samples without a sample row, and an Ecogroup name for outgroups. (SampleDatabase)
18. A metadata source for the remaining samples without Ecogroup, Sex and LabReared (the Malinsky supplement); `Subgroup` and `Category` for new samples. (SampleDatabase)
19. ERR14734399: stay excluded or fetch from ENA. (02)

## Appendix A. For the methods

Reads of the samples of Blumer et al. and Malinsky et al. that were not already in the lab's database (334 runs) were retrieved from NCBI's copy of the Sequence Read Archive (public AWS bucket) because transfers from ENA were about 0.3 MB/s. Each run was validated (`vdb-validate`), converted with `fasterq-dump`, and compared with the ENA record: the number of read pairs agreed with ENA's `read_count` within 5% (SRA held 0.01% to 1.71% more spots than ENA in every run), and the sequence and base qualities of the first 1,000 reads were identical to ENA's FASTQ in 323 of 327 runs checked (the remaining 4 were not comparable at the time and are to be rechecked). One run, ERR14734399, whose first reads did not match ENA's at all, was excluded. Species and sex were taken from Blumer et al.'s Table S1 and inversion genotypes from their Table S7 for the samples they describe; the species names recorded in SRA are kept separately.
(Wording to be confirmed once the 4 unverified runs are rechecked and the final run counts are fixed.)
