# CichlidSRSequencing: specification

Status: v5 (2026-11-09). Owner: McGrath lab. This file lives at the repository root and is the source of truth for what each part of the repository is for, what it consumes, what it produces, and when it is done. It extends v4 (2026-11-08), which merged the downloads, the alignment redesign and the SampleDatabase of v3 with the genotyping version. The log of decisions is `docs/decisions.md` (section 11).
Where a statement is not yet verified it is marked **(unverified)**; where a decision is still to be made it is marked **(open)**.

## 0. What changed since v4

* **Module 03 is built and validated on real data.** On one sample the old pipeline had aligned, the new pipeline gave `all.bam` and the chimeric BAM identical to the old ones, equal statistics and a current panel JSON, in 0.5 hours. Its first real batch (the 8 unaligned Deep_Benthic samples) is running on yellowhead; its upload, ledger write and cleanup are not yet exercised (section 6, module 03; decisions 46 to 52).
* **Decisions of 2026-11-09 for module 03:** no quality thresholds in the ledger; reruns realign everything together, replace the BAMs and the row, and delete the stale GVCF and old splits; disjoint selections instead of claims; a pinned environment; three known defects to fix before the next batch.
* **The PCA basis panel is proposed** (module 04): a second method for the inversion states of LG2, 9, 10, 11, 13 and 20, and the check of a sample's ecogroup, built from the large polymorphism file of the Malinsky species. Nothing is built; the first step is to describe the file and run a pilot on one inversion (decision 54).
* **Three samples flagged `Error` or `Exclude`** were looked at (section 7): their species agree in SRA, Blumer and the sheet, and they cluster with their species in Blumer's tree (decision 53).
* **The decision log has entries 46 to 54.**

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

### Where the code is today

Three branches of the one GitHub repository (`ptmcgrat/CichlidSRSequencing`) hold the code:

| Branch | Layout | Holds |
|---|---|---|
| `YH_SD_Paper` | original (`cichlid_sr_sequencing/`, `YH_SD_Paper/`) | the downloader and its tools, the SampleDatabase tools, the Google retry, the alignment audit and registration (pushed 2026-11-07) |
| `mcgrath_dev` | original, older (last commit 2026-09-21) | the old alignment pipeline as run on utaka |
| `reorg_spec` | the module folders above | module 04 (`sd_genotypes`, the panel) and the work on it; `docs/decisions.md` |

| Piece | Location |
|---|---|
| Downloader (ENA and SRA sources, `--manifest`, `--add-pending [--bulk]`, `--max-runs`) | `YH_SD_Paper`: `cichlid_sr_sequencing/downloadShortReadData_v2.py` with `helper_modules/ena_fetch.py`, `sra_fetch.py`, `ena_scheduler.py`, `grab_ena_v2.py` |
| DownloadDatabase tool (import, status, reconcile, release) | `downloadManifest.py`, `helper_modules/download_manifest.py` |
| Retry of every Google call | `helper_modules/gspread_retry.py`, installed by `helper_modules/file_manager.py` at import |
| File manager patcher (download-database methods; `--retry-only`) | `unit_scripts/patch_file_manager.py` |
| uBAM recheck against ENA | `unit_scripts/verify_ubam.py` |
| Alignment audit (read-only) and registration of uploaded-but-unregistered alignments | `unit_scripts/audit_alignments.py`, `unit_scripts/register_alignments.py` |
| SampleDatabase checks and the Blumer update | `unit_scripts/sampledb.py`, input tables in `data/blumer/` |
| Network benchmark; tests (offline); documentation of these tools | `unit_scripts/ena_bench.py`; `tests/`; `docs/DOWNLOAD_AND_DATABASE_TOOLS.md` |
| Alignment pipeline, old: `alignFastQReads.py`, `helper_modules/alignment_worker.py`, with the download queue of patches v3/v4 (`helper_modules/download_queue.py`) | `./`, `helper_modules/` |
| Alignment pipeline, v6 (the base of `align_pipeline`; offline tests): `alignFastQReads_v2.py` with `seq_pipeline.py`, `bam_split.py`, `bam_slices.py`, `bam_stats.py`, `alignment_db.py`, `disk_ledger.py`, `download_manager.py`, `slice_regions/`, `unit_scripts/validate_against_old.py`; described in `PIPELINE_V2_README.md` | `./`, `helper_modules/` |
| The panel (gsdf, te, lg5inv, bchrom), retrospective measurement, parity with sexloci | `reorg_spec`: `04_genotype_sd_systems/sd_genotypes/`, `docs/USAGE.md` |
| Analyses of the paper (centromeres, inversion breakpoints, sample table, BC genotyping, callsets, figures) | `YH_SD_Paper/` (package `yhsd`) |

The tools pushed to `YH_SD_Paper` on 2026-11-07 (`sampledb`, `register_alignments`, `audit_alignments`, the hardened `gspread_retry`, `data/blumer/`) are also in the reorganised tree prepared the same day; whether `reorg_spec` already holds them has to be checked **(open)**. The utaka checkout runs `mcgrath_dev`, with the old alignment pipeline and an older `file_manager.py` (410 lines, no download-database methods, no retry).

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
| 03 | `Bamfiles/<genome>/<sample>/`: `<sample>.all.bam` + `<sample>.all.bai` (the names are part of the contract), `<sample>.chimeric.bam` + index, slices of LG2, LG10 and LG13, `qc.json`, `<sample>.panel.json`; one row per sample in `AlignmentDatabase` | 04 (retrospective measurement streams `all.bam` by name), 05, 07 |
| 04 | `<sample>.panel.json` (raw measurements, with the definition key of every module); call tables (`gsdf_calls.tsv`, `te_calls.tsv`, `inversion_calls.tsv`, `bchrom_calls.tsv`); calls written back to `SampleDatabase` columns | 05 (cohort definitions), 06, 07, 09, 10 |
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

**Constraint (resolved 2026-11-08).** The alignment run on utaka used an older `file_manager.py` that can overwrite another machine's sheet changes, so `common/` first re-exported the existing file manager. That run has ended (decision 43); `common/` still re-exports whatever file manager the machine has until every machine is updated, and the Google retry is installed by the file manager at import.

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


**Proposed finish plan (2026-11-08, awaiting confirmation).**

*What the literature and the published tables say.* Kumar et al. report six rearrangements: single inversions on chromosomes 2, 9, 10 and 13 and **double (tandem) inversions on 11 and 20**; on 20 the second part (20b, 4.1 Mb) is absent in four *Rhamphochromis*, which carry only the first (20a). Blumer et al. also describe two adjacent inversions on chromosome 20 (fig. S16, text S1), but their Table S7 covers only chromosomes 2, 9, 10, 11 and 13, so **LG20 has no labelled samples**. The published centromere table (Chromosome Research 2026, Table 3) annotates all 22 chromosomes of eight Lake Malawi assemblies (fMayZeb221, fLabFul222, fLabTre303, fLetChi2, fAulStu2, fDipLim2, fRhaChi3, fAstCal68; GenBank accessions in its Table 1) and *A. burtoni*; for *M. zebra* it gives LG20 at 19.58 to 20.36 Mb (not T2T) and LG13 at 1.90 to 3.78 Mb (1.88 Mb wide, ungapped, not T2T). The *Genome Biology and Evolution* paper in the project notes that the orientation of LG13, LG14 and LG20 differs between *M. zebra* and *A. calliptera*.

*Proposals.*
1. **LG20 is in scope as two rows, 20a and 20b,** sharing a system id; LG11's double inversion is recorded the same way, so tandem inversions have one representation.
2. **`reference_state` and `derived_state` columns** in `inversions.tsv`, polarity from Kumar's outgroup analysis (*O. niloticus*, *P. nyererei*). The mbuna are 0/0 in Table S7 for 2, 9, 10, 11 and 13, which suggests GT3 carries the non-inverted state there; LG5 (`reference_is_standard` unset) is settled from the assemblies.
3. **LG9 is declared a misassembly only if the mbuna assemblies** (fMayZeb221, fLabFul222, fLabTre303) agree with each other but not with GT3; a true polymorphism would not show in species that share the reference's arrangement.
4. **LG13 centromere:** report the array scan of this module and the published interval, with the published width as the uncertainty; no further work.
5. **The published Malawi assemblies are the additional assemblies:** fAulStu2 for LG2, 10 and 13, fDipLim2 and fRhaChi3 for LG20, expected from Kumar's species distributions (inversions 2, 10, 13 only in deep benthic *Aulonocara*; 9, 11, 20 mainly in benthics and *Diplotaxodon*) and to be confirmed by this module's own alignments. Kumar's own *Aulonocara* assemblies have accession "TBD" in the text at hand and are not relied on.
6. **Schemas and validators** for `centromeres.tsv`, `inversions.tsv`, `inversion_breakpoint_windows.bed` and `known_misassemblies.tsv`; the windows are each breakpoint's uncertainty interval plus a fixed pad, chosen once the uncertainties are known. 03 and 04 read the windows (04's goal 6 is deferred until they exist).

*Not part of 01:* sample genotypes (04), read-level refinement of breakpoints, figures (10), LG20 sample labels.

*Order of work:* register the published assemblies and align them on utaka (the long step); re-run breakpoint calling for LG2, 5, 9, 10, 11, 13 and 20; lift the published centromeres and compare with the array scans; apply the rules of 3 and 4; write schemas, validators and tests; produce the final run folder; update this section.

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

### 03_align_pipeline: align, mark duplicates, split, slice, QC, measure, upload

**Goal.** Turn uBAMs into aligned, duplicate-marked BAMs for a registered genome, with the one split and the slices later modules need and the per-sample measurements of 04 made in one pass while the BAM is on local disk, one sample at a time, with a ledger (the AlignmentDatabase) that says what exists, how it was made and from which reads.

**Inputs.** uBAMs (02), the reference and its minimap2 index, the slice definitions (LG2, LG10, LG13; contig names from the genome registry, later from 01's tables), the panel definitions (04; `panel check`), the SampleDatabase (selection and metadata).

**Outputs.** In `Bamfiles/<genome>/<sample>/`: `<sample>.all.bam` and `<sample>.all.bai` (**these names are part of the contract**: 04's retrospective measurement streams the BAM from Dropbox by name), `<sample>.chimeric.bam` with its index, the slices of LG2, LG10 and LG13 (BAM and index), the duplication metrics, `qc.json`, `<sample>.panel.json`; one ledger row; no uBAM kept after use. The five other splits (unmapped, discordant, inversion, duplication, clipped) are no longer produced.

**Where it is built.** `03_align_pipeline/` of the reorganised layout on branch `reorg_spec`, in the same checkout as `04_genotype_sd_systems/sd_genotypes` (decision 43). The code of today (the old script and v6) is in the original layout (section 3).

**Interface.** `python -m align_pipeline run <genome> [-s -c -p -e -u] [--include-excluded] [--batch-size N] [--samples-in-flight N] [--keep-local] [--no-upload] [--no-slices] [--no-panel] [--panel-expect KEY...] [--disk-reserve-factor F] [--list] [--print-commands SAMPLE]` and `python -m align_pipeline setup-database <genome>`, after `source tools/env.sh` (it puts every module on `PYTHONPATH`; the panel subprocess inherits it). `-r` selects samples that already have a ledger row (v6's meaning; the rerun rule below replaces it). The old scripts stay on `YH_SD_Paper` and `mcgrath_dev`.

#### What exists

* **The old pipeline** (`alignFastQReads.py`, `alignment_worker.py`; utaka and yellowhead, the latter with the download queue of patches v3/v4): every stage is a barrier over all samples, alignment is one sample at a time with `minimap2 -t 96`, HaplotypeCaller GVCFs are made for every sample, and the upload and the database rows come last, one row at a time through `FileManager.addAlignmentRow`.
* **`align_pipeline`** (this module, `reorg_spec`; v6 ported, with the decisions of 2026-11-08 and 2026-11-09 applied; **validated on one real sample on 2026-11-09 and running its first batch**). Each sample moves on its own (download, align, mark duplicates, chimeric split / statistics / slices / panel, qc, verified upload, ledger row, cleanup); several samples in flight with a shared core and RAM budget; space-aware continuous downloads with a disk ledger; a separate upload pool; stage checkpoints so a rerun resumes; no GVCFs; `samtools stats` and `idxstats` instead of seven `gatk CountReads`; **one row appended atomically with read-back and retry** (`AppendingDatabaseWriter`); a panel preflight before anything is aligned; time limits on transfers, the quick commands and the panel; status lines every 10 minutes; samples with Ecogroup `Error` or `Exclude` are skipped unless `--include-excluded`; the commit of the checkout is read from the `.git` folder when there is no git program. 140 offline tests (158 where the audit and registration tools are present: `reorg_spec` does not hold them), with fake tools and a few real subprocesses running stand-in programs. Deployment and validation steps: `03_align_pipeline/DEPLOY_AND_VALIDATE.md`.
* **Tools** (`unit_scripts/` on `YH_SD_Paper`): `audit_alignments.py` (read-only: the ledger against Dropbox against the local machine; classes ok, needs_row, partial_upload, needs_upload, local_incomplete, row_only, not_started) and `register_alignments.py` (adds rows for uploaded alignments: GATK counts as the old loop did, in parallel, one read and one write of the sheet, backup, read-back; `--check-existing` compares every computed column with existing rows).

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

#### Validated on real data (2026-11-09, yellowhead)

One sample the old pipeline had aligned, `SAMEA1904833` (one run, 49,908,521 records, 5.2x, 6.1 GB of uBAM), was run with `--no-upload` and compared with the old BAMs on Dropbox, in the environment `align` (minimap2 2.31, gatk4 4.6.2.0, pysam 0.24.1: the versions of the environment that made the old BAMs).

| Check | Result |
|---|---|
| `all.bam` and `chimeric.bam` against the old pipeline's | **identical**, record by record (49,908,521 and 408,251 records) |
| `TotalReads` (`samtools idxstats`) and `Coverage` against the old ledger row (GATK `CountReads`) | equal: 49,908,521 and 5.1871 |
| panel preflight, four keys; panel measurement | passed (`5d71ee8d6c2b`, `2af781e90119`, `80bbea57e10b`, `75367d0dfba8`); `panel.identity` says current |
| stages (minutes) | download 5.0 (20 MB/s), align 16.3, markdup 8.9, split 1.7, stats 1.3, slices 0.3, panel 0.2, qc 0.0: 0.5 h in all |
| resources | peak 22 of 88 cores, 27 of 1,099 GB RAM; peak local disk 4.6x the uBAM (the reservation was 3.5x) |
| slices | LG2 238 MB, LG10 181 MB, LG13 179 MB: 598 MB beside a 5.2 GB BAM (11.5%) |
| QC warning | `fully unmapped reads differ`: Picard's `UNMAPPED_READS` (1,896,658) is **all** unmapped primary reads, equal to samtools' `reads unmapped`, and the check compared it with the 582,860 reads without a position; **a defect of the check**, to be fixed |
| not exercised | the upload, the ledger write, the cleanup of local files and several samples in flight (`--no-upload`; the first batch does) |

The local folder of the sample still held the old pipeline's files (per-contig GVCF fragments, `unmapped.bam`): a realigned sample must not upload them (decision 48).


#### Decided (2026-11-07 to 2026-11-09)

1. **No GVCF stage.** v6 has none; the old script made them. Downstream code builds them on demand: `callset.py` and the BC code call `ensure_gvcf` (a cloud GVCF if there is one, HaplotypeCaller from the BAM otherwise; `paper.yaml`: "only used if a GVCF is missing"). The yellowhead batch was stopped in its GVCF stage, uploaded without GVCFs and registered (decision 35).
2. **The ledger has 22 columns** (confirmed 2026-11-09; decision 37; `alignment_db.LEDGER_COLUMNS`): SampleID, GenomeVersion, Species, RunIDs, NRuns, InputFingerprint, Coverage, UsableCoverage, TotalReads, MappedFraction, DuplicationRate, BamSize, Instrument, ReadLength, QCWarnings, PanelStatus, AlignedOn, Machine, PipelineVersion and the three tool versions. About 20 further statistics stay in `<sample>.qc.json`. `DuplicationRate` is Picard's PERCENT_DUPLICATION pooled over libraries, not the old `DuplicationReads` (the `.duplication.bam` split). `RunIDs` keeps the double-comma separator of the existing rows. `Species` is a copy made at alignment time and is **not** refreshed automatically (the audit flags a mismatch). The six split-read fractions are gone. **No quality thresholds** (coverage, mapping, duplication): the applications need different amounts of data, so `QCWarnings` flags pipeline faults only (decision 47).
3. **The rerun rule** (decided 2026-11-09, decision 48; **not built**). Reruns are rare (a sample with less data than expected gets more runs), and everything is realigned together. A sample is selected when the runs of its ledger row are not all the runs now in the DNAReadsDatabase (compared as sets; the `InputFingerprint` catches a uBAM re-made under the same RunID). A run in the database whose uBAM is not on Dropbox makes the sample wait and be reported. The sample is realigned from ALL its runs; its BAMs are replaced and **its row is replaced in place, never appended**; after the verified upload and the row, the **stale GVCF and the five old splits are deleted** from its Dropbox folder (each deletion logged), and only the files the pipeline wrote are uploaded (never the old files of the local folder). The old row's six fraction columns are cleared with it.
4. **Ledger writes.** v6's append with read-back is the base and is kept. Still to build: replacement of a row (item 3) and a local journal replayed when Google is unreachable (v6 raises after four attempts and leaves an uploaded sample without a row, which the audit and `register_alignments` repair). **Claims are not planned:** two machines get disjoint selections (decision 49).
5. **Supervision.** Built: a time limit on every transfer (it must make 5 MB/s, with a floor of 30 minutes), one hour on the quick commands, half an hour on the panel; a stop request kills every process group. Selection reports the Error and Exclude samples it skips. **(unverified: samples missing from the SampleDatabase.)**
6. **The panel stage (decisions 40, 45).** Once per run, before anything is aligned: `python -m sd_genotypes panel check --modules all`; the run does not start if it fails, nor, with `--panel-expect`, if the definition keys differ from the pilot's frozen ones. Per sample, after duplicate marking and while the BAM is on disk: `python -m sd_genotypes panel measure --bam <sample>.all.bam --sample <ID> --out-json <ID>.panel.json --modules all`, as a subprocess. The JSON (schema 2) records the BAM it was measured on (`bam`, `bam_size_bytes`); the pipeline checks it after every measurement: another BAM is a failed panel and the JSON is removed, no identity (an old `sd_genotypes`) is accepted with a QC warning. **A failed panel never fails the sample:** the BAM is uploaded, the row says `PanelStatus` failed with the reason, and module 04's own commands (`panel measure`, `retro`) measure it later; there is no genotype-only mode in 03 (decision 49). The panel JSON lives in the sample's Bamfiles folder beside the BAM, with the other files the calls need (`qc.json`, the Picard and samtools output). The recipe stays `minimap2 -x sr` (the gsdf calibration was made on those alignments, decision 17).
7. **Splits (decision 41).** Only the chimeric BAM is kept. The old splits and GVCFs of the existing sample folders are never used: they are deleted when a sample is realigned (item 3) and cleaned from the rest, together with a backfill of the new slices and the corrected duplication rates, by a separate tool (**(open)**).
8. **Slices (decision 42).** LG2, LG10 and LG13 as whole chromosomes (`slice_regions/Mzebra_GT3_NCBI.tsv`; `END` is the chromosome's length from the reference index; checked before a run starts); LG9, LG11 and LG20 are not sliced. **(open: whole chromosomes or 01's windows.)**
9. **One checkout, one environment (decision 43).** `align` (from `03_align_pipeline/env.yml`) pins minimap2 2.31, gatk4 4.6.2.0, pysam 0.24.1 and samtools 1.24, the versions of yellowhead's `genomics`, so that a cohort is aligned with one set of tools. Git lives in `genomics` only, and is needed only to clone and pull; the pipeline does not need it.
10. **The pilot.** The first real batch is the 8 unaligned `Deep_Benthic` samples (148 GB of uBAMs; IDs from the Blumer and Malinsky data), started on yellowhead on 2026-11-09 with `--samples-in-flight 4 --disk-reserve-factor 5.5`. The mbuna cohort follows: it calibrates the LG5 and B thresholds of the panel. Different machines get disjoint selections (decision 51).

**Known defects, to fix before the next batch:** the unmapped-reads check (compare Picard's `UNMAPPED_READS` with `reads unmapped`, tested on the numbers above); the default disk reservation (3.5x, observed 4.6x: 5.5x); uploading only the files the pipeline wrote.

**Work.**
1. **Done:** the port, its tests, the validation on one sample (above).
2. **Running:** the Deep_Benthic batch on yellowhead; it exercises the upload, the ledger write, the cleanup and four samples at once. Its result is not yet known.
3. The three fixes above; then the rerun logic of item 3, with the cleanup of stale files.
4. The journal for a Google outage.
5. The cleanup and backfill tool for the existing sample folders (item 7), and the port of `audit_alignments.py` and `register_alignments.py` to `reorg_spec` with the audit's new list of usual files.
6. After the first batches: the mbuna cohort, with every BAM measured while it is local; retire the old script.

**Done when.** One sample gives `all.bam` and a chimeric BAM identical to the old pipeline's (**met**); a batch runs on one machine without intervention, uploading, writing its rows and cleaning up (the Deep_Benthic batch); a killed run resumes without redoing finished samples; every sample has a ledger row and a panel JSON, or a recorded reason for the missing JSON; a sample with a new run is selected for rerun, replaced in place and its stale files deleted, and one without is not; a Google outage delays the ledger but loses nothing; tests pass offline.

### 04_genotype_sd_systems: gsdf duplication, transposon insertions, inversions, B chromosome

**Goal.** Per-sample genotypes of the SD systems, using user-provided definitions, normally produced by 03 while the BAM is local, and re-runnable on samples that missed it.

**Inputs.** A local indexed BAM (or a stored slice); definitions: the gsdf region and calibration, the transposon insertion definitions for gsdf and chr10 (sexloci's `known_sites.py`, copied and md5-pinned), inversion breakpoint windows (01), B-chromosome core blocks (its own discovery step).

**Outputs.** `<sample>.panel.json`; call tables; columns in `SampleDatabase` (for example `gsdf number`, `TE1`, `TE2`, `Inversion*`, B status) with a source and flag column for each.

**Interface.** `sd_genotypes panel measure --bam B --sample S --out-json J [--modules gsdf,te,lg5inv,bchrom]`; `sd_genotypes gsdf call|validate`; `python -m sd_genotypes.lg5inv.call --panel DIR --out T`; `python -m sd_genotypes.bchrom.score --panel DIR --out T`; `python -m sd_genotypes.te.call --panel DIR --out T [--gsdf-calls T]`; `python -m sd_genotypes.te.definition --ref GT3.fna --out-dir D`; `python -m sd_genotypes.retro --remote R --samples S --modules gsdf,te --out-dir D`; `python -m sd_genotypes.parity --sexloci sv_genotypes.tsv ...`; `sd_genotypes bchrom prepare|call|compare`; `python -m sd_genotypes.lg5inv prepare|status|clean` (pool reads for the LG5 markers); `sd_genotypes sheet-update`.

**The six goals, and how each is measured.** A genotype is made in two separate steps: a *measurement* (raw counts from reads, stored in `<sample>.panel.json` with the hash of its definition) and a *call* (thresholds and flags applied later to the stored counts, so a threshold can change without touching a BAM; measurements made under different definitions are never mixed). Each measurement runs in one or both of two modes. **Disk:** while the sample's BAM is local (03's pipeline calls `panel measure`). **Retro:** a retrospective measurement for a sample whose BAM is no longer local, which reads only the regions the measurement needs, from a stored slice or from the BAM on Dropbox without downloading all of it. (the mechanism is decided, decision 19; 03 stores slices of LG2, LG10 and LG13, decision 42)

| # | System | Evidence (the measurement) | Mode | Definition | State |
|---|---|---|---|---|---|
| 1 | gsdf duplication | depth ratio of the 21.7 kb tandem unit to its flanks (5 kb bins, region NC_135173.1:50.70-50.83 Mb) and unique-read (MAPQ >= 30) junction coverage | disk and retro | `config/gsdf.yaml` (calibration from the MC pedigree; minimap2 alignments) | rules aligned with sexloci 0.1.26; tests with fakes; MC family: 46 of 46 unflagged calls equal the sheet; **junction cross-check to be recalibrated; agreement over all 910 sheet values pending** |
| 2 | gsdf transposon insertion (chr7ins, directly upstream of gsdf; two reference copies) | read pairs crossing the junctions of a 603 bp junction allele at NC_135173.1:50,765,050 and :50,786,787, summed into one call | disk and retro | `sd_genotypes/te/sites.py` (copied from sexloci, md5 8b7b5cf1d8c7), `te_definition.json` (reference flanks, built once), `config/tepanel.yaml` | measurement, call and tests written (`te`); definition built from GT3a (md5 34f6c2b5dd25); chr7ins matches Munby's calls for 47 fish (3/3, 12/13, 24/24, 7/7; AUC 1.000); 18 samples reproduce VARIANTS.md; the sheet has no chr7ins column |
| 3 | two candidate transposon insertions on chr10 (TE1 +1,139 bp at NC_135176.1:24,921,178; TE2 +9,868 bp at :24,921,394) | read pairs crossing the junctions of each full allele | disk and retro | as 2 (md5 1957b81b6617 and d2b75a3b30e2) | TE1 133 of 133 and TE2 130 of 130 calls equal the sheet (146 samples; every non-reference sheet value included) |
| 4 | LG5 (orange blotch / ZW) inversion | ref and alt reads at 954 SNV markers inside the inversion (NC_135171.1:19.68-26.73 Mb) and reads containing the junction sequences at the right breakpoint | disk; retro only if 03 slices the inversion region (19.40-26.98 Mb) | `sd_genotypes/lg5inv/data/`, `config/lg5panel.yaml` | measurement and call written; measurement validated on real pool reads and run end to end on one real BAM; thresholds provisional |
| 5 | B chromosome | depth over the 874 core blocks (836 tier A) divided by the genome-wide depth, estimated from 12,000 sampled 50 bp bins | **disk only** (blocks and bins are genome-wide) | `sd_genotypes/bchrom/data/`, `config/bchrompanel.yaml` | measurement and call written, tests with fakes; scale estimate validated on the 21 tracks; `bedcov` agrees with mosdepth on a real BAM (17x: score ratio 1.0033) |
| 6 | inversions on LG2, 9, 10, 11, 13 and 20 | breakpoint-read evidence at the windows from 01 | disk; retro from 03's stored slices of LG2, LG10 and LG13 (LG9, LG11 and LG20 on disk only) | 01's windows | **deferred until the mbuna pilot is complete** |
| 7 | ecogroup placement, and a second method for 6 (proposed) | reference and alternative read counts at a basis set of common SNPs, projected onto fixed principal components | disk (retro to be examined) | `sd_genotypes/pca/` and `basis/<version>/` (proposed) | **proposed 2026-11-09; nothing built** |

The pilot starts with the 8 unaligned Deep_Benthic samples (running since 2026-11-09) and continues with the mbuna cohort. For the pilot, 4 and 5 and 1 to 3 are measured while the BAM is local; 6 is not started before it is complete. 7 is proposed (below).

**Current state, per goal.**
- *1. gsdf duplication* (`gsdf`, `panel`; tests with fakes): raw 5 kb bins and junction coverage stored, copy number called with the calibration of the MC pedigree (0.530 and 0.268). Aligned with sexloci 0.1.26, the code that the calibration and the existing `gsdf number` column came from: a 5 kb bin belongs to a window when it *starts* inside it (unit 7 bins, flanks 8 + 4; an earlier version of 04 took only bins entirely inside); the junction ratio divides by the first flank only; the reported number is the rounded estimate clipped to 1-6 and named CNn outside -/-, dup/-, dup/dup; a call is flagged only when the unrounded estimate is below 0.6 or above 5.4 (the 2-4 range of an earlier draft was an assumption, and wrong: five copies exist in the sheet). The junction cross-check is kept in a separate `check_flags` column and does not change whether a call counts as unflagged. `fit_anchors` and a theory calibration are available for another aligner. The calibration was made on **minimap2 -x sr** alignments (confirmed by the lab; the `bwa` in `VARIANTS.md` is wrong and should be corrected there). In sexloci: MC family 51 samples, all 46 unflagged calls right; agreement with Munby et al. Table S5 98.9% of 629 fish. **Real-data results (D05, samples streamed with `retro`): in the MC family 46 of 46 unflagged calls equal the sheet's `gsdf number` (females 25 of 25 at 2 copies, males 21 of 21 at 3 copies, classed male-specific); the five MC samples flagged here (MC_049, MC_060, MC_061, MC_075, MC_092) are among those the sheet leaves blank. Depth ratios do not overlap: 2 copies median 0.535 (0.463-0.584), 3 copies 0.796 (0.752-0.885). The junction cross-check expected 0.5 and 1.0 for 3 and 4 copies and is wrong in magnitude: the unique-read junction ratio is at most 0.054 in the 25 two-copy samples and has median 0.243 (0 to 0.442, zero in at least 10%) in the 21 three-copy samples, at a median flank depth of only 5.8x (minimum 2.5x), so the old rule would flag most true three-copy males. It only ever set `check_flags`, never the call, and has been replaced (decisions.md 22) by a one-sided rule: a unique-read junction ratio above 0.15 in a sample called with two copies (all 25 two-copy MC samples were at or below 0.054) is reported as `junction_in_two_copy`; a missing junction in a sample called with three or more copies is not checked until well-covered samples are measured, because true three-copy samples are often at 0 at 5-6x. SAMEA117806076 is called 4 copies (ratio 0.99, estimate 3.7) with no junction reads at about 18x; VARIANTS.md says three copies (source not established; the sheet has no gsdf value for the Munby fish). (unverified: the 98% agreement over all 910 sheet values.)**
- *2 and 3. Transposon insertions* (`te`; tests with fakes): the lab's IndelReadClassifier and the built-in sites of sexloci 0.1.26 are copied, byte-identical and md5-pinned, into `sd_genotypes/te/` (decision: copy, not depend; sexloci stays the reference and keeps the sheet writer and the retrospective runs it already does). For each site, the read pairs in the window anchor +- (250 + 1 + 200) bp (unmapped mates stored at their mate's position included; duplicate, secondary and supplementary excluded) are aligned (Smith-Waterman) to a local reference allele and insertion allele; a pair counts for an allele when a read clearly prefers it (score margin 20) and crosses a junction with 15 bp on each side. The REF and ALT pair counts are stored per site; the call is sexloci's binomial (P(alt) 0.02, 0.5, 0.98; fewer than 5 informative pairs gives ./.; GQ >= 10 counts as confident). The reference flanks are stored in `te_definition.json` (built once from GT3a by `te.definition`, which also checks every anchor base and the pinned fingerprints), so measuring a BAM needs no FASTA. chr7ins is called from the summed pairs of both reference copies; a `dilution_possible` flag marks a 0/0 with insertion pairs in a fish with 3 or more gsdf copies (the copy-number-aware model of VARIANTS.md 5.7 is not written). TE1 and TE2: how the original records were called is NOT RECORDED, the transposon family is unassigned, and phasing is not attempted. With the numpy aligner a read costs about 28 ms against the TE2 allele (parasail, optional, is about 100 times faster). **Real-data results (D05, `retro`, numpy backend): 18 samples reproduced the numbers quoted in VARIANTS.md (SAMEA117806076 chr7ins 63 reference : 6 insertion pairs `0/0`; AG_567_1 TE1 2 : 5 `0/1` GQ 13.3; YH_004_m, YH_009_m, YH_014_m without TE pairs) and agree with every sheet value (TE1 14 of 14, TE2 13 of 13, blanks included). chr7ins against Munby's calls for 47 fish: 3 of 3 homozygotes `1/1`, 12 of 13 heterozygotes `0/1`, 24 of 24 non-carriers `0/0`, 7 of 7 fish without a Munby call `0/0`, AUC 1.000 (16 carriers, 24 non-carriers), the sexloci result. The one discordance, SAMEA117806076, carries `dilution_possible` (4 gsdf copies): the flag works on real data. Targeted run on D05 (146 samples: every sample with a non-reference TE1 or TE2 value in the sheet, 93 and 62 carriers, plus 40 random reference samples): TE1 133 of 133 and TE2 130 of 130 confident calls equal the sheet, including 85 TE1 and 52 TE2 heterozygotes and 8 TE1 and 10 TE2 homozygotes, and the samples without a confident call are exactly those the sheet leaves blank (13 and 16). This is parity with sexloci (same code, same BAMs), not independent truth; the independent evidence is sexloci's pedigree association in YH (TE1 in 21 of 24 males and 0 of 17 females, TE2 in 22 of 26 and 0 of 20) and the pair counts quoted in VARIANTS.md, which were reproduced. (unverified: the remaining reference-state sheet values, 731 TE1 and 751 TE2 `0/0`, of which 40 and 68 were tested.)**
- *Retrospective mode and parity.* `sd_genotypes.retro` streams only the regions the modules need (their `regions(cfg)`: for gsdf and te, the 130 kb gsdf region, which contains both chr7ins windows, and one 1.1 kb TE1/TE2 window) from the BAM on Dropbox through a read-only `rclone serve http` and one `samtools view` per sample, sorts and indexes the small BAM, runs the panel on it and writes the usual `<sample>.panel.json` (resumable; failures in `retro_failures.tsv`; disk-only modules are refused before any network step). Verified on real data: the Dropbox layout is `Bamfiles/Mzebra_GT3_NCBI/<sample>/<sample>.all.bam` with `<sample>.all.bai` (1,337 sample folders; each also holds the chimeric, clipped, discordant, duplication, inversion and unmapped split BAMs, which the panel does not use), and 18 samples streamed without a failure. `sd_genotypes.parity` compares the calls with sexloci's `sv_genotypes.tsv` for the same samples (pair counts, genotypes, depth ratios, copy numbers).
- *One real BAM, all four modules* (D05; SAMEA117806270, 17x, 9.7 GB, `panel measure --modules gsdf,te,lg5inv,bchrom`): gsdf ratio 0.513, 2 copies, unflagged; chr7ins 0 ref : 12 alt `1/1` (a Kingiri homozygote), TE1 20 : 0 and TE2 17 : 0 `0/0` (parasail backend); LG5 median site fraction 0.000, junction 0/0/8, `std/std`; B score 1.075, `B-`. The disk measurement equals the streamed (retro) one: identical chr7ins, TE1 and TE2 pair counts and the same gsdf ratio. Cost with the BAM on local disk: gsdf 0.3 s, te 0.5 s, lg5inv 1.0 s, bchrom 4.3 s, each under 80 MB of memory; the four measurements together are 154 KB of JSON (bchrom 138 KB).

**Calling the panel from another pipeline (for example 03).** Preflight once per run: `python -m sd_genotypes panel check [--modules all]` loads every module's config and definition files, prints the definition key of each and exits with 1 if any is missing or broken (so a run stops before anything is aligned). Per sample, after mark-duplicates and while the BAM is on disk: `python -m sd_genotypes panel measure --bam <sample>.all.bam --sample <ID> --out-json <ID>.panel.json --modules all`, or from Python `sd_genotypes.panel.measure.measure(bam, sample, modules='all')` and `write_json`. Input: a coordinate-sorted BAM aligned to Mzebra_GT3_NCBI with its index next to it (`.bai`, also `<name>.bai`); the contig lengths are checked and a BAM of another reference is refused. Output: one JSON (schema 2), written atomically, with `bam` (the file name) and `bam_size_bytes` of the BAM measured, the definition key of every module and the repository version; a measurement belongs to a BAM only if both equal the file's name and size (a realigned sample usually keeps its name, so the size tells old from new): `python -m sd_genotypes.panel.identity --json J --bam B` exits 0 (current), 1 (stale, re-measure) or 2 (no size recorded), and `identity.same_bam` does the same from Python; schema 1 files written before this are still read; every call table carries `bam` and `bam_size_bytes`; for samples measured with `retro` they are those of the SOURCE BAM on Dropbox (an HTTP HEAD request), and `retro` skips a sample that already has a JSON, so after a realignment use `--refresh` or remove the JSON; a non-zero exit code and no JSON if a module fails (the caller decides whether to retry or skip). Needs only samtools and Python with numpy, pandas and pyyaml (environment `sdpanel`, `envs/panel.yml`; parasail optional); no Project, Dropbox or FileManager. The calls are made later, on a folder of collected JSONs, by `sd_genotypes gsdf call`, `sd_genotypes.te.call`, `sd_genotypes.lg5inv.call` and `sd_genotypes.bchrom.score`, which refuse to mix measurements made under different definitions. Samples whose BAM is not on disk are measured by `sd_genotypes.retro` (gsdf and te only). **Definition keys at the start of the pilot** (`panel check --modules all`, 2026-11-08): gsdf `5d71ee8d6c2b`, te `2af781e90119`, lg5inv `80bbea57e10b`, bchrom `75367d0dfba8`. A measurement carries its key; changing a measurement key in a `config/*.yaml` or a file under `sd_genotypes/*/data/` changes the key and splits the cohort, so these are not edited during the pilot (the thresholds of the calls are outside the key and can change at any time). The guide for running all of this is `04_genotype_sd_systems/docs/USAGE.md`.
- *4. LG5 inversion.* GT3 carries the standard arrangement; LatrZW and MezeOBm carry the inverted one, which has no standard haplotype in either assembly. At the right breakpoint (GT3 26,726,680) the junction is unique sequence and genotypable (inverted junctions a and b and the standard sequence C); the left breakpoint lies in a repeat and is not. 954 SNV markers (alleles shared by both inverted assemblies, clean in BB pools and Luwino, at about the inverted fraction in the OB pools of both genera) were chosen from five pools of PRJNA808791 aligned with the lab recipe (`lg5inv prepare`; D05). In pools the markers give 0.50-0.51 in the OB pools and 0.001-0.004 in the BB pools (selected on them), 0.53 and 0.57 in the OB pool of the other genus (held out); the held-out BB-type pools are not clean (0.06-0.20: the alt alleles also occur on standard chromosomes of other species), so a diploid fish is called on the median per-site alt fraction. The measurement reproduces the pool counts exactly through a second route. Calls: 0, 1 or 2 inverted chromosomes from the median fraction (cuts 0.25 and 0.75, **provisional**), with junction and aggregate cross-checks as flags. **(unverified: diploid fish; Luwino shows a uniform ~3% component that is either one inverted chromosome in 40 or contamination.)**
- *5. B chromosome* (`bchrom`): the discovery used Clark et al.'s 12 carrier females and 9 B-free pools (21 runs, aligned on D07 with the lab recipe). Core set on GT3: 874 blocks, 1.40 Mb (Clark: 619 blocks, 1.37 Mb), of which 836 (1.36 Mb) are tier A; 825 blocks lie on named LGs. Carriers score 4.6-8.2 (leave-one-out 4.3-8.0), B-free samples 1.08-1.09 (tier A: carriers 4.7-8.4, B-free 1.10-1.11); all three automatic checks pass. The measurement (`bchrom.measure`, `bchrom.score`) stores block depths and the sampled bins; the score is the median scaled coverage over tier A. The sampled-bin scale reproduces the exact genome-wide scale on all 21 tracks (largest deviation 0.8% of the scale and of the score; the exact route reproduces the discovery's own report). On a real BAM (SAMEA117806270, 17x, D05) the panel's score is 1.0751 against mosdepth's 1.0716 (ratio 1.0033; median block ratio 1.0037; 98.5% of blocks within 5%; genome-wide depth 15.66 in both) and the call is `B-`. **(unverified: thresholds 2.5 and 1.6 are provisional; mbuna males as negatives; this B-free single fish had 18% of its tier A blocks at 2x or more against 8-10% in the B-free pools, so `fraction_disagrees` (limit 0.30) may fire on poorer data.)**
- *6. Inversions on LG2, 9, 10, 11, 13 and 20*: breakpoints exist for LG2, 5, 9, 10, 11, 13 (01, run 20261105) and LG20 is not covered. Blumer et al. (Table S7) give chr2, 9, 10, 11, 13 for 1,375 samples, 1,364 of them in the SampleDatabase (columns `Inversion2` to `Inversion13`, section 7); all 440 mbuna among them are 0/0 on all five, so for the mbuna these are validation labels that mostly confirm absence, not an input. LG20 has no labels. The LG5 result is the guide: the first step for each breakpoint is a sequence-only check of how unique it is (one LG5 breakpoint was clean, the other lay in a repeat). The genotyper is not written. **Deferred until the mbuna pilot is complete.**

**Work.**
1. Pilot readiness: **done** (decisions 15 and 26-28). Measurement, call and definitions of goals 1 to 5 are committed, registered in `panel` and pass `panel check` (keys above); the B scale estimate and `bedcov` are validated; one real BAM went through all four modules on D05 with the disk measurement equal to the streamed one. The alignment machine (utaka) is reported ready by the lab (not verified here). Remaining for the pilot itself: refit the LG5 and B thresholds on the cohort's distributions, check the sex pattern of the calls (B absent from males; LG5 female-biased in the OB populations), and the agreement of gsdf with all 910 sheet values (not yet run; for te every non-reference value agrees).
2. 1: parity with sexloci on about 40 samples (raw depth ratios and copy numbers), then independent validation (Munby et al. Table S5; the MC family); run `retro` on the existing alignments.
3. 2 and 3: `te_definition.json` is built and committed (md5 34f6c2b5dd25); parity with sexloci's genotypes is shown for every non-reference sheet value; chr7ins is to be run on the pilot mbuna. Still open (VARIANTS.md section 7): phasing TE1 and TE2, a copy-number-aware chr7ins model, TE family, chr19-ins.
4. Slices are not needed for 1 to 3 (retro streams the BAM). 03 stores slices of LG2, LG10 and LG13 for goal 6 (decision 42); it would slice the LG5 inversion region only if goal 4 were ever to be measured retrospectively.
5. One `panel` JSON schema for all modules (done: gsdf, te, lg5inv, bchrom). The sheet write-back stays with sexloci (`--write-sample-database`: guarded, cell by cell, RAW values) unless 04 takes it over **(open)**.
6. After the pilot: 6.

**Done when.** 1: parity with sexloci (equal copy numbers on samples unflagged in both; depth ratios equal to rounding); unflagged calls agree with the existing column for at least 98% of the samples and females are 2 and males 3 in the MC family. 2 and 3: equal REF/ALT pair counts and agreeing genotypes with sexloci on the same samples (any difference listed and explained); chr7ins: all 16 of Munby's carriers show the junction and none of his 24 non-carriers, as in sexloci; agreement with the old `TE1`/`TE2` columns where they exist. 4: BB pools about 0 and OB pools about 0.5 inverted on the markers (reached in the pools); on the pilot the median site fraction forms clusters near 0, 0.5 and 1, the thresholds are refit on them, junction and marker evidence agree where the junction is informative, and calls are female-biased in the OB populations. 5: carriers score above the threshold leave-one-out (reached), mbuna males below it; the sampled-scale score within 2% of the exact score on the 21 tracks, and the first pilot BAM's measurement within 2% of its mosdepth score. 6: at least 98% agreement with Blumer on labelled samples, with hets reported separately. Every measurement carries its definition hash. **Reached so far:** gsdf, MC family 46 of 46 unflagged calls equal the sheet; chr7ins against Munby's 47 fish 3/3, 12/13, 24/24, 7/7 with AUC 1.000; LG5 pools 0.00-0.004 (BB) and 0.50 (OB), held-out OB pools 0.53 and 0.57; B carriers 4.3-8.0 leave-one-out against 1.1 for B-free samples, the sampled scale within 0.8% of the exact one and `bedcov` within 0.3% of mosdepth. TE1 and TE2: 133 of 133 and 130 of 130 calls equal the sheet, every non-reference sheet value included. **Not yet:** the gsdf agreement over all 910 sheet values and the reference-state TE values beyond the 40 and 68 tested; thresholds on the pilot cohort.

**Proposed: the PCA basis panel (goal 7, and a second method for goal 6) (2026-11-09; nothing is built).**

*Why.* Two questions are answered by the same kind of measurement: **are the inversion states of LG2, 9, 10, 11, 13 and 20 what the breakpoint reads say** (goal 6), and **does a sample sit in the ecogroup of its species** (goal 7). A principal component analysis of common SNPs does both: inside an inversion the first component separates the karyotypes 0/0, 0/1 and 1/1 (heterozygotes in the middle), and over the rest of the genome it separates the ecogroups. It is also the check that was made by hand for three samples whose Ecogroup is `Error` or `Exclude` (section 7).

*The reference data.* A polymorphism file of about 100 GB for many of the Malinsky species is available. **(unverified: its format, individuals, reference and content.)** The Malinsky et al. dataset has 134 individuals of 73 species across the seven main eco-morphological groups, 102 at about 15x and 32 at about 6x, with 30.6 million variants called against *Metriaclima zebra* assembly version 1.1, **not** against GT3: the sites have to be lifted.

*Two basis sets.*

| | Ecogroup placement | Inversion state (one set per inversion) |
|---|---|---|
| Sites | genome-wide autosomal SNPs outside the inversions, thinned for linkage (a few tens of thousands) | common SNPs inside the inversion's breakpoints (01's tables), thinned |
| Selection | biallelic; high minor-allele frequency across the lake and present in several ecogroups; a high call rate; no excess heterozygosity (paralogs); no indel nearby; unique, non-repeat sequence in GT3; mostly neutral, with a share chosen for differences between groups | as for the ecogroup set, within the interval |
| Reference labels | the ecogroup of each reference individual (the sheet, Blumer's clade) | Blumer's Table S7 for chr2, 9, 10, 11 and 13; **none for LG20** (its tandem inversions 20a and 20b need Kumar's labels or the assemblies to name the clusters) |
| Why separate | the inversions are polymorphic within ecogroups and would blur the group structure | the karyotype is the only structure that matters there |

*Building a basis (offline, once, on a machine with the file).* (1) Stream the file once per chromosome and write per-site statistics (allele frequency overall and per ecogroup, call rate, heterozygosity excess). (2) Select the sites by the criteria above. (3) Lift them to GT3 by mapping each site's flanking sequence with minimap2 and keeping only the sites whose reference allele matches; then drop those in repeats, near indels or of low mappability in GT3. (4) Write the reference genotypes of the reference individuals at the final sites and fix the principal components. (5) Store the result under `basis/<version>/` (sites, reference genotypes, loadings, group labels) with a definition file whose hash enters `panel check`, as for the other modules.

*Per sample.* A new panel module `pca` of `sd_genotypes` counts the reference and alternative reads at the basis sites (base and mapping quality limits in its config) while the BAM is on disk, and stores the counts in `<sample>.panel.json` beside the BAM. The alignment pipeline needs no change: it runs `panel measure --modules all`. **(unverified: the cost per sample; the sites are scattered over the genome, so the measurement is disk only unless `retro` can fetch them.)**

*Offline calls.* A probabilistic projection onto the fixed components, which tolerates low coverage and missing sites, gives a position, the nearest ecogroup and a confidence; a local projection per inversion gives a karyotype with thresholds taken from the labelled references. A sample whose ecogroup or karyotype disagrees with the sheet or with Table S7 is flagged, not changed. Hybrids and the intermediate species of the Malinsky paper fall between groups: that is a result, not an error. **(open: which sheet columns hold the calls.)**

*Validation.* Held-out species and Blumer's 1,375 samples (clade labels, and Table S7 for five of the six inversions) that are not in the basis file; reference individuals downsampled to 5x, since the pipeline's samples range from 5x to 48x; a comparison with goal 6's breakpoint reads where both exist. Choosing sites for group differences inflates the separation of the reference individuals, so only held-out samples count.

*Pilot.* One inversion first, to learn the number of usable sites, the cost and the separation before building the genome-wide set. **(open: which; one that is common in several ecogroups, so that all three karyotypes occur, for example LG9 or LG11.)**

*Dependencies.* 01 (the inversion intervals, and the alignment of assembly v1.1 to GT3 if the lift uses it), 03 (BAMs), the SampleDatabase (ecogroups) and Table S7 (labels).

*Done when (proposed).* On held-out samples the projection puts at least 95% in their own ecogroup (the intermediate species excepted and listed); for the pilot inversion the calls agree with Table S7 for at least 98% of the labelled samples; every measurement carries its definition hash; tests pass offline.

**Requests from 03 to 04 (2026-11-08).** (1) The panel JSON should record the identity of the BAM it was measured on (file name and size at least), so that a re-alignment makes the old measurement detectable as stale; to be asked of the genotyping side. (2) 03 treats a non-zero exit of `panel measure` as "no JSON, sample not failed"; the exit codes and the atomic write already satisfy this. (3) 03 needs `sd_genotypes` and its dependencies installable from the same checkout into the environment of the alignment machine (item 10 of module 03).

**Reply from 04 (2026-11-09).** (1) Done: the panel JSON (schema 2) records `bam` (the file name) and `bam_size_bytes`; for samples measured with `retro` they are those of the source BAM on Dropbox; `python -m sd_genotypes.panel.identity --json J --bam B` exits 0 (current), 1 (stale: re-measure) or 2 (no size recorded); schema 1 files are still read (decision 45). 03 decides what to do with a stale measurement; the measurements already made by `retro` are schema 1 and report as unknown. (2) Confirmed: `panel measure` exits non-zero and writes no JSON if a module fails, and the JSON is written atomically. (3) The panel needs only samtools and Python with numpy, pandas and pyyaml (`envs/panel.yml`, environment `sdpanel`) and runs from a checkout after `source tools/env.sh`; `sd_genotypes` is not packaged for pip.


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

**Three flagged samples (2026-11-09).** `SAMEA2661230` (*Ctenopharynx intermedius*), `SAMEA3388870` (*Dimidiochromis strigatus*) and `SAMEA117806170` (*Otopharynx argyrosoma*) carry an `Error` or `Exclude` marker, which the lab put there because they did not look like the right species in a PCA. Their species is the same in SRA, in Blumer's Table S1 and in the sheet, and in Blumer's tree (`Data_S1.nwk`) each sits with its own species or genus (10 of the 10 nearest tips conspecific for *O. argyrosoma*, 4 of 4 for *D. strigatus*, the two other *C. intermedius* as nearest tips); in Blumer's table all three are in the Shallow clade at 16.6 to 17.7x. Two of them have only `all.bam` and its index in their Bamfiles folder (the partial uploads of the 2026-11-07 audit); the third has a complete old output. The marker is in different columns for different samples (SAMEA2661230 reads Ecogroup `Tilapiine` with Subgroup `Error`), and the pipeline skips on **Ecogroup** only. A redownload of their reads is on hold; whether their uBAMs are intact is not yet checked. The PCA basis panel (module 04) would test such cases systematically.

**Module 04 and the SampleDatabase.** The genotype columns (`gsdf number`, `TE1`, `TE2`, `Inversion*`, B status) are written by sexloci's guarded sheet writer (`--write-sample-database`, cell by cell, raw values) unless 04 takes the write-back over **(open)**. Any writer follows the rule at the top of this section.

## 8. Machines and environments

| Machine | Role | Notes |
|---|---|---|
| utaka | alignment and heavy work (96 cores); the machine for the pilot (reported ready by the lab, not verified here) | branch `mcgrath_dev`, Python 3.7 environment `CichlidSRSequencing`; the old file manager (no download-database methods, no retry); the alignment audit of 2026-11-07 ran here |
| mzebra | downloads, conversions, sheet maintenance (96 cores, 503 GB RAM) | login `bio-mcgrath-dropbox`; environment `enadl` (modules 02 and 04, decision 23); the download campaign ran here and the checks of this repository run here |
| D07 (workstation) | the B-chromosome discovery (32 cores, 62 GB RAM) | Aspera blocked; SRA about 11 MB/s with 8 connections (2.7 MB/s with one) |
| D05 (workstation) | module 04's real-data checks of 2026-11-08 (retrospective runs, one BAM through four modules) | rclone from snap, its configuration outside `~/.config/rclone` (decision 23) |
| yellowhead | validation of the pipeline (student machine); the first real batch | the 45-sample batch of 2026-11-07 was uploaded and registered; the new pipeline was validated here on 2026-11-09 (environment `align`, from `03_align_pipeline/env.yml`; git only in `genomics`) and runs the Deep_Benthic batch; its old local BAMs stay as a reference until the validation is repeated if needed |
| laptop | development, tests, publishing | tests run offline; the `syri_env` environment has an installed package called `tests` that shadows the repository's, so run `python3 -m unittest discover -s tests`, not `tests.<module>`; `psutil` and `pysam` are missing there and on mzebra, so three older tests of the old alignment code cannot import |

Each module has an `env.yml`. `envs/download.yml` (`enadl`) serves modules 02 and 04 and now also holds matplotlib and parasail-python (`tools/check_env.sh` verifies it); `sdpanel` (`envs/panel.yml`) is the minimal panel environment. A known conflict: `aria2` and `sra-tools` cannot be solved in one conda environment on macOS; on Linux they can.

## 9. Order of work

1. **Done:** module 02 for both Run Info files; the Blumer update of the SampleDatabase; the yellowhead batch uploaded and registered (AlignmentDatabase 1,296 rows); the panel's goals 1 to 5 built and validated on real data, their definitions frozen for the pilot (decision 28); **module 03 built and validated on one real sample** (identical BAMs, equal statistics, the panel current; section 6).
2. **Now:** the first real batch of the new pipeline, the 8 unaligned Deep_Benthic samples, running on yellowhead (upload, ledger write, cleanup and four samples at once are its test); finish 01 (LG20, schemas, the breakpoint windows that goal 6 and 03's slices need); decide the 186 samples without a sample row and the three flagged samples.
3. **Next (module 03):** the three fixes (the unmapped-reads check, the disk reservation, uploading only what the pipeline wrote), then the rerun logic with the cleanup of stale files, then the journal for a Google outage.
4. **The mbuna cohort** (about 477 samples), every BAM measured by the panel while local: refit the LG5 and B thresholds on the cohort's distributions, check the sex pattern of the calls and the agreement of gsdf with all 910 sheet values.
5. **The PCA basis panel** (module 04): describe the 100 GB file, then a pilot on one inversion, then the genome-wide set; it verifies the inversion states of goal 6 and places samples in their ecogroups.
6. **After the pilots:** goal 6 on the slices and on disk; then 05 and 06, 08, 07; 09 and 10 read everything.

Critical path to the mbuna cohort: the Deep_Benthic batch clean, the three fixes, 01's windows (for the slices and goal 6), the Ecogroup of every sample set.

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

**Rules for the move.** It happens in one change, after the pilots and while no download or alignment is running on any machine; the downloader launches its workers by file path, and the alignment run imports the old modules. Thin wrapper scripts under today's names stay for one release. The test suites are the safety net: they import by name and fail loudly if anything was missed. The package name `yhsd` is retired in favour of the module packages above. This rule governs retiring the original layout from the machines; the reorganised layout already exists on `reorg_spec`, and new code (module 03) is built there before the pilot (decision 43).

**Status (2026-11-08).** The layout of section 3 exists on branch `reorg_spec` (module 04 was built there); the original layout stays on `YH_SD_Paper` and `mcgrath_dev`. **Module 03 now moves with the others** (decision 43): the constraint behind decision 11, an alignment run that used the old file manager, has lapsed as far as is known (the utaka run ended with its rows registered; the yellowhead run was stopped), to be confirmed that no other alignment job runs. The old scripts stay usable through wrappers and move to `legacy/` after the pilot. The tools of 2026-11-07 that exist only on `YH_SD_Paper` have to be ported (section 3).

## 11. Decisions and measurements

The log is `docs/decisions.md` (newest last, entries 1 to 54; an entry that is refined later says so, and the later entry names it). The entries that shape this document:

* Downloads: 1 (NCBI's copy when ENA is slow), 2 and 31 (every SRA run checked; counts), 4 (binned qualities), 6 (the DownloadDatabase), 30 (bulk row addition; whole-sheet rewrites only with a single writer), 32 (the campaign is complete), 33 (supervision; the retry cannot break its program).
* Alignment: 34 and 35 (the utaka crash and its recovery; the yellowhead batch), 36 (no GVCF), 37 (the ledger), 40 to 43 (the panel stage, the splits, the slices, the location of the code), 45 (the panel JSON records its BAM), 46 to 52 (validation of the new pipeline, no QC thresholds, the rerun rule, disjoint selections, the environment, the pilot, the defects found).
* SampleDatabase: 38 (Blumer is the correct species; `SRASpecies`), 53 (the three flagged samples).
* Genotyping: 13 to 29 (measurement and call are separate; the panel contract; definitions frozen for the pilot; real-data results), 54 (the PCA basis panel, proposed).
* Layout: 10 to 12, 43 and 44.

## 12. Open decisions

1. LG20 and the other proposals of the module 01 finish plan: confirm or change them. (01)
2. The reference-state convention for LG5 (`reference_is_standard` unset). (01)
3. TE1 and TE2: phasing, the transposon family, a copy-number-aware chr7ins model, chr19-ins (VARIANTS.md section 7). (04)
4. Who writes the genotype columns back to the SampleDatabase: sexloci or 04. (04)
5. Variant caller, joint or per-family calling, filters, region scope. (05)
6. Species list and sample-size rules for SD discovery. (06)
7. QTL software. (07)
8. Annotation and mappability sources for ASE. (08)
9. Web server hosting, audience and data store. (09)
10. The figure and table list. (10)
11. When the reorganised layout is merged into the main branch, and who updates each machine.
12. Whether `SPEC.md` is edited directly in the repository, or changes go through the decision log.
13. Slices: whole chromosomes or 01's windows; LG9, LG11 and LG20 are not sliced. (03)
14. Whether runs can be excluded individually, or only whole samples. (03)
15. Whether utaka's environment can host the panel or needs a new one. (03)
16. The cleanup and backfill tool for the existing sample folders: new slices, corrected duplication rates, deleting the old splits and GVCFs. (03)
17. The port of `audit_alignments.py` and `register_alignments.py` to `reorg_spec`, with the audit's new list of usual files (GVCFs and the five old splits are no longer expected). (layout)
18. The 186 samples without a sample row, and an Ecogroup name for outgroups. (SampleDatabase)
19. A metadata source for the remaining samples without Ecogroup, Sex and LabReared (the Malinsky supplement); `Subgroup` and `Category` for new samples. (SampleDatabase)
20. The three flagged samples (section 7): whether to check their uBAMs, whether to realign them with `--include-excluded`, and what the `Error` marker means in the Subgroup column. (SampleDatabase)
21. ERR14734399: stay excluded or fetch from ENA; the 4 runs accepted without the ENA comparison. (02)
22. The PCA basis panel (04): the 100 GB file (format, individuals, reference, genotypes or frequencies); where the module lives (a panel module of `sd_genotypes` is proposed); the site-selection thresholds; the inversion for the pilot; the reference individuals and their ecogroup labels; the source of LG20 labels; which sheet columns hold the calls.

## Appendix A. For the methods

Reads of the samples of Blumer et al. and Malinsky et al. that were not already in the lab's database (334 runs) were retrieved from NCBI's copy of the Sequence Read Archive (public AWS bucket) because transfers from ENA were about 0.3 MB/s. Each run was validated (`vdb-validate`), converted with `fasterq-dump`, and compared with the ENA record: the number of read pairs agreed with ENA's `read_count` within 5% (SRA held 0.01% to 1.71% more spots than ENA in every run), and the sequence and base qualities of the first 1,000 reads were identical to ENA's FASTQ in 323 of 327 runs checked (the remaining 4 were not comparable at the time and are to be rechecked). One run, ERR14734399, whose first reads did not match ENA's at all, was excluded. Species and sex were taken from Blumer et al.'s Table S1 and inversion genotypes from their Table S7 for the samples they describe; the species names recorded in SRA are kept separately.
(Wording to be confirmed once the 4 unverified runs are rechecked and the final run counts are fixed.)
