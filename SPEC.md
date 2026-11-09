# CichlidSRSequencing: specification

Status: v4 (2026-11-08). Owner: McGrath lab. This file lives at the repository root and is the source of truth for what each part of the repository is for, what it consumes, what it produces, and when it is done. It merges two working versions: v3 (2026-11-07: downloads, the alignment redesign, the SampleDatabase) and the genotyping version (draft v2 of 2026-11-06, with module 04 as it stood on 2026-11-08). The log of decisions is `docs/decisions.md` (section 11).
Where a statement is not yet verified it is marked **(unverified)**; where a decision is still to be made it is marked **(open)**.

## 0. What changed since v3

* **The two versions are merged.** Module 04 is the genotyping version (the six goals, the pipeline-facing panel contract, retrospective measurement, the real-data results of 2026-11-08); modules 02 and 03 and the SampleDatabase section are v3's, updated.
* **Module 03 has its decisions of 2026-11-08:** a panel stage that never fails a sample, only the chimeric split kept, slices of LG2, LG10 and LG13, and the new alignment code built in the reorganised layout next to `sd_genotypes` (section 6, module 03; decisions 40 to 43).
* **Module 01 has a proposed finish plan** (section 6, module 01), awaiting confirmation: LG20 as two tandem inversions, a reference-state column, rules for the LG9 and LG13 questions, the published Malawi assemblies as the additional assemblies, schemas.
* **The reorganised layout exists** on branch `reorg_spec`; `sd_genotypes` was built there. Section 3 says which branch holds what.
* **The decision log is merged** into `docs/decisions.md`: entries 1 to 29 (layout, downloads, genotyping) and 30 to 44 (downloads finished, the alignment crash and its recovery, the ledger, the SampleDatabase, the decisions of 2026-11-08).

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
| Alignment pipeline, v6 (offline tests, never run on real data): `alignFastQReads_v2.py` with `seq_pipeline.py`, `bam_split.py`, `bam_slices.py`, `bam_stats.py`, `alignment_db.py`, `disk_ledger.py`, `download_manager.py`, `slice_regions/`, `unit_scripts/validate_against_old.py`; described in `PIPELINE_V2_README.md` | `./`, `helper_modules/` |
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

**Interface.** Today: `alignFastQReads.py <genome> -e|-u|-p|-c|-s ...` (old) and `alignFastQReads_v2.py <genome> [--samples-in-flight N --align-threads N --uploads N --list --no-upload --setup-database --include-excluded]` (v6). Target: `align_pipeline run <genome> [-e ecogroup] [-u subgroups] [--samples ...] [--list] [--force]`, `align_pipeline setup-database`, `align_pipeline genotype-only <samples>`.

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

#### Decided on 2026-11-07 and 2026-11-08, set against v6

1. **No GVCF stage.** v6 has none; the old script still makes them. Downstream code builds them on demand: `callset.py` and the BC code call `ensure_gvcf` (a cloud GVCF if there is one, HaplotypeCaller from the BAM otherwise; `paper.yaml`: "only used if a GVCF is missing"). The yellowhead batch was stopped in its GVCF stage, uploaded without GVCFs and registered (decision 35).
2. **The ledger keeps** SampleID, GenomeVersion, RunIDs, Coverage, TotalReads, the version columns (minimap2, gatk, pysam) and BamSize; **adds** DuplicationRate and Species; **drops** the six split-read fractions. The discussion was held without reference to v6's own new columns, which are many (below), so the reconciliation is **open**.

| Column | Agreed or proposed | Relation to v6 |
|---|---|---|
| DuplicationRate | agreed | v6 `DuplicateFraction` (Picard's definition, pooled over libraries); not the old `DuplicationReads`, which counted the `.duplication.bam` split (structural-duplication evidence) |
| Species | agreed | new: a copy of `SampleDatabase.Species` at alignment time; the audit flags a mismatch, a refresh command updates it |
| MappedFraction, QC, PipelineVersion | proposed, welcomed | v6 has `MappedFraction`, `QCWarnings` and `PipelineVersion` (the text `v2`; to become the git commit) |
| Usable coverage | proposed | v6 measures `CoverageNonDuplicate` (and `CoverageMapped`) instead of a derived formula |
| NRuns, InputFingerprint, AlignedOn, Machine, Instrument, ReadLength | proposed | new (the fingerprint is a hash of the sorted (RunID, uBAM size) pairs) |
| about 20 other v6 statistics (`PrimaryReads`, `UnmappedFraction`, `UnmappedNoPositionFraction`, `UnmappedOrphanFraction`, `DuplicateFractionFromBam`, `OpticalDuplicateFraction`, `MeanReadLength`, `MappedBases`, `ProperlyPairedFraction`, `MQ0Fraction`, `InsertSizeMean`, `InsertSizeSD`, `ErrorRate`, `UBamReads`) | **(open)** | proposal: keep them in the per-sample `qc.json` beside the BAM and keep the sheet lean; `QCWarnings` stays in the sheet |

   Dropping the six fractions means changing v6's required legacy columns (`LEGACY_KEYS`) and the stage that counts the splits. Of the splits only the chimeric BAM is still produced (item 7). Panel and genotype results do not belong in the ledger: a changed panel definition must never trigger a realignment.
3. **The rerun rule (new in either pipeline).** A sample needs (re)alignment if it has no row; or its run set now differs from the row's (runs added or removed; the old RunIDs separator is a double comma and is parsed as a set); or its InputFingerprint differs (a uBAM re-created under the same RunID, for example re-downloaded from SRA instead of ENA); `--force` overrides. A run that is in the database but whose uBAM is not yet on Dropbox makes the sample wait and be reported, never be aligned from a subset. A rerun replaces: upload the new set, verify, update the row in place (v6 only appends). A change of tool versions or reference is recorded but triggers nothing unless asked. The panel JSON carries the alignment fingerprint so stale genotype results are detectable.
4. **Ledger writes.** v6's append with read-back is the base. Still to add: update in place for reruns; a local journal replayed when Google is unreachable (v6 raises after four attempts and leaves an uploaded sample without a row, which the audit and `register_alignments` then repair); claims so two machines do not align the same sample (nothing in v6 prevents it). **(open: ledger only, or also a queue with status and claims, as the DownloadDatabase has.)**
5. **Supervision and reported skips.** Every subprocess has a timeout or stall limit and dies with its parent; selection lists the samples it leaves out and why (no SampleDatabase row, blank Ecogroup, missing uBAM). **(unverified: how far v6 already does this; `--list` previews the selection.)**

6. **The panel stage (decision 40).** Once per run, before anything is aligned: `python -m sd_genotypes panel check --modules all` (loads every module's config and definition, prints the four definition keys, exits 1 if any is missing; the keys are logged with the run). Per sample, after duplicate marking and while the BAM is on disk: `python -m sd_genotypes panel measure --bam <sample>.all.bam --sample <ID> --out-json <ID>.panel.json --modules all`, as a subprocess; no panel logic is copied into 03. Input: the coordinate-sorted, duplicate-marked BAM aligned to `Mzebra_GT3_NCBI` with its index next to it (gsdf and te exclude duplicate, secondary and supplementary reads, so the duplicate flags must be set). Output: one atomic JSON (schema 1) carrying the definition key of every module, or a non-zero exit and no JSON. Measured cost per sample: 0.3 s (gsdf), 0.5 s (te), 1.0 s (lg5inv), 4.3 s (bchrom), under 80 MB each, 154 KB of JSON. **A panel failure does not fail the sample:** the BAM is the primary product, so the sample is uploaded and registered with a `panel_failed` warning (`QCWarnings`) and `genotype-only` fills it in later. The alignment recipe must stay `minimap2 -x sr`: the gsdf calibration was made on those alignments (decision 17). The JSON does not yet record which BAM it was measured on, which the rerun rule needs to detect a stale panel after a re-alignment: **requested of 04 (open)**; until then 03 writes a sidecar with the BAM's size and the InputFingerprint.
7. **Splits (decision 41).** Only the chimeric BAM is kept; the panel does not use the other five. For v6 this means changing the required legacy columns (`LEGACY_KEYS`), the single pass that writes the splits, and the stage that counts them; for the audit, a new list of usual files. The 1,337 existing sample folders on Dropbox keep their old splits; whether to delete them is **(open)**.
8. **Slices (decision 42).** 03 stores slices of LG2, LG10 and LG13 for goal 6 of 04, as v6 already does for LG10 (whole chromosomes; **(open: whole chromosomes or 01's windows)**); contig names come from the genome registry, not from an assumed numbering. LG9, LG11 and LG20 are not sliced: goal 6 on them is measured on disk only (or after a re-download). Three of 22 chromosomes are about an eighth of `all.bam` per sample (estimate), which matters for Dropbox space.
9. **Genotype-only mode.** For samples whose BAM is not local: gsdf and te through `sd_genotypes.retro` (it streams only the regions they need from the BAM on Dropbox through a read-only `rclone serve http` and `samtools view`; no BAM download); LG5 and the B chromosome read megabases or the whole BAM, so the BAM is downloaded, measured and deleted.
10. **One checkout, one environment.** 03 is built in the reorganised layout beside `sd_genotypes` (decision 43), so a machine that aligns checks out one branch. One environment file for 03 holds the alignment tools (minimap2, samtools, the GATK/Picard tools used for MergeBamAlignment and MarkDuplicates) and the panel's dependencies (numpy, pandas, pyyaml; parasail optional); `sdpanel` stays the minimal environment for panel-only use. The interpreter that runs `sd_genotypes` is a per-machine setting. **(unverified: whether utaka's Python 3.7 environment can host it; a new environment there is the likely answer.)**

**Work.**
1. Validate v6 on real data: the BAMs of the yellowhead run are the reference (`validate_against_old.py`: the seven BAMs of a sample must be identical); keep them until then. With the new splits policy the comparison covers `all.bam` and the chimeric BAM.
2. Port v6 into `03_align_pipeline/align_pipeline/` with one `cli.py`, tests, and an `env.yml`.
3. The panel stage (item 6), tested with a fake `sd_genotypes` that follows the contract's exit codes and file names; then the genotype-only mode (item 9).
4. The chimeric-only split and the slice list (items 7 and 8); change `LEGACY_KEYS`, the header setup and the audit's file list.
5. The ledger columns, the rerun rule, update in place, the journal, claims (items 2 to 4).
6. Supervision and reporting of skipped samples (item 5).
7. A two-machine procedure for the mbuna samples (about 477 samples, 7 TB of uBAMs), documented and rehearsed; the pilot is the first mbuna alignments, with every BAM measured while it is local.
8. Retire the old script after the pilot.

**Done when.** One sample gives `all.bam` and a chimeric BAM identical to the old pipeline's; a batch of mbuna samples runs on two machines without intervention; every sample has a ledger row and a panel JSON, or a recorded reason for the missing JSON; a killed run resumes without redoing finished samples; the ledger row of every sample states its runs, its fingerprint and its pipeline version; a sample with a new run is selected for rerun and one without is not; a Google outage delays the ledger but loses nothing; tests pass offline.

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

The pilot is the first mbuna alignments. For the pilot, 4 and 5 and 1 to 3 are measured while the BAM is local; 6 is not started before it is complete.

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

**Module 04 and the SampleDatabase.** The genotype columns (`gsdf number`, `TE1`, `TE2`, `Inversion*`, B status) are written by sexloci's guarded sheet writer (`--write-sample-database`, cell by cell, raw values) unless 04 takes the write-back over **(open)**. Any writer follows the rule at the top of this section.

## 8. Machines and environments

| Machine | Role | Notes |
|---|---|---|
| utaka | alignment and heavy work (96 cores); the machine for the pilot (reported ready by the lab, not verified here) | branch `mcgrath_dev`, Python 3.7 environment `CichlidSRSequencing`; the old file manager (no download-database methods, no retry); the alignment audit of 2026-11-07 ran here |
| mzebra | downloads, conversions, sheet maintenance (96 cores, 503 GB RAM) | login `bio-mcgrath-dropbox`; environment `enadl` (modules 02 and 04, decision 23); the download campaign ran here and the checks of this repository run here |
| D07 (workstation) | the B-chromosome discovery (32 cores, 62 GB RAM) | Aspera blocked; SRA about 11 MB/s with 8 connections (2.7 MB/s with one) |
| D05 (workstation) | module 04's real-data checks of 2026-11-08 (retrospective runs, one BAM through four modules) | rclone from snap, its configuration outside `~/.config/rclone` (decision 23) |
| yellowhead | validation of the pipeline (student machine) | the batch of about 45 samples was stopped in its GVCF stage on 2026-11-07 and uploaded and registered; its local BAMs stay as the reference for validating v6 |
| laptop | development, tests, publishing | tests run offline; the `syri_env` environment has an installed package called `tests` that shadows the repository's, so run `python3 -m unittest discover -s tests`, not `tests.<module>`; `psutil` and `pysam` are missing there and on mzebra, so three older tests of the old alignment code cannot import |

Each module has an `env.yml`. `envs/download.yml` (`enadl`) serves modules 02 and 04 and now also holds matplotlib and parasail-python (`tools/check_env.sh` verifies it); `sdpanel` (`envs/panel.yml`) is the minimal panel environment. A known conflict: `aria2` and `sra-tools` cannot be solved in one conda environment on macOS; on Linux they can.

## 9. Order of work

1. **Done:** module 02 for both Run Info files; the Blumer update of the SampleDatabase; the yellowhead batch uploaded and registered (AlignmentDatabase 1,296 rows); the panel's goals 1 to 5 built and validated on real data, their definitions frozen for the pilot (decision 28).
2. **Now:** finish 01 (LG20, schemas, the breakpoint windows that goal 6 and 03's slices need); decide the 186 samples without a sample row; recheck the 4 unverified runs; check what `reorg_spec` holds of the tools of 2026-11-07 and port the rest.
3. **Module 03, built in the reorganised layout:** validate v6 against the yellowhead BAMs; then the panel stage, the chimeric-only split, the slices of LG2, LG10 and LG13, the ledger columns and the rerun rule. This gates the pilot.
4. **The pilot:** the first mbuna alignments (about 477 samples), every BAM measured by the panel while local; refit the LG5 and B thresholds on the cohort's distributions, check the sex pattern of the calls, and the agreement of gsdf with all 910 sheet values.
5. **After the pilot:** goal 6 (inversions on LG2, 9, 10, 11, 13, 20) on the slices and on disk; then 05 and 06, 08, 07; 09 and 10 read everything.

Critical path to the pilot: 03 validated and extended, 01's windows (for the slices and goal 6), the Ecogroup of every sample set.

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

The log is `docs/decisions.md` (newest last, entries 1 to 44; an entry that is refined later says so, and the later entry names it). The entries that shape this document:

* Downloads: 1 (NCBI's copy when ENA is slow), 2 and 31 (every SRA run checked; counts), 4 (binned qualities), 6 (the DownloadDatabase), 30 (bulk row addition; whole-sheet rewrites only with a single writer), 32 (the campaign is complete), 33 (supervision; the retry cannot break its program).
* Alignment: 34 and 35 (the utaka crash and its recovery; the yellowhead batch), 36 (no GVCF), 37 (the ledger and the rerun rule), 40 to 43 (the panel stage, the splits, the slices, the location of the code).
* SampleDatabase: 38 (Blumer is the correct species; `SRASpecies`).
* Genotyping: 13 to 29 (measurement and call are separate; the panel contract; definitions frozen for the pilot; real-data results).
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
13. The ledger columns: reconcile the agreed set with v6's ~20 statistics columns (sheet or `qc.json`). (03)
14. The BAM identity in the panel JSON (requested of 04), and what the sidecar holds until then. (03, 04)
15. Slices: whole chromosomes or 01's windows; LG9, LG11 and LG20 are not sliced. (03)
16. Whether to delete the old splits of the 1,337 existing sample folders on Dropbox. (03)
17. Whether runs can be excluded individually, or only whole samples. (03)
18. The ledger as a result sheet only, or also a queue with claims. (03)
19. How a stale `Species` in the ledger is refreshed. (03)
20. Whether utaka's environment can host the panel or needs a new one. (03)
21. Whether `reorg_spec` holds the tools of 2026-11-07, and the porting of the rest. (layout)
22. The 186 samples without a sample row, and an Ecogroup name for outgroups. (SampleDatabase)
23. A metadata source for the remaining samples without Ecogroup, Sex and LabReared (the Malinsky supplement); `Subgroup` and `Category` for new samples. (SampleDatabase)
24. ERR14734399: stay excluded or fetch from ENA; the 4 runs accepted without the ENA comparison. (02)

## Appendix A. For the methods

Reads of the samples of Blumer et al. and Malinsky et al. that were not already in the lab's database (334 runs) were retrieved from NCBI's copy of the Sequence Read Archive (public AWS bucket) because transfers from ENA were about 0.3 MB/s. Each run was validated (`vdb-validate`), converted with `fasterq-dump`, and compared with the ENA record: the number of read pairs agreed with ENA's `read_count` within 5% (SRA held 0.01% to 1.71% more spots than ENA in every run), and the sequence and base qualities of the first 1,000 reads were identical to ENA's FASTQ in 323 of 327 runs checked (the remaining 4 were not comparable at the time and are to be rechecked). One run, ERR14734399, whose first reads did not match ENA's at all, was excluded. Species and sex were taken from Blumer et al.'s Table S1 and inversion genotypes from their Table S7 for the samples they describe; the species names recorded in SRA are kept separately.
(Wording to be confirmed once the 4 unverified runs are rechecked and the final run counts are fixed.)
