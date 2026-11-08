# CichlidSRSequencing: specification

Status: draft v2 (2026-11-06): the layout of section 3 now exists (see MIGRATE.md); module 04 downloads the B-chromosome reads through the SRA source. Owner: McGrath lab. This file lives at the repository root and is the source of truth for what each folder is for, what it consumes, what it produces, and when it is done.
Where a statement is not yet verified it is marked **(unverified)**; where it is a decision still to be made it is marked **(open)**.

---

## 1. Purpose

Take public and in-house short-read data for Lake Malawi cichlids from a read archive to a set of results about sex determination (SD): where the SD loci are, what they are made of (inversions, transposon insertions, a gene duplication, a B chromosome), how they behave in backcross families, and which variants are candidate causal changes. Everything that goes into the paper must be reproducible from versioned code and recorded inputs.

The work is ten modules. Each is one folder. They communicate only through **data contracts** (section 4), never by importing each other's internals, so any module can be re-run, replaced or tested alone. A module may import the package of an upstream module it declares a dependency on (04 uses `download_reads.sources`; 06 and 07 use `variant_calling.callset`; all use `cichlid_common`), never its internals.

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

**Inputs.** A local indexed BAM (or a stored slice); definitions: the gsdf region and calibration, the transposon insertion definitions for gsdf and chr10 (sexloci's `known_sites.py`, copied and md5-pinned), inversion breakpoint windows (01), B-chromosome core blocks (its own discovery step).

**Outputs.** `<sample>.panel.json`; call tables; columns in `SampleDatabase` (for example `gsdf number`, `TE1`, `TE2`, `Inversion*`, B status) with a source and flag column for each.

**Interface.** `sd_genotypes panel measure --bam B --sample S --out-json J [--modules gsdf,te,lg5inv,bchrom]`; `sd_genotypes gsdf call|validate`; `python -m sd_genotypes.lg5inv.call --panel DIR --out T`; `python -m sd_genotypes.bchrom.score --panel DIR --out T`; `python -m sd_genotypes.te.call --panel DIR --out T [--gsdf-calls T]`; `python -m sd_genotypes.te.definition --ref GT3.fna --out-dir D`; `python -m sd_genotypes.retro --remote R --samples S --modules gsdf,te --out-dir D`; `python -m sd_genotypes.parity --sexloci sv_genotypes.tsv ...`; `sd_genotypes bchrom prepare|call|compare`; `python -m sd_genotypes.lg5inv prepare|status|clean` (pool reads for the LG5 markers); `sd_genotypes sheet-update`.

**The six goals, and how each is measured.** A genotype is made in two separate steps: a *measurement* (raw counts from reads, stored in `<sample>.panel.json` with the hash of its definition) and a *call* (thresholds and flags applied later to the stored counts, so a threshold can change without touching a BAM; measurements made under different definitions are never mixed). Each measurement runs in one or both of two modes. **Disk:** while the sample's BAM is local (03's pipeline calls `panel measure`). **Retro:** a retrospective measurement for a sample whose BAM is no longer local, which reads only the regions the measurement needs, from a stored slice or from the BAM on Dropbox without downloading all of it. **(open: the mechanism of reading regions from Dropbox, and which regions 03 slices; see 10.)**

| # | System | Evidence (the measurement) | Mode | Definition | State |
|---|---|---|---|---|---|
| 1 | gsdf duplication | depth ratio of the 21.7 kb tandem unit to its flanks (5 kb bins, region NC_135173.1:50.70-50.83 Mb) and unique-read (MAPQ >= 30) junction coverage | disk and retro | `config/gsdf.yaml` (calibration from the MC pedigree; minimap2 alignments) | rules aligned with sexloci 0.1.26; tests with fakes; first real-data run (18 samples) reproduced the MC_049 and MC_060 results; **junction cross-check to be recalibrated; parity with sexloci's table pending** |
| 2 | gsdf transposon insertion (chr7ins, directly upstream of gsdf; two reference copies) | read pairs crossing the junctions of a 603 bp junction allele at NC_135173.1:50,765,050 and :50,786,787, summed into one call | disk and retro | `sd_genotypes/te/sites.py` (copied from sexloci, md5 8b7b5cf1d8c7), `te_definition.json` (reference flanks, built once), `config/tepanel.yaml` | measurement, call and tests written (`te`); definition built from GT3a (md5 34f6c2b5dd25); first real-data run (18 samples) reproduced the numbers quoted in VARIANTS.md; **table-level parity with sexloci pending** |
| 3 | two candidate transposon insertions on chr10 (TE1 +1,139 bp at NC_135176.1:24,921,178; TE2 +9,868 bp at :24,921,394) | read pairs crossing the junctions of each full allele | disk and retro | as 2 (md5 1957b81b6617 and d2b75a3b30e2) | as 2 |
| 4 | LG5 (orange blotch / ZW) inversion | ref and alt reads at 954 SNV markers inside the inversion (NC_135171.1:19.68-26.73 Mb) and reads containing the junction sequences at the right breakpoint | disk; retro only if 03 slices the inversion region (19.40-26.98 Mb) | `sd_genotypes/lg5inv/data/`, `config/lg5panel.yaml` | measurement and call written; measurement validated on real pool reads; thresholds provisional |
| 5 | B chromosome | depth over the 874 core blocks (836 tier A) divided by the genome-wide depth, estimated from 12,000 sampled 50 bp bins | **disk only** (blocks and bins are genome-wide) | `sd_genotypes/bchrom/data/`, `config/bchrompanel.yaml` | measurement and call written, tests with fakes; scale estimate still to be checked (see below) |
| 6 | inversions on LG2, 9, 10, 11, 13 and 20 | breakpoint-read evidence at the windows from 01 | disk; retro from the stored breakpoint-window slices | 01's windows | **deferred until the mbuna pilot is complete** |

The pilot is the first mbuna alignments. For the pilot, 4 and 5 and 1 to 3 are measured while the BAM is local; 6 is not started before it is complete.

**Current state, per goal.**
- *1. gsdf duplication* (`gsdf`, `panel`; tests with fakes): raw 5 kb bins and junction coverage stored, copy number called with the calibration of the MC pedigree (0.530 and 0.268). Aligned with sexloci 0.1.26, the code that the calibration and the existing `gsdf number` column came from: a 5 kb bin belongs to a window when it *starts* inside it (unit 7 bins, flanks 8 + 4; an earlier version of 04 took only bins entirely inside); the junction ratio divides by the first flank only; the reported number is the rounded estimate clipped to 1-6 and named CNn outside -/-, dup/-, dup/dup; a call is flagged only when the unrounded estimate is below 0.6 or above 5.4 (the 2-4 range of an earlier draft was an assumption, and wrong: five copies exist in the sheet). The junction cross-check is kept in a separate `check_flags` column and does not change whether a call counts as unflagged. `fit_anchors` and a theory calibration are available for another aligner. The calibration was made on **minimap2 -x sr** alignments (confirmed by the lab; the `bwa` in `VARIANTS.md` is wrong and should be corrected there). In sexloci: MC family 51 samples, all 46 unflagged calls right; agreement with Munby et al. Table S5 98.9% of 629 fish. **First real-data run (18 samples, streamed with `retro`, numpy backend): MC_049 called 2 copies and MC_060 called 4, both flagged, exactly as in the MC family result of VARIANTS.md. Two findings: (a) the junction cross-check assumed a junction ratio of 0.5 and 1.0 for 3 and 4 copies, but a known duplication carrier (AG_563_1) gives 0.18, the 2-copy median is 0.02, and SAMEA117806076 (4 copies by depth, ratio 0.99) gives 0.00; reads wholly inside the 500 bp junction window map ambiguously, so only spanning reads count and the expected ratio is far below 0.5. The check's magnitude is wrong and is to be recalibrated on the MC family (it only affects `check_flags`, never the call). (b) SAMEA117806076 is called 4 copies (EstCopies 3.7) where VARIANTS.md says three: to be compared with sexloci's own value for that sample. (unverified: parity with sexloci's depth numbers over a whole table.)**
- *2 and 3. Transposon insertions* (`te`; tests with fakes): the lab's IndelReadClassifier and the built-in sites of sexloci 0.1.26 are copied, byte-identical and md5-pinned, into `sd_genotypes/te/` (decision: copy, not depend; sexloci stays the reference and keeps the sheet writer and the retrospective runs it already does). For each site, the read pairs in the window anchor +- (250 + 1 + 200) bp (unmapped mates stored at their mate's position included; duplicate, secondary and supplementary excluded) are aligned (Smith-Waterman) to a local reference allele and insertion allele; a pair counts for an allele when a read clearly prefers it (score margin 20) and crosses a junction with 15 bp on each side. The REF and ALT pair counts are stored per site; the call is sexloci's binomial (P(alt) 0.02, 0.5, 0.98; fewer than 5 informative pairs gives ./.; GQ >= 10 counts as confident). The reference flanks are stored in `te_definition.json` (built once from GT3a by `te.definition`, which also checks every anchor base and the pinned fingerprints), so measuring a BAM needs no FASTA. chr7ins is called from the summed pairs of both reference copies; a `dilution_possible` flag marks a 0/0 with insertion pairs in a fish with 3 or more gsdf copies (the copy-number-aware model of VARIANTS.md 5.7 is not written). TE1 and TE2: how the original records were called is NOT RECORDED, the transposon family is unassigned, and phasing is not attempted. With the numpy aligner a read costs about 28 ms against the TE2 allele (parasail, optional, is about 100 times faster). **First real-data run (18 samples, streamed with `retro`, numpy backend) reproduced the numbers quoted in VARIANTS.md: SAMEA117806076 chr7ins 63 reference : 6 insertion pairs `0/0`; AG_567_1 TE1 2 : 5 `0/1` with GQ 13.3; the two Kingiri carriers `1/1` and `0/1`; the eight non-carriers `0/0` or `./.`; YH_004_m, YH_009_m and YH_014_m without any TE1 or TE2 insertion pair; AG_563_1 TE1 `0/1`. (unverified: the dilution flag on real data (the run was made without `--gsdf-calls`), parasail against numpy, parity with sexloci's table.)**
- *Retrospective mode and parity.* `sd_genotypes.retro` streams only the regions the modules need (their `regions(cfg)`: for gsdf and te, the 130 kb gsdf region, which contains both chr7ins windows, and one 1.1 kb TE1/TE2 window) from the BAM on Dropbox through a read-only `rclone serve http` and one `samtools view` per sample, sorts and indexes the small BAM, runs the panel on it and writes the usual `<sample>.panel.json` (resumable; failures in `retro_failures.tsv`; disk-only modules are refused before any network step). Verified on real data: the Dropbox layout is `Bamfiles/Mzebra_GT3_NCBI/<sample>/<sample>.all.bam` with `<sample>.all.bai` (1,337 sample folders; each also holds the chimeric, clipped, discordant, duplication, inversion and unmapped split BAMs, which the panel does not use), and 18 samples streamed without a failure. `sd_genotypes.parity` compares the calls with sexloci's `sv_genotypes.tsv` for the same samples (pair counts, genotypes, depth ratios, copy numbers).
- *4. LG5 inversion.* GT3 carries the standard arrangement; LatrZW and MezeOBm carry the inverted one, which has no standard haplotype in either assembly. At the right breakpoint (GT3 26,726,680) the junction is unique sequence and genotypable (inverted junctions a and b and the standard sequence C); the left breakpoint lies in a repeat and is not. 954 SNV markers (alleles shared by both inverted assemblies, clean in BB pools and Luwino, at about the inverted fraction in the OB pools of both genera) were chosen from five pools of PRJNA808791 aligned with the lab recipe (`lg5inv prepare`; D05). In pools the markers give 0.50-0.51 in the OB pools and 0.001-0.004 in the BB pools (selected on them), 0.53 and 0.57 in the OB pool of the other genus (held out); the held-out BB-type pools are not clean (0.06-0.20: the alt alleles also occur on standard chromosomes of other species), so a diploid fish is called on the median per-site alt fraction. The measurement reproduces the pool counts exactly through a second route. Calls: 0, 1 or 2 inverted chromosomes from the median fraction (cuts 0.25 and 0.75, **provisional**), with junction and aggregate cross-checks as flags. **(unverified: diploid fish; Luwino shows a uniform ~3% component that is either one inverted chromosome in 40 or contamination.)**
- *5. B chromosome* (`bchrom`): the discovery used Clark et al.'s 12 carrier females and 9 B-free pools (21 runs, aligned on D07 with the lab recipe). Core set on GT3: 874 blocks, 1.40 Mb (Clark: 619 blocks, 1.37 Mb), of which 836 (1.36 Mb) are tier A; 825 blocks lie on named LGs. Carriers score 4.6-8.2 (leave-one-out 4.3-8.0), B-free samples 1.08-1.09 (tier A: carriers 4.7-8.4, B-free 1.10-1.11); all three automatic checks pass. The measurement (`bchrom.measure`, `bchrom.score`) stores block depths and the sampled bins; the score is the median scaled coverage over tier A. **(unverified: the sampled-bin scale against the exact genome-wide median on the 21 tracks (`bchrom.validate_scale`); `samtools bedcov` against the mosdepth tracks on a real BAM; thresholds 2.5 and 1.6 are provisional; mbuna males as negatives.)**
- *6. Inversions on LG2, 9, 10, 11, 13 and 20*: breakpoints exist for LG2, 5, 9, 10, 11, 13 (01, run 20261105) and LG20 is not covered. Blumer et al. (Table S7) give chr2, 9, 10, 11, 13 for 1,101 samples; all 440 mbuna among them are 0/0 on all five, so for the mbuna these are validation labels that mostly confirm absence, not an input. LG20 has no labels. The LG5 result is the guide: the first step for each breakpoint is a sequence-only check of how unique it is (one LG5 breakpoint was clean, the other lay in a repeat). The genotyper is not written. **Deferred until the mbuna pilot is complete.**

**Work.**
1. Pilot readiness: 4 and 5 complete (measurement, call, definitions committed, both registered in `panel`); run the scale validation on the 21 tracks; on the first pilot BAMs compare the B measurement with a mosdepth track of the same BAM, run the LG5 measurement, and refit the thresholds on the cohort's own distributions.
2. 1: parity with sexloci on about 40 samples (raw depth ratios and copy numbers), then independent validation (Munby et al. Table S5; the MC family); run `retro` on the existing alignments.
3. 2 and 3: build `te_definition.json` from GT3a (on D05) and commit it; parity of pair counts and genotypes with sexloci on the same samples; chr7ins in the pilot mbuna. Still open (VARIANTS.md section 7): phasing TE1 and TE2, a copy-number-aware chr7ins model, TE family, chr19-ins.
4. Slices are not needed for 1 to 3 (retro streams the BAM). 03 would need to slice the LG5 inversion region only if goal 4 were ever to be measured retrospectively.
5. One `panel` JSON schema for all modules (done: gsdf, te, lg5inv, bchrom). The sheet write-back stays with sexloci (`--write-sample-database`: guarded, cell by cell, RAW values) unless 04 takes it over **(open)**.
6. After the pilot: 6.

**Done when.** 1: parity with sexloci (equal copy numbers on samples unflagged in both; depth ratios equal to rounding); unflagged calls agree with the existing column for at least 98% of the samples and females are 2 and males 3 in the MC family. 2 and 3: equal REF/ALT pair counts and agreeing genotypes with sexloci on the same samples (any difference listed and explained); chr7ins: all 16 of Munby's carriers show the junction and none of his 24 non-carriers, as in sexloci; agreement with the old `TE1`/`TE2` columns where they exist. 4: BB pools about 0 and OB pools about 0.5 inverted on the markers (reached in the pools); on the pilot the median site fraction forms clusters near 0, 0.5 and 1, the thresholds are refit on them, junction and marker evidence agree where the junction is informative, and calls are female-biased in the OB populations. 5: carriers score above the threshold leave-one-out (reached), mbuna males below it; the sampled-scale score within 2% of the exact score on the 21 tracks, and the first pilot BAM's measurement within 2% of its mosdepth score. 6: at least 98% agreement with Blumer on labelled samples, with hets reported separately. Every measurement carries its definition hash.

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

One environment, `enadl` (`envs/download.yml`), covers modules 02 and 04 on Linux: downloads, alignment of the pool reads, the panel (measure, call, retro, parity) and the code copied from sexloci (with parasail as its fast aligner). New machine: `conda env create -f envs/download.yml`; existing: `conda env update -n enadl -f envs/download.yml --prune`; check with `bash tools/check_env.sh`. rclone is not in it: each machine already has one with the `ptm_dropbox` remote (on D05 the snap build keeps its configuration apart from a conda rclone's, which would not see the remote). A known conflict: `aria2` and `sra-tools` cannot be solved in one conda environment on macOS; on Linux they can.

## 8. Order of work

1. **Now, in parallel:** finish 02 (pilots, campaign); validate 03 on real data; finish 04's definitions before any BAM is downloaded again.
2. **Before the mbuna pilot:** 04's goals 4 and 5 complete (LG5 and B-chromosome measurement, call and definitions committed and registered in `panel`) so every pilot BAM is measured once; goals 1 to 3 (gsdf, transposon insertions): measurement and call written and tested with fakes, the definition of `te` still to be built from GT3a, parity with sexloci and a first real-data check still to be done. **The inversions on LG2, 9, 10, 11, 13 and 20 (goal 6) are deferred until the pilot is complete**; 01 finished for LG20 and the open items is needed before them, not before the pilot.
3. **Then:** 05, then 06 and 08 (independent of each other), then 07.
4. **Last:** 09 and 10, which read everything.

Critical path to the mbuna run: 02 complete, 03 validated, the 04 panel defined for goals 4 and 5 (and 1 to 3 when their definitions exist).

## 9. Migration from the old layout (done on branch `reorg_spec`; see MIGRATE.md)

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
2. ~~The transposon insertion definitions and the evidence to use.~~ **Resolved**: sexloci's `known_sites.py` and the junction-pair method, copied into 04 (see 04 and decisions.md). (04)
3. Variant caller, joint or per-family calling, filters, region scope. (05)
4. Species list and sample-size rules for SD discovery. (06)
5. QTL software. (07)
6. Annotation and mappability sources for ASE. (08)
7. Web server hosting, audience and data store. (09)
8. The figure and table list. (10)
9. When to run the migration, and who updates each machine.
10. Whether `SPEC.md` is edited directly in the repo, or changes go through a decision log in `docs/`.
11. ~~Retro mode for 1 to 3: the mechanism.~~ **Resolved**: `sd_genotypes.retro` (stream the needed regions, as sexloci does). Still open: whether 03 should also store slices of the gsdf, transposon and LG5 regions. (03, 04)
12. Final thresholds of the LG5 and B-chromosome calls, set on the pilot cohort's own distributions. (04)
13. Who writes `gsdf number`, `TE1`, `TE2`, `chr7ins` into the SampleDatabase: sexloci (guarded writer, exists) or 04. (04)
