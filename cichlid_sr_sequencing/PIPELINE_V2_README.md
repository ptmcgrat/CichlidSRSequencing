# Pipelined alignment pipeline (v2): `alignFastQReads_v2.py`

New files only; the existing `alignFastQReads.py` and `alignment_worker.py` are untouched and remain the fallback.
Requires the download patch (`helper_modules/download_queue.py`, lab_repo_patch v3/v4), `samtools` >= 1.10, `pysam`, `gatk`, `minimap2`, `rclone` in the conda environment.
This package (lab_repo_patch_v6) replaces v5: apply v6 instead of v5, not on top of it.

## What changed and why

| Old (`alignFastQReads.py`) | New (`alignFastQReads_v2.py`) |
|---|---|
| every stage is a barrier over ALL samples; all uBAMs downloaded first | each sample moves on its own: download -> align -> mark duplicates -> split / stats / slice -> qc -> verified upload -> database row -> cleanup |
| alignment: one sample at a time with `minimap2 -t 96` | `--samples-in-flight` samples at once, `--align-threads` (default 16) each |
| `max_processes` = number of cores, as if every job used one core | one shared core + RAM budget; every command declares what it needs and waits for it; `nice`; **later stages are served before new alignments** (finished work drains, uploads start early) |
| downloads for the whole batch at once, no limit on disk | **continuous, space-aware downloads** (below): never ahead of the disk; integrity checks run in the background; a continuous queue keeps every download slot busy |
| uploads only at the very end, one at a time | **a separate upload pool**: a finished sample is uploaded while other samples align, mark and split; a slow upload never holds a compute slot; at most `--max-pending-uploads` finished samples wait locally |
| nothing written to the database until the whole batch finished; a late failure repeated everything | database row as soon as the sample is verified on Dropbox; stage checkpoints, so a rerun resumes; a failure stops only that sample |
| only the LAST process of the alignment pipe was checked | `bash -o pipefail`: a crash of SamToFastq, minimap2 or MergeBamAlignment fails the sample |
| whole-sheet rewrite of the database for every sample (two machines lose rows) | **atomic append of one row** with read-back verification (below) |
| ~64 python processes per sample + 6 `MergeSamFiles` runs to split the BAM | one pass over the BAM writes the six split BAMs (decision logic identical to the old script, tested on 60,000 random reads) |
| HaplotypeCaller GVCFs (about 40 core-hours per sample) | none (known sites are genotyped from the BAMs or the slice later) |
| `pdb.set_trace()` on errors; uBAMs deleted right after alignment | no `pdb`; **uBAMs deleted as soon as the duplicate-marked BAM is verified and the read count is done**; local files deleted after the upload is verified (`--keep-local` to keep) |
| 7 x `gatk CountReads` per sample | `samtools stats` (one pass) + `samtools idxstats` (from the index, seconds); the uBAM read count runs alongside the alignment |
| no region slices | `slice_regions/<Genome>.tsv`; default: **all of chromosome 10** (`NC_135176.1`, about 3.6% of each BAM) |

Unchanged on purpose: the SamToFastq | minimap2 `-x sr` | MergeBamAlignment recipe and its options, MarkDuplicates and its options, all file names, the six split-BAM definitions.
Changed JVM settings only: explicit `-Xmx` and `-XX:ParallelGCThreads=2`.

## Disk safety: the download queue cannot run ahead of the disk

Two gates must both be open before ANY file of a sample is downloaded:
1. **Disk slots:** at most `--samples-in-flight + --prefetch` samples have uBAMs on disk at once.
2. **The disk ledger:** every admitted sample reserves `--disk-reserve-factor` (3.5) x its uBAM size (the uBAMs, sorted BAM, final BAM, scratch and split files at the peak). A new sample is admitted only if
   `free disk now - (what samples already in flight still have to write) >= its own reservation + --min-free-gb (200)`.

If there is no room the queue simply waits; a finished and cleaned-up sample frees space. A reservation shrinks to the sample's measured size when its large writes are over, and disappears when its local files are deleted or it fails.
If nothing is in flight that could ever free space (the disk is simply too full), the next sample is skipped with `not enough free disk ...` after 30 minutes instead of waiting for ever.
The uBAM sizes come from the sample database (`FileSize`); 20 GB per file is assumed if a size is missing.
The summary and the `STATUS` lines report the **peak local disk use actually observed** as a multiple of the uBAM size: if it ever approaches 3.5x, raise `--disk-reserve-factor`.
The 200 GB margin absorbs estimation error: in the tests, samples that outgrew their reservation by 14% still never filled the disk.

## The database (two machines): atomic append

Every finished sample is APPENDED as one row with a single atomic Google Sheets call, then the SampleID column is read back: the sample must be there exactly once (a lost row is appended again; a duplicate is reported).
Two machines can therefore finish samples at the same moment without losing rows. The old pipeline rewrote the whole sheet for every sample, which does lose rows when two writes overlap (reproduced in `tests/test_alignment_db.py`).
**Once, from ONE machine, before the first real run:**
```bash
python3 alignFastQReads_v2.py Mzebra_GT3_NCBI --setup-database       # adds the new statistics columns to the worksheet header; idempotent
```
A real run checks the header first and stops before any work if columns are missing. Do not let the OLD pipeline's database step (the BC run's final rewrite of the sheet) overlap with these runs.

## The statistics (new columns; the old ones keep their old meaning)

The old `UnmappedReads` and `DuplicationReads` are NOT what their names suggest:
* `UnmappedReads` = records in `unmapped.bam` / all records. `unmapped.bam` holds only unmapped reads whose MATE is mapped (not duplicates, mean base quality > 25). Pairs with both reads unmapped are never in it.
* `DuplicationReads` = records in `duplication.bam` / all records: read pairs that are evidence of a tandem DUPLICATION (orientation signal). It is not the PCR duplicate rate; that was never recorded.
* all split-BAM fractions divide by every record, duplicates included, although the split BAMs exclude duplicates.

| new column | meaning |
|---|---|
| `PrimaryReads` | reads sequenced (R1 + R2), excluding secondary and supplementary alignments |
| `MappedFraction`, `UnmappedFraction` | mapped / unmapped reads over `PrimaryReads` |
| `UnmappedNoPositionFraction` | both mates unmapped (no position in the BAM; not in `unmapped.bam`) |
| `UnmappedOrphanFraction` | unmapped, placed at the mate's position (the pool `unmapped.bam` draws from). The two add up to `UnmappedFraction` |
| `DuplicateFraction` | PCR/optical duplicate rate by Picard's definition, pooled over libraries: (unpaired dups + 2 x pair dups) / (unpaired examined + 2 x pairs examined) |
| `DuplicateFractionFromBam` | the same from the duplicate flags in the BAM (`samtools stats`); the two should agree |
| `OpticalDuplicateFraction` | optical duplicates (Picard) over reads examined |
| `MeanReadLength`, `MappedBases`, `CoverageMapped`, `CoverageNonDuplicate` | measured read length; bases mapped (CIGAR); mean depth with and without duplicates |
| `ProperlyPairedFraction`, `MQ0Fraction`, `InsertSizeMean`, `InsertSizeSD`, `ErrorRate` | from `samtools stats` |
| `UBamReads` | reads counted in the input uBAMs; should equal `PrimaryReads` |
| `QCWarnings` | filled when independent measurements disagree (duplicate rate by more than 1 point, fully-unmapped reads, uBAM vs BAM read count, ...) |
| `PipelineVersion` | `v2` |

## Options that matter

| option | default | meaning |
|---|---|---|
| `-e`, `-u`, `-p`, `-c`, `-s` | | select samples by Ecogroup, Subgroup, Project, Species or ID (as in the old script). Samples with Ecogroup `Error` / `Exclude` are skipped unless `--include-excluded` |
| `--max-cores` | all but 8 | the only thing that matters for sharing a machine |
| `--samples-in-flight` | 6 | samples computed at the same time |
| `--align-threads` | 16 | minimap2 threads per sample |
| `-d`, `--MaxDownloads` | 4 | files downloaded at the same time (try 8: a single stream was ~10 MB/s) |
| `--uploads` | 2 | samples uploaded at the same time |
| `--prefetch` | 2 | samples whose uBAMs may be downloaded beyond those in flight |
| `--min-free-gb`, `--disk-reserve-factor` | 200, 3.5 | the disk gate above |
| `--markdup-heap-gb` | 16 | MarkDuplicates heap; try 32 on the test run |
| `--tmp-dir` | sample temp folder | base folder for the GATK/Picard scratch files (a fast local disk if you have one) |
| `--status-minutes` | 10 | status line: what every sample is doing, transfers, free disk, cores |
| `--batch-size N` | all | only the first N selected samples |
| `--list`, `--print-commands SAMPLE` | | preview the selection (with the uBAM GB) or the commands; nothing runs |
| `--no-upload` | | test mode: everything local, no upload, no database row, nothing deleted |

## Install

```bash
# laptop, in the repository root
git apply --check lab_pipeline_v2.patch && git apply lab_pipeline_v2.patch
cd cichlid_sr_sequencing && python3 -m unittest discover -s tests          # all tests, no external tools needed
git add -A && git commit -m "pipelined alignment pipeline v2 (v6)" && git push
# the machines that run it: git pull
```

## First use: validate against the old pipeline on one test sample

The test samples are already in the database, so use `-r`; `--no-upload` writes no database row, uploads nothing and deletes nothing (and needs no `--setup-database`).
```bash
G=Mzebra_GT3_NCBI; B=/Data/$USER/Temp/CichlidSequencingData/Bamfiles/$G
cp -r $B/SAMEA117806260 $B/SAMEA117806260.old                          # a new run writes into the same folder name
python3 alignFastQReads_v2.py $G -s SAMEA117806260 -r --no-upload --samples-in-flight 1 | tee ~/v2_test.log
python3 unit_scripts/validate_against_old.py $B/SAMEA117806260.old $B/SAMEA117806260 SAMEA117806260
cat $B/SAMEA117806260/SAMEA117806260.qc.json
```
Expect IDENTICAL for all seven BAMs. Then inspect `qc.json`: `QCWarnings` should be empty. The slice is `$B/SAMEA117806260/slices/SAMEA117806260.slice.LG10.bam`.

## The mbuna run (two machines)

```bash
tmux new -s align
# once, from one machine:   python3 alignFastQReads_v2.py Mzebra_GT3_NCBI --setup-database
python3 alignFastQReads_v2.py Mzebra_GT3_NCBI -e Mbuna -u Maylandia Labeotropheus --batch-size 8 --max-cores 88 -d 8 2>&1 | tee ~/v2_batch1.log     # machine 1; machine 2: other subgroups
```
Rerun the same command to continue: finished samples are skipped (database row), unfinished ones resume at their last completed stage. Ctrl-C kills every command it started; finished stages are kept.
Run the two machines on DIFFERENT subgroups; two machines on the same samples would align them twice.
Log lines show every stage per sample, the per-sample download and upload rates, and periodic `STATUS` lines. The summary lists failures with the stage and the tail of the error log.
uBAMs of a failed sample stay on disk (a rerun does not download them again); delete them by hand if the failure is final.

## Things I could not verify (no GATK, samtools, pysam, rclone or Google access where this was written)
* The commands are those of the old pipeline plus `samtools` for statistics, slices and indexes; syntax was checked by reading, not by running. `--print-commands` shows them for review.
* The Google Sheets calls (`row_values`, `col_values`, `append_rows`, `add_cols`, `update_cell`) are written against the gspread API and tested against a fake worksheet only.
* That `MergeBamAlignment` keeps every read of the uBAM, so that `PrimaryReads` equals the uBAM read count. If `QCWarnings` reports a difference on every sample, this assumption is wrong; it is a warning, not a failure.
* `samtools stats` field names and the Picard metrics columns were taken from their documented format.
* The default speed settings (6 samples in flight, 16 alignment threads, 4 downloads, 2 uploads, 3.5x disk reservation) are guesses: the summary prints per-stage times, transfer rates, peak cores and the peak disk use actually observed.
* Reads without a position (both mates unmapped) are, as before, not in any split BAM; they stay in `all.bam` and are counted in `UnmappedNoPositionFraction`.
* Kept from the old split script on purpose (worth a look): the adapter filter `not read.has_tag('XM')` is probably inert under minimap2; only the clipped-read rule excludes MAPQ 0.
