# lab_repo_patch_v9o_ena: downloadShortReadData_v2 / grab_ena_v2 + the download database + the SRA source (supersedes v7a-v8b)

New files only (nothing of the old scripts is touched; `downloadShortReadData.py` and `grabENA.py` keep working, and its `--Local` mode stays there).

    files/cichlid_sr_sequencing/helper_modules/sra_fetch.py        the SRA source: aria2c from NCBI's copy on AWS -> vdb-validate -> fasterq-dump -> checks
    files/cichlid_sr_sequencing/downloadManifest.py                import / status / release for the DownloadDatabase spreadsheet
    files/cichlid_sr_sequencing/helper_modules/download_manifest.py   the manifest: tabs, claims, cell-level batched writes, journal (Google or in memory)
    files/cichlid_sr_sequencing/unit_scripts/patch_file_manager.py    adds downloadManifest / downloadTabs / readDownloadTab + self.dl_ID to file_manager.py (diff preview, backup, idempotent)
    files/cichlid_sr_sequencing/downloadShortReadData_v2.py        orchestrator
    files/cichlid_sr_sequencing/helper_modules/ena_fetch.py        ENA lookups, file selection, md5-verified downloads (ascp + https), FASTQ cleaner
    files/cichlid_sr_sequencing/helper_modules/ena_scheduler.py    two-pool scheduler with an adaptive download limit
    files/cichlid_sr_sequencing/helper_modules/grab_ena_v2.py      one run: fetch stage and convert stage
    files/cichlid_sr_sequencing/unit_scripts/ena_bench.py          measure how many parallel downloads help on your network
    files/cichlid_sr_sequencing/tests/                             131 tests (incl. your MalinskyPaperData.csv)

## Install
    unzip lab_repo_patch_v7_ena_*.zip -d /tmp && cd ~/CichlidSRSequencing/cichlid_sr_sequencing        # your clone
    rsync -a --ignore-existing /tmp/lab_repo_patch_7_ena/files/cichlid_sr_sequencing/ ./              # only adds files
    python3 -m unittest discover -s tests                                                              # 131 tests, ~1 s

## Use
    python3 downloadShortReadData_v2.py runs.csv --dry-run                         # classify every run, write the plan, download nothing
    python3 downloadShortReadData_v2.py --samples-file blumer_missing_from_sheet.tsv --dry-run
    python3 downloadShortReadData_v2.py runs.csv --dl-workers 3 --workers 4 --ascp-rate 300m
    python3 downloadShortReadData_v2.py --retry download_status_<stamp>.tsv        # retry what failed in a retryable way
Run list from a Run Selector / ENA CSV (no `Library Name` column needed: ENA's, else the run id), from ENA accessions (`--ena PRJEB20078`), or from a sample list (`--samples-file`, Biosample [+ Bioproject] column).
Outputs in the temp folder: `ena_plan_<stamp>.tsv` (every run and its class: ok / in_database / no_fastq / single_end / missing_mate / bad_listing / not_on_ena),
`download_status_<stamp>.tsv` (every run, state, stage, failure kind, attempts, error), per-run `<run>.fetch.log`, `<run>.convert.log`. A run goes into the sample database only after its uBAM is seen on the cloud.

## The SRA source (when ENA is slow)
`--source ena | sra | auto` (default auto). With auto the 20 s Aspera probe decides: ENA first if Aspera moves data, SRA first otherwise; the other source is the fallback for a run that fails (attempts rotate between the sources).
SRA route per run: `aria2c` (8 connections, resumable) from `https://sra-pub-run-odp.s3.amazonaws.com/sra/<run>/<run>` -> `vdb-validate` -> `fasterq-dump --split-3` -> `<run>_1.fastq` / `<run>_2.fastq` into the SAME conversion stage
(FastqToSam, adapter marking, upload) -> .sra, FASTQ and scratch deleted. Needs `conda install -c conda-forge -c bioconda sra-tools aria2`; the downloader checks for the tools before it starts. Peak disk is reserved at 8x the ENA FASTQ size
(.sra + uncompressed FASTQ + scratch + uBAM) whenever SRA is among the sources. The run list and metadata still come from ENA's portal API.
Every SRA run is CHECKED: (1) vdb-validate; (2) fasterq-dump's `spots read` must agree with ENA's `read_count`, which has NO single meaning (measured: ERR13389718 has 109,512,218 = 2 x 54,756,109 spots, i.e. reads; ERR14734229 has 40,613,306
against 40,650,256 spots, i.e. about pairs): both readings are tried and the closer is used; fewer spots than ENA by more than 0.1 % is rejected, up to 5 % more is accepted and written to `ReadCountCheck` (measured on 11 runs: +0.04 % to +1.7 %, never fewer) (with the number of unpaired spots
fasterq-dump put in a separate file, which the pipeline does not use); (3) the first 1000 reads are compared with ENA's FASTQ (sequence and quality): identical -> `identical-to-ENA`; different -> the run is rejected (`source_unavailable`: the next source is tried), because NCBI also
distributes SRA Lite files with simplified qualities. If ENA does not answer in time the check is recorded as `skipped (ENA unreachable)` (`--sra-require-check` refuses such runs instead), and a run with only one or two distinct quality
symbols is rejected as Lite-like. Binned qualities alone are NOT a sign of Lite: the original HiSeq 2500 data of ERR13389718 have 6 distinct symbols at ENA and at NCBI alike. The pilot showed identical first reads (sequence and quality string).
The sheet gets `Source` (SRA-aws), `QualityCheck`, `ReadsInUBAM` (the pairs fasterq-dump wrote) and `MBPerSec`. `--max-runs N` processes only the N smallest pending runs: use it for the first real test. `--sra-connections`, `--fq-threads`.

## Retry of every Google call, and the recheck of unverified runs (v9f)
`helper_modules/gspread_retry.py` wraps the one method every gspread call goes through and retries connection resets, timeouts and HTTP 429 / 500 / 502 / 503 / 504 (6 tries, about 5-60 s apart, a line printed for each retry). A request that
may already have reached Google is NOT repeated if it is not idempotent (an append would add a row twice). It is installed by file_manager.py at import: `python3 unit_scripts/patch_file_manager.py --dry-run`, then without `--dry-run` (it also
upgrades the download-database methods; a backup is written; idempotent). Motivation: an alignment run on utaka died in `addAlignmentRow` with `Connection reset by peer`, because FileManager retried the write of the alignment sheet but not its read.
`unit_scripts/verify_ubam.py --skipped [--write]` compares the first reads of finished uBAMs (streamed from Dropbox) with ENA's FASTQ for the runs whose SRA quality check was skipped because ENA did not answer.

## Adding the pending rows in bulk (v9g)
`--add-pending FILE` adds run by run: every run re-reads both sheets and rewrites both whole sheets, which for a few hundred runs exceeds Google's 60 reads per minute (HTTP 429) and takes hours. `--add-pending FILE --bulk` makes the same
rows (the rules of FileManager.addDNAReadRow: a run is added only if it is not in DNAReadsDatabase and its uBAM is on the cloud; a sample only if it is not in SampleDatabase) with ONE read and ONE write of the two sheets: both sheets are
saved to CSV in the temp folder first, nothing is written if any line of the pending file is malformed, the grid of a worksheet is enlarged if it has too few rows (set_with_dataframe does not), and the result is read back and checked.
Existence and size of every uBAM come from ONE `rclone lsjson` per project folder (FileManager.returnCloudFiles is an uncached `rclone lsf` of the whole folder and returnFileSize an `rclone size`, one of each per run: that, not Google, made the
per-run version take hours); a dry run asks for no sizes. Only when nothing else writes the two sheets meanwhile. `--dry-run` (with or without --bulk) reads and reports, and never writes (v9h still reconciled the DownloadDatabase tab in a dry run). Like every write of FileManager it rewrites the whole sheet from values, so cell formulas become values;
Google Sheets' version history has the old sheet.

## Auditing the alignments, and a minimal patch for an older file manager (v9j)
`unit_scripts/audit_alignments.py --genome G` is READ-ONLY: it compares the AlignmentDatabase rows of one genome with `Bamfiles/<genome>/` on Dropbox (one recursive rclone lsjson) and on this machine, and classifies every sample: ok, needs_row (uploaded,
row never added), partial_upload (some files, an empty file, or a size that differs from the local copy), needs_upload (complete here only), local_incomplete, row_only (a row but files missing on Dropbox), not_started. The usual files are not assumed:
they are the endings that at least 80% of the finished samples have. Written for Python 3.7 (utaka). `unit_scripts/patch_file_manager.py --retry-only` installs only the retry of Google calls, for a file manager that predates the download
database (the one of the alignment pipeline on utaka).

## Registering alignments that are already uploaded (v9k)
`unit_scripts/register_alignments.py` adds the AlignmentDatabase rows of samples whose BAMs are on Dropbox but whose row is missing (audit class `needs_row`: what is left after a crash in alignment_worker.uploadAndUpdateDatabase, which uploads
every sample first and only then counts and adds the rows one by one). It builds each row as that loop does (gatk CountReads of the seven BAMs, coverage = reads * ReadLength / 2 / reference size, fractions, RunIDs joined with ',,', tool versions, BamSize),
counts in parallel (`--jobs`), and adds ALL rows with one read and one write of the sheet after a CSV backup, with the grid enlarged if needed and a read-back check. `--counter gatk|samtools|pysam`: gatk is what alignment_worker used (one thread per BAM, tens of minutes for an all.bam); samtools / pysam run `samtools view -c` with several threads (many times faster, the same count for well-formed BAMs; pysam uses the samtools built into pysam, nothing to install);
each BAM is its own task and a progress line is printed per sample. `--check-existing N` recounts N registered samples and compares EVERY computed column with the sheet (the counts, the fractions, Coverage and RunIDs)
(validates the counting and the formulas); `--limit` and `--dry-run` try a few without writing. It needs the BAMs on this machine; samples without them are reported and skipped, and nothing already in the sheet is changed. A sample with reads but no row in the SampleDatabase is skipped, as alignment_worker skips it (its merged table has no ReadLength for it); `--allow-unlisted` registers it anyway with the ReadLength of the DNAReadsDatabase.

## The SampleDatabase: checks and the Blumer update (v9o)
`unit_scripts/sampledb.py check [--s1 Table_S1.tsv --s7 Table_S7.tsv]` is READ-ONLY: duplicates, orphans between SampleDatabase / DNAReadsDatabase / DownloadDatabase, codings, blank key fields, taxa written in several ways, and how the sheet compares with Blumer.
`unit_scripts/sampledb.py blumer --s1 ... --s7 ... [--write]` adds Blumer et al.'s information to the samples whose SampleID is a Blumer BioSample: Species (Blumer's name overwrites; the SRA name is kept in the new column SRASpecies, taken from the DownloadDatabase
for EVERY sample that has runs there), Sex and Ecogroup (filled if blank, a differing value is kept and listed), LabReared 0, BlumerSequenceID, BlumerSpecimenID, Inversion2/9/10/11/13. Without --write it only writes reports (every change, cell by cell, the
species pairs, unmatched samples). With --write: CSV backup, grid enlarged, ONLY the SampleDatabase tab written once, read back cell by cell. Measured on the export of 1,820 rows: 1,101 samples match S1, Species differs for 324 (103 pairs), Sex and Inversion10 agree for all,
Ecogroup differs for 4 (Error x2, Exclude, one Utaka: manual curation, kept).

## The download database (Google Sheet `DownloadDatabase`, one TAB per Run Info file)
Instead of CSV files, the runs and what happened to them live in a spreadsheet of their own (`FileManager.dl_ID` = 1gNFAd_B8tWxQZRJeEmuhVl1qsqJf1KJtKYpehoEMfrQ), shared as Editor with the service account. Tabs: Blumer, Malinsky, ...
Columns: the 14 Run Info columns, then `DownloadStatus` (pending / claimed / downloading / converting / uploaded / in_database / failed / skipped), `InDatabase`, `DownloadedOn` (UTC), `Source` (the site: ENA-aspera, ENA-https,
legacy for runs that were already there; SRA-aws when that source exists), `Machine`, `ClaimedBy`, `ClaimedAt`, `UpdatedAt`, `Attempts` (added up over invocations), `FailureKind` (permanent / transient / throttle), `LastError`,
`Verified` (md5 matched), `ReadsSource` (ENA's read count), `ReadsInUBAM` (reserved: filled by the SRA route's read-count check, empty for ENA runs), `UBamBytes`, `UBamPath`, `DownloadSeconds`, `MBPerSec`, `PipelineVersion`,
`Notes` (yours), `Skip` (TRUE = never download this run).

    python3 unit_scripts/patch_file_manager.py --dry-run              # once: shows the 17 lines it adds to helper_modules/file_manager.py; then without --dry-run (a backup is written)
    python3 downloadManifest.py import Blumer BlumerPaperData.csv     # creates / extends the tab; runs already in DNAReadsDatabase start as in_database, uBAMs on Dropbox without a row as uploaded
    python3 downloadManifest.py import Malinsky MalinskyPaperData.csv
    python3 downloadManifest.py status
    python3 downloadShortReadData_v2.py --manifest Blumer Malinsky --dry-run
    python3 downloadShortReadData_v2.py --manifest Blumer Malinsky --dl-workers 3 --workers 4 [--no-db]
    python3 downloadManifest.py reconcile Blumer Malinsky             # after --add-pending: runs now in the database become in_database (also: --add-pending FILE --manifest Blumer does it)
    python3 downloadManifest.py release Blumer --hours 6              # after a crash: give stale claims back

How it writes (NOT like setSampleDatabase, which rewrites whole worksheets and flattens formulas): only the CELLS of changed rows, in one batched call every 15 s; the row of a run is found by RunID at every flush (sorting or
inserting rows in the sheet is harmless); a run is CLAIMED (re-read to confirm) before it is downloaded so two machines do not take the same run, claims expire after `--claim-hours` (6); if Google cannot be reached the updates go to a
local journal (`download_manifest_<tab>.journal.jsonl` in the temp folder) and are replayed at the next flush, and claims are skipped so that an outage never stops a download. Re-importing never touches existing rows, so your Notes, Skip
flags and manual status edits survive. `--retry-failed no|retryable|all` (default retryable). A dry run only reads the sheet.
Verified here: all of the logic, against an in-memory sheet and a fake gspread client (131 tests). NOT verified: the real Google API (the backend uses worksheet.get_all_values / col_values / row_values / batch_update /
append_rows / add_worksheet / add_cols / freeze; if a call differs in your gspread version, `import` fails with the traceback and it is a one-line fix).

## When Aspera moves no data (UDP blocked)
Symptom: the benchmark shows ~0 MB/s at every setting, the free disk does not shrink, no session dies. ascp's TCP control connection (port 33001) works but its UDP data channel is dropped by a firewall, so the session
stays open and silent. Defences built in: (1) before the downloads start a 20 s PROBE decides Aspera vs https (`--no-probe` skips it, `--probe-seconds`); (2) an Aspera transfer that receives no data for `--stall-seconds`
(90) is killed and the run falls back to https; (3) after two such stalls the remaining runs go straight to https. https is md5-verified like Aspera. To measure https: `python3 unit_scripts/ena_bench.py <runs> --https`.
Permanent fix: ask IT to allow outbound UDP 33001 from the machine to fasp.sra.ebi.ac.uk (193.62.193.164).

## Run Info files: what is downloaded
    python3 downloadShortReadData_v2.py blumer.csv malinsky.csv --dry-run          # any number of files, which must be in .../Reads/ReadDownloadFiles/ on Dropbox (as for the old script)
    python3 downloadShortReadData_v2.py blumer.csv malinsky.csv --dl-workers 3 --workers 4
For every run the plan says: `in_database` (a row in DNAReadsDatabase: skipped), `on_cloud` (its uBAM is already on Dropbox but it has no row: not downloaded again; `--adopt-cloud` adds the row),
`ok` (to download), or why it cannot be (no_fastq / single_end / missing_mate / bad_listing / not_on_ena). Printed: the counts, `to download: N runs, X GB of FASTQ`, and per SAMPLE how many are complete / partial / none in
the database (`ena_plan_samples_<stamp>.tsv`: a sample with 9 runs is complete only when all 9 are in). The check is per run, as in the old script. Each ENA accession is asked once, however many files list it.
Put the Run Info files BEFORE options that take several values (`--ena`).

## Two machines, or a running job whose file manager can overwrite the sample database: `--no-db`
Downloading and uploading uBAMs does not touch the sample database; only the database ROW does. With `--no-db` a run is downloaded, converted, uploaded and verified on the cloud exactly as before, but its row is appended to
`pending_db_rows.jsonl` in the temp folder instead of being written. When nothing else writes the database, add the rows from ONE machine:

    python3 downloadShortReadData_v2.py --samples-file blumer_missing_from_sheet.tsv --no-db
    python3 downloadShortReadData_v2.py --add-pending ~/Temp/<...>/pending_db_rows.jsonl      # only rows whose uBAM is on the cloud and whose run is not in the database yet; the rest is reported

## Setting up another machine (read from helper_modules/file_manager.py)
* The local folder is `$HOME/Temp/CichlidSequencingData/` for every host that is not utaka / yellowhead / burtoni: nothing to edit. Temp files go to `.../Temp/`.
* Credentials: `FileManager()` downloads `cichlidsrsequencing_api_creds.json` itself through rclone, so a working `ptm_dropbox:` remote is all that is needed.
* Python packages: pandas, openpyxl, gspread, gspread-dataframe, google-auth (the sample database is a Google Sheet), plus gatk4 / samtools / curl for the conversion.
* `addDNAReadRow` re-reads both sheets (SampleDatabase, DNAReadsDatabase) and rewrites them COMPLETELY: do not run another writer of that spreadsheet at the same time (a manual edit, `sample_table --write-sheet`). The alignment database is a different spreadsheet.
* If the database write fails for any reason (including the assertion that the run is already in the database) the run keeps its row in `pending_db_rows.jsonl` and is reported as `done (row pending)`; add them with `--add-pending`.
* Upload failures (`rclone check` raises) are retried, not treated as permanent.

## What was fixed (line numbers of your uploads)
| old | problem | now |
|---|---|---|
| grabENA 113, 125 | `sys.exit()` = status 0 after a failed conversion: the parent recorded the run | exit 10 / 11 / 12 and a status JSON; the parent adds a row only after `done` AND the cloud listing shows the uBAM |
| dSRD 108-109 | mate files taken by position; ENA lists an unpaired file first for some runs | chosen by name (`_1` / `_2`), orphans ignored and noted |
| dSRD 56, 99 | `layout`, `data.RunID` undefined: one single-end / no-FTP run aborts the batch | classified before anything is downloaded; recorded, never fatal |
| dSRD 89-94 | bare retries without a pause; the third unprotected | back-off, and a failure is a recorded status |
| grabENA 83 | `zip()` of the two mates hides a truncated file; no md5 | ENA md5 + size verified; unequal read counts raise |
| grabENA 135, 137 | cleanup only on success | `finally`, after failures too |
| dSRD 144, 123 | the next run's errors written to the wrong file; logs per sample (50 samples of your CSV have 3-9 runs) | logs per run and stage |
| dSRD 125 | busy loop | sleeps |
| grabENA 118-130 | non-Illumina platforms fail at upload | the uBAM is used as it is |
| grabENA 72-99 | Python rewrite of both FASTQs for every run (slow, ~3x the disk) | only as a fallback when FastqToSam rejects the files |
| dSRD 44 | `'Instrument:'` | fixed |

## Speed
**What limits it.** One Aspera (FASP) transfer is built to fill your link by itself; the target rate `-l` is a ceiling, not a share. N transfers at `-l 1000m` ask the network and ENA for N Gbit/s, so beyond the link capacity they only
compete: UDP loss, retransmissions, aborted sessions, and whatever per-IP session limit ENA applies (I do not know ENA's current limits; the error texts in your logs will tell). The old workers also did the CPU-heavy conversion, so keeping the CPUs busy meant more simultaneous downloads.

**What changed.**
* Two pools. `--dl-workers` (default 3) downloads; `--workers` (default 4) conversions, which start as soon as a run's FASTQs are on disk, so the CPUs are busy while the link is busy. `--max-ready` bounds the runs waiting for a conversion slot (disk).
* The download limit adapts (AIMD): a refused / aborted / timed-out session halves it (not below 1) and cools down 60 s; three clean downloads add one back (not above `--dl-workers`).
* Starts are staggered (`--stagger 8` s), failed runs wait (`--backoff`) and retry (`--max-attempts 3`), and the scheduler prints a status line a minute: downloads, ready, converting, waiting, done, failed, average MB/s, free disk.
* ascp: `-k 1` (a retry RESUMES a partial file; the old call restarted from zero), both mates in ONE session (one handshake; `--no-one-session` to undo), the paths from ENA's `fastq_aspera` field, `-l` per transfer from `--ascp-rate`, https with md5 as the fallback when a session is refused or UDP is blocked.
* Disk: a run is admitted only if its peak (4 x its FASTQ size, `--disk-factor`) fits next to what the runs in flight still need, keeping `--min-free-gb`.

**Find your numbers (10 minutes).** Pick 4+ big runs and run
    python3 unit_scripts/ena_bench.py ERR715541 ERR14733391 ERR1749458 ERR13389714 --concurrency 1 2 3 4 --seconds 60 --rate 300m
It prints the aggregate MB/s for 1..4 parallel transfers and how many sessions died; use the smallest setting within 10 % of the best as `--dl-workers`, and `--ascp-rate` of about link capacity / that number.
Repeat at another time of day: ENA's load varies.

## What I could not verify (no network, no ascp, no FileManager here)
* `fastq_aspera` as the ENA field name and its path format (if it is missing the path is derived from the ftp path as the old script did).
* The ascp error patterns in `classify_ascp_error` are my reading of typical ascp messages. If a throttled session is classed 'transient', or a missing file 'throttle', adjust the two regexes in `ena_fetch.py`; the logs show the real text.
* That ascp accepts two sources in one call for ENA, and that `-k 1` resumes against ENA as it does elsewhere. If either fails the https fallback runs; `--no-one-session` removes the first.
* Your FileManager: I assumed `returnCloudFiles(dir)` returns file names, `uploadData(path)` uploads and `returnFileSize(path)` works after the local copy is deleted (as in your script).
* Not included: `--Local`; downloading a CRAM / BAM when ENA has no FASTQ (such runs are listed as `no_fastq` with the submitted file name).
