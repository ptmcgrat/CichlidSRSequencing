# align_pipeline

One sample at a time, from uBAMs to a verified upload and a ledger row:

    download -> align (+merge) -> mark duplicates -> [chimeric split | statistics | slices | panel | count uBAM reads] -> qc -> verified upload -> ledger row -> cleanup

It is v6 (`alignFastQReads_v2.py`: pipelined, space-aware, resumable, no GVCFs) in the module layout, with the decisions of 2026-11-08 applied (docs/decisions.md 40 to 43). **It has not run on real data yet**:
its 128 tests (v6's 76, adapted, three of them removed with the old whole-sheet writer, and 55 new) all pass offline, with fake tools and a handful of real subprocesses running stand-in programs.

## Use

    source tools/env.sh                                                      # puts every module on PYTHONPATH (the panel subprocess inherits it); run from the repository root
    python -m align_pipeline setup-database Mzebra_GT3_NCBI                  # ONCE, from one machine: adds the ledger's columns to the worksheet header
    python -m align_pipeline run Mzebra_GT3_NCBI -e Mbuna --list             # what would run, and how much data
    python -m align_pipeline run Mzebra_GT3_NCBI -s SAMPLE --print-commands SAMPLE
    python -m align_pipeline run Mzebra_GT3_NCBI -e Mbuna --panel-expect 5d71ee8d6c2b 2af781e90119 80bbea57e10b 75367d0dfba8

`run` takes the selection options of the old script (`-s -c -p -e -u`), `--include-excluded`, `--batch-size`, the resource options (`--max-cores --max-ram-gb --samples-in-flight --align-threads --uploads ...`),
`--keep-local`, `--no-upload` (test mode: nothing uploaded, no row, nothing deleted; the genome folder is still synchronised, as in the old pipeline), `--no-slices`, `--no-panel`, `--panel-expect KEY...` (the run does not
start unless `panel check` prints these definition keys: the pilot's frozen keys). `-r` selects samples that already have a ledger row (v6's meaning; the rerun rule that replaces it is not built yet).
`python -m align_pipeline run --help` lists all of them. The exit code is 1 when any sample failed. Two machines: give each a disjoint selection.

## What a finished sample leaves in `Bamfiles/<genome>/<sample>/` (all of it is uploaded with the folder)

`<sample>.all.bam` and `.all.bai`; `.all.bam.duplication_metrics.txt`, `.all.bam.stats.txt`, `.all.bam.idxstats.txt` (the Picard and samtools output the statistics are computed from); `<sample>.chimeric.bam` and `.bai`
with `<sample>.split_summary.json`; `slices/<sample>.slice.LG2.bam`, `.LG10`, `.LG13` and their `.bai`; `<sample>.qc.json` (every statistic, about 20); `<sample>.panel.json` (the panel's raw measurements, which
record the BAM they were measured on: calls are made offline from them). The uBAMs, the temp folder and the scratch files are deleted after the verified upload.

## What changed from v6

| | v6 | here |
|---|---|---|
| splits | six BAMs in one pass | **only the chimeric BAM** (the panel uses no other); `route()` still knows all six |
| slices | whole LG10 | **LG2, LG10, LG13**, whole chromosomes (`slice_regions/<Genome>.tsv`, `END` = the length from the reference index); checked against the index before the run starts |
| panel | none | **a stage**: `panel check` once before anything is aligned (the run does not start if it fails), `panel measure` per sample after duplicate marking; **a failure never fails the sample** |
| ledger row | the old columns plus about 20 statistics | **22 columns** (`alignment_db.LEDGER_COLUMNS`); everything else stays in `<sample>.qc.json` |
| time limits | none (a hung `rclone` would hang the sample, as one did for 88 days) | uploads must make 5 MB/s (`--upload-min-mb-s`, at least 30 minutes), the quick commands get an hour, the panel half an hour |
| dependencies | `download_queue.py` (and psutil) of the old layout | none: `procs.py` holds `bam_check`, `kill_tree`, `configure_rclone_env`; psutil is optional |
| writer | `DatabaseWriter` (whole sheet) and the appender | the appender only |
| entry point | a top-level script | `python -m align_pipeline run|setup-database` |

## The ledger row (the sheet)

`SampleID, GenomeVersion, Species, RunIDs, NRuns, InputFingerprint, Coverage, UsableCoverage, TotalReads, MappedFraction, DuplicationRate, BamSize, Instrument, ReadLength, QCWarnings, PanelStatus, AlignedOn,
Machine, PipelineVersion, minimap2_version, gatk_version, pysam_version`. `DuplicationRate` is Picard's PERCENT_DUPLICATION pooled over libraries (not the old `DuplicationReads`, the `.duplication.bam` split);
`UsableCoverage` is the measured non-duplicate coverage; `RunIDs` is joined with `,,` as in the existing rows; `InputFingerprint` hashes the sorted (RunID, uBAM size) pairs (the input of the rerun rule, which compares it with the runs now in the DNAReadsDatabase: not built yet);
`PanelStatus` is `ok`, `failed` or `off`; `PipelineVersion` is `align_pipeline <version> <commit>` (the commit is read from the .git folder when there is no git program); `QCWarnings` is empty when nothing is wrong. Old rows keep their old columns. **Which columns the sheet carries is an open decision (SPEC.md 12.13): it is one constant.**

## The panel stage

Exactly the contract of SPEC.md module 04: `python -m sd_genotypes panel check --modules all` and `python -m sd_genotypes panel measure --bam <sample>.all.bam --sample <ID> --out-json <ID>.panel.json --modules all`,
as subprocesses of `--panel-python` (default: this interpreter). After every measurement the pipeline checks that the JSON names this BAM (`bam` and `bam_size_bytes`, schema 2, decision 45; `panel.identity_status` follows
`python -m sd_genotypes.panel.identity`): another BAM is a failed panel and the JSON is removed; no identity (an old `sd_genotypes`) is accepted with a QC warning. A failed panel never fails the sample: the BAM is uploaded,
the row says `PanelStatus` = failed with the reason, and module 04's own commands (`panel measure`, `retro`) measure it later.

## Before the first real run

Read `DEPLOY_AND_VALIDATE.md`: how to put this build in the repository from the laptop, check yellowhead's environment, and validate on one small sample against the old BAMs on Dropbox.

## Not built yet

* **The rerun rule** (decided 2026-11-09): a sample whose ledger row does not hold all the runs now in the DNAReadsDatabase (compared as sets, and by `InputFingerprint`) is realigned from ALL its runs together; its
  BAMs are replaced, its row is replaced in place (never appended), and after the verified upload the stale GVCF and the five old splits are deleted from its Dropbox folder (logged). Until then `-r` is v6's.
* A local journal for when Google is unreachable (v6 raises after four attempts and leaves an uploaded sample without a row, which `audit_alignments.py` and `register_alignments.py` repair).
* Not planned: claims (two machines get disjoint selections) and a genotype-only mode (module 04's commands do that).

## Tests

`python tools/run_tests.py 03` from the repository root (all offline). `test_real_exec.py` runs real subprocesses (stand-in `gatk`, `minimap2`, `samtools`, `rclone` and a stand-in `python -m sd_genotypes`).
