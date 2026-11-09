# align_pipeline

One sample at a time, from uBAMs to a verified upload and a ledger row:

    download -> align (+merge) -> mark duplicates -> [chimeric split | statistics | slices | panel | count uBAM reads] -> qc -> verified upload -> ledger row -> cleanup

It is v6 (`alignFastQReads_v2.py`: pipelined, space-aware, resumable, no GVCFs) in the module layout, with the decisions of 2026-11-08 applied (docs/decisions.md 40 to 43). **It has not run on real data yet**:
its 128 tests (v6's 76, adapted, three of them removed with the old whole-sheet writer, and 55 new) all pass offline, with fake tools and a handful of real subprocesses running stand-in programs.

## Use

    cd 03_align_pipeline
    export PYTHONPATH=$PWD:$PWD/../common                                  # or install the packages
    python -m align_pipeline setup-database Mzebra_GT3_NCBI                  # ONCE, from one machine: adds the ledger's columns to the worksheet header
    python -m align_pipeline run Mzebra_GT3_NCBI -e Mbuna --list             # what would run, and how much data
    python -m align_pipeline run Mzebra_GT3_NCBI -s SAMPLE --print-commands SAMPLE
    python -m align_pipeline run Mzebra_GT3_NCBI -e Mbuna --panel-pythonpath ../04_genotype_sd_systems

`run` takes the selection options of the old script (`-s -c -p -e -u`), `--include-excluded`, `--batch-size`, the resource options (`--max-cores --max-ram-gb --samples-in-flight --align-threads --uploads ...`),
`--keep-local`, `--no-upload` (test mode: nothing uploaded, no row, nothing deleted), `--no-slices`, `--no-panel`. `python -m align_pipeline run --help` lists all of them. The exit code is 1 when any sample failed.

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
`UsableCoverage` is the measured non-duplicate coverage; `RunIDs` is joined with `,,` as in the existing rows; `InputFingerprint` hashes the sorted (RunID, uBAM size) pairs (the input of the rerun rule, not built yet);
`PanelStatus` is `ok`, `failed` or `off`; `QCWarnings` is empty when nothing is wrong. Old rows keep their old columns. **Which columns the sheet carries is an open decision (SPEC.md 12.13): it is one constant.**

## The panel stage

Exactly the contract of SPEC.md module 04: `python -m sd_genotypes panel check --modules all` and `python -m sd_genotypes panel measure --bam <sample>.all.bam --sample <ID> --out-json <ID>.panel.json --modules all`,
as subprocesses of `--panel-python` (default: this interpreter). `<sample>.panel.meta.json` next to the JSON records which BAM (name, size) and which runs (fingerprint) it was measured from, so that a re-alignment
makes the old measurement detectable as stale (`panel.sidecar_is_current`); module 04's JSON does not record this yet (requested).

## Before the first real run

1. `--list` and `--print-commands` on a sample: read the commands.
2. **Validate against the old pipeline on one sample.** A run writes into the sample's own folder, the same one the old pipeline wrote to: **copy the old folder first**, then
   `python -m align_pipeline run Mzebra_GT3_NCBI -s SAMPLE --no-upload --keep-local --no-panel`, then
   `python validate_against_old.py <copy of the old folder> <the new folder> SAMPLE` (compares `all` and `chimeric`, record by record). They must be identical.
3. Then with the panel (`--panel-pythonpath`), still `--no-upload`; then a small real batch.

## Not built yet

The rerun rule (a sample whose runs changed is realigned; the fingerprint and the sidecar are in place), updating a row in place, a local journal for when Google is unreachable, claims so two machines do not align
the same sample, `genotype-only`, and the supervision of what the pipeline does not time-limit.

## Tests

`python tools/run_tests.py 03` from the repository root (all offline). `test_real_exec.py` runs real subprocesses (stand-in `gatk`, `minimap2`, `samtools`, `rclone` and a stand-in `python -m sd_genotypes`).
