# 02 ENA / SRA download into uBAMs on Dropbox, tracked in the DownloadDatabase sheet

Specification: [`SPEC.md`](../SPEC.md), module 02. Package: `download_reads` (the number is only for ordering folders).

```
python -m download_reads run --manifest Blumer Malinsky [--source ena|sra|auto] [--max-runs N] [--no-db] [--dry-run]
    python -m download_reads manifest import|status|reconcile|release|tabs
    python -m download_reads bench | patch-file-manager
```

Tests: `python tools/run_tests.py 02` (or `python -m download_reads selftest`). Configuration: `config/`. Environment: `../envs/`.

Current state: v9d: two sources (ENA, SRA on AWS) with per-run checks, claims and cell-level writes to the sheet; one run completed end to end; pilots on utaka and mzebra and the campaign are outstanding.
