# 03 Align, mark duplicates, split, slice, QC, upload

Specification: [`SPEC.md`](../SPEC.md), module 03. Package: `align_pipeline` (the number is only for ordering folders).

```
python alignFastQReads_v2.py <genome> ...    (still in cichlid_sr_sequencing/: see ../MIGRATE.md)
```

Tests: `python tools/run_tests.py 03` (or `python -m align_pipeline selftest`). Configuration: `config/`. Environment: `../envs/`.

Current state: lab_repo_patch_v6 (tests with fakes); real-data validation (validate_against_old.py) not done; not moved yet because an alignment run is using the old file manager.
