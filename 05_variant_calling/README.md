# 05 VCFs for families and species

Specification: [`SPEC.md`](../SPEC.md), module 05. Package: `variant_calling` (the number is only for ordering folders).

```
python -m variant_calling callset build <name> [--dry-run]
```

Tests: `python tools/run_tests.py 05` (or `python -m variant_calling selftest`). Configuration: `config/`. Environment: `../envs/`.

Current state: the callset builder of the BC work (bc_family, yh_broods, mc_broods); caller choice and cohort definitions are open.
