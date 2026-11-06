# 01 Centromeres and inversion locations (LG2, 5, 9, 10, 11, 13, 20)

Specification: [`SPEC.md`](../SPEC.md), module 01. Package: `landmarks` (the number is only for ordering folders).

```
python -m landmarks centromeres arrays | scan | lift --table <xlsx>
    python -m landmarks inversion_breakpoints align | call | figure | all
```

Tests: `python tools/run_tests.py 01` (or `python -m landmarks selftest`). Configuration: `config/`. Environment: `../envs/`.

Current state: centromeres and inversion_breakpoints exist for LG2, 5, 9, 10, 11, 13 (tested on synthetic data and run on the cluster). LG20 is not covered; the LG9 suspected misassembly, the LG13 ambiguity and a second assembly for LG2, 10, 13 are open.
