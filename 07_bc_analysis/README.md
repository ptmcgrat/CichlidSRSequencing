# 07 Backcross genotyping, recombination, TD, QTL mapping, fine mapping

Specification: [`SPEC.md`](../SPEC.md), module 07. Package: `bc_analysis` (the number is only for ordering folders).

```
python -m bc_analysis bc_genotyping provenance|flagstat|original|v2-sites|v2-counts|v2-hmm|validate|compare
```

Tests: `python tools/run_tests.py 07` (or `python -m bc_analysis selftest`). Configuration: `config/`. Environment: `../envs/`.

Current state: bc_genotyping (about 266 fish, 9,336 markers) is the part to fix; the analyses after it are not written.
