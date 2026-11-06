# 06 SD regions in many species

Specification: [`SPEC.md`](../SPEC.md), module 06. Package: `sd_regions` (the number is only for ordering folders).

```
python -m sd_regions parental_sex_region callset|depth|scan|summarize [--species YH,MC]
    python -m sd_regions sex_association fisher [--table T --out D]
```

Tests: `python tools/run_tests.py 06` (or `python -m sd_regions selftest`). Configuration: `config/`. Environment: `../envs/`.

Current state: parental_sex_region and sex_association exist (tested on synthetic data); not run on real data.
