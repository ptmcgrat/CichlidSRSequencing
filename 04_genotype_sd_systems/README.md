# 04 gsdf duplication, transposon insertions, inversions, B chromosome

Specification: [`SPEC.md`](../SPEC.md), module 04. Package: `sd_genotypes` (the number is only for ordering folders).

```
python -m sd_genotypes panel measure --bam B --sample S --out-json J
    python -m sd_genotypes gsdf call|validate ...
    python -m sd_genotypes bchrom prepare|call|all|status|clean|compare
    python -m sd_genotypes sample_table inversions --blumer-dir D
```

Tests: `python tools/run_tests.py 04` (or `python -m sd_genotypes selftest`). Configuration: `config/`. Environment: `../envs/`.

Current state: gsdf and panel (tested with fakes, not yet on real BAMs); bchrom downloads through the SRA source with ENA as the fallback (download.source in config/bchrom.yaml) and has not run on real data; LG5 genotyper, transposon genotyper and the B score are not written.
