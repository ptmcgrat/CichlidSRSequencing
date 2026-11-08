# 04 gsdf duplication, transposon insertions, inversions, B chromosome

Specification: [`SPEC.md`](../SPEC.md), module 04. Package: `sd_genotypes` (the number is only for ordering folders). **How to use it: [`docs/USAGE.md`](docs/USAGE.md).**

```
python -m sd_genotypes panel check --modules all                                  # preflight: definition files and keys
python -m sd_genotypes panel measure --bam B --sample S --out-json J --modules all    # one BAM on disk -> raw counts (gsdf, te, lg5inv, bchrom)
python -m sd_genotypes.retro --remote R --samples-file F --modules gsdf,te --out-dir D  # samples without a local BAM (stream the regions)
python -m sd_genotypes gsdf call|validate | python -m sd_genotypes.te.call | .lg5inv.call | .bchrom.score    # the calls, on a folder of panel JSONs
python -m sd_genotypes.parity --sexloci sv_genotypes.tsv ...                          # compare with sexloci's table
python -m sd_genotypes bchrom prepare|call|all|status|clean|compare                   # discovery of the B chromosome
python -m sd_genotypes.lg5inv prepare|status|clean                                    # pool reads for the LG5 markers
```

Tests: `python tools/run_tests.py 04` (or `python -m sd_genotypes selftest`). Configuration: `config/`. Environments: `../envs/` (`enadl` for everything on Linux, `sdpanel` to run the panel only).

Current state (2026-11-08): measurement, call and definitions for gsdf, the transposon insertions (TE1, TE2, chr7ins), the LG5 inversion and the B chromosome are committed and registered in `panel`; `panel check` passes. Validated on real data:
gsdf (MC family 46 of 46 against the sheet), chr7ins (Munby's 47 fish, AUC 1.000), LG5 (pools of known composition), B (scale estimate on the 21 discovery tracks, `bedcov` against mosdepth on a real BAM). Thresholds of the LG5 and B calls are
provisional until the mbuna pilot. The inversions of LG2, 9, 10, 11, 13 and 20 are deferred until the pilot is complete.
