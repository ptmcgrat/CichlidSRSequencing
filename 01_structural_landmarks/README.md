# 01 Centromeres and inversion locations (LG2, 5, 9, 10, 11, 13, 20)

Specification: [`SPEC.md`](../SPEC.md), module 01 (decision 70). Package: `landmarks` (the number is only for ordering folders).

```
python -m landmarks centromeres arrays | scan | lift --table <xlsx>
python -m landmarks inversion_breakpoints align | call | figure | all          # the breakpoint CALLER: no longer developed, see below
from landmarks.inversions_table import load_inversions, breakpoint_windows       # the breakpoints this project uses
```

Tests: `python tools/run_tests.py 01` (or `python -m landmarks selftest`). Configuration: `config/`. Environment: `../envs/`.

## The inversion breakpoints of this project: the published table (decision 70, 2026-11-11)

The breakpoints are **Kumar et al. 2025 (eLife, Supplementary File 4)** for the six inversions they report: LG2, LG9, LG10, LG11 (as 11a and 11b), LG13 and LG20 (as 20a and 20b). They are recorded in
`data/published_inversions_Mzebra_GT3_NCBI.tsv`: per inversion the chromosome accession, the interval that contains the left breakpoint (`bp1_start`, `bp1_end`) and the one that contains the right breakpoint
(`bp2_start`, `bp2_end`), the reference state (the reference haplotype of GT3a is the ancestral, non-inverted one for all six: Kumar et al.), a status and the source. A tandem (double) inversion is two rows whose shared breakpoint is
the identical interval. Refining the breakpoints further (more assemblies, the LG9 misassembly question, a LG20 resolver, read-level work) was judged too large and too dependent on other references for this project, and belongs to another one.

`landmarks/inversions_table.py` reads the table, **validates it** (accessions match the linkage groups, intervals are ordered, tandem pairs share their breakpoint exactly, every row has a source) and derives the interval middles, the length range and
`breakpoint_windows(inv, pad)`: one window per breakpoint, the shared tandem breakpoint once, for 03's slices and 04's breakpoint-read evidence (the pad is an open decision). It refuses a table with a `pending` row unless asked.
`tests/test_inversions_table.py` holds a frozen copy of the published numbers, so an accidental edit fails.

**LG5 is not in Kumar et al.** (it is the W / orange-blotch inversion), so its row is this pipeline's own call (LatrZW and MezeOBm), taken from the high-confidence consensus row of run `20261110-183552_aa6580ae` (decision 73): `bp1_start` and `bp1_end` are its
`left_min` and `left_max`, `bp2_start` and `bp2_end` its `right_min` and `right_max` (the same columns that reproduce Kumar's LG13 intervals exactly). Both assemblies agree, so the intervals are narrow: **NC_135171.1 19,657,762 to 19,657,848 and 26,726,676 to 26,726,685,
7.07 Mb** (86 bp and 9 bp wide). The other six rows of that consensus table (methods B and C, low or medium confidence, one or two calls) are not used. Behrens et al. 2025 describe a simple inversion of about 7 Mb that LatrZW and MezeOBm share with the same breakpoints
(6.89 Mb in the UMD2a reference, a different reference), which agrees in size and structure; the coordinates are not comparable. `reference_is_standard` stays unset for LG5 (open decision 2 of the spec).

## What our own calls found, for the record (2026-11-10)

`data/comparison_kumar_vs_pipeline_20261110.tsv` compares the calls of `inversion_breakpoints` (run `20261110-183552_aa6580ae`) with Kumar's intervals for 12 breakpoints on LG9, 10, 11, 13 and 20: 3 identical, 5 with Kumar's interval inside ours, 1 with ours inside
Kumar's, 2 overlapping, 1 disjoint (the LG20 internal breakpoint, 20.1 kb). The exact agreements (LG13 both ends, LG10 right) most likely come from shared alignments, so they are consistency and not independent confirmation. The largest disagreement is the
LG9 left breakpoint (ours 14,621,505, 291 kb to the left of Kumar's interval). On LG11 Kumar's internal breakpoint is a 38 bp interval at 18,425,177 inside the 640 kb gap between our two segments. LG2 was not compared (the current run's LG2 table was not at hand).
For LG20 the 20b part is 3.47 to 3.57 Mb from Kumar's own intervals, against 4.1 Mb in their text; our call gave 3.59 to 3.68 Mb. These numbers are not used by the project: the published table is.

## Centromeres and the breakpoint tool

`centromeres` (array scan and the lift of the published table, Chromosome Research 2026, Table 3) is unchanged by this decision; the LG13 centromere ambiguity is open. The `inversion_breakpoint` tool (`align`, `call`, `figure`) stays in the repository, tested on
synthetic data and run on utaka for LG2, 5, 9, 10, 11, 13 (and LG20 on 2026-11-10), but it is **no longer developed**; its only role was LG5's row, which is recorded. The whole-genome alignments it makes are kept only in the `tmp/` folder of the machine that ran it (utaka).
