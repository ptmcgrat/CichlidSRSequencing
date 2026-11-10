# 01 Centromeres and inversion locations (LG2, 5, 9, 10, 11, 13, 20)

Specification: [`SPEC.md`](../SPEC.md), module 01. Package: `landmarks` (the number is only for ordering folders).

```
python -m landmarks centromeres arrays | scan | lift --table <xlsx>
    python -m landmarks inversion_breakpoints align | call | figure | all
```

Tests: `python tools/run_tests.py 01` (or `python -m landmarks selftest`). Configuration: `config/`. Environment: `../envs/`.

Current state: centromeres and inversion_breakpoints exist for LG2, 5, 9, 10, 11, 13 (tested on synthetic data and run on the cluster). **LG20 is configured (below) but has not been run on the real alignments.**
The LG9 suspected misassembly, the LG13 ambiguity and a second assembly for LG2, 10, 13 are open.

## LG20: the double inversion 20a + 20b (and LG11 as 11a + 11b)

**What is known** (Kumar et al. 2025, eLife). Chromosome 20 carries two tandem inversions: 20a (the left arm) formed first, 20b (the right one, 4.1 Mb) after it, in that background. The benthics, Diplotaxodon and utaka carry both,
*Rhamphochromis* only 20a, the mbuna and *A. calliptera* neither, so GT3 (*M. zebra*) is the standard orientation (`reference_is_standard: true`). The Darwin Tree of Life assembly of *Rhamphochromis* sp. 'chilingali' (**fRhaChi2.1,
GCA_963969265.1**; the accession is from NCBI's assembly hub page and is verified by the md5 check of the fetch) is the one that carries 20a alone: Kumar report that its alignment to GT3a shows a single inversion on chr20. **Not settled:**
Blumer et al. 2025 (supplement) call chr20 an inversion differentially fixed between Diplotaxodon and *all* benthics, which sits uneasily with benthics carrying both parts. The benthic assemblies (YH, fAulStu2) are in the run to find out and are deliberately not in `expect`.

**How it is identified, the same way as LG11.** `align` and `call` give, per assembly, the strand runs and the inversion with its **segments**: a tandem inversion A B' C' D, each part inverted in place, is one inverted run along which the reference coordinate
jumps back up, so it is called as ONE inversion with TWO segments and an internal breakpoint (`paf.inversion_segments`). What LG20 adds is an assembly that carries only the first part: it is an ordinary single inversion whose right breakpoint is that
internal breakpoint, and the consensus (which needs both breakpoints to agree) lists it as a different inversion. `parts.py` puts them together, for every chromosome that declares `parts: [a, b]` (names left to right in reference coordinates) in
`config/inversion_breakpoints.yaml`: the structure comes from the assemblies that resolve two segments; a single-part carrier is recognised by one outer breakpoint plus the internal one and sharpens the internal breakpoint (otherwise it is only the gap
between the two segments, hundreds of kb wide); every assembly gets a state; `expect` is compared with the result.

```
python -m cichlid_common genomes fetch --only fRhaChi2          # downloads and md5-verifies the assembly (the registry is common/config/genomes.yaml)
python -m landmarks inversion_breakpoints align                  # the whole-genome alignments of the other assemblies are cached; only fRhaChi2 is aligned (asm5); the targets now include LG20
python -m landmarks inversion_breakpoints call
```

New files next to the existing ones in `<master>/YH_SD_Paper/<genome>/inversion_breakpoints/call/<run>/` (also uploaded by the runner): `LG20_parts.tsv` (each part: breakpoints with their intervals, length, carriers, the evidence),
`LG20_assembly_states.tsv` (each assembly: both / `20a only` / `20b only` / both with the parts not resolved / another inversion / none, with its expectation), `LG20_parts_breakpoint_regions.bed` (the four breakpoint intervals), and the polarity JSON
and `summary.json` carry the same. The call step prints the table and **DIFFERS FROM THE LITERATURE** when an assembly does not match `expect`. Read the dotplots (`LG20_<assembly>_dotplot.png`) and `LG20_strand_runs.tsv` before trusting a result.
LG11 gets the same files (`LG11_parts.tsv`, as 11a + 11b): all three of its assemblies are expected to carry both parts, a regression check on the real data.

If fRhaChi2 is `20b only`, or a part is placed where the literature does not put it, check first that part 1 is really the left part in GT3 coordinates (the names follow the reference, Kumar's left and right follow their figures).

**What is tested and what is not.** `tests/test_double_inversion_parts.py` (19 tests) builds synthetic alignments like the real LG11 case: two-segment carriers (also reverse-complemented), a first-part-only carrier, a second-part-only carrier, a collinear one, an
inversion with other breakpoints, no double carrier at all, and the same data through the real `call_step`. **Nothing here has been run on the real LG20 alignments**: the breakpoints, the carriers and the verdicts are for the first real run to say.
The time of the fRhaChi2 alignment is not measured.
