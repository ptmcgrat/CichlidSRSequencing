# Detection and validation scripts

These stand-alone scripts did the detection and validation behind `docs/VARIANTS.md`. They are **not part of the installed package**, are not covered by the test suite beyond a check that each file parses and matches the checksum below,
and were written for the lab laptop: they expect to be run from `~/sexloci` in the environment where `sexloci` works, reuse its Dropbox streaming (needs `samtools` and `rclone` with the lab's `ptm_dropbox:` remote) and write under `~/Temp/SexDeterminationAnalysis`. Each script's docstring has its options.

| Script | Variant | What it does | Main input | Main output | md5 (first 12) |
|---|---|---|---|---|---|
| `selfdot.py` | gsdf | Self dot plot of one genome region; reports the repeat period and both copies of a tandem duplication. | a GT3a FASTA, a chromosome and a window | `<out>.png`, `<out>_segments.tsv`, a printed summary | `4b5e46b48c7e` |
| `gsdf_cn.py` | gsdf | The original stand-alone copy-number genotyper (depth over the unit relative to its flanks). Superseded by `sexloci sv-genotype`, which uses the same constants. | a depth table and optional sex table | `gsdf_cn.tsv` and a plot | `157c0d5ea315` |
| `gsdf_dup_check.py` | gsdf | Is the tandem duplication present in a genome assembly? Builds a probe from GT3a, then tests assemblies. | GT3a FASTA; assemblies | probe FASTA; a presence/absence report per assembly | `94ed0761d80d` |
| `gsdf_primer_tools.py` | gsdf | Designs three primer assays for the duplication (junction PCR and qPCR) and checks them in silico on any assembly. **Not tested at the bench.** | GT3a FASTA or one fish's BAM | primer tables (see `data/gsdf_primers_MC050.tsv`) | `b3f85f304720` |
| `clip_check.py` | TE1, TE2 | Counts soft-clipped reads at TE1/TE2 in the region BAMs `sv-genotype` keeps; separates a clean reference, the expected insertion, and a different insertion at the same site. | region BAMs from `sv-genotype` | a table per sample and site | `3a263bef7141` |
| `cov_check.py` | all | Small coverage QC helper (no docstring). |  |  | `811b564a516e` |
| `alignment_check.py` | chr7ins | Which carriers have alignments, and against which genome (alignment sheet and Dropbox listing, read-only). | a carriers table | a report | `4961e6688e73` |
| `find_insertion.py` | chr7ins | Locates an insertion absent from the reference by carrier-versus-control soft-clipped junction reads, and assembles the two ends. | carrier and control BAM regions streamed from the cloud | `candidates.tsv`, `report.txt`, end FASTAs, `anomalous_mate_peaks.tsv` (in `insertion_scan/`) | `97a5938f6ac9` |
| `mates_vs_repeatmasker.py` | chr7ins | Identifies the transposon family from the RepeatMasker annotation at the positions of the anomalous reads' mates. | region BAMs; the GT3a RepeatMasker `.out` | an enrichment table of repeat families, carriers versus controls | `11c5446bde73` |
| `element_from_annotation.py` | chr7ins | Reconstructs the element from intact genomic copies: pairs the two end families, lists the body families, extracts and ranks copies. | RepeatMasker `.out`; GT3a FASTA | `element_template.fasta`, `element_modal_copies.fasta`, `element_candidates.fasta`, `element_best_per_class.fasta` | `d595c728f868` |
| `template_check_with_mates.py` | chr7ins | Aligns the carriers' mates to a template, ranks the modal copies, reports where the carriers differ and writes the polished estimate. | carriers' BAMs; the template and modal copies | `element_final_estimate.fasta` (3,841 bp) | `b3c15cedb1da` |
| `build_chr7ins_alleles.py` | chr7ins | Builds the genotyping alleles (anchor + first 300 bp + last 300 bp + TG) at both reference copies. The result is built into sexloci as set `chr7_gsdf_ins`; rerun only if the element estimate changes. | `element_final_estimate.fasta`; GT3a FASTA | `chr7ins_alleles.vcf` | `0c4415a066ab` |
| `chr7ins_combine.py` | chr7ins | Sums the informative pairs of the two copies and calls one genotype per sample. Now done inside `sv-genotype` (`chr7ins_*` columns). | `sv_genotypes.tsv` | `chr7ins_combined.tsv` | `71e50b0b9c5a` |
| `chr7ins_concordance.py` | chr7ins | Compares the combined (or candidate) genotypes with Munby's calls for 47 fish: AUC, per-group tables, genotype agreement. | genotypes; `data/chr7ins_validation_calls.tsv` | `chr7ins_concordance.tsv` and a printed report | `7638908b570f` |

## Order for chr7ins

1. `alignment_check.py` (optional): which carriers have GT3a alignments.
2. `find_insertion.py`: carrier-versus-control junction scan; candidate #1 is the insertion.
3. `mates_vs_repeatmasker.py`: the family of the element, from the mates.
4. `element_from_annotation.py`: intact copies and a template.
5. `template_check_with_mates.py`: check and polish the template with the carriers' reads; `element_final_estimate.fasta`.
6. `build_chr7ins_alleles.py`: the two genotyping sites (now built in as `chr7_gsdf_ins`).
7. `sexloci sv-genotype`: genotype every sample (the combination of the two copies is now built in).
8. `chr7ins_concordance.py`: validate against Munby's calls.

## Order for the gsdf duplication

1. `selfdot.py` on the GT3a region, then `minimap2 -X -c` of the region against itself for exact boundaries.
2. `sexloci sv-genotype` (or the older `gsdf_cn.py`) for copy number; `gsdf_dup_check.py` for assemblies.
3. `gsdf_primer_tools.py` for the bench assays (untested).

## Data

* `data/chr7ins_validation_calls.tsv`: Munby's chr7-ins genotypes, read counts, location and link basis for 47 fish (input to `chr7ins_concordance.py`).
* `data/chr7ins_validation_set.tsv`: the 28-fish validation set (11 carriers, 17 non-carriers) with ENA run accessions.
* `data/gsdf_primers_MC050.tsv`: primers designed for the gsdf duplication from fish MC_050.
