# Using the 04 panel (gsdf, transposon insertions, LG5 inversion, B chromosome)

What it does: a BAM aligned to Mzebra_GT3_NCBI is reduced to a small JSON of raw counts (a *measurement*), for four systems at once. Genotypes come later, in a separate *call* step that reads a folder of
those JSONs. A threshold can therefore change without touching a BAM, and every JSON carries the definition key of what was measured, so measurements made under different definitions are never mixed.
Design and results: `SPEC.md` (module 04), `docs/decisions.md`, `docs/lg5inv_design.md`.

## 1. Setup
```bash
cd ~/CichlidSRSequencing && git checkout reorg_spec && git pull
conda activate enadl            # or sdpanel (envs/panel.yml) on a machine that only runs the panel: python, numpy, pandas, pyyaml, samtools; parasail optional
source tools/env.sh             # puts every module on PYTHONPATH; run commands from the repo root
bash tools/check_env.sh         # enadl only: tools, packages, the ptm_dropbox remote
```
Retrospective measurements (section 5) also need `rclone` with the `ptm_dropbox:` remote. Everything else needs only samtools and Python.

## 2. Before a long run: the preflight
```bash
python -m sd_genotypes panel check --modules all      # exit code 0 = every module has its config and definition files
```
It prints one line per module with its definition key. Keys at the start of the pilot (2026-11-08): `gsdf 5d71ee8d6c2b`, `te 2af781e90119`, `lg5inv 80bbea57e10b`, `bchrom 75367d0dfba8`. A
different key means a measurement setting or a definition file changed: stop and find out why before measuring more samples. `MISSING` names the command that builds the file.

## 3. Measure one BAM (the BAM is on local disk)
```bash
python -m sd_genotypes panel measure --bam <sample>.all.bam --sample <SampleID> --out-json <SampleID>.panel.json --modules all
```
- Input: a coordinate-sorted BAM aligned to Mzebra_GT3_NCBI, duplicates marked, with its index next to it (`<name>.bam.bai` or `<name>.bai`). A BAM of another reference is refused (contig lengths are checked).
- Output: one JSON (about 154 KB for all four modules), written atomically. It starts with `schema` (2), `sample`, **`bam` (the BAM's file name) and `bam_size_bytes` (its size)**, `panel_version`, `measured` and `samtools`, then one block per module with its definition key. On failure: a traceback, exit code 1, no JSON. The default is `--modules gsdf`; `all` means every registered module; a list such as `gsdf,te` works.
- Cost on a 17x BAM: gsdf 0.3 s, te 0.5 s (with parasail), lg5inv 1.0 s, bchrom 4.3 s; under 80 MB of memory.
- From another program, as a subprocess (preferred: no logic is copied):
```bash
measure_panel() {   # measure_panel <SampleID> <bam> <out_dir>
  python -m sd_genotypes panel measure --bam "$2" --sample "$1" --out-json "$3/$1.panel.json" --modules all || { echo "panel failed for $1" >&2; return 1; }
}
```
- Or from Python: `from sd_genotypes.panel import measure as pm; d = pm.measure(bam, sample, modules='all'); pm.write_json(path, d)`.
- Where the JSON lives is decided by the calling pipeline; module 03's plan is `Bamfiles/<genome>/<sample>/<sample>.panel.json`. Measuring again overwrites it.

### Which BAM was this measured on? (realigned samples)
A realigned sample usually keeps its file name (`<sample>.all.bam`) but not its size, so a stored measurement belongs to a BAM only if BOTH `bam` and `bam_size_bytes` match the file. To ask:
```bash
python -m sd_genotypes.panel.identity --json <sample>.panel.json --bam <sample>.all.bam    # exit 0 CURRENT, 1 STALE (re-measure), 2 UNKNOWN (no size recorded)
```
or in Python `from sd_genotypes.panel import identity as ID; status, reason = ID.same_bam(path_or_dict, bam_path)` with status `'same'`, `'different'` or `'unknown'`. Files written before the identity was stored (schema 1) are still read:
their `bam` is reduced to the file name and the size is taken from the old `bam_size`; those of streamed samples have no usable size (`unknown`). A name and a size are not a checksum: two different alignments of exactly the same size would not be told apart.
Every call table (`gsdf_calls.tsv`, `te_calls.tsv`, `lg5_calls.tsv`, `bchrom_scores.tsv`) carries `bam` and `bam_size_bytes`, so a called genotype can always be traced to the BAM it came from.

## 4. Measure many BAMs that are already on disk
```bash
mkdir -p panel && for b in /path/to/*.all.bam; do s=$(basename $b .all.bam); [ -f panel/$s.panel.json ] || python -m sd_genotypes panel measure --bam $b --sample $s --out-json panel/$s.panel.json --modules all || echo "FAILED $s" >> panel/failed.txt; done
```

## 5. Measure samples whose BAM is not on disk (gsdf and te only)
```bash
python -m sd_genotypes.retro --remote 'ptm_dropbox:/COS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/Bamfiles/Mzebra_GT3_NCBI/' \
    --samples-file ids.txt --modules gsdf,te --out-dir retro_panel --jobs 4
```
`ids.txt` has one SampleID per line (no spaces). The JSON records the identity of the SOURCE BAM on Dropbox (name and size, read first with an HTTP HEAD request: a missing BAM is reported as `no_bam` before anything is fetched), not of the temporary regional BAM that is measured. Only the needed regions are streamed (a read-only `rclone serve http` and one `samtools view` per sample); nothing is written to Dropbox; samples already in
`--out-dir` are skipped whatever BAM they were measured on (after a realignment use `--refresh` or delete the old JSON; `panel.identity` tells which JSONs are stale); failures go to `retro_failures.tsv` (`no_bam`, `no_index`, `error`) and never stop the others. The JSONs are ordinary panel JSONs. `lg5inv` and
`bchrom` are refused: they read megabases or the whole BAM and are measured on disk.

## 6. Make the calls (on a folder of `*.panel.json`)
| System | Command | Output |
|---|---|---|
| gsdf duplication | `python -m sd_genotypes gsdf call --measurements DIR --out-dir OUT` | `OUT/gsdf_calls.tsv`, `sheet_update.tsv`, `summary.txt` |
| transposon insertions (TE1, TE2, chr7ins) | `python -m sd_genotypes.te.call --panel DIR --out te_calls.tsv [--gsdf-calls OUT/gsdf_calls.tsv]` | one TSV |
| LG5 inversion | `python -m sd_genotypes.lg5inv.call --panel DIR --out lg5_calls.tsv` | one TSV and a printed summary |
| B chromosome | `python -m sd_genotypes.bchrom.score --panel DIR --out bchrom_scores.tsv` | one TSV and a printed summary |

Each refuses measurements made under another definition than the current one, and a sample measured twice. Panel files without the module are skipped and counted.

## 7. Reading the calls
- **gsdf** (`gsdf_calls.tsv`): `gsdf_number` is the TOTAL number of copies of the 21.7 kb unit (2 = none duplicated, 3 = one chromosome, 4 = both; 5 and 6 occur); `genotype` names -/-, dup/-, dup/dup, or CNn. Depth
  cannot say WHICH chromosome carries the extra copy: that needs sex or pedigree. `call_flags` (the call is still made): `low_flank_depth` (flank under 3x), `non_integer` (estimate more than 0.35 from an integer),
  `outside_calibrated_range` (below 0.6 or above 5.4 copies), `no_estimate`. Only unflagged calls are written to the sheet. `check_flags`: `junction_in_two_copy` (unique reads across the junction in a sample called with two copies).
- **Transposon insertions** (`te_calls.tsv`): per site `<site>_GT` (0/0 absent, 0/1 heterozygous, 1/1 homozygous, ./. = fewer than 5 informative read pairs), `_ref` and `_alt` (pair counts), `_GQ`, `_status`; `<site>_call`
  holds the genotype only when GQ >= 10 (this is what the sheet gets). chr7ins is called from the pairs of both gsdf copies together (`chr7ins_*`). Only junction pairs inform, so no-calls are common below about 10x. With
  `--gsdf-calls`, `chr7ins_flags = dilution_possible` marks a 0/0 with insertion pairs in a fish with 3 or more gsdf copies (extra reference-like copies dilute the insertion). It cannot say whether TE1 and TE2 are on the same chromosome.
- **LG5 inversion** (`lg5_calls.tsv`): `genotype` std/std, std/inv, inv/inv (`n_inverted` 0, 1, 2) from `median_site_fraction`, the median alt-allele fraction over the 954 markers (cuts at 0.25 and 0.75). std is the
  arrangement of the reference (GT3), inv the OB/W inversion. `junction_A/B/C` and `junction_fraction` are an independent check at the right breakpoint. Flags: `few_markers`, `low_depth`, `near_threshold`,
  `aggregate_disagrees`, `junction_uninformative`, `junction_disagrees`. **The cuts are provisional**: refit them on the pilot cohort's own distribution (it should form clusters near 0, 0.5 and 1).
- **B chromosome** (`bchrom_scores.tsv`): `score` is the median depth over the 836 tier A blocks relative to the genome-wide depth: about 1 without a B, 4-8 with one. `bchrom` B+ (score >= 2.5), B- (<= 1.6), B? in
  between. Flags: `low_depth`, `few_scale_bins`, `few_blocks`, `fraction_disagrees`. **The cuts are provisional**; mbuna males are the expected negatives.

## 8. Where the settings are
| Setting | File | Part of the definition key? |
|---|---|---|
| gsdf region, bins, windows, MAPQ | `config/gsdf.yaml` (region ... mapq_unique) | yes |
| gsdf calibration, flags, labels | `config/gsdf.yaml` (calibration, thresholds, labels) | no |
| te flanks, aligner scores, MAPQ, alleles | `config/tepanel.yaml`, `sd_genotypes/te/data/te_definition.json` | yes |
| te call model | `config/tepanel.yaml` (call) | no |
| LG5 markers, junction probes, region | `config/lg5panel.yaml`, `sd_genotypes/lg5inv/data/` | yes |
| LG5 cuts, flags | `config/lg5panel.yaml` (thresholds) | no |
| B blocks, scale bins | `config/bchrompanel.yaml`, `sd_genotypes/bchrom/data/` | yes |
| B cuts, flags, block set | `config/bchrompanel.yaml` (thresholds, blocks) | no |

Anything in the "yes" rows changes the key and splits the cohort: do not edit those during a run; re-measure instead. Anything in the "no" rows can change at any time (the stored counts do not depend on it).
Rebuilding a definition (only if a "yes" input changes): `python -m sd_genotypes.te.definition --ref <GT3 .fna> --out-dir sd_genotypes/te/data`; `python -m sd_genotypes.lg5inv.definition --ref <GT3 .fna> --markers <markers.tsv> --out-dir sd_genotypes/lg5inv/data`;
`python -m sd_genotypes.bchrom.definition --core-blocks ... --lengths-from ... --out-dir sd_genotypes/bchrom/data` (arguments in `sd_genotypes/*/data/README.md`).

## 9. Checks against known data
- gsdf against the sheet (species and sex tables as well): `python -m sd_genotypes gsdf validate --table sheet.tsv --calls OUT/gsdf_calls.tsv --out-dir V`.
- gsdf and TE against sexloci's table: `python -m sd_genotypes.parity --sexloci sv_genotypes.tsv --te-calls te_calls.tsv --gsdf-calls OUT/gsdf_calls.tsv --out-dir P` (pair counts, genotypes, depth ratios).
- B scale estimate against the 21 discovery tracks: `python -m sd_genotypes.bchrom.validate_scale --tracks DIR --summary sample_summary.tsv`; against mosdepth on one BAM (`mosdepth -n --by 50 PREFIX bam`):
  `python -m sd_genotypes.bchrom.validate_scale --panel S.panel.json --mosdepth PREFIX.regions.bed.gz` (expect a score ratio within 0.02 of 1).

## 10. When something goes wrong
| Message | Meaning and what to do |
|---|---|
| `... has no index` | the BAM has no `.bai`/`.csi` next to it: `samtools index` it |
| `... is N bp in the header, the project reference has M` / `no contig ...` | the BAM was not aligned to Mzebra_GT3_NCBI |
| `MISSING` in `panel check`, or `... is missing: build the definition with ...` | the definition file is not in the repo; run the command in the message |
| `measurements made under definition(s) [...] the current definition is ...` | the JSONs were made with another config or definition: re-measure them, or call them separately |
| `sample measured more than once` | two JSONs for one sample in the folder |
| retro: `no_bam` / `no_index` | no `<sample>/<sample>.all.bam` (or index) in the Bamfiles folder |
| retro: `rclone serve exited immediately` | the `--remote` path is wrong: test it with `rclone lsf '<remote>' --dirs-only` |
| retro: `sample names must not contain whitespace` | the id file has spaces or tabs pasted as spaces |
| te slow (seconds per sample at 20x) | parasail is not installed; install it (`pip install parasail`); the JSON records the backend used |
| gsdf `low_flank_depth` or lots of `./.` | the sample has too little coverage; the call is reported but should not be used |

## 11. What it does not do
It does not write to the SampleDatabase (sexloci's guarded writer does, for `gsdf number`, `TE1`, `TE2`, `chr7ins`); it does not phase TE1 and TE2; it has no copy-number-aware chr7ins model (only the flag); it does not
genotype the inversions of LG2, 9, 10, 11, 13 and 20 (deferred until the mbuna pilot is complete); it does not fit thresholds (that is done on the pilot cohort's distributions).
