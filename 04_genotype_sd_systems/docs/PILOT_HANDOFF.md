# Module 04 panel, PILOT hand-off for the one-pass writer (AlignmentDatabase and SampleDatabase)

Status: **pilot**. Every value below is provisional and will change as the definitions and thresholds improve; the call settings carry `status: pilot` and a `version`, and every flags cell written carries the tag
`pilot-v<version>`. Measurements and calls are separate: calls are recomputed from the stored measurements, never from a BAM. Longer guide: `USAGE.md`; design and results: `SPEC.md` (module 04), `decisions.md`.

## 1. Measuring (one local, indexed BAM aligned to Mzebra_GT3_NCBI)
```python
from sd_genotypes.panel import measure as pm
d = pm.measure(bam, sample, modules=('gsdf', 'te', 'lg5inv', 'bchrom'))     # or modules='all'; no Project, Dropbox or FileManager
pm.write_json(f'{sample}.panel.json', d)                                     # atomic; or: python -m sd_genotypes panel measure --bam B --sample S --out-json J --modules all
```
- Registered in `panel/measure.py` MODULES: `gsdf`, `te` (the three transposon sites: the one upstream of gsdf = `chr7ins`, `TE1`, `TE2`: one module because they share one classifier and one definition), `lg5inv`, `bchrom`.
- Needs only samtools and Python with numpy, pandas and pyyaml (parasail optional: about 100x faster transposon alignment). Run `python -m sd_genotypes panel check --modules all` first: it prints the definition key of every module and exits 1 if a definition file is missing.
- The JSON (schema 2) holds `bam` (file name) and `bam_size_bytes` (to tell a measurement of an old BAM from one of a realigned BAM), `measured`, `panel_version`, and one block per module with its **definition key**.
  The key is a hash of the measurement settings and the CONTENT of the definition files (not of the whole config): changing a threshold never changes it, changing what is measured always does. Keys at the start of the pilot:
  gsdf `5d71ee8d6c2b`, te `2af781e90119`, lg5inv `80bbea57e10b`, bchrom `75367d0dfba8`.

## 2. Calling (pure: stored measurement + config, no BAM)
```python
from sd_genotypes.panel import calls
r = calls.call_panel(d)                       # or a path; r['calls'] is a list of calls (te gives three); also r['bam'], r['bam_size_bytes']
values, flags = calls.sheet_row(r)            # {column: value} for the CONFIDENT calls, {flags column: 'pilot-v1;flag;flag'} always
```
A call: `column`, `flags_column`, `value` (what the data say; None if there is no estimate), `flags` (list), `confident`, `value_to_write` (value if confident else None), `call_tag` (`pilot-v1`), `definition`, `details`.
`confident` = there is a value and none of the module's `blocking_flags` is set. Thresholds, `status`, `version`, the sheet columns and the blocking flags are in `config/{gsdf,tepanel,lg5panel,bchrompanel}.yaml`; change them and recompute.
A measurement made under another definition than the repository's current one is refused (`DefinitionMismatch`): use the commit that holds its definition, or measure again.

## 3. Sheet columns and the flags convention (proposal)
| Value column | Flags column | Value written | Exists in the sheet |
|---|---|---|---|
| `gsdf number` | `gsdf flags` | integer number of copies (2 = none duplicated) | yes (910 values, written by sexloci) |
| `chr7ins` | `chr7ins flags` | `0/0`, `0/1`, `1/1` (insertion directly upstream of gsdf, both reference copies together) | **no: new** |
| `TE1` | `TE1 flags` | `0/0`, `0/1`, `1/1` | value column yes (824), flags new |
| `TE2` | `TE2 flags` | `0/0`, `0/1`, `1/1` | value column yes (813), flags new |
| `Inversion5` | `Inversion5 flags` | 0, 1 or 2 inverted (OB/W) chromosomes, like `Inversion10`; 0 = the arrangement of the reference | **no: new** |
| `B score` | `B flags` | float: median scaled depth over the 836 tier A blocks (about 1 without a B, 4-8 in the discovery carriers); no B+/B- call | **no: new** |

Flags cell = `;`-separated tokens, **the first always the call tag** (`pilot-v1`), then the flags. A cell with a panel tag is recognisable as written by the panel (a hand-entered value has none), and can be
recomputed when the version is bumped. Blocking flags (the value is NOT written, the flags cell is): gsdf `low_flank_depth`, `non_integer`, `outside_calibrated_range`, `no_estimate`; transposons `no_call` (fewer than 5 informative
read pairs), `low_gq` (GQ < 10), `site_error`; LG5 `few_markers`, `low_depth`, `no_estimate`; B `low_depth`, `few_scale_bins`, `few_blocks`, `no_estimate`. Information-only flags (the value is written): gsdf `junction_in_two_copy`;
chr7ins `dilution_possible` (a 0/0 with insertion pairs in a fish with 3 or more gsdf copies); LG5 `near_threshold`, `aggregate_disagrees`, `junction_uninformative`, `junction_disagrees`; B `fraction_disagrees`.

## 4. Validation so far (provisional) and what is unverified
- **gsdf**: in the MC family 46 of 46 unflagged calls equal the sheet's `gsdf number` (females 25 of 25 at 2 copies, males 21 of 21 at 3); 13 other samples with sheet values agree 13 of 13; the five MC samples flagged here are
  among the seven the sheet leaves blank. Depth ratios: 2 copies median 0.535 (0.463-0.584), 3 copies 0.796 (0.752-0.885). The calibration (0.530, +0.268 per copy) was fitted on minimap2 alignments, like the pilot's.
  On 188 more samples (every sample of the sheet export with 3 or more copies, 148, plus 40 random 2-copy samples) 188 of 188 unflagged calls equal the sheet (2: 40, 3: 131, 4: 14, 5: 3); depth ratios do not overlap (2: 0.449-0.559, 3: 0.712-0.887, 4: 0.999-1.110, 5: 1.261-1.311). Parity with sexloci, not independent truth. *Unverified*: the other 2-copy values; any mbuna value (the sheet has none).
- **TE1 / TE2**: 146 samples (every non-reference sheet value plus 40 random `0/0`): TE1 133 of 133 and TE2 130 of 130 confident calls equal the sheet (TE1 85 `0/1`, 8 `1/1`, 40 `0/0`; TE2 52 `0/1`, 10 `1/1`, 68 `0/0`), and the samples
  without a confident call are exactly those the sheet leaves blank (13 and 16). This is parity with sexloci (same code, same BAMs), not independent truth. *Unverified*: the other reference-state values (731 TE1, 751 TE2 `0/0`).
- **chr7ins** (no sheet column): against Munby's calls for 47 fish: 3 of 3 homozygotes `1/1`, 12 of 13 heterozygotes `0/1`, 24 of 24 non-carriers `0/0`, the 7 without a call `0/0`; AUC 1.000. The one discordance
  (SAMEA117806076, 63 ref : 6 alt, 4 gsdf copies) carries `dilution_possible`. *Unverified*: a copy-number-aware model; calls in fish with more than 2 copies.
- **LG5 inversion**: 954 SNV markers + the right-breakpoint junction, tested only on pools of known composition (5 pools): OB pools 0.50-0.51 of alt reads at the markers (selected on them), held-out OB pools of the other genus 0.53 and 0.57;
  BB pools 0.001-0.004 (selected on them) but held-out BB-type pools 0.06-0.20 (the alt alleles also occur on standard chromosomes of other species: hence a call from the MEDIAN per-site fraction). Cuts 0.25 and 0.75 are the midpoints of
  0 / 0.5 / 1, set without any diploid fish. One real individual measured (an *A. calliptera*, median 0.000, `std/std`): no individual of known LG5 genotype has been validated. Luwino pool: possibly one inverted chromosome in 40.
  **Treat every LG5 value as unvalidated**; the pilot cohort should show clusters near 0, 0.5 and 1, and the cuts are to be refit on it.
- **B chromosome**: discovery on 12 carrier females and 9 B-free pools (874 core blocks, 836 tier A): carriers 4.7-8.4, B-free 1.10-1.11 (tier A, in sample; leave-one-out carriers 4.3-8.0). The sampled-bin genome-wide depth reproduces the exact one
  within 0.8% on all 21 discovery tracks, and `samtools bedcov` reproduces mosdepth on a real BAM (score 1.0751 against 1.0716, 98.5% of blocks within 5%). One B-free diploid BAM scored 1.075, with 18% of its blocks at 2x or more (pools: 8-10%): the flag
  `fraction_disagrees` (limit 0.30) may fire on poorer data. *Unverified*: any individual carrier through the panel; mbuna males as negatives; the cut-offs 2.5 and 1.6 (they only drive that flag).
- **Not verified here at all**: samtools 1.24 on mzebra (everything was run with the conda samtools of D05/D07); runtime without parasail on a full BAM; BAMs other than the one 9.7 GB BAM timed.

## 5. Runtime and where the code is
- One 9.7 GB (about 17x) BAM on local disk (D05), each module alone as its own process, Python start-up included: gsdf 0.33 s, te 0.46 s (parasail), lg5inv 1.03 s, bchrom 4.26 s; under 80 MB of memory; 154 KB of JSON for the four (bchrom 138 KB).
  The four together are about 6 s. Without parasail the transposon step is slower (the numpy aligner takes about 28 ms per read against the 9.9 kb TE2 allele): my estimate is seconds to about 20 s per sample, not measured on a full BAM.
- Code: github.com/ptmcgrat/CichlidSRSequencing, branch `reorg_spec`, folder `04_genotype_sd_systems/sd_genotypes` (`panel/measure.py`, `panel/calls.py`, `panel/identity.py`, `te/`, `gsdf/`, `lg5inv/`, `bchrom/`), commit `8c831dbb` and the commit that adds
  `panel/calls.py` and the `status`/`version`/`sheet` keys of the four configs. From a checkout: `source tools/env.sh` (puts every module on PYTHONPATH), then `import sd_genotypes...`. Tests: `python tools/run_tests.py 04`.
- Check on a new machine (the numbers for SAMEA117806270, 9.7 GB, measured on D05): gsdf ratio 0.513, 2 copies; chr7ins 0 ref : 12 alt `1/1`; TE1 20 : 0 and TE2 17 : 0 `0/0`; LG5 median site fraction 0.000; B score 1.075.

## 6. Verified on utaka and the first mbuna (2026-11-17)
- **utaka** (a copy of the repository, samtools 1.25, Python 3.11, parasail; data directory `/Data/mcgrath-lab/Temp`, since `/home` is 98% full): `panel check` prints the four keys above; the 9.7 GB BAM of SAMEA117806270 gives the calls of section 5 (gsdf 2, chr7ins `1/1`, TE1 and TE2 `0/0`, LG5 0, B score 1.0751) in 0.38 s (gsdf), 0.54 s (te), 2.00 s (lg5inv) and 8.49 s (bchrom), under 85 MB each.
- **The mbuna**: the sheet export has 573 (72 species), all with a BAM in `Bamfiles/Mzebra_GT3_NCBI` (5.8 TB, median 9.3 GB, largest 34.0 GB: 16 are over 20 GB) and none with a `gsdf number`, `TE1` or `TE2`: there is nothing to compare or conflict with, so the first values are new. The pass must delete each BAM once it is measured (the data directory has 5.0 TB free).
- **Six mbuna through the whole pass** (download 81-94 MB/s, measurement 6-8 s, BAM deleted afterwards; no error): LG5 1 in both *Labeotropheus* females, 0 in both males and in both *Cynotilapia* (median site fractions exactly 0.5 and 0.0); gsdf 4, 4, 4, 5 in *Labeotropheus* (depth ratios 1.08-1.32, junction ratios 0.60-0.85, as for the sheet's 4- and 5-copy samples) and 3 (male) / 2 (female) in *Cynotilapia*; chr7ins, TE1 and TE2 `0/0` wherever called (the 6.5x *Cynotilapia* female has no TE call); B scores 0.985-1.107.
- **What to write for the mbuna**: `gsdf number`, `chr7ins`, `TE1`, `TE2` now. `Inversion5` and `B score` stay in the JSONs until the distribution over about 50 mbuna shows the expected clusters (LG5 near 0, 0.5 and 1; B near 1 against 4-8): six fish cannot set a threshold.
