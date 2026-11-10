# Changelog

## 0.1.26
- **Handover documentation.** `docs/VARIANTS.md` gives, for the gsdf duplication, TE1, TE2 and chr7ins: what is known about the genetic change, how each was detected, and how to genotype it against GT3a
  (methods, constants, validation numbers, limits, open items). Anything the project never recorded (how TE1/TE2 were first called, their element family) is marked NOT RECORDED.
- `scripts/detection/`: the stand-alone scripts that did the detection and validation (dot plot, copy number, junction scan, family lookup, element reconstruction, allele build, concordance),
  with an index, checksums and the small data tables they read. A test checks that each file parses and matches its checksum.
- **Correction (docstring only):** the `chr7_gsdf_ins` note in `known_sites.py` (0.1.25) named a specific GT3a copy as the element's template and implied it was the closest. That was not supported: the coordinate was only one of
  twelve candidate copies listed by an early script run, and the carriers differ from all 25 modal GT3a copies by 14-16 bases per kb. The note now says so. No code or allele changed.
- Two tests tie `docs/VARIANTS.md` to the code (built-in coordinates, allele fingerprints, gsdf constants), so the document cannot silently drift.

## 0.1.25
- **The chr7 insertion upstream of gsdf is built in** (built-in set `chr7_gsdf_ins`, `known_sites.py`): the ~3.8 kb CACTA-type transposon that Munby et al. 2021
  report as a Y allele in *A. calliptera* and some mbuna, at both GT3a copies of the gsdf unit (sites `chr7ins_c1`, NC_135173.1:50,765,050, and `chr7ins_c2`,
  :50,786,787; 603 bp alleles, md5 8b7b5cf1d8c7). Validated against Munby's genotypes for 47 fish: all 16 carriers show the junction, none of 24 non-carriers.
- **Combined calls over sites that are one locus** (`SITE_GROUPS`): `sv-genotype` sums the informative pairs of the two copies and writes one genotype per
  sample, `chr7ins_GT` / `_ref` / `_alt` / `_GQ` / `_reads` / `_status` in `sv_genotypes.tsv`, a summary table in the log, and one `chr7ins` column in the
  Sample Database (the per-copy columns stay in the TSV). The combined call is made when all members are genotyped.
- `--builtin-sites` takes several sets and now defaults to **all** (`deep_benthic_TE` + `chr7_gsdf_ins`), so a plain `sv-genotype` adds the chr7 sites;
  `--builtin-sites deep_benthic_TE` gives the previous behaviour. `write-builtin-vcf --set all` writes every site. `load_builtin_sites()` without
  arguments still returns TE1 and TE2.
- The genotype model is one function, `longindel.call_from_counts(n_ref, n_alt)`, used by the classifier and the combined calls (tested equal on a 22 x 22 grid).
- `cli.build_parser()` exposes the parser (for tests). 8 new tests.

## 0.1.24
- `genome_overview.png`: the known sex-determination loci are marked on every panel and labelled above the first: gsdf (LG7) and the two LG10 transposon
  insertions (TE1/TE2, 216 bp apart, so one shared label). They are placed by GT3a accession or, for VCFs using UMD2a chromosome names, by linkage group;
  a locus beyond the analysed windows is skipped. Positions are in `knownloci.py`. Config: `show_known_loci: false` turns it off.

## 0.1.23
- `sv-genotype` finds the genome FASTA in an NCBI download folder: rna.fna, cds_from_genomic.fna, protein and similar files are ignored, and *_genomic.fna
  is preferred when several whole-genome files remain. `--reference-fasta` still overrides.

## 0.1.22
- The two candidate transposon insertions (NC_135176.1:24,921,178 and :24,921,394) are built into sexloci (`known_sites.py`), so `sv-genotype` needs no
  VCF: `--indel-vcf` is now optional (default: built-in set `deep_benthic_TE`, sites TE1 and TE2). `--no-indel` skips them. The log prints the length and MD5
  of each inserted allele for comparison with the original file. `write-builtin-vcf OUT` writes them as a VCF.

## 0.1.21
- New command **`sv-genotype`**: genotypes the gsdf duplication AND long insertions (transposons) for every sample in the database, streaming only
  the needed BAM regions (one `samtools view` request per sample: the 130 kb gsdf region plus a ~1 kb window per insertion site). Insertion
  genotypes use the lab's `IndelReadClassifier` (new module `longindel.py`: read pairs are scored against a local reference and insertion allele,
  binomial genotype model, phred GQ). Sites come from a VCF-style file (`--indel-vcf`, `--indel-names`); the VCF REF is checked against the
  reference FASTA. Each locus has its own cache, so a new site fetches only its own window, and gsdf results from `gsdf-genotype` are imported.
- `--write-sample-database` now writes one column per locus in a single guarded batch (`gsdf number` plus one genotype column per site).
- The sheet writer sends values as RAW. (USER_ENTERED would let Google Sheets turn a genotype such as 1/1 or 0/1 into a date.) Text values are
  supported next to numbers.
- `sv-genotype --refresh` measures only the SELECTED samples again (use it after adding sequencing to existing samples); results of the others are kept.
- Aligner: parasail (fast, `pip install parasail` or `pip install -e ".[indel]"`) if installed, else a numpy Smith-Waterman (same scores, slower).

## 0.1.20
- `gsdf-genotype --write-sample-database [--yes] [--overwrite] [--include-flagged] [--column-name]`: writes the copy number into a new column
  'gsdf number' of the SampleDatabase worksheet. Dry run by default; touches only that column; matches by SampleID; snapshot, plan log, change
  guard and read-back check. All other sheet access remains read-only.

## 0.1.19
- gsdf-genotype writes a Pedigree column (offspring / parent) from the Parents entries of the sample sheet; gsdf-summary counts independent fish
  only by default (`--include-pedigree` to keep relatives), because one family can inflate the sex comparison.

## 0.1.18
- gsdf-genotype carries the sample sheet's Subgroup column (and any `--extra-columns`) into the table and prints genotype counts by it.
- `sexloci gsdf-summary`: carrier frequency by sex, with a Fisher test and a pattern label, for any column of the table.

## 0.1.17
- `sexloci gsdf-genotype`: copy-number genotype of the gsdf-region tandem duplication for every sample with a GT3a alignment, reading only the
  130 kb around the locus from the cloud BAMs; writes a table with SampleID, Species, Ecogroup, Sex and genotype. Resumable, parallel, with
  subset filters, calibration options and flags. Config keys gsdf_r2 / gsdf_step.

## 0.1.16
- XY/ZW: the pedigree decides when it can (the parent whose inherited haplotype predicts sex is the heterogametic one). The sign of dHet and the
  diagnostic-SNP counts only decide without a pedigree. With one sire a negative dHet (alleles private to the sire's X) is normal for XY and
  had produced a spurious "conflicting"/"ZW" call.
- A run deletes region plots and pedigree outputs left in the folder by an earlier run, so old figures can no longer be mistaken for the current ones.

## 0.1.15
- Broods that share a parent (one sire crossed to several dams) are independent broods by default (`family_independence: pair`); `parent` keeps the
  strict view. Regions whose supporting broods share a parent say so (`units_share_parent`, note in the report).

## 0.1.14
- The pedigree summary lists every parent pair with its number of offspring, and the `shared_parents` note names the parent that links broods
  (broods that share a sire or dam count as one independent family), so a wrong sire/dam in the sample sheet is easy to spot.

## 0.1.13
- Fixed: a local VCF was never re-downloaded once it existed, so a corrected file uploaded under the same name was silently ignored. Each run now
  compares with the cloud copy (size + modification time) and updates it; `folder_summary.tsv` shows downloaded / updated / unchanged.

## 0.1.12
- Chromosome labels (LG1 ...) for the GT3a accessions are read from the annotation's chromosome records, in addition to the Assignment3 table for the
  UMD2a accessions.
- Gene lists work for VCFs called on GT3 sequence under UMD2a chromosome names (YH): chromosomes are matched to the GT3a annotation by linkage group
  (`annotation_alias_by_linkage_group`). Previously the `genes` column was empty for YH. New notes `annotation_matched_by_linkage_group` / `no_annotation`.

## 0.1.11
- `sexloci filter-audit --vcf X`: where along the genome the site filter removes sites (per-window table and whole-genome PNG/PDF with
  the share flagged by each reason), to judge whether a filter, e.g. mapping quality, is thinning a candidate region.

## 0.1.10
- Broods with a missing parent are handled explicitly: a parent that is listed but not in the VCF, or a Parents entry with only one parent
  (`sire,` / `sire,?` / `,dam` / a single known sample ID), excludes that brood from the transmission test only. Its offspring stay together
  as one family in the sex contrasts (previously a one-parent entry put them in the unrelated pool) and count as an independent unit there.
- A parent that is in the VCF but whose brood cannot use the transmission test is no longer treated as an unrelated animal in the contrasts.
- `parent_not_in_vcf` and the new `parent_unknown` say which brood, how many offspring and who is missing (why / how to fix, `fix_list.tsv`);
  `few_independent_families` states how many families can use the transmission test; the pedigree summary counts one-parent entries.

## 0.1.9
- **Differentiated core** interval per region (`core_start`/`core_end`, windows where male-female heterozygosity is an outlier),
  reported next to the transmission plateau, which overshoots because linkage fades slowly. On the YH data the LG10 core is 11.75-30.00 Mb
  (plateau 11.75-31.00 Mb).
- The XY/ZW decision only uses evidence that fired in the region (a window-wide dHet average no longer outvotes diagnostic SNPs); this removes
  the spurious "conflicting" call on a 3-SNP region.
- Fst / dHet / depth outliers need >= 30 markers in the window (`min_evidence_markers`): a 5-marker window on the repetitive LG3 produced a tier C lead.
- `narrow_spike` flag for regions of a few SNPs within < 20 kb with no shoulders.

## 0.1.8
- The pedigree panel (genome overview and region plots) now shows the effect size, mean |r| between inherited haplotype and sex
  (0-1), instead of -log10 permutation p, which saturated at 1/(n_perm+1) and drew every strongly linked window at one height.
  Significance is shown as circled points. README explains association versus pedigree and why a flat top at |r| = 1 is expected.

## 0.1.7
- **Whole-genome plots** (PNG + PDF): `genome_overview` now labels chromosomes as linkage groups (LG1 ... LG23, from the
  Assignment3 table); new `genome_haplotypes` (inherited parental haplotype of every offspring along the genome, sorted by sex,
  per family and parent) and `genome_by_chromosome` (strongest haplotype-sex association per chromosome, sire vs dam).
- `filter-compare` writes `retention_genome` (per-window retention of an existing filtered VCF, depleted windows highlighted).
- `ignore_filter_labels` config option: do not honour chosen FILTER labels (e.g. `allele_freq`, `min_DP`) to measure what each
  filter costs.
- Silenced the harmless "Mean of empty slice" warnings from windows with no depth data.

## 0.1.6
- **Low coverage no longer penalises the pedigree analysis.** The transmission test now uses its own genotype thresholds
  (offspring depth >= 2, no GQ filter; parents depth >= 6 / GQ >= 20) and is exempt from the call-rate filter. In simulation it
  recovers XY and ZW loci with offspring at 3x and 2x (tier A, all families supporting); before, the pipeline collapsed at 8x
  (424 of 9,902 sites) and failed at 5x and 3x.
- Sex contrasts keep strict genotypes, but samples that keep < 50% of their genotypes are left out of them (listed, with
  `sample_qc.tsv` role) instead of removing sites for all samples. All-excluded + pedigree -> transmission-only; all-excluded
  and no pedigree -> `too_few_usable_samples` block.
- Default `min_gq` is now 0 (the lab script filters on DP only); GQ >= 20 was an unrequested filter that discards most genotypes below ~7x.
- Coverage warnings are role-aware: family offspring are warned only below 2x; low coverage of offspring that are left out of
  the contrasts is a note. `expected_site_loss` removed (superseded by the exclusion above).
- Fixed a false alarm: the parentage (Mendelian) check treated low-depth heterozygote under-calls as errors and capped the
  tier; it now uses only sites where both parents are homozygous for the same allele.

## 0.1.5
- Low depth is now quantified instead of silent: the low-coverage warning says what share of a sample's genotypes survive
  `min_dp` (Poisson estimate) and lists every sample in `fix_list.tsv`; a new `expected_site_loss` warning predicts, from the
  alignment-table coverage, how many sites will reach the call-rate filter.
- The run reports how many sites the call-rate filter removed (`call_rate_filter`, plus `heavy_call_rate_filtering` above 30%)
  and lists it in `variant_filter_summary.tsv`. Previously those sites were dropped without being counted.

## 0.1.4
- Blocking messages now say what is wrong (counts, percentage, breakdown by value, affected IDs compressed into ranges),
  why it blocks and how to fix it; near-miss sample IDs are suggested (`U_M1 -> U-M1`).
- New `fix_list.tsv` per VCF: one row per affected sample (problem, sample, column, current value, needed value);
  removed automatically once the problem is fixed. `findings.tsv` / `folder_findings.tsv` gain `why` and `fix` columns.
- A pedigree-completeness note is always reported (also when blocked): how many samples list `sire,dam`, how many
  parent pairs, and how many listed parents are in the VCF.
- Fixed: empty spreadsheet cells (pandas NA) were not recognised as blank.

## 0.1.3
- `run-folder` now prints the blocking and warning findings of every VCF (code + message, with sample IDs) and writes
  `folder_findings.tsv`; `--check-only` rows carry `block_codes` (previously the reason for a block was only in each VCF's
  `preflight_report.md`).

## 0.1.2
- The lab FileManager is no longer needed: its read-only parts are packaged as `sexloci.labdata.LabData`
  (rclone list/download, SampleDatabase / DNAReadsDatabase / alignment sheets, genome GFF). It has no write path,
  authenticates lazily, requests read-only Google scopes and raises on rclone errors instead of ignoring them.
  The `--filemanager` option is gone; the lab sheets are used by default (`--sample-table` etc. for offline use).
- Default output location `~/Temp/SexDeterminationAnalysis/` (config key `output_root`); single runs write to
  `<root>/<vcf name>/`, `run-folder` to `<root>/<vcf name>/` per VCF plus `<root>/folder_summary.tsv`.
- New dependencies: gspread, gspread-dataframe, google-auth (rclone must be installed and configured).

## 0.1.1
- Zip files are now versioned (`sexloci-<version>.zip`).
- Default genome is `Mzebra_GT3_NCBI`; `run-folder` downloads its GFF so regions list genes.
- `acknowledged_codes` config option to show reviewed warnings as notes.
- `run-folder`, `list-vcfs`, `describe` and `filter-compare` commands; variant-level filtering of bad sites.

## 0.1.0
- Strict preflight (BLOCK/WARN/INFO), per-stratum association, pedigree transmission test,
  independent-unit replication tiers, simulation test-suite.
