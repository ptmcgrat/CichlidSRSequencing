# sexloci

Map sex-linked loci in a cichlid species from a **single-species VCF** plus the lab sample / alignment databases.
Combines pedigree transmission mapping, family-aware association and Fst/heterozygosity contrasts, calls **XY vs ZW**,
never assumes an inversion, and refuses to run (or warns loudly) when the data cannot support the analysis.

## Install (from a downloaded zip, into your home directory)
```bash
cd ~
[ -d sexloci ] && mv sexloci sexloci.bak_$(date +%Y%m%d_%H%M%S)   # keep the previous version, if any
unzip ~/Downloads/sexloci-0.1.26.zip -d ~                            # creates ~/sexloci
cd ~/sexloci && pip install -e ".[fast]"                            # [fast] adds scikit-allel (recommended for real VCFs)
sexloci --help && python -c "import sexloci; print(sexloci.__version__)"
```
Needs `rclone` with the lab remote (`ptm_dropbox:`) already configured. The install is editable: it always runs the code in
`~/sexloci`; to upgrade, repeat the commands with the new zip name.

## Downloads stay in step with the cloud
Every run compares each VCF (and its index) with the cloud copy by size and modification time and fetches it again if it differs, so a
corrected file uploaded under the same name is picked up automatically. `folder_summary.tsv` has a `download` column (downloaded / updated /
unchanged) and `vcf_modified`; `run_summary.json` records the VCF's size and modification time. `--refresh` forces a fresh download.

## Where things go
| what | default location |
|---|---|
| **results** | `~/Temp/SexDeterminationAnalysis/<vcf name>/` (`output_root` in the config; `--outdir` overrides) |
| downloaded VCFs, GFF, credentials | `~/Temp/CichlidSequencingData/...` (mirrors the cloud folders, as the lab's other code does) |
| default cloud folder for `run-folder` | `Outputs/PedigreeVCFs/` under the lab master directory |

Re-running a VCF overwrites its own results folder only.

## Lab data access (built in, read-only)
`sexloci.labdata.LabData` packages the parts of the lab FileManager that are needed: `rclone lsf` / `rclone copy` between the
cloud and local master directories, the *SampleDatabase*, *DNAReadsDatabase* and genome-version alignment sheets, and the
genome GFF. It cannot write: no sheet is modified, Google is asked for read-only scopes (`sheets_read_only: true`), the
credential file is fetched from the cloud only the first time a sheet is read, and rclone errors are raised, not ignored.
Sheet IDs, remote and directories are config keys (`config/example_config.yaml`). Offline alternative: export the sheets and pass
`--sample-table`, `--alignment-table`, `--reads-table`.

## Run
```bash
sexloci list-vcfs                                    # what is in Outputs/PedigreeVCFs/
sexloci run-folder --check-only                      # download + preflight every VCF (no analysis)
sexloci run-folder --contig-regex '^NC_'             # analyse every VCF, chromosomes only
sexloci run-folder --include 'YH*' --manifest manifest.tsv

sexloci check --vcf ~/Temp/CichlidSequencingData/Outputs/PedigreeVCFs/saulosi_pedigree.vcf.gz   # one local VCF
sexloci run   --vcf <local vcf> --n-perm 1000
```
`run-folder` downloads each VCF (+ `.tbi`), runs each as its own analysis unit and records every outcome (OK / BLOCKED / ERROR /
SKIPPED) in `folder_summary.tsv`; a blocked or failing file never stops the batch. `--genome-version` defaults to
**Mzebra_GT3_NCBI** (it selects the alignment sheet and the GFF so regions list their genes); a manifest
(`examples/manifest_example.tsv`: `vcf`, `genome_version`, `skip`, `note`) overrides it per file. RefSeq chromosomes are `NC_*`
(unplaced scaffolds `NW_*`). Because this reference is a mbuna, every non-mbuna species gets a `reference_divergence` warning;
once reviewed, `acknowledged_codes: [reference_divergence]` shows it as a note.

## Coverage: each analysis gets the depth it needs
| analysis | genotypes used | what low coverage does |
|---|---|---|
| **transmission (pedigree)** | offspring: depth >= 2, no GQ filter; parents: depth >= 6 and GQ >= 20; sites are not subject to a call-rate filter | **tolerant**: it reads which parental haplotype each offspring inherited from hundreds of informative markers per window, recombination is rare, and low depth only turns heterozygotes into homozygotes. Works at ~2-3x in simulation. Warned only below 2x; parents should be >= 10x |
| **sex contrasts** (association, Fst, heterozygosity, diagnostic SNPs) | depth >= `min_dp` (5), no GQ filter (`min_gq: 0`, as in the lab script); sites need `min_call_rate` among the contrast samples | needs confident genotypes: a sample that keeps < 50% of its genotypes is **left out of the contrasts** (listed in `samples_excluded_from_contrasts`, and in `sample_qc.tsv`) instead of removing sites for everyone. Family offspring left out this way still count for the transmission test |

If every contrast sample is excluded and a pedigree exists, the result rests on the transmission test alone (replication across
independent families is still assessed); with no pedigree that is a block (`too_few_usable_samples`). Raise or lower
`min_dp`, `ped_min_dp`, `ped_parent_min_dp` and `max_sample_missing_contrast` in the config if needed.

## Bad-variant filtering
Before any statistic: a FILTER label other than PASS/`.` removes a site, and so does any hard filter that fires (QUAL<30, QD<2,
FS>40, MQ<50, SOR>3, MQRankSum<-12.5, ReadPosRankSum<-8; skipped when the INFO field is absent) or INFO/DP above the 99th
percentile (collapsed repeats). The same logic handles raw, FILTER-annotated and already-filtered files; what was removed
and why is in `variant_filter_summary.tsv`, with a warning if the VCF carries nothing to filter on or >60% is removed.
`ignore_filter_labels: [allele_freq, min_DP]` (config) stops chosen FILTER labels from being honoured and `variant_filters: []` switches
off the hard filters, so you can see what any single filter costs.
**Deliberately not applied:** allele-frequency, `InbreedingCoeff`/`ExcessHet` and a lower depth cutoff, because in an XY region all
males are heterozygous and hemizygous X/Y sequence has low depth in one sex: these filters would delete the sex-linked sites.
```bash
sexloci filter-audit --vcf MC.vcf.gz                                       # WHERE along the genome does the filter remove sites?
sexloci describe --vcf a.vcf.gz b.vcf.gz                                  # raw / annotated / already filtered?
sexloci filter-compare --raw YH_raw.vcf.gz --filtered YH_filtered.vcf.gz   # what did the filtered file lose?
```

## Whole-genome plots (every run, PNG + vector PDF)
| file | what it shows |
|---|---|
| `genome_overview` | every chromosome (LG1, LG2, ...) in one view: association (max chi2), pedigree transmission (mean \|r\| between the inherited sire/dam haplotype and offspring sex, 0-1; circled = significant), Fst, male-female heterozygosity, depth ratio per window; gold = called regions; dashed/dotted lines = the known sex-determination loci: gsdf (LG7, purple) and the two transposon insertions TE1/TE2 (LG10, green), placed by accession or linkage group, so UMD2a-named VCFs get them too (switch off with `show_known_loci: false` in the config) |
| `genome_haplotypes` | **inheritance along the genome**: one panel per family and parent; rows = offspring (females red / males blue in the bar at right), columns = windows, colour = which of the parent's two haplotypes was inherited. Blocks are haplotypes, a single row switching colour is a recombinant, and a sex-linked region shows as a block where the colours split by sex for one parent only |
| `genome_by_chromosome` | per chromosome and family, the strongest association between inherited haplotype and sex (sire vs dam): a sex chromosome stands out for one parent in every family |
`sexloci filter-compare` additionally writes `retention_genome` (sites kept per window by an existing filtered VCF, depleted windows
in red). Chromosome names use the linkage-group table from the lab's Assignment3 script (`contig_labels` in the config; check it fits the
assembly). Colours in `genome_haplotypes` are arbitrary per chromosome (the two parental haplotypes), constant along it.

## Rules (decided with the lab)
**BLOCK** (the run stops; `preflight_report.md` lists the offending sample IDs)

| code | rule |
|---|---|
| `sex_not_mf` | every VCF sample must have Sex exactly `M` or `F` (`U`, `F?`, `M?`, blank all block) |
| `species_mismatch` / `species_missing` | all VCF samples must have the identical Species (parents included; cross families are not run) |
| `sample_not_in_table` | every VCF sample ID must exist in the sample table (exact match) |
| `genome_version_mismatch` | alignment records must come from `--genome-version` |
| `too_few_per_sex` | at least 3 analysed males and 3 females |

**WARN** (some cap the confidence tier at B): low sex counts (<6), small families (<12), fewer than 3 independent families, no
usable pedigree, parents missing from the VCF / wrong sex / malformed `Parents` (with suggested ID matches, never applied
silently), low coverage (family offspring <2x; other animals <5x, <3x unusable; parents <10x), coverage spread (note when a pedigree drives the analysis), alignment-QC outliers, reference species far from
the samples, sex confounded with batch/project/lab-reared/brood (**cap B**), sex-depth difference (**cap B**), heterozygosity
tracking depth, high sample missingness, VCF-vs-table coverage mismatch, Mendelian errors against listed parents (**cap B**),
no/over-aggressive variant filtering. `Inversion10` and `Category` are never read.

## Outputs (per VCF)
A blocked run writes `preflight_report.md` (what is wrong, why it blocks, how to fix, affected IDs), `findings.tsv` and `fix_list.tsv`.
For the full set after a successful run:
`report.md` (findings + regions), `sex_linked_regions.tsv`, `window_stats.tsv`, `sex_associated_markers.tsv`,
`stratum_window_stats.tsv`, `family_transmission_stats.tsv`, `variant_filter_summary.tsv`, `pedigree_chromosome_summary.tsv`,
`sample_qc.tsv`, `findings.tsv` (with why / fix), `fix_list.tsv` (blocked runs: one row per sample to correct), `genome_overview.png`, `region_*.png`, `run_summary.json`.

## Association versus pedigree (the two panels)
**Pedigree (transmission)** asks, per window, which of a parent's two haplotypes each offspring inherited (read from markers where one
parent is heterozygous and the other homozygous) and whether that predicts offspring sex. It needs genotyped parents, says whether
the sire (XY) or dam (ZW) carries it, and works even when X and Y barely differ in sequence. **Association** ignores parentage: per SNP,
does the genotype (alt dosage or heterozygosity) differ between males and females? It needs no parents and is sharp in unrelated animals.
Within one family both use the same offspring information, so they never count as independent confirmation.
In the pedigree panel a **flat top at |r| = 1** is expected for a sex-linked region without recombination in the heterogametic parent
(e.g. an inversion): every window in the block carries the same haplotype-sex information. The slopes on either side are where
recombinant offspring break the link. The panel shows the effect size (|r|) rather than -log10 p because permutation p-values saturate
at 1/(n_perm+1) and would draw every strongly linked window at the same height.

## Chromosome names and gene lists
VCFs on this lab's assemblies use two naming systems: the GT3a RefSeq accessions (`NC_1351xx`) and, for older calls such as the YH pedigree,
GT3 sequence under the UMD2a accessions (`NC_0367xx`). Plots label chromosomes as linkage groups (LG1 ... LG23): the `NC_0367xx` names come from the
table in the lab's Assignment3 script (`contig_labels`), the `NC_1351xx` names from the chromosome records of the annotation (GFF). Genes inside a
called region are listed from the GT3a annotation; for a VCF whose chromosome names are not in the annotation, the chromosomes are matched by linkage
group (`annotation_alias_by_linkage_group`, on by default; reported as `annotation_matched_by_linkage_group`). That assumes the VCF coordinates
are those of the annotated GT3a sequence, which holds for YH; switch it off for a VCF called against a different assembly.

## Genotyping the gsdf-region duplication across the sample database
`sexloci gsdf-genotype` genotypes the tandem duplication of the 21.7 kb unit around gsdf (NC_135173.1:50,751,322-50,794,777 in GT3a) for every
sample that has an alignment to the GT3a reference, and writes `gsdf_genotypes.tsv`: SampleID, Species, Ecogroup, Sex, Genotype, CopyNumber,
EstCopies, DepthRatio, FlankDepth, JunctionJ, Flag, Coverage, Status. It needs `samtools` and `rclone` on the PATH (nothing is written to the cloud).

    sexloci gsdf-genotype --limit 10                      # trial run on 10 samples first
    sexloci gsdf-genotype --jobs 4                        # the whole database; safe to interrupt and rerun (it resumes)
    sexloci gsdf-genotype --ecogroup Mbuna --species zebra --samples 'MZ_*'     # subsets

For each sample only the 130 kb around the locus is read from the cloud BAM, through a read-only rclone web server on
`Bamfiles/<genome version>/` (a few hundred kB per sample, kept in `region_bams/`; `--no-keep-bams` deletes them). The reference carries two tandem
copies, so depth over the unit relative to the flanks is (total copies)/4. Genotype notation: `-/-` two copies in total, `dup/-` three (one
chromosome carries the tandem duplication, e.g. a male whose Y does), `dup/dup` four (both do, like the reference); `CN1`, `CN5`, `CN6` otherwise.
Depth alone cannot say which chromosome carries the duplication. `JunctionJ` is the uniquely mapped coverage where the two copies join, relative to
the flank: about 0 in samples without the tandem arrangement, 0.2-0.4 in those with one duplicated chromosome; it is specific but not sensitive,
so read it as confirmation only.
A `Subgroup` column from the sample sheet (any spelling) is carried into the table automatically; `--extra-columns NAME ...` adds others.
The table also has a `Pedigree` column (offspring / parent, from the sheet's Parents entries). Relatives share haplotypes, so one family can make
a result look far stronger than it is: `gsdf-summary` leaves pedigree fish out by default (`--include-pedigree` keeps them).
`sexloci gsdf-summary --by Subgroup Ecogroup Species --min-n 5` summarises the table by any column: samples, carriers (at least one
duplicated chromosome), carriers among males and females with a Fisher test, dup/dup, 5+ copies, flagged, and a pattern label (absent,
male-limited, male-biased, present not sex-associated, no sex information; the rules are in the function docstring). Rerunning
`gsdf-genotype` with the same output folder reuses the cached measurements, so adding a column needs no new downloads.
Calibration: by default the MC family (2 copies -> ratio 0.530, +0.268 per extra copy; bwa/GT3a, 4-8x). For data from another pipeline give
`--anchors anchors.tsv` (SampleID, copies of a few samples of known state, at least two different copy numbers), or `--theory` (ratio = copies/4);
the raw `DepthRatio` is always in the table so you can recalibrate. Flags: `low depth` (flank depth < 3x), `uncertain` (between two classes),
`outside calibrated range`. Samples without a GT3a alignment, with coverage below `--min-coverage`, or without a BAM/index in the cloud are listed
with that Status and no genotype. Species in which the locus is rearranged can look odd in depth alone; check the flags.

## Writing the result into the Sample Database
`sexloci gsdf-genotype --write-sample-database` adds a column **gsdf number** (the integer total number of gsdf copies in the diploid genome: 2 = no
duplicated chromosome, 3 = one, 4 = two; the same scale as the 'gsdf copies' of Munby et al.) to the SampleDatabase worksheet. It is the only
thing in sexloci that can modify a sheet, and it is built to be hard to misuse:

    sexloci gsdf-genotype --jobs 4 --write-sample-database          # DRY RUN: prints and saves the plan, writes nothing
    sexloci gsdf-genotype --jobs 4 --write-sample-database --yes    # really writes

 * Only that one column is touched, and only the cells that change (one batch request); rows are matched by SampleID, never by position.
 * The column is created at the end of the header if it does not exist (any capitalisation of 'gsdf number' is recognised; two columns with
   that name stop the run). A note on the header cell says how the values were made.
 * Existing values are kept unless `--overwrite` (use it after recalibrating); conflicts are listed in the plan file. `--column-name` renames it.
 * Flagged calls (low depth, uncertain, out of calibrated range) are left blank unless `--include-flagged`; samples without a call stay blank.
 * Before writing it saves a full snapshot of the sheet and the plan (old and new values) in the output folder, refuses to write if the SampleID
   column changed while it ran (someone sorted or inserted rows), and reads the column back to check it.
 * The service account in the lab credential file needs EDIT access to the spreadsheet; reading elsewhere in sexloci stays read-only.

### Genotyping transposon insertions together with gsdf: `sv-genotype`

`sexloci sv-genotype` streams only the parts of each BAM that are needed (one request per sample) and genotypes the gsdf duplication and the
insertions built into sexloci (`known_sites.py`) at once: the two candidate transposon insertions on NC_135176.1 (LG10: TE1 at 24,921,178 and TE2
at 24,921,394; set `deep_benthic_TE`) and the transposon insertion directly upstream of gsdf on NC_135173.1 (set `chr7_gsdf_ins`: the chr7-ins Y
allele of Munby et al. 2021, see below). No VCF is needed.

    sexloci sv-genotype --limit 10                                     # trial run on 10 samples
    sexloci sv-genotype --jobs 4                                       # the whole database; resumable
    sexloci sv-genotype --jobs 4 --write-sample-database               # DRY RUN: shows the plan, writes nothing
    sexloci sv-genotype --jobs 4 --write-sample-database --yes         # writes
    sexloci sv-genotype --no-indel                                     # gsdf only
    sexloci sv-genotype --builtin-sites chr7_gsdf_ins --no-gsdf        # only the chr7 insertion upstream of gsdf (default: all built-in sets)
    sexloci sv-genotype --indel-vcf my_sites.vcf --indel-names A B     # other insertions, from a VCF-style file (CHROM POS ID REF ALT ...)
    sexloci write-builtin-vcf transposons.vcf --set all                # the built-in sites as a VCF file (--set deep_benthic_TE | chr7_gsdf_ins | all)

The log prints the length and an MD5 checksum of each inserted allele (`ALT 9,869 bp, md5 d2b75a3b30e2`). To confirm they equal your original VCF:

    python3 -c "import sys,hashlib; [print(l.split('\t')[1], len(l.split('\t')[4]), hashlib.md5(l.split('\t')[4].upper().encode()).hexdigest()[:12]) for l in open(sys.argv[1]) if l.strip() and l[0]!='#']" original.vcf

* **Reads.** Per sample ONE `samtools view` request through the read-only rclone web server fetches the 130 kb gsdf region and a window of
  (flank + 200) bp on each side of every insertion (overlapping windows are merged). Unmapped mates stored at their mate's position are
  included, since they carry the inserted sequence. Duplicates, secondary and supplementary alignments are skipped.
* **Genotype.** Every read pair is aligned to a local reference allele and a local insertion allele (the same flank on both sides) and counted
  as supporting REF, ALT, equal or uninformative; only reads crossing a junction with at least 15 bp on both sides are informative. A
  binomial model (P(alt) = 0.02, 0.5, 0.98 for 0/0, 0/1, 1/1) gives the genotype and a phred-scaled quality. Fewer than
  `--min-informative-pairs` (default 5) pairs gives `./.`.
* **The reference.** The FASTA is only read for the ~500 bp around each site. `--reference-fasta` names it (default: the only FASTA in the genome
  folder). The VCF REF base is checked against it; a mismatch stops the run and says whether the position is off by one.
* **Caches.** `gsdf_measurements.tsv` and `te_measurements.tsv` (one row per sample and site). A later run fetches only what is missing: a new
  site fetches its own window, not the gsdf region again. Finished gsdf measurements from `gsdf-genotype` are imported on the first run.
* **New sequencing on existing samples.** Their cached results are stale. `--refresh --samples 'MC_*'` (or `--species`, `--ecogroup`) measures just the
  selected samples again and keeps everyone else's results.
* **Output.** `sv_genotypes.tsv`: gsdf call plus, per site, `<name>_GT`, `_ref`, `_alt` (read pairs), `_GQ`, `_reads`, `_status`.
* **Sample Database.** `--write-sample-database` adds `gsdf number`, one column per site (named like the site: `TE1`, `TE2`) and the combined `chr7ins`
  column (below), holding `0/0`, `0/1` or `1/1`.
  Only calls with GQ >= `--min-gq` (default 10) are written; no-calls and failed samples stay blank. Everything else is as for gsdf: dry run
  unless `--yes`, only the new columns are touched, SampleID matching, change guard, snapshot, read-back check. Values are sent RAW, so a
  genotype can never turn into a date.
* **Aligner.** parasail if installed (`pip install parasail`), otherwise a numpy implementation with the same scoring that is several times slower.


For the full methods of all four variants (the gsdf duplication, TE1, TE2 and chr7ins: what is known, how each was detected, and how to genotype it against GT3a), see
[`docs/VARIANTS.md`](docs/VARIANTS.md); the detection scripts are in [`scripts/detection/`](scripts/detection/README.md).

#### The chr7 insertion upstream of gsdf (`chr7ins`)

The Y allele that Munby et al. (2021) found in *A. calliptera* and some mbuna is a ~3.8 kb CACTA-type DNA transposon (class DNA/CMC-EnSpm, 13 bp
terminal inverted repeat, 2 bp target-site duplication TG) inserted 2,564 bp upstream of the first coding exon of gsdf. GT3a holds two copies of the gsdf unit, and a carrier
has one copy per chromosome, so its reads split between the two reference copies. `chr7_gsdf_ins` therefore holds the two junction alleles
(`chr7ins_c1` at NC_135173.1:50,765,050 and `chr7ins_c2` at :50,786,787; each is the anchor base, the first and last 300 bp of the element and the
duplicated TG), and sv-genotype writes **one combined call**, `chr7ins`: the informative pairs of both copies are summed and called with the usual
model (the same function the single sites use). The per-copy results stay in `sv_genotypes.tsv`; the Sample Database gets only `chr7ins`. The combined
call needs both members, so renaming them (`--indel-names`) or genotyping one alone gives per-copy calls only.

* **Validation** (47 fish with Munby's calls): all 16 of his carriers show the junction and none of his 24 non-carriers does (4 of them from the
  carriers' own population, Lake Kingiri); 15 of 16 carriers get his genotype (3 of 3 homozygotes, 12 of 13 heterozygotes), 24 of 24 non-carriers and
  7 of 7 fish without a call are `0/0`; across 1,043 aligned fish, all 874 confident calls in the other ecogroups (benthics, *Diplotaxodon*, *Rhamphochromis*, utaka, hybrids) are `0/0`.
* **Copy number.** The model expects half of the informative pairs from the insertion in a heterozygote. A fish with three or four gsdf copies has
  extra reference-like copies of this locus, which dilute the insertion allele (about a third at three copies). The one discordant fish
  (SAMEA117806076: three gsdf copies, 6 insertion pairs of 69 = 9%, Munby's own ratio 20%) was called `0/0`; the dilution probably contributes, but its
  insertion fraction is below even a third, so it stays uncertain. Treat a `0/0` with a few insertion pairs (`chr7ins_alt` in the TSV) in a fish with more
  than two gsdf copies (`gsdf number` in the Sample Database) as uncertain.
* **Genome version.** The alleles are for Mzebra_GT3_NCBI (NC_135173.1); the built-in anchor bases are checked against the FASTA like any other site.

## Broods that share a parent
One male crossed with several females is normal in cichlids. By default (`family_independence: pair`, or `--family-independence pair`) every distinct
sire x dam pair is an **independent brood**: separate meioses with different partners. The shared parent's haplotypes are common to both broods,
so agreement between them shows consistency across partners, not that the locus is shared across unrelated parents; the run says so (`shared_parents`
note, and a line under any region whose supporting broods share a parent, column `units_share_parent`). `family_independence: parent` is the strict
view: broods that share a sire or dam count as one family, which caps the confidence tier at B.

## Broods with a missing parent
A brood whose sire or dam is **not in the VCF** (or lists only one parent) is left out of the transmission (pedigree) test, because informative
markers need both parents. Its offspring are *not* dropped: they stay together as one family in the sex contrasts (association, Fst,
heterozygosity), still count as an independent unit there, and a parent that is in the VCF is kept out of the contrasts (it is a relative, not an
unrelated animal). `parent_not_in_vcf` / `parent_unknown` say which brood, how many offspring and who is missing, and `fix_list.tsv` lists the
parents to add. Accepted ways of writing a missing parent in `Parents`: `sire,` / `sire,?` / `,dam`, or a single known sample ID (its recorded
sex decides sire or dam). If no brood has both parents the transmission test is switched off and the sex contrasts run alone.

## Plateau versus core, and spikes
A called region reports two intervals. The **transmission plateau** (`start-end`) is where the sire's (or dam's) inherited haplotype
predicts sex; with few recombinants it fades slowly, so it overshoots. The **differentiated core** (`core_start-core_end`) is where male-female
heterozygosity is an outlier (the heterogametic signature) and is the better estimate of the sex-differentiated segment's edges. A region
whose few significant SNPs span < 20 kb with no shoulders in the neighbouring windows is flagged `narrow_spike`: a linked locus in a family
raises neighbouring windows too, so a spike usually means Y/W-linked sequence mapping elsewhere, a mapping artifact or a small structural
variant. Population-genetic outliers (Fst, dHet, depth) only count from windows with >= 30 markers.

## Reading the confidence tier
Replication across **independent units** (families sharing no parent, and the pool of non-family animals) drives confidence:
**A** = combined test significant and >= 2 independent units significant; **B** = combined test significant only; **C** =
population-genetic outlier only. A single pool or a single family can reach B at most. In small pedigrees linkage spans much of a
chromosome: use `start-end` (plateau), `linked_start-linked_end` and `sig_marker_span` together. See `docs/METHODS.md`.

## Tests
`pytest`: simulated cohorts with a known XY/ZW locus, every BLOCK/WARN rule, null datasets, variant filtering, and the lab-data
layer (real path/error logic over a fake rclone and fake Google modules). Not exercised here: the real rclone/Dropbox, the real Google
sheets, and the `scikit-allel` reader.
