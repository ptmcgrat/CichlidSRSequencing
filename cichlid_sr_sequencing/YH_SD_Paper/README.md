# YH_SD_Paper: sex determination in YH and MC (analysis code, Mzebra_GT3_NCBI)

Run everything on **utaka.biosci.gatech.edu** from `~/CichlidSRSequencing/cichlid_sr_sequencing/YH_SD_Paper/`
(`helper_modules/` is found automatically in the parent folders). Data are NEVER stored in git: **Dropbox is the system of record**.

## Layout
    YH_SD_Paper/
      config/paper.yaml                 shared settings (genome, cores, pedigree, BC pipeline parameters)
      config/callsets.yaml              every VCF in the paper = a named callset defined here
      config/parental_sex_region.yaml   parameters of the parental-species analysis
      yhsd/                             library: common, genome, runner (manifests + Dropbox), callset, inheritance, popgen, HMM
      yhsd/analyses/bc_genotyping/      BC fish: original vs v2 genotyping, validation, comparison
      yhsd/analyses/parental_sex_region/  YH and MC broods: inheritance, Fst, heterozygosity, depth ratio -> candidate regions
Each analysis folder gets its own README (question, inputs, run command, outputs, figure it feeds). Add new analyses as
`yhsd/analyses/<name>/steps.py` + `config/<name>.yaml`; build any new VCF through `callsets.yaml`.

## Dropbox record keeping (automatic)
Every step writes an immutable run folder  `.../YH_SD_Paper/<genome>/<analysis>/<step>/<run_id>/`  (run_id = date_time + git commit)
containing the outputs and `manifest.json` (command, git commit/dirty flag, config, tool versions, input/output files with sizes+md5).
The runner uploads it to the same path on Dropbox (rclone) and verifies it with `rclone check`; a `LATEST` pointer names the current run.
Later steps read the latest run, pulling it from Dropbox if it is not on local disk. Directories of many small files are tarred.
`--no-upload` disables this; `--test` runs write to `YH_SD_Paper_test/` and never upload. BC working directories (`v2/`, `original/`, `qc/`)
are published as snapshots after each step.

## A. Workflow: the laptop does git, utaka does `git pull`
Run commands on utaka from `~/CichlidSRSequencing/cichlid_sr_sequencing/YH_SD_Paper/`. Never commit on utaka; it only pulls.

Laptop, for every new zip (use the exact new file name each time):
    cd ~/CichlidSRSequencing/cichlid_sr_sequencing
    git checkout YH_SD_Paper
    rm -rf YH_SD_Paper && unzip ~/Downloads/<zip name>.zip
    git add YH_SD_Paper && git commit -m "<message>" && git push -u origin YH_SD_Paper
utaka, to receive it:
    cd ~/CichlidSRSequencing && git pull --ff-only origin YH_SD_Paper
    cat cichlid_sr_sequencing/YH_SD_Paper/VERSION.txt
To pick up lab-wide changes to helper_modules, merge on the LAPTOP (`git checkout YH_SD_Paper && git merge mcgrath_dev && git push`), then pull on utaka.
Alternative without git on utaka: `bash YH_SD_Paper/scripts/deploy_to_utaka.sh` on the laptop (rsync + selftest). Do not mix the two on the same checkout.

### One-time: make utaka's branch an exact mirror of origin (if it still has the old unzipped folder or an upstream-less local branch)
    cd ~/CichlidSRSequencing
    git status --short | head; git branch -vv; git log origin/YH_SD_Paper..YH_SD_Paper --oneline     # local-only commits? back them up:
    git branch backup_utaka_YH_SD_Paper                          # only if the last command listed commits
    mv cichlid_sr_sequencing/YH_SD_Paper ~/YH_SD_Paper_old_unzipped      # untracked copy would block the checkout
    git fetch origin && git checkout -B YH_SD_Paper origin/YH_SD_Paper
    git branch --set-upstream-to=origin/YH_SD_Paper YH_SD_Paper

## B. Environment and tests
    conda env create -f environment.yml && conda activate yh_sd_paper
    python -m yhsd.cli selftest        # synthetic tests, no data needed; must end with ALL SELFTESTS PASSED
    python -m yhsd.cli chroms          # RefSeq <-> LG table for Mzebra_GT3_NCBI; CHECK the LG labels
    python -m yhsd.cli bc_genotyping provenance

### Trial on real data (writes to YH_SD_Paper_test/, never uploads)
    T="--test --chroms LG10"
    python -m yhsd.cli callset build yh_broods --dry-run        # lists the samples that would be genotyped
    python -m yhsd.cli $T callset build yh_broods
    python -m yhsd.cli $T parental_sex_region scan --species YH
    python -m yhsd.cli $T parental_sex_region summarize --species YH
Send me any traceback. Untested here (no GATK/bcftools/rclone/mosdepth in the build environment): `callset`, `depth`, VCF loading in
`scan`, and all BAM-level BC steps. Tested by `selftest`: HMM, inheritance/Fst/heterozygosity/depth-ratio statistics at your real depths,
region calling, sample selection, genotype-class logic, and the runner (manifest, upload, verify, restore) against a fake rclone.

## C. Parental-species analysis (YH and MC broods)
    python -m yhsd.cli parental_sex_region callset      # GATK callsets yh_broods, mc_broods
    python -m yhsd.cli parental_sex_region depth        # mosdepth windows -> male/female depth ratio
    python -m yhsd.cli parental_sex_region scan
    python -m yhsd.cli parental_sex_region summarize    # candidate regions (>= 2 independent signals) + genome-scan figures
Interpretation notes: 24 meioses per family localise an XY locus only to roughly +/-5-10 Mb at ~1 cM/Mb (simulated at your depths;
depth is not the limit, the number of meioses is). In YH the X carries the chr10 inversion, so expect a plateau across the inversion.
Fst/heterozygosity are computed within a small number of families: treat family replication as the main evidence.

## C2. Sample table: inversion genotypes (Blumer 2025 + BC fish inferred from haplotype calls)
Put Table_S1.tsv, Table_S7.tsv, Data_S1.nwk and Data_S2.nwk in one folder on utaka (e.g. `~/blumer2025/`). The BC genotype grid is the
`OutputGenotypes.csv` from the first analysis (Mconophoros_GT1 coordinates; use `--bc-legacy-mc`) or, later, `v2/OutputPosteriors.csv` /
`original/OutputGenotypes_original.csv` from the Mzebra_GT3_NCBI run (drop `--bc-legacy-mc`).
    python -m yhsd.cli sample_table inversions --blumer-dir ~/blumer2025 \
        --bc-genotypes legacy --bc-legacy-mc        # dry run, sheet untouched; downloads 2026_Nikesh/OutputGenotypes.csv via the FileManager if needed
    python -m yhsd.cli sample_table inversions --blumer-dir ~/blumer2025 --bc-genotypes legacy --bc-legacy-mc --write-sheet  # then append to the sheet
BC fish: the F1 sires are heterozygous across the chr2 and chr13 inversions, so each fish carries one F1 haplotype across the recombination-free block.
With MC non-inverted (Kumar 2025 Table 1; Blumer Table S7) the YH haplotype is the inverted one: YH haplotype -> Inversion = 1, MC haplotype -> 0.
Fish whose haplotype is not constant across the block (consistency < 0.95) are left blank. Config: `bc_inversions` in config/paper.yaml.
Pedigree expectations for BC fish (chr9 = 2, chr10 = 0, chr11 = 2) are filled too, labelled `pedigree_expected:` in InversionSource, and only if the data show
no recombination-suppressed block on that chromosome (the F1 sires are then not heterozygous). A cold block alone is not proof of an inversion: LG3 and LG7
also have cold blocks, so the chr2/chr13 inference relies on the known inversions there (Kumar 2025, Blumer 2025).
Sheet access: the code uses the FileManager's authenticated gspread client (`proj.fm.gc`) and sheet ID (`proj.fm.s_ID`) but NOT its write methods. With
`--write-sheet` it (1) appends the 7 new columns, and (2) fills only BLANK cells of the existing Inversion10 column (never overwrites). A backup of the sheet is
saved in the run folder first, and without `--write-sheet` nothing is written.
The run folder (Dropbox) holds the annotated table, join report, unmatched samples, `bc_inversion_calls.tsv`, `bc_inversion_blocks.tsv`, a full backup
of the sheet and copies of the external inputs. The writer only appends columns (BlumerSequenceID, BlumerSpecimenID, Inversion2/9/11/13,
InversionSource) and refuses to run if they exist or the target cells are not empty, so run it BEFORE adding lab PCA genotypes by hand, then fill the
remaining rows (set InversionSource to e.g. `lab_PCA`). Do not use FileManager.setSampleDatabase (it rewrites the whole sheet).
Table S7 coding: 0 non-inverted, 1 het, 2 inverted.

## C3. Sex-association tests (Fisher): is each allele / inversion a sex determiner in each clade?
Uses SampleDatabase columns: gsdf copy number (`gsdf number`), the gsdf TE insertion genotype, `Inversion9/10/11` and the LG5 genotype. Edit the column names and the
carrier definitions in `config/sex_association.yaml` (two lines are marked CHECK). For every marker x clade it builds carrier/non-carrier x expected-sex/other-sex and
reports counts, odds ratio, Fisher p, Benjamini-Hochberg q and a STATUS: 'not testable' (no variation / too few individuals), 'underpowered' (cannot reach the FDR even
with perfect association, usually too few females) or 'tested'. Lab-reared and wild samples are tested separately, offspring of known broods are left out, uncertain sexes
('M?', 'F?') and pedigree/BC-inferred inversion genotypes are excluded. A hit means association with sex, not mechanism: check relatedness in lab cohorts.
    python -m yhsd.cli sex_association fisher                                   # from the Google sheet (FileManager); uploads a run folder to Dropbox
    python -m yhsd.cli sex_association fisher --table SampleDatabase.tsv --out ./sex_association_out   # offline on any tab-separated table
    python -m yhsd.cli sex_association fisher --group-by species                # one test per species instead of per clade
Outputs: fisher_results.tsv, sample_flow.tsv (why samples were dropped), summary_matrix.tsv, sex_association_matrix.svg/png (rows in the order of the intro figure).
Intro figure script: `figures/fig_intro_sex_determination.py` (literature only; `--with-this-study` adds this-study markers).

## C4. Genomes folder, repeat annotation, inversion breakpoints (LG2, 5, 10, 13) and centromeres
All genomes of the project live in ONE folder, `<CichlidSequencingData>/YH_SD_Paper/genomes/<key>/` (and on Dropbox), registered in `config/genomes.yaml` with a manifest per genome
(accession, md5, contig stats). Linked from the FileManager (Dropbox, never copied): Mzebra_GT3_NCBI (reference), Mconophoros_GT1, YH (NCBI GCF_057382095.1, male GT1.1, inverted
homozygous on chr2/10/13, X haplotype of chr10; md5sum.txt checked) and the hifiasm chr10 haplotypes YH_X_hifiasm (h1tg000012l, 28,074,903 bp, reverse-complemented) and
YH_Y_hifiasm (h2tg000014l, 27,689,011 bp), whose md5s are in the config (for the separate X-versus-Y analysis, not used for breakpoints). Downloaded from NCBI and md5-verified:
LatrZW (ASM4766397v1) and MezeOBm (ASM4771593v1). Repeat annotations sit next to their genome in `genomes/<key>/annotation/`.
    conda install -n yh_sd_paper -c bioconda -c conda-forge minimap2         # once
    python -m yhsd.cli genomes fetch                                          # registers/downloads all genomes; lists any that still need input
    python -m yhsd.cli genomes annotate --only Mzebra_GT3_NCBI                # copies the GT3a RepeatMasker files, maps names, records version/library, uploads to Dropbox
    python -m yhsd.cli inversion_breakpoints align                            # minimap2 (asm5) of each registered assembly to the reference
    python -m yhsd.cli inversion_breakpoints call                             # breakpoints per chromosome, consensus, polarity, dotplots, BED
    python -m yhsd.cli centromeres scan                                       # satellite-density centromere calls for every genome with a registered annotation
    python -m yhsd.cli inversion_breakpoints figure                           # overview table + draft diagram (inversions and centromeres)
Annotation import (GT3a now; YH and the orange-blotch genomes when the student's RepeatMasker runs finish): fill `annotation:` of the genome in config/genomes.yaml
(`source_dir`, `files`, `main`, `lengths_from`, `compare_to: Mzebra_GT3_NCBI`) and run `genomes annotate --only <key>`. The importer maps annotation names to the genome (identity, LG label,
or unique length: scaffold01 -> NW_...), refuses the annotation if a chromosome is missing or has a different length, tolerates unmapped scaffolds and the mitochondrion (reported), and
reads RepeatMasker version, rmblastn version and library name from the .tbl, warning if they differ from GT3a (4.1.2-p1, rmblastn 2.14.1+, cichlid_merged_lib.fa).
Chromosome-to-LG labels for a genome come from `chrom_to_lg` in config/paper.yaml, else the NCBI assembly report, else the FASTA headers ("linkage group LG10"); `python -m yhsd.cli chroms` prints the table.
Inversions: LG5 from LatrZW + MezeOBm; LG2, LG10 and LG13 from the NCBI YH chromosomes. Calls do not depend on the strand of a contig: a contig and its reverse complement give the same
breakpoints. Methods (reference coordinates): A = flank / inverted segment / flank inside one contig; B = windows covered on both strands by different haplotigs (needs two haplotypes);
A-half and C = the contig ends at a breakpoint or is only the inverted segment (low confidence, flagged strand_dependent because they assume the reference orientation). `reference_is_standard`
in the config says what the inverted strand means (true for LG2/10/13, unknown for LG5). Breakpoints are intervals; inv_q_start / inv_q_end give the segment's coordinates in the assembly.
Check every `*_dotplot.png` by eye before using a coordinate. Assemblies that are not registered yet are skipped with a message.
Centromeres: satellite-repeat density per 50 kb window from the registered annotation, smoothed; one call per chromosome with a confidence class (high / medium / low / none; 'low' and
'none' mean not located) and its evidence. Non-reference genomes (YH) are lifted to the reference by the nearest alignable flanks. Cichlid centromeres are often poorly assembled: an
unassembled array is reported as not located, never guessed.

## D. BC analysis
    python -m yhsd.cli callset build bc_family                  # 7-sample GATK family VCF
    python -m yhsd.cli bc_genotyping flagstat                   # downloads BAMs, flagstat, MAPQ by chromosome
    python -m yhsd.cli bc_genotyping original-vcfs && python -m yhsd.cli bc_genotyping original
    python -m yhsd.cli bc_genotyping v2-sites                   # bcftools cross-check + extra-MC-female recount
    python -m yhsd.cli bc_genotyping v2-counts && python -m yhsd.cli bc_genotyping v2-hmm
    python -m yhsd.cli bc_genotyping v2-hmm --rho 5e-9 --tag _rho5e-9      # sensitivity (also 2e-8, 5e-8)
    python -m yhsd.cli bc_genotyping validate && python -m yhsd.cli bc_genotyping compare
Use `tmux` for long runs.

## E. Submission
Tag the submission commit on the laptop (`git tag paper-v1 && git push origin paper-v1`) and archive it to Zenodo.
Run records store the git commit: on utaka it comes from the pulled checkout (clean tree), else from DEPLOYED_FROM.txt, else VERSION.txt.
