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
    python -m yhsd selftest        # synthetic tests, no data needed; must end with ALL SELFTESTS PASSED
    python -m yhsd chroms          # RefSeq <-> LG table for Mzebra_GT3_NCBI; CHECK the LG labels
    python -m yhsd bc_genotyping provenance

### Trial on real data (writes to YH_SD_Paper_test/, never uploads)
    T="--test --chroms LG10"
    python -m yhsd callset build yh_broods --dry-run        # lists the samples that would be genotyped
    python -m yhsd $T callset build yh_broods
    python -m yhsd $T parental_sex_region scan --species YH
    python -m yhsd $T parental_sex_region summarize --species YH
Send me any traceback. Untested here (no GATK/bcftools/rclone/mosdepth in the build environment): `callset`, `depth`, VCF loading in
`scan`, and all BAM-level BC steps. Tested by `selftest`: HMM, inheritance/Fst/heterozygosity/depth-ratio statistics at your real depths,
region calling, sample selection, genotype-class logic, and the runner (manifest, upload, verify, restore) against a fake rclone.

## C. Parental-species analysis (YH and MC broods)
    python -m yhsd parental_sex_region callset      # GATK callsets yh_broods, mc_broods
    python -m yhsd parental_sex_region depth        # mosdepth windows -> male/female depth ratio
    python -m yhsd parental_sex_region scan
    python -m yhsd parental_sex_region summarize    # candidate regions (>= 2 independent signals) + genome-scan figures
Interpretation notes: 24 meioses per family localise an XY locus only to roughly +/-5-10 Mb at ~1 cM/Mb (simulated at your depths;
depth is not the limit, the number of meioses is). In YH the X carries the chr10 inversion, so expect a plateau across the inversion.
Fst/heterozygosity are computed within a small number of families: treat family replication as the main evidence.

## C2. Sample table: inversion genotypes (Blumer 2025 + BC fish inferred from haplotype calls)
Put Table_S1.tsv, Table_S7.tsv, Data_S1.nwk and Data_S2.nwk in one folder on utaka (e.g. `~/blumer2025/`). The BC genotype grid is the
`OutputGenotypes.csv` from the first analysis (Mconophoros_GT1 coordinates; use `--bc-legacy-mc`) or, later, `v2/OutputPosteriors.csv` /
`original/OutputGenotypes_original.csv` from the Mzebra_GT3_NCBI run (drop `--bc-legacy-mc`).
    python -m yhsd sample_table inversions --blumer-dir ~/blumer2025 \
        --bc-genotypes legacy --bc-legacy-mc        # dry run, sheet untouched; downloads 2026_Nikesh/OutputGenotypes.csv via the FileManager if needed
    python -m yhsd sample_table inversions --blumer-dir ~/blumer2025 --bc-genotypes legacy --bc-legacy-mc --write-sheet  # then append to the sheet
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

## D. BC analysis
    python -m yhsd callset build bc_family                  # 7-sample GATK family VCF
    python -m yhsd bc_genotyping flagstat                   # downloads BAMs, flagstat, MAPQ by chromosome
    python -m yhsd bc_genotyping original-vcfs && python -m yhsd bc_genotyping original
    python -m yhsd bc_genotyping v2-sites                   # bcftools cross-check + extra-MC-female recount
    python -m yhsd bc_genotyping v2-counts && python -m yhsd bc_genotyping v2-hmm
    python -m yhsd bc_genotyping v2-hmm --rho 5e-9 --tag _rho5e-9      # sensitivity (also 2e-8, 5e-8)
    python -m yhsd bc_genotyping validate && python -m yhsd bc_genotyping compare
Use `tmux` for long runs.

## E. Submission
Tag the submission commit on the laptop (`git tag paper-v1 && git push origin paper-v1`) and archive it to Zenodo.
Run records store the git commit: on utaka it comes from the pulled checkout (clean tree), else from DEPLOYED_FROM.txt, else VERSION.txt.
