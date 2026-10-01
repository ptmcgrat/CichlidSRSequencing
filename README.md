# CichlidSRSequencing (Kumar_eLife branch)

This branch holds the code used for the short-read (Illumina) analyses in our eLife paper:

> Kumar et al. **Large inversions in Lake Malawi cichlids are associated with habitat preference, lineage, and sex determination.** *eLife* (2025). https://doi.org/10.7554/eLife.104923

The branch was trimmed down from the main repo so it only has the scripts we used for the paper. The Bionano optical mapping and the PacBio HiFi genome assemblies were done outside of this repo. The code here picks up once we have the genomes and the Illumina reads, and covers aligning reads to the M_zebra_GT3a reference, calling variants, genotyping the inversions with PCA, building phylogenies, calculating genetic distance, testing for sex association, and whole-genome alignments.

## A quick heads-up before running anything

These scripts were written to run on our lab server (Utaka), and they pull data from and push data to the lab Dropbox using `rclone`. All of that is handled by the `FileManager` class in `helper_modules/file_manager.py`. Because of this, they won't run out of the box on another machine without changing paths and the rclone remote. We're posting them mainly so you can see exactly what we ran and with what parameters.

Some other things to know:

- Sample metadata (ecogroup, organism, sex, which samples go in which analysis, etc.) is read from `SampleDatabase_v2.xlsx` (the `SampleLevel` sheet). Columns like `CorePCA`, `PCAFigure`, `PhylogenyFigure`, `BionanoData`, and `YHPedigree` are Yes/No flags that decide which samples go into each analysis. See the paper's supplementary tables for the sample list and accessions.
- The GT3a reference keeps the old UMD2a chromosome names (e.g. `NC_036780.1` = LG1), so the linkage group names hardcoded in the scripts are the NCBI UMD2a names.
- Most scripts are run from inside the `cichlid_sr_sequencing/` directory, since they import from `helper_modules/` and look for some files relative to where you run them.

## Pipeline overview

The scripts are run in this order:

| Script | What it does |
|--------|--------------|
| `downloadShortReadData.py` | Generates uBAM files from FASTQ per GATK guidelines before alignment |
| `alignFastQ.py` | Aligns Illumina reads to M_zebra_GT3a and makes a GVCF for each sample |
| `callVariants.py` | Joint genotypes all samples into a single cohort VCF |
| `process_vcf.py` | Concatenates VCFs and performs compression and filtering to generate the pass_variants VCF file |
| `pca_maker.py` | Runs PCA on the whole genome and on each inversion region to genotype the inversions |
| `buildPhylogenies.py` | Builds maximum likelihood trees for the whole genome and each inversion |
| `analyzePiXY.py` | Calculates pairwise genetic distance (dXY) between ecogroups for each inversion |
| `analyzeSexAssociation.py` | Scans for Fst and heterozygosity differences between males and females in the lab broods |
| `alignGenomes.py` | Runs whole-genome alignments of the new assemblies and the outgroups to GT3a |

## Dependencies

All dependencies were installed and managed using conda. Please refer to the Methods of the *eLife* paper for details on software versions.

## Citation

If you use any of this code, please cite:

> Kumar et al. Large inversions in Lake Malawi cichlids are associated with habitat, lineage, and sex determination. *eLife* (2025). https://doi.org/10.7554/eLife.104923

## License

This code is released under the MIT License. See [LICENSE](LICENSE) for details.
