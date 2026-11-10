# 06 Sex-linked regions from cohort VCFs (sexloci)

Specification: [`SPEC.md`](../SPEC.md), module 06 (decision 58). Packages: `sexloci` (vendored, unchanged) and `sd_regions` (the module's own command line and the older scripts).

```
python -m sd_regions sexloci list-vcfs | run-folder | run | check | describe | filter-audit | filter-compare | simulate     # sexloci 0.1.26, every option as in its README
python -m sd_regions sexloci run-folder --check-only                  # download and preflight every VCF in Outputs/PedigreeVCFs/, no analysis
python -m sd_regions sexloci run --vcf <local vcf> --n-perm 1000
python -m sd_regions parental_sex_region callset|depth|scan|summarize [--species YH,MC]      # the older scripts, until parity (SPEC module 06, Work 4)
python -m sd_regions sex_association fisher [--table T --out D]
```

Environment: `../envs/sex_mapping.yml` (`conda env create -f envs/sex_mapping.yml`, then `source tools/env.sh`). Tests: `python tools/run_tests.py 06`.

## sexloci, vendored

`sexloci/` is the package, `sexloci_dist/` the rest of sexloci 0.1.26 (README, METHODS and VARIANTS documents, CHANGELOG, configuration, examples, detection scripts and its own 105 pytest tests). **Both are copies, unchanged**:
`SEXLOCI_MANIFEST.tsv` lists the size and md5 of every file, and `tests/test_vendored_sexloci.py` fails when one changes or when a file is added. Change sexloci in its own project; to take a new release run
`python 06_sd_region_discovery/vendor_sexloci.py <sexloci-x_y_z.zip>` (it replaces both folders and rewrites the manifest), then the tests, and read `sexloci_dist/CHANGELOG.md`.

Read first: `sexloci_dist/README.md` (the rules, the outputs, the confidence tiers) and `sexloci_dist/docs/METHODS.md`. One VCF is one species is one analysis; module 05 builds the VCFs and checks a cohort against the same rules first.

Its genotyping commands (`gsdf-genotype`, `sv-genotype`, `gsdf-summary`, `write-builtin-vcf`, and the guarded sheet writer) come along unchanged but belong to module 04, which holds their parity-tested equivalents: use
`python -m sd_genotypes` for them (SPEC open decision 4 is who writes the genotype columns back).

## Tests

`tests/test_vendored_sexloci.py`: the manifest; the import and version; the command-line delegation; sexloci's simulator (an XY locus planted on chr2 is recovered at tier A by all 3 families; a wrong genome version is
refused by the preflight); and, **when pytest is installed**, sexloci's own 105 tests (skipped otherwise: they were not run when this module was built, because the sandbox had no pytest).
Not exercised, as in sexloci's own README: the scikit-allel backend, the live rclone and the live Google sheets.

## Current state

Copied in and tested as above. Not done yet (SPEC module 06, Work): parity on the positive controls (the chr10 Y inversion in the YH pedigree, gsdf dosage in *M. conophoros*, the LG5 locus in the orange-blotch
species), a configuration for bcftools VCFs, retiring the two older scripts, replacing sexloci's own sheet access (`LabData`) by the shared file manager, and a location for the results on Dropbox.
