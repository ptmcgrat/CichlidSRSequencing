# Blumer et al. supplementary tables

Input of `unit_scripts/sampledb.py` (the Blumer update of the SampleDatabase). Copied unchanged from the paper's supplementary tables.

* `Table_S1.tsv`: one row per sequenced individual (1,375): `sequence_id`, `specimen_id`, `seq_depth`, `genus`, `species`, `sex`, `clade`, `location`, `sublocation`, `Biosample`, `Bioproject`.
  The SampleDatabase's SampleID of these samples is the `Biosample` accession.
* `Table_S7.tsv`: the inversion genotypes (`chr2`, `chr9`, `chr10`, `chr11`, `chr13`; 0, 1 or 2 inverted alleles). `Table_S7.sample_id` is `Table_S1.sequence_id`.

`python3 unit_scripts/sampledb.py check | blumer [--write]` reads them from here by default (`--s1` / `--s7` take other files).
