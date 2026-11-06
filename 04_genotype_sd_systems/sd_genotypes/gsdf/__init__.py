"""gsdf copy-number genotyping from short reads (see config/gsdf.yaml for the method).

  python -m yhsd.cli gsdf depth      stream NC_135173.1:50.70-50.83 Mb of every selected sample's BAM from Dropbox (or read a local BAM) and keep the raw bin depths
  python -m yhsd.cli gsdf call       copy number per sample (ratio -> copies), flags, junction check
  python -m yhsd.cli gsdf validate   compare with the existing `gsdf number` column, and sex x copy number per species
  python -m yhsd.cli gsdf all        depth + call (+ validate when a sample table is available)

`depth` stores the raw bins, so thresholds and calibration can change without streaming again."""
