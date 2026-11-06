"""lg5inv: reads for the LG5 (orange blotch / ZW) inversion, from pools of known composition.

  python -m sd_genotypes.lg5inv prepare [--only RUN,...] [--runs-in-flight 1] [--threads 12]   download -> align -> mark duplicates -> keep the LG5 inversion region and the unmapped reads
  python -m sd_genotypes.lg5inv status                                                           what is done, what is on Dropbox, disk space
  python -m sd_genotypes.lg5inv clean [--force]                                                  delete the work folder once every finished run is on Dropbox

Independent of bchrom (its own config/lg5inv.yaml and work folder), but it reuses bchrom's tested per-run processing (the lab alignment recipe, SRA/ENA download, disk ledger, checkpoints).
Each finished pool is uploaded as its own run folder `lg5inv/pool_<run>/` (summary, depth track, slice BAMs) as soon as it is done, so a crash loses at most the run in flight.
Marker discovery and genotyping are separate steps that read these slices."""
