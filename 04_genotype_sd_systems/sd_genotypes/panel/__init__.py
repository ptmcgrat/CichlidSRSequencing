"""panel: every per-sample measurement that needs a BAM, taken in ONE pass over the BAM while it is on local disk.

  python -m yhsd.cli panel measure --bam <sample.all.bam> --sample <ID> --out-json <ID>.panel.json [--modules gsdf]

Why one panel: the BAMs of all samples will be downloaded once more for the LG5 and B-chromosome genotypes, so everything that is measured from a BAM (gsdf copy number now; the LG5 inversion
and the B chromosome later) is defined before that, and measured together. A measurement module exposes measure_bam(bam, cfg, run) and definition_key(cfg); its raw numbers go into the JSON,
and calling (copy number, genotype, flags) is done later from the JSON by a separate, offline step, so changing a threshold never needs a BAM. The lab pipeline calls this command per sample
(a subprocess), so the logic lives only in this repository."""
