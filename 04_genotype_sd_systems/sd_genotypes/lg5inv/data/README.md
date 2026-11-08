# LG5 inversion definition (generated, committed)

`lg5_markers.tsv` (chrom, pos, ref, alt) and `lg5_junctions.json` (the junction probe sequences) are built from the pool analysis and the reference by

    python -m sd_genotypes.lg5inv.definition --ref <Mzebra_GT3_NCBI .fna> --markers <markers.tsv from `markers select`> --out-dir sd_genotypes/lg5inv/data

They are part of the DEFINITION of the measurement (`sd_genotypes.lg5inv.measure.definition_key` hashes them): do not edit by hand; rebuild and commit.
