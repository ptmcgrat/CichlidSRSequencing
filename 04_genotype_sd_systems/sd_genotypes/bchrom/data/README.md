# B-chromosome measurement definition (generated, committed)

`bchrom_blocks.tsv`, `bchrom_scale_bins.bed` and `bchrom_meta.json` are built from the discovery run and the reference by

    python -m sd_genotypes.bchrom.definition --core-blocks core_blocks.tsv --block-signal block_signal.tsv --sample-summary sample_summary.tsv \
        --lengths-from <a depth track .json or the reference .fai> --out-dir sd_genotypes/bchrom/data

(`core_blocks.tsv`, `block_signal.tsv` and `sample_summary.tsv` are in the `bchrom/call/<run>/` folder on Dropbox; the depth tracks are in `bchrom/depth_tracks_20261107/`.)
They are part of the DEFINITION of the measurement (`sd_genotypes.bchrom.measure.definition_key` hashes them): do not edit by hand; rebuild and commit.
