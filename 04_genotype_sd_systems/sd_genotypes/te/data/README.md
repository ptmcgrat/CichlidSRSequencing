# Transposon-insertion definition (generated, committed)

`te_definition.json` holds, for each site of `sd_genotypes/te/sites.py` (TE1, TE2, chr7ins_c1, chr7ins_c2): the anchored record (chrom, 1-based pos, ref, alt) and the reference flanks (250 bp each side) read from
Mzebra_GT3_NCBI, the contig lengths used to check that a BAM is aligned to it, the combined-call groups, and the allele fingerprints. Build it once, on a machine that has the reference FASTA:

    python -m sd_genotypes.te.definition --ref <GCF_041146795.1_Mzebra_GT3a_genomic.fna> --out-dir sd_genotypes/te/data

It stops if the anchor base of any site is not the base of the reference at that position (wrong reference or off-by-one) or if a fingerprint differs from the pinned one. The file is part of the DEFINITION of
the measurement (`sd_genotypes.te.measure.definition_key` hashes it): do not edit by hand; rebuild and commit.
