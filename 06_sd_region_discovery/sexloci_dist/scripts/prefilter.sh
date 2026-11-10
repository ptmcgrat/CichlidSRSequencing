#!/usr/bin/env bash
# Optional: thin a large cohort VCF (biallelic SNPs on chromosome-scale contigs) before running sexloci.
# usage: prefilter.sh in.vcf.gz out.vcf.gz "chr1,chr2,..."
set -euo pipefail
bcftools view -r "$3" -m2 -M2 -v snps "$1" -Oz -o "$2" && tabix -p vcf "$2"
