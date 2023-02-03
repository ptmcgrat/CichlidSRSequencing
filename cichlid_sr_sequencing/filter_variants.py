import subprocess as sp, argparse, shlex

parser = argparse.ArgumentParser(usage = 'Pipeline for filtering an input VCF file by various INFO and FORMAT fields followed by sample filtering, and PCA plot generation.')
# parser.add_argument('input', help = 'Input VCF file for Filtering')
# parser.add_argument('filters', help = ' ')
args = parser.parse_args()

#### Local file paths
# vcf = '/Users/kmnike/Data/CichlidSequencingData/Outputs/small_master_file.vcf.gz'
# ref = '/Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
# out = '/Users/kmnike/Data/CichlidSequencingData/Outputs/test_out.vcf.gz'

#### Mzebra server File paths
vcf = '/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/master_file.vcf.gz'
ref = '/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
out = '/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/master_file_filtered.vcf.gz'

sp.run(shlex.split(f"gatk VariantFiltration \
   -R {ref} \
   -V {vcf} \
   -O {out} \
   --filter-name 'allele_freq' \
   --filter-expression 'AF < 0.001916' \
   --filter-name 'inbreeding_test' \
   --filter-expression 'InbreedingCoeff < -0.6' \
   --filter-name 'depth_Qual' \
   --filter-expression 'QD < 2.0' \
   --filter-name 'max_DP' \
   --filter-expression 'QD > 3000.0' \
   --filter-name 'min_DP' \
   --filter-expression 'QD < 1000.0' \
   --filter-name 'strand_bias' \
   --filter-expression 'FS > 40.0' \
   --filter-name 'mapping_quality' \
   --filter-expression 'MQ < 40.0' \
   --filter-name 'no_calls' \
   --filter-expression 'NCC > 125.0'"))

