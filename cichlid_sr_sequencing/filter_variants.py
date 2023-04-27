import subprocess as sp, argparse, shlex
# /Users/kmnike/anaconda3/envs/mcgrath/bin/python /Users/kmnike/McGrath/genomics/CichlidSRSequencing/cichlid_sr_sequencing/filter_variants.py
parser = argparse.ArgumentParser(usage = 'Pipeline for filtering an input VCF file by various INFO and FORMAT fields followed by sample filtering, and PCA plot generation.')
# parser.add_argument('input', help = 'Input VCF file for Filtering')
# parser.add_argument('filters', help = ' ')
args = parser.parse_args()

### Local file paths
# vcf = '/Users/kmnike/Data/CichlidSequencingData/Outputs/NaN_troubleshooting/small_lg1-22_master_file.vcf'
# ref = '/Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
# out = '/Users/kmnike/Data/CichlidSequencingData/Outputs/NaN_troubleshooting/test_out.vcf'
# gunzip_out = '/Users/kmnike/Data/CichlidSequencingData/Outputs/NaN_troubleshooting/test_out.vcf'

# sp.run(shlex.split(f"gatk VariantFiltration \
#    -R {ref} \
#    -V {vcf} \
#    -O {out} \
#    --filter-name 'allele_freq' \
#    --filter-expression 'AF < 0.001916' \
#    --filter-name 'inbreeding_test' \
#    --filter-expression 'InbreedingCoeff < -0.6' \
#    --filter-name 'depth_Qual' \
#    --filter-expression 'QD < 2.0' \
#    --filter-name 'max_DP' \
#    --filter-expression 'DP > 12000' \
#    --filter-name 'min_DP' \
#    --filter-expression 'DP < 3000' \
#    --filter-name 'strand_bias' \
#    --filter-expression 'FS > 40.0' \
#    --filter-name 'mapping_quality' \
#    --filter-expression 'MQ < 40.0' \
#    --filter-name 'no_calls' \
#    --filter-expression 'NCC > 125' \
#    --verbosity ERROR"))
# sp.run(shlex.split(f'gunzip {out}'))

### Mzebra server File paths
vcf = '/Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/vcf_concat_output/576_cohort_master_file.vcf'
ref = '/Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
out = '/Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/vcf_concat_output/filtered_576_cohort_master_file.vcf'
# gunzip_out = '/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/master_file_filtered.vcf'
sp.run(shlex.split(f"gatk VariantFiltration \
   -R {ref} \
   -V {vcf} \
   -O {out} \
   --filter-name 'allele_freq' \
   --filter-expression 'AF < 0.000958' \
   --filter-name 'inbreeding_test' \
   --filter-expression 'InbreedingCoeff < -0.6' \
   --filter-name 'depth_Qual' \
   --filter-expression 'QD < 2.0' \
   --filter-name 'max_DP' \
   --filter-expression 'DP > 11000' \
   --filter-name 'min_DP' \
   --filter-expression 'DP < 8000' \
   --filter-name 'strand_bias' \
   --filter-expression 'FS > 40.0' \
   --filter-name 'mapping_quality' \
   --filter-expression 'MQ < 50.0' \
   --filter-name 'no_calls' \
   --filter-expression 'NCC > 125.0' \
   --verbosity ERROR"))

# sp.run(shlex.split(f'gunzip {out}'))

#### test local filter file for number of filters called in each variant
# filters = {'allele_freq':0, 'inbreeding_test':0, 'depth_Qual':0, 'max_DP':0, 'min_DP':0, 'strand_bias':0, 'mapping_quality':0, 'no_calls':0, 'PASS':0}

# with open(gunzip_out, 'r') as fh:
#    for line in fh:
#       if line.strip().split()[0].startswith('##') or line.strip().split()[0].startswith('#'):
#          continue
#       filter = line.split('\t')[6]
#       filter = filter.split(';')
#       for criteria in filter:
#          filters[criteria] += 1
# print(filters)



#### see how many filters were triggered for variants in the whole file on the server

# filters = {'allele_freq':0, 'inbreeding_test':0, 'depth_Qual':0, 'max_DP':0, 'min_DP':0, 'strand_bias':0, 'mapping_quality':0, 'no_calls':0, 'PASS':0}
# with open(gunzip_out, 'r') as fh:
#    for line in fh:
#       if line.strip().split()[0].startswith('##') or line.strip().split()[0].startswith('#'):
#          continue
#       filter = line.split('\t')[6]
#       filter = filter.split(';')
#       for criteria in filter:
#          filters[criteria] += 1
# print(filters)

# gatk VariantFiltration -R ~/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V small_lg1-22_master_file.vcf.gz -O pass.vcf --filter-name 'allele_freq' --filter-expression 'AF < 0.000958' --filter-expression 'DP > 11000' --filter-name 'min_DP'  --filter-not-in-mask