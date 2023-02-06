import subprocess as sp, argparse, shlex

parser = argparse.ArgumentParser(usage = 'Pipeline for filtering an input VCF file by various INFO and FORMAT fields followed by sample filtering, and PCA plot generation.')
# parser.add_argument('input', help = 'Input VCF file for Filtering')
# parser.add_argument('filters', help = ' ')
args = parser.parse_args()

#### Local file paths
# vcf = '/Users/kmnike/Data/CichlidSequencingData/Outputs/small_master_file.vcf.gz'
# ref = '/Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
# out = '/Users/kmnike/Data/CichlidSequencingData/Outputs/test_out.vcf.gz'

# sp.run(shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk VariantFiltration \
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
#    --filter-expression 'DP > 3000' \
#    --filter-name 'min_DP' \
#    --filter-expression 'DP < 1000' \
#    --filter-name 'strand_bias' \
#    --filter-expression 'FS > 40.0' \
#    --filter-name 'mapping_quality' \
#    --filter-expression 'MQ < 40.0' \
#    --filter-name 'no_calls' \
#    --filter-expression 'NCC > 125.0' \
#    --verbosity ERROR"))
# sp.run(shlex.split(f'gunzip {out}'))

### Mzebra server File paths
vcf = '/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/original_data/master_file.vcf.gz'
ref = '/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
out = '/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/master_file_correctly_filtered.vcf.gz'

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
   --filter-expression 'DP > 3000.0' \
   --filter-name 'min_DP' \
   --filter-expression 'DP < 1000.0' \
   --filter-name 'strand_bias' \
   --filter-expression 'FS > 40.0' \
   --filter-name 'mapping_quality' \
   --filter-expression 'MQ < 40.0' \
   --filter-name 'no_calls' \
   --filter-expression 'NCC > 125.0' \
   --verbosity ERROR"))
sp.run(shlex.split(f'gunzip {out}'))

#### test local filter file for number of filters called in each variant
# filters = {'allele_freq':0, 'inbreeding_test':0, 'depth_Qual':0, 'max_DP':0, 'min_DP':0, 'strand_bias':0, 'mapping_quality':0, 'no_calls':0, 'PASS':0}

# with open('/Users/kmnike/Data/CichlidSequencingData/Outputs/test_out.vcf', 'r') as fh:
#    for line in fh:
#       if line.strip().split()[0].startswith('##') or line.strip().split()[0].startswith('#'):
#          continue
#       filter = line.split('\t')[6]
#       filter = filter.split(';')
#       for criteria in filter:
#          filters[criteria] += 1
# print(filters)



#### see how many filters were triggered for variants in the whole file. filters = {'allele_freq':0, 'inbreeding_test':0, 'depth_Qual':0, 'max_DP':0, 'min_DP':0, 'strand_bias':0, 'mapping_quality':0, 'no_calls':0}

filters = {'allele_freq':0, 'inbreeding_test':0, 'depth_Qual':0, 'max_DP':0, 'min_DP':0, 'strand_bias':0, 'mapping_quality':0, 'no_calls':0, 'PASS':0}
with open('/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/master_file_correctly_filtered.vcf', 'r') as fh:
   for line in fh:
      if line.strip().split()[0].startswith('##') or line.strip().split()[0].startswith('#'):
         continue
      filter = line.split('\t')[6]
      filter = filter.split(';')
      for criteria in filter:
         filters[criteria] += 1
   
print(filters)