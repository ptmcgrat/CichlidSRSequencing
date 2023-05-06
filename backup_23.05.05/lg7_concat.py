import os, subprocess as sp, shlex

#### Be sure to have the master file and the lg7 file in the same dir

#### Local code testing; should add 765 lines to small_master_file.vcf (26419 + 765 = 27184)
# lg7_file = '/Users/kmnike/Data/CichlidSequencingData/Outputs/vcf_sort_test/small_lg7.vcf'
# master_file = '/Users/kmnike/Data/CichlidSequencingData/Outputs/vcf_sort_test/small_master_file.vcf'

# with open(lg7_file, 'r') as f1:
#     for line in f1:
#         if not line.startswith('#'):
#             with open(master_file, 'a') as f2:
#                 f2.write(f"{line.strip()}\n")
# bcftools sort small_master_file.vcf > sorted_small_master_file.vcf

#### Run code on Mzebra
# IMPORTANT: run a word count before on master_file.vcf, whole NC_036786.1_output.vcf, and only variant lines of NC_036786.1_output.vcf:
# master_file.vcf lines: 138381918
# NC_036786.1_output.vcf all lines: 12843606
# NC_036786.1_output.vcf only variant lines: 12841864
lg7_file = '/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/original_data/NC_036786.1_output.vcf'
master_file = '/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/original_data/master_file.vcf'

with open(lg7_file, 'r') as f1:
    for line in f1:
        if not line.startswith('#'):
            with open(master_file, 'a') as f2:
                f2.write(f"{line.strip()}\n")

# bcftools sort master_file.vcf > lg1-22_master_file.vcf
# bgzip -c lg1-22_master_file.vcf > lg1-22_master_file.vcf.gz
# tabix -p vcf small_master_file.vcf.gz