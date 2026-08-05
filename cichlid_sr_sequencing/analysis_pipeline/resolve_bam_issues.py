import pdb
import pandas as pd
import subprocess as sp


temp = sp.run(['ls', '/Data/mcgrath-lab/Temp/CichlidSequencingData/Temp'], capture_output=True, encoding='utf-8')
temp_samples = temp.stdout.split('\n')[:-1]
mc_bam_samples = sp.run(['ls', '/Data/mcgrath-lab/Temp/CichlidSequencingData/Bamfiles/Mconophoros_GT1'], capture_output=True, encoding='utf-8').stdout.split('\n')[:-1]
s_df = pd.read_csv('pm262728_run_samples.csv')
seq_samples = s_df.samplename.to_list()

# get samples that did not get .all.bam files written:
missing_bam_samples = []
with open('/home/mcgrath-lab/Patrick/CichlidSRSequencing/cichlid_sr_sequencing/nikesh_logs/error_retry_failed_lc_samples.txt', 'r' ) as f:
    for line in f:
        if line.startswith('['):
            missing_bam_samples.append(line.split('/')[-2])


hc = ['MCYHBC1-564-1', 'MCYHBC1-579-1', 'MCYHBC1-602-1', 'MCYHBC1-614-1', 'MCYHBC1-629-1', 'MCYHBC1-642-1', 'MCYHBC1-645-1', 'MCYHBC1-655-1', 'MCYHBC1-661-1', 'MCYHBC1-737-1', 'MCYHBC1-804-1', 'MCYHBC1-809-1', 'MCYHBC1-825-1', 'MCYHBC1-834-1', 'MCYHBC1-840-1', 'MCYHBC1-886-1', 'MCYHBC1-894-1', 'MCYHBC1-905-1']

hc_dict = {}
for sample in hc:
    # Code figures out the size of the *.all.bam files created for the high coverage samples located in /Data/mcgrath-lab/Temp/CichlidSequencingData/Bamfiles/Mconophoros_GT1
    # Can confirm that the *.all.bam files were not rewritten and they are large compared to bam files for the LC samples.
    bam_file_size = sp.check_output(['rclone', 'size', '/Data/mcgrath-lab/Temp/CichlidSequencingData/Bamfiles/Mconophoros_GT1/' + sample + "_HC", '--include', '*.all.bam'], encoding='utf-8')
    file_size = bam_file_size.strip().split(':')[2].split(' ')[1]
    hc_dict[sample] = file_size
# Print out the HC samples's *.all.bam file size
# for key in hc_dict:
#     print(f"{key}\t{hc_dict[key]}")


# for key in hc_dict:
#     # Change names for firectories of high coverage BC1 samples for now so that they can be rerun with their corresponding low coverage data.
#     # Note that only 17 of the 18 HC samples were incldued for LC sequencing. Sample MCYHBC1-840-1 was not included for LC sequencing due to an error on my part. 
#     sp.run(['mkdir', '/home/mcgrath-lab/nkumar317/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/temp/' + key])
#     sp.run(['mv', '/Data/mcgrath-lab/Temp/CichlidSequencingData/Temp/' + key, '/Data/mcgrath-lab/Temp/CichlidSequencingData/Temp/' + key + '_HC'])
#     sp.run(['mv', '/Data/mcgrath-lab/Temp/CichlidSequencingData/Bamfiles/Mconophoros_GT1/' + key, '/Data/mcgrath-lab/Temp/CichlidSequencingData/Bamfiles/Mconophoros_GT1/' + key + '_HC'])
#     pass

# confirm that no HC sample was in the missing sample list
for sample in hc:
    if sample in missing_bam_samples:
        print(sample + ' is also in the missing samples')
    else:
        print(sample + " is not in the missing bam samples")

# concat the missing samples and hc sample lists since the LC samples need to be run alongside 
all = hc + missing_bam_samples
all_samples_to_rerun = []
for sample in all:
    if sample == 'MCYHBC1-840-1':
        continue
    else:
        all_samples_to_rerun.append(sample)
pdb.set_trace()


"""
MISSING: 
MISSING: 
MISSING: 
MISSING: 
MISSING: 
MISSING: 
MISSING: 
MISSING: 
MISSING: 
MISSING: 
MISSING: 
MISSING: 
MISSING: 
"""

# MCYHBC1-673-1 MCYHBC1-704-1 MCYHBC1-729-1 MCYHBC1-752-1 MCYHBC1-761-1 MCYHBC1-765-1 MCYHBC1-796-1 MCYHBC1-807-1 MCYHBC1-817-1 MCYHBC1-823-1 MCYHBC1-843-1 MCYHBC1-931-1 MCYHBC1-784-1

