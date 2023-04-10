import subprocess as sp
import pandas as pd
from helper_modules.nikesh_file_manager import FileManager as FM
import pathlib, pdb

"""
The script will be used to clean up the various files in the Bamfiles folder on Utaka. These folders should only contain the *.g.vcf.gz files and the *.g.vcf.gz.tbi files.
All other files need to be removed. 
Looks like rclone will upload the new version and overwrite the previous if the destination dir cpntains that original file 

For the error files:
- A new GVCF file has been regenerated. This will need to be uploaded to replace the old one on Dropbox.
- The "rerun" files and their indexes need to be removed

For the BrainDiversity_S1 & BigBrain Cohort:
- the *.g.vcf.gz & *.g.vcf.gz.tbi files need to be reuploaded to dropbox. 

All Directories:
- for all Directories, the final code to run will be to remove any bam and bai files 


"""
fm_obj = FM('Mzebra_UMD2a')
projectIDs = ['BrainDiversity_s1', 'BigBrain']
alignment_df = pd.read_csv(fm_obj.localAlignmentFile)
filtered_df = alignment_df[alignment_df['ProjectID'].isin(projectIDs)]
sampleIDs = filtered_df['SampleID'].tolist()
error_files = ['SAMEA4033321', 'SAMEA2661241', 'SAMEA4032067', 'SAMEA4033318', 'SAMEA4032108', 'SAMEA4033341', 'SAMEA4033314', 'SAMEA4032129', 'SAMEA4033252', 'SAMEA4032105', 'SAMEA3388855', 'SAMEA4032131', 'SAMEA4033248', 'SAMEA2661359', 'SAMEA1877414', 'SAMN06685761', 'SAMEA4032096', 'SAMEA2661294', 'SAMEA4032053', 'SAMEA1920091', 'SAMEA4032046', 'SAMEA4033254', 'SAMEA2661396', 'SAMEA4033283', 'SAMEA4032136', 'SAMEA2661367', 'SAMEA3388862', 'SAMEA4033301', 'SAMEA2661292', 'SAMEA4033271', 'SAMEA4032127', 'SAMEA2661414', 'SAMEA4033317', 'SAMEA2661221', 'SAMEA2661310', 'SAMEA1904322', 'SAMEA4032051', 'SAMEA4032033', 'SAMEA4033277', 'SAMEA4033307', 'SAMEA4032038', 'SAMEA4032042', 'SAMEA2661239', 'SAMN08051112', 'SAMEA4032125', 'SAMEA2661258', 'SAMEA2661306', 'SAMEA4033322', 'SAMEA1920095', 'SAMEA4033284', 'SAMEA4033304', 'SAMEA4033305', 'SAMEA4032137', 'SAMEA1920096', 'SAMEA3388871', 'SAMEA4033278', 'SAMEA4033286', 'SAMEA4032049', 'SAMEA1920092', 'SAMEA4033289', 'SAMEA4033324', 'SAMEA1904832', 'SAMEA2661339', 'SAMEA1877499', 'SAMEA2661277', 'SAMEA2661248', 'SAMEA2661280', 'SAMEA2661381', 'SAMEA1904328', 'SAMEA4033261', 'SAMEA4033287', 'SAMEA4032088', 'SAMEA4032070', 'SAMEA4032069']
test_samples = ['CK_1003_p']

def upload_data(samples):
    for sample in samples:
        fm_obj.createSampleFiles(sample)
        fm_obj.uploadData(fm_obj.localGVCFFile)

def removeExtraFiles(samples):
    for sample in samples:
        fm_obj.createSampleFiles(sample)
        extra_files = [fm_obj.localBamFile]

# upload_data(error_files)
removeExtraFiles(error_files)

print('DONE')



"""
LocalTesting:

/Users/kmnike/anaconda3/envs/pipeline/bin/python3 upload_new_GVCFs_and_BAMs.py

"""