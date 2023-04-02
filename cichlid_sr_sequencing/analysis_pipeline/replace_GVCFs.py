import os, pdb
import subprocess as sp
import pandas as pd
from helper_modules.nikesh_file_manager import FileManager as FM

fm_obj = FM("Mzebra_UMD2a")
projectIDs = ['BrainDiversity_s1', 'BigBrain']
alignment_df = pd.read_csv(fm_obj.localAlignmentFile)
filtered_df = alignment_df[alignment_df['ProjectID'].isin(projectIDs)]
sampleIDs = filtered_df['SampleID'].tolist()

local_sample_IDs = ['MC_1_m', 'SAMEA4032100']

for sample in local_sample_IDs:
    fm_obj.createSampleFiles(sample)
    sp.run(f'mv {fm_obj.testRedo_GVCFFile} {fm_obj.testGVCFFile}', shell=True)
    sp.run(f'mv {fm_obj.testRedoGVCFIndexFile} {fm_obj.testGVCFIndexFile}', shell=True)


for sample in sampleIDs:
    fm_obj.createSampleFiles(sample)
    sp.run(f'mv {fm_obj.localRedoGVCFFile} {fm_obj.localGVCFFile}', shell=True)
    sp.run(f"mv {fm_obj.localRedoGVCFFile}.tbi {fm_obj.localGVCFFile}.tbi", shell=True)


"""
/Users/kmnike/anaconda3/envs/variant/bin/python3 replace_GVCFs.py

"""
