from helper_modules.file_manager import FileManager as FM
import pandas as pd 
import pdb, subprocess, sys

# Create fm_obj and grab sample file
fm_obj = FM(genome_version = 'Mzebra_UMD2a')

fm_obj.createSampleFiles('MZ_1_m')

for raw_bam in fm_obj.localRawBamFiles:
	fm_obj.downloadData(raw_bam)

for raw_bam in fm_obj.localRawBamFiles:
	subprocess.run(['samtools','fastq',raw_bam], stdout = open(fm_obj.localReadsDir + 'MZ_1_m.fq', 'a'))

"""
YH_1_m
CV_1_m
CV_2_f
MZ_1_m
MZ_2_f
MC_1_m
MC_2_f
"""