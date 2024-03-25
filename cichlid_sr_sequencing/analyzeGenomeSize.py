from helper_modules.file_manager import FileManager as FM
import pandas as pd 
import pdb, subprocess, sys

def process_sample(sampleID, fm_obj):
	fm_obj.createSampleFiles(sampleID)

	for raw_bam in fm_obj.localRawBamFiles:
		fm_obj.downloadData(raw_bam)

	for raw_bam in fm_obj.localRawBamFiles:
		subprocess.run(['samtools','fastq',raw_bam], stdout = open(fm_obj.localReadsDir + sampleID + '.fq', 'a'))

	subprocess.run(['fastqc', fm_obj.localReadsDir + sampleID + '.fq'])
	subprocess.run(['jellyfish','count','-m','21','-s','6G', '-C', '-t', '96', '-U','100000', '-o', fm_obj.localReadsDir + sampleID + '_counts.jf', fm_obj.localReadsDir + sampleID + '.fq'])
	subprocess.run(['jellyfish','histo','-t','96',fm_obj.localReadsDir + sampleID + '_counts.jf'], stdout = open(fm_obj.localReadsDir + sampleID + '.histo', 'w'))
	subprocess.run(['rm','-f', fm_obj.localReadsDir + sampleID + '.fq', fm_obj.localReadsDir + sampleID + '_counts.jf'] + fm_obj.localRawBamFiles)

# Create fm_obj and grab sample file
fm_obj = FM(genome_version = 'Mzebra_HybridScaffold')
fm_obj.downloadData(fm_obj.localGenomeFile)

repeatmodelerpath = '/home/ad.gatech.edu/bio-mcgrath-dropbox/anaconda3/envs/GenomeSize/share/RepeatModeler/'

subprocess.run([repeatmodelerpath + 'BuildDatabase','-name','HybridScaffold',fm_obj.localGenomeFile])
subprocess.run([repeatmodelerpath + 'RepeatModeler','-database','HybridScaffold','-pa','12','-LTRStruct'])



#for sample in ['YH_1_m','CV_1_m','CV_2_f','MZ_1_m','MZ_2_f','MC_1_m','MC_2_f']:
#	process_sample(sample, fm_obj)


"""
YH_1_m
CV_1_m
CV_2_f
MZ_1_m
MZ_2_f
MC_1_m
MC_2_f
"""