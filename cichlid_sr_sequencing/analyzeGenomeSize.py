from helper_modules.file_manager import FileManager as FM
import pandas as pd 
import pdb, subprocess, sys
import gffpandas.gffpandas as gffpd
import pyfaidx

inversions = {'LG2':('NC_036781.1',19743639,20311160,43254805,43658853),'LG9':('NC_036789.1',14453796,15649299,32255605,33496468),'LG10':('NC_036790.1',11674905,11855817,29898615,29898615),
				'LG11':('NC_036791.1',8302039,8309764,30371888,30459686),'LG13':('NC_036792.1',2249541,2453698,23046928,23131968),'LG20':('NC_036799.1',19614379,19689710,32872827,33764042)}



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
genome_version = 'Mzebra_GT3'
fm_obj = FM(genome_version = genome_version)
#fm_obj.downloadData(fm_obj.localGenomeFile)

#print(fm_obj.localGenomeFile)
dna = pyfaidx.Fasta(fm_obj.localGenomeFile + '.masked')

# Problem running this line using the () in the Dropbox name
#subprocess.run(['BuildDatabase','-name',fm_obj.localGenomeDir + genome_version,fm_obj.localGenomeFile])
#subprocess.run(['RepeatModeler','-database',genome_version,'-threads','12','-LTRStruct'])
#rep_ele = gffpd.read_gff3(fm_obj.localGenomeFile + '.out.gff')
#r_dt = rep_ele.df
for lg,loc in inversions.items():
	
	#rep_left = r_dt[(r_dt.seq_id == loc[0]) & (r_dt.end > loc[1]) & (r_dt.start < loc[2])]
	dna_left = dna[loc[0]][loc[1]:loc[2]]
	dna_right = dna[loc[0]][loc[3]:loc[4]]
	if len(dna_left) + len(dna_right) < 1000000:
		with open(lg + '_left.fa', 'w') as f:
			print('>'+lg + '_Left', file = f)
			print(dna[loc[0]][loc[1]:loc[2]], file = f)
		with open(lg + '_right.fa', 'w') as f:
			print('>'+lg + '_Right', file = f)
			print(dna[loc[0]][loc[3]:loc[4]], file = f)
	

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