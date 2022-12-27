import argparse, pdb, subprocess as sp, pysam
from helper_modules.file_manager import FileManager as FM
from helper_modules.Timer import Timer
import shlex
import pandas as pd

parser = argparse.ArgumentParser(usage = 'This script will download fastq data the McGrath lab dropbox and align it to the Genome version of choice')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Make sure genome version is a valid option
if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))

fm_obj.downloadData(fm_obj.localAlignmentFile)

a_dt = pd.read_csv(fm_obj.localAlignmentFile)

a_dt = a_dt[a_dt.GenomeVersion == args.Genome]

sampleIDs = set(a_dt.SampleID)

gvcffiles = []

# Create timer object
# timer = Timer()

# Download the reference genome files and build the proper directry structure for them locally.
# fm_obj.downloadData(fm_obj.localGenomeDir)
fasta_obj = pysam.FastaFile(fm_obj.localGenomeFile)
# below line defines the LG names for the first 22 LGs in the genome.
contigs = fasta_obj.references[0:22]


# Download the GVCF and GVCF.idx files for each sample in the AlignmentDatabase
# for sampleID in sampleIDs:
# 	if sampleID in []:
# 		continue
# 	print(sampleID)
# 	fm_obj.createSampleFiles(sampleID)
# 	fm_obj.downloadData(fm_obj.localGVCFFile)
# 	fm_obj.downloadData(fm_obj.localGVCFFile + '.tbi')
# 	gvcffiles.append(fm_obj.localGVCFFile)
	# pdb.set_trace() # used for troubleshooting

processes = []
processes2 = []
# First gatk command takes in chromosome names and a tab delimited cohort of samples for which to generate a genomicsdb workspace. The location of the workspace, per chromosome, must be specified using an absolute filepath.
# The loop is parallelized to run each chromosome in parallel on 4 cores.
# for contig in contigs:
#     p = sp.Popen(shlex.split(f"gatk GenomicsDBImport --genomicsdb-workspace-path {'/Data/mcgrath-lab/Data/CichlidSequencingData/Databases/' + contig + '_database'} --intervals {contig} --sample-name-map sample_map_utaka.txt --reader-threads 4"))
#     processes.append(p)

#     if len(processes) == 22:
#         for p in processes:
#             p.communicate()
#         processes = []
# This second gatk command takes in the reference genome and the path to the genomicsdb workspace and outputs all variants per chromosome per sample included in the cohort per chromosome. 
# The -V flag specifying the genomicsdb workspace location must start with 'gendb://' but it will look in the curret dir for the workspace, so a new relative path to the correct workspace must be appended to the 'gendb://'
# An output location and filename must also be provided, per chromosome. The loop is also parallelized to generate a VCF file for all samples in the cohort, per chromosome. 
for contig in contigs:
	p2 = sp.Popen(shlex.split(f"gatk GenotypeGVCFs -R /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {'gendb://../../../../../Data/mcgrath-lab/Data/CichlidSequencingData/Databases/' + contig + '_database/'}  -O {'/Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/' + contig + '_output.vcf'} --heterozygosity 0.0012"))
	processes2.append(p2)

	if len(processes2) == 22:
		for p in processes2:
			p.communicate()
		processes2=[]

print('Pipeline Completed')