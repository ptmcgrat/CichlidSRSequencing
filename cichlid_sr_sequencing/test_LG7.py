# for testing on the mzebra server only
import subprocess as sp, argparse, pandas as pd, shlex, pysam, os
from helper_modules.file_manager import FileManager as FM


parser = argparse.ArgumentParser(usage = 'Testing gatk GenomicsDBImport parallelization optimization')
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

lg7 = 'NC_036786.1'
# sp.run(shlex.split(f"gatk GenomicsDBImport --genomicsdb-workspace-path {'/Data/mcgrath-lab/Data/CichlidSequencingData/TestingDatabases/' + lg7 + '_database'} --intervals lg7.interval_list --sample-name-map sample_map_utaka.txt --max-num-intervals-to-import-in-parallel 4"))

processes = []
processes2 = []
#### First gatk command takes in chromosome names and a tab delimited cohort of samples for which to generate a genomicsdb workspace. The location of the workspace, per chromosome, must be specified using an absolute filepath.
#### The loop is parallelized to run each chromosome in parallel on 4 cores.
#### Below is the new code to use after splitting the contigs into intervals and running by importing 4 intervals at a time. 100kbp took about 2.41 mins
# for contig in test_contigs:
#     p = sp.Popen(shlex.split(f"gatk GenomicsDBImport --genomicsdb-workspace-path {'/Data/mcgrath-lab/Data/CichlidSequencingData/TestingDatabases/' + contig + '_database'} --intervals small_contig.interval_list --sample-name-map sample_map_utaka.txt --max-num-intervals-to-import-in-parallel 4"))
#     processes.append(p)

#     if len(processes) == 3:
#         for p in processes:
#             p.communicate()
#         processes = []


#### here's the orignal code. 
# for contig in test_contigs:
#     p = sp.Popen(shlex.split(f"gatk GenomicsDBImport --genomicsdb-workspace-path {'/Data/mcgrath-lab/Data/CichlidSequencingData/TestingDatabases/' + contig + '_database2'} --intervals small_contig.interval_list --sample-name-map sample_map_utaka.txt --reader-threads 4"))
#     processes.append(p)

#     if len(processes) == 3:
#         for p in processes:
#             p.communicate()
#         processes = []

#### This second gatk command takes in the reference genome and the path to the genomicsdb workspace and outputs all variants per chromosome per sample included in the cohort per chromosome. 
#### The -V flag specifying the genomicsdb workspace location must start with 'gendb://' but it will look in the curret dir for the workspace, so a new relative path to the correct workspace must be appended to the 'gendb://'
#### An output location and filename must also be provided, per chromosome. The loop is also parallelized to generate a VCF file for all samples in the cohort, per chromosome.
# for contig in test_contigs:
# 	p2 = sp.Popen(shlex.split(f"gatk GenotypeGVCFs -R /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {'gendb://../../../../../Data/mcgrath-lab/Data/CichlidSequencingData/TestingDatabases/' + contig + '_database/'}  -O {'/Data/mcgrath-lab/Data/CichlidSequencingData/TestingOutputs/' + contig + '_output.vcf'} --heterozygosity 0.0012"))
# 	processes2.append(p2)

# 	if len(processes2) == 3:
# 		for p in processes2:
# 			p.communicate()
# 		processes2=[]

print('Pipeline Completed')
# START BY ADDING THE CODE TO DOWNLOAD THE DATA TO AUTOMATICALLY CREATE THE DIR STRUCTURE 