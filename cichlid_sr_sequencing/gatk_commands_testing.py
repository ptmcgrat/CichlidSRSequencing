import argparse
import subprocess as sp
import shlex
import pysam
import pandas as pd
from helper_modules.file_manager import FileManager as FM


parser = argparse.ArgumentParser(usage = 'This script will download fastq data the McGrath lab dropbox and align it to the Genome version of choice')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
args = parser.parse_args()

# call the File Manager Object
fm_obj = FM(args.Genome)
# Use pysam to store the contig names into fasta_obj
fasta_obj = pysam.FastaFile(fm_obj.localGenomeFile)
# read in the AlignmentFile.csv file conatiing all sample names for which bam files were generated and store as a pandas dataframe
a_dt = pd.read_csv(fm_obj.localAlignmentFile)
# rename the genomeversion to match the reference used in this analysis
a_dt = a_dt[a_dt.GenomeVersion == args.Genome]
# Use the dataframe to define a list of the sample names as sampleIDs
sampleIDs = set(a_dt.SampleID)
# define the contigs we want to call variants for
contigs = fasta_obj.references[0:2]

processes = []
processes2 = []
# for contig in contigs:
# 	print(shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk GenomicsDBImport --genomicsdb-workspace-path {'/Users/kmnike/Data/CichlidSequencingData/Databases' + contig + '_database'} --intervals {contig} --sample-name-map test_sample_map.txt --reader-threads 4"))
# UNCOMMENT BELOW FOR LOOPS AND TIME THE RUN FOR 2 SAMPLES ON 2 LGS USING PARALLELIZATION AND WITHOUT PARALLELIZATION. CREATE THE DATABASES DIR AND DIRECT THE DATABASES THERE FOR NOW WHILE EDITING THE FILE MANAGER SCRIPT IS FIGURED OUT
# for contig in contigs:
#     # Change the absolute path to gatk to simply 'gatk' on the server
#     p = sp.Popen(shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk GenomicsDBImport --genomicsdb-workspace-path {'/Users/kmnike/Data/CichlidSequencingData/Databases/' + contig + '_database'} --intervals {contig} --sample-name-map test_sample_map.txt --reader-threads 4"))
#     processes.append(p)

#     if len(processes) == 2:
#         for p in processes:
#             p.communicate
#         processes = []


for contig in contigs:
	# print(shlex.split(f"gatk GenotypeGVCFs -R /Users/kmnike/Data/CichlidSequencingData/Genome/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {'gendb://Users/kmnike/Data/CichlidSequencingData/Databases' + contig + '_database'}  -O {contig + '_output.vcf'}"))
	p2 = sp.Popen(shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk GenotypeGVCFs -R /Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {'gendb://../../../../../Data/CichlidSequencingData/Databases/' + contig + '_database'}  -O {'/Users/kmnike/Data/CichlidSequencingData/Outputs/' + contig + '_output.vcf'}"))
	processes2.append(p2)

	if len(processes2) == 2:
		for p in processes2:
			p.communicate
		processes2=[]

"""
# command testing for gendb workspace
gatk GenotypeGVCFs -R /Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V /Users/kmnike/Data/CichlidSequencingData/Databases/NC_036780.1_database  -O /Users/kmnike/Data/CichlidSequencingData/Outputs/NC_036780_output.vcf
gatk GenotypeGVCFs -R /Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V /Users/kmnike/Data/CichlidSequencingData/Databases/NC_036780.1_database/callset.json -O /Users/kmnike/Data/CichlidSequencingData/Outputs/NC_036780_output.vcf
gatk GenotypeGVCFs -R /Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V gendb://../../../../../Data/CichlidSequencingData/Databases/NC_036780.1_database/  -O /Users/kmnike/Data/CichlidSequencingData/Outputs/NC_036780_output.vcf
-V gendb://test
"""

# command2 = shlex.split('gatk GenotypeGVCFs -R /Users/kmnike/Data/CichlidSequencingData/Genome/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V gendb://my_database -O test_output.vcf')
# for contig in contigs:
# 	print(shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk GenomicsDBImport --genomicsdb-workspace-path {'/Users/kmnike/Data/CichlidSequencingData/Databases' + contig + '_database'} --intervals {contig} --sample-name-map test_sample_map.txt --reader-threads 4"))
	# sp.Popen(shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk GenomicsDBImport --genomicsdb-workspace-path {'/Users/kmnike/Data/CichlidSequencingData/Databases/' + contig + '_database'} --intervals {contig} --sample-name-map test_sample_map.txt --reader-threads 4"))
	
# Possibly Useful Code:
# command1 = shlex.split('gatk GenomicsDBImport --genomicsdb-workspace-path my_database_{contigs} --intervals {contigs} --sample-name-map test_sample_map.txt --reader-threads 4')
# for contig in contigs:
#     print(shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk GenomicsDBImport --genomicsdb-workspace-path {contig + '_database'} --intervals {contig} --sample-name-map test_sample_map.txt --reader-threads 4"))
"""
commands = [shlex.split('python3 timer.py 3'), shlex.split('python3 timer.py 6'), shlex.split('python3 timer.py 2')]
processes = []
for command in commands:
	p = sp.Popen(command)
	processes.append(p)

	if len(processes) == 3:
		for p in processes:
			p.communicate()
		processes = []
"""