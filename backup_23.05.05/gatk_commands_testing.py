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
# First gatk command takes in chromosome names and a tab delimited cohort of samples for which to generate a genomicsdb workspace. The location of the workspace, per chromosome, must be specified using an absolute filepath. 
# The loop is parallelized to run each chromosome in parallel on 4 cores.
for contig in contigs:
    p = sp.Popen(shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk GenomicsDBImport --genomicsdb-workspace-path {'/Users/kmnike/Data/CichlidSequencingData/Databases/' + contig + '_database'} --intervals {contig} --sample-name-map test_sample_map.txt --reader-threads 4"))
    processes.append(p)

    if len(processes) == 2:
        for p in processes:
            p.communicate()
        processes = []
# This second gatk command takes in the reference genome and the path to the genomicsdb workspace and outputs all variants per chromosome per sample included in the cohort per chromosome. 
# The -V flag specifying the genomicsdb workspace location must start with 'gendb://' but it will look in the curret dir for the workspace, so a new relative path to the correct workspace must be appended to the 'gendb://'
# An output location and filename must also be provided, per chromosome. The loop is also parallelized to generate a VCF file for all samples in the cohort, per chromosome. 
for contig in contigs:
	p2 = sp.Popen(shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk GenotypeGVCFs -R /Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {'gendb://../../../../../Data/CichlidSequencingData/Databases/' + contig + '_database'}  -O {'/Users/kmnike/Data/CichlidSequencingData/Outputs/' + contig + '_output.vcf'}"))
	processes2.append(p2)

	if len(processes2) == 2:
		for p in processes2:
			p.communicate()
		processes2=[]







"""
Possibly Useful Code:
command1 = shlex.split('gatk GenomicsDBImport --genomicsdb-workspace-path my_database_{contigs} --intervals {contigs} --sample-name-map test_sample_map.txt --reader-threads 4')
for contig in contigs:
    print(shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk GenomicsDBImport --genomicsdb-workspace-path {contig + '_database'} --intervals {contig} --sample-name-map test_sample_map.txt --reader-threads 4"))

command testing for gendb workspace:
gatk GenotypeGVCFs -R /Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V /Users/kmnike/Data/CichlidSequencingData/Databases/NC_036780.1_database  -O /Users/kmnike/Data/CichlidSequencingData/Outputs/NC_036780_output.vcf
gatk GenotypeGVCFs -R /Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V /Users/kmnike/Data/CichlidSequencingData/Databases/NC_036780.1_database/callset.json -O /Users/kmnike/Data/CichlidSequencingData/Outputs/NC_036780_output.vcf
gatk GenotypeGVCFs -R /Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V gendb://../../../../../Data/CichlidSequencingData/Databases/NC_036780.1_database/  -O /Users/kmnike/Data/CichlidSequencingData/Outputs/NC_036780_output.vcf
-V gendb://test


# for contig in contigs:
# 	print(shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk GenomicsDBImport --genomicsdb-workspace-path {'/Users/kmnike/Data/CichlidSequencingData/Databases' + contig + '_database'} --intervals {contig} --sample-name-map test_sample_map.txt --reader-threads 4"))

"""