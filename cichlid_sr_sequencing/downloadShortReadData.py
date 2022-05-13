import argparse, subprocess, os, urllib, shutil, contextlib, datetime, sys, pdb
import pandas as pd
import numpy as np
from helper_modules.file_manager import FileManager as FM

# --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${TMP_DIR} -Dsamjdk.compression_level=5"

parser = argparse.ArgumentParser(usage = 'This script will download fastq data from the ENA database and place it into the McGrath Apps Sequencing folder\n \
Data is stored as uBam format to follow best practices:\n \
https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM\n \
https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently\n \
Data in the Run Info File should be Run,AvgSpotLen,Bases,BioProject,BioSample,Experiment,Instrument,LibraryName,LibraryLayout,LibrarySelect,LibrarySource,Organism,Platform,SRA Study')
parser.add_argument('Run_Info_File', type = str, help = 'File containing information on each run')

args = parser.parse_args()

# Download and open master sample database file and read it in
fm_obj = FM()
master_sample_data = fm_obj.localSampleFile
fm_obj.downloadData(master_sample_data)
sample_dt = pd.read_csv(master_sample_data)

# Download and open run info file that contains new data to include
fm_obj.downloadData(fm_obj.localReadDownloadDir + args.Run_Info_File)
new_dt = pd.read_csv(fm_obj.localReadDownloadDir + args.Run_Info_File)

# Basic Error checking
if len(set(new_dt.Run)) != len(new_dt):
	raise Exception('Each line of Run_Info_File should have unique Run data')

# Rename columns to be consistent with master sample Database
remapper = {'Run':'RunID', 'BioSample':'SampleID','AvgSpotLen':'ReadLength','Bases':'TotalBases','BioProject':'ProjectID','Instrument:':'Instrument','Library Name':'LibraryID','LibraryLayout':'LibraryLayout','LibrarySource':'LibrarySource','Organism':'Organism','Platform':'Platform'}
new_dt = new_dt.rename(columns = remapper)[remapper.values()]

# Add columns to hold info about read group info and file location
new_dt['File'] = ''

# Loop through runs and download data and convert to uBam
processes = []
for index, row in new_dt.iterrows():
	run_id, library_id, sample_id, platform, layout = row['RunID'], row['LibraryID'], row['SampleID'], row['Platform'], row['LibraryLayout']
	output_bamfile = fm_obj.localReadsDir + row['ProjectID'] + '/' + run_id + '.unmapped_marked_adapters.bam'

	# Make sure we are analyzing paired end reads
	if layout != 'PAIRED':
		print('Error on ' + row.RunID + ': Can only handle paired end data. Library layout is: ' + layout, file = sys.stderr)
		continue

	# Make sure we this run hasn't already been added to the sample database
	if run_id in sample_dt['RunID']:
		print('Error on ' + row.RunID + ': Run already added to sample database', file = sys.stderr)
		continue

	# Create directories for temp and final data to be stored in
	os.makedirs(fm_obj.localReadsDir + row['ProjectID'], exist_ok = True)
	os.makedirs(fm_obj.localTempDir, exist_ok = True)

	# Download ENA data to determine ftp site of fastq files
	try:
		ena_dt = pd.read_csv('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=' + row['RunID'] + '&result=read_run&fields=fastq_ftp&format=tsv&limit=0', sep = '\t')
	except:
		try:
			ena_dt = pd.read_csv('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=' + row['RunID'] + '&result=read_run&fields=fastq_ftp&format=tsv&limit=0', sep = '\t')
		except:
			ena_dt = pd.read_csv('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=' + row['RunID'] + '&result=read_run&fields=fastq_ftp&format=tsv&limit=0', sep = '\t')

	# If ftp site doesn't exist it is None
	if ena_dt.fastq_ftp[0] != ena_dt.fastq_ftp[0]:
		print('Error on ' + row.RunID + ': Cant find ftp site locations', file = sys.stderr)
		continue 

	# Store file locations for remote and local fq files
	ftps = ena_dt.fastq_ftp[0].split(';')
	ena_fq1 = ftps[0]
	ena_fq2 = ftps[1]

	# Asynchronously download fastq files (up to 12 at a time)
	command = [str(x) for x in ['python3', 'unit_scripts/grabENA.py', run_id, ena_fq1, ena_fq2, output_bamfile, fm_obj.localTempDir, sample_id, library_id, platform]]
	processes.append(subprocess.Popen(command))

	row.File = row['ProjectID'] + '/' + run_id + '.unmapped_marked_adapters.bam'
	sample_dt = sample_dt.append(row)

	if len(processes) == 12:
		print('  Waiting for processes to complete')
		for p in processes:
			p.communicate()
		pdb.set_trace()
		sample_dt.to_csv(master_sample_data, index = False)
		fm_obj.uploadData(master_sample_data)
		print('Database uploaded')
		processes = []

if len(processes) != 0:
	for p in processes:
		p.communicate()
	pdb.set_trace()
	print('Waiting for processes to complete')
	sample_dt.to_csv(master_sample_data, index = False)
	fm_obj.uploadData(master_sample_data)
	print('Database uploaded')
	processes = []
