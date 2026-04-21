import argparse, subprocess, os, urllib, shutil, contextlib, datetime, sys, pdb
import pandas as pd
import numpy as np
from helper_modules.file_manager import FileManager as FM
from types import SimpleNamespace

# --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${TMP_DIR} -Dsamjdk.compression_level=5"

parser = argparse.ArgumentParser(usage = 'This script will download fastq data from the ENA database and place it into the McGrath Apps Sequencing folder\n \
Data is stored as uBam format to follow best practices:\n \
https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM\n \
https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently\n \
Data in the Run Info File should be Run,AvgSpotLen,Bases,BioProject,BioSample,Experiment,Instrument,Library Name,LibraryLayout,LibrarySelection,LibrarySource,Organism,Platform,SRA Study')
parser.add_argument('Run_Info_File', type = str, help = 'File containing information on each run')
parser.add_argument('-t', '--TestData', action = 'store_true', help = 'Use this flag if you want to create a small test file (1000 reads) instead of the entire read set')
parser.add_argument('-l', '--Local', action = 'store_true', help = 'Use this flag if the data is local. The Run Info File should include a FileLocations column that lists the absolute or relative path to the Reads files (split by ,,)')
args = parser.parse_args()

# Download and open master sample database file and read it in
fm_obj = FM()
fm_obj.readSampleDatabase()
# Download and open run info file that contains new data to include
fm_obj.downloadData(fm_obj.localReadDownloadDir + args.Run_Info_File)
new_dt = pd.read_csv(fm_obj.localReadDownloadDir + args.Run_Info_File)

# Basic Error checking
if len(set(new_dt.Run)) != len(new_dt):
	raise Exception('Each line of Run_Info_File should have unique Run data')

# Check all necessary data is included in the run info file
columns = ['Run','AvgSpotLen','Bases','BioProject','BioSample','Experiment','Instrument','Library Name','LibraryLayout','LibrarySelection','LibrarySource','Organism','Platform','SRA Study']
if args.Local:
	columns += ['FileLocations']

bad_data_flag = False
for c in columns:
	if c not in new_dt.columns:
		print(c + ' must be in the Run Info File')
		bad_data_flag = True
if bad_data_flag:
	sys.exit()

# Rename columns to be consistent with master sample Database
remapper = {'Run':'RunID', 'BioSample':'SampleID','AvgSpotLen':'ReadLength','Bases':'TotalBases','BioProject':'ProjectID','Instrument:':'Instrument','Library Name':'LibraryID','LibraryLayout':'LibraryLayout','LibrarySource':'LibrarySource','Organism':'Organism','Platform':'Platform'}
if args.Local:
	remapper['FileLocations'] = 'FileLocations'
new_dt = new_dt.rename(columns = remapper)[remapper.values()]
#new_dt['FileLocations'] = ''


# Loop through runs and download data and convert to uBam
commands = []
for index, row in new_dt.iterrows():
	if row['LibraryLayout'] != 'PAIRED':
		print('Error on ' + row.RunID + ': Can only handle paired end data. Library layout is: ' + layout, file = sys.stderr)
		continue
	
	data = SimpleNamespace(runID = row['RunID'], libraryID = row['LibraryID'], sampleID = row['SampleID'],
		platform = row['Platform'], layout = row['LibraryLayout'], projectID = row['ProjectID'],
		readLength = row['ReadLength'], totalBases = row['TotalBases'], instrument = row['Instrument'],
		libraryLayout = row['LibraryLayout'], librarySource = row['LibrarySource'],
		organism = row['Organism'])
	
	data.outputBamfile = fm_obj.localReadsDir + data.projectID + '/' + data.runID + '.unmapped_marked_adapters.bam'

	# Make sure we this run hasn't already been added to the sample database
	if data.runID in set(fm_obj.reads_dt['RunID']):
		print('Error on ' + data.runID + ': Run already added to sample database', file = sys.stderr)
		continue

	existing_bamfiles = set([x.split('.')[0] for x in fm_obj.returnCloudFiles(fm_obj.localReadsDir + row['ProjectID'] + '/')])
	if data.runID in existing_bamfiles:
		print('Warning on ' + data.runID + ': Run data on cloud but not in Sample Database. Rerunning...', file = sys.stderr)
		#row.FileLocations = row['ProjectID'] + '/' + run_id + '.unmapped_marked_adapters.bam'
		
	# Create directories for temp and final data to be stored in
	os.makedirs(fm_obj.localReadsDir + row['ProjectID'], exist_ok = True)
	os.makedirs(fm_obj.localTempDir, exist_ok = True)

	# Download ENA data to determine ftp site of fastq files
	if args.Local:
		fq1,fq2 = [fm_obj.localSeqCoreDataDir + x for x in row['FileLocations'].split(',,')]

	else:
		try:
			ena_dt = pd.read_csv('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=' + row['RunID'] + '&result=read_run&fields=fastq_ftp&format=tsv&limit=0', sep = '\t')
		except:
			try:
				ena_dt = pd.read_csv('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=' + row['RunID'] + '&result=read_run&fields=fastq_ftp&format=tsv&limit=0', sep = '\t')
			except:
				ena_dt = pd.read_csv('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=' + row['RunID'] + '&result=read_run&fields=fastq_ftp&format=tsv&limit=0', sep = '\t')

		# If ftp site doesn't exist it is None
		if ena_dt.fastq_ftp[0] != ena_dt.fastq_ftp[0]:
			print('Error on ' + data.RunID + ': Cant find ftp site locations', file = sys.stderr)
			continue 

		# Store file locations for remote and local fq files
		ftps = ena_dt.fastq_ftp[0].split(';')
		fq1 = ftps[0]
		fq2 = ftps[1]

	# Asynchronously download fastq files (up to 12 at a time)
	command = [str(x) for x in ['python3', 'helper_modules/grabENA.py', data.runID, fq1, fq2, data.outputBamfile, fm_obj.localTempDir, data.sampleID, data.libraryID, data.platform, data.layout]]
	#print(command)
	if args.Local:
		command += ['--Local']
	data.command = command
	data.fileLocations = data.projectID + '/' + data.runID + '.unmapped_marked_adapters.bam'
	commands.append(data)
	
for i,data in enumerate(commands):
	data.process = None
	if i < 4:
		data.process = subprocess.Popen(data.command, stderr = open(fm_obj.localTempDir + data.sampleID + '_errors.txt', 'w'))
		print('Starting analysis of ' + data.sampleID)
while commands:
	finished_processes = [x for x in commands if x.process is not None and x.process.poll() is not None]
	for data in finished_processes:
		if data.process.returncode != 0:
			print(data.sampleID + ' did not complete properly. Something went wrong')
			pdb.set_trace()
		else:
			subprocess.run(['rm',fm_obj.localTempDir + data.sampleID + '_errors.txt'])
			read_data = {'SampleID':data.sampleID,'ProjectID':data.projectID,'RunID':data.runID,
			'ReadLength':data.readLength,'TotalBases':data.totalBases,'Instrument':data.instrument,
			'LibraryID':data.libraryID,'LibraryLayout':data.libraryLayout,'LibrarySource':data.librarySource,
			'Platform':data.platform,'FileLocations':data.fileLocations}
			read_data['FileSize'] = fm_obj.returnFileSize(fm_obj.localReadsDir + data.fileLocations)
			sample_data = {'SampleID':data.sampleID,'Sex':'','Species':data.organism,'DoB':'','BroodID':'','Parents':'','Ecogroup':'','LabReared':'','Inversion10':''}
			fm_obj.addDNAReadRow(read_data, sample_data)
			print('Finished analysis of ' + data.sampleID)
		commands.remove(data)  # Remove finished process from monitoring list
		next_command = next((x for x in commands if x.process is None),None)
		if next_command is not None:
			next_command.process = subprocess.Popen(next_command.command, stderr = open(fm_obj.localTempDir + data.sampleID + '_errors.txt', 'w'))
			print('Starting analysis of ' + next_command.sampleID)


