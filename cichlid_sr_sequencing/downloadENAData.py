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
Data in the Run Info File should be run_accession,nominal_length,base_count,study_accession,sample_accession,experiment_accession,instrument_model,library_name,library_layout,library_selection,library_source,scientific_name,instrument_platform,fastq_ftp,fastq_aspera,fastq_md5')
parser.add_argument('Run_Info_File', type = str, help = 'File containing information on each run')
args = parser.parse_args()

# Download and open master sample database file and read it in
fm_obj = FM()
fm_obj.readSampleDatabase()
# Download and open run info file that contains new data to include
fm_obj.downloadData(fm_obj.localReadDownloadDir + args.Run_Info_File)
dt = pd.read_csv(fm_obj.localReadDownloadDir + args.Run_Info_File)

# Basic Error checking
if len(set(dt.run_accession)) != len(dt):
	raise Exception('Each line of Run_Info_File should have unique Run data')

# Check all necessary data is included in the run info file
#columns = ['Run','AvgSpotLen','Bases','BioProject','BioSample','Experiment','Instrument','Library Name','LibraryLayout','LibrarySelection','LibrarySource','Organism','Platform','SRA Study']
columns = ['run_accession','nominal_length','base_count','study_accession','sample_accession','experiment_accession','instrument_model','library_name','library_layout','library_selection','library_source','scientific_name','instrument_platform','fastq_ftp','fastq_aspera','fastq_md5']

missing_columns = [x for x in columns if x not in dt.columns]
if missing_columns != []:
	[print(x + ' must be in the Run Info File') for x in missing_columns]
	sys.exit()

# Rename columns to be consistent with master sample Database
remapper = {'run_accession':'RunID', 'sample_accession':'SampleID','experiment_accession':'ExperimentID','nominal_length':'ReadLength','base_count':'TotalBases','study_accession':'ProjectID','instrument_model':'Instrument','library_name':'LibraryID','library_layout':'LibraryLayout','library_source':'LibrarySource','scientific_name':'Organism','instrument_platform':'Platform','fastq_ftp':'FastqFtp','fastq_aspera':'FastqAspera','fastq_md5':'FastqMd5'}
dt = dt.rename(columns = remapper)[remapper.values()]
dt.to_csv(fm_obj.localTempDir + args.Run_Info_File)
# Loop through runs and download data and convert to uBam
commands = []
for index, row in dt.iterrows():
	if row['LibraryLayout'] != 'PAIRED':
		print('Error on ' + row.RunID + ': Can only handle paired end data. Library layout is: ' + layout, file = sys.stderr)
		continue
	
	data = SimpleNamespace(runID = row['RunID'], libraryID = row['LibraryID'], sampleID = row['SampleID'],
		platform = row['Platform'], layout = row['LibraryLayout'], projectID = row['ProjectID'],
		readLength = row['ReadLength'], totalBases = row['TotalBases'], instrument = row['Instrument'],
		libraryLayout = row['LibraryLayout'], librarySource = row['LibrarySource'],
		organism = row['Organism'], fq1 = row.FastqAspera.split(';')[0], fq2 = row.FastqAspera.split(';'))
	
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


	# Asynchronously download fastq files (up to 12 at a time)
	command = [str(x) for x in ['python3', 'helper_modules/grabENA_2.py', fm_obj.localTempDir + args.Run_Info_File, data.runID, data.outputBamfile, fm_obj.localTempDir]]
	print(command)
	data.command = command
	data.fileLocations = data.projectID + '/' + data.runID + '.unmapped_marked_adapters.bam'
	commands.append(data)
	
for i,data in enumerate(commands):
	data.process = None
	data.process = subprocess.run(data.command, capture_output = True)
	continue
	if i < 4:
		
		data.process = subprocess.Popen(data.command, stderr = open(fm_obj.localTempDir + data.sampleID + '_errors.txt', 'w'))
		print('Starting analysis of ' + data.sampleID)
"""
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
"""

