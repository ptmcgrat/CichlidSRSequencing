import argparse, subprocess, os, urllib, shutil, contextlib, datetime
import pandas as pd
from helper_modules.file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage = 'This script will download fastq data from the SRA database and place it into the McGrath Apps Sequencing folder\n \
										  Data in the Run Info File should be Run,AvgSpotLen,Bases,BioProject,BioSample,Experiment,Instrument,LibraryName,LibraryLayout,LibrarySelect,LibrarySource,Organism,Platform,SRA Study')
parser.add_argument('Run_Info_File', type = str, help = 'File containing information on each run')

args = parser.parse_args()

fm_obj = FM()
master_sample_data = fm_obj.localSampleFile
fm_obj.downloadData(master_sample_data)

new_dt = pd.read_csv(args.Run_Info_File)
sample_dt = pd.read_csv(master_sample_data)

if len(set(new_dt.Run)) != len(new_dt):
	raise Exception('Each line of Run_Info_File should have unique Run data')

remapper = {'Run':'RunID', 'BioSample':'SampleID','AvgSpotLen':'ReadLength','Bases':'TotalBases','BioProject':'ProjectID','Instrument:':'Instrument','Library Name':'LibraryID','LibraryLayout':'LibraryLayout','LibrarySource':'LibrarySource','Organism':'Organism','Platform':'Platform'}

new_dt = new_dt.rename(columns = remapper)[remapper.values()]

new_dt['ReadGroup'] = ''
new_dt['Files'] = ''

processes = []
for index, row in new_dt.iterrows():
	#print(' Grabbing file locations for: ' + row['RunID'] + ', Time:' + str(datetime.datetime.now()))
	if not os.path.exists(fm_obj.localReadsDir + row['ProjectID']):
		os.makedirs(fm_obj.localReadsDir + row['ProjectID'])
	rg = '@RG\tID:' + row['RunID'] + '\tLB:' + row['LibraryID'] + '\tSM:' + row['SampleID'] + '\tPL:' + row['Platform']
	fqs = row['ProjectID'] + '/' + row['RunID'] + '_1.fastq.gz,,' + row['ProjectID'] + '/' + row['RunID'] + '_2.fastq.gz'

	try:
		ena_dt = pd.read_csv('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=' + row['RunID'] + '&result=read_run&fields=fastq_ftp&format=tsv&limit=0', sep = '\t')
	except:
		try:
			ena_dt = pd.read_csv('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=' + row['RunID'] + '&result=read_run&fields=fastq_ftp&format=tsv&limit=0', sep = '\t')
		except:
			ena_dt = pd.read_csv('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=' + row['RunID'] + '&result=read_run&fields=fastq_ftp&format=tsv&limit=0', sep = '\t')


	if ena_dt.fastq_ftp[0] is np.nan:
		print('<<<<<<----------Cant find data for ' + row.RunID + ': ' + row.Organism)
		continue 
	ftps = ena_dt.fastq_ftp[0].split(';')

	ena_fq1 = ftps[0]
	ena_fq2 = ftps[1]
	local_fq1 = fm_obj.localReadsDir + row['ProjectID'] + '/' + row['RunID'] + '_1.fastq.gz'
	local_fq2 = fm_obj.localReadsDir + row['ProjectID'] + '/' + row['RunID'] + '_2.fastq.gz'

	if fm_obj.checkCloudFile(local_fq1) and fm_obj.checkCloudFile(local_fq2):
		# data already downloaded
		print('Already downloaded ' + row['RunID'])
		continue

	processes.append(subprocess.Popen(['python3', 'unit_scripts/grabENA.py', row['RunID'], ena_fq1, ena_fq2, local_fq1, local_fq2]))

	row.ReadGroup = rg
	row.Files = fqs
	sample_dt = sample_dt.append(row)

	if len(processes) == 12:
		print('Waiting for processes to complete')
		for p in processes:
			p.communicate()
		sample_dt.to_csv(master_sample_data, index = False)
		fm_obj.uploadData(master_sample_data)
		print('Database uploaded')
		processes = []

if len(processes) != 0:
	for p in processes:
		p.communicate()
	print('Waiting for processes to complete')
	sample_dt.to_csv(master_sample_data, index = False)
	fm_obj.uploadData(master_sample_data)
	print('Database uploaded')
	processes = []
