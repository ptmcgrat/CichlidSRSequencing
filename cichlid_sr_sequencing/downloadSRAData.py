import argparse, subprocess, os
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

for index, row in new_dt.iterrows():
	if not os.path.exists(fm_obj.localReadsDir + row['ProjectID']):
		os.makedirs(fm_obj.localReadsDir + row['ProjectID'])
	rg = '@RG\tID:' + row['RunID'] + '\tLB:' + row['LibraryID'] + '\tSM:' + row['SampleID'] + '\tPL:' + row['Platform']
	subprocess.run(['prefetch', row['RunID']])
	subprocess.run(['fastq-dump', os.getenv('HOME') + '/ncbi/public/sra/' + row['RunID'] + '.sra', '--split-files'])
	fqs = row['ProjectID'] + '/' + row['RunID'] + '_1.fastq,,' + row['ProjectID'] + '/' + row['RunID'] + '_2.fastq'
	subprocess.run(['mv', row['RunID'] + '_1.fastq', row['RunID'] + '_2.fastq', fm_obj.localReadsDir + row['ProjectID']])
	fm_obj.uploadData(fm_obj.localReadsDir + row['ProjectID'] + '/' + row['RunID'] + '_1.fastq')
	fm_obj.uploadData(fm_obj.localReadsDir + row['ProjectID'] + '/' + row['RunID'] + '_2.fastq')

	subprocess.run(['rm', os.getenv('HOME') + '/ncbi/public/sra/' + row['RunID'] + '.sra'])
	subprocess.run(['rm', fm_obj.localReadsDir + row['ProjectID'] + '/' + row['RunID'] + '_2.fastq'])
	subprocess.run(['rm', fm_obj.localReadsDir + row['ProjectID'] + '/' + row['RunID'] + '_1.fastq'])

	row.ReadGroup = rg
	row.Files = fqs
	sample_dt = sample_dt.append(row)
	sample_dt.to_csv(master_sample_data)
	fm_obj.uploadData(master_sample_data)
	break