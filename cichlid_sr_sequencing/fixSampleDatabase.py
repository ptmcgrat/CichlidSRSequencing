import argparse, pdb
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

sample_species_dt = new_dt.groupby(['BioSample','ArrayExpress-SPECIES']).count()['Organism'].reset_index()[['BioSample','ArrayExpress-SPECIES']]

pd.merge(sample_dt, sample_species_dt)
new_dt = pd.merge(sample_species_dt, sample_dt, left_on = 'BioSample', right_on = 'SampleID')
new_dt.to_csv(fm_obj.localSampleFile, index = False)
fm_obj.uploadData(fm_obj.localSampleFile)

pdb.set_trace()
