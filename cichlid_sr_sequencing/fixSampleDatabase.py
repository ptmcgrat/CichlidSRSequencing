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

pdb.set_trace()
