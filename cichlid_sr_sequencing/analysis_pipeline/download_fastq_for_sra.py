import pdb
import subprocess
import argparse
import multiprocessing
import pandas as pd
from helper_modules.nikesh_file_manager import FileManager as FM

"""
python download_fastq_for_sra.py Mzebra_GT3a SRA_metadata_for_short_read_data.xlsx --local_test
"""

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('genome', help='name to make FileManager work')
parser.add_argument('sra_data', help='path to an excel sheet containing the sampleIDs and associated information for the samples you want to download data for')
parser.add_argument('--local_test', help='for testing using 2024_KumarInversions directory', action='store_true')
args = parser.parse_args()

# FileManager setup
fm = FM(args.genome)

# Load sample data from SRA metadata Excel
s_df = pd.read_excel(args.sra_data, sheet_name='SRA_data')
samples = s_df['sample_name'].to_list()

# Test samples for local test
test_samples = [
    'Figure1_ReferenceImprovement', 'Figure2_Inversions', 'Figure3-Distribution',
    'Figure4_EvolutionaryHistory', 'Figure5-SexDetermination', 'Figure6_Model'
]

# Fetch all fastq paths from Dropbox
if args.local_test:
    all_fastq_paths = subprocess.check_output(
        ['rclone', 'lsf', 'ptm_dropbox:/CoS-EXT/BioSci-Ext/BioSci-McGrath-Ext/2024_KumarInversions/', '-R', '-F', 'p', '--include', '*.png'],
        encoding='utf-8'
    )
    remote_path = 'ptm_dropbox:/CoS-EXT/BioSci-Ext/BioSci-McGrath-Ext/2024_KumarInversions/'
else:
    all_fastq_paths = subprocess.check_output(
        ['rclone', 'lsf', 'ptm_dropbox:/CoS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/Reads/', '-R', '-F', 'p', '--include', '*.fastq.gz'],
        encoding='utf-8'
    )
    remote_path = 'ptm_dropbox:/CoS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/Reads/'

# Prepare fastq paths
all_fastq_paths = all_fastq_paths.strip().split('\n')
all_fastq_paths = [remote_path + x for x in all_fastq_paths]

# Function to check and match sample names (dashed/underscored)
def find_fastq_paths(sample_name, fastq_paths):
    dashed_sample = sample_name.replace('_', '-')
    underscore_sample = sample_name.replace('-', '_')
    
    # Check if the dashed or underscored version exists in the list of paths
    matched_paths = [path for path in fastq_paths if dashed_sample in path or underscore_sample in path]
    return matched_paths

# Main function to download fastqs
def download_fastqs(sample):
    matched_files = find_fastq_paths(sample, all_fastq_paths)
    
    # Define the target directory based on testing or not
    if args.local_test:
        target_dir = '/Users/kmnike/Data/CichlidSequencingData/Reads/SRA'
    else:
        target_dir = '/Data/mcgrath-lab/Data/CichlidSequencingData/Reads/SRA'
    
    if matched_files:
        print(f"Matched files for {sample}:")
        for file in matched_files:
            subprocess.run(['rclone', 'copy', file, target_dir, '-P'])
    else:
        print(f"No files found for {sample}")
        with open('error_samples.txt', 'a') as error_file:
            error_file.write(f"{sample}\n")

# Multiprocessing setup
def multiprocess(samples):
    if args.local_test:
        concurrent_processes = 6
    else:
        concurrent_processes = 48
    
    try:
        with multiprocessing.Pool(processes=concurrent_processes) as pool:
            pool.map(download_fastqs, samples)
    except Exception as e:
        print(f"Error occurred during multiprocessing: {e}")

# Entry point
if __name__ == "__main__":
    if args.local_test:
        multiprocess(test_samples)
    else:
        multiprocess(samples)