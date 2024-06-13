"""
2024 June 03
NOTE:
These intervals are generated to breakdown the main LGs into as many equal pieces as possible to maximize parallelization when running.
The goal is to run all 96 cores at once and give each core as equal a portion of the 933Mbp of sequence as possible that is contained in the 22 LGs.
that comes out to processing roughly 9.17875Mbp per core. 
To keep things as simple as possible, I will split each LG into intervals that are as close to this 9Mbp value as possible.
This will breakdown large LGs into more intervals and allow for faster processing through GenomicsDBImport
This will also generate 96 Databases of roughly equal size that can each be passed into GenotypeGVCFs whcih will drastically speed up how quickly the VCF generation goes compared to before


The header for each interval_list file comes from the Mzebra_GT3.dict file. Note that the @HD header line is chnaged as follows:
@HD	VN:1.6 
@HD	VN:1.6	SO:coordinate			
Note the extra tabs in the updated header


Since GenotypeGVCFs takes in one DB per process, I'll need to generate 96 Databases
Therefore, each interval_list will have just one interval in it
the --max-intervals-to-import flag will thus become irrelevant since each interval will itself determine parallelization and instead of processing one LG at a time, we'll split into 96 "equal" regions of the genome

TODO:
per LG:
- take in the length of the LG using pyfaidx
- divide by 9.17875 and take the floor value = num_intervals
- divide length by num_intervals and take the floor = equal_bases_for_current_LG
- write num_intervals # of files with the equal_bases_for_current_LG
- ensure that the last inerval goes to the end value equal to the length of the LG
- repeaaattttt 
"""
import os
from pyfaidx import Fasta
import math

# Load the FASTA file using pyfaidx
fasta_file = '/Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_GT3/Mzebra_GT3a.fasta'
fasta = Fasta(fasta_file)

# Filter contigs based on the specified criteria
contigs = {contig: len(fasta[contig]) for contig in fasta.keys() if contig.startswith('LG') and contig != 'Mzebra_UMD2a_mitochrondrial_genome'} # update if contig names are changed 

# Calculate the total number of base pairs
total_base_pairs = sum(contigs.values())

# Determine the approximate interval size
approx_interval_size = total_base_pairs // 96

# Create the output directory if it doesn't exist
output_dir = '/Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/GT3_intervals' #update when rerunning 
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def write_header(file):
    if not os.path.exists(file) or os.path.getsize(file) == 0:
        with open(file, 'w') as f1:
            with open('/Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/GT3_intervals/generate_intervals/header_GT3.tsv', 'r') as f2: # update the file path if needed. Update the header if a new dict file is created
                for line in f2:
                    f1.write(line)
                f1.write('\n')

# Generate the interval list files
total_files_created = 0
file_num = 1

for contig, length in contigs.items():
    # Use rounding to determine the number of intervals for each contig
    num_intervals = max(1, round(length / approx_interval_size))
    # print(f"{contig}: {num_intervals} intervals")  # Print the number of intervals for each contig
    start = 1
    interval_size = math.ceil(length / num_intervals)
    
    for i in range(num_intervals):
        end = min(start + interval_size - 1, length)
        file_path = f'{output_dir}/{file_num}.interval_list'
        write_header(file_path)
        with open(file_path, 'a') as f:
            f.write(f"{contig}\t{start}\t{end}\t+\t.\t\n")
        start = end + 1
        file_num += 1

    total_files_created += num_intervals

print(f"Total files created: {total_files_created}")
