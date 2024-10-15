import pdb, multiprocessing, subprocess, argparse, os
from helper_modules.nikesh_file_manager import FileManager as FM

"""
# for testing the commands needed to get a read depth value:

minimap2 -a -x map-hifi /Output/Genomes/Mzebra_GT3a/Mzebra_GT3a.fasta test.fastq.gz -t 16 -K 150G > test.sam
samtools sort test.sam -o sorted_test.bam -m 10G -@16
samtools coverage sorted_test.bam | awk 'NR > 1 { total_depth+=$7; total_length+=$3-$2+1 } END { print total_depth/NR "X" }'

alignment takes in 2 parameters and an output:
- ref
- fastq
- sam_out


sam_to_bam takes in 1 parameter and an output 
- sam
- sorted_output_bam

calculate_depth takes in 1 parameter and an output 
- sorted_output_bam
- depth.txt

"""

parser = argparse.ArgumentParser()
parser.add_argument('--local_test', action='store_true')
args = parser.parse_args()

fm = FM('Mzebra_GT3')
read_ref_dict = {'/Output/Reads/M_Met_zebra_Male.hifi_reads.fastq.gz': '/Output/Genomes/Mzebra_GT3a/Mzebra_GT3a.fasta',
               '/Output/Reads/N_Met_zebra_Female.hifi_reads.fastq.gz':'/Output/Genomes/MZ_GT3.2/N_Met_zebra_Female_error_corrected_contigs_hs_with_MZ7.1m_molecules_all_scaffolds.fasta',
               '/Output/Reads/combined_MZ4f_reads.fastq.gz': '/Output/Genomes/MZ_GT3.3/MZ4f_ptm_error_corrected_contigs_hs_with_MZ-010-f_molecules_all_scaffolds.fasta',
               '/Output/Reads/G_Aulon_yelhead_Male.hifi_reads.fastq.gz': '/Output/Genomes/YH_GT1.1/G_Aulon_yelhead_Male_error_corrected_contigs_hs_with_kocher_1m_molecules_all_scaffolds.fasta',
               '/Output/Reads/H_Aulon_yelhead_Female.hifi_reads.fastq.gz': '/Output/Genomes/YH_GT1.2/kocher_YH_female_mabs_error_corrected_hs_with_YH15.fasta',
               '/Output/Reads/m84053_231214_031015_s1.hifi_reads.bc2029.fastq.gz': '/Output/Genomes/YH_GT1.3/YH7f_ptm_error_corrected_contigs_hs_with_kocher_2f_molecules_all_scaffolds.fasta'}
genome_names = ['Mzebra_GT3', 'MZ_GT3.2', 'MZ_GT3.3', 'YH_GT1.1', 'YH_GT1.2', 'YH_GT1.3']
# mkdir Mzebra_GT3 MZ_GT3.2 MZ_GT3.3 YH_GT1.1 YH_GT1.2 YH_GT1.3

out_dir = '/Output/Depth'
if args.local_test:
    read_ref_dict = {'/Output/Reads/test.fastq.gz':'/Output/Genomes/Mzebra_GT3a/Mzebra_GT3a.fasta'}
    genome_names = ['test']

def align_reads(data):
    reference_key, genome = data  # Unpack the tuple
    output_sam = f"{out_dir}/{genome}/{genome}.sam"  # Path to output SAM file

    # Run minimap2 and write to SAM file
    with open(output_sam, 'w') as output:
        subprocess.run([
            'minimap2', '-a', '-x', 'map-hifi',
            read_ref_dict[reference_key], reference_key,
            '-t', '16'
        ], stdout=output, check=True)

# Function to convert SAM to BAM (dummy function for now)
def sam_to_bam(data):
    reference_key, genome = data  # Unpack the tuple
    input_sam = f"{out_dir}/{genome}/{genome}.sam"
    output_bam = f"{out_dir}/{genome}/sorted_{genome}.bam"

    # Run samtools sort to convert SAM to BAM
    subprocess.run([
        'samtools', 'sort', input_sam,
        '-o', output_bam, '-m', '10G', '-@16'
    ], check=True)

# Function to calculate depth using samtools coverage
def calculate_depth(data):
    reference_key, genome = data  # Unpack the tuple
    bam_file = f"{out_dir}/{genome}/sorted_{genome}.bam"

    # Run samtools coverage and calculate mean depth
    depth_output = subprocess.run(
        ['samtools', 'coverage', bam_file],
        stdout=subprocess.PIPE, text=True
    )

    # Use awk-like logic to calculate total depth and length from the coverage output
    total_depth = 0
    total_length = 0
    num_records = 0
    for line in depth_output.stdout.splitlines():
        if line.startswith('#'):
            continue  # Skip header lines
        fields = line.split()
        startpos = int(fields[1])
        endpos = int(fields[2])
        depth = float(fields[6])
        num_records += 1
        total_depth += depth
        total_length += endpos - startpos + 1

    # Calculate and print the mean depth (coverage) in "X"
    mean_depth = total_depth / num_records
    output_file = f"{out_dir}/{genome}/depth_summary.txt"

    # Ensure the directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Write the mean depth to the file
    with open(output_file, 'w') as depth_file:
        depth_file.write(f"{genome} mean depth: {mean_depth:.2f}X\n")

# Modified multiprocess function to handle tuples of (reference_key, genome)
def multiprocess(function):
    inputs = [(key, genome_names[idx]) for idx, key in enumerate(read_ref_dict.keys())]
    concurrent_processes = 6
    try:
        with multiprocessing.Pool(processes=concurrent_processes) as pool:
            pool.map(function, inputs)
    except Exception as e:
        print(f"Error occurred during multiprocessing: {e}")

if __name__ == "__main__":
    # Run the multiprocessing for each step
    multiprocess(align_reads)
    multiprocess(sam_to_bam)
    multiprocess(calculate_depth)
