import os
from Bio import SeqIO

def get_fasta_info(fasta_file):
    """Returns a dictionary with header as key and sequence length as value."""
    fasta_info = {}
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fasta_info[record.id] = len(record.seq)
    return fasta_info

def process_directory(directory):
    """Processes all FASTA files in a directory and returns a dictionary of filename to header/length pairs."""
    fasta_data = {}
    for filename in os.listdir(directory):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            file_path = os.path.join(directory, filename)
            fasta_data[filename] = get_fasta_info(file_path)
    return fasta_data

def write_combined_fasta_info(original_data, ncbi_data, output_dir, filename_mapping):
    """Writes combined information for each corresponding FASTA file."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for original_filename, original_info in original_data.items():
        # Extract corresponding NCBI filename based on the mapping
        ncbi_filename = filename_mapping.get(original_filename)
        if ncbi_filename and ncbi_filename in ncbi_data:
            output_file = os.path.join(output_dir, f"{original_filename}_vs_{ncbi_filename}.txt")
            with open(output_file, "w") as out:
                for original_header, length in original_info.items():
                    # Find the corresponding NCBI header that has the same length
                    ncbi_header = next((header for header, ncbi_length in ncbi_data[ncbi_filename].items() if ncbi_length == length), None)
                    if ncbi_header:
                        out.write(f"{original_header}\t{ncbi_header}\t{length}\n")

# Directories for the original and NCBI genomes
original_dir = "/Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/S3_scaffolded_genomes"
ncbi_dir = "/Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/ncbi_genomes"
output_dir = "/Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/scaffolds"

# Placeholder for filename mapping
filename_mapping = {
    'N_Met_zebra_Female_error_corrected_contigs_hs_with_MZ7.1m_molecules_all_scaffolds.fasta': 'MZ_GT3.2.fasta',
    'MZ4f_ptm_error_corrected_contigs_hs_with_MZ-010-f_molecules_all_scaffolds.fasta': 'MZ_GT3.3.fasta',
    'G_Aulon_yelhead_Male_error_corrected_contigs_hs_with_kocher_1m_molecules_all_scaffolds.fasta': 'YH_GT1.1.fasta',
    'kocher_YH_female_mabs_error_corrected_hs_with_YH15.fasta': 'YH_GT1.2.fasta',
    'YH7f_ptm_error_corrected_contigs_hs_with_kocher_2f_molecules_all_scaffolds.fasta': 'YH_GT1.3.fasta',
    'PunNye_v2_hybrid_scaffold_genome.fasta': 'renamed_PunNye_v2_hybrid_scaffold_genome.fasta'
    }

# Process both directories
original_fasta_data = process_directory(original_dir)
ncbi_fasta_data = process_directory(ncbi_dir)

# Write the combined output files
write_combined_fasta_info(original_fasta_data, ncbi_fasta_data, output_dir, filename_mapping)

print("Combined FASTA info files created successfully.")
