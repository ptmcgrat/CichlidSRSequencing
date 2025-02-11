import pandas as pd
import re, pdb

# Load the metadata Excel file
metadata_path = '/Users/kmnike/Desktop/SRA_metadata_for_short_read_data.xlsx'
fastq_path = '/Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/fastq.txt'
metadata_df = pd.read_excel(metadata_path, sheet_name="SRA_data")

# Read the FASTQ file list
with open(fastq_path, 'r') as file:
    fastq_files = [line.strip() for line in file.readlines()]

# Initialize columns for filename mapping if they don't already exist
filename_cols = ['filename', 'filename2', 'filename3', 'filename4']
for col in filename_cols:
    if col not in metadata_df.columns:
        metadata_df[col] = None

# Mapping fastq files to each sample in metadata
for sample_name in metadata_df['sample_name'].unique():
    # Filter relevant FASTQ files for this sample
    sample_fastqs = [f for f in fastq_files if sample_name in f]
    
    # Sort based on read (R1, R2) and lane (L001, L002)
    pdb.set_trace()
    sample_fastqs_sorted = sorted(sample_fastqs, key=lambda x: (re.search(r'_L00(\d)', x).group(1), re.search(r'_R(\d)', x).group(1)))

    # Assign to metadata columns
    for idx, fastq_file in enumerate(sample_fastqs_sorted[:4]):
        col_name = filename_cols[idx]
        # Fill in the matched row in metadata with the correct file
        metadata_df.loc[(metadata_df['sample_name'] == sample_name) & (metadata_df[col_name].isna()), col_name] = fastq_file

# Save the modified metadata to a new Excel file
output_path = '/Users/kmnike/Desktop/SRA_metadata_with_filenames.xlsx'
metadata_df.to_excel(output_path, sheet_name="SRA_data", index=False)
print(f"Updated metadata saved to {output_path}")
