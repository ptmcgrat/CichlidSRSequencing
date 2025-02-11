from Bio import SeqIO
import os, pathlib, pdb

"""
ToDo: 
for each genome the following needs to be done: 
- Each LG needs to be read in and written into a new FASTA file following the conventions written in Table S3. Need to macth the metrics from that with the genomes.
- Each unmapped contig (starting with 'NW') needs to be written as "scaffold01", "scaffold02", etc
- The mito genome needs to be written as "Mzebra_UMD2a_mitochondrial_genome"
- Any Unmapped contig needs to be further writtten as "scaffoldXX" at the end of the file

"""
fasta_dir = '/Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/S3_scaffolded_genomes'
out_dir = '/Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/ncbi_genomes/'
# fasta_dir = '/Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/p_nyererei_v2_genome_for_ncbi'
# out_dir = '/Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/p_nyererei_v2_genome_for_ncbi'
YH_files = ['G_Aulon_yelhead_Male_error_corrected_contigs_hs_with_kocher_1m_molecules_all_scaffolds.fasta', 'YH7f_ptm_error_corrected_contigs_hs_with_kocher_2f_molecules_all_scaffolds.fasta', 'kocher_YH_female_mabs_error_corrected_hs_with_YH15.fasta']
MZ_files = ['MZ4f_ptm_anchorMZ4f_ptm_error_corrected_contigs_hs_with_MZ-010-f_molecules_all_scaffolds.fasta', 'N_Met_zebra_Female_error_corrected_contigs_hs_with_MZ7.1m_molecules_all_scaffolds.fasta']
PN_file = ['p_nyererei_v2_genome_for_ncbi']

def rename_fasta_headers(input_fasta, output_fasta):
    scaffold_count = 1
    contig_count = 1

    with open(output_fasta, 'w') as output_handle:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            if record.__len__() >= 200:
                if record.id.startswith('Super-Scaffold'):
                    new_header = f"scaffold{scaffold_count:03d}"
                    scaffold_count += 1
                elif record.id.startswith('ptg'):
                    new_header = f"contig{contig_count:03d}"
                    contig_count += 1
                else:
                    new_header = record.id  # Leave unchanged if it doesn't match the expected patterns

                # Write the record with the new header
                output_handle.write(f">{new_header}\n{record.seq}\n")

def verify_sequences(input_fasta, output_fasta):
    """Verify that all sequences are identical between the input and output files."""
    input_records = {record.id: str(record.seq) for record in SeqIO.parse(input_fasta, 'fasta') if record.__len__() >= 200}
    output_records = {record.id: str(record.seq) for record in SeqIO.parse(output_fasta, 'fasta') if record.__len__() >= 200}
    # input_records = {record.id: str(record.seq) for record in SeqIO.parse(input_fasta, 'fasta')} # if record.__len__() >= 200}
    # output_records = {record.id: str(record.seq) for record in SeqIO.parse(output_fasta, 'fasta')} # if record.__len__() >= 200}

    if len(input_records) != len(output_records):
        print(f"Mismatch in number of records between {input_fasta} and {output_fasta}\n")
        return False

    for (input_id, input_seq), (output_id, output_seq) in zip(input_records.items(), output_records.items()):
        if input_seq != output_seq:
            print(f"Sequence mismatch for scaffold/contig: {input_id} vs {output_id}")
            return False

    print(f"Verification passed: Sequences are identical between {input_fasta} and {output_fasta}\n")
    return True

if __name__ == "__main__":
    # List your input fasta files here
    input_files = os.listdir(fasta_dir)

    for input_fasta in input_files:
        if input_fasta.endswith('.fasta'):
            # pdb.set_trace()
            # Define output file names based on the input file names
            input_fasta = os.path.join(fasta_dir, input_fasta)
            output_fasta = os.path.join(out_dir, f"renamed_{os.path.basename(input_fasta)}")
            # Rename headers and write to new output file
            rename_fasta_headers(input_fasta, output_fasta)

            # Verify the integrity of sequences
            if verify_sequences(input_fasta, output_fasta):
                print(f"Integrity check passed for {input_fasta}")
            else:
                print(f"Integrity check failed for {input_fasta}")

        print(f"Processed {input_fasta} and saved as {output_fasta}")
