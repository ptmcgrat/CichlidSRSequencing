import os
import pathlib
from pyfaidx import Fasta
from helper_modules.nikesh_file_manager import FileManager as FM

def compute_genome_stats(fasta_path, output_csv_path):
    pyfaidx_obj = Fasta(str(fasta_path))

    with open(output_csv_path, 'w') as fh:
        fh.write('Contig_Name,Contig_Length,Num_Ns,%N_Content,%Repetitive_DNA\n')
        for contig in pyfaidx_obj.keys():
            print(f'Writing stats for {contig} in {fasta_path.name}')
            seq = str(pyfaidx_obj[contig][:])
            scaffold_length = len(seq)
            num_Ns = seq.count('N')
            percent_N = (num_Ns / scaffold_length) * 100 if scaffold_length else 0
            percent_repeat = (sum(base.islower() for base in seq) / scaffold_length) * 100 if scaffold_length else 0
            fh.write(f'{contig},{scaffold_length},{num_Ns},{percent_N:.2f},{percent_repeat:.2f}\n')

def find_fasta_files(input_dir):
    fasta_exts = {'.fasta', '.fa', '.fna'}
    return [f for f in pathlib.Path(input_dir).rglob('*') if f.suffix.lower() in fasta_exts]

def main(input_root_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    fasta_files = find_fasta_files(input_root_dir)

    if not fasta_files:
        print("No FASTA files found.")
        return

    for fasta_path in fasta_files:
        sample_name = fasta_path.stem
        output_csv = pathlib.Path(output_dir) / f"{sample_name}_genome_stats.csv"
        compute_genome_stats(fasta_path, output_csv)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Batch genome stats calculator for FASTA files.")
    parser.add_argument('input_directory', help="Path to the top-level input directory containing FASTA files.")
    parser.add_argument('output_directory', help="Path to the output directory for CSV files.")
    args = parser.parse_args()

    main(args.input_directory, args.output_directory)
