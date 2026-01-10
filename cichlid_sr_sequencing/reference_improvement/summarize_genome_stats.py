import pandas as pd
import pathlib

def summarize_genome_stats(csv_path):
    df = pd.read_csv(csv_path)

    total_contigs = df.shape[0]
    total_length = df['Contig_Length'].sum()
    total_Ns = df['Num_Ns'].sum()

    percent_N = (total_Ns / total_length) * 100 if total_length else 0
    percent_repeat = (df['%Repetitive_DNA'] * df['Contig_Length']).sum() / total_length if total_length else 0

    return {
        'Genome': csv_path.stem.replace('_genome_stats', ''),
        'Num_Contigs': total_contigs,
        'Total_Length': total_length,
        'Total_Ns': total_Ns,
        '%N_Content': round(percent_N, 2),
        '%Repetitive_DNA': round(percent_repeat, 2),
    }

def load_contig_stats(csv_path):
    df = pd.read_csv(csv_path)
    genome_name = csv_path.stem.replace('_genome_stats', '')
    df = df.set_index('Contig_Name')
    df.columns = [f"{genome_name}_{col}" for col in df.columns]
    return df

def main(input_dir, output_excel_path):
    input_dir = pathlib.Path(input_dir)
    csv_files = list(input_dir.glob('*.csv'))

    if not csv_files:
        print("No CSV files found.")
        return

    # Sheet 1: Per-genome summary
    summary_rows = []
    for csv_file in csv_files:
        print(f"Summarizing {csv_file.name}")
        summary = summarize_genome_stats(csv_file)
        summary_rows.append(summary)
    genome_summary_df = pd.DataFrame(summary_rows)

    # Sheet 2: Per-contig wide-format summary
    contig_dfs = []
    for csv_file in csv_files:
        print(f"Processing contigs in {csv_file.name}")
        contig_dfs.append(load_contig_stats(csv_file))

    contig_summary_df = pd.concat(contig_dfs, axis=1)
    contig_summary_df.index.name = 'Contig_Name'

    # Write to Excel
    with pd.ExcelWriter(output_excel_path) as writer:
        genome_summary_df.to_excel(writer, sheet_name='Genome_Summary', index=False)
        contig_summary_df.to_excel(writer, sheet_name='Contig_Summary')

    print(f"\nExcel workbook written to: {output_excel_path}")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Summarize genome stats and save to Excel.")
    parser.add_argument('csv_directory', help="Directory containing per-genome CSVs.")
    parser.add_argument('output_excel', help="Path to output Excel file.")
    args = parser.parse_args()

    main(args.csv_directory, args.output_excel)
