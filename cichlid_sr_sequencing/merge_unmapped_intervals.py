import subprocess as sp, argparse, shlex, pysam, os
from helper_modules.file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage = 'This script will concatenate VCF ouputs from the CallSmallSNVs.py pipeline into a master VCF file')
# parser.add_argument('Genome', type = str, help = 'Version of the genome use')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM('Mzebra_UMD2a')

#### get unmapped contig names
fasta_obj = pysam.FastaFile(fm_obj.localGenomeFile)
contigs = fasta_obj.references[22:]

unmapped_contigs_dir = os.getcwd() + '/' + 'unmapped_contig_intervals/'
out_dir = '/Users/kmnike/McGrath/genomics/CichlidSRSequencing/cichlid_sr_sequencing/intervals_unmapped_contigs/'
file_list = ['0000-scattered.interval_list', '0001-scattered.interval_list', '0002-scattered.interval_list', '0003-scattered.interval_list']


for file in os.listdir(unmapped_contigs_dir):
    current_contig = unmapped_contigs_dir + file + '/'
    file_name = out_dir + file.removesuffix('_intervals') + '.interval_list'
    with open(file_name, 'w') as f1:
        for interval in file_list:
            with open(current_contig + interval, 'r') as f2:
                if os.path.getsize(file_name) == 0:
                    print(f"{current_contig + interval}")
                    f1.write(sp.check_output(shlex.split(f"cat {current_contig + interval}"), encoding='utf-8'))
                else:
                    f1.write(sp.check_output(shlex.split(f"tail -1 {current_contig + interval}"), encoding='utf-8'))
