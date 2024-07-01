#!/Users/kmnike/anaconda3/envs/variant/bin/python3
import pdb
from pyfaidx import Fasta
from helper_modules.nikesh_file_manager import FileManager as FM

contigs = ['NC_036780.1', 'NC_036781.1', 'NC_036782.1', 'NC_036783.1', 'NC_036784.1', 'NC_036785.1', 'NC_036786.1', 'NC_036787.1', 'NC_036788.1', 'NC_036789.1', 'NC_036790.1', 'NC_036791.1', 'NC_036792.1', 'NC_036793.1', 'NC_036794.1', 'NC_036795.1', 'NC_036796.1', 'NC_036797.1', 'NC_036798.1', 'NC_036799.1', 'NC_036800.1', 'NC_036801.1']
for contig in contigs:
    with open ('4_per_LG.interval_list', 'r') as f1:
        with open(contig + '.interval_list', 'w') as f2:
            for line in f1:
                if line.startswith('@') or line.startswith(contig):
                    f2.write(line)
            