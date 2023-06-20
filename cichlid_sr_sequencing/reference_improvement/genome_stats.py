import pdb, pathlib, re
from pyfaidx import Fasta
from helper_modules.nikesh_file_manager import FileManager as FM

"""
stuff to write...
contig name
    pyfaidx_obj.keys()
contig length
Num Ns
%N content
%Repetitive DNA
"""

fm = FM('Mzebra_GT1')
genome = fm.localGenomeFile
pyfaidx_obj = Fasta(genome)

with open('GT1_v2_genome_stats.csv', 'w') as fh:
    fh.write('Contig_Name,Contig_Length,Num_Ns,%N_Content,%Repetitive_DNA\n')
    for contig in pyfaidx_obj.keys():
        print(f'writing stats for {contig}')
        name = contig
        scaffold_length = len(pyfaidx_obj[contig])
        Num_Ns = pyfaidx_obj[contig][:].__str__().count('N')
        percent_N = (Num_Ns / scaffold_length) * 100
        percent_repeat = (sum(base.islower() for base in pyfaidx_obj[contig][:].__str__()) / scaffold_length) * 100
        fh.write(f'{name},{scaffold_length},{Num_Ns},{percent_N:.2f},{percent_repeat:.2f}\n')

"""
/Users/kmnike/miniforge3/envs/pipeline_x86/bin/python3 genome_stats.py

"""