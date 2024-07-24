import pdb, pathlib, re
from pyfaidx import Fasta # youll need to conda install this into your conda env
from helper_modules.nikesh_file_manager import FileManager as FM # Keep the helper_modules folder in the same directory as this script, but be sure that the "nikesh_file_manager.py" script is within helper_moduels.
# Here's what the structure should look like. Assume you just put these files into a directory called "genome_analysis" on your laptop:
# genome_analysis
    # genome_stats.py
    # csv_output_files_when_you_run_thge_script.csv
    # helper_modules
        # nikesh_file_manager.py
# Filemanager may also have some dependencies u may have to install. Conda installs or pip installs should solve those issues if they come up! 
# Let me know if something else doesn't work! The script should run for each genome in under a minute or 2 

"""
stuff to write...
contig name
    pyfaidx_obj.keys()
contig length
Num Ns
%N content
%Repetitive DNA
"""

# fm = FM('Mzebra_GT2')
# genome = fm.localGenomeFile
# pyfaidx_obj = Fasta(genome)

pyfaidx_obj = Fasta('/Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/reference_improvement/Mzebra_GT3a.fasta') # replace the filepath with the path to each dgenies genome when u run the script 

with open('/Users/kmnike/Desktop/Mzebra_GT3a_stats.csv', 'w') as fh: 
    fh.write('Contig_Name,Contig_Length,Num_Ns,%N_Content,%Repetitive_DNA\n')
    for contig in pyfaidx_obj.keys():
        print(f'writing stats for {contig}')
        name = contig
        scaffold_length = len(pyfaidx_obj[contig])
        Num_Ns = pyfaidx_obj[contig][:].__str__().count('N')
        percent_N = (Num_Ns / scaffold_length) * 100
        percent_repeat = (sum(base.islower() for base in pyfaidx_obj[contig][:].__str__()) / scaffold_length) * 100
        fh.write(f'{name},{scaffold_length},{Num_Ns},{percent_N:.2f},{percent_repeat:.2f}\n')



