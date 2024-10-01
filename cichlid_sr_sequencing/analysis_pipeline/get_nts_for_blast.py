from pyfaidx import Fasta
from helper_modules.nikesh_file_manager import FileManager as FM
import pdb

fm = FM('Mzebra_GT3')
file = Fasta(fm.localGenomeFile)

with open('lg7_inv.fasta', 'w') as f:
    f.write(file['NC_036786.1'][4752576:5129585].__str__())

