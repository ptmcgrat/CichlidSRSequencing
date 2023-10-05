import pdb
import subprocess as sp
from pyfaidx import Fasta
# conda install -c bioconda pyfaidx

# replace the path below to location of the Mzebra_UMD2a reference genome
fasta = Fasta('/Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna')

# for contig in fasta.keys():
#     # you'll need an if statement to exclude the first 22 chromosomes. Hint: they all start with the same 2 letters. 
#     # use the linux touch command to create a file with the name of teh contig you're iterating through and store the file in a trash directory
#     print('hi')
#     pdb.set_trace()
#     print('hi again')

p = sp.Popen(['ls', '-l'])
pdb.set_trace()
p.wait()

