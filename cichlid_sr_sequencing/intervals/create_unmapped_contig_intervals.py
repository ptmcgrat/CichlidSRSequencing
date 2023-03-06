import subprocess as sp, argparse, pandas as pd, shlex, pysam, os
from helper_modules.file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage = 'This script will generate a master intervals list for all of the unmapped contigs in Mzebra_UMD2a')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

fasta_obj = pysam.FastaFile(fm_obj.localGenomeFile)
contigs = fasta_obj.references[22:]

#### temp filepath to the location in which the intervals will be generated. 
try:
    os.makedirs(os.getcwd() + '/Temp/')
except FileExistsError:
    pass
temp_dir = os.getcwd() + '/Temp/'


#### use SplitIntervals to generate 4 intervals per unmapped contig. the tool spits out 4 files with the same name by default, so the -O flag is used to pass each set of files into a new directory with the name of the contig used. the INTERNAL SUBDIVISION flag splits the intervals as evenly as possible per contig. 
for contig in contigs:
    sp.run(shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk SplitIntervals -R /Users/kmnike/McGrath/genomics/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna --scatter-count 4 -O {temp_dir + contig + '_intervals'}/ --subdivision-mode INTERVAL_SUBDIVISION -L {contig}"))

#### write intervals to a master interval list file. Navigate to the dir containing all the contigs' directories which have the intervals.
#### check if the master.intervals_list file is empty. If so, write the whole first file (to get the header). Otherwise, use tail -n1 to get the last line which has the interval for that contig. The files are not read in the correct order, so we need to sort the intervals afterwards
with open('master.intervals_list', 'w+') as fh:
    for contig in os.listdir(temp_dir):
            for interval in os.listdir(f"{temp_dir +'/' + contig}"):
                if os.path.getsize(f"{os.getcwd() + '/master.intervals_list'}") == 0:
                    fh.write(sp.check_output(shlex.split(f"cat {temp_dir + '/' + contig + '/' + interval}"), encoding = 'utf-8'))
                else:
                    fh.write(sp.check_output(shlex.split(f"tail -n1 {temp_dir + '/' + contig + '/' + interval}"), encoding='utf-8'))

# sorting the interval file using GATK tools
sp.run(shlex.split('/Users/kmnike/bin/gatk-4.2.6.1/gatk IntervalListTools -I master.intervals_list --SORT True -O master.intervals_list'))