import subprocess as sp, argparse, pandas as pd, shlex, pysam
from helper_modules.file_manager import FileManager as FM
from os import listdir
from os.path import isfile, join

parser = argparse.ArgumentParser(usage = 'Testing gatk GenomicsDBImport parallelization optimization')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Make sure genome version is a valid option
if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))

a_dt = pd.read_csv(fm_obj.localAlignmentFile)
a_dt = a_dt[a_dt.GenomeVersion == args.Genome]
sampleIDs = set(a_dt.SampleID)

fasta_obj = pysam.FastaFile(fm_obj.localGenomeFile)

test_contigs = ['NW_020192349.1', 'NW_020192340.1', 'NW_020192348.1']
for contig in test_contigs:
    path = f"/home/ad.gatech.edu/bio-mcgrath-dropbox/interval_testing/{contig + '_intervals'}/"
    # sp.run(shlex.split(f"gatk SplitIntervals -R /home/ad.gatech.edu/bio-mcgrath-dropbox/interval_testing/genome/GCF_000238955.4_M_zebra_UMD2a_genomic.fna --scatter-count 4 -O {path} --subdivision-mode INTERVAL_SUBDIVISION -L {contig}"))
    intervals = [[(filename, sp.check_output(['tail', '-1', path + filename]).decode('utf-8'))] for filename in [f for f in listdir(path) if isfile(join(path, f)) and f.endswith('.interval_list')]]
    with open(f"{path + 'test.intervals'}", 'w') as f:
        f.write(intervals)



"""
path = "{your_path}"

# last_lines = [subprocess.check_output(['tail', '-1', path + filename]) for filename in [f for f in listdir(path) if isfile(join(path, f)) and f.endswith('.txt')]]

filename_last_lines = [[(filename, subprocess.check_output(['tail', '-1', path + filename]))] for filename in [f for f in listdir(path) if isfile(join(path, f)) and f.endswith('.txt')]]

# print(last_lines)
print(filename_last_lines)
"""