import subprocess as sp, argparse, pandas as pd, shlex, pysam, os
from helper_modules.file_manager import FileManager as FM


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
    # sp.run(shlex.split(f"gatk SplitIntervals -R /home/ad.gatech.edu/bio-mcgrath-dropbox/interval_testing/genome/GCF_000238955.4_M_zebra_UMD2a_genomic.fna --scatter-count 4 -O /home/ad.gatech.edu/bio-mcgrath-dropbox/interval_testing/{contig + '_intervals'}/ --subdivision-mode INTERVAL_SUBDIVISION -L {contig}"))
    with open(f"{path + 'master.intervals_list'}", 'w') as fh:
        for dir in os.listdir():
            for file in dir:
                fh.write(sp.run(shlex.split(f"tail -n1 {file}")))
