import subprocess as sp, argparse, pandas as pd, shlex, pysam
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

test_contig = 'NW_020193136.1'
for contig in fasta_obj.references:
    sp.run(shlex.split(f'gatk SplitIntervals -O /home/ad.gatech.edu/bio-mcgrath-dropbox/interval_testing -R /home/ad.gatech.edu/bio-mcgrath-dropbox/interval_testing/genome/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -L {contig} -scatter 4'))
# sp.run(shlex.split('/Users/kmnike/bin/gatk-4.2.6.1/gatk ScatterIntervalsByNs -O /Users/kmnike/McGrath/genomics/intervals_testing/intervals_no_Ns.interval_list -R /Users/kmnike/McGrath/genomics/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna --OUTPUT_TYPE ACGT'))