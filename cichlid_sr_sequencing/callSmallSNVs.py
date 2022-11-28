# test line
import argparse, pdb, subprocess, pysam
from helper_modules.file_manager import FileManager as FM
from helper_modules.Timer import Timer

import pandas as pd

parser = argparse.ArgumentParser(usage = 'This script will download fastq data the McGrath lab dropbox and align it to the Genome version of choice')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Make sure genome version is a valid option
if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))

fm_obj.downloadData(fm_obj.localAlignmentFile)

a_dt = pd.read_csv(fm_obj.localAlignmentFile)

a_dt = a_dt[a_dt.GenomeVersion == args.Genome]

sampleIDs = set(a_dt.SampleID)

gvcffiles = []

# Create timer object
timer = Timer()

#fm_obj.downloadData(fm_obj.localGenomeDir)
#fasta_obj = pysam.FastaFile(fm_obj.localGenomeFile)

for sampleID in sampleIDs:
	if sampleID in []:
		continue
	print(sampleID)
	fm_obj.createSampleFiles(sampleID)
	fm_obj.downloadData(fm_obj.localGVCFFile)
	fm_obj.downloadData(fm_obj.localGVCFFile + '.tbi')
	gvcffiles.append(fm_obj.localGVCFFile)
	pdb.set_trace()

processes = []
for contig in fasta_obj.references:
	#timer.start('Calling SNVs for ' + str(len(bamfiles)) + ' bamfiles.')

	# p1 = subprocess.Popen(['bcftools', 'mpileup', '-r', contig, '-C', '50', '-pm2', '-F', '0.2', '-f', fm_obj.localGenomeFile] + bamfiles, stdout = subprocess.PIPE)
	# processes.append(subprocess.Popen(['bcftools', 'call', '-vmO', 'v', '-f', 'GQ', '-o', contig + '.vcf'], stdin = p1.stdout))

	if len(processes) > 23:
		for p in processes:	
			p.communicate()
		processes = []
	#timer.stop()


