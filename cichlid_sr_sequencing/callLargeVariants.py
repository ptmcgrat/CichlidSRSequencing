import argparse, pdb
from helper_modules.McGrathChimericVariantDetector import ChimericCaller
from helper_modules.file_manager import FileManager as FM
import pandas as pd

parser = argparse.ArgumentParser(usage = 'This script will download fastq data the McGrath lab dropbox and align it to the Genome version of choice')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
parser.add_argument('-p', '--ProjectIDs', nargs = '+', help = 'Restrict analysis to a specific ProjectID')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Make sure genome version is a valid option
if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))

fm_obj.downloadData(fm_obj.localAlignmentFile)
fm_obj.downloadData(fm_obj.localGenomeDir)
a_dt = pd.read_csv(fm_obj.localAlignmentFile)
a_dt = a_dt[(a_dt.GenomeVersion == args.Genome) & (a_dt.ProjectID.isin(args.ProjectIDs))]

sampleIDs = set(a_dt.SampleID)

chimeric_bamfiles = []
discordant_bamfiles = []
geno_bamfiles = []


for sampleID in sampleIDs:
	print(sampleID)
	fm_obj.createSampleFiles(sampleID)
	fm_obj.downloadData(fm_obj.localChimericBamFile)
	fm_obj.downloadData(fm_obj.localChimericBamFile.replace('.bam','.bai'))
	fm_obj.downloadData(fm_obj.localDiscordantBamFile)
	fm_obj.downloadData(fm_obj.localDiscordantBamFile.replace('.bam','.bai'))
	chimeric_bamfiles.append(fm_obj.localChimericBamFile)
	discordant_bamfiles.append(fm_obj.localDiscordantBamFile)

	geno_bamfiles.append(fm_obj.localBamFile)

cc_obj = ChimericCaller(fm_obj.localGenomeFile, 'TestLargeVariantFile.vcf')
#cc_obj.identifyChimericVariants(chimeric_bamfiles)
cc_obj.identifyLargeInsertions(discordant_bamfiles)

