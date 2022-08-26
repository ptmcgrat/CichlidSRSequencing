import argparse, pdb
from helper_modules.McGrathChimericVariantDetector import ChimericCaller
from helper_modules.file_manager import FileManager as FM
import pandas as pd

parser = argparse.ArgumentParser(usage = 'This script will download fastq data the McGrath lab dropbox and align it to the Genome version of choice')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
parser.add_argument('AnalysisID', type = str, help = 'Name of Analysis ID holding ')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Make sure genome version is a valid option
if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))

fm_obj.createAnalysisIDFiles(args.AnalysisID)
fm_obj.downloadData(fm_obj.localAnalysisFile)

a_dt = pd.read_csv(fm_obj.localAnalysisFile)

sampleIDs = set(a_dt[a_dt.Ecogroup != 'Riverine'].SampleID)

bamfiles = []
count = 0

pdb.set_trace()

for sampleID in sampleIDs:
	print(sampleID)
	fm_obj.createBamFiles(sampleID)
	fm_obj.downloadData(fm_obj.localChimericBamFile)
	fm_obj.downloadData(fm_obj.localChimericBamFile + '.bai')
	bamfiles.append(fm_obj.localChimericBamFile)
	if count > 10:
		break
	count += 1


cc_obj = ChimericCaller(fm_obj.localGenomeFile)
cc_obj.identifyChimericLocations(bamfiles)

pdb.set_trace()
