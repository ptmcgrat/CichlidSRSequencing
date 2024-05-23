import argparse, os
from helper_modules.file_manager import FileManager as FM
from helper_modules.alignment_worker import AlignmentWorker as AW
from helper_modules.Timer import Timer

import pandas as pd


import argparse, os, pdb, subprocess, sys, datetime
from collections import defaultdict
from multiprocessing import cpu_count


# Need to make SampleIDs and ProjectIDs mutually exclusive
parser = argparse.ArgumentParser(usage = 'This script will download fastq data the McGrath lab dropbox and align it to the Genome version of choice. It will also create gvcf files')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
parser.add_argument('-s', '--SampleIDs', nargs = '+', help = 'Restrict analysis to the following sampleIDs')
parser.add_argument('-p', '--ProjectID', type = str, help = 'Restrict analysis to a specific ProjectID')
parser.add_argument('-e', '--Ecogroups', nargs = '+', type = str, help = 'Restrict analysis to specific Ecogroups')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Create timer object to keep track of time
timer = Timer()

# Make sure genome version is a valid option
if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))

# Download master sample database and read it in (this contains info on valid sampleIDs, projectIDs, and file locations)
fm_obj.downloadData(fm_obj.localSampleFile_v2)
s_dt = pd.read_excel(fm_obj.localSampleFile_v2, sheet_name = 'SampleLevel')

# If running on projectID, make sure it is valid and subset sample database to those with the right projectID
if args.ProjectID is not None:
	if args.ProjectID not in set(s_dt.ProjectID):
		raise argparse.ArgumentTypeError('ProjectID ' + args.ProjectID + ' does not exist. Options are: ' + ','.join(set(s_dt.ProjectID)))
	s_dt = s_dt[s_dt.ProjectID == args.ProjectID]
	good_samples = set(s_dt.SampleID)

# If running on sampleIDs, make sure they are valid
if args.SampleIDs is not None:
	bad_samples = []
	for sample in args.SampleIDs:
		if sample not in list(s_dt.SampleID):
			bad_samples.append(sample)

	if len(bad_samples) > 0:
		raise argparse.ArgumentTypeError('The following samples were not found in sample database: ' + ','.join(bad_samples))
	good_samples = set(args.SampleIDs)

# If running on ecogroup, make sure they are valid
if args.Ecogroups is not None:
	good_samples = set()
	for ecogroup in args.Ecogroups:
		if ecogroup not in set(s_dt.Ecogroup_PTM):
			raise argparse.ArgumentTypeError('Ecogroup ' + ecogroup + ' does not exist. Options are: ' + ','.join(set(s_dt.Ecogroup_PTM)))
		sub_dt = s_dt[s_dt.Ecogroup_PTM == ecogroup]
		good_samples.update(sub_dt.SampleID)

# If no filtering options are given, run on all samples
if args.ProjectID is None and args.SampleIDs is None and args.Ecogroups is None:
	good_samples = set(s_dt.SampleID)

# Download master alignment database to keep track of samples that have been aligned
fm_obj.downloadData(fm_obj.localAlignmentFile)
a_dt = pd.read_csv(fm_obj.localAlignmentFile)


# Make directories necessary for analysis
os.makedirs(fm_obj.localMasterDir, exist_ok = True)
os.makedirs(fm_obj.localTempDir, exist_ok = True)
os.makedirs(fm_obj.localBamRefDir, exist_ok = True)

# Download genome data necessary for analysis
timer.start('Downloading genome')		
fm_obj.downloadData(fm_obj.localGenomeDir)
timer.stop()
pdb.set_trace()

# Loop through each sample, determine if it needs to be rerun, and align it to genome
for sample in good_samples:
	# Manually exclude samples that are problematic until debugging can be completed
	# Also SAMEA1904330 'SAMEA1904323', 'SAMEA4032094', 'SAMEA1904322', 'SAMEA4032090', 'SAMEA1904329', 'SAMEA1904328', 'SAMEA4032091', 'SAMEA1920092'
	if sample in ['SAMEA2661255', 'SAMEA2661406']:
		print('Skipping ' + sample + '. Problematic.')
		continue

	# Determine if sample has already been aligned to genome version
	sub_a_dt = a_dt[(a_dt.SampleID == sample) & (a_dt.GenomeVersion == args.Genome)]
	if len(sub_a_dt) != 0:
		print(sample + ' already aligned to ' + args.Genome + '. Skipping...')
		continue

	# Make directories and appropriate files
	print(' Processing sample: ' + sample + '; Start time: ' + str(datetime.datetime.now()))

	# Loop through all of the runs for a sample
	aw_obj = AW(fm_obj, s_dt, sample)
	timer.start('  Downloading uBam files for Sample: ' + sample)
	fm_obj.downloadData(fm_obj.localSampleBamDir)
	os.makedirs(fm_obj.localTempDir, exist_ok = True)

	#aw_obj.downloadReadData()
	timer.stop()
	timer.start('  Aligning fastq files for Sample: ' + sample)
	#aw_obj.alignData()
	timer.stop()
	timer.start('  Marking duplicates for Sample: ' + sample)
	#aw_obj.markDuplicates()
	timer.stop()
	timer.start('  Splitting reads based upon their alignment for Sample: ' + sample)
	#aw_obj.splitBamfiles()
	timer.stop()
	timer.start('  Creating GVCF file for Sample: ' + sample)
	aw_obj.createGVCF()
	timer.stop()
	timer.start('  Uploading data for Sample: ' + sample)
	fm_obj.uploadData(fm_obj.localSampleBamDir)
	stats = aw_obj.calculateStats()

	read_length = s_dt[s_dt['SampleID'] == sample]['ReadLength'].values[0]/2
	reference_size = sum(pysam.FastaFile(fm_obj.localGenomeFile).lengths)
	coverage = stats['all'] * read_length / reference_size

	sample_data = {'SampleID':sample, 'Organism':s_dt.Organism.values[0], 'GenomeVersion': args.Genome, 'RunIDs':',,'.join(list(s_dt[s_dt['SampleID'] == sample].RunID)), 'ProjectID':s_dt[s_dt['SampleID'] == sample]['ProjectID'].values[0], 
				   'Coverage':coverage, 'TotalReads':stats['all'], 'UnmappedReads':stats['unmapped'], 'DiscordantReads':stats['discordant'], 'InversionReads':stats['inversion'],
				   'DuplicationReads':stats['duplication'], 'ClippedReads':stats['clipped'], 'ChimericReads':stats['chimeric']}

	output = subprocess.run(['conda', 'list'], capture_output = True)
	sample_data['bwa_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('bwa')][0]
	sample_data['gatk_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('gatk4')][0]
	sample_data['pysam_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('pysam')][0]
	sample_data['BamSize'] = os.path.getsize(fm_obj.localBamFile)

	# Upload data and delete
	subprocess.run(['rm','-rf', fm_obj.localSampleBamDir])
	subprocess.run(['rm','-rf', fm_obj.localTempDir])
	fm_obj.downloadData(fm_obj.localAlignmentFile)
	a_dt = pd.read_csv(fm_obj.localAlignmentFile)
	a_dt = pd.concat([a_dt, pd.DataFrame.from_records([sample_data])])
	a_dt.to_csv(fm_obj.localAlignmentFile, index = False)
	fm_obj.uploadData(fm_obj.localAlignmentFile, upload_async = True)
	timer.stop()
	print(' Finished with sample ' + sample + ': ' + str(datetime.datetime.now()))
	print()

