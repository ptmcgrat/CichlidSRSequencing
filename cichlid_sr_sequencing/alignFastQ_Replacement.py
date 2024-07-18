import argparse, os, math, pysam
from helper_modules.file_manager_Replacement import FileManager as FM
from helper_modules.alignment_worker_Replacement import AlignmentWorker as AW
#from helper_modules.Timer import Timer

import pandas as pd


import argparse, os, pdb, subprocess, sys, datetime
from collections import defaultdict
from multiprocessing import cpu_count


# Need to make SampleIDs and ProjectIDs mutually exclusive
parser = argparse.ArgumentParser(usage = 'This script will download fastq data the McGrath lab dropbox and align it to the Genome version of choice. It will also create gvcf files')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
parser.add_argument('-n', '--NumberParallel', type = int, default = 48, help = 'Specify the number of samples run in parallel. Default is 96')
parser.add_argument('-s', '--SampleIDs', nargs = '+', type = str, help = 'Restrict analysis to the following sampleIDs')
parser.add_argument('-p', '--ProjectIDs', nargs = '+', type = str, help = 'Restrict analysis to a specific ProjectIDs')
parser.add_argument('-e', '--Ecogroups', nargs = '+', type = str, help = 'Restrict analysis to specific Ecogroups')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Create timer object to keep track of time
#timer = Timer()

# This command identifies all the samples that will need to be run based upon user input and stores in self.samples and self.s_dt
fm_obj.setSamples(projectIDs = args.ProjectIDs, sampleIDs = args.SampleIDs, ecogroupIDs = args.Ecogroups)

# Download genome data necessary for analysis
timer.start('Downloading genome')		
fm_obj.downloadData(fm_obj.localGenomeDir)
timer.stop()

# Create alignment worker object:
aw_obj = AW(args.Genome, fm_obj)

print('The following ' + str(len(fm_obj.samples)) + ' samples will be analyzed:')
print(','.join(fm_obj.samples))

timer.start('  Parallel Downloading uBams files')
aw_obj.downloadReadData('Popen')
timer.stop()

print('  Aligning reads to create sorted Bam files')
aw_obj.alignData()
#aw_obj.alignData(linked=True)

print('  Marking duplicates for bamfiles')
#aw_obj.markDuplicates()
aw_obj.markDuplicates(parallel = True)

timer.start('  Splitting reads based upon their alignment')
aw_obj.splitBamfiles()
timer.stop()

print('  Calling haplotypes to create gvcf files')
aw_obj.createGVCF(parallel = True)

processes = []
for sample in fm_obj.samples:
	fm_obj.createSampleFiles(sample)
	fm_obj.uploadData(fm_obj.localSampleBamDir)

	stats = aw_obj.calculateStats(sample)

	read_length = s_dt[s_dt['SampleID'] == sample]['ReadLength'].values[0]/2
	reference_size = sum(pysam.FastaFile(fm_obj.localGenomeFile).lengths)
	coverage = stats['all'] * read_length / reference_size

	sample_data = {'SampleID':sample, 'Organism':s_dt[s_dt['SampleID'] == sample].Organism.values[0], 'GenomeVersion': args.Genome, 'RunIDs':',,'.join(list(s_dt[s_dt['SampleID'] == sample].RunID)), 'ProjectID':s_dt[s_dt['SampleID'] == sample]['ProjectID'].values[0], 
			   'Coverage':coverage, 'TotalReads':stats['all'], 'UnmappedReads':stats['unmapped'], 'DiscordantReads':stats['discordant'], 'InversionReads':stats['inversion'],
			   'DuplicationReads':stats['duplication'], 'ClippedReads':stats['clipped'], 'ChimericReads':stats['chimeric']}

	output = subprocess.run(['conda', 'list'], capture_output = True)
	sample_data['bwa_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('bwa')][0]
	sample_data['gatk_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('gatk4')][0]
	sample_data['pysam_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('pysam')][0]
	sample_data['BamSize'] = os.path.getsize(fm_obj.localBamFile)

	# Upload data and delete
	#subprocess.run(['rm','-rf', fm_obj.localSampleBamDir])
	#subprocess.run(['rm','-rf', fm_obj.localTempDir])
	fm_obj.downloadData(fm_obj.localAlignmentFile)
	a_dt = pd.read_csv(fm_obj.localAlignmentFile)
	a_dt = pd.concat([a_dt, pd.DataFrame.from_records([sample_data])])
	a_dt.to_csv(fm_obj.localAlignmentFile, index = False)
	fm_obj.uploadData(fm_obj.localAlignmentFile)
	timer.stop()
	print(' Finished with sample ' + sample + ': ' + str(datetime.datetime.now()))
	print()

