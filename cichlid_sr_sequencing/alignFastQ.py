import argparse, os, pysam, pdb, subprocess, sys, datetime
from helper_modules.file_manager import FileManager as FM
from collections import defaultdict
from multiprocessing import cpu_count
import pandas as pd


# Need to make SampleIDs and ProjectIDs mutually exclusive
parser = argparse.ArgumentParser(usage = 'This script will download fastq data the McGrath lab dropbox and align it to the Genome version of choice')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
parser.add_argument('-s', '--SampleIDs', nargs = '+', help = 'Restrict analysis to the following sampleIDs')
parser.add_argument('-p', '--ProjectID', type = str, help = 'Restrict analysis to a specific ProjectID')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Make sure genome version is a valid option
if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))

# Download master sample database and read it in (this contains info on valid sampleIDs, projectIDs, and file locations)
fm_obj.downloadData(fm_obj.localSampleFile)
s_dt = pd.read_csv(fm_obj.localSampleFile)

# If running on projectID, make sure it is valid and subset sample database to those with the right projectID
if args.ProjectID is not None:
	if args.ProjectID not in s_dt.ProjectID:
		raise argparse.ArgumentTypeError('ProjectID does not exist. Options are: ' + ','.join(set(args.ProjectID)))
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

# Download master alignment database to keep track of samples that have been aligned
fm_obj.downloadData(fm_obj.localAlignmentFile)
a_dt = pd.read_csv(fm_obj.localAlignmentFile)

# Make directories necessary for analysis
os.makedirs(fm_obj.localMasterDir, exist_ok = True)
os.makedirs(fm_obj.localTempDir, exist_ok = True)
os.makedirs(fm_obj.localBamRefDir, exist_ok = True)

# Download genome data necessary for analysis		
print('Downloading genome: ' + str(datetime.datetime.now()))
fm_obj.downloadData(fm_obj.localGenomeDir)

# Loop through each sample, determine if it needs to be rerun, and align it to genome
for sample in good_samples:
	# Manually exclude samples that are problematic until debugging can be completed
	if sample in ['SAMEA1904330', 'SAMEA1904323', 'SAMEA4032094', 'SAMEA1904322', 'SAMEA4032090', 'SAMEA1904329', 'SAMEA1904328', 'SAMEA4032091', 'SAMEA1920092']:
		print('Skipping ' + sample + '. Problematic.')
		continue

	# Determine if sample has already been aligned to genome version
	sub_a_dt = a_dt[(a_dt.SampleID == sample) & (a_dt.GenomeVersion == args.Genome)]
	if len(sub_a_dt) != 0:
		print(sample + ' already aligned to ' + args.Genome + '. Skipping...')
		continue

	# Make directories and appropriate files
	print('Processing sample: ' + sample + ': ' + str(datetime.datetime.now()))
	fm_obj.createBamFiles(sample)
	os.makedirs(fm_obj.localSampleBamDir, exist_ok = True)
	sorted_bam = fm_obj.localTempDir + sample + '.sorted.bam'
	sample_dt = s_dt[s_dt.SampleID == sample] # <- dataframe of all runs that match the same sample

	# Loop through all of the runs for a sample
	for i, (index,row) in enumerate(sample_dt.iterrows()):
		print('Downloading uBam files for Run: ' + row['RunID'] + ' :' + str(datetime.datetime.now()))

		# Download unmapped bam file
		uBam_file = fm_obj.localReadsDir + row.File
		fm_obj.downloadData(uBam_file)

		# Create temporary outputfile
		t_bam = fm_obj.localTempDir + sample + '.' + str(i) + '.sorted.bam'

		# Align unmapped bam file following best practices
		# https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently
		print('Aligning fastq files for Run: ' + row['RunID'] + ': ' + str(datetime.datetime.now()))
		# Align fastq files and sort them
		# First command coverts unmapped bam to fastq file, clipping out illumina adapter sequence by setting quality score to #
		command1 = ['gatk', 'SamToFastq', '-I', uBam_file, '--FASTQ', '/dev/stdout', '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2']
		command1 += ['--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', fm_obj.localTempDir]

		# Debugging - useful for ensuring command is working properly, saving intermediate files instead of piping into each other
		#command1 = ['gatk', 'SamToFastq', '-I', uBam_file, '--FASTQ', fm_obj.localTempDir + 'testing.fq', '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2']
		#command1 += ['--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', fm_obj.localTempDir]
		#subprocess.run(command1)
		#pdb.set_trace()

		# Second command aligns fastq data to reference
		command2 = ['bwa', 'mem', '-t', str(cpu_count()), '-M', '-p', fm_obj.localGenomeFile, '/dev/stdin']

		# Debugging - useful for ensuring command is working properly, saving intermediate files instead of piping into each other
		#command2 = ['bwa', 'mem', '-t', str(cpu_count()), '-M', '-p', fm_obj.localGenomeFile, fm_obj.localTempDir + 'testing.fq', '-o', fm_obj.localTempDir + 'testing.sam']
		#subprocess.run(command2)
		#pdb.set_trace()

		# Final command reads read group information to aligned bam file and sorts it
		# Figure out how to keep hard clipping
		command3 = ['gatk', 'MergeBamAlignment', '-R', fm_obj.localGenomeFile, '--UNMAPPED_BAM', uBam_file, '--ALIGNED_BAM', '/dev/stdin']
		command3 += ['-O', t_bam, '--ADD_MATE_CIGAR', 'true', '--CLIP_ADAPTERS', 'false', '--CLIP_OVERLAPPING_READS', 'true']
		command3 += ['--INCLUDE_SCEONDARY_ALIGNMENTS', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1', '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant']
		command3 += ['--ATTRIBUTES_TO_RETAIN', 'XS', '--TMP_DIR', fm_obj.localTempDir]

		# Debugging - useful for ensuring command is working properly, saving intermediate files instead of piping into each other
		command3 = ['gatk', 'MergeBamAlignment', '-R', fm_obj.localGenomeFile, '--UNMAPPED_BAM', uBam_file, '--ALIGNED_BAM', fm_obj.localTempDir + 'testing.sam']
		command3 += ['-O', t_bam, '--ADD_MATE_CIGAR', 'true', '--CLIP_ADAPTERS', 'false', '--CLIP_OVERLAPPING_READS', 'true']
		command3 += ['--INCLUDE_SCEONDARY_ALIGNMENTS', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1', '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant']
		command3 += ['--ATTRIBUTES_TO_RETAIN', 'XS', '--TMP_DIR', fm_obj.localTempDir]
		subprocess.run(command3)
		pdb.set_trace()

		# Figure out how to pipe 3 commands together

		# Remove unmapped reads
		subprocess.run(['rm', '-f', uBam_file])

	print('Merging bam files if necessary... ' + row['RunID'] + ': ' + str(datetime.datetime.now()))
	if i == 0:
		subprocess.run(['mv', t_bam, sorted_bam])
	else:
		inputs = []
		ind_files = [fm_obj.localTempDir + sample + '.' + str(x) + '.sorted.bam' for x in range(i+1)]
		for ind_file in ind_files:
			inputs = inputs + ['-I', ind_file]
		subprocess.run(['gatk', 'MergeSamFiles', '--TMP_DIR', fm_obj.localTempDir] + inputs + ['-O', sorted_bam], stderr = open('TempErrors.txt', 'a'))
		subprocess.run(['rm','-f'] + ind_files)

	print('Marking duplicates... ' + row['RunID'] + ': ' + str(datetime.datetime.now()))
	subprocess.run(['gatk', 'MarkDuplicates', '-I', sorted_bam, '-O', fm_obj.localBamFile, '--tmp-dir', fm_obj.localTempDir, '-OBI', '--spark-runner', 'LOCAL'], stderr = open('TempErrors.txt', 'a'))

	# Remove remaining files
	subprocess.run(['rm','-f',sorted_bam])

	align_file = pysam.AlignmentFile(fm_obj.localBamFile) 
	unmapped = pysam.AlignmentFile(fm_obj.localUnmappedBamFile, mode = 'wb', template = align_file)
	discordant = pysam.AlignmentFile(fm_obj.localDiscordantBamFile, mode = 'wb', template = align_file)
	inversion = pysam.AlignmentFile(fm_obj.localInversionBamFile, mode = 'wb', template = align_file)
	duplication = pysam.AlignmentFile(fm_obj.localDuplicationBamFile, mode = 'wb', template = align_file)
	clipped = pysam.AlignmentFile(fm_obj.localClippedBamFile, mode = 'wb', template = align_file)
	chimeric = pysam.AlignmentFile(fm_obj.localChimericBamFile, mode = 'wb', template = align_file)

	# Go through all reads and process them into appropriate categories
	print('Splitting reads based upon their alignment: ' + str(datetime.datetime.now()))
	read_data = defaultdict(int)
	for read in align_file.fetch(until_eof=True):
		read_data['TotalReads'] += 1
		if read.is_paired:
			if not read.is_unmapped and not read.is_duplicate:
				read_data['MappedReads'] += 1
			if read.is_duplicate:
				read_data['DuplicatedReads'] += 1
			# Both reads are unmapped
			elif read.is_unmapped and read.mate_is_unmapped:
				unmapped.write(read)
				read_data['UnmappedReads'] += 1
			# One read is unmapped
			elif read.is_unmapped or read.mate_is_unmapped:
				discordant.write(read)
				read_data['DiscordantReads'] += 1             
			# Chromosome fusion
			elif read.reference_id!=read.next_reference_id:
				discordant.write(read)
				read_data['DiscordantReads'] += 1
			# Inversion
			elif read.is_reverse == read.mate_is_reverse:
				inversion.write(read)
				read_data['InversionReads'] += 1
			# Duplication
			elif ((read.pos < read.mpos and read.is_reverse) or (read.pos > read.mpos and read.mate_is_reverse)) and abs(read.isize) > 102:
				duplication.write(read)
				read_data['DuplicationReads'] += 1
		else:
			if read.is_unmapped:
				unmapped.write(read)
				read_data['UnmappedReads'] += 1

		# Clipped
		if read.cigarstring is not None:
			for pair in read.cigartuples:
				if pair[0] == 4 and pair[1] > 4:
					clipped.write(read)
					read_data['ClippedReads'] += 1
					break
				elif pair[0] == 5 and pair[1] > 4:
					clipped.write(read)
					read_data['ClippedReads'] += 1
					break
		# Chimeric
		if read.has_tag('SA'):
			chimeric.write(read)
			read_data['ChimericReads'] += 1


	coverage = read_data['MappedReads'] / sum(align_file.lengths) * len(read.seq)

	align_file.close()
	unmapped.close()
	discordant.close()
	inversion.close()
	duplication.close()
	clipped.close()
	chimeric.close()

	pysam.index(fm_obj.localUnmappedBamFile)
	pysam.index(fm_obj.localDiscordantBamFile)
	pysam.index(fm_obj.localInversionBamFile)
	pysam.index(fm_obj.localDuplicationBamFile)
	pysam.index(fm_obj.localClippedBamFile)
	pysam.index(fm_obj.localChimericBamFile)


	sample_data = {'SampleID':sample, 'Organism':sample_dt.Organism.values[0], 'GenomeVersion': args.Genome, 'RunIDs':',,'.join(list(sample_dt.RunID)), 'Coverage':coverage, 'TotalReads': read_data['TotalReads'], 
				   'MappedReads': read_data['MappedReads'], 'UnmappedReads': read_data['UnmappedReads'], 'DiscordantReads': read_data['DiscordantReads'], 'InversionReads': read_data['InversionReads'], 
				   'DuplicationReads': read_data['DuplicationReads'], 'ClippedReads': read_data['ClippedReads'], 'ChimericReads': read_data['ChimericReads'], 'DuplicatedReads': read_data['DuplicatedReads']}

	output = subprocess.run(['conda', 'list'], capture_output = True)
	sample_data['bwa_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('bwa')][0]
	sample_data['gatk_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('gatk4')][0]
	sample_data['pysam_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('pysam')][0]
	sample_data['BamSize'] = os.path.getsize(fm_obj.localBamFile)

	# Upload data and delete
	print('Uploading data: ' + str(datetime.datetime.now()))
	fm_obj.uploadData(fm_obj.localSampleBamDir)
	subprocess.run(['rm','-rf', fm_obj.localSampleBamDir])

	fm_obj.downloadData(fm_obj.localAlignmentFile)
	a_dt = pd.read_csv(fm_obj.localAlignmentFile)
	a_dt = pd.concat([a_dt, pd.DataFrame.from_records([sample_data])])
	a_dt.to_csv(fm_obj.localAlignmentFile, index = False)
	fm_obj.uploadData(fm_obj.localAlignmentFile, upload_async = True)

	print('Finished with sample ' + sample + ': ' + str(datetime.datetime.now()))

