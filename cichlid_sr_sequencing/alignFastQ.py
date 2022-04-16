import argparse, os, pysam, pdb, subprocess, sys, datetime
from helper_modules.file_manager import FileManager as FM
from collections import defaultdict
from multiprocessing import cpu_count
import pandas as pd


parser = argparse.ArgumentParser(usage = 'This script will download fastq data the McGrath lab dropbox and align it to the Genome version of choice')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
parser.add_argument('-s', '--SampleIDs', nargs = '+', help = 'Restrict analysis to the following sampleIDs')
parser.add_argument('-p', '--ProjectID', type = str, help = 'Restrict analysis to a specific ProjectID')
args = parser.parse_args()

fm_obj = FM(args.Genome)

if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))


fm_obj.downloadData(fm_obj.localSampleFile)
s_dt = pd.read_csv(fm_obj.localSampleFile)

if args.ProjectID is not None:
	s_dt = s_dt[s_dt.ProjectID == args.ProjectID]

fm_obj.downloadData(fm_obj.localAlignmentFile)
a_dt = pd.read_csv(fm_obj.localAlignmentFile)

if args.SampleIDs is not None:
	bad_samples = []
	for sample in args.SampleIDs:
		if sample not in list(s_dt.SampleID):
			bad_samples.append(sample)

	if len(bad_samples) > 0:
		raise argparse.ArgumentTypeError('The following samples were not found: ' + ','.join(bad_samples))

	good_samples = set(args.SampleIDs)

else:
	good_samples = set(s_dt.SampleID)

# Make directories necessary for analysis
os.makedirs(fm_obj.localMasterDir, exist_ok = True)
os.makedirs(fm_obj.localTempDir, exist_ok = True)
os.makedirs(fm_obj.localBamRefDir, exist_ok = True)

# Download genome data necessary for analysis		
print('Downloading genome and sample file: ' + str(datetime.datetime.now()))
fm_obj.downloadData(fm_obj.localGenomeDir)

# Get list of bam files already run
existing_bams = list(a_dt.SampleID)


for sample in good_samples:
	if sample in existing_bams:
		print(sample + ' already analyzed. Skipping...')
		continue

	print('Processing sample: ' + sample + ': ' + str(datetime.datetime.now()))
	fm_obj.createBamFiles(sample)
	os.makedirs(fm_obj.localSampleBamDir, exist_ok = True)

	unsorted_sam = fm_obj.localTempDir + sample + '.unsorted.sam'

	sample_dt = s_dt[s_dt.SampleID == sample]

	# Loop through all of the fastq files


	for i, (index,row) in enumerate(sample_dt.iterrows()):

		print('Downloading fastq files for Run: ' + row['RunID'] + ' :' + str(datetime.datetime.now()))
		# Download fastq files
		fq1 = fm_obj.localReadsDir + row.Files.split(',,')[0]
		fq2 = fm_obj.localReadsDir + row.Files.split(',,')[1]
		fm_obj.downloadData(fq1)
		fm_obj.downloadData(fq2)

		print('Aligning fastq files for Run: ' + row['RunID'] + ': ' + str(datetime.datetime.now()))
		# Align fastq files and sort them
		t_sam = fm_obj.localTempDir + sample + '.' + str(i) + '.unsorted.sam'
		subprocess.run(['bwa', 'mem', '-t', str(cpu_count()), '-R', row.ReadGroup.replace('\t','\\t'), '-M', fm_obj.localGenomeFile, fq1, fq2], stdout = open(t_sam, 'w'), stderr = open('TempErrors.txt', 'a'))
		subprocess.run(['rm', '-f', fq1, fq2])

	if i == 0:
		subprocess.run(['mv', t_sam, unsorted_sam])
	else:
		inputs = []
		ind_files = [fm_obj.localTempDir + sample + '.' + str(x) + '.unsorted.sam' for x in range(i+1)]
		for ind_file in ind_files:
			inputs = inputs + ['-I', ind_file]
		subprocess.run(['gatk', 'MergeSamFiles'] + inputs + ['-O', unsorted_sam], stderr = open('TempErrors.txt', 'a'))
		subprocess.run(['rm','-f'] + ind_files)

	print('Marking duplicates and sorting... ' + row['RunID'] + ': ' + str(datetime.datetime.now()))
	subprocess.run(['gatk', 'MarkDuplicatesSpark', '-I', unsorted_sam, '-O', fm_obj.localBamFile, '--tmp-dir', fm_obj.localTempDir, '-OBI', '--spark-runner', 'LOCAL'], stderr = open('TempErrors.txt', 'a'))

	# Remove remaining files
	subprocess.run(['rm','-f',unsorted_sam])

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

