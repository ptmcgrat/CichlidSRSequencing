import subprocess, os, pdb, psutil
import pysam
import pandas as pd
from helper_modules.file_manager_Replacement import FileManager as FM
from helper_modules.Timer import Timer
from multiprocessing import cpu_count

class AlignmentWorker():
	def __init__(self, genome, sample_dt, samples):
		self.samples = samples
		self.fileManagers = {}
		self.sample_dt = sample_dt[sample_dt.SampleID.isin(samples)]
		for sample in samples:
			self.fileManagers[sample] = FM(genome)
			self.fileManagers[sample].createSampleFiles(sample)

	def monitorProcess(self,command,base_text):
		
		timer = Timer()
		
		data_file = open(self.fm_obj.localSampleTempDir + base_text + '_resources.txt', 'w')
		error_file = open(self.fm_obj.localSampleTempDir + base_text + '_errors.txt', 'w')

		timer.start('    ' + base_text)
		print('cpu,threads,memory', file = data_file)

		p1 = subprocess.Popen(command, stderr = error_file, stdout = subprocess.DEVNULL)
		proc = psutil.Process(pid = p1.pid)
		print(','.join([str(x) for x in [proc.cpu_percent(interval = 1), proc.num_threads(), proc.memory_info().rss/1000000000]]), file = data_file)

		while p1.poll() is None:
			try:
				print(','.join([str(x) for x in [proc.cpu_percent(interval = 60), proc.num_threads(), proc.memory_info().rss/1000000000]]), file = data_file)
			except psutil.ZombieProcess:
				break
			except psutil.NoSuchProcess:
				break
			data_file.flush()

		p1.communicate()

		data_file.close()
		dt = pd.read_csv(self.fm_obj.localSampleTempDir + base_text + '_resources.txt')
		mean = dt.mean()
		max_usage = dt.max()
		#print(' CPU_avg,max: ' + str(round(mean.cpu,1)) + ',' + str(round(max_usage.cpu,1)) + ' RAM_avg,max: ' + str(round(mean.memory,1)) + ',' + str(round(max_usage.memory,1)) + ' Threads: ' + str(mean.threads) + ' ')
		print(' CPU_avg,max: {:0.1f},{:0.1f} RAM_avg,max: {:0.1f},{:0.1f} Threads: {}.... '.format(mean.cpu, max_usage.cpu, mean.memory, max_usage.memory, mean.threads), end = '')
		timer.stop()

	def monitorProcesses(self,commands,base_text):
		

		timer = Timer()
		timer.start('   ' + base_text)

		proc = psutil.Process(pid = os.getpid())

		data_file = open(self.fm_obj.localSampleTempDir + base_text + '_resources.txt', 'w')
		print('cpu,threads,memory', file = data_file)

		error_files = []
		processes = []
		for i,command in enumerate(commands):
			self.fm_obj = self.fileManagers[self.samples[i]]
			error_file = open(self.fm_obj.localSampleTempDir + base_text + '_errors.txt', 'w')
			processes.append(subprocess.Popen(command, stderr = error_file, stdout = subprocess.DEVNULL))

		print(','.join([str(x) for x in [proc.cpu_percent(interval = 1), proc.num_threads(), proc.memory_info().rss/1000000000]]), file = data_file)
		while processes[0].poll() is None:
			try:
				print(','.join([str(x) for x in [proc.cpu_percent(interval = 60), proc.num_threads(), proc.memory_info().rss/1000000000]]), file = data_file)
			except ZombieProcess:
				break
			data_file.flush()

		for p in processes:
			p.communicate()

		timer.stop()
		data_file.close()
		dt = pd.read_csv(self.fm_obj.localSampleTempDir + base_text + '_resources.txt')


	def downloadReadData(self, approach = 'Normal'):
		self.downloaded_files = []
		processes = []

		# Loop through all of the runs for a sample
		for i, (index,row) in enumerate(self.sample_dt.iterrows()):
			fm_obj = self.fileManagers[row.SampleID]
			uBam_file = fm_obj.localReadsDir + row.FileLocations
			if approach == 'Normal':
				fm_obj.downloadData(uBam_file)
			elif approach == 'Popen':
				processes.append(fm_obj.downloadData(uBam_file, parallel = True))
			elif approach == 'rclone_more_threads':
				fm_obj.downloadData(uBam_file, rclone = True)
			self.downloaded_files.append(uBam_file)

		if approach == 'Popen':
			for p in processes:
				p.communicate()


	def delete_read_data(self):
		subprocess.run(['rm', '-f'] + self.downloaded_files)

	def alignData(self, linked = False):
		# Loop through all of the runs for a sample
		timer = Timer()
		for sample in self.samples:

			s_dt = self.sample_dt[self.sample_dt.SampleID == sample]
			self.fm_obj = self.fileManagers[sample]

			sorted_bam = self.fm_obj.localTempSortedBamFile
			if os.path.isfile(sorted_bam):
				print(sample + ' already run')
				continue

			if linked:
				timer.start('  Aligning reads to create sorted Bam files')


			for i, (index,row) in enumerate(s_dt.iterrows()):

				# Download unmapped bam file
				uBam_file = self.fm_obj.localReadsDir + row.FileLocations

				# Create temporary outputfile
				t_bam = self.fm_obj.localSampleTempDir + sample + '.' + str(i) + '.sorted.bam'

				if linked:
					command1 = ['gatk', 'SamToFastq', '-I', uBam_file, '--FASTQ', '/dev/stdout', '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2']
					command1 += ['--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', self.fm_obj.localSampleTempDir]
					command2 = ['bwa', 'mem', '-t', str(cpu_count()), '-M', '-p', self.fm_obj.localGenomeFile, '/dev/stdin']
					command3 = ['gatk', 'MergeBamAlignment', '-R', self.fm_obj.localGenomeFile, '--UNMAPPED_BAM', uBam_file, '--ALIGNED_BAM', '/dev/stdin']
					command3 += ['-O', t_bam, '--ADD_MATE_CIGAR', 'true', '--CLIP_ADAPTERS', 'false', '--CLIP_OVERLAPPING_READS', 'true']
					command3 += ['--INCLUDE_SECONDARY_ALIGNMENTS', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1', '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant']
					command3 += ['--ATTRIBUTES_TO_RETAIN', 'XS', '--TMP_DIR', self.fm_obj.localSampleTempDir]

					error_file = open(self.fm_obj.localSampleTempDir + 'Alignment_errors.txt', 'w')
					p1 = subprocess.Popen(command1, stdout=subprocess.PIPE, stderr = error_file)
					p2 = subprocess.Popen(command2, stdin = p1.stdout, stdout = subprocess.PIPE, stderr = error_file)
					p1.stdout.close()
					p3 = subprocess.Popen(command3, stdin = p2.stdout, stderr = error_file, stdout = subprocess.DEVNULL)
					p2.stdout.close()
					output = p3.communicate()

				else:
					# Align unmapped bam file following best practices
					# https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently
					# Align fastq files and sort them
					# First command coverts unmapped bam to fastq file, clipping out illumina adapter sequence by setting quality score to #
					# Debugging - useful for ensuring command is working properly, saving intermediate files instead of piping into each other
					command1 = ['gatk', 'SamToFastq', '-I', uBam_file, '--FASTQ', self.fm_obj.localSampleTempDir + 'testing.fq', '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2']
					command1 += ['--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', self.fm_obj.localSampleTempDir]
					self.monitorProcess(command1, 'SamToFastq_' + sample + '_' +str(i))
					#subprocess.run(command1)

					# Second command aligns fastq data to reference
					command2 = ['bwa', 'mem', '-t', str(cpu_count()), '-M', '-o', self.fm_obj.localSampleTempDir + 'testing.sam', '-p', self.fm_obj.localGenomeFile, self.fm_obj.localSampleTempDir + 'testing.fq']
					#print(command2)
					self.monitorProcess(command2, 'BWA_' + sample + '_' + str(i))
					subprocess.run(['rm', '-f', self.fm_obj.localSampleTempDir + 'testing.fq'])


					# Final command reads read group information to aligned bam file and sorts it
					# Figure out how to keep hard clipping
					

					# Debugging - useful for ensuring command is working properly, saving intermediate files instead of piping into each other
					command3 = ['gatk', 'MergeBamAlignment', '-R', self.fm_obj.localGenomeFile, '--UNMAPPED_BAM', uBam_file, '--ALIGNED_BAM', self.fm_obj.localSampleTempDir + 'testing.sam']
					command3 += ['-O', t_bam, '--ADD_MATE_CIGAR', 'true', '--CLIP_ADAPTERS', 'false', '--CLIP_OVERLAPPING_READS', 'true']
					command3 += ['--INCLUDE_SECONDARY_ALIGNMENTS', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1', '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant']
					command3 += ['--ATTRIBUTES_TO_RETAIN', 'XS', '--TMP_DIR', self.fm_obj.localSampleTempDir]
					self.monitorProcess(command3, 'SortMergeBam' + sample + '_' + str(i))
					subprocess.run(['rm', '-f', self.fm_obj.localSampleTempDir + 'testing.sam'])

				# Remove unmapped reads
				subprocess.run(['rm', '-f', uBam_file])
				
			if i == 0:
				subprocess.run(['mv', t_bam, sorted_bam])
			else:
				inputs = []
				ind_files = [self.fm_obj.localSampleTempDir + sample + '.' + str(x) + '.sorted.bam' for x in range(i+1)]
				for ind_file in ind_files:
					inputs = inputs + ['-I', ind_file]
				output = subprocess.run(['gatk', 'MergeSamFiles', '--TMP_DIR', self.fm_obj.localSampleTempDir] + inputs + ['-O', sorted_bam], stderr = open(self.fm_obj.localSampleTempDir + 'MergeSamFiles_errors.txt', 'w'), stdout = subprocess.DEVNULL)
				subprocess.run(['rm','-f'] + ind_files)
			
			if linked:
				timer.stop()
	def markDuplicates(self, parallel = False):
		commands = []
		del_files = []
		for sample in self.samples:
			s_dt = self.sample_dt[self.sample_dt.SampleID == sample]
			self.fm_obj = self.fileManagers[sample]
			unsorted_file = self.fm_obj.localTempSortedBamFile + 'unsorted_dedup.bam'

			command = ['gatk', 'MarkDuplicates', '-I', self.fm_obj.localTempSortedBamFile, '-O', self.fm_obj.localBamFile, '-M', self.fm_obj.localBamFile + '.duplication_metrics.txt', '--TMP_DIR', self.fm_obj.localSampleTempDir, '--CREATE_INDEX']
			#command = ['gatk', 'MarkDuplicatesSpark', '--create-output-bam-index',  '-I', self.fm_obj.localTempSortedBamFile, '-O', self.fm_obj.localBamFile, '-M', self.fm_obj.localBamFile + '.duplication_metrics.txt', '--tmp-dir', self.fm_obj.localSampleTempDir]
			
			#command2 = ['gatk', 'SortSam', '-I', unsorted_file, '-O', self.fm_obj.localBamFile, '-S', 'coordinate', '--TMP_DIR', self.fm_obj.localSampleTempDir, '--CREATE_INDEX']
			commands.append(command)
			if not parallel:
				self.monitorProcess(command, 'MarkDuplicates_' + sample)
				#self.monitorProcess(command2,'SortBam_' + sample)
				#subprocess.run(['rm', '-f', self.fm_obj.localTempSortedBamFile])
			else:
				del_files.append(self.fm_obj.localTempSortedBamFile)

		if parallel:
			self.monitorProcesses(commands, 'MarkDuplicates_' + str(len(self.samples)))
			#for del_file in del_files:
			#	subprocess.run(['rm','-f',del_file])

			command = ['gatk', 'MarkDuplicates', '--CREATE_INDEX',  '-I', self.fm_obj.localTempSortedBamFile, '-O', self.fm_obj.localBamFile, '-M', self.fm_obj.localBamFile + '.duplication_metrics.txt', '--TMP_DIR', self.fm_obj.localSampleTempDir]

	def splitBamfiles(self):
		for sample in self.samples:
			self.fm_obj = self.fileManagers[sample]
			print('  Splitting sample ' + sample)
			# Get contigs
			try:
				bam_obj = pysam.AlignmentFile(self.fm_obj.localBamFile)
			except OSError:
				print( '.........ERROR WITH THIS bam file. Probably truncated ' + sample)
				continue
			contigs = bam_obj.references  
			
			processes = []
			for contig in contigs:
				processes.append(subprocess.Popen(['python3', 'unit_scripts/split_bamfile_by_contig.py', self.fm_obj.localBamFile, contig]))

				if len(processes) == cpu_count():
					for p1 in processes:
						p1.communicate()
					processes = []
			for p1 in processes:
				p1.communicate()

			for bam_type in ['unmapped', 'discordant', 'inversion', 'duplication', 'clipped', 'chimeric']:
				bam_files = [self.fm_obj.localBamFile.replace('bam', x + '.' + bam_type + '.bam') for x in contigs]
				command = ['gatk', 'MergeSamFiles']
				for bam_file in bam_files:
					command += ['-I', bam_file]
				command += ['-O', self.fm_obj.localBamFile.replace('all.bam', bam_type + '.bam'), '--CREATE_INDEX']
				output = subprocess.run(command, capture_output = True)
				if output.returncode != 0:
					pdb.set_trace()
				for bam_file in bam_files:
					subprocess.run(['rm', bam_file])


	def createGVCF(self, parallel = False):

		commands = []
		for sample in self.samples:
			s_dt = self.sample_dt[self.sample_dt.SampleID == sample]
			self.fm_obj = self.fileManagers[sample]
			
			command = ['gatk', 'HaplotypeCaller', '-R', self.fm_obj.localGenomeFile, '-I', self.fm_obj.localBamFile, '-ERC', 'GVCF', '-O', self.fm_obj.localGVCFFile]
			commands.append(command)

			if not parallel:
				self.monitorProcess(command, 'HaplotypeCaller_' + sample)

		if parallel:
			self.monitorProcesses(commands, 'HaplotypeCaller_' + self.samples[0] + '_' + str(len(self.samples)) + 'OtherSamples')


	def calculateStats(self, sample):
		stats = {}
		self.fileManager = self.fileManagers[sample]
		for filename in [self.fileManager.localBamFile, self.fileManager.localUnmappedBamFile, self.fileManager.localDiscordantBamFile, self.fileManager.localInversionBamFile, self.fileManager.localDuplicationBamFile, self.fileManager.localClippedBamFile, self.fileManager.localChimericBamFile]:
			output = subprocess.run(['gatk', 'CountReads', '-I', filename], capture_output = True, encoding = 'utf-8')
			stats[filename.split('.')[-2]] = int(output.stdout.split('\n')[1])
		return stats

