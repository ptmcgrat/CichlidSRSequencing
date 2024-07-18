import subprocess, os, pdb, psutil, shutil
import pysam
import pandas as pd
from helper_modules.file_manager_Replacement import FileManager as FM
#from helper_modules.Timer import Timer
from multiprocessing import cpu_count

class AlignmentWorker():
	def __init__(self, genome, fm_obj):
		self.fm_obj = fm_obj
		self.genome = genome
		self.s_dt = fm_obj.s_dt

		self.fileManagers = {}
		
		self.uBam_files = {}

		sizes = {}
		for sample in fm_obj.samples:
			# Create sample file manager (need to keep them all in memory for parallelization)
			self.fileManagers[sample] = FM(genome)
			self.fileManagers[sample].createSampleFiles(sample)

			sub_dt = fm_obj.s_dt[fm_obj.s_dt.SampleID == sample]

			self.uBam_files[sample] = [self.fileManagers[sample].localReadsDir + x for x in sub_dt.FileLocations]
			sizes[sample] = sum([fm_obj.returnFileSize(x) for x in self.uBam_files[sample]])


		# Make sure there is enough
		total_sample_size = sum(sizes.values())
		free_memory = shutil.disk_usage(fm_obj.localMasterDir).free
		if 3*total_sample_size > free_memory:
			raise Exception('Need more space to run this analysis')

		self.samples = list({k: v for k, v in sorted(sizes.items(), key=lambda item: item[1], reverse = True)}.keys())
		print('The order of analysis based on size will be: ' + ',' + ','.join(self.samples))

	def monitorProcess(self,command,base_text,resource_file,error_file):
		
		timer = Timer()
		
		fm_obj = self.fileManagers[self.samples[i]]

		resource_fp = open(resource_file, 'w')
		error_fp= open(error_file, 'w')
		timer.start('    ' + base_text)
		print('cpu,threads,memory', file = resource_fp)

		p1 = subprocess.Popen(command, stderr = error_fp, stdout = subprocess.DEVNULL)
		proc = psutil.Process(pid = p1.pid)
		print(','.join([str(x) for x in [proc.cpu_percent(interval = 1), proc.num_threads(), proc.memory_info().rss/1000000000]]), file = resource_fp)

		while p1.poll() is None:
			try:
				print(','.join([str(x) for x in [proc.cpu_percent(interval = 60), proc.num_threads(), proc.memory_info().rss/1000000000]]), file = resource_fp)
			except psutil.ZombieProcess:
				break
			except psutil.NoSuchProcess:
				break
			data_file.flush()

		p1.communicate()

		data_file.close()
		dt = pd.read_csv(resource_file)
		mean = dt.mean()
		max_usage = dt.max()
		#print(' CPU_avg,max: ' + str(round(mean.cpu,1)) + ',' + str(round(max_usage.cpu,1)) + ' RAM_avg,max: ' + str(round(mean.memory,1)) + ',' + str(round(max_usage.memory,1)) + ' Threads: ' + str(mean.threads) + ' ')
		print(' CPU_avg,max: {:0.1f},{:0.1f} RAM_avg,max: {:0.1f},{:0.1f} Threads: {}.... '.format(mean.cpu, max_usage.cpu, mean.memory, max_usage.memory, mean.threads), end = '')
		timer.stop()

	def monitorProcesses(self, command_dict, base_text, num_parallel):
		
		resource_fp = open(self.fm_obj.localProcessesFile,'w')

		proc = psutil.Process(pid = os.getpid())

		print('cpu,threads,memory', file = resource_fp)

		error_files = []
		processes = []
		for strain,command in command_dict.items():
			fm_obj = self.fileManagers[strain]
			error_file = open(fm_obj.localSampleErrorDir + base_text + '_errors.txt', 'w')
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
		processes = []
		# Loop through all of the runs for a sample
		for sample in self.samples:
			fm_obj = self.fileManagers[sample]
			for uBam_file in self.uBam_files[sample]:
				if approach == 'Normal':
					fm_obj.downloadData(uBam_file)
				elif approach == 'Popen':
					processes.append(fm_obj.downloadData(uBam_file, parallel = True))
				elif approach == 'rclone_more_threads':
					fm_obj.downloadData(uBam_file, rclone = True)

		if approach == 'Popen':
			for p in processes:
				p.communicate()


	def delete_read_data(self):
		subprocess.run(['rm', '-f'] + self.downloaded_files)

	def alignData(self, linked = False):
		# Loop through all of the runs for a sample
		timer = Timer()
		for sample in self.samples:
			fm_obj = self.fileManagers[sample]

			sorted_bam = fm_obj.localTempSortedBamFile
			if os.path.isfile(sorted_bam):
				print(sample + ' already run')
				continue

			if linked:
				timer.start('  Aligning reads to create sorted Bam files')

			for i,uBam_file in enumerage(self.uBam_files[sample]):

				# Create temporary outputfile
				t_bam = fm_obj.localSampleTempDir + sample + '.' + str(i) + '.sorted.bam'

				if linked:
					command1 = ['gatk', 'SamToFastq', '-I', uBam_file, '--FASTQ', '/dev/stdout', '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2']
					command1 += ['--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', fm_obj.localSampleTempDir]
					command2 = ['bwa', 'mem', '-t', str(cpu_count()), '-M', '-p', fm_obj.localGenomeFile, '/dev/stdin']
					command3 = ['gatk', 'MergeBamAlignment', '-R', fm_obj.localGenomeFile, '--UNMAPPED_BAM', uBam_file, '--ALIGNED_BAM', '/dev/stdin']
					command3 += ['-O', t_bam, '--ADD_MATE_CIGAR', 'true', '--CLIP_ADAPTERS', 'false', '--CLIP_OVERLAPPING_READS', 'true']
					command3 += ['--INCLUDE_SECONDARY_ALIGNMENTS', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1', '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant']
					command3 += ['--ATTRIBUTES_TO_RETAIN', 'XS', '--TMP_DIR', fm_obj.localSampleTempDir]

					error_file = open(fm_obj.localSampleTempDir + 'Alignment_errors.txt', 'w')
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
					command1 = ['gatk', 'SamToFastq', '-I', uBam_file, '--FASTQ', fm_obj.localSampleTempDir + 'testing.fq', '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2']
					command1 += ['--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', fm_obj.localSampleTempDir]
					self.monitorProcesses({strain + str(i):command1}, 'SamToFastq_')
					#subprocess.run(command1)

					# Second command aligns fastq data to reference
					command2 = ['bwa', 'mem', '-t', str(cpu_count()), '-M', '-o', fm_obj.localSampleTempDir + 'testing.sam', '-p', fm_obj.localGenomeFile, fm_obj.localSampleTempDir + 'testing.fq']
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

