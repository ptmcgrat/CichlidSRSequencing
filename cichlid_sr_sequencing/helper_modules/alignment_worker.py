import subprocess, os, pdb, psutil, shutil
import pysam
import pandas as pd
from helper_modules.file_manager import FileManager as FM
from helper_modules.Timer import Timer
from multiprocessing import cpu_count
from types import SimpleNamespace

class AlignmentWorker():
	def __init__(self, genome, fm_obj, check_size = True):
		self.fm_obj = fm_obj
		self.fm_obj.readSampleDatabase()
		self.fm_obj.readAlignmentDatabase()
		
		self.genome = genome

		self.uBam_files = {}
		self.setMaxProcesses()

		for sampleID in fm_obj.samples:
			# Create sample file manager (need to keep them all in memory for parallelization)
			self.fm_obj.createSampleFiles(sampleID)
			self.uBam_files[sampleID] = self.fm_obj.localRawBamFiles
		self.samples = fm_obj.samples

	def setMaxProcesses(self):
		cpus = cpu_count()
		ram_units = int(psutil.virtual_memory().available/2500000000)
		if cpus > ram_units:
			self.max_processes = ram_units
			self.ram_unit = '2g'
		else:
			self.max_processes = cpus
			self.ram_unit = str(int(psutil.virtual_memory().available/cpus/1000000000)) + 'g'
		print(f'  Analysis using {self.max_processes} cores and {self.ram_unit} of RAM per core',)
	
	def checkSize(self):
		sizes = {}
		for sampleID in fm_obj.samples:
			try:
				sizes[sampleID] = fm_obj.merged_dt[fm_obj.merged_dt.SampleID == sampleID].FileSize.sum()
			except ValueError:
				pdb.set_trace()
		total_sample_size = sum(sizes.values())
		free_memory = shutil.disk_usage(fm_obj.localMasterDir).free
		if 3*total_sample_size > free_memory:
			print('Total_sample_size: ' + str(3*total_sample_size) + ', Free memory: ' + str(free_memory))
			raise Exception('Need more space to run this analysis')
		self.samples = list({k: v for k, v in sorted(sizes.items(), key=lambda item: item[1], reverse = True)}.keys())
		print('The order of analysis based on size will be: ' + ',' + ','.join(self.samples))
		
	def monitorProcesses(self, command_list, num_parallel = None):
		# command_dict is a dictionary where the key is the sample name that allows
		# the user to keep track of each separate command they want run
		# the values are a SimpleNamespace object containing command: a command as a list
		# and the location for an error file to be printed to (which is deleted if the
		# command runs sucessfully)

		bad_samples = []
		if num_parallel is None:
			num_parallel = self.max_processes
		for i,data in enumerate(command_list):
			data.error_fp = open(data.error_file, 'w')
			data.process = None
			if i < num_parallel:
				data.process = subprocess.Popen(data.command, stderr = data.error_fp, stdout = subprocess.DEVNULL)

		while command_list:
			finished_processes = [x for x in command_list if x.process is not None and x.process.poll() is not None]
			for data in finished_processes:
				print(f'..{data.sampleID} complete..', end = '')
				data.error_fp.close()
				if data.process.returncode != 0:
					bad_samples.append(data.sampleID)
				else:
					subprocess.run(['rm',data.error_file])
				command_list.remove(data)  # Remove finished process from monitoring list
				next_command = next((x for x in command_list if x.process is None),None)
				if next_command is not None:
					next_command.process = subprocess.Popen(next_command.command, stderr = next_command.error_fp, stdout = subprocess.DEVNULL)

		return bad_samples

	def downloadReadData(self):
		processes = {}
		bad_samples = []
		# Loop through all of the runs for a sample
		for sample in self.samples:
			fm_obj = self.fm_obj
			fm_obj.createSampleFiles(sample)
			for uBam_file in self.uBam_files[sample]:
				processes[sample] = fm_obj.downloadData(uBam_file, parallel = True)
				
		for sample, p in processes.items():
			p.communicate()
			if p.returncode != 0:
				bad_samples.append(sample)
		return bad_samples

	def delete_read_data(self):
		subprocess.run(['rm', '-f'] + self.downloaded_files)

	def alignData(self):
		# Loop through all of the runs for a sample
		timer = Timer()
		print('  Starting alignment of individual samples')
		bad_samples = []
		for sample in self.samples:
			fm_obj = self.fm_obj
			fm_obj.createSampleFiles(sample)
			os.makedirs(fm_obj.localSampleTempDir, exist_ok = True)
			os.makedirs(fm_obj.localSampleBamDir, exist_ok = True)

			sorted_bam = fm_obj.localTempSortedBamFile
			if os.path.isfile(sorted_bam):
				print(sample + ' already run')
				continue
			timer.start('  Aligning reads to create sorted Bam files for: ' + sample)

			for i,uBam_file in enumerate(self.uBam_files[sample]):

				# Create temporary outputfile
				t_bam = fm_obj.localSampleTempDir + sample + '.' + str(i) + '.sorted.bam'

				command1 = ['gatk', 'SamToFastq', '-I', uBam_file, '--FASTQ', '/dev/stdout', '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2']
				command1 += ['--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', fm_obj.localSampleTempDir]
				#command2 = ['bwa', 'mem', '-t', str(cpu_count()), '-M', '-p', fm_obj.localGenomeFile, '/dev/stdin']
				try:
					fm_obj.localMinimapGenomeFile
					command2 = ['minimap2', '-x','sr', '-a', '-t', str(cpu_count()), fm_obj.localMinimapGenomeFile, '-']
				except AttributeError:
					command2 = ['minimap2', '-x','sr', '-a', '-t', str(cpu_count()), fm_obj.localGenomeFile, '-']

				command3 = ['gatk', 'MergeBamAlignment', '-R', fm_obj.localGenomeFile, '--UNMAPPED_BAM', uBam_file, '--ALIGNED_BAM', '/dev/stdin']
				command3 += ['-O', t_bam, '--ADD_MATE_CIGAR', 'true', '--CLIP_ADAPTERS', 'false', '--CLIP_OVERLAPPING_READS', 'true']
				command3 += ['--INCLUDE_SECONDARY_ALIGNMENTS', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1', '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant']
				command3 += ['--ATTRIBUTES_TO_RETAIN', 'XS', '--TMP_DIR', fm_obj.localSampleTempDir]

				error_file_StF = open(fm_obj.localSampleTempDir + 'Alignment_errors_StF.txt', 'w')
				error_file_MM = open(fm_obj.localSampleTempDir + 'Alignment_errors_MM.txt', 'w')
				error_file_Merge = open(fm_obj.localSampleTempDir + 'Alignment_errors_Merge.txt', 'w')

				p1 = subprocess.Popen(command1, stdout=subprocess.PIPE, stderr = error_file_StF)
				p2 = subprocess.Popen(command2, stdin = p1.stdout, stdout = subprocess.PIPE, stderr = error_file_MM)
				p1.stdout.close()
				p3 = subprocess.Popen(command3, stdin = p2.stdout, stderr = error_file_Merge, stdout = subprocess.DEVNULL)
				p2.stdout.close()
				output = p3.communicate()
				if p3.returncode != 0:
					print('  Error aligning with sample ' + sample)
					bad_samples.append(sample)
					timer.stop()
					continue
				else:
					subprocess.run(['rm', '-f', uBam_file])
				
			if i == 0:
				subprocess.run(['mv', t_bam, sorted_bam])
				timer.stop()
				continue
			else:
				inputs = []
				ind_files = [fm_obj.localSampleTempDir + sample + '.' + str(x) + '.sorted.bam' for x in range(i+1)]
				for ind_file in ind_files:
					inputs = inputs + ['-I', ind_file]
				output = subprocess.run(['gatk', 'MergeSamFiles', '--TMP_DIR', fm_obj.localSampleTempDir] + inputs + ['-O', sorted_bam], stderr = open(fm_obj.localSampleTempDir + 'MergeSamFiles_errors.txt', 'w'), stdout = subprocess.DEVNULL)
				if output.returncode == 0 and os.path.exists(sorted_bam):
					subprocess.run(['rm','-f'] + ind_files)
				else:
					bad_samples.append(sample)
					
			timer.stop()
		return bad_samples

	def markDuplicates(self):
		commands = []

		for sample in self.samples:
			self.fm_obj.createSampleFiles(sample)
			error_file = self.fm_obj.localErrorsDir + 'MarkDuplicates_' + sample + '_errors.txt'
			command = ['gatk', '--java-options', '-Xmx'+self.ram_unit, 'MarkDuplicates', '-I', self.fm_obj.localTempSortedBamFile, '-O', self.fm_obj.localBamFile, '-M', self.fm_obj.localBamFile + '.duplication_metrics.txt', '--TMP_DIR', self.fm_obj.localSampleTempDir, '--SORTING_COLLECTION_SIZE_RATIO', '.2','--CREATE_INDEX']
			commands.append(SimpleNamespace(sampleID=sample, error_file = error_file, command = command))

		bad_samples = self.monitorProcesses(commands)
		
		for sample in self.samples:
			self.fm_obj.createSampleFiles(sample)
			if sample not in bad_samples:
				subprocess.run(['rm','-f',self.fm_obj.localTempSortedBamFile])
		return bad_samples

	def splitBamfiles(self):
		bad_samples = []
		commands = []
		for sample in self.samples:
			self.fm_obj.createSampleFiles(sample)
			try:
				bam_obj = pysam.AlignmentFile(self.fm_obj.localBamFile)
			except OSError:
				print( '.........ERROR WITH THIS bam file. Probably truncated ' + sample)
				bad_samples.append(sample)
				continue
			contigs = bam_obj.references  
			
			for contig in contigs:
				error_file = self.fm_obj.localErrorsDir + 'SplitBamFiles_' + sample +'__' + contig + '_errors.txt'
				command = ['python3', 'unit_scripts/split_bamfile_by_contig.py', self.fm_obj.localBamFile, contig]
				commands.append(SimpleNamespace(sampleID=sample + '__' + contig, error_file = error_file, command = command))

		tb_samples = self.monitorProcesses(commands)
		bad_samples = list(set(bad_samples + [x.split('__')[0] for x in tb_samples]))

		for sample in self.samples:
			if sample in bad_samples:
				continue
			for bam_type in ['unmapped', 'discordant', 'inversion', 'duplication', 'clipped', 'chimeric']:
				bam_files = [self.fm_obj.localBamFile.replace('bam', x + '.' + bam_type + '.bam') for x in contigs]
				command = ['gatk', 'MergeSamFiles']
				for bam_file in bam_files:
					command += ['-I', bam_file]
				command += ['-O', self.fm_obj.localBamFile.replace('all.bam', bam_type + '.bam'), '--CREATE_INDEX']
				output = subprocess.run(command, capture_output = True)
				if output.returncode != 0:
					bad_samples.append(sample)
				for bam_file in bam_files:
					subprocess.run(['rm', bam_file])
		return bad_samples

	def createGVCF(self):
		fasta_obj = pysam.FastaFile(self.fm_obj.localGenomeFile)
		contigs = fasta_obj.references
		commands = []
		bad_samples = []
		for sample in self.samples:
			self.fm_obj.createSampleFiles(sample)
			for contig in contigs:
				contig_vcf = self.fm_obj.localGVCFFile.replace('.g.vcf.gz','_' + contig + '.g.vcf.gz')
				error_file = self.fm_obj.localErrorsDir + 'HaplotypeCaller_' + sample + '__' + contig + '_errors.txt'
				command = ['gatk', '--java-options', '-Xmx' + self.ram_unit, 'HaplotypeCaller', '-R', self.fm_obj.localGenomeFile, '-I', self.fm_obj.localBamFile, '-L', contig, '-ERC', 'GVCF', '-O', contig_vcf]
				commands.append(SimpleNamespace(sampleID=sample + '__' + contig, error_file = error_file, command = command))
				
		tb_samples = self.monitorProcesses(commands)
		bad_samples = list(set(bad_samples + [x.split('__')[0] for x in tb_samples]))

		for sample in self.samples:
			if sample in bad_samples:
				continue
			self.fm_obj.createSampleFiles(sample)
			command = ['gatk','GatherVcfs']
			for contig in contigs:
				contig_vcf = self.fm_obj.localGVCFFile.replace('.g.vcf.gz','_' + contig + '.g.vcf.gz')
				command += ['-I', contig_vcf]
			command += ['-O', self.fm_obj.localGVCFFile]
			output = subprocess.run(command, capture_output = True)
			if output.returncode != 0:
				bad_samples.append(sample)
				continue
			for vcf_file in [self.fm_obj.localGVCFFile.replace('.g.vcf.gz','_' + x + '.g.vcf.gz') for x in contigs]:
				subprocess.run(['rm', vcf_file])

		return bad_samples

	def uploadAndUpdateDatabase(self, upload = True):
		commands = []

		for sample in self.samples:
			self.fm_obj.createSampleFiles(sample)
			error_file = self.fm_obj.localErrorsDir + 'UploadData_' + sample + '_errors.txt'
			command = self.fm_obj.uploadData(self.fm_obj.localSampleBamDir, parallel = True)
			commands.append(SimpleNamespace(sampleID=sample, error_file = error_file, command = command))

		bad_samples = self.monitorProcesses(commands, 1)

		for sample in self.samples:
			if sample in bad_samples:
				continue
			stats = {}
			
			for filename in [self.fm_obj.localBamFile, self.fm_obj.localUnmappedBamFile, self.fm_obj.localDiscordantBamFile, self.fm_obj.localInversionBamFile, self.fm_obj.localDuplicationBamFile, self.fm_obj.localClippedBamFile, self.fm_obj.localChimericBamFile]:
				output = subprocess.run(['gatk', 'CountReads', '-I', filename], capture_output = True, encoding = 'utf-8')
				try:
					stats[filename.split('.')[-2]] = int(output.stdout.split('\n')[1])
				except IndexError:
					pdb.set_trace()
	
			try:
				read_length = self.fm_obj.merged_dt[self.fm_obj.merged_dt['SampleID'] == sample]['ReadLength'].values[0]/2
			except IndexError:
				print('Weird Error. Somehow this Sample is not in the DNAReadsDatabase. Skipping...')
				continue
			reference_size = sum(pysam.FastaFile(self.fm_obj.localGenomeFile).lengths)
			coverage = stats['all'] * read_length / reference_size
			stats = {k:v/stats['all'] if k!= 'all' else v for k,v in stats.items()}
			
			sample_data = {'SampleID':sample, 'GenomeVersion': self.fm_obj.genome_version, 'RunIDs':',,'.join(list(self.fm_obj.reads_dt[self.fm_obj.reads_dt['SampleID'] == sample].RunID)), 
			   'Coverage':coverage, 'TotalReads':stats['all'], 'UnmappedReads':stats['unmapped'], 'DiscordantReads':stats['discordant'], 'InversionReads':stats['inversion'],
			   'DuplicationReads':stats['duplication'], 'ClippedReads':stats['clipped'], 'ChimericReads':stats['chimeric']}

			output = subprocess.run(['conda', 'list'], capture_output = True)
			sample_data['minimap2_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('minimap2')][0]
			sample_data['gatk_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('gatk4')][0]
			sample_data['pysam_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('pysam')][0]
			sample_data['BamSize'] = os.path.getsize(self.fm_obj.localBamFile)
	
			self.fm_obj.addAlignmentRow(sample_data)
			self.fm_obj.setAlignmentDatabase()

			#subprocess.run(['rm','-rf', self.fm_obj.localSampleBamDir])
			#subprocess.run(['rm','-rf', self.fm_obj.localSampleTempDir])
		

