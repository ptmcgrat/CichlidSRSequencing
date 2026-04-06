import subprocess, os, pdb, psutil, shutil
import pysam
import pandas as pd
from helper_modules.file_manager import FileManager as FM
from helper_modules.Timer import Timer
from multiprocessing import cpu_count

class AlignmentWorker():
	def __init__(self, genome, fm_obj, check_size = True):
		self.fm_obj = fm_obj
		self.fm_obj.readSampleDatabase()
		self.fm_obj.readAlignmentDatabase()
		
		self.genome = genome

		self.uBam_files = {}

		self.max_processes = int(psutil.virtual_memory().available/4000000000)
		sizes = {}

		for sampleID in fm_obj.samples:
			# Create sample file manager (need to keep them all in memory for parallelization)
			self.fm_obj.createSampleFiles(sampleID)

			#sub_dt = fm_obj.s_dt[fm_obj.s_dt.SampleID == sampleID]
			self.uBam_files[sampleID] = self.fm_obj.localRawBamFiles
			if check_size:
				try:
					sizes[sampleID] = fm_obj.merged_dt[fm_obj.merged_dt.SampleID == sampleID].FileSize.sum()
				except ValueError:
					pdb.set_trace()
		if check_size:
			# Make sure there is enough room
			total_sample_size = sum(sizes.values())
			free_memory = shutil.disk_usage(fm_obj.localMasterDir).free
			if 3*total_sample_size > free_memory:
				print('Total_sample_size: ' + str(3*total_sample_size) + ', Free memory: ' + str(free_memory))
				raise Exception('Need more space to run this analysis')
			self.samples = list({k: v for k, v in sorted(sizes.items(), key=lambda item: item[1], reverse = True)}.keys())
			print('The order of analysis based on size will be: ' + ',' + ','.join(self.samples))
		else:
			self.samples = fm_obj.samples

	def monitorProcess(self,command,base_text,resource_file,error_file):
		
		timer = Timer()
		
		fm_obj = self.fm_obj
		fm_obj.createSampleFiles(self.samples[i])

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
		timer = Timer()
		timer.start('    ' + base_text)

		error_files = {}
		processes = {}
		current_samples = []
		error_samples = []
		for sample,command in command_dict.items():
			fm_obj = self.fm_obj
			fm_obj.createSampleFiles(strain)
			error_files[sample] = fm_obj.localErrorsDir + base_text + '_' + sample + '_errors.txt'
			error_fp = open(error_files[sample], 'w')
			processes[sample] = subprocess.Popen(command, stderr = error_fp, stdout = subprocess.DEVNULL)
			current_samples.append(sample)
			if len(current_samples) == num_parallel:

				for sample, p in processes.items():
					p.communicate()
					if p.returncode != 0:
						print('  Failure of command for sample: ' + sample)
						error_samples.append(sample)
					else:
						subprocess.run(['rm',error_files[sample]])

				current_samples = []

		for sample, p in processes.items():
			p.communicate()
			if p.returncode != 0:
				print('  Failure of command for sample: ' + sample)
				error_samples.append(sample)
			else:
				subprocess.run(['rm',error_files[sample]])

		timer.stop()
		return error_samples

	def downloadReadData(self):
		processes = []
		# Loop through all of the runs for a sample
		for sample in self.samples:
			fm_obj = self.fm_obj
			fm_obj.createSampleFiles(sample)
			for uBam_file in self.uBam_files[sample]:
				processes.append(fm_obj.downloadData(uBam_file, parallel = True))
				
		for p in processes:
			p.communicate()

	def delete_read_data(self):
		subprocess.run(['rm', '-f'] + self.downloaded_files)

	def alignData(self):
		# Loop through all of the runs for a sample
		timer = Timer()
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
				if output.returncode != 0:
					pass
					# Remove unmapped reads
					#subprocess.run(['rm', '-f', uBam_file])
				else:
					pdb.set_trace()

			if i == 0:
				subprocess.run(['mv', t_bam, sorted_bam])
			else:
				inputs = []
				ind_files = [fm_obj.localSampleTempDir + sample + '.' + str(x) + '.sorted.bam' for x in range(i+1)]
				for ind_file in ind_files:
					inputs = inputs + ['-I', ind_file]
				output = subprocess.run(['gatk', 'MergeSamFiles', '--TMP_DIR', fm_obj.localSampleTempDir] + inputs + ['-O', sorted_bam], stderr = open(fm_obj.localSampleTempDir + 'MergeSamFiles_errors.txt', 'w'), stdout = subprocess.DEVNULL)
				subprocess.run(['rm','-f'] + ind_files)
			
			timer.stop()

	def markDuplicates(self):
		commands = {}
		del_files = {}
		for sample in self.samples:
			fm_obj = self.fm_obj
			fm_obj.createSampleFiles(sample)

			command = ['gatk', '--java-options', '-Xmx2g','MarkDuplicates', '-I', fm_obj.localTempSortedBamFile, '-O', fm_obj.localBamFile, '-M', fm_obj.localBamFile + '.duplication_metrics.txt', '--TMP_DIR', fm_obj.localSampleTempDir, '--SORTING_COLLECTION_SIZE_RATIO', '.2','--CREATE_INDEX']
			commands[sample] = command
			del_files[sample] = fm_obj.localTempSortedBamFile

		bad_samples = self.monitorProcesses(commands, 'MarkDuplicates_' + str(len(self.samples)), self.max_processes)
		for sample in self.samples:
			if sample not in bad_samples:
				subprocess.run(['rm','-f',del_files[sample]])

	def splitBamfiles(self):
		for sample in self.samples:
			fm_obj = self.fm_obj
			fm_obj.createSampleFiles(sample)
			#print('  Splitting sample ' + sample)
			# Get contigs
			try:
				bam_obj = pysam.AlignmentFile(fm_obj.localBamFile)
			except OSError:
				print( '.........ERROR WITH THIS bam file. Probably truncated ' + sample)
				continue
			contigs = bam_obj.references  
			
			processes = []
			for contig in contigs:
				processes.append(subprocess.Popen(['python3', 'unit_scripts/split_bamfile_by_contig.py', fm_obj.localBamFile, contig]))

				if len(processes) == cpu_count():
					for p1 in processes:
						p1.communicate()
					processes = []
			for p1 in processes:
				p1.communicate()

			for bam_type in ['unmapped', 'discordant', 'inversion', 'duplication', 'clipped', 'chimeric']:
				bam_files = [fm_obj.localBamFile.replace('bam', x + '.' + bam_type + '.bam') for x in contigs]
				command = ['gatk', 'MergeSamFiles']
				for bam_file in bam_files:
					command += ['-I', bam_file]
				command += ['-O', fm_obj.localBamFile.replace('all.bam', bam_type + '.bam'), '--CREATE_INDEX']
				output = subprocess.run(command, capture_output = True)
				if output.returncode != 0:
					pdb.set_trace()
				for bam_file in bam_files:
					subprocess.run(['rm', bam_file])


	def createGVCF(self):
		fasta_obj = pysam.FastaFile(self.fm_obj.localGenomeFile)
		chromosomes = fasta_obj.references
		commands = {}
		if len(self.samples) < 20:
			split = True
		
		for sample in self.samples:
			fm_obj = self.fm_obj
			fm_obj.createSampleFiles(sample)
			if split:
				processes = []
				vcfs = []
				for chrom in chromosomes:
					contig_vcf = fm_obj.localGVCFFile.replace('.g.vcf.gz','_' + chrom + '.g.vcf.gz')
					command = ['gatk', '--java-options', '-Xmx2g', 'HaplotypeCaller', '-R', fm_obj.localGenomeFile, '-I', fm_obj.localBamFile, '-L', chrom, '-ERC', 'GVCF', '-O', contig_vcf]
					processes.append(subprocess.Popen(command, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL))
					vcfs.append(contig_vcf)
				for p1 in processes:
					p1.communicate()
				command = ['gatk','GatherVcfs']
				for chrom in chromosomes:
					contig_vcf = fm_obj.localGVCFFile.replace('.g.vcf.gz','_' + chrom + '.g.vcf.gz')
					command += ['-I', contig_vcf]
				command += ['-O', fm_obj.localGVCFFile]
				output = subprocess.run(command, capture_output = True)
				for vcf_file in vcfs:
					subprocess.run(['rm', vcf_file])

			else:
				commands[strain] = ['gatk', '--java-options', '-Xmx2g', 'HaplotypeCaller', '-R', fm_obj.localGenomeFile, '-I', fm_obj.localBamFile, '-ERC', 'GVCF', '-O', fm_obj.localGVCFFile]
		if not split:
			self.monitorProcesses(commands, 'HaplotypeCaller_' + str(len(self.samples)) + 'Samples', self.max_processes)

	def uploadAndUpdateDatabase(self, upload = True, sample_override = False):
		if sample_override:
			self.samples = [sample_override]
		for sample in self.samples:
			fm_obj = self.fm_obj
			fm_obj.createSampleFiles(sample)
			
			if not os.path.exists(fm_obj.localGVCFFile):
				print(sample + ' did not complete. You need to rerun.')
				continue

			if upload:
				fm_obj.uploadData(fm_obj.localSampleBamDir)

			# Calcultate stats
			stats = {}
			
			for filename in [fm_obj.localBamFile, fm_obj.localUnmappedBamFile, fm_obj.localDiscordantBamFile, fm_obj.localInversionBamFile, fm_obj.localDuplicationBamFile, fm_obj.localClippedBamFile, fm_obj.localChimericBamFile]:
				output = subprocess.run(['gatk', 'CountReads', '-I', filename], capture_output = True, encoding = 'utf-8')
				stats[filename.split('.')[-2]] = int(output.stdout.split('\n')[1])
		
			try:
				read_length = fm_obj.merged_dt[fm_obj.merged_dt['SampleID'] == sample]['ReadLength'].values[0]/2
			except IndexError:
				print('Weird Error. Somehow this Sample is not in the DNAReadsDatabase. Skipping...')
				continue
			reference_size = sum(pysam.FastaFile(fm_obj.localGenomeFile).lengths)
			coverage = stats['all'] * read_length / reference_size
			stats = {k:v/stats['all'] if k!= 'all' else v for k,v in stats.items()}
			
			sample_data = {'SampleID':sample, 'GenomeVersion': fm_obj.genome_version, 'RunIDs':',,'.join(list(fm_obj.reads_dt[fm_obj.reads_dt['SampleID'] == sample].RunID)), 
			   'Coverage':coverage, 'TotalReads':stats['all'], 'UnmappedReads':stats['unmapped'], 'DiscordantReads':stats['discordant'], 'InversionReads':stats['inversion'],
			   'DuplicationReads':stats['duplication'], 'ClippedReads':stats['clipped'], 'ChimericReads':stats['chimeric']}

			output = subprocess.run(['conda', 'list'], capture_output = True)
			sample_data['minimap2_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('minimap2')][0]
			sample_data['gatk_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('gatk4')][0]
			sample_data['pysam_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('pysam')][0]
			sample_data['BamSize'] = os.path.getsize(fm_obj.localBamFile)
	
			fm_obj.addAlignmentRow(sample_data)
			fm_obj.setAlignmentDatabase()

			subprocess.run(['rm','-rf', fm_obj.localSampleBamDir])
			subprocess.run(['rm','-rf', fm_obj.localSampleTempDir])
		

