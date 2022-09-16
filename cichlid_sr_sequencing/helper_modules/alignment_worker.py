import subprocess, pysam, os
from multiprocessing import cpu_count

class AlignmentWorker():
	def __init__(self, fileManager, sample_dt, sampleID):
		self.fileManager = fileManager
		self.sample_dt = sample_dt[sample_dt.SampleID == sampleID]
		self.sampleID = sampleID
		self.fileManager.createSampleFiles(self.sampleID)
		os.makedirs(self.fileManager.localSampleBamDir, exist_ok = True)

	def downloadReadData(self):
		# Loop through all of the runs for a sample
		for i, (index,row) in enumerate(self.sample_dt.iterrows()):
			# Download unmapped bam file
			uBam_file = self.fileManager.localReadsDir + row.File
			self.fileManager.downloadData(uBam_file)

	def alignData(self):
		# Loop through all of the runs for a sample
		sorted_bam = self.fileManager.localTempDir + self.sampleID + '.sorted.bam'

		for i, (index,row) in enumerate(self.sample_dt.iterrows()):

			# Download unmapped bam file
			uBam_file = self.fileManager.localReadsDir + row.File
			self.fileManager.downloadData(uBam_file)

			# Create temporary outputfile
			t_bam = self.fileManager.localTempDir + self.sampleID + '.' + str(i) + '.sorted.bam'

			# Align unmapped bam file following best practices
			# https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently
			# Align fastq files and sort them
			# First command coverts unmapped bam to fastq file, clipping out illumina adapter sequence by setting quality score to #
			command1 = ['gatk', 'SamToFastq', '-I', uBam_file, '--FASTQ', '/dev/stdout', '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2']
			command1 += ['--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', self.fileManager.localTempDir]

			# Debugging - useful for ensuring command is working properly, saving intermediate files instead of piping into each other
			#command1 = ['gatk', 'SamToFastq', '-I', uBam_file, '--FASTQ', fm_obj.localTempDir + 'testing.fq', '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2']
			#command1 += ['--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', fm_obj.localTempDir]
			#subprocess.run(command1)
			#pdb.set_trace()

			# Second command aligns fastq data to reference
			command2 = ['bwa', 'mem', '-t', str(cpu_count()), '-M', '-p', self.fileManager.localGenomeFile, '/dev/stdin']

			# Debugging - useful for ensuring command is working properly, saving intermediate files instead of piping into each other
			#command2 = ['bwa', 'mem', '-t', str(cpu_count()), '-M', '-p', fm_obj.localGenomeFile, fm_obj.localTempDir + 'testing.fq', '-o', fm_obj.localTempDir + 'testing.sam']
			#subprocess.run(command2)
			#pdb.set_trace()

			# Final command reads read group information to aligned bam file and sorts it
			# Figure out how to keep hard clipping
			command3 = ['gatk', 'MergeBamAlignment', '-R', self.fileManager.localGenomeFile, '--UNMAPPED_BAM', uBam_file, '--ALIGNED_BAM', '/dev/stdin']
			command3 += ['-O', t_bam, '--ADD_MATE_CIGAR', 'true', '--CLIP_ADAPTERS', 'false', '--CLIP_OVERLAPPING_READS', 'true']
			command3 += ['--INCLUDE_SECONDARY_ALIGNMENTS', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1', '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant']
			command3 += ['--ATTRIBUTES_TO_RETAIN', 'XS', '--TMP_DIR', self.fileManager.localTempDir]

			# Debugging - useful for ensuring command is working properly, saving intermediate files instead of piping into each other
			#command3 = ['gatk', 'MergeBamAlignment', '-R', fm_obj.localGenomeFile, '--UNMAPPED_BAM', uBam_file, '--ALIGNED_BAM', fm_obj.localTempDir + 'testing.sam']
			#command3 += ['-O', t_bam, '--ADD_MATE_CIGAR', 'true', '--CLIP_ADAPTERS', 'false', '--CLIP_OVERLAPPING_READS', 'true']
			#command3 += ['--INCLUDE_SECONDARY_ALIGNMENTS', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1', '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant']
			#command3 += ['--ATTRIBUTES_TO_RETAIN', 'XS', '--TMP_DIR', fm_obj.localTempDir]
			#subprocess.run(command3)
			#pdb.set_trace()

			# Figure out how to pipe 3 commands together
			p1 = subprocess.Popen(command1, stdout=subprocess.PIPE, stderr = subprocess.DEVNULL)
			p2 = subprocess.Popen(command2, stdin = p1.stdout, stdout = subprocess.PIPE, stderr = subprocess.DEVNULL)
			p1.stdout.close()
			p3 = subprocess.Popen(command3, stdin = p2.stdout, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
			p2.stdout.close()
			output = p3.communicate()
			# Remove unmapped reads
			subprocess.run(['rm', '-f', uBam_file])

		if i == 0:
			subprocess.run(['mv', t_bam, sorted_bam])
		else:
			inputs = []
			ind_files = [self.fileManager.localTempDir + self.sampleID + '.' + str(x) + '.sorted.bam' for x in range(i+1)]
			for ind_file in ind_files:
				inputs = inputs + ['-I', ind_file]
			output = subprocess.run(['gatk', 'MergeSamFiles', '--TMP_DIR', self.fileManager.localTempDir] + inputs + ['-O', sorted_bam], stderr = open('TempErrors.txt', 'a'), stdout = subprocess.DEVNULL)
			subprocess.run(['rm','-f'] + ind_files)
		

	def markDuplicates(self):
		sorted_bam = self.fileManager.localTempDir + self.sampleID + '.sorted.bam'
		output = subprocess.run(['gatk', 'MarkDuplicates', '-I', sorted_bam, '-O', self.fileManager.localBamFile, '-M', self.fileManager.localBamFile + '.duplication_metrics.txt', '--TMP_DIR', self.fileManager.localTempDir, '--CREATE_INDEX', 'true'], stdout = subprocess.DEVNULL, stderr = open('TempErrors.txt', 'a'))
		# Remove remaining files
		subprocess.run(['rm','-f',sorted_bam])

	def splitBamfiles(self):

		# Get contigs
		bam_obj = pysam.AlignmentFile(self.fileManager.localBamFile)
		contigs = bam_obj.references  
		
		processes = []
		for contig in contigs:
			processes.append(subprocess.Popen(['python3', 'unit_scripts/split_bamfile_by_contig.py', self.fileManager.localBamFile, contig]))

			if len(processes) == cpu_count():
				for p1 in processes:
					p1.communicate()
				processes = []
				pdb.set_trace()
		for p1 in processes:
			p1.communicate()



		for bam_type in ['unmapped', 'discordant', 'inversion', 'duplication', 'clipped', 'chimeric']:
			bam_files = [self.fileManager.localBamFile.replace('bam', x + bam_type + '.bam') for x in contigs]

			command = ['gatk', 'MergeSamFiles']
			for bam_file in bam_files:
				command += ['-I', bam_file]
			command += ['-O', self.fileManager.localBamFile.replace('bam', bam_type + '.bam'), '--CREATE_INDEX']
			subprocess.run(command)


	def createGVCF(self):
		# Get contigs
		bam_obj = pysam.AlignmentFile(self.fileManager.localBamFile)
		contigs = bam_obj.references  

		processes = []
		vcf_files = []
		for contig in contigs:
			vcf_files.append(self.fileManager.localTempDir + self.sampleID + '_' + contig + '.g.vcf')
			processes.append(subprocess.Popen(['gatk', '-R', self.fileManager.localGenomeFile, '-I', self.fileManager.localBamFile, '-ERC', 'GCVF', '-L', contig, '-O', vcf_files[-1]]))

			if len(processes) == cpu_count():
				for p1 in processes:
					p1.communicate()
				processes = []

		for p1 in processes:
			p1.communicate()

		command = ['gatk', 'MergeVcfs']
		for vcf_file in vcf_files:
			command += ['-I', vcf_file]
		command += ['-O', self.fileManager.localGVCFFile]

		subprocess.run(command)

