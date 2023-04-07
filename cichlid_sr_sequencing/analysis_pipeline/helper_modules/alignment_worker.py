import subprocess, pysam, os, pdb
from multiprocessing import cpu_count

"""
impletement a platform argument. I can add this to the argparse block in alignFastQ.py and have it get passed here. 
"""



class AlignmentWorker():
	def __init__(self, fileManager, sample_dt, sampleID, platform):
		self.fileManager = fileManager
		self.sample_dt = sample_dt[sample_dt.SampleID == sampleID]
		self.sampleID = sampleID
		self.fileManager.createSampleFiles(self.sampleID)
		self.platform = platform
		os.makedirs(self.fileManager.localSampleBamDir, exist_ok = True)

	def downloadReadData(self):
		# Loop through all of the runs for a sample
		for i, (index,row) in enumerate(self.sample_dt.iterrows()):
			# Download unmapped bam file
			uBam_file = self.fileManager.localReadsDir + row.FileLocations
			self.fileManager.downloadData(uBam_file)

	def alignData(self):
		# Loop through all of the runs for a sample
		sorted_bam = self.fileManager.localTempDir + self.sampleID + '.sorted.bam'

		for i, (index,row) in enumerate(self.sample_dt.iterrows()):

			# Download unmapped bam file
			uBam_file = self.fileManager.localReadsDir + row.FileLocations
			self.fileManager.downloadData(uBam_file)

			# Create temporary outputfile
			t_bam = self.fileManager.localTempDir + self.sampleID + '.' + str(i) + '.sorted.bam'
			if self.platform == 'illumina':
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
			elif self.platform == 'pacbio':
				"""
				For Pacbio Reads:
				Seems like minimap2 only works on FASTA sequences too so I won't be able to run it on UBAMs just like for illumna reads
				Step 1: UBAM to Fastq 
				command1: gatk SamToFastq -I GGBC5251_MZ-003-m_l1.unmapped_marked_adapters.bam --FASTQ test.fastq --INCLUDE_NON_PF_READS true

				Unsure if we need the following extra commands: 
				--CLIPPING_ATTRIBUTE XT
				The attribute that stores the position at which the SAM record should be clipped  Default value: null. 
				--CLIPPING_ACTION 2
				The action that should be taken with clipped reads: 'X' means the reads and qualities should be trimmed at the clipped position;
				'N' means the bases should be changed to Ns in the clipped region; and any integer means that the base qualities should be set to that value in the clipped region.  Default value: null. 
				(So the attribute is where the SAM record should be clipped ad this is replaces by a str(2) value... I think this may have specifically been used due to the nature of the L1 and L2 reads in the BAMs from Illumina data)
				--Interleave true
				Will generate an interleaved fastq if paired, each line will have /1 or /2 to describe which end it came from  Default value: false. Possible values: {true, false}
				Set to False since we don't have paired reads.
				--NON_PF true
				not foudn in version 4.3.0.0 for some reason?
				--INCLUDE_NON_PF_READS true
				Use this flag instead to get the non_filtered reads to be included in the fastq
				-- TMP_DIR leave as fm_obj.lcoalTempDir

				command2: minimap2 -asm20 /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_GT1/Mzebra_GT1_v1.fna GGBC5251_MZ-003-m_l1.unmapped_marked_adapters.bam > test.sam

				-asm20 is needed to indicate the reads are PacBio HiFi/CCS reads. 

				"""
				# command for splitting UBAM to fastq file. Many options used for the illumina reads are exlcuded here. 
				command1 = ['gatk', 'SamToFastq', '-I', uBam_file, '--FASTQ', '/dev/stdout', '--INCLUDE_NON_PF_READS', 'true']
				# command for generating a sam file using minimap2 (v2.18 or earlier):
				# be sure that the outpit is piped into the next command instead of into
				command2 = ['minimap2', '-ax', 'asm20', self.fileManager.localGenomeFile, '/dev/stdin']
				command3 = ['samtools' 'view' ]
				command4 = ['gatk', 'MergeBamAlignment', '-R', self.fileManager.localGenomeFile, '--UNMAPPED_BAM', uBam_file, '--ALIGNED_BAM', '/dev/stdin']
				command4 += ['-O', t_bam, '--ADD_MATE_CIGAR', 'true', '--CLIP_ADAPTERS', 'false', '--CLIP_OVERLAPPING_READS', 'true']
				command4 += ['--INCLUDE_SECONDARY_ALIGNMENTS', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1', '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant']
				command4 += ['--ATTRIBUTES_TO_RETAIN', 'XS', '--TMP_DIR', self.fileManager.localTempDir]

				# the chaining commands part:
				p1 = subprocess.Popen(command1, stdout=subprocess.PIPE, stderr = subprocess.DEVNULL)
				p2 = subprocess.Popen(command2, stdin = p1.stdout, stdout = subprocess.PIPE, stderr = subprocess.DEVNULL)
				p1.stdout.close()
				p3 = subprocess.Popen(command3, stdin = p2.stdout, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
				p2.stdout.close()
				p4 = subprocess.Popen(command4, stdin = p3.stdout, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
				output = p4.communicate()
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
		for p1 in processes:
			p1.communicate()

		for bam_type in ['unmapped', 'discordant', 'inversion', 'duplication', 'clipped', 'chimeric']:
			bam_files = [self.fileManager.localBamFile.replace('bam', x + '.' + bam_type + '.bam') for x in contigs]
			command = ['gatk', 'MergeSamFiles']
			for bam_file in bam_files:
				command += ['-I', bam_file]
			command += ['-O', self.fileManager.localBamFile.replace('all.bam', bam_type + '.bam'), '--CREATE_INDEX']
			output = subprocess.run(command, capture_output = True)
			if output.returncode != 0:
				pdb.set_trace()
			for bam_file in bam_files:
				subprocess.run(['rm', bam_file])


	def createGVCF(self):
		# Get contigs
		bam_obj = pysam.AlignmentFile(self.fileManager.localBamFile)
		contigs = bam_obj.references  

		processes = []
		vcf_files = []
		for contig in contigs:

			command = ['gatk', 'HaplotypeCaller', '-R', self.fileManager.localGenomeFile, '-I', self.fileManager.localBamFile, '-ERC', 'GVCF']
			command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
			command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
			command += ['-A', 'DepthPerSampleHC', '-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']

			##contig=<ID=NC_036786.1,length=64916660>
			if contig == 'NC_036786.1':
				vcf_files.append(self.fileManager.localTempDir + self.sampleID + '_' + contig + '_1.g.vcf')
				vcf_files.append(self.fileManager.localTempDir + self.sampleID + '_' + contig + '_2.g.vcf')

				command1 = command + ['-L', contig + ':1-32400000', '-O', vcf_files[-2]]
				command2 = command + ['-L', contig + ':32400000-64916660', '-O', vcf_files[-1]]

				processes.append(subprocess.Popen(command1, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL))
				processes.append(subprocess.Popen(command2, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL))

			else:	

				vcf_files.append(self.fileManager.localTempDir + self.sampleID + '_' + contig + '.g.vcf')
				command = command + ['-L', contig , '-O', vcf_files[-1]]
				processes.append(subprocess.Popen(command, stderr = subprocess.PIPE, stdout = subprocess.PIPE))

			if len(processes) == int(cpu_count()/4):
				for p1 in processes:
					pdb.set_trace()
					p1.communicate()
				processes = []

		for p1 in processes:
			p1.communicate()

		command = ['gatk', 'MergeVcfs']
		for vcf_file in vcf_files:
			command += ['-I', vcf_file]
		command += ['-O', self.fileManager.localGVCFFile]
		output = subprocess.run(command, capture_output = True)
		if output.returncode != 0:
			pdb.set_trace()


	def calculateStats(self):
		stats = {}
		for filename in [self.fileManager.localBamFile, self.fileManager.localUnmappedBamFile, self.fileManager.localDiscordantBamFile, self.fileManager.localInversionBamFile, self.fileManager.localDuplicationBamFile, self.fileManager.localClippedBamFile, self.fileManager.localChimericBamFile]:
			output = subprocess.run(['gatk', 'CountReads', '-I', filename], capture_output = True, encoding = 'utf-8')
			stats[filename.split('.')[-2]] = int(output.stdout.split('\n')[1])
		return stats