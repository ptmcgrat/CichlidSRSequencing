import os, subprocess, pdb, random, datetime
import pandas as pd
import numpy as np
from xml.etree import ElementTree as ET
from multiprocessing import cpu_count
#from pysam import VariantFile
from collections import defaultdict

class FileManager():
	def __init__(self, rcloneRemote = 'cichlidVideo:', masterDir = 'McGrath/Apps/CichlidSequencingData/'):

		self.localMasterDir = os.getenv('HOME').rstrip('/') + '/' + 'Temp/CichlidSequencingData/' #Master directory for local data

		# Identify cloud directory for rclone
		self.rcloneRemote = rcloneRemote
		# On some computers, the first directory is McGrath, on others it's BioSci-McGrath. Use rclone to figure out which
		output = subprocess.run(['rclone', 'lsf', self.rcloneRemote + masterDir], capture_output = True, encoding = 'utf-8')
		if output.stderr == '':
			self.cloudMasterDir = self.rcloneRemote + masterDir
		else:
			output = subprocess.run(['rclone', 'lsf', self.rcloneRemote + 'BioSci-' + masterDir], capture_output = True, encoding = 'utf-8')
			if output.stderr == '':
				self.cloudMasterDir = self.rcloneRemote + 'BioSci-' + masterDir
			else:
				raise Exception('Cant find master directory (' + masterDir + ') in rclone remote (' + rcloneRemote + '')

		"""self.linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
							  'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
							  'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
							  'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

		"""
		self._createMasterDirs()

	def _createMasterDirs(self,version = 'Mzebra_UMD2a'):
		self.localGenomesDir = self.localMasterDir + 'Genomes/'
		self.localPolymorphismsDir = self.localMasterDir + 'Polymorphisms/'		
		self.localReadsDir = self.localMasterDir + 'Reads/'		
		self.localSeqCoreDataDir = self.localMasterDir + 'SeqCoreData/'
		self.localBamfilesDir = self.localMasterDir + 'Bamfiles/'
		self.localTempDir = self.localMasterDir + 'Temp/'
		self.localBamRefDir = self.localBamfilesDir + version + '/'
		self.localGenomeDir = self.localGenomesDir + version + '/'
		self.localGenomeFile = self.localGenomeDir + 'UMD2a_LG_only.fna'
		self.localSampleFile = self.localReadsDir + 'SampleDatabase.csv'
		#self.localSampleFile = self.localReadsDir + 'MCs_to_add.csv'


	def _createSeqCoreFiles(self, coreID):
		self.coreID = coreID
		self.localCoreDir = self.localSeqCoreDataDir + coreID + '/'
#		self.downloadData(self.localCoreDir)

		core_directories = self.returnCloudDirs(self.localCoreDir)
		for d in core_directories:
			c_files = self.returnCloudFiles(self.localCoreDir + d)
			if 'RunInfo.xml' in c_files and 'SampleSheet.csv' in c_files and 'RunParameters.xml' in c_files:
				self.localRunInfoFile = self.localCoreDir + d + '/RunInfo.xml'
				self.localSampleInfoFile = self.localCoreDir + d + '/SampleSheet.csv'
				self.localParametersInfo = self.localCoreDir + d + '/RunParameters.xml'
			c_dirs = self.returnCloudDirs(self.localCoreDir + d)
			if len(c_dirs) == 1 and 'FASTQ' in c_dirs[0] and 'Unindexed' not in d:
				self.fastqDir = self.localCoreDir + d + '/' + c_dirs[0]

		self.downloadData(self.localRunInfoFile)
		self.downloadData(self.localSampleInfoFile)
		self.downloadData(self.localParametersInfo)

	def _addMCData(self):
		data = subprocess.run(['rclone', 'lsf', self.localSeqCoreDataDir.replace(self.localMasterDir, self.cloudMasterDir) + 'bs_fastq'], capture_output = True, encoding = 'utf-8').stdout
		fastq_dirs = [x for x in data.split('\n') if '_ds' in x]
		
		samples = sorted(list(set([x.split('_')[1] for x in fastq_dirs])))

		out_dt = pd.DataFrame(columns = ['SampleID','Datatype','Date','Paired','RG','Files'])
		for sample_small in samples:
			sample = 'MC_' + sample_small

			print(sample)
			for fastq_dir in [x for x in fastq_dirs if sample_small in x]:
				data = subprocess.run(['rclone', 'size', self.localSeqCoreDataDir.replace(self.localMasterDir, self.cloudMasterDir) + 'bs_fastq/' + fastq_dir], capture_output = True, encoding = 'utf-8').stdout
				try:
					data_size = int(data.split(' Bytes')[0].split(' (')[-1])
				except ValueError:
					print(fastq_dir)
					continue
				if data_size > 10000000:
					data = subprocess.run(['rclone', 'lsf', self.localSeqCoreDataDir.replace(self.localMasterDir, self.cloudMasterDir) + 'bs_fastq/' + fastq_dir], capture_output = True, encoding = 'utf-8').stdout
					fq1, fq2 = [x for x in data.split('\n') if x.endswith('fastq.gz')]

					date = str(datetime.datetime.now().date())

					read_group = '@RG\\tID:' + fastq_dir.split('_')[1] + '.' + fastq_dir.split('_')[2] + '.' + sample + '\\tLB:' + sample + '\\tSM:' + sample + '\\tPL:ILLUMINA'
					out_dt = out_dt.append({'SampleID':sample,'Datatype': 'GenomicDNA', 'Date':date,'Paired':'True','RG':read_group, 'Files': 'MC/' + fq1 + ',,' + 'MC/' + fq2}, ignore_index=True)
					
					command = ['rclone', 'copy', self.localSeqCoreDataDir.replace(self.localMasterDir, self.cloudMasterDir) + 'bs_fastq/' + fastq_dir + fq1, self.localReadsDir.replace(self.localMasterDir, self.cloudMasterDir) + 'MC/']
					print(command)
					subprocess.run(command)
					command = ['rclone', 'copy', self.localSeqCoreDataDir.replace(self.localMasterDir, self.cloudMasterDir) + 'bs_fastq/' + fastq_dir + fq2, self.localReadsDir.replace(self.localMasterDir, self.cloudMasterDir) + 'MC/']
					print(command)
					subprocess.run(command)

		out_dt.to_csv(self.localReadsDir + 'MCs_to_add.csv')
		self.uploadData(self.localReadsDir + 'MCs_to_add.csv')
		pdb.set_trace()

	def _alignQTLData(self):
		# Make directories necessary for analysis
		subprocess.run(['mkdir',self.localMasterDir])

		subprocess.run(['mkdir',self.localTempDir])
		subprocess.run(['mkdir',self.localBamfilesDir])
		subprocess.run(['mkdir',self.localBamRefDir])

		# Download data necessary for analysis		

		print('Downloading genome and sample file')
		self.downloadData(self.localGenomeDir)
		subprocess.run(['bwa', 'index', self.localGenomeFile])
		self.downloadData(self.localSampleFile)
		
		dt = pd.read_csv(self.localSampleFile)
		for sample in dt.SampleID:
			existing_bams = subprocess.run(['rclone', 'lsf', self.localBamRefDir.replace(self.localMasterDir, self.cloudMasterDir)], capture_output = True, encoding = 'utf-8').stdout.split()
			if sample + '.bam' in existing_bams:
				print(sample + ' already analyzed.')
				continue

			print('Processing sample: ' + sample)
			outbam = self.localBamRefDir + sample + '.bam'

			sample_dt = dt[dt.SampleID == sample]
			tfiles = []

			# Loop through all of the fastq files
			for index,row in sample_dt.iterrows():
				tfile1 = self.localTempDir + str(random.randint(10000000,99999999)) + '.sam'
				tfile2 = self.localTempDir + str(random.randint(10000000,99999999)) + '.bam'

				# Download fastq files
				fq1 = self.localReadsDir + row.Files.split(',,')[0]
				fq2 = self.localReadsDir + row.Files.split(',,')[1]
				self.downloadData(fq1)
				self.downloadData(fq2)

				# Align fastq files and sort them
				subprocess.run(['bwa', 'mem', '-t', str(cpu_count()), '-R', row.RG.replace('\t','\\t'), '-M', self.localGenomeFile, fq1, fq2], stdout = open(tfile1, 'w'),stderr = open('TempErrors.txt', 'a'))
				subprocess.run(['picard', 'SortSam', 'I='+tfile1, 'O='+tfile2, 'SORT_ORDER=coordinate', 'TMP_DIR=~/Temp/', 'VALIDATION_STRINGENCY=LENIENT'], stderr = open('TempErrors.txt', 'a'))
				
				# Remove files and append output
				subprocess.run(['rm','-f', tfile1, fq1, fq2])
				tfiles.append(tfile2)
			
			# Merge bam files if necessary
			tfile3 = self.localTempDir + str(random.randint(10000000,99999999)) + '.bam'
			if len(tfiles) > 1:
				subprocess.call(['samtools', 'merge', '-f', tfile3] + tfiles)
			else:
				subprocess.call(['mv', tfiles[0], tfile3])

			# Remove reads that aren't perfectly aligned			
			merged_bam = pysam.AlignmentFile(tfile3)
			good_bam = pysam.AlignmentFile(outbam, mode = 'wb', template = merged_bam)

			for read in merged_bam.fetch(until_eof=True):
				if not read.is_paired: 
					continue
				elif read.is_unmapped or read.mate_is_unmapped:
					continue
				elif not read.is_proper_pair:
					continue
				elif abs(read.isize) > 3000:
					continue
				elif read.mapq < 40:
					continue
				elif len(read.cigartuples) != 1:
					continue
				elif read.cigartuples[0][0] != 0:
					continue
				good_bam.write(read)
			
			good_bam.close()
	
			subprocess.call(['samtools', 'index', outbam])

			# Remove remailing files
			subprocess.run(['rm','-f'] + tfiles)
			subprocess.run(['rm','-f', tfile3])

			# Upload data and delete
			self.uploadData(outbam)
			self.uploadData(outbam + '.bai')
			subprocess.run(['rm','-f', outbam, outbam + '.bai'])

	def _filterMainFiles(self):
		mainPolymorphismFile = self.localPolymorphismsDir + 'UMD2a_all_annotated.vcf'
		mainGenomeFile = self.localGenomeDir + 'GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
		self.downloadData(mainPolymorphismFile)
		self.downloadData(mainGenomeFile)
		
		fasta_in = pysam.FastaFile(mainGenomeFile)
		
		fasta_out = open(self.localGenomeFile, 'w')

		for ref in fasta_in.references:
			if ref in self.linkageGroups:
				print('>' + ref, file = fasta_out)
				print(fasta_in[ref], file = fasta_out)

		fasta_in.close()
		self.uploadData(self.localGenomeDir + 'UMD2a_LG_only.fna')
		pass

		vcf_in = pysam.VariantFile(mainPolymorphismFile)  # auto-detect input format

		vcf_out = open(self.localPolymorphismsDir + 'UMD2a_StartingSNVs.bed', 'w')
		#print('##fileformat=VCFv4.2', file = vcf_out)
		#for lg in self.linkageGroups:
		#	print('##contig=<ID=' + lg + ',length=' + str(len(fasta_in[lg])), file = vcf_out)
		#print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT', file = vcf_out)
		for i,rec in enumerate(vcf_in.fetch()):

			if i % 100000 == 0: 
				print('Processing record ' + str(i))

			if rec.contig not in self.linkageGroups:
				continue

			if len(rec.alts)!=1 or len(rec.ref) != 1 or len(rec.alts[0]) != 1:
				continue

			if rec.qual < 300:
				continue

			CV = rec.samples.values()[7].alleles
			MC = rec.samples.values()[15].alleles
			TI = rec.samples.values()[26].alleles

			if CV[0] is None or MC[0] is None or TI[0] is None:
				continue
			if CV[0] != CV[1]:
				continue
			elif MC[0] != MC[1]:
				continue
			elif TI[0] != TI[1]:
				continue
			if CV == MC and TI == MC:
				continue
			#print('\t'.join([str(x) for x in [rec.contig,rec.pos,'.',rec.ref, rec.alts[0], rec.qual,'.','.' ,'.']]), file = vcf_out)
			print(rec.contig + '\t' + str(rec.pos-1) + '\t' + str(rec.pos), file = vcf_out)
		vcf_out.close()
		#vcftools --min-meanDP 20 --max-meanDP 60 --indv CV_all --indv MC_all --indv TI_all --remove-indels --minQ 300 --vcf UMD2a_all_annotated__CV_MC_TI_unique.vcf --out RILSNVs --recode

		self.uploadData(self.localPolymorphismsDir + 'UMD2a_StartingSNVs.bed')

	def _parallelMPileup(self, bedfile, referencefile, bamfiles, outfile):
		dt = pd.read_csv(bedfile, sep = '\t', header = None)
		
		processes = []
		for contig in sorted(set(dt[0])):
			if contig == 'NC_036795.1':
				continue
			contig_dt = dt[dt[0] == contig]
			contig_dt.to_csv(self.localPolymorphismsDir + contig + '.bed', index = False, header = False, sep = '\t')
			
			command = ['bcftools', 'mpileup', '-R', self.localPolymorphismsDir + contig + '.bed', '--fasta-ref', self.localGenomeFile, '-a', 'FORMAT/AD', '-o', self.localPolymorphismsDir + contig + '.vcf'] + bamfiles
			processes.append(subprocess.Popen(command))
		
		for process in processes:
			process.communicate()
		

		for contig in sorted(set(dt[0])):
			if contig == 'NC_036795.1':
				continue
			vcf_in = VariantFile(self.localPolymorphismsDir + contig + '.vcf')
			for i,rec in enumerate(vcf_in.fetch()):
				try:
					print(rec, file = vcf_out, end = '')
				except NameError:
					vcf_out = open(outfile, 'w')
					print(vcf_in.header, file = vcf_out, end = '')
					print(rec, file = vcf_out, end = '')

			vcf_in.close()
			#subprocess.call(['rm', '-f', self.localPolymorphismsDir + contig + '.vcf', self.localPolymorphismsDir + contig + '.bed'])
		vcf_out.close()

	def _identifyFixedDifferences(self):

		mainPolymorphismFile = self.localPolymorphismsDir + 'UMD2a_StartingSNVs.bed'
		refPolymorphismFile = self.localPolymorphismsDir + 'UMD2a_ReferenceSamples.vcf'

		self.downloadData(mainPolymorphismFile)
		self.downloadData(self.localGenomeFile)
		self.downloadData(self.localBamRefDir)

		sample_dict = {}
		sample_dict['CV'] = ['PM17_CV1','PM17_CV2','CVfem1', 'CVfem2','CVmale1']
		sample_dict['MC'] = ['PM17_MC1', 'PM17_MC2', 'MC_1B11', 'MC_1B18', 'MC_1B4', 'MC_1B5', 'MC_1C11', 'MC_1C17', 'MC_1C4', 'MC_1C5', 'MC_2B13', 'MC_2B17', 'MC_2B19', 'MC_2C10', 'MC_2C19', 'MC_3B2', 'MC_3B6', 'MC_3B7', 'MC_3B9', 'MC_3C2', 'MC_3C6', 'MC_3C7', 'MC_3C9', 'MC_4B12', 'MC_4B14', 'MC_4B25', 'MC_4C12', 'MC_4C13', 'MC_4C25', 'MC_5B22', 'MC_5B23', 'MC_5B24', 'MC_5B26', 'MC_5C22', 'MC_5C23', 'MC_5C24', 'MC_5C26']
		sample_dict['TI'] = ['PM17_TI1', 'PM17_TI2']

		bamfiles = [self.localBamRefDir + x + '.bam' for x in sample_dict['CV'] + sample_dict['MC'] + sample_dict['TI']]

		self._parallelMPileup(mainPolymorphismFile, self.localGenomeFile, bamfiles, refPolymorphismFile)

		vcf_in = VariantFile(refPolymorphismFile)  # auto-detect input format
		vcf_out_CV = open(self.localPolymorphismsDir + 'UMD2a_CVMC_SNVs.bed', 'w')
		vcf_out_TI = open(self.localPolymorphismsDir + 'UMD2a_TIMC_SNVs.bed', 'w')

		for i,rec in enumerate(vcf_in.fetch()):
			af = self._alleleFrequencies(rec, sample_dict)

			if af['TI'][1] < 4 or af['CV'][1] < 8 or af['MC'][1] < 20:
				continue

			cv_diff = abs(af['MC'][0]/af['MC'][1] - af['CV'][0]/af['CV'][1])
			ti_diff = abs(af['MC'][0]/af['MC'][1] - af['TI'][0]/af['TI'][1])

			if cv_diff == 1:
				print(rec.contig + '\t' + str(rec.pos-1) + '\t' + str(rec.pos), file = vcf_out_CV)
			if ti_diff == 1:
				print(rec.contig + '\t' + str(rec.pos-1) + '\t' + str(rec.pos), file = vcf_out_TI)

		vcf_out_CV.close()
		vcf_out_TI.close()

	def _alleleFrequencies(self, rec, sample_dict):
		out_af = {}
		geno_converter = {(0,0):np.array([0,2]), (0,1):np.array([1,2]), (1,1):np.array([2,2]), (None,):np.array([0,0])}

		for population,samples in sample_dict.items():
			out_af[population] = sum([self._sampleGenotyper(rec.samples[x]['AD'][0:2]) for x in samples])

			#out_af[population] = sum([geno_converter[rec.samples[x]['GT']] for x in samples])
		
		return out_af

	def _sampleGenotyper(self, data):
		reads = sum(data)
		ref_reads = data[0]
		alt_reads = data[1]
		if reads == 0:
			return np.array([0,0])
		elif ref_reads == 0:
			return np.array([2,2])
		elif alt_reads == 0:
			return np.array([0,2])
		else:
			return np.array([1,2])


	def _genotypeRILs(self):
		version = 'Mzebra_UMD2a'

		#self.downloadData(self.localPolymorphismFile)
		self.downloadData(self.localGenomeFile)
		self.downloadData(self.localBamRefDir)

		cv_bamfiles = [self.localBamRefDir + x for x in os.listdir(self.localBamRefDir) if x.endswith('.bam') and 'xCV_F2' in x] + [self.localBamRefDir + 'PM17_MC1.bam'] + [self.localBamRefDir + 'PM17_MC2.bam']
		ti_bamfiles = [self.localBamRefDir + x for x in os.listdir(self.localBamRefDir) if x.endswith('.bam') and 'TIxMC' in x] + [self.localBamRefDir + 'PM17_MC1.bam'] + [self.localBamRefDir + 'PM17_MC2.bam']

		#subprocess.run(['freebayes'] + [val for pair in zip(['--bam']*len(cv_bamfiles),cv_bamfiles) for val in pair] + ['--variant-input', self.localPolymorphismsDir + 'UMD2a_CVMC_SNVs.vcf', '--fasta-reference', self.localGenomeFile, '--vcf', self.localPolymorphismsDir + 'UMD2a_CVMC_SNVs_RILs.vcf', '--only-use-input-alleles'])
		#subprocess.run(['freebayes'] + [val for pair in zip(['--bam']*len(ti_bamfiles),ti_bamfiles) for val in pair] + ['--variant-input', self.localPolymorphismsDir + 'UMD2a_TIMC_SNVs2.vcf', '--fasta-reference', self.localGenomeFile, '--vcf', self.localPolymorphismsDir + 'UMD2a_TIMC_SNVs_RILs2.vcf', '--only-use-input-alleles'])

		self._parallelMPileup(self.localPolymorphismsDir + 'UMD2a_CVMC_SNVs.bed', self.localGenomeFile, cv_bamfiles, self.localPolymorphismsDir + 'UMD2a_CVMC_SNVs_RILs.vcf')
		self._parallelMPileup(self.localPolymorphismsDir + 'UMD2a_TIMC_SNVs.bed', self.localGenomeFile, ti_bamfiles, self.localPolymorphismsDir + 'UMD2a_TIMC_SNVs_RILs.vcf')


	def _transformReads(self, x):
		if x[0]/(x[0] + x[1]) <= 0.2:
			return 'MC_' + str(x)
		elif x[0]/(x[0] + x[1]) >= 0.8:
			return 'TI_' + str(x)
		else:
			return 'Het_' + str(x)


	def _analyzeVCF(self):
		
		
		vcf_in = VariantFile(self.localPolymorphismsDir + 'UMD2a_TIMC_SNVs_RILs.vcf')  # auto-detect input format

		converter_MCref = {(0,0):0,(0,1):1,(1,1):2,(None,):np.nan}
		converter_MCmut = {(0,0):2,(0,1):1,(1,1):0,(None,):np.nan}

		sample_dict = {'MC': ['PM17_MC1', 'PM17_MC2'], 'All':[x for x in vcf_in.header.samples]}

		data = defaultdict(list)
		count = 0
		for i,rec in enumerate(vcf_in.fetch()):
			af = self._alleleFrequencies(rec, sample_dict)
			dp = sum([np.array(rec.samples[x]['AD'][0:2]) for x in sample_dict['All']]).sum()

			if af['MC'][1] != 4:
				continue
			if af['MC'][0] == 4:
				reads = [np.array([rec.samples[x]['AD'][0],rec.samples[x]['AD'][1]]) for x in rec.samples.keys()]
			elif af['MC'][0] == 0:
				reads = [np.array([rec.samples[x]['AD'][1],rec.samples[x]['AD'][0]]) for x in rec.samples.keys()]
			else:
				continue

			#num_GT = np.count_nonzero(~np.isnan(genos))
			if (af['All'][0]/af['All'][1]) > 0.7 or (af['All'][0]/af['All'][1]) < 0.3:
				continue
			if dp < 100:
				continue
			if af['All'][1]/2 < 30:
				continue
			if rec.contig not in self.linkageGroups:
				continue

			data['Bin'].append(int(count/20))
			data['Contig'].append(self.linkageGroups[rec.contig])
			data['Loc'].append(rec.pos)
			data['DP'].append(dp)
			data['AF'].append(af['All'][0]/af['All'][1])
			data['num_GT'].append(af['All'][1]/2)
			for i,sample in enumerate([x for x in rec.samples.keys() if x != 'MC']):
				data[sample].append(reads[i])
			count += 1

			"""
			try:
				corr = np.ma.corrcoef(np.ma.masked_invalid(genos), np.ma.masked_invalid(prev_genos)).data[1,0]
			except NameError:
				corr = 1
			data['Corr'].append(corr)
			prev_genos = genos
			count += 1
			"""
		dt = pd.DataFrame(data)

		agg_dict = {x: 'sum' for x in dt.columns if 'x' in x}
		agg_dict['Loc'] = 'first'
		agg_dict['AF'] = 'mean'

		out_dt = dt.groupby(['Contig','Bin']).agg(agg_dict)
		for col in [x for x in out_dt.columns if x not in  ['Loc','AF']]:
			out_dt[col] = out_dt[col].apply(self._transformReads)

		dt = pd.merge(out_dt.reset_index()[['Contig','Bin','AF']], dt, on = ['Contig','Bin'])
		dt['AF_diff'] = np.abs(dt['AF_x'] - dt['AF_y'])
		dt = dt.rename(columns = {'AF_y':'AF'})
		dt = dt[dt.AF_diff < 0.1]

		dt.to_csv(self.localPolymorphismsDir + 'TIMCIndividualData.csv')
		self.uploadData(self.localPolymorphismsDir + 'TIMCIndividualData.csv')

		out_dt = dt.groupby(['Contig','Bin']).agg(agg_dict)
		for col in [x for x in out_dt.columns if x not in  ['Loc','AF']]:
			out_dt[col] = out_dt[col].apply(self._transformReads)


		out_dt.to_csv(self.localPolymorphismsDir + 'TIMCAggregateData.csv')
		self.uploadData(self.localPolymorphismsDir + 'TIMCAggregateData.csv')
		pdb.set_trace()

		replace_nones = {(None,None):(0,0)}
		vcf_in = VariantFile(self.localPolymorphismsDir + 'UMD2a_genotyped_secondfilterCV.recode.vcf')  # auto-detect input format
		with open(self.localPolymorphismsDir + 'UMD2a_genotyped_secondfilterCV.recode.txt', 'w') as f:
			print('Contig\tLocation\tBinnedLocation\tQuality\tAlleleFrequency\t' + '\t'.join(CV_specific_strains), file = f)
			for i,rec in enumerate(vcf_in.fetch()):
				if rec.samples['PM17_CV1']['GT'] == (0,0):
					genos = [(rec.samples[x]['RO'],rec.samples[x]['AO'][0]) for x in CV_specific_strains]
					genos = [(0,0) if x == (None,None) else x for x in genos]
				elif rec.samples['PM17_CV1']['GT'] == (1,1):
					genos = [(rec.samples[x]['AO'][0],rec.samples[x]['RO']) for x in CV_specific_strains]
					genos = [(0,0) if x == (None,None) else x for x in genos]
				else:
					print('Error!')
				af = sum([x[0] for x in genos]) / (sum([x[0] for x in genos]) + sum([x[1] for x in genos]))
				print(rec.chrom + '\t' + str(rec.start) + '\t' + str(500000*int(rec.start/500000)) + '\t' + str(rec.qual) + '\t' + str(af) + '\t' + '\t'.join([str(x) for x in genos]), file = f)

		vcf_in = VariantFile(self.localPolymorphismsDir + 'UMD2a_genotyped_secondfilterTI.recode.vcf')  # auto-detect input format
		with open(self.localPolymorphismsDir + 'UMD2a_genotyped_secondfilterTI.recode.csv', 'w') as f:

			for i,rec in enumerate(vcf_in.fetch()):
				genos = [allele_name1[rec.samples[x]['GT']] for x in TI_specific_strains]
				print(rec.chrom + ',' + str(rec.start) + ',' + ','.join(genos), file = f)

		dt = pd.read_csv(self.localPolymorphismsDir + 'UMD2a_genotyped_secondfilterCV.recode.txt', sep = '\t')

		for cn in [x for x in dt.columns if 'PM17' in x]:
			dt[cn] = dt[cn].apply(eval).apply(np.array)

		grouped = dt.groupby(['Contig','BinnedLocation'])

		agg_functions = {}

		for head in [x for x in dt.columns if 'PM17' in x]:
			agg_functions[head] = sum

		final_data = grouped.agg(agg_functions)
		
		for head in [x for x in dt.columns if 'PM17' in x]:
			final_data[head] = final_data[head].apply(lambda x: x/x.sum()).apply(lambda x: 'CV' if x[0] > 0.85 else ('Het' if x[0] > 0.15 else 'MC')) + final_data[head].astype(str)
	
		pdb.set_trace()




	def _addSeqCoreData(self, coreID, datatype, run_info_file = 'RunInfo.xml', sample_info_file = 'SampleSheet.csv', run_parameters_info = 'RunParameters.xml'):
		self.localCoreDir = self.localSeqCoreDataDir + coreID.rstrip('/') + '/'

		self.downloadData(self.localCoreDir)
		subprocess.run(['mkdir',self.localMasterDir])
		subprocess.run(['mkdir',self.localReadsDir])
		subprocess.run(['mkdir',self.localReadsDir + coreID])

		for root,dirs,files in os.walk(self.localCoreDir):
			if run_info_file in files:
				self.localRunInfoDir = root + '/' + run_info_file
			if sample_info_file in files:
				self.localSampleInfoDir = root + '/' + sample_info_file
			if run_parameters_info in files:
				self.localParametersInfoDir = root + '/' + run_parameters_info


		tree = ET.parse(self.localRunInfoDir).getroot()
		run_info = tree.findall('Run')[0]
		flowcell = run_info.findall('Flowcell')[0].text
		instrument = run_info.findall('Instrument')[0].text
		date = run_info.findall('Date')[0].text

		dt = pd.read_csv(self.localSampleInfoDir, header = 19)


		out_dt = pd.DataFrame(columns = ['SampleID','Datatype','Date','Paired','RG','Files'])
		for sample in dt.Sample_ID:
			for root,dirs,files in os.walk(self.localCoreDir):
				for d in dirs:
					if sample + '_L' in d:
						fqs = [x for x in os.listdir(root + '/' + d) if '.fastq.gz' in x]
						try:
							lane = [x for x in fqs[0].split('_') if x[0] == 'L'][0]
						except IndexError:
							print(d)
							continue
							#pdb.set_trace()
						read_group = '@RG\\tID:' + flowcell + '.' + lane + '.' + sample + '\\tLB:' + sample + '\\tSM:' + sample + '\\tPL:ILLUMINA'
						#print(read_group)

						try:
							out_dt = out_dt.append({'SampleID':sample,'Datatype': 'GenomicDNA', 'Date':date,'Paired':'True','RG':read_group, 'Files':coreID + '/' + fqs[0] + ',,' + coreID + '/' + fqs[1]}, ignore_index=True)
						except IndexError:
							pdb.set_trace()
						#print(sample + '\t' + str(date) + '\t' + 'True' + '\t' + read_group + '\t' + fqs[0] + ',,' + fqs[1])
						subprocess.run(['mv', root + '/' + d + '/' + fqs[0], self.localReadsDir + coreID])
						subprocess.run(['mv', root + '/' + d + '/' + fqs[1], self.localReadsDir + coreID])
						
						try:
							print('Uploading ' + fqs[0])
							self.uploadData(self.localReadsDir + coreID + '/' + fqs[0])
							self.uploadData(self.localReadsDir + coreID + '/' + fqs[1])
						except Exception:
							print('Cant upload ' + str(fqs[0]))
			if sample not in list(out_dt.SampleID):
				print('No sample data for: ' + sample)
		out_dt.to_csv(self.localReadsDir + coreID + '_ToBeAdded.csv')
		self.uploadData(self.localReadsDir + coreID + '_ToBeAdded.csv')

	def downloadData(self, local_data, tarred = False, tarred_subdirs = False):

		relative_name = local_data.rstrip('/').split('/')[-1] + '.tar' if tarred else local_data.rstrip('/').split('/')[-1]
		local_path = local_data.split(local_data.rstrip('/').split('/')[-1])[0]
		cloud_path = local_path.replace(self.localMasterDir, self.cloudMasterDir)

		cloud_objects = subprocess.run(['rclone', 'lsf', cloud_path], capture_output = True, encoding = 'utf-8').stdout.split()

		if relative_name + '/' in cloud_objects: #directory
			output = subprocess.run(['rclone', 'copy', cloud_path + relative_name, local_path + relative_name], capture_output = True, encoding = 'utf-8')
		elif relative_name in cloud_objects: #file
			output = subprocess.run(['rclone', 'copy', cloud_path + relative_name, local_path], capture_output = True, encoding = 'utf-8')
		else:
			raise FileNotFoundError('Cant find file for download: ' + cloud_path + relative_name)

		if not os.path.exists(local_path + relative_name):
			raise FileNotFoundError('Error downloading: ' + local_path + relative_name)

		if tarred:
			# Untar directory
			output = subprocess.run(['tar', '-xvf', local_path + relative_name, '-C', local_path], capture_output = True, encoding = 'utf-8')
			output = subprocess.run(['rm', '-f', local_path + relative_name], capture_output = True, encoding = 'utf-8')

		if tarred_subdirs:
			for d in [x for x in os.listdir(local_data) if '.tar' in x]:
				output = subprocess.run(['tar', '-xvf', local_data + d, '-C', local_data, '--strip-components', '1'], capture_output = True, encoding = 'utf-8')
				os.remove(local_data + d)

	def uploadData(self, local_data, tarred = False):

		relative_name = local_data.rstrip('/').split('/')[-1]
		local_path = local_data.split(relative_name)[0]
		cloud_path = local_path.replace(self.localMasterDir, self.cloudMasterDir)

		if tarred:
			output = subprocess.run(['tar', '-cvf', local_path + relative_name + '.tar', '-C', local_path, relative_name], capture_output = True, encoding = 'utf-8')
			if output.returncode != 0:
				print(output.stderr)
				raise Exception('Error in tarring ' + local_data)
			relative_name += '.tar'

		if os.path.isdir(local_path + relative_name):
			output = subprocess.run(['rclone', 'copy', local_path + relative_name, cloud_path + relative_name], capture_output = True, encoding = 'utf-8')
			#subprocess.run(['rclone', 'check', local_path + relative_name, cloud_path + relative_name], check = True)

		elif os.path.isfile(local_path + relative_name):
			#print(['rclone', 'copy', local_path + relative_name, cloud_path])
			output = subprocess.run(['rclone', 'copy', local_path + relative_name, cloud_path], capture_output = True, encoding = 'utf-8')
			output = subprocess.run(['rclone', 'check', local_path + relative_name, cloud_path], check = True, capture_output = True, encoding = 'utf-8')
		else:
			raise Exception(local_data + ' does not exist for upload')

		if output.returncode != 0:
			pdb.set_trace()
			raise Exception('Error in uploading file: ' + output.stderr)

	def returnCloudDirs(self, local_data):
		output = subprocess.run(['rclone', 'lsf', local_data.replace(self.localMasterDir, self.cloudMasterDir)], capture_output = True, encoding = 'utf-8')
		return [x.rstrip('/') for x in output.stdout.split('\n') if x.endswith('/') ]

	def returnCloudFiles(self, local_data):
		output = subprocess.run(['rclone', 'lsf', local_data.replace(self.localMasterDir, self.cloudMasterDir)], capture_output = True, encoding = 'utf-8')
		return [x.rstrip('/') for x in output.stdout.split('\n') if not x.endswith('/') ]


fm_obj = FileManager()
#fm_obj._addMCData()
#fm_obj._addSeqCoreData('PM17', 'GenomicDNA')
#fm_obj._alignQTLData()
#fm_obj._filterMainFiles()
#fm_obj._identifyFixedDifferences()
#fm_obj._runRILData()
#fm_obj._filterVCF()
#fm_obj._genotypeRILs()
fm_obj._analyzeVCF()
