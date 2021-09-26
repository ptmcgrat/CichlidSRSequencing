import os, subprocess, pdb, random, pysam, datetime
import pandas as pd
import numpy as np
from xml.etree import ElementTree as ET
from multiprocessing import cpu_count
from pysam import VariantFile
from collections import defaultdict

class FileManager():
	def __init__(self, rcloneRemote = 'cichlidVideo:', masterDir = 'McGrath/Apps/CichlidSequencingData/'):

		self.localMasterDir = os.getenv('HOME').rstrip('/') + '/' + 'Temp/CichlidSequencingData/' #Master directory for local data

		# Identify cloud directory for rclone
		self.rcloneRemote = rcloneRemote
		# On some computers, the first directory is McGrath, on others it's BioSci-McGrath. Use rclone to figure out which
		self.cloudMasterDir = self.rcloneRemote + masterDir
		"""
		output = subprocess.run(['rclone', 'lsf', self.rcloneRemote + masterDir], capture_output = True, encoding = 'utf-8')
		if output.stderr == '':
			self.cloudMasterDir = self.rcloneRemote + masterDir
		else:
			output = subprocess.run(['rclone', 'lsf', self.rcloneRemote + 'BioSci-' + masterDir], capture_output = True, encoding = 'utf-8')
			if output.stderr == '':
				self.cloudMasterDir = self.rcloneRemote + 'BioSci-' + masterDir
			else:
				raise Exception('Cant find master directory (' + masterDir + ') in rclone remote (' + rcloneRemote + '')
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
		self.localGenomeFile = self.localGenomeDir + 'GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
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

	def _filterMainVCF(self):
		self.localPolymorphismFile = self.localPolymorphismsDir + 'UMD2a_all_annotated.vcf'
		self.downloadData(self.localPolymorphismFile)

		vcf_in = VariantFile(self.localPolymorphismFile)  # auto-detect input format
		vcf_out = VariantFile(self.localPolymorphismsDir + 'UMD2a_all_annotated__CV_MC_TI_unique.vcf', 'w', header=vcf_in.header)
		for i,rec in enumerate(vcf_in.fetch()):
			if i % 10000 == 0: 
				print('Processing record ' + str(i))
			if rec.alts[0] not in 'ACGT' or len(rec.alts)!=1:
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
			vcf_out.write(rec)

		vcf_out.close()
		#vcftools --min-meanDP 20 --max-meanDP 60 --indv CV_all --indv MC_all --indv TI_all --remove-indels --minQ 300 --vcf UMD2a_all_annotated__CV_MC_TI_unique.vcf --out RILSNVs --recode

		self.uploadData(self.localPolymorphismsDir + 'UMD2a_all_annotated__CV_MC_TI_unique.vcf')

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

	def _identifyFixedDifferences(self):

		self.localPolymorphismFile = self.localPolymorphismsDir + 'UMD2a_filtered_CV_MC_TI.vcf'

		self.downloadData(self.localPolymorphismFile)
		self.downloadData(self.localGenomeFile)
		self.downloadData(self.localBamRefDir)

		sample_dict = {}
		sample_dict['CV'] = ['CV','PM17_CV1','PM17_CV2','CVfem1', 'CVfem2','CVmale1']
		sample_dict['MC'] = ['MC','PM17_MC1', 'PM17_MC2', 'MC_1B11', 'MC_1B18', 'MC_1B4', 'MC_1B5', 'MC_1C11', 'MC_1C17', 'MC_1C4', 'MC_1C5', 'MC_2B13', 'MC_2B17', 'MC_2B19', 'MC_2C10', 'MC_2C19', 'MC_3B2', 'MC_3B6', 'MC_3B7', 'MC_3B9', 'MC_3C2', 'MC_3C6', 'MC_3C7', 'MC_3C9', 'MC_4B12', 'MC_4B14', 'MC_4B25', 'MC_4C12', 'MC_4C13', 'MC_4C25', 'MC_5B22', 'MC_5B23', 'MC_5B24', 'MC_5B26', 'MC_5C22', 'MC_5C23', 'MC_5C24', 'MC_5C26']
		sample_dict['TI'] = ['TI','PM17_TI1', 'PM17_TI2']

		bamfiles = [self.localBamRefDir + x + '.bam' for x in sample_dict['CV'] + sample_dict['MC'] + sample_dict['TI']]

		subprocess.run(['freebayes'] + [val for pair in zip(['--bam']*len(bamfiles),bamfiles) for val in pair] + ['--variant-input', self.localPolymorphismFile, '--fasta-reference', self.localGenomeFile, '--vcf', self.localPolymorphismsDir + 'UMD2a_genotypedReferenceStrains.vcf', '--only-use-input-alleles'])


	def _alleleFrequencies(self, rec, sample_dict):
		out_af = {}
		geno_converter = {(0,0):np.array([0,2]), (0,1):np.array([1,2]), (1,1):np.array([2,2]), (None,):np.array([0,0])}

		for population,samples in sample_dict.values(): 
			out_af[population] = sum([geno_converter[rec.samples[x]['GT']] for x in samples])
		
		return out_af

	def _genotypeRILs(self):
		version = 'Mzebra_UMD2a'
		self.localBamRefDir = self.localBamfilesDir + version + '/'
		self.localGenomeDir = self.localGenomesDir + version + '/'
		self.localGenomeFile = self.localGenomeDir + 'GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
		self.localPolymorphismFile = self.localPolymorphismsDir + 'UMD2a_filtered_CV_MC_TI.vcf'

		self.downloadData(self.localPolymorphismFile)
		self.downloadData(self.localGenomeFile)
		self.downloadData(self.localBamRefDir)

		bamfiles = [self.localBamRefDir + x for x in os.listdir(self.localBamRefDir) if x.endswith('.bam')]

		subprocess.run(['freebayes'] + [val for pair in zip(['--bam']*len(bamfiles),bamfiles) for val in pair] + ['--variant-input', self.localPolymorphismFile, '--fasta-reference', self.localGenomeFile, '--vcf', self.localPolymorphismsDir + 'UMD2a_genotypedRILsAll.vcf', '--only-use-input-alleles'])
		#subprocess.run(['gatk','HaplotypeCaller'] + [val for pair in zip(['--input']*len(bamfiles),bamfiles) for val in pair] +  ['--output', self.localPolymorphismsDir + 'UMD2a_genotypedRILsAll.vcf', '--reference', self.localGenomeFile, '--genotyping-mode', 'GENOTYPE_GIVEN_ALLELES', '--alleles', self.localPolymorphismFile])
		self.uploadData(self.localPolymorphismsDir + 'UMD2a_genotypedRILsAll.vcf')

	def _play(self):
		from pysam import VariantFile
		vcf_in = VariantFile(self.localPolymorphismsDir + 'UMD2a_genotypedRILsAllFixedHeaders.vcf')
		vcf_MC_in = VariantFile(self.localPolymorphismsDir + 'all_no_cr.vcf.gz')

		ref_strains = ['PM17_CV1', 'PM17_CV2', 'PM17_MC1', 'PM17_MC2', 'PM17_TI1', 'PM17_TI2']
		translate = {(0,0):-1, (0,1):0, (1,1):1, (None,):0}
		for i,rec in enumerate(vcf_in.fetch()):
			pdb.set_trace()
			print([rec.samples[x]['GT'] for x in ref_strains])

		
	def _analyzeVCF(self):
		from pysam import VariantFile
		from collections import defaultdict
		
		reference_strains = ['PM17_CV1','PM17_CV2','CV', 'CVfem1', 'CVfem2','CVmale1','PM17_MC1','PM17_MC2','MC','PM17_TI1','PM17_TI2']

		vcf_in = VariantFile(self.localPolymorphismsDir + 'UMD2a_genotypedRILsAllFixedHeaders.vcf')  # auto-detect input format
		vcf_out_CV = VariantFile(self.localPolymorphismsDir + 'UMD2a_genotyped_firstfilterCV.vcf', 'w', header=vcf_in.header)
		vcf_out_TI = VariantFile(self.localPolymorphismsDir + 'UMD2a_genotyped_firstfilterTI.vcf', 'w', header=vcf_in.header)


		for i,rec in enumerate(vcf_in.fetch()):
			genos = [rec.samples[x]['GT'] for x in reference_strains]
			CV_genos = [x for x in genos[:6] if x in [(0,0),(1,1)]]
			MC_genos = [x for x in genos[6:9] if x in [(0,0),(1,1)]]
			TI_genos = [x for x in genos[9:] if x in [(0,0),(1,1)]]

			if len(CV_genos) != 2 or len(MC_genos) != 2 or len(TI_genos) != 2:
				#print('Bad: ' + str(genos))
				continue

			CV_set = set(CV_genos)
			MC_set = set(MC_genos)
			TI_set = set(TI_genos)

			if len(CV_set) != 1 or len(MC_set) != 1 or len(TI_set) != 1:
				#print('Bad: ' + str(genos))
				continue

			if CV_set == MC_set and TI_set == MC_set:
				#print('Bad: ' + str(genos))
				continue
			if CV_set != MC_set:
				vcf_out_CV.write(rec)
			if TI_set != MC_set:
				vcf_out_TI.write(rec)

		CV_specific_strains = sorted([x for x in rec.samples.keys() if 'CV' in x])
		TI_specific_strains = sorted([x for x in rec.samples.keys() if 'TI' in x])

		command = ['vcftools', '--remove-indels', '--recode']
		for x in TI_specific_strains:
			command += ['--remove-indv', x]
		command += ['--maf', '0.3333', '--max-maf', '0.66666']
		command += ['--vcf', self.localPolymorphismsDir + 'UMD2a_genotyped_firstfilterCV.vcf']
		command += ['--out', self.localPolymorphismsDir + 'UMD2a_genotyped_secondfilterCV']

		subprocess.run(command)
		command = ['vcftools', '--remove-indels', '--recode']
		command += ['--maf', '0.3333', '--max-maf', '0.66666']

		for x in CV_specific_strains:
			command += ['--remove-indv', x]
		command += ['--vcf', self.localPolymorphismsDir + 'UMD2a_genotyped_firstfilterTI.vcf']
		command += ['--out', self.localPolymorphismsDir + 'UMD2a_genotyped_secondfilterTI']

		subprocess.run(command)

		allele_name1 = {(0,0):'CV', (0,1):'Het', (1,1):'MC', (None,):'NA'}
		allele_name2 = {(0,0):'MC', (0,1):'Het', (1,1):'CV', (None,):'NA'}

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
#fm_obj._play()
#fm_obj._alignQTLData()
fm_obj._identifyFixedDifferences()
#fm_obj._runRILData()
#fm_obj._filterVCF()
#fm_obj._genotypeRILs()
#fm_obj._analyzeVCF()
