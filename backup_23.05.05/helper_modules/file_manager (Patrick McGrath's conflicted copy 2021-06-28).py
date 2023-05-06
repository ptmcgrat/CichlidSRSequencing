import os, subprocess, pdb, random
import pandas as pd
from xml.etree import ElementTree as ET
from multiprocessing import cpu_count


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

		self._createMasterDirs()

	def _createMasterDirs(self):
		self.localGenomesDir = self.localMasterDir + 'Genomes/'
		self.localPolymorphismsDir = self.localMasterDir + 'Polymorphisms/'		
		self.localReadsDir = self.localMasterDir + 'Reads/'		
		self.localSeqCoreDataDir = self.localMasterDir + 'SeqCoreData/'
		self.localBamfilesDir = self.localMasterDir + 'Bamfiles/'


		self.localTempDir = self.localMasterDir + 'Temp/'

	def _runRILData(self):
		subprocess.run(['mkdir', self.localMasterDir])
		subprocess.run(['mkdir',self.localTempDir])

		version = 'Mzebra_UMD2a'
		self.localBamRefDir = self.localBamfilesDir + version + '/'
		self.localGenomeDir = self.localGenomesDir + version + '/'
		self.localGenomeFile = self.localGenomeDir + 'GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz'
		#self.localPolymorphismFile = self.localPolymorphismsDir + 'UMD2a_all_annotated.vcf'
		self.localSampleFile = self.localReadsDir + 'SampleDatabase.csv'
		self.downloadData(self.localGenomeFile)
		#self.downloadData(self.localPolymorphismFile)
		self.downloadData(self.localSampleFile)
		subprocess.run(['bwa', 'index', self.localGenomeFile])
		subprocess.run(['mkdir',self.localBamfilesDir])
		subprocess.run(['mkdir',self.localBamRefDir])

		dt = pd.read_csv(self.localSampleFile)
		samples = set(dt.SampleID)
		for sample in samples:
			outbam = self.localBamRefDir + sample + '.bam'

			sample_dt = dt[dt.SampleID == sample]
			tfiles = []

			for index,row in sample_dt.iterrows():
				tfile1 = self.localTempDir + str(random.randint(10000000,99999999)) + '.sam'
				tfile2 = self.localTempDir + str(random.randint(10000000,99999999)) + '.bam'

				fq1 = self.localReadsDir + row.Files.split(',,')[0]
				fq2 = self.localReadsDir + row.Files.split(',,')[1]
				self.downloadData(fq1)
				self.downloadData(fq2)
				subprocess.run(['bwa', 'mem', '-t', str(cpu_count()), '-R', row.RG.replace('\t','\\t'), '-M', self.localGenomeFile, fq1, fq2], stdout = open(tfile1, 'w'))
				subprocess.run(['rm','-f', fq1, fq2])

				subprocess.run(['picard', 'SortSam', 'I='+tfile1, 'O='+tfile2, 'SORT_ORDER=coordinate', 'TMP_DIR=~/Temp/'])

				#p1 = subprocess.Popen(['samtools', 'fixmate', '-m', '-O', 'BAM',tfile1, '-'], stdout=subprocess.PIPE)
				#p2 = subprocess.Popen(['samtools', 'sort','-o', tfile2, '-@', str(cpu_count()), '-'], stdin = p1.stdout)
				#p2.communicate()
				subprocess.run(['rm','-f', tfile1])
				tfiles.append(tfile2)
			tfile3 = self.localTempDir + str(random.randint(10000000,99999999)) + '.bam'

			if len(tfiles) > 1:
				subprocess.call(['samtools', 'merge', '-f', tfile3] + tfiles)
			else:
				subprocess.call(['mv', tfiles[0], tfile3])

			subprocess.run(['rm','-f'] + tfiles)
			subprocess.run(['picard', 'MarkDuplicates', 'I='+tfile3, 'O=' + outbam, 'METRICS_FILE=' + outbam + '.metrics', 'TMP_DIR=~/Temp/'])

			#subprocess.call(['samtools', 'markdup', '-s', tfile3, outbam])
			subprocess.call(['samtools', 'index', outbam])
			subprocess.run(['rm','-f', tfile3])
			self.uploadData(outbam)
			self.uploadData(outbam + '.bai')
			subprocess.run(['rm','-f', outbam, outbam + '.bai'])



	def _genotypeRILs(self):
		version = 'Mzebra_UMD2a'
		self.localBamRefDir = self.localBamfilesDir + version + '/'
		self.localGenomeDir = self.localGenomesDir + version + '/'
		self.localGenomeFile = self.localGenomeDir + 'GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz'
		self.localPolymorphismFile = self.localPolymorphismsDir + 'UMD2a_all_annotated__CV_MC_TI_unique.vcf'
		self.downloadData(self.localPolymorphismFile)
		self.downloadData(self.localGenomeFile)
		self.downloadData(self.localBamRefDir)

		bamfiles = [self.localBamRefDir + x for x in os.listdir(self.localBamRefDir) if x.endswith('.bam')]
		subprocess.run(['gatk','HaplotypeCaller', '--input'] + bamfiles +  ['--output', self.localPolymorphismsDir + 'UMD2a_genotypedRILsAll.vcf', '--reference', self.localGenomeFile, '--genotyping-mode', 'GENOTYPE_GIVEN_ALLELES', '--alleles', self.localPolymorphismFile])
		self.uploadData(self.localPolymorphismsDir + 'UMD2a_genotypedRILsAll.vcf')


	def _filterVCF(self):
		from pysam import VariantFile
		self.localPolymorphismFile = self.localPolymorphismsDir + 'UMD2a_all_annotated.vcf'
		self.downloadData(self.localPolymorphismFile)
		vcf_in = VariantFile(self.localPolymorphismFile)  # auto-detect input format
		vcf_out = VariantFile(self.localPolymorphismsDir + 'UMD2a_all_annotated__CV_MC_TI_unique.vcf', 'w', header=vcf_in.header)
		from collections import defaultdict
		counts = defaultdict(int)
		for rec in vcf_in.fetch():
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
			counts[rec.chrom]+=1
			vcf_out.write(rec)
		self.uploadData(self.localPolymorphismsDir + 'UMD2a_all_annotated__CV_MC_TI_unique.vcf')

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
						read_group = '@RG\\tID:' + flowcell + '.' + lane + '\\tLB:' + sample + '\\tSM:' + sample + '\\tPL:ILLUMINA'
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


fm_obj = FileManager()
fm_obj._runRILData()
#fm_obj._filterVCF()
#fm_obj._genotypeRILs()
