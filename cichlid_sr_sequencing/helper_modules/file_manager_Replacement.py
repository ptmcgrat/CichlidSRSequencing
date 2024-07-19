import os, subprocess, pdb, random, datetime, platform
import pandas as pd
#import numpy as np
from xml.etree import ElementTree as ET
from multiprocessing import cpu_count
#from pysam import VariantFile
from collections import defaultdict

class FileManager():
	def __init__(self, genome_version = '', rcloneRemote = 'ptm_dropbox:/', masterDir = 'COS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/'):

		self.genome_version = genome_version

		if platform.node() == 'ebb-utaka.biosci.gatech.edu' or platform.node() == 'utaka.biosci.gatech.edu' or 'utaka' in platform.node():
			self.localMasterDir = '/Data/' + os.getenv('USER') + '/Temp/CichlidSequencingData/'
			self.localStorageDir = '/Output/'

		else:
			self.localMasterDir = os.getenv('HOME').rstrip('/') + '/' + 'Temp/CichlidSequencingData/' #Master directory for local data
			self.localStorageDir = '/Output/'

		# Identify cloud directory for rclone
		self.rcloneRemote = rcloneRemote
		self.cloudMasterDir = self.rcloneRemote + masterDir
	
		"""self.linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
							  'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
							  'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
							  'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

		"""
		self._createMasterDirs()
		

	def _createMasterDirs(self):
		self.localGenomesDir = self.localMasterDir + 'Genomes/'
		self.localGenomesComparisonDir = self.localMasterDir + 'Genomes/Comparisons/'

		self.localPolymorphismsDir = self.localMasterDir + 'Polymorphisms/'	
		self.localPileupDir = self.localMasterDir + '/Pileups/'	+ self.genome_version 	
	
		self.localReadsDir = self.localMasterDir + 'Reads/'		
		self.localSeqCoreDataDir = self.localMasterDir + 'SeqCoreData/'
		self.localBamfilesDir = self.localMasterDir + 'Bamfiles/'
		self.localTempDir = self.localMasterDir + 'Temp/'
		self.localAnalysisDir = self.localMasterDir + 'Analysis/'

		self.localBamRefDir = self.localBamfilesDir + self.genome_version + '/'
		self.localGenomeDir = self.localGenomesDir + self.genome_version + '/'
		if self.genome_version == 'Mzebra_UMD2a':
			self.localGenomeFile = self.localGenomeDir + 'GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
		elif self.genome_version == 'Mzebra_GT1':
			self.localGenomeFile = self.localGenomeDir + 'Mzebra_GT1_v1.fna'
		elif self.genome_version == 'Mzebra_HybridScaffold':
			self.localGenomeFile = self.localGenomeDir + 'HYBRID_SCAFFOLD.fasta'
		elif self.genome_version == 'Mzebra_GT3':
			self.localGenomeFile = self.localGenomeDir + 'Mzebra_GT3.fasta'
		elif self.genome_version == 'kocher_N_Met_zebra_Female':
			self.localGenomeFile = self.localGenomeDir + 'kocher_N_Met_zebra_Female_anchored_assembly.fasta.gz'

		elif self.genome_version == 'MZ4f_ptm':
			self.localGenomeFile = self.localGenomeDir + 'MZ4f_ptm_anchored_assembly.fasta.gz'
		elif self.genome_version == 'kocher_H_Aulon_yelhead_Female':
			self.localGenomeFile = self.localGenomeDir + 'A_spYH_GT1.fasta'
		elif self.genome_version == 'kocher_G_Aulon_yelhead_Male':
			self.localGenomeFile = self.localGenomeDir + 'kocher_G_Aulon_yelhead_Male_anchored_assembly.fasta.gz'
		elif self.genome_version == 'YH7f_ptm':
			self.localGenomeFile = self.localGenomeDir + 'YH7f_ptm_anchored_assembly.fasta.gz'

		elif self.genome_version == 'P_nyererei_v2':
			self.localGenomeFile = self.localGenomeDir + 'PunNye_v2_hyrbid_scaffold/PunNye_v2_hybrid_scaffold_genome.fasta'
		elif self.genome_version == 'O_niloticus_UMD_NMBU':
			self.localGenomeFile = self.localGenomeDir + 'GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna'

		elif self.genome_version == 'Rhamp_chilingali':
			self.localGenomeFile = self.localGenomeDir + 'GCA_963969265.1_fRhaChi2.1_genomic.fna'

		else:
			raise FileNotFoundError(self.genome_version + ' not an option')
		self.localSampleFile = self.localReadsDir + 'SampleDatabase.csv'
		self.localSampleFile_v2 = self.localReadsDir + 'SampleDatabase_v2.xlsx'

		self.localAlignmentFile = self.localBamfilesDir + 'AlignmentDatabase.csv'
		self.localReadDownloadDir = self.localReadsDir + 'ReadDownloadFiles/'

		self.localProcessesFile = self.localTempDir + 'ProcessInfo.csv'
		self.localErrorsDir = self.localMasterDir + 'Errors/'
		os.makedirs(self.localMasterDir, exist_ok = True)
		os.makedirs(self.localTempDir, exist_ok = True)
		os.makedirs(self.localBamRefDir, exist_ok = True)
		os.makedirs(self.localErrorsDir, exist_ok = True)


		#self.localSampleFile = self.localReadsDir + 'MCs_to_add.csv'

	def setSamples(self, projectIDs = None, sampleIDs = None, ecogroupIDs = None):

		self.downloadData(self.localSampleFile_v2)
		s_dt = pd.read_excel(self.localSampleFile_v2, sheet_name = 'SampleLevel')

		if projectIDs is not None:
			bad_projects = []
			for projectID in projectIDs:
				if projectID not in set(s_dt.ProjectID_PTM):
					bad_projects.append(projectID)
			if len(bad_projects) > 0:
				raise FileNotFoundError('The following projects were not found in sample database: ' + ','.join(bad_projects))

			s_dt = s_dt[s_dt.ProjectID_PTM.isin(projectIDs)]

		elif sampleIDs is not None:
			bad_samples = []
			for sample in sampleIDs:
				if sample not in list(s_dt.SampleID):
					bad_samples.append(sample)

			if len(bad_samples) > 0:
				raise FileNotFoundError('The following samples were not found in sample database: ' + ','.join(bad_samples))

			s_dt = s_dt[s_dt.ProjectID_PTM.isin(projectIDs)]

		elif ecogroupIDs is not None:
			bad_ecogroups = []
			for ecogroup in ecogroupIDs:
				if ecogroup not in set(s_dt.Ecogroup_PTM):
					bad_ecogroups.add(ecogroup)

			if len(bad_ecogroups) > 0:
				raise FileNotFoundError('Ecogroup ' + ecogroup + ' does not exist. Options are: ' + ','.join(set(s_dt.Ecogroup_PTM)))

			s_dt = s_dt[s_dt.Ecogroup_PTM.isin(ecogroupIDs)]

		else: 
			s_dt = s_dt

		# Download master alignment database to keep track of samples that have been aligned
		self.downloadData(self.localAlignmentFile)
		a_dt = pd.read_csv(self.localAlignmentFile)
		a_dt = a_dt[(a_dt.GenomeVersion == self.genome_version)]

		self.samples = set()
		already_run_samples = []
		for sample in set(s_dt.SampleID):
			if sample in set(a_dt.SampleID):
				already_run_samples.append(sample)
			else:
				self.samples.add(sample)

		if len(already_run_samples) > 0:
			print('The following samples have already been aligned to the genome and will not be rerun:')
			print(','.join(sorted(already_run_samples)))

		print('The following samples will be run:')
		print(','.join(sorted(self.samples)))

		self.downloadData(self.localSampleFile)
		s_dt = pd.read_csv(self.localSampleFile)

		self.s_dt = s_dt[s_dt.SampleID.isin(self.samples)]


	def createSampleFiles(self, sampleID):
		if not os.path.isfile(self.localSampleFile):
			print('DownloadingSampleFile')
			self.downloadData(self.localSampleFile)

		dt = pd.read_csv(self.localSampleFile)
		s_dt = dt[dt.SampleID == sampleID]
		projectID = s_dt.ProjectID.iloc[0]

		self.localRawBamFiles = [self.localReadsDir + projectID + '/' + x +'.unmapped_marked_adapters.bam' for x in s_dt.RunID.to_list()]

		self.sampleID = sampleID
		self.localSampleBamDir = self.localBamRefDir + sampleID + '/'
		self.localSampleTempDir = self.localTempDir + sampleID + '/'
		self.localTempSortedBamFile = self.localSampleTempDir + self.sampleID + '.sorted.bam'

		self.localBamFile = self.localSampleBamDir + sampleID + '.all.bam'
		self.localUnmappedBamFile = self.localSampleBamDir + sampleID + '.unmapped.bam'
		self.localDiscordantBamFile = self.localSampleBamDir + sampleID + '.discordant.bam'
		self.localInversionBamFile = self.localSampleBamDir + sampleID + '.inversion.bam'
		self.localDuplicationBamFile = self.localSampleBamDir + sampleID + '.duplication.bam'
		self.localClippedBamFile = self.localSampleBamDir + sampleID + '.clipped.bam'
		self.localChimericBamFile = self.localSampleBamDir + sampleID + '.chimeric.bam'
		self.localGVCFFile = self.localSampleBamDir + sampleID + '.g.vcf.gz'

		os.makedirs(self.localSampleBamDir, exist_ok = True)
		os.makedirs(self.localSampleTempDir, exist_ok = True)		

	def returnTempGVCFFile(self, contig):
		return self.localTempDir + contig + '_' + sampleID + '.g.vcf.gz'

	def returnTempBamFiles(self, contig):
		return [self.localTempDir + contig + '_' + sampleID + '.all.bam', self.localTempDir + contig + '_' + sampleID + '.unmapped.bam', self.localTempDir + contig + '_' + sampleID + '.discordant.bam', 
				self.localTempDir + contig + '_' + sampleID + '.inversion.bam', self.localTempDir + contig + '_' + sampleID + '.duplication.bam', self.localTempDir + contig + '_' + sampleID + '.clipped.bam',
				self.localTempDir + contig + '_' + sampleID + '.chimeric.bam']

	def createPileupFiles(self, sampleID):
		self.localSamplePileupDir = self.localPileupDir + sampleID + '/'
		self.localSampleSAMPileupFile = self.localSamplePileupDir + sampleID + '.mpileup'
		self.localSampleGVCFFile = self.localSamplePileupDir + sampleID + '.gvcf'

	def createAnalysisIDFiles(self, analysisID):
		self.localAnalysisFile = self.localAnalysisDir + analysisID + '.csv'

	def returnGenomeVersions(self):
		return self.returnCloudDirs(self.localGenomesDir)

	def returnSamples(self):
		if not os.path.isfile(self.localSampleFile):
			print('DownloadingSampleFile')
			self.downloadData(self.localSampleFile)
		dt = pd.read_csv(self.localSampleFile)
		return set(dt.SampleID)

	def downloadData(self, local_data, tarred = False, tarred_subdirs = False, parallel = False, rclone=False):

		relative_name = local_data.rstrip('/').split('/')[-1] + '.tar' if tarred else local_data.rstrip('/').split('/')[-1]
		local_path = local_data.split(local_data.rstrip('/').split('/')[-1])[0]
		cloud_path = local_path.replace(self.localMasterDir, self.cloudMasterDir)

		cloud_objects = subprocess.run(['rclone', 'lsf', cloud_path], capture_output = True, encoding = 'utf-8').stdout.split()

		if relative_name + '/' in cloud_objects: #directory
			output = subprocess.run(['rclone', 'copy', cloud_path + relative_name, local_path + relative_name], capture_output = True, encoding = 'utf-8')
		elif relative_name in cloud_objects: #file
			if parallel:
				process = subprocess.Popen(['rclone', 'copy', cloud_path + relative_name, local_path])
				return process
			elif rclone:
				output = subprocess.run(['rclone', 'copy', '--multi-thread-streams', '96', '--multi-thread-cutoff','100Mi', cloud_path + relative_name, local_path], capture_output = True, encoding = 'utf-8')
			else:
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

	def uploadData(self, local_data, tarred = False, upload_async = False):

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
			if upload_async:
				subprocess.Popen(['rclone', 'copy', local_path + relative_name, cloud_path + relative_name])
			else:
				output = subprocess.run(['rclone', 'copy', local_path + relative_name, cloud_path + relative_name], capture_output = True, encoding = 'utf-8')
			#subprocess.run(['rclone', 'check', local_path + relative_name, cloud_path + relative_name], check = True)

		elif os.path.isfile(local_path + relative_name):
			#print(['rclone', 'copy', local_path + relative_name, cloud_path])
			if upload_async:
				subprocess.Popen(['rclone', 'copy', local_path + relative_name, cloud_path])
			else:
				output = subprocess.run(['rclone', 'copy', local_path + relative_name, cloud_path], capture_output = True, encoding = 'utf-8')
				output = subprocess.run(['rclone', 'check', local_path + relative_name, cloud_path], check = True, capture_output = True, encoding = 'utf-8')
		else:
			raise Exception(local_data + ' does not exist for upload')

		if not upload_async:
			if output.returncode != 0:
				pdb.set_trace()
				raise Exception('Error in uploading file: ' + output.stderr)

	def returnFileSize(self, local_data):
		output = subprocess.run(['rclone', 'size', local_data.replace(self.localMasterDir, self.cloudMasterDir)], capture_output = True, encoding = 'utf-8')
		return int(output.stdout.split(' Byte)')[0].split('(')[-1])

	def returnCloudDirs(self, local_data):
		output = subprocess.run(['rclone', 'lsf', local_data.replace(self.localMasterDir, self.cloudMasterDir)], capture_output = True, encoding = 'utf-8')
		return [x.rstrip('/') for x in output.stdout.split('\n') if x.endswith('/') ]

	def returnCloudFiles(self, local_data):
		output = subprocess.run(['rclone', 'lsf', local_data.replace(self.localMasterDir, self.cloudMasterDir)], capture_output = True, encoding = 'utf-8')
		return [x.rstrip('/') for x in output.stdout.split('\n') if not x.endswith('/') ]

	def checkCloudFile(self, local_data):

		relative_name = local_data.rstrip('/').split('/')[-1]
		local_path = local_data.split(relative_name)[0]
		cloud_path = local_path.replace(self.localMasterDir, self.cloudMasterDir)

		uploadedFiles = self.returnCloudFiles(local_path)

		if relative_name in uploadedFiles:
			return True 
		else:
			return False







