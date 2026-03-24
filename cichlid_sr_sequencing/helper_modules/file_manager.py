import platform, os, pdb, gspread, subprocess
import pandas as pd

from gspread_dataframe import get_as_dataframe
from google.oauth2.service_account import Credentials


class FileManager():
	def __init__(self, genome_version = None, sampleID = None, rcloneRemote = 'ptm_dropbox:/', masterDir = 'COS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/'):

		self.genome_version = genome_version

		if platform.node() == 'ebb-utaka.biosci.gatech.edu' or platform.node() == 'utaka.biosci.gatech.edu' or 'utaka' in platform.node():
			self.localMasterDir = '/Data/' + os.getenv('USER') + '/Temp/CichlidSequencingData/'
		else:
			self.localMasterDir = os.getenv('HOME').rstrip('/') + '/' + 'Temp/CichlidSequencingData/' #Master directory for local data

		self.rcloneRemote = rcloneRemote
		self.cloudMasterDir = self.rcloneRemote + masterDir
	
		"""self.linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
							  'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
							  'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
							  'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

		"""
		self._createMasterDirs()
		self._readDatabases()
		if genome_version is not None:
			self.setGenome(genome_version)
			if sampleID is not None:
				self.createSampleFiles(sampleID)

	def _createMasterDirs(self):

		self.localPolymorphismsDir = self.localMasterDir + 'Polymorphisms/'	
		self.localReadsDir = self.localMasterDir + 'Reads/'		
		self.localSeqCoreDataDir = self.localMasterDir + 'SeqCoreData/'
		self.localBamfilesDir = self.localMasterDir + 'Bamfiles/'
		self.localGenomesDir = self.localMasterDir + 'Genomes/'
		self.localTempDir = self.localMasterDir + 'Temp/'
		self.localReadDownloadDir = self.localReadsDir + 'ReadDownloadFiles/'
		
		self.localCredentialFile = self.localMasterDir + 'cichlidsrsequencing_api_creds.json'

	def _readDatabases(self):
		g_ID = '1NmgB_TWoO01Qz2ufvECuZFkxXayhUsyu8wQGStVB_8k'
		self.downloadData(self.localCredentialFile)

		scopes = ['https://www.googleapis.com/auth/spreadsheets', 'https://www.googleapis.com/auth/drive']
		credentials = Credentials.from_service_account_file(self.localCredentialFile, scopes=scopes)
		gc = gspread.authorize(credentials)

		spreadsheet = gc.open_by_key(g_ID) # Or use open('Spreadsheet Name')
		
		worksheet = spreadsheet.worksheet('GenomeDatabase') # Access a specific sheet tab
		self.g_dt = get_as_dataframe(worksheet, evaluate_formulas=True)
		worksheet = spreadsheet.worksheet('SampleDatabase') # Access a specific sheet tab
		s_dt = get_as_dataframe(worksheet, evaluate_formulas=True)
		worksheet = spreadsheet.worksheet('DNAReads') # Access a specific sheet tab
		d_dt = get_as_dataframe(worksheet, evaluate_formulas=True)
		self.s_dt = pd.merge(s_dt,d_dt, on = 'SampleID')
		worksheet = spreadsheet.worksheet('AlignmentDatabase') # Access a specific sheet tab
		self.a_dt = get_as_dataframe(worksheet, evaluate_formulas=True)

	def _createGenomeFiles(self):
		self.localBamRefDir = self.localBamfilesDir + self.genome_version + '/'
		self.localGenomeDir = self.localGenomesDir + self.genome_version + '/'
		if self.genome_version == 'Mzebra_UMD2a':
			self.localGenomeFile = self.localGenomeDir + 'GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
		elif self.genome_version == 'Mzebra_GT3':
			self.localGenomeFile = self.localGenomeDir + 'Mzebra_GT3.fasta'
		elif self.genome_version == 'Mconophoros_GT1':
			self.localGenomeFile = self.localGenomeDir + 'anchored_kocher_E_Mchenga_conof_Male_contigs_hs_with_kocher_MC_female_molecules_mito_corrected.fasta'
		else:
			raise FileNotFoundError(self.genome_version + ' not an option')

	def returnOptions(self, datatype):
		if datatype == 'Genomes':
			return self.g_dt.GenomeID.to_list()
		if datatype == 'Samples':
			return self.s_dt.SampleID.unique().tolist()
		if datatype == 'Species':
			return self.s_dt.Species.unique().tolist()
		if datatype == 'ProjectIDs':
			return self.s_dt.ProjectID.unique().tolist()
			
	def setGenome(self, genome_version):
		self.genome_version = genome_version
		self._createGenomeFiles()
		
	def setSamples(self, projectIDs, sampleIDs, species, rerun):

		if projectIDs is not None:
			self.s_dt = self.s_dt[self.s_dt.ProjectID.isin(projectIDs)]

		if sampleIDs is not None:
			self.s_dt = self.s_dt[self.s_dt.SampleID.isin(sampleIDs)]

		if species is not None:
			self.s_dt = self.s_dt[self.s_dt.Ecogroup_PTM.isin(ecogroupIDs)]

		# Filter alignment database for requested genome version
		self.a_dt = self.a_dt[(self.a_dt.GenomeVersion == self.genome_version)]

		# Identify already run samples
		filter_set = set(self.a_dt.SampleID)
		already_run_samples = [x for x in set(self.s_dt.SampleID) if x in filter_set]
		samples = [for x in set(self.s_dt.SampleID) if x not in filter_set]
		
		if not rerun:
			if len(already_run_samples) > 0:
				print('The following samples have already been aligned to the genome and will not be rerun:')
				print(','.join(sorted(already_run_samples)))
			self.samples = samples
		else:
			if len(already_run_samples) > 0:
				print('The following samples have already been aligned to the genome and will be overwritten:')
				print(','.join(sorted(already_run_samples)))
			self.samples = already_run_samples + samples

		print('The following samples will be run:')
		print(','.join(sorted(self.samples)))

		self.s_dt = self.s_dt[self.s_dt.SampleID.isin(self.samples)]

	def createSampleFiles(self, sampleID):
		self.s_dt = self.s_dt[self.s_dt.SampleID == sampleID]

		self.localRawBamFiles = [self.localReadsDir + x + for x in self.s_dt.FileLocations.to_list()]

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
				if '.DS_Store' not in output.stderr:
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
				if '.DS_Store' not in output.stderr:
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







