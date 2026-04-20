import platform, os, pdb, gspread, subprocess
import pandas as pd

from gspread_dataframe import get_as_dataframe
from gspread_dataframe import set_with_dataframe
from google.oauth2.service_account import Credentials


class FileManager():
	def __init__(self, genome_version = None, sampleID = None, rcloneRemote = 'ptm_dropbox:/', masterDir = 'COS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/'):
		
		# Set master directories/remotes
		if platform.node() == 'ebb-utaka.biosci.gatech.edu' or platform.node() == 'utaka.biosci.gatech.edu' or 'utaka' in platform.node():
			self.localMasterDir = '/Data/' + os.getenv('USER') + '/Temp/CichlidSequencingData/'
		else:
			self.localMasterDir = os.getenv('HOME').rstrip('/') + '/' + 'Temp/CichlidSequencingData/' #Master directory for local data
		self.rcloneRemote = rcloneRemote
		self.cloudMasterDir = self.rcloneRemote + masterDir
	
		self.s_ID = '1NmgB_TWoO01Qz2ufvECuZFkxXayhUsyu8wQGStVB_8k' # Sample Database
		self.a_ID = '1vlA_eJ09RSKCaM0foFMc0n5tL95XK4deB3dr7kAFx3E' # Alignment Database

		# Create file structure for data
		self._createMasterDirs()

		# Authenticate for databases
		self._authenticateGS()

		# If genome and sample file is given then store and create file structures
		if genome_version is not None:
			self.setGenome(genome_version)
		if sampleID is not None:
			if genome_version is None:
				raise Exception('Cant set sampleID without setting genome version')
			self.createSampleFiles(sampleID)

	def setGenome(self, genome_version):
		self.genome_version = genome_version
		try:
			self.g_dt
		except AttributeError:
			self.readGenomeDatabase()
		self._createGenomeFiles()

	def createSampleFiles(self, sampleID):
		try:
			self.sample_dt
		except AttributeError:
			self.readSampleDatabase()

		self.sampleID = sampleID

		self.localRawBamFiles = [self.localReadsDir + x for x in self.reads_dt[self.reads_dt.SampleID == sampleID].FileLocations.to_list()]

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

		#os.makedirs(self.localSampleBamDir, exist_ok = True)
		#os.makedirs(self.localSampleTempDir, exist_ok = True)	

	def readGenomeDatabase(self):
		spreadsheet = self.gc.open_by_key(self.s_ID) # Or use open('Spreadsheet Name')
		worksheet = spreadsheet.worksheet('GenomeDatabase') # Access a specific sheet tab
		self.g_dt = get_as_dataframe(worksheet, evaluate_formulas=True)
	
	def setGenomeDatabase(self):
		for i in range(3):
			try:
				spreadsheet = self.gc.open_by_key(self.s_ID) # Or use open('Spreadsheet Name')
				set_with_dataframe(spreadsheet.worksheet('GenomeDatabase'), self.g_dt) # df is your DataFrame
				return True
			except Exception as e:
				print('Gspread exception: ' + e)
		return False

	def readSampleDatabase(self):
		spreadsheet = self.gc.open_by_key(self.s_ID) # Or use open('Spreadsheet Name')
		worksheet = spreadsheet.worksheet('SampleDatabase') # Access a specific sheet tab
		self.sample_dt = get_as_dataframe(worksheet, evaluate_formulas=True)
		worksheet = spreadsheet.worksheet('DNAReadsDatabase') # Access a specific sheet tab
		self.reads_dt = get_as_dataframe(worksheet, evaluate_formulas=True)
		self.merged_dt = pd.merge(self.sample_dt,self.reads_dt, on = 'SampleID')

	def addDNAReadRow(self, dna_dict, sample_dict):
		self.readSampleDatabase()
		assert set(dna_dict.keys()) == set(['SampleID','ProjectID','RunID','ReadLength','TotalBases','Instrument','LibraryID','LibraryLayout','LibrarySource','Platform','FileLocations','FileSize'])
		assert dna_dict['RunID'] not in self.reads_dt.RunID.to_list()
		self.reads_dt = pd.concat([self.reads_dt,pd.DataFrame([dna_dict])], ignore_index = True)
		assert set(sample_dict.keys()) == set(['SampleID','Sex','Species','DoB','BroodID','Parents','Ecogroup','LabReared','Inversion10'])
		if sample_dict['SampleID'] in self.sample_dt.SampleID.to_list():
			print('Warning: ' + sample_dict['SampleID'] + ' already in SampleDatabase')
		else:
			self.sample_dt = pd.concat([self.sample_dt,pd.DataFrame([sample_dict])], ignore_index = True)	
		self.setSampleDatabase()
	
	def addSampleRow(self, row_dict):
		self.readSampleDatabase()
		assert set(row_dict.keys()) == set(['SampleID','Sex','Species','DoB','BroodID','Parents','Ecogroup','LabReared','Inversion10'])
		assert row_dict['SampleID'] not in self.sample_dt.SampleID.to_list()
		self.sample_dt = pd.concat([self.sample_dt,pd.DataFrame([row_dict])], ignore_index = True)
		self.setSampleDatabase()
	
	def setSampleDatabase(self):
		for i in range(3):
			try:
				spreadsheet = self.gc.open_by_key(self.s_ID) # Or use open('Spreadsheet Name')
				set_with_dataframe(spreadsheet.worksheet('SampleDatabase'), self.sample_dt) # df is your DataFrame
				set_with_dataframe(spreadsheet.worksheet('DNAReadsDatabase'), self.reads_dt) # df is your DataFrame
				return True
			except Exception as e:
				print('Gspread exception: ' + str(e))
		return False

	def readAlignmentDatabase(self):
		assert self.genome_version

		spreadsheet = self.gc.open_by_key(self.a_ID) # Or use open('Spreadsheet Name')
		worksheet = spreadsheet.worksheet(self.genome_version) # Access a specific sheet tab
		self.alignment_dt = get_as_dataframe(worksheet, evaluate_formulas=True)

	def addAlignmentRow(self, row_dict):
		assert set(row_dict.keys()) == set(['SampleID','GenomeVersion','RunIDs','Coverage','TotalReads','UnmappedReads','DiscordantReads','InversionReads','DuplicationReads','ClippedReads','ChimericReads','minimap2_version','gatk_version','pysam_version','BamSize'])
		self.readAlignmentDatabase()
		assert row_dict['SampleID'] not in self.alignment_dt.SampleID.to_list()
		self.alignment_dt = pd.concat([self.alignment_dt,pd.DataFrame([row_dict])], ignore_index = True)
		self.setAlignmentDatabase()
		
	def setAlignmentDatabase(self):
		assert self.genome_version

		for i in range(3):
			try:
				spreadsheet = self.gc.open_by_key(self.a_ID) # Or use open('Spreadsheet Name')
				set_with_dataframe(spreadsheet.worksheet(self.genome_version), self.alignment_dt) # df is your DataFrame
				return True
			except Exception as e:
				print('Gspread exception: ' + str(e))
		return False


	def setSamples(self, projectIDs, sampleIDs, species, ecogroups, subgroups, rerun):
		assert self.genome_version

		try:
			self.sample_dt
		except AttributeError:
			self.readSampleDatabase()

		try:
			self.alignment_dt
		except AttributeError:
			self.readAlignmentDatabase()

		temp_dt = self.merged_dt
		if projectIDs is not None:
			temp_dt = temp_dt[temp_dt.ProjectID.isin(projectIDs)]
		if sampleIDs is not None:
			temp_dt = temp_dt[temp_dt.SampleID.isin(sampleIDs)]
		if species is not None:
			temp_dt = temp_dt[temp_dt.Species.isin(species)]
		if ecogroups is not None:
			temp_dt = temp_dt[temp_dt.Ecogroup.isin(ecogroups)]
		if subgroups is not None:
			temp_dt = temp_dt[temp_dt.Subgroup.isin(subgroups)]

		# Filter alignment database for requested genome version
		a_dt = self.alignment_dt[(self.alignment_dt.GenomeVersion == self.genome_version)]

		# Identify already run samples
		filter_set = set(a_dt.SampleID)
		already_run_samples = [x for x in set(temp_dt.SampleID) if x in filter_set]
		samples = [x for x in set(temp_dt.SampleID) if x not in filter_set]
		
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

	def removeSamples(self, bad_samples):
		self.samples = [x for x in self.samples if x not in bad_samples]

	def returnOptions(self, datatype):
		if datatype == 'Genomes':
			try:
				return self.g_dt.GenomeID.to_list()
			except AttributeError:
				self.readGenomeDatabase()
				return self.g_dt.GenomeID.to_list()
		if datatype == 'Samples':
			try:
				return self.reads_dt.SampleID.unique().tolist()
			except AttributeError:
				self.readSampleDatabase()
				return self.reads_dt.SampleID.unique().tolist()
		if datatype == 'Species':
			try:
				return self.merged_dt.Species.unique().tolist()
			except AttributeError:
				self.readSampleDatabase()
				return self.merged_dt.Species.unique().tolist()
		if datatype == 'ProjectIDs':
			try:
				return self.merged_dt.ProjectID.unique().tolist()
			except AttributeError:
				self.readSampleDatabase()
				return self.merged_dt.ProjectID.unique().tolist()
		if datatype == 'Ecogroups':
			try:
				return self.merged_dt.Ecogroup.dropna().unique().tolist()
			except AttributeError:
				self.readSampleDatabase()
				return self.merged_dt.Ecogroup.dropna().unique().tolist()
		if datatype == 'Subgroups':
			try:
				return self.merged_dt.Subgroup.dropna().unique().tolist()
			except AttributeError:
				self.readSampleDatabase()
				return self.merged_dt.Subgroup.dropna().unique().tolist()

	def _authenticateGS(self):
		self.downloadData(self.localCredentialFile)
		scopes = ['https://www.googleapis.com/auth/spreadsheets', 'https://www.googleapis.com/auth/drive']
		credentials = Credentials.from_service_account_file(self.localCredentialFile, scopes=scopes)
		self.gc = gspread.authorize(credentials)

	def _createMasterDirs(self):

		self.localPolymorphismsDir = self.localMasterDir + 'Polymorphisms/'	
		self.localReadsDir = self.localMasterDir + 'Reads/'		
		self.localSeqCoreDataDir = self.localMasterDir + 'SeqCoreData/'
		self.localBamfilesDir = self.localMasterDir + 'Bamfiles/'
		self.localGenomesDir = self.localMasterDir + 'Genomes/'
		self.localTempDir = self.localMasterDir + 'Temp/'
		self.localReadDownloadDir = self.localReadsDir + 'ReadDownloadFiles/'
		
		self.localCredentialFile = self.localMasterDir + 'cichlidsrsequencing_api_creds.json'
		self.localProcessesFile = self.localTempDir + 'ProcessInfo.csv'
		self.localErrorsDir = self.localMasterDir + 'Errors/'
		os.makedirs(self.localErrorsDir, exist_ok = True)



	def _createGenomeFiles(self):
		self.localBamRefDir = self.localBamfilesDir + self.genome_version + '/'
		self.localGenomeDir = self.localGenomesDir + self.genome_version + '/'
		if self.genome_version == 'Mzebra_UMD2a':
			self.localGenomeFile = self.localGenomeDir + 'GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
		elif self.genome_version == 'Mzebra_GT3':
			self.localGenomeFile = self.localGenomeDir + 'Mzebra_GT3.fasta'
		elif self.genome_version == 'Mzebra_GT3_NCBI':
			self.localGenomeFile = self.localGenomeDir + 'GCF_041146795.1_Mzebra_GT3a_genomic.fna'
			self.localMinimapGenomeFile = self.localGenomeDir + 'GCF_041146795.1_Mzebra_GT3a_genomic.mmi'
		elif self.genome_version == 'Mconophoros_GT1':
			self.localGenomeFile = self.localGenomeDir + 'anchored_kocher_E_Mchenga_conof_Male_contigs_hs_with_kocher_MC_female_molecules_mito_corrected.fasta'
		else:
			raise FileNotFoundError(self.genome_version + ' not an option')


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

	def uploadData(self, local_data, tarred = False, parallel = False):

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
			if parallel:
				command = ['rclone', 'copy', local_path + relative_name, cloud_path + relative_name]
				return command
			else:
				output = subprocess.run(['rclone', 'copy', local_path + relative_name, cloud_path + relative_name], capture_output = True, encoding = 'utf-8')
			#subprocess.run(['rclone', 'check', local_path + relative_name, cloud_path + relative_name], check = True)

		elif os.path.isfile(local_path + relative_name):
			#print(['rclone', 'copy', local_path + relative_name, cloud_path])
			if parallel:
				command = ['rclone', 'copy', local_path + relative_name, cloud_path]
				return command
			else:
				output = subprocess.run(['rclone', 'copy', local_path + relative_name, cloud_path], capture_output = True, encoding = 'utf-8')
				output = subprocess.run(['rclone', 'check', local_path + relative_name, cloud_path], check = True, capture_output = True, encoding = 'utf-8')
		else:
			raise Exception(local_data + ' does not exist for upload')

		if not parallel:
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







