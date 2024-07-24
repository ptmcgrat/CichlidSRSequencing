import os, subprocess, pdb, random, datetime, platform
import pandas as pd
import numpy as np
from xml.etree import ElementTree as ET
from multiprocessing import cpu_count
#from pysam import VariantFile
from collections import defaultdict

class FileManager():
	def __init__(self, genome_version = '', rcloneRemote = 'ptm_dropbox:/', masterDir = 'CoS/BioSci/McGrath/Apps/CichlidSequencingData/'):

		self.genome_version = genome_version
		
		if platform.node() == 'ebb-utaka.biosci.gatech.edu' or platform.node() == 'utaka.biosci.gatech.edu' or platform.node() == 'utaka':
			# basically ignore the below line because we'll work on mzebra and not have to deal with the messy setup of utaka server
			self.localMasterDir = '/Data/' + os.getenv('USER') + '/Data/CichlidSequencingData/'
			# In Oct 2023, Curtis & I mounted a new RAID array to Utaka server that I then mounted into /Output. The below line will allow the scripts to access data stroed in this directory. 
			self.localStorageDir = '/Output/'
		else:
			#Master directory for local data. This is where the data structure will be setup on the server. It will be '/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/'
			self.localMasterDir = os.getenv('HOME').rstrip('/') + '/' + 'Data/CichlidSequencingData/'
			self.localStorageDir = '/Output/'
		# Identify cloud directory for rclone
		self.rcloneRemote = rcloneRemote
		
		# On some computers, the first directory is McGrath, on others it's BioSci-McGrath. Use rclone to figure out which
		output = subprocess.run(['rclone', 'lsf', self.rcloneRemote + masterDir], capture_output = True, encoding = 'utf-8')
		if output.stderr == '':
			self.cloudMasterDir = self.rcloneRemote + masterDir
		else:
			masterDir = 'CoS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/'
			output = subprocess.run(['rclone', 'lsf', self.rcloneRemote + masterDir], capture_output = True, encoding = 'utf-8')
			if output.stderr == '':
				self.cloudMasterDir = self.rcloneRemote + masterDir
			else:
				raise Exception('Cant find master directory (' + masterDir + ') in rclone remote (' + rcloneRemote + '')

		self._createMasterDirs()

	def _createMasterDirs(self):
		# localGenomesDir = /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Genomes
		# Same logic applies to all of the below in creating the local directory structure.
		self.localGenomesDir = self.localMasterDir + 'Genomes/'
		self.localPolymorphismsDir = self.localMasterDir + 'Polymorphisms/'	
		self.localPileupDir = self.localMasterDir + '/Pileups/'	+ self.genome_version
		
		self.localReadsDir = self.localMasterDir + 'Reads/'
		self.localSeqCoreDataDir = self.localMasterDir + 'SeqCoreData/'
		self.localBamfilesDir = self.localMasterDir + 'Bamfiles/'
		self.localTempDir = self.localMasterDir + 'Temp/'
		self.localAnalysisDir = self.localMasterDir + 'Analysis/'

		self.localBamRefDir = self.localBamfilesDir + self.genome_version + '/'
		self.localGenomeDir = self.localGenomesDir + self.genome_version + '/'
		self.localPBGenomeFile = self.localGenomeDir + self.genome_version + '_v1.fa'
		self.localPBIndex = self.localGenomeDir + self.genome_version + '_v1.mmi'
		if self.genome_version == 'Mzebra_UMD2a':
			self.localGenomeFile = self.localGenomeDir + 'GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
		elif self.genome_version == 'Mzebra_GT2':
			self.localGenomeFile = self.localGenomeDir + 'Mzebra_GT2.fna'
		elif self.genome_version == 'Mzebra_GT3':
			self.localGenomeFile = self.localGenomeDir + 'Mzebra_GT3.fasta'
		

		self.localSampleFile = self.localReadsDir + 'SampleDatabase.csv'
		self.localSampleFile_v2 = self.localReadsDir + 'SampleDatabase_v2.xlsx'
		self.localSampleFile_for_grant = self.localReadsDir + 'SampleDatabase_v2_2024_feb_grant.xlsx'
		self.localAlignmentFile = self.localBamfilesDir + 'AlignmentDatabase.csv'
		self.localReadDownloadDir = self.localReadsDir + 'ReadDownloadFiles/'
		self.localDatabasesDir = self.localMasterDir + 'Databases/'
		self.localOutputDir = self.localMasterDir + 'Outputs/'
		self.localPCADir = self.localOutputDir + 'pca_outputs'

		# Below block is to map file structures in /Output on the Utaka server:
		self.StorageBamfilesDir = self.localStorageDir + 'Bamfiles/'
		self.StorageBamRefDir = self.StorageBamfilesDir + self.genome_version + '/'
		self.StorageOutputDir = self.localStorageDir + 'Outputs/'

	def createSampleFiles(self, sampleID):
		self.sampleID = sampleID
		self.localSampleBamDir = self.localBamRefDir + sampleID + '/'
		self.localSampleTestingDir = self.localBamfilesDir + 'nikesh_local_testing/' + sampleID + '/' 
		self.localBamFile = self.localSampleBamDir + sampleID + '.all.bam'
		self.localTestBamFile = self.localSampleTestingDir + sampleID + '_0.0001_subset.all.bam'
		self.localBamIndex = self.localSampleBamDir + sampleID + '.all.bai'
		self.localTestBamIndex = self.localSampleTestingDir + sampleID + '_0.0001_subset.all.bai'
		self.localUnmappedBamFile = self.localSampleBamDir + sampleID + '.unmapped.bam'
		self.localDiscordantBamFile = self.localSampleBamDir + sampleID + '.discordant.bam'
		self.localInversionBamFile = self.localSampleBamDir + sampleID + '.inversion.bam'
		self.localDuplicationBamFile = self.localSampleBamDir + sampleID + '.duplication.bam'
		self.localClippedBamFile = self.localSampleBamDir + sampleID + '.clipped.bam'
		self.localChimericBamFile = self.localSampleBamDir + sampleID + '.chimeric.bam'
		self.localGVCFFile = self.localSampleBamDir + sampleID + '.g.vcf.gz'
		self.localTestGVCFFile = self.localSampleTestingDir + sampleID + '_0.0001_subset.g.vcf.gz'
		self.localGVCFIndex = self.localSampleBamDir + sampleID + '.g.vcf.gz.tbi'
		self.localTestGVCFIndex = self.localSampleTestingDir + sampleID + '_0.01_subset.g.vcf.gz.tbi'

		# Below block is to access files in the /Output storage directory
		self.StorageSampleBamDir = self.StorageBamRefDir + sampleID + '/'
		self.StorageGVCFFile = self.StorageSampleBamDir + sampleID + '.g.vcf.gz'
		self.StorageGVCFIndex = self.StorageSampleBamDir + sampleID + '.g.vcf.gz.tbi'

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
		self.downloadData(self.localSampleFile)
		dt = pd.read_csv(self.localSampleFile)
		return set(dt.SampleID)

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
