import os, subprocess, pdb, random, datetime, platform
# Google Sheets access (SampleData + Alignment Database). This file only READS the sheets: set_with_dataframe is
# deliberately not imported, so nothing here can overwrite the shared lab databases.
# Server env needs: pip install gspread gspread-dataframe google-auth
import gspread
from gspread_dataframe import get_as_dataframe
from google.oauth2.service_account import Credentials
import pandas as pd
import numpy as np
from xml.etree import ElementTree as ET
from multiprocessing import cpu_count
#from pysam import VariantFile
from collections import defaultdict

class FileManager():
	def __init__(self, genome_version = '', rcloneRemote = 'ptm_dropbox:/', masterDir = 'COS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/'): # same Dropbox path as the PI's file_manager

		self.genome_version = genome_version
		self.s_ID = '1NmgB_TWoO01Qz2ufvECuZFkxXayhUsyu8wQGStVB_8k' # SampleData sheet
		self.a_ID = '1vlA_eJ09RSKCaM0foFMc0n5tL95XK4deB3dr7kAFx3E' # Alignment Database sheet
		# genome_version must match an Alignment Database TAB name (e.g. 'Mzebra_GT3_NCBI', 'Mconophoros_GT1_NCBI'),
		# because readAlignmentDatabase() opens the tab with this name.

		# Local master directory: same rule as the PI's file_manager.py (lines 13-16), so everything is written where the lab's
		# other pipelines expect it.
		#   lab servers (utaka, yellowhead, burtoni): /Data/<USER>/Temp/CichlidSequencingData/
		#       e.g. Outputs -> /Data/student/Temp/CichlidSequencingData/Outputs/
		#   any other machine (e.g. a laptop):        ~/Temp/CichlidSequencingData/
		node = platform.node().lower()
		if 'utaka' in node or 'yellowhead' in node or 'burtoni' in node:
			self.localMasterDir = '/Data/' + os.getenv('USER') + '/Temp/CichlidSequencingData/'
		else:
			self.localMasterDir = os.getenv('HOME').rstrip('/') + '/' + 'Temp/CichlidSequencingData/'
		# In Oct 2023, Curtis & I mounted a new RAID array to Utaka server that I then mounted into /Output (used by --Output in callVariants.py)
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
		# Google Sheets sign-in is NOT done here (the PI's file_manager does it at this point). Instead it happens the first time
		# a read*Database() method needs the sheets (see _sheetsClient), so scripts/runs that never touch the sheets
		# (e.g. --local_test, uploads) don't need the credential file or network access.

	def _createMasterDirs(self):
		# e.g. on a lab server: localGenomesDir = /Data/<USER>/Temp/CichlidSequencingData/Genomes/
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
		elif self.genome_version == 'kocher_E_Mchenga_conof_Male':
			self.localGenomeFile = self.localGenomeDir + 'anchored_kocher_E_Mchenga_conof_Male_contigs_hs_with_kocher_MC_female_molecules_mito_corrected.fasta'
		elif self.genome_version == 'Mconophoros_GT1':
			self.localGenomeFile = self.localGenomeDir + 'anchored_kocher_E_Mchenga_conof_Male_contigs_hs_with_kocher_MC_female_molecules_mito_corrected.fasta'
		elif self.genome_version == 'Mzebra_GT3_NCBI':
			self.localGenomeFile = self.localGenomeDir + 'GCF_041146795.1_Mzebra_GT3a_genomic.fna'
		elif self.genome_version == 'Mconophoros_GT1_NCBI':
			self.localGenomeFile = self.localGenomeDir + 'GCF_057382185.1_M_conophoros_male_GT1.1_genomic.fna'

		# Service-account credentials for the Google Sheets (same file and same local location as the PI's file_manager uses);
		# downloadData() maps it to <cloudMasterDir>/cichlidsrsequencing_api_creds.json.
		self.localCredentialFile = self.localMasterDir + 'cichlidsrsequencing_api_creds.json'
		self.localReadDownloadDir = self.localReadsDir + 'ReadDownloadFiles/'
		# Roots shared by all analyses. Each analysis gets its own <genome>/<analysisID>/ subfolders under these - see createAnalysisIDFiles().
		self.localDatabasesDir = self.localMasterDir + 'Databases/'
		self.localOutputDir = self.localMasterDir + 'Outputs/'
		self.localPCADir = self.localOutputDir + 'pca_outputs'
		self.localBEDDir = self.localOutputDir + 'BEDFiles/'

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
		self.localTestGVCFIndex = self.localSampleTestingDir + sampleID + '_0.0001_subset.g.vcf.gz.tbi'

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
		# Defines (does NOT create) the folders/files for one analysis on self.genome_version. callVariants.py builds its
		# analysis folders from these paths, and any later script (PCA, filtering, ...) can find an analysis's outputs again
		# from just the genome version + analysis ID, e.g.:
		#   fm = FileManager('Mconophoros_GT1_NCBI'); fm.createAnalysisIDFiles('CaitVCF'); print(fm.localAnalysisVCFFile)
		# Layout:
		#   Outputs/<genome>/<analysisID>/                       frozen sample list, sample map, command log
		#   Outputs/<genome>/<analysisID>/interval_vcfs/         per-interval VCFs from GenotypeGVCFs
		#   Outputs/<genome>/<analysisID>/vcf_concat_output/     final VCF: <analysisID>.vcf.gz
		#   Databases/<genome>/<analysisID>/                     GenomicsDB workspaces
		self.analysisID = analysisID
		self.localAnalysisFile = self.localAnalysisDir + analysisID + '.csv'
		self.localAnalysisOutputDir = self.localOutputDir + self.genome_version + '/' + analysisID + '/'
		self.localAnalysisIntervalVCFDir = self.localAnalysisOutputDir + 'interval_vcfs/'
		self.localAnalysisConcatDir = self.localAnalysisOutputDir + 'vcf_concat_output/'
		self.localAnalysisDatabasesDir = self.localDatabasesDir + self.genome_version + '/' + analysisID + '/'
		self.localSampleManifestFile = self.localAnalysisOutputDir + analysisID + '_samples.txt' # frozen sample list
		self.localSampleMapFile = self.localAnalysisOutputDir + 'sample_map.txt'
		self.localCompressionListFile = self.localAnalysisOutputDir + 'compression_file_list.txt'
		self.localCommandsLogFile = self.localAnalysisOutputDir + 'commands_log.txt'
		self.localMissingGVCFsFile = self.localAnalysisOutputDir + 'missing_gvcfs.txt'
		self.localAnalysisVCFFile = self.localAnalysisConcatDir + analysisID + '.vcf.gz' # a re-run of -c writes <analysisID>_2.vcf.gz etc. next to it

	def analysisIDExists(self, analysisID):
		# True if either folder of this analysis already exists for self.genome_version
		return os.path.exists(self.localOutputDir + self.genome_version + '/' + analysisID) or os.path.exists(self.localDatabasesDir + self.genome_version + '/' + analysisID)

	def returnAnalysisIDs(self):
		# Analyses that already have an output folder for self.genome_version
		genome_output_dir = self.localOutputDir + self.genome_version + '/'
		return sorted(os.listdir(genome_output_dir)) if os.path.isdir(genome_output_dir) else []

	def returnGenomeVersions(self):
		return self.returnCloudDirs(self.localGenomesDir)

	# -----------------------------------------------------------------------------------------------------------------
	# Google Sheets - READ ONLY. Adapted from file_manager_PI.py: _authenticateGS (PI lines 243-247),
	# readSampleDatabase (PI 90-96), readAlignmentDatabase (PI 128-133), setSamples (PI 155-200 -> selectSamples).
	# Differences from the PI's versions:
	#   - sign-in happens lazily (_sheetsClient) instead of in __init__
	#   - read-only Google scopes, and no write methods (addSampleRow, setSampleDatabase, ... are intentionally not copied)
	#   - empty rows and rows without a SampleID are dropped (get_as_dataframe returns hundreds of blank rows, and the
	#     Mzebra_GT3_NCBI alignment tab has summary rows at the bottom with no SampleID)
	# -----------------------------------------------------------------------------------------------------------------
	def _authenticateGS(self):
		self.downloadData(self.localCredentialFile)
		# Read-only scopes: even a bug in this file can't modify the shared sheets. (The PI's file_manager uses the read/write
		# scopes 'https://www.googleapis.com/auth/spreadsheets' and '.../auth/drive'.)
		scopes = ['https://www.googleapis.com/auth/spreadsheets.readonly', 'https://www.googleapis.com/auth/drive.readonly']
		credentials = Credentials.from_service_account_file(self.localCredentialFile, scopes = scopes)
		self.gc = gspread.authorize(credentials) # gc = "gspread client": the signed-in connection to Google Sheets

	def _sheetsClient(self):
		# Lazy sign-in: connect the first time the sheets are needed, then reuse the same connection
		if not hasattr(self, 'gc'):
			self._authenticateGS()
		return self.gc

	def _readWorksheet(self, sheet_key, tab_name):
		# Read one tab into a DataFrame, dropping blank rows, unnamed blank columns and rows without a SampleID.
		# Named columns are kept even if empty (e.g. a new AnalysisSets column that hasn't been filled in yet).
		worksheet = self._sheetsClient().open_by_key(sheet_key).worksheet(tab_name)
		df = get_as_dataframe(worksheet, evaluate_formulas = True)
		df = df.dropna(how = 'all')
		df = df.loc[:, [not str(column).startswith('Unnamed') for column in df.columns]]
		df.columns = [str(column).strip() for column in df.columns]
		if 'SampleID' in df.columns:
			df = df.dropna(subset = ['SampleID'])
			df['SampleID'] = df['SampleID'].astype(str).str.strip()
		return df.reset_index(drop = True)

	def __getstate__(self):
		# Python calls this when the object is pickled (copied byte-for-byte), which multiprocessing.Pool does for every task it
		# sends to a worker process: callVariants passes bound methods like self.RunGenomicsDBImport, and those carry this
		# FileManager along. The Google Sheets connection (self.gc) holds network sessions/credentials that may not pickle, and the
		# workers never use it, so it's left out of the copy. If a worker ever did read the sheets, _sheetsClient would reconnect.
		state = self.__dict__.copy()
		state.pop('gc', None)
		return state

	def readSampleDatabase(self):
		self.sample_dt = self._readWorksheet(self.s_ID, 'SampleDatabase')
		self.reads_dt = self._readWorksheet(self.s_ID, 'DNAReadsDatabase')
		self.merged_dt = pd.merge(self.sample_dt, self.reads_dt, on = 'SampleID') # one row per sequencing RUN, so de-duplicate SampleIDs after filtering

	def readAlignmentDatabase(self):
		assert self.genome_version, 'genome_version (= Alignment Database tab name) must be set before reading the Alignment Database'
		self.alignment_dt = self._readWorksheet(self.a_ID, self.genome_version)

	def readAnalysisSets(self):
		# AnalysisSets tab of the SampleData sheet: one row per SampleID, one column per analysis (e.g. PM30, CaitVCF, CORE_PCA).
		self.analysis_sets_dt = self._readWorksheet(self.s_ID, 'AnalysisSets')

	def returnAnalysisSetNames(self):
		# Analysis columns = AnalysisSets columns that aren't SampleDatabase metadata columns
		if not hasattr(self, 'analysis_sets_dt'):
			self.readAnalysisSets()
		if not hasattr(self, 'sample_dt'):
			self.readSampleDatabase()
		return [column for column in self.analysis_sets_dt.columns if column not in self.sample_dt.columns]

	def returnAnalysisSetSamples(self, analysis_set):
		# SampleIDs marked for one analysis column. 'Yes' marks a sample (case-insensitive; 'Y' and TRUE from a checkbox also count).
		# Blank / No / FALSE = not in the analysis. Any other value is ignored but printed, so a typo can't silently drop a sample.
		if not hasattr(self, 'analysis_sets_dt'):
			self.readAnalysisSets()
		dt = self.analysis_sets_dt
		if analysis_set not in dt.columns:
			raise Exception(f"'{analysis_set}' is not a column in the AnalysisSets tab. Analysis columns found: {self.returnAnalysisSetNames()}")
		values = dt[analysis_set].fillna('').astype(str).str.strip().str.upper()
		selected = values.isin(['YES', 'Y', 'TRUE'])
		unrecognized = values[~selected & ~values.isin(['', 'NO', 'N', 'FALSE'])]
		if len(unrecognized) > 0:
			print(f"WARNING: AnalysisSets column {analysis_set} has values that aren't Yes/blank; these rows were NOT selected:")
			print('  ' + ', '.join(f"{sample}='{dt.loc[i, analysis_set]}'" for i, sample in dt.loc[unrecognized.index, 'SampleID'].items()))
		duplicates = sorted(set(dt.SampleID[dt.SampleID.duplicated()]))
		if duplicates:
			print(f"WARNING: SampleIDs listed more than once in the AnalysisSets tab: {','.join(duplicates)}")
		return set(dt.loc[selected, 'SampleID'])

	def selectSamples(self, analysis_set = None, projectIDs = None, species = None, ecogroups = None, subgroups = None,
					  extra_samples = None, exclude_samples = None, whole_alignment_tab = False):
		# Returns the sorted list of SampleIDs for a variant-calling analysis on self.genome_version (also stored in self.samples).
		# Adapted from the PI's setSamples, with one important difference: setSamples is for the ALIGNMENT pipeline, so it keeps
		# samples that are NOT yet in the alignment tab. Variant calling needs the opposite - only samples that ARE in the alignment
		# tab for this genome have BAMs/GVCFs.
		# How the selectors combine:
		#   1. union of: the AnalysisSets column (analysis_set)
		#              | the metadata filter (projectIDs / species / ecogroups / subgroups; combined with AND, like setSamples)
		#              | every sample in this genome's alignment tab (whole_alignment_tab)
		#              | extra_samples
		#   2. minus exclude_samples
		#   3. minus anything not aligned to self.genome_version (printed, never silently dropped)
		assert self.genome_version
		self.readSampleDatabase()
		self.readAlignmentDatabase()
		aligned = set(self.alignment_dt.SampleID)
		known = set(self.sample_dt.SampleID)
		candidates = set()
		summary = []

		if analysis_set:
			self.readAnalysisSets()
			set_samples = self.returnAnalysisSetSamples(analysis_set)
			not_in_sample_db = sorted(set_samples - known)
			if not_in_sample_db:
				print(f"WARNING: {len(not_in_sample_db)} SampleIDs marked in AnalysisSets column {analysis_set} are not in the SampleDatabase tab (renamed or removed?): {','.join(not_in_sample_db)}")
			candidates |= set_samples
			summary.append(f"AnalysisSets column {analysis_set}: {len(set_samples)}")

		if projectIDs or species or ecogroups or subgroups:
			temp_dt = self.merged_dt
			if projectIDs:
				temp_dt = temp_dt[temp_dt.ProjectID.isin(projectIDs)]
			if species:
				temp_dt = temp_dt[temp_dt.Species.isin(species)]
			if ecogroups:
				temp_dt = temp_dt[temp_dt.Ecogroup.isin(ecogroups)]
			if subgroups:
				temp_dt = temp_dt[temp_dt.Subgroup.isin(subgroups)]
			metadata_samples = set(temp_dt.SampleID)
			candidates |= metadata_samples
			summary.append(f"metadata filter (ProjectID={projectIDs}, Species={species}, Ecogroup={ecogroups}, Subgroup={subgroups}): {len(metadata_samples)}")

		if whole_alignment_tab:
			candidates |= aligned
			summary.append(f"all samples in Alignment Database tab {self.genome_version}: {len(aligned)}")

		if extra_samples:
			typos = sorted(set(extra_samples) - known - aligned)
			if typos:
				raise Exception(f"extra_samples not found in the SampleDatabase tab or the {self.genome_version} alignment tab (typo?): {','.join(typos)}")
			candidates |= set(extra_samples)
			summary.append(f"extra_samples: {len(set(extra_samples))}")

		if exclude_samples:
			not_selected_anyway = sorted(set(exclude_samples) - candidates)
			if not_selected_anyway:
				print(f"Note: these exclude_samples weren't in the selection anyway: {','.join(not_selected_anyway)}")
			candidates -= set(exclude_samples)
			summary.append(f"excluded: {len(set(exclude_samples) - set(not_selected_anyway))}")

		not_aligned = sorted(candidates - aligned)
		if not_aligned:
			print(f"WARNING: {len(not_aligned)} selected samples are NOT aligned to {self.genome_version} (no row in that Alignment Database tab) and are left out:")
			print('  ' + ','.join(not_aligned))
		self.samples = sorted(candidates & aligned)

		flagged = self.sample_dt[self.sample_dt.SampleID.isin(self.samples) & self.sample_dt.Ecogroup.isin(['Error', 'Exclude'])]
		if len(flagged) > 0:
			print(f"WARNING: selected samples with Ecogroup Error/Exclude in the SampleDatabase tab: {','.join(sorted(flagged.SampleID))}")

		print(f"Sample selection for {self.genome_version}: " + '; '.join(summary))
		print(f"{len(self.samples)} samples selected ({len(not_aligned)} dropped as not aligned to {self.genome_version})")
		return self.samples

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