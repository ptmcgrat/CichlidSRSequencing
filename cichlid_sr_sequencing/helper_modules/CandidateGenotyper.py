from helper_modules.smallVariantGenotyper import normalize_sites
from helper_modules.file_manager import FileManager as FM
import pandas as pd
import os,subprocess,pdb
from types import SimpleNamespace

class CandidateGenotyper:
	def __init__(self, genome_version, candidate_QTNs_ID, ecogroups):
		self.genome_version = genome_version
		self.fm_obj = FM(genome_version)
		self.QTNs_ID = candidate_QTNs_ID
		self.ecogroups = ecogroups
		
		# Set samples
		self.fm_obj.readSampleDatabase()
		self.fm_obj.readAlignmentDatabase()
		all_samples = self.fm_obj.sample_dt[self.fm_obj.sample_dt.Ecogroup.isin(self.ecogroups)].SampleID.to_list()
		self.samples = self.fm_obj.alignment_dt[self.fm_obj.alignment_dt.SampleID.isin(all_samples)].SampleID.to_list()

		self._createMasterFiles()
		self._downloadData()	
		self._createMasterVCFs()

	def _createMasterFiles(self):
		self.masterDir = self.fm_obj.localCandidateQTNDir + self.QTNs_ID + '/'
		self.masterSampleVCFDir = self.masterDir + 'SampleVCFs/'
		self.masterTSV = self.masterDir + self.QTNs_ID + '.tsv'
		self.masterSV_VCF = self.masterDir + self.QTNs_ID + '_sv.vcf'
		self.masterSV_Norm_VCF = self.masterDir + self.QTNs_ID + '_sv.norm.vcf.gz'
		self.masterLV_Norm_VCF = self.masterDir + self.QTNs_ID + '_lv.csv'

	def _downloadData(self):
		os.makedirs(self.fm_obj.localErrorsDir, exist_ok = True)
		self.fm_obj.downloadData(self.fm_obj.localGenomeFile)
		self.fm_obj.downloadData(self.masterDir)

	def _createMasterVCFs(self, threshold = 10):
		# Master tsv file of QTN candidates
		dt = pd.read_csv(self.masterTSV, sep = '\t')

		# Create normalized small SV vcf
		small_dt = dt[(dt.Alt.str.len() <= threshold) & (dt.Reference.str.len() <= threshold)]
		with open(self.masterSV_VCF,'w') as fp:
			print('##fileformat=VCFv4.2', file = fp)
			print('##source=MSAUniqueVariantCaller', file = fp)
			print('##sample=Y_reg', file = fp)
			print('##contig=<ID=NC_135176.1>', file = fp)
			print('##INFO=<ID=TYPE,Number=1,Type=String,Description="SNV, INS, or DEL">', file = fp)
			print('##INFO=<ID=ALN_COL,Number=1,Type=Integer,Description="1-based alignment column where the variant starts">', file = fp)
			print('##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length of inserted or deleted sequence">', file = fp)
			print('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']), file = fp)
			
			for i,row in small_dt.iterrows():
				print('\t'.join([row.Chromosome,str(row.Position),row.Name,row.Reference.upper(), row.Alt.upper(), str(row.Q), 'PASS', row.Info]), file = fp)

		normalize_sites(self.masterSV_VCF, self.fm_obj.localGenomeFile, self.masterSV_Norm_VCF)

		# Create large LV csv
		dt[(dt.Alt.str.len() > threshold) | (dt.Reference.str.len() > threshold)].to_csv(self.masterLV_Norm_VCF)

	def genotypeSamples(self, num_parallel = 48):

		commands = []
		bad_samples = []
		for sampleID in self.samples[0:2]:
			out_vcf = self.masterSampleVCFDir + sampleID + '_' + self.QTNs_ID + '.vcf.gz'
			command = ['python','-m', 'unit_scripts.genotypeCandidates', self.masterSV_Norm_VCF, self.masterLV_Norm_VCF, out_vcf, self.genome_version, sampleID]
		error_file = self.fm_obj.localErrorsDir + 'QTGFinder_' + sampleID + '_errors.txt'
		commands.append(SimpleNamespace(sampleID=sampleID, command = command, error_file = error_file))
	
		for i,data in enumerate(commands):
	
			if i < num_parallel:
				print(data.sampleID)
				data.error_fp = open(data.error_file, 'w')
				data.process = subprocess.Popen(data.command, stderr = data.error_fp, stdout = subprocess.DEVNULL)
			else:
				data.process = None
	
		while commands:
			finished_processes = [x for x in commands if x.process is not None and x.process.poll() is not None]
			for data in finished_processes:
				data.error_fp.close()
				print(f'..{data.sampleID} complete..', end = '')
				if data.process.returncode != 0:
					bad_samples.append(data.sampleID)
				else:
					subprocess.run(['rm',data.error_file])
				commands.remove(data)  # Remove finished process from monitoring list
				next_command = next((x for x in commands if x.process is None),None)
				if next_command is not None:
					next_command.error_fp = open(next_command.error_file, 'w')
					next_command.process = subprocess.Popen(next_command.command, stderr = next_command.error_fp, stdout = subprocess.DEVNULL)

		print(bad_samples)