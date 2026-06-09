import pdb, vcf
from helper_modules.file_manager import FileManager as FM

class GenotypeGroup:
	def __init__(self, name, samples):
		self.name = name
		self.samples = samples

	def addSample(self, vcf_file):
		pass

fm_obj = FM(genome_version = 'Mzebra_GT3_NCBI')
fm_obj.readSampleDatabase()
fm_obj.readAlignmentDatabase()

all_samples = fm_obj.sample_dt[fm_obj.sample_dt.Ecogroup.isin(['Deep_Benthic','Shallow_Benthic'])].SampleID.to_list()
aligned_samples = fm_obj.alignment_dt[fm_obj.alignment_dt.SampleID.isin(all_samples)].SampleID.to_list()

for sampleID in aligned_samples:
	out_vcf = fm_obj.localNikeshDir + 'QTG_Candidates/' + sampleID + '_candidate_QTNs.vcf.gz'
	fm_obj.downloadData(out_vcf)
	vcf_obj = vcf.VCFReader(out_vcf)
	for record in vcf_obj:
		pdb.set_trace()
