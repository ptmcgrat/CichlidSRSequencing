import pysam, pdb, math

from collections import defaultdict

from helper_modules.file_manager import FileManager as FM

insertions = [24869679,24918628,24921178,24921394]

fm_obj = FM('Mzebra_GT3','YH_005_m')
print(fm_obj.localSampleBamDir)

fm_obj = FM(genome_version = 'Mzebra_GT3')
samples = fm_obj.a_dt[fm_obj.a_dt.GenomeVersion == 'Mzebra_GT3'].SampleID.to_list()

for sample in samples:
	try:
		organism = fm_obj.sample_dt[fm_obj.sample_dt.SampleID == sample].Species.values[0]
	except IndexError:
		organism = 'UnknownSpecies'

	out_dt = pd.DataFrame(columns = ['InsertionLocation','Organism','Sample','AlnChrom','AlnStart','AlnEnd','PairChrom','PairStart','TemplateLength','Duplicate'])

	fm_obj.createSampleFiles(sample)
	fm_obj.downloadData(fm_obj.localDiscordantBamFile)
	pysam.index(fm_obj.localDiscordantBamFile)
	bam_obj = pysam.AlignmentFile(fm_obj.localDiscordantBamFile)

	for insertion in insertions:
		for aln in bam_obj.fetch('NC_036789.1',insertion - 2000,insertion+2000):
			if not aln.is_proper_pair:
				out = [str(insertion), organism, sample, aln.reference_name, aln.reference_start, aln.reference_end, aln.next_reference_name, aln.next_reference_start, aln.template_length, aln.is_duplicate]
				out_dt.loc[len(out_dt)] = out


out_dt.to_csv('InsertionReadsInAllSamples.csv')