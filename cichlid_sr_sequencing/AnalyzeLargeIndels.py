import pysam, pdb
import pandas as pd
from helper_modules.file_manager import FileManager as FM
fm_obj = FM(genome_version = 'Mzebra_GT3')

#aln_obj = pysam.AlignmentFile('/Users/pmcgrath7/Temp/CichlidSequencingData/Bamfiles/Mzebra_GT3/YH_005_m/YH_005_m.all.bam')

insertions = [24869679,24918628,24921178,24921394]
male_yhs = ['YH_011_m','YH_008_m','YH_1_m','YH_005_m','YH_003_m','YH_025','YH_026','YH_027','YH_028','YH_029','YH_030','YH_037','YH_038','YH_039','YH_040','YH_041','YH_042','YH_046','YH_047','YH_048','YH_049','YH_050','YH_051','YH_052','YH_053','YH_054','YH_068','YH_069','YH_070','YH_012_m','YH_014_m','Kocher_YH1m','YH_015_m']

out_dt = pd.DataFrame(columns = ['InsertionLocation','Sample','AlnChrom','AlnStart','AlnEnd','PairChrom','PairStart','TemplateLength','Duplicate'])

for sample in male_yhs:
	fm_obj.createSampleFiles(sample)
	fm_obj.downloadData(fm_obj.localDiscordantBamFile)
	pysam.index(fm_obj.localDiscordantBamFile)
	bam_obj = pysam.AlignmentFile(fm_obj.localDiscordantBamFile)

	for insertion in insertions:
		for aln in bam_obj.fetch('NC_036789.1',insertion - 2000,insertion+2000):
			if not aln.is_proper_pair:
				out = [str(insertion), sample, aln.reference_name, aln.reference_start, aln.reference_end, aln.next_reference_name, aln.next_reference_start, aln.template_length, aln.is_duplicate]
				out_dt.loc[len(out_dt)] = out

out_dt.to_csv('InsertionReadsInYHMales.csv')