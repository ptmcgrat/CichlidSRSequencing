import pysam, pdb, math

from collections import defaultdict

from helper_modules.file_manager import FileManager as FM
from helper_modules.ChimericData import ChimericRead

#LG10_inversion = ('NC_036789.1',11817927,11841666,29881282,29902558)
LG10_inversion = ('NC_036789.1',11817927,24677955,25067688,29902558)
female_yhs = ['Kocher_YH2f','YH_007_f','YH_006_f','YH_013_f','YH_016','YH_017','YH_018','YH_019','YH_020','YH_021','YH_022','YH_023','YH_024','YH_031','YH_032','YH_033','YH_055','YH_056','YH_057','YH_058','YH_059','YH_060','YH_061','YH_062','YH_063','YH_064','YH_065','YH_066','YH_010_f']
male_yhs = ['YH_011_m','YH_009_m','YH_008_m','YH_1_m','YH_005_m','YH_003_m','YH_004_m','YH_025','YH_026','YH_027','YH_028','YH_029','YH_030','YH_037','YH_038','YH_039','YH_040','YH_041','YH_042','YH_046','YH_047','YH_048','YH_049','YH_050','YH_051','YH_052','YH_053','YH_054','YH_068','YH_069','YH_070','YH_012_m','YH_014_m','Kocher_YH1m','YH_015_m']
mcs = ['MC-1G2G-f','MC-1O7O-f','MC-1R6R-f','MC-2B4B-f','MC-2G8G-f','MC-2P4P-f','MC-3P6P-f','MC-3R6R-f','MC-4O5O-f','MC-5B6B-f','MC-5G11G-f']

minMapQ = 30

fm_obj = FM(genome_version = 'Mzebra_GT3')
samples = fm_obj.a_dt[fm_obj.a_dt.GenomeVersion == 'Mzebra_GT3'].SampleID.to_list()

for sample in samples:
	try:
		organism = fm_obj.sample_dt[fm_obj.sample_dt.SampleID == sample].Species.values[0]
	except IndexError:
		organism = 'UnknownSpecies'
	fm_obj.createSampleFiles(sample)
	fm_obj.downloadData(fm_obj.localChimericBamFile)
	pysam.index(fm_obj.localChimericBamFile)
	discoveryChimeras = defaultdict(int)	
	bam_obj = pysam.AlignmentFile(fm_obj.localChimericBamFile)
	for read in bam_obj.fetch(LG10_inversion[0],LG10_inversion[2],LG10_inversion[3]):
		if not read.is_secondary and read.mapq > minMapQ:
			try:
				SA = read.get_tag('SA')
			except KeyError:
				continue
			if SA.split(',')[0]!=LG10_inversion[0]:
				continue
			read.set_tag('SA',SA.replace(SA.split(',')[0], '9'))

			try:
				newRead = ChimericRead(read)
			except TypeError:
				continue
			discoveryChimeras[newRead.data] += 1
	for ps,num in discoveryChimeras.items():
		if ps[1] > 24870500 and ps[4] < 24872500 and ps[7] == 'del':
			print(sample + '\t' + organism + '\t' + 'NC_036789.1\t' + str(ps[1]) + '\t' + str(ps[4]) + '\t' + ps[7] + '\t' + str(num))

		if ps[1] > 24986100 and ps[4] < 24987000 and ps[7] == 'del':
			print(sample + '\t' + organism + '\t' + 'NC_036789.1\t' + str(ps[1]) + '\t' + str(ps[4]) + '\t' + ps[7] + '\t' + str(num))

		if ps[1] > 24997000 and ps[4] < 24998200 and ps[7] == 'del':
			print(sample + '\t' + organism + '\t' + 'NC_036789.1\t' + str(ps[1]) + '\t' + str(ps[4]) + '\t' + ps[7] + '\t' + str(num))

