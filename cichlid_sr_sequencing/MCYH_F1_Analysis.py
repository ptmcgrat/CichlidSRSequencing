import pysam, pdb, math

from collections import defaultdict

from helper_modules.file_manager_Replacement import FileManager as FM
from helper_modules.ChimericData import ChimericRead

LG10_inversion = ('NC_036789.1',11817927,11841666,29881282,29902558)
parents = ['YH_007_f','YH_006_f','YH_010_f','YH_013_f','YH_005_m','MC-4O5O-f','MC-3P6P-f','MC-2B4B-f','MC-2G8G-f','MC-5B6B-f']
minMapQ = 30

fm_obj = FM(genome_version = 'Mzebra_GT3')
fm_obj.downloadData(fm_obj.localGenomeFile)
refObj = pysam.FastaFile(fm_obj.localGenomeFile)

discoveryChimeras = {}

for parent in parents:
	fm_obj.createSampleFiles(parent)
	discoveryChimeras[parent] = defaultdict(int)
	fm_obj.downloadData(fm_obj.localChimericBamFile)
	pysam.index(fm_obj.localChimericBamFile)
	
for parent in parents:
	fm_obj.createSampleFiles(parent)
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
			discoveryChimeras[parent][newRead.data] += 1
	
for ps,num in discoveryChimeras['YH_005_m'].items():
	if num <5:
		continue
	length = ps[4] - ps[1]
	if length > 80 and length < 200 and ps[6] == '' and ps[7] == 'del':
		nums = [discoveryChimeras[x][ps] for x in parents]
		if sum(nums[5:]) == 0 and sum(nums[:4]) == 0:
			print('NC_036789.1\t' + str(ps[1]) + '\t' + str(ps[4]) + '\t' + refObj['NC_036789.1'][ps[1] - 400:ps[1]] + '[' + refObj['NC_036789.1'][ps[1]:ps[4]] + ']' + refObj['NC_036789.1'][ps[4]:ps[4]+400])


pdb.set_trace()