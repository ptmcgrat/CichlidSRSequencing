import pysam, pdb, math

from collections import defaultdict

from helper_modules.file_manager_Replacement import FileManager as FM
from helper_modules.ChimericData import ChimericRead

LG10_inversion = ('NC_036789.1',11817927,11841666,29881282,29902558)
parents = ['YH_007_f','YH_006_f','YH_010_f','YH_013_f','MC-4O5O-f','MC-3P6P-f','MC-2B4B-f','MC-2G8G-f','MC-5B6B-f']
minMapQ = 30

fm_obj = FM(genome_version = 'Mzebra_GT3')
fm_obj.downloadData(fm_obj.localGenomeFile)
refObj = pysam.FastaFile(fm_obj.localGenomeFile)

discoveryChimeras = defaultdict(int)
exclusionChimeras = defaultdict(int)

for parent in parents:
	fm_obj.createSampleFiles(parent)
	fm_obj.downloadData(fm_obj.localChimericBamFile)
	pysam.index(fm_obj.localChimericBamFile)
	
for parent in parents:
	fm_obj.createSampleFiles(parent)
	bam_obj = pysam.AlignmentFile(fm_obj.localChimericBamFile)
	#for read in bam_obj.fetch(LG10_inversion[0],LG10_inversion[2],LG10_inversion[3]):
	for read in bam_obj.fetch():

		if not read.is_secondary and read.mapq > minMapQ:
			try:
				SA = read.get_tag('SA')
			except KeyError:
				continue
			#if SA.split(',')[0]!=LG10_inversion[0]:
			#	continue
			read.set_tag('SA',SA.replace(SA.split(',')[0], str(refObj.references.index(SA.split(',')[0]))))

			try:
				newRead = ChimericRead(read)
			except TypeError:
				continue
			if 'YH' in parent:
				discoveryChimeras[newRead.data] += 1
			else:
				exclusionChimeras[newRead.data] += 1

candidates = [x for x in discoveryChimeras if discoveryChimeras[x] > 10 and exclusionChimeras[x] == 0]
for candidate in candidates:
	print('\t'.join([str(x) for x in candidate]) + '\t' + str(discoveryChimeras[candidate])
#			print('NC_036789.1\t' + str(ps[1]) + '\t' + str(ps[4]) + '\t' + refObj['NC_036789.1'][ps[1] - 400:ps[1]] + '[' + refObj['NC_036789.1'][ps[1]:ps[4]] + ']' + refObj['NC_036789.1'][ps[4]:ps[4]+400])


