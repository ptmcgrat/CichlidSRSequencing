from helper_modules.file_manager_Replacement import FileManager as FM
import pdb

fm_obj = FM('Mzebra_UMD2a')

outfiles = fm_obj.returnCloudFiles(fm_obj.localReadsDir + 'BigBrain/')

for i,x in enumerate([x for x in outfiles if x.endswith('.bam')]):
	runID = x.split('-')[1].split('.')[0]
	sampleID = runID.split('_L')[0]
	#print(','.join([sampleID,runID,'300','0','BigBrain','Illumina NovaSeq',sampleID + '_library','PAIRED','GENOMIC','ILLUMINA','BigBrain/' + x,'Mchenga conophoros','Shallow_Benthic']))
	print('BigBrain/' + x)
	#if i %2 == 0:
	#	print(sampleID)
##pdb.set_trace()