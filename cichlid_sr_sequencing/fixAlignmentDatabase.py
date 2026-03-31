from helper_modules.file_manager import FileManager as FM
from helper_modules.alignment_worker import AlignmentWorker as AW
import datetime, pdb

for genome_version in ['Mzebra_GT3','Mconophoros_GT1']:
	fm_obj = FM(genome_version = 'Mzebra_GT3')
	fm_obj.setSamples(None,None,None,False)
	aw_obj = AW(genome_version, fm_obj, check_size = False)
	samples = fm_obj.returnCloudDirs(fm_obj.localBamRefDir)
	for sample in samples:
		print(sample + ', ' + str(datetime.datetime.now()))
		fm_obj.createSampleFiles(sample)
		fm_obj.downloadData(fm_obj.localSampleBamDir)
		aw_obj.uploadAndUpdateDatabase(upload = False, sample_override = sample)

		if sample not in fm_obj.sample_dt.SampleID.to_list():
			print('Adding ' + sample + ' to SampleDatbase')
			#sample_row = {'SampleID':sample,'Sex':'','Species':'','DoB':'','BroodID':'','Parents':'','Ecogroup':'','LabReared':''}
			#fm_obj.sample_dt = fm_obj.sample_dt.append(pd.Series(sample_row), ignore_index = True)
			#fm_obj._setDatabase('SampleDatabase', fm_obj.sample_dt)

