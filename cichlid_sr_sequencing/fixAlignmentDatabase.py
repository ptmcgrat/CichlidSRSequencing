from helper_modules.file_manager import FileManager as FM
from helper_modules.alignment_worker import AlignmentWorker as AW
import datetime, pdb
import pandas as pd
fm_obj = FM('Mzebra_GT3')
fm_obj.readSampleDatabase()
fm_obj.createSampleFiles('MCYHBC1-809-1')
for uBam in fm_obj.localRawBamFiles:
	fm_obj.downloadData(uBam)
pdb.set_trace()

for i,row in fm_obj.reads_dt.iterrows():
	if pd.isna(fm_obj.reads_dt.loc[i,'FileSize']):
		size = 0
		for ubam in row.FileLocations.split(','):
			try:
				size+=fm_obj.returnFileSize(fm_obj.localReadsDir + ubam)
			except ValueError:
				print('Bad size for ' + row.SampleID)
		fm_obj.reads_dt.loc[fm_obj.reads_dt.FileLocations == row.FileLocations,'FileSize'] = size
	if row.SampleID not in fm_obj.sample_dt.SampleID.tolist():
		print(row.SampleID + ' not in sample database')
		fm_obj.sample_dt.loc[len(fm_obj.sample_dt)] = [row.SampleID] + ['']*8
fm_obj.setSampleDatabase()
	

pdb.set_trace()


"""
for genome_version in ['Mconophoros_GT1']:
	fm_obj = FM(genome_version = genome_version)
	fm_obj.setSamples(None,None,None,False)
	aw_obj = AW(genome_version, fm_obj, check_size = False)
	samples = fm_obj.returnCloudDirs(fm_obj.localBamRefDir)
	for sample in samples:
		if sample in fm_obj.alignment_dt[(fm_obj.alignment_dt.GenomeVersion == genome_version)].SampleID.to_list():
			print(sample + ' already run. Skipping...')
			continue
		print(sample + ', ' + str(datetime.datetime.now()))
		fm_obj.createSampleFiles(sample)
		fm_obj.downloadData(fm_obj.localSampleBamDir)
		aw_obj.uploadAndUpdateDatabase(upload = False, sample_override = sample)

		if sample not in fm_obj.sample_dt.SampleID.to_list():
			print('Adding ' + sample + ' to SampleDatbase')
			#sample_row = {'SampleID':sample,'Sex':'','Species':'','DoB':'','BroodID':'','Parents':'','Ecogroup':'','LabReared':''}
			#fm_obj.sample_dt = fm_obj.sample_dt.append(pd.Series(sample_row), ignore_index = True)
			#fm_obj._setDatabase('SampleDatabase', fm_obj.sample_dt)
"""
