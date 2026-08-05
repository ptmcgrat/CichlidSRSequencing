from helper_modules.file_manager import FileManager as FM
import pandas as pd
from xml.etree import ElementTree as ET
import pdb

class SampleWorker():
	def __init__(self, fileManager):
		self.fileManager = fileManager

	def _addSeqCoreData(self, coreID, datatype):
		self.coreID = coreID
		self.datatype = datatype
		self.fileManager._createSeqCoreFiles(coreID)
		pdb.set_trace()


		tree = ET.parse(self.fileManager.localRunInfoFile).getroot()
		run_info = tree.findall('Run')[0]
		flowcell = run_info.findall('Flowcell')[0].text
		instrument = run_info.findall('Instrument')[0].text
		date = run_info.findall('Date')[0].text

		dt = pd.read_csv(self.fileManager.localSampleInfoFile, header = 19)


		out_dt = pd.DataFrame(columns = ['SampleID','Datatype','Date','Paired','RG','Files'])
		for sample in dt.Sample_ID:
			for root,dirs,files in os.walk(self.localCoreDir):
				for d in dirs:
					if sample + '_L' in d:
						fqs = [x for x in os.listdir(root + '/' + d) if '.fastq.gz' in x]
						try:
							lane = [x for x in fqs[0].split('_') if x[0] == 'L'][0]
						except IndexError:
							print(d)
							continue
							#pdb.set_trace()
						read_group = '@RG\\tID:' + flowcell + '.' + lane + '.' + sample + '\\tLB:' + sample + '\\tSM:' + sample + '\\tPL:ILLUMINA'
						#print(read_group)

						try:
							out_dt = out_dt.append({'SampleID':sample,'Datatype': 'GenomicDNA', 'Date':date,'Paired':'True','RG':read_group, 'Files':coreID + '/' + fqs[0] + ',,' + coreID + '/' + fqs[1]}, ignore_index=True)
						except IndexError:
							pdb.set_trace()
						#print(sample + '\t' + str(date) + '\t' + 'True' + '\t' + read_group + '\t' + fqs[0] + ',,' + fqs[1])
						subprocess.run(['mv', root + '/' + d + '/' + fqs[0], self.localReadsDir + coreID])
						subprocess.run(['mv', root + '/' + d + '/' + fqs[1], self.localReadsDir + coreID])
						
						try:
							print('Uploading ' + fqs[0])
							self.uploadData(self.localReadsDir + coreID + '/' + fqs[0])
							self.uploadData(self.localReadsDir + coreID + '/' + fqs[1])
						except Exception:
							print('Cant upload ' + str(fqs[0]))
			if sample not in list(out_dt.SampleID):
				print('No sample data for: ' + sample)
		out_dt.to_csv(self.localReadsDir + coreID + '_ToBeAdded.csv')
		self.uploadData(self.localReadsDir + coreID + '_ToBeAdded.csv')

fm_obj = FM()
sp_obj = SampleWorker(fm_obj)
sp_obj._addSeqCoreData('PM17', 'GenomicDNA')