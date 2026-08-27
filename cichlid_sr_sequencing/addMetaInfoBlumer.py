import argparse, subprocess, os, urllib, shutil, contextlib, datetime, sys, pdb
import pandas as pd
import numpy as np
from helper_modules.file_manager import FileManager as FM



# Download and open master sample database file and read it in
fm_obj = FM()

fm_obj.readSampleDatabase()
fm_obj.downloadData(fm_obj.localReadDownloadDir + 'BlumerPaperData.csv')

dt1 = pd.read_csv(fm_obj.localReadDownloadDir + 'BlumerPaperData.csv')
dt2 = pd.read_csv('../TempBioProjectData/Table_S1.tsv', sep = '\t')
dt3 = pd.merge(dt1, dt2[['Biosample','clade']], left_on = 'BioSample', right_on='Biosample', how = 'left')
dt3.to_csv(fm_obj.localReadDownloadDir + 'BlumerPaperData_withClades.csv')
pdb.set_trace()
# Download and open run info file that contains new data to include
#fm_obj.downloadData(fm_obj.localReadDownloadDir)

dt = pd.read_csv('../TempBioProjectData/Table_S1.tsv', sep = '\t')
dt2 = pd.read_csv('../TempBioProjectData/Table_S7.tsv', sep = '\t')
dt3 = pd.merge(dt,dt2,left_on='sequence_id', right_on='sample_id')
for idx, row in dt3.iterrows():
	if row.Biosample in fm_obj.sample_dt.SampleID.tolist():
		if not pd.isna(row.chr10):
			fm_obj.sample_dt.loc[fm_obj.sample_dt.SampleID == row.Biosample,'Inversion10'] = int(row.chr10)

fm_obj.setSampleDatabase()