import argparse, subprocess, os, urllib, shutil, contextlib, datetime, sys, pdb
import pandas as pd
import numpy as np
from helper_modules.file_manager import FileManager as FM



# Download and open master sample database file and read it in
fm_obj = FM()

fm_obj.readSampleDatabase()
# Download and open run info file that contains new data to include
#fm_obj.downloadData(fm_obj.localReadDownloadDir)

dt = pd.read_csv('../TempBioProjectData/Table_S1.tsv', sep = '\t')
dt2 = pd.read_csv('../TempBioProjectData/Table_S7.tsv', sep = '\t')
dt3 = pd.merge(dt,dt2,left_on='sequence_id', right_on='sample_id')
for idx, row in dt3.iterrows():
	if row.Biosample in fm_obj.sample_dt.SampleID.tolist():
		if not pd.isna(row.chr10):
			fm_obj.sample_dt.loc[fm_obj.sample_dt.SampleID == row.Biosample,'Inversion10'] = str(row.chr10)

fm_obj.setSampleDatabase()