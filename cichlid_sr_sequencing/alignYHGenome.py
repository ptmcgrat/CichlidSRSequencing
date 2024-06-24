from helper_modules.file_manager_Replacement import FileManager as FM
from Bio import AlignIO
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import subprocess, pdb

#https://dgenies.toulouse.inra.fr/

inversions = {'LG2':('NC_036781.1',19705000,19748000,43254805,43658853),'LG9':('NC_036789.1',14453796,15649299,32255605,33496468),'LG10':('NC_036790.1',11674905,11855817,29898615,29898615),
				'LG11':('NC_036791.1',8302039,8309764,30371888,30459686),'LG12':('NC_036792.1',2249541,2453698,23046928,23131968),'LG20':('NC_036799.1',19614379,19689710,32872827,33764042)}

linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
							  'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
							  'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
							  'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

fm_obj_mz = FM(genome_version = 'Mzebra_GT3')
fm_obj_yh = FM(genome_version = 'kocher_YH_female')

#fm_obj_mz.downloadData(fm_obj_mz.localGenomeFile)
#fm_obj_yh.downloadData(fm_obj_yh.localGenomeFile)

#subprocess.run(['GSAlign','-dp','-i',fm_obj_mz.localGenomeFile,'-q',fm_obj_yh.localGenomeFile, '-o', fm_obj_mz.localGenomesDir + 'MZ_YH_Alignment'])

processes = []
for contig, lg in linkageGroups.items():
	print('Running :' + lg)
	subprocess.run(['faidx', fm_obj_mz.localGenomeFile, contig, '-o', fm_obj_mz.localGenomesDir + lg + '_MZ.fa'])
	subprocess.run(['faidx', fm_obj_yh.localGenomeFile, contig, '-o', fm_obj_mz.localGenomesDir + lg + '_YH.fa'])

	command = ['minimap2', fm_obj_mz.localGenomesDir + lg + '_MZ.fa', fm_obj_mz.localGenomesDir + lg + '_YH.fa', '-o', fm_obj_mz.localGenomesDir + lg + '.paf']
	#command = ['lastz', fm_obj_mz.localGenomesDir + inversions[inv_contig][0] + '_MZ.fa', fm_obj_mz.localGenomesDir + inversions[inv_contig][0] + '_YH.fa']
	#command += ['--strand=both','--nochain','--gap=600,150', '--gappedthresh=3000', '--masking=254', '--hspthresh=4500']
	#command += ['--allocate:traceback=1.99G', '--output=' + fm_obj_mz.localGenomesDir + inv_contig + '.maf', '--format=maf']
	subprocess.run(command)
	#processes.append(subprocess.Popen(command))
#for p in processes:
#	p.communicate()

for lg in linkageGroups.values():
	dt = pd.read_csv(fm_obj_mz.localGenomesDir + lg + '.paf', sep = '\t', 
	names = ['Q_Name','Q_Size','Q_Start','Q_Stop','Q_Strand','R_Name','R_Size','R_Start','R_Stop','ResiduesMatch','AlignmentLength','MappingQuality','AlignmentType','c','d','e','f','g'])
	dt = dt[['R_Name','R_Start','R_Stop','Q_Name','Q_Start','Q_Stop','Q_Strand','ResiduesMatch','AlignmentLength','MappingQuality','AlignmentType']].sort_values('R_Start')
	dt['PercentMatch'] = dt.ResiduesMatch/dt.AlignmentLength
	dt['NewStart'] = np.where(dt.Q_Strand == '-', dt.Q_Stop,dt.Q_Start)
	dt['NewStop'] = np.where(dt.Q_Strand == '-', dt.Q_Start,dt.Q_Stop)
	g_dt = dt[(dt.PercentMatch > 0.5) & (dt.AlignmentType == 'tp:A:P') & (dt.AlignmentLength > 10000)]
	a_dt = g_dt[['ResiduesMatch','R_Start','R_Stop']].melt(id_vars = 'ResiduesMatch', var_name = 'AlnLoc', value_name = 'Position')
	a_dt['Position2'] = g_dt[['ResiduesMatch','NewStart','NewStop']].melt(id_vars = 'ResiduesMatch', var_name = 'AlnLoc', value_name = 'Position')['Position']
	sns.lineplot(data = a_dt, x = 'Position', y = 'Position2', hue = 'ResiduesMatch').set(title = lg)
	plt.show()



aln = AlignIO.parse(fm_obj_mz.localGenomesDir + 'MZ_YH_Alignment.maf', "maf")

dt = pd.read_csv('map_A_spYH_GT1_to_Mzebra_GT3.paf', sep = '\t', names = ['QueryContig','QueryContigSize','QueryStart','QueryStop','QueryStrand','RefContig','RefContigSize','RefStart','RefStop','ResiduesMatch','AlignmentBlockLength','MappingQuality','b','c','d','e','f','g'])
dt = dt[['RefContig','RefStart','RefStop','QueryContig','QueryStart','QueryStop','QueryStrand','ResiduesMatch','AlignmentBlockLength','MappingQuality']]
LG = dt[(dt.QueryContig == 'NC_036781.1') & (dt.RefContig == 'NC_036781.1')].sort_values('RefStart')
pdb.set_trace()




dt = dt[['RefContig','RefStart','RefStop','QueryContig','QueryStart','QueryStop','QueryStrand']]
mz = []
yh = []
dt = pd.DataFrame(columns = ['RefContig','RefStart','RefStop','RefSize','QueryContig','QueryStart','QueryStop','QuerySize'])
for ma in aln:
	ref_contig = ma[0].id 
	ref_start = ma[0].annotations['start']
	ref_length = ma[0].annotations['size']
	ref_stop = ref_start + ma[0].annotations['strand'] * ref_lengthLG
	qry_contig = ma[1].id 
	qry_start = ma[1].annotations['start']
	qry_length = ma[1].annotations['size']
	qry_stop = qry_start + ma[1].annotations['strand'] * qry_length

	dt.loc[len(dt)] = [ref_contig,ref_start,ref_stop,ref_length,qry_contig,qry_start,qry_stop,qry_length]

LG = dt[(dt.QueryContig == 'NC_036781.1') & (dt.TargetContig == 'NC_036781.1')].sort_values('QueryStart')
pdb.set_trace()

plt.scatter(mz,yh)
plt.show()
