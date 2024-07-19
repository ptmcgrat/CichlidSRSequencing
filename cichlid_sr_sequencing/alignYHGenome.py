from helper_modules.file_manager_Replacement import FileManager as FM
from matplotlib.backends.backend_pdf import PdfPages

from Bio import AlignIO
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import subprocess, pdb

#https://dgenies.toulouse.inra.fr/

def align_genomes_contigbycontig(genome_version1,genome_version2,contig_mapping):
	fm_obj_1 = FM(genome_version = genome_version1)
	fm_obj_2 = FM(genome_version = genome_version2)

	all_dt = pd.DataFrame(columns = ['R_Name','R_Start','R_Stop','Q_Name','Q_Start','Q_Stop','Q_Strand','ResiduesMatch','AlignmentLength','MappingQuality','AlignmentType','PercentMatch','NewStart','NewStop'])
	for contig1,contig2 in contig_mapping.items():
		#subprocess.run(['faidx', fm_obj_1.localGenomeFile, contig1, '-o', fm_obj_1.localGenomeDir + contig1 + '.fa'])
		#subprocess.run(['faidx', fm_obj_2.localGenomeFile, contig2, '-o', fm_obj_2.localGenomeDir + contig2 + '.fa'])
		#subprocess.run(['minimap2', fm_obj_1.localGenomeDir + contig1 + '.fa', fm_obj_2.localGenomeDir + contig2 + '.fa'], stdout = open(fm_obj_1.localTempDir + contig1 + '_' + genome_version1 + '_' + genome_version2 + '.paf', 'w'))

		dt = pd.read_csv(fm_obj_1.localTempDir + contig1 + '_' + genome_version1 + '_' + genome_version2 + '.paf', sep = '\t', 
			names = ['Q_Name','Q_Size','Q_Start','Q_Stop','Q_Strand','R_Name','R_Size','R_Start','R_Stop','ResiduesMatch','AlignmentLength','MappingQuality','AlignmentType','c','d','e','f','g'])
		dt = dt[['R_Name','R_Start','R_Stop','Q_Name','Q_Start','Q_Stop','Q_Strand','ResiduesMatch','AlignmentLength','MappingQuality','AlignmentType']].sort_values('R_Start')
		dt['PercentMatch'] = dt.ResiduesMatch/dt.AlignmentLength
		dt['NewStart'] = np.where(dt.Q_Strand == '-', dt.Q_Stop,dt.Q_Start)
		dt['NewStop'] = np.where(dt.Q_Strand == '-', dt.Q_Start,dt.Q_Stop)
		all_dt = pd.concat([all_dt,dt])
		filter_dt = all_dt[(all_dt.PercentMatch > 0.3) & (all_dt.AlignmentType == 'tp:A:P') & (all_dt.AlignmentLength > 300)]
		
	all_dt.to_csv(fm_obj_1.localGenomesComparisonDir + genome_version1 + '_' + genome_version2 + '_contig_to_contig.csv')
	with PdfPages(fm_obj_1.localGenomesComparisonDir + genome_version1 + '_' + genome_version2 + '_contig_to_contig.pdf') as pdf_pages:
		for i,contig in enumerate(linkageGroups.keys()):
			lg = linkageGroups[contig]
			contig_dt = filter_dt[filter_dt.R_Name == contig][['R_Start','R_Stop']].melt(var_name = 'AlnLoc', value_name = 'Position', ignore_index = False)
			contig_dt['Position2'] = filter_dt[filter_dt.R_Name == contig][['NewStart','NewStop']].melt(var_name = 'AlnLoc', value_name = 'Position', ignore_index = False)['Position']
			figu = plt.figure(i)
			lineplot = sns.lineplot(data = contig_dt.reset_index(), x = 'Position', y = 'Position2', hue = 'index', hue_norm = (-255,0)).set(title = lg)
			pdf_pages.savefig(figu)

def align_genomes(genome_version1,genome_version2):
	fm_obj_1 = FM(genome_version = genome_version1)
	fm_obj_2 = FM(genome_version = genome_version2)

	out_paf = fm_obj_1.localGenomesComparisonDir + genome_version1 + '_' + genome_version2 + '.paf',
	out_csv = fm_obj_1.localGenomesComparisonDir + genome_version1 + '_' + genome_version2 + '_whole_genome.csv'
	out_pdf = fm_obj_1.localGenomesComparisonDir + genome_version1 + '_' + genome_version2 + '_whole_genome.pdf'


	subprocess.run(['minimap2', fm_obj_1.localGenomeFile, fm_obj_2.localGenomeFile], stdout = open(out_paf,'w'))

	all_dt = pd.read_csv(out_paf, sep = '\t', 
		names = ['Q_Name','Q_Size','Q_Start','Q_Stop','Q_Strand','R_Name','R_Size','R_Start','R_Stop','ResiduesMatch','AlignmentLength','MappingQuality','AlignmentType','c','d','e','f','g'])
	all_dt = all_dt[['R_Name','R_Start','R_Stop','Q_Name','Q_Start','Q_Stop','Q_Strand','ResiduesMatch','AlignmentLength','MappingQuality','AlignmentType']].sort_values(['R_Name','R_Start'])

	all_dt['PercentMatch'] = all_dt.ResiduesMatch/all_dt.AlignmentLength
	all_dt['NewStart'] = np.where(all_dt.Q_Strand == '-', all_dt.Q_Stop,all_dt.Q_Start)
	all_dt['NewStop'] = np.where(all_dt.Q_Strand == '-', all_dt.Q_Start,all_dt.Q_Stop)

	filter_dt = all_dt[(all_dt.PercentMatch > 0.3) & (all_dt.AlignmentType == 'tp:A:P') & (all_dt.AlignmentLength > 300)]

	all_dt.to_csv(out_csv)

	with PdfPages(out_pdf) as pdf_pages:
		for i,contig in enumerate(linkageGroups.keys()):
			lg = linkageGroups[contig]
			t_dt = filter_dt[filter_dt.R_Name == contig][['R_Start','R_Stop','Q_Name','NewStart','NewStop','AlignmentLength']]
			tophits = t_dt.groupby('Q_Name').sum()['AlignmentLength'].sort_values(ascending = False).head(3)
			t_dt = t_dt[t_dt.Q_Name.isin(tophits.index)]

			contig_dt = t_dt[['R_Start','R_Stop','Q_Name']].melt(id_vars = ['Q_Name'], var_name = 'AlnLoc', value_name = 'Position', ignore_index = False)
			contig_dt['Position2'] = t_dt[['NewStart','NewStop']].melt(var_name = 'AlnLoc', value_name = 'Position', ignore_index = False)['Position']

			figu = plt.figure(i)
			lineplot = sns.lineplot(data = contig_dt.reset_index(), x = 'Position', y = 'Position2', hue = 'Q_Name', units = 'index', estimator = None).set(title = lg)
			#lineplot = sns.lineplot(data = contig_dt.reset_index(), x = 'Position', y = 'Position2', hue = 'index', hue_norm = (-255,0)).set(title = lg)
			pdf_pages.savefig(figu)


#inversions = {'LG2':('NC_036781.1',19705000,19748000,43254805,43658853),'LG9':('NC_036789.1',14453796,15649299,32255605,33496468),'LG10':('NC_036790.1',11674905,11855817,29898615,29898615),
				#'LG11':('NC_036791.1',8302039,8309764,30371888,30459686),'LG12':('NC_036792.1',2249541,2453698,23046928,23131968),'LG20':('NC_036799.1',19614379,19689710,32872827,33764042)}

inversions = {'LG2':('NC_036781.1',19745291,19746383,43556654,43556721),'LG9':('NC_036788.1',14503645,14972797,32763228,32765988),'LG10':('NC_036789.1',11817927,11841666,29881282,29902558),
				'LG11':('NC_036790.1',8304758,8303033,30372595,30385320),'LG13':('NC_036792.1',2398101,2405661,23071814,23124402),'LG20':('NC_036799.1',19616371,19616393,33288065,33254389)}

linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
							  'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
							  'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
							  'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

LG_MZtoON = {'NC_036780.1':'NC_031965.2', 'NC_036781.1':'NC_031966.2', 'NC_036782.1':'NC_031967.2', 'NC_036783.1':'NC_031969.2', 'NC_036784.1':'NC_031970.2', 'NC_036785.1':'NC_031971.2', 
							  'NC_036786.1':'NC_031972.2', 'NC_036787.1':'NC_031973.2', 'NC_036788.1':'NC_031974.2', 'NC_036789.1':'NC_031975.2', 'NC_036790.1':'NC_031976.2',
							  'NC_036791.1':'NC_031977.2', 'NC_036792.1':'NC_031978.2', 'NC_036793.1':'NC_031979.2', 'NC_036794.1':'NC_031980.2', 'NC_036795.1':'NC_031987.2', 'NC_036796.1':'NC_031981.2',
							  'NC_036797.1':'NC_031982.2', 'NC_036798.1':'NC_031983.2', 'NC_036799.1':'NC_031984.2', 'NC_036800.1':'NC_031985.2', 'NC_036801.1':'NC_031986.2'}

LG_MZtoYH = {x:x for x in linkageGroups.keys()}

genome_versions = ['Mzebra_GT3','kocher_Mzebra_female','kocher_YH_female','kocher_YH_male','YH_3','O_niloticus_UMD_NMBU','P_nyererei_v2','Rhamp_chilingali']
fm_objs = {}

for gv in genome_versions:
	fm_objs[gv] = FM(genome_version = gv)
	fm_objs[gv].downloadData(fm_objs[gv].localGenomeFile)
		
for gv in genome_versions:
	if gv == 'Mzebra_GT3':
		continue
	align_genomes('Mzebra_GT3',gv)



#fm_obj_mz.downloadData(fm_obj_mz.localGenomeFile)
#fm_obj_pn.downloadData(fm_obj_pn.localGenomeFile)

#fm_obj_mz.createSampleFiles('YH_1_m')
fm_obj_mz.createSampleFiles('MZ_1_m')
#fm_obj_mz.downloadData(fm_obj_mz.localSampleBamDir)
#fm_obj_yh.downloadData(fm_obj_yh.localGenomeFile)
#fm_obj_on.downloadData(fm_obj_on.localGenomeFile)
#fm_obj_yhhf.downloadData(fm_obj_yhhf.localGenomeFile)
#align_genomes_contigbycontig('Mzebra_GT3','O_niloticus_UMD_NMBU',LG_MZtoON)
#align_genomes_contigbycontig('Mzebra_GT3','kocher_YH_female',LG_MZtoYH)
#align_genomes('Mzebra_GT3','P_nyererei_v2')
align_genomes('Mzebra_GT3','O_niloticus_UMD_NMBU')
#subprocess.run(['GSAlign','-dp','-i',fm_obj_mz.localGenomeFile,'-q',fm_obj_yh.localGenomeFile, '-o', fm_obj_mz.localGenomesDir + 'MZ_YH_Alignment'])

fm_obj_mz.uploadData(fm_obj_mz.localGenomesComparisonDir)

"""with PdfPages(fm_obj_mz.localGenomesDir + 'Comparisons/MZ_GT3vsON.pdf') as pdf_pages:
	for i,(contig,lg) in enumerate(linkageGroups.items()):
		a_dt = g_dt[g_dt.R_Name == contig]
		a_dt = g_dt[g_dt.R_Name == contig][['R_Start','R_Stop']].melt(var_name = 'AlnLoc', value_name = 'Position', ignore_index = False)
		a_dt['Position2'] = g_dt[g_dt.R_Name == contig][['NewStart','NewStop']].melt(var_name = 'AlnLoc', value_name = 'Position', ignore_index = False)['Position']
		figu = plt.figure(i)
		lineplot = sns.lineplot(data = a_dt.reset_index(), x = 'Position', y = 'Position2', hue = 'index', hue_norm = (-255,0)).set(title = lg)
		pdf_pages.savefig(figu)
"""


# Old code

"""
processes = []

for contig, lg in linkageGroups.items():
	print('Running :' + lg)
	#subprocess.run(['faidx', fm_obj_mz.localGenomeFile, contig, '-o', fm_obj_mz.localGenomesDir + lg + '_MZ.fa'])
	#subprocess.run(['faidx', fm_obj_yh.localGenomeFile, contig, '-o', fm_obj_mz.localGenomesDir + lg + '_YH.fa'])

	command = ['minimap2', fm_obj_mz.localGenomesDir + lg + '_MZ.fa', fm_obj_mz.localGenomesDir + lg + '_YH.fa', '-o', fm_obj_mz.localGenomesDir + lg + '.paf']
	#command = ['lastz', fm_obj_mz.localGenomesDir + inversions[inv_contig][0] + '_MZ.fa', fm_obj_mz.localGenomesDir + inversions[inv_contig][0] + '_YH.fa']
	#command += ['--strand=both','--nochain','--gap=600,150', '--gappedthresh=3000', '--masking=254', '--hspthresh=4500']
	#command += ['--allocate:traceback=1.99G', '--output=' + fm_obj_mz.localGenomesDir + inv_contig + '.maf', '--format=maf']
	#subprocess.run(command)
	#processes.append(subprocess.Popen(command))
#for p in processes:
#	p.communicate()
with PdfPages('Test.pdf') as pdf_pages:
	for i,(contig,lg) in enumerate(linkageGroups.items()):

		#dt = pd.read_csv(fm_obj_mz.localGenomesDir + lg + '.paf', sep = '\t', 
		#names = ['Q_Name','Q_Size','Q_Start','Q_Stop','Q_Strand','R_Name','R_Size','R_Start','R_Stop','ResiduesMatch','AlignmentLength','MappingQuality','AlignmentType','c','d','e','f','g'])
		#dt = dt[['R_Name','R_Start','R_Stop','Q_Name','Q_Start','Q_Stop','Q_Strand','ResiduesMatch','AlignmentLength','MappingQuality','AlignmentType']].sort_values('R_Start')
		#dt['PercentMatch'] = dt.ResiduesMatch/dt.AlignmentLength
		#dt['NewStart'] = np.where(dt.Q_Strand == '-', dt.Q_Stop,dt.Q_Start)
		#dt['NewStop'] = np.where(dt.Q_Strand == '-', dt.Q_Start,dt.Q_Stop)
		#g_dt = dt[(dt.PercentMatch > 0.5) & (dt.AlignmentType == 'tp:A:P') & (dt.AlignmentLength > 10000)]
		#a_dt = g_dt[['ResiduesMatch','R_Start','R_Stop']].melt(id_vars = 'ResiduesMatch', var_name = 'AlnLoc', value_name = 'Position')
		#a_dt['Position2'] = g_dt[['ResiduesMatch','NewStart','NewStop']].melt(id_vars = 'ResiduesMatch', var_name = 'AlnLoc', value_name = 'Position')['Position']
		pdb.set_trace()


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
"""