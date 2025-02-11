from helper_modules.file_manager_Replacement import FileManager as FM
from matplotlib.backends.backend_pdf import PdfPages

from Bio import AlignIO
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import subprocess, pdb, os
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

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

def align_genomes(genome_version1,genome_version2, hs_flag = False):
	fm_obj_1 = FM(genome_version = genome_version1)
	fm_obj_2 = FM(genome_version = genome_version2)

	if hs_flag:
		out_paf = fm_obj_1.localGenomesComparisonDir + genome_version1 + '_' + genome_version2 + '.hs.paf'
		out_csv = fm_obj_1.localGenomesComparisonDir + genome_version1 + '_' + genome_version2 + '_whole_genome.hs.csv'
		out_pdf = fm_obj_1.localGenomesComparisonDir + genome_version1 + '_' + genome_version2 + '_whole_genome.hs.pdf'
	else:
		out_paf = fm_obj_1.localGenomesComparisonDir + genome_version1 + '_' + genome_version2 + '.paf'
		out_csv = fm_obj_1.localGenomesComparisonDir + genome_version1 + '_' + genome_version2 + '_whole_genome.csv'
		out_pdf = fm_obj_1.localGenomesComparisonDir + genome_version1 + '_' + genome_version2 + '_whole_genome.pdf'

	if not os.path.isfile(out_paf):
		if hs_flag:
			subprocess.run(['minimap2', fm_obj_1.localGenomeFile, fm_obj_2.localHybridScaffoldFile, '-t','12'], stdout = open(out_paf,'w'))
		else:
			subprocess.run(['minimap2', fm_obj_1.localGenomeFile, fm_obj_2.localGenomeFile, '-t','12'], stdout = open(out_paf,'w'))

	all_dt = pd.read_csv(out_paf, sep = '\t', 
		names = ['Q_Name','Q_Size','Q_Start','Q_Stop','Q_Strand','R_Name','R_Size','R_Start','R_Stop','ResiduesMatch','AlignmentLength','MappingQuality','AlignmentType','c','d','e','f','g'])
	all_dt = all_dt[['R_Name','R_Start','R_Stop','Q_Name','Q_Start','Q_Stop','Q_Strand','ResiduesMatch','AlignmentLength','MappingQuality','AlignmentType']].sort_values(['R_Name','R_Start'])

	all_dt['PercentMatch'] = all_dt.ResiduesMatch/all_dt.AlignmentLength
	all_dt['NewStart'] = np.where(all_dt.Q_Strand == '-', all_dt.Q_Stop,all_dt.Q_Start)
	all_dt['NewStop'] = np.where(all_dt.Q_Strand == '-', all_dt.Q_Start,all_dt.Q_Stop)

	filter_dt = all_dt[(all_dt.PercentMatch > 0.3) & (all_dt.AlignmentType == 'tp:A:P') & (all_dt.AlignmentLength > 10000)]

	filter_dt.to_csv(out_csv)
	chr_dt = filter_dt.groupby(['R_Name','Q_Name']).sum()['AlignmentLength'].reset_index()
	chr_dt = chr_dt[chr_dt.AlignmentLength > 2000000]

	with PdfPages(out_pdf) as pdf_pages:
		for i,contig in enumerate(linkageGroups.keys()):
			lg = linkageGroups[contig]
			t_dt = filter_dt[(filter_dt.R_Name == contig)][['R_Start','R_Stop','Q_Name','NewStart','NewStop','AlignmentLength']]
			#tophits = t_dt.groupby('Q_Name').sum()['AlignmentLength'].sort_values(ascending = False).head(3)
			#(filter_dt.Q_Name.isin(chr_dt[chr_dt.R_Name == filter_dt.R_Name].Q_Name)
			#t_dt = t_dt[t_dt.Q_Name.isin(tophits.index)]
			t_dt = t_dt[t_dt.Q_Name.isin(chr_dt[chr_dt.R_Name == contig].Q_Name)]
			contig_dt = t_dt[['R_Start','R_Stop','Q_Name']].melt(id_vars = ['Q_Name'], var_name = 'AlnLoc', value_name = 'Position', ignore_index = False)
			contig_dt['Position2'] = t_dt[['NewStart','NewStop']].melt(var_name = 'AlnLoc', value_name = 'Position', ignore_index = False)['Position']

			figu = plt.figure(i)
			lineplot = sns.lineplot(data = contig_dt.reset_index(), x = 'Position', y = 'Position2', hue = 'Q_Name', units = 'index', estimator = None).set(title = lg)
			#lineplot = sns.lineplot(data = contig_dt.reset_index(), x = 'Position', y = 'Position2', hue = 'index', hue_norm = (-255,0)).set(title = lg)
			pdf_pages.savefig(figu)
			figu.clf()


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

genome_versions = ['kocher_N_Met_zebra_Female','MZ4f_ptm','kocher_H_Aulon_yelhead_Female',
					'kocher_G_Aulon_yelhead_Male','YH7f_ptm','Rhamp_chilingali','O_niloticus_UMD_NMBU','P_nyererei_v2']
['O_niloticus_UMD_NMBU','P_nyererei_v2','Rhamp_chilingali']
fm_objs = {}

for gv in genome_versions:
	print(gv)
	fm_objs[gv] = FM(genome_version = gv)
	#fm_objs[gv].downloadData(fm_objs[gv].localHybridScaffoldFile)

fm_objs[gv].downloadData(fm_objs[gv].localGenomesComparisonDir)
	

for gv in genome_versions:
	if gv == 'Mzebra_GT3':
		continue
	#align_genomes('Mzebra_GT3',gv,hs_flag = True)
	"""if gv == 'O_niloticus_UMD_NMBU':
		align_genomes('Mzebra_GT3',gv,LG_MZtoON)
	elif gv == 'P_nyrerei_v2':
		align_genomes('Mzebra_GT3',gv,LG_MZtoPN)
	elif gv == 'Rhamp_chilingali':
		align_genomes('Mzebra_GT3',gv,LG_MZtoRC)
	else:
		align_genomes('Mzebra_GT3',gv,LG_MZtoYH)
	"""

for gv in ['kocher_N_Met_zebra_Female','MZ4f_ptm','kocher_H_Aulon_yelhead_Female',
					'kocher_G_Aulon_yelhead_Male','YH7f_ptm','Rhamp_chilingali','O_niloticus_UMD_NMBU','P_nyererei_v2']:

	if gv in ['Rhamp_chilingali','P_nyererei_v2','O_niloticus_UMD_NMBU']:
		a_dt = pd.read_csv(fm_objs[gv].localGenomesComparisonDir + 'Mzebra_GT3'+ '_' + gv + '_whole_genome.csv')
	else:
		a_dt = pd.read_csv(fm_objs[gv].localGenomesComparisonDir + 'Mzebra_GT3'+ '_' + gv + '_whole_genome.hs.csv')
	a_dt['GenomeVersion'] = gv

	chr_dt = a_dt.groupby(['R_Name','Q_Name']).sum()['AlignmentLength'].reset_index()
	chr_dt = chr_dt[chr_dt.AlignmentLength > 2500000]

	a_dt = pd.merge(left = a_dt,right = chr_dt[['R_Name','Q_Name']], on = ['R_Name','Q_Name'])

	try:
		total_dt = pd.concat([total_dt,a_dt])
	except:
		total_dt = a_dt

total_dt.to_csv(fm_objs[gv].localGenomesComparisonDir + 'YH_Reference_Comparisions_whole_genome.hs.csv')
pdb_out = fm_objs[gv].localGenomesComparisonDir + 'YH_Reference_Comparisions_whole_genome.hs.pdf'
with PdfPages(pdb_out) as pdf_pages:
	for i,contig in enumerate(linkageGroups.keys()):
		lg = linkageGroups[contig]
		t_dt = total_dt[(total_dt.R_Name == contig)][['R_Start','R_Stop','Q_Name','NewStart','NewStop','GenomeVersion']]
		contig_dt = t_dt[['R_Start','R_Stop','Q_Name']].melt(id_vars = ['Q_Name'], var_name = 'AlnLoc', value_name = 'Position', ignore_index = False)
		contig_dt['Position2'] = t_dt[['NewStart','NewStop']].melt(var_name = 'AlnLoc', value_name = 'Position', ignore_index = False)['Position']

		figu = plt.figure(i)
		lineplot = sns.lineplot(data = contig_dt.reset_index(), x = 'Position', y = 'Position2', hue = 'Q_Name', units = 'index', estimator = None).set(title = lg)
		#plt.legend(loc = "upper left", bbox_to_anchor=(1, 1))

		#lineplot = sns.lineplot(data = contig_dt.reset_index(), x = 'Position', y = 'Position2', hue = 'index', hue_norm = (-255,0)).set(title = lg)
		pdf_pages.savefig(figu)
		figu.clf()

pdb_out = fm_objs[gv].localGenomesComparisonDir + 'InversionSummary.hs.pdf'
with PdfPages(pdb_out) as pdf_pages:
	gv_names = ['MZf_2','MZf_3','YHf_1','YHm_2','YHf_3','RC']
	for lg,(contig,start1,start2,stop1,stop2) in inversions.items():
		fig,axes = plt.subplots(3,2, figsize = (6,9))
		fig.subplots_adjust(hspace=0.6)
		fig.subplots_adjust(wspace=0.6)
		fig.suptitle(lg)
		current_axes = [axes[0,0],axes[1,0],axes[0,1],axes[1,1],axes[2,1],axes[2,0]]

		t_dt = total_dt[(total_dt.R_Name == contig)][['R_Start','R_Stop','Q_Name','NewStart','NewStop','GenomeVersion']]
		for i,gv in enumerate(['kocher_N_Met_zebra_Female','MZ4f_ptm','kocher_H_Aulon_yelhead_Female','kocher_G_Aulon_yelhead_Male','YH7f_ptm', 'Rhamp_chilingali']):
			contig_dt = t_dt[t_dt.GenomeVersion == gv][['R_Start','R_Stop','Q_Name']].melt(id_vars = ['Q_Name'], var_name = 'AlnLoc', value_name = 'Position 1 (Mb)', ignore_index = False)
			contig_dt['Position 2 (Mb)'] = t_dt[t_dt.GenomeVersion == gv][['NewStart','NewStop']].melt(var_name = 'AlnLoc', value_name = 'Position 1 (Mb)', ignore_index = False)['Position 1 (Mb)']
			contig_dt['Position 1 (Mb)'] = contig_dt['Position 1 (Mb)']/1000000
			contig_dt['Position 2 (Mb)'] = contig_dt['Position 2 (Mb)']/1000000
			lineplot = sns.lineplot(data = contig_dt.reset_index(), x = 'Position 1 (Mb)', y = 'Position 2 (Mb)', hue = 'Q_Name', units = 'index', estimator = None, ax = current_axes[i], legend = True, palette = 'Dark2').set(title = gv_names[i])
			current_axes[i].add_patch(plt.Rectangle((start1/1000000, -10), (stop1-start1)/1000000, 100, edgecolor='grey', fill=True, alpha = 0.3, facecolor = 'grey'))
			if lg == 'LG11':
				current_axes[i].add_patch(plt.Rectangle((18425216/1000000, -10), 10000/1000000, 100, edgecolor='black', fill=True, alpha = 0.3, facecolor = 'grey'))
			if lg == 'LG20':
				current_axes[i].add_patch(plt.Rectangle((29716509/1000000, -10), 10000/1000000, 100, edgecolor='black', fill=True, alpha = 0.3, facecolor = 'grey'))

			#axes[2,0].set_axis_off()
		pdf_pages.savefig(fig)
		fig.clf()


pdb_out = fm_objs[gv].localGenomesComparisonDir + 'InversionAncestry.pdf'
with PdfPages(pdb_out) as pdf_pages:
	fig,axes = plt.subplots(3,2, figsize = (8.5,11))
	fig.subplots_adjust(hspace=0.6)
	fig.subplots_adjust(wspace=0.6)
	fig.suptitle('Outgroup alignments')
	current_axes = [axes[0,0],axes[0,1],axes[1,0],axes[1,1],axes[2,0],axes[2,1]]
	for i,(lg,(contig,start1,start2,stop1,stop2)) in enumerate(inversions.items()):
		t_dt = total_dt[(total_dt.R_Name == contig)][['R_Start','R_Stop','Q_Name','NewStart','NewStop','GenomeVersion']]
		short_gv = ['O_niloticus','P_nyererei']
		palettes = ['Reds_d','Blues_d']
		for j,gv in enumerate(['O_niloticus_UMD_NMBU','P_nyererei_v2']):
			contig_dt = t_dt[t_dt.GenomeVersion == gv][['R_Start','R_Stop','Q_Name']].melt(id_vars = ['Q_Name'], var_name = 'AlnLoc', value_name = 'Position 1 (Mb)', ignore_index = False)
			contig_dt['Position 2 (Mb)'] = t_dt[t_dt.GenomeVersion == gv][['NewStart','NewStop']].melt(var_name = 'AlnLoc', value_name = 'Position 1 (Mb)', ignore_index = False)['Position 1 (Mb)']
			contig_dt['Position 1 (Mb)'] = contig_dt['Position 1 (Mb)']/1000000
			contig_dt['Position 2 (Mb)'] = contig_dt['Position 2 (Mb)']/1000000
			lineplot = sns.lineplot(data = contig_dt.reset_index(), x = 'Position 1 (Mb)', y = 'Position 2 (Mb)', hue = 'Q_Name', units = 'index', estimator = None, ax = current_axes[i], legend = True, palette = palettes[j]).set(title = lg)
		current_axes[i].add_patch(plt.Rectangle((start1/1000000, -10), (stop1-start1)/1000000, 100, edgecolor='grey', fill=True, alpha = 0.3, facecolor = 'grey'))
		if lg == 'LG11':
			current_axes[i].add_patch(plt.Rectangle((18425216/1000000, -10), 10000/1000000, 100, edgecolor='black', fill=True, alpha = 0.3, facecolor = 'grey'))
		if lg == 'LG20':
			current_axes[i].add_patch(plt.Rectangle((29716509/1000000, -10), 10000/1000000, 100, edgecolor='black', fill=True, alpha = 0.3, facecolor = 'grey'))

	pdf_pages.savefig(fig)

fm_objs[gv].uploadData(fm_objs[gv].localGenomesComparisonDir)

# For 2
#faidx kocher_G_Aulon_yelhead_Male_anchored_assembly.fasta NC_036781.1:29631300-29631500 > YH_G_Inv2_L.fa
#NC_036781.1:29631300-29631500   NC_036781.1     100.000 102     0       0       100     201     19925153        19925052        1.53e-46        189
#NC_036781.1:29631300-29631500   NC_036781.1     100.000 100     0       0       1       100     43556556        43556655        1.98e-45        185

#faidx kocher_G_Aulon_yelhead_Male_anchored_assembly.fasta NC_036781.1:5101060-5102938 > YH_G_Inv2_R.fa
#blastn -db ../Mzebra_GT3/Mzebra_GT3.fasta -query YH_G_Inv2_R.fa -outfmt 6 | more                      
#NC_036781.1:5101060-5102938     NC_036781.1     99.138  812     5       2       1       812     43559853        43559044        0.0     1459

#faidx kocher_G_Aulon_yelhead_Male_anchored_assembly.fasta NC_036781.1:5110060-5125938 > YH_G_Inv2_R.fa
#(genomes) ➜  kocher_G_Aulon_yelhead_Male blastn -db ../Mzebra_GT3/Mzebra_GT3.fasta -query YH_G_Inv2_R.fa -outfmt 6 | more                      
#NC_036781.1:5110060-5125938     NC_036781.1     99.309  6079    38      4       6170    12247   19787402        19793477        0.0     10990
#NC_036781.1:5110060-5125938     NC_036781.1     99.531  4904    23      0       7344    12247   20018261        20023164        0.0     8929

# For 9


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