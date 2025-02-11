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



inversions = {'LG2':('NC_036781.1',19745291,19746383,43556654,43556721),'LG9':('NC_036788.1',14503645,14972797,32763228,32765988),'LG10':('NC_036789.1',11817927,11841666,29881282,29902558),
				'LG11':('NC_036790.1',8304758,8303033,30372595,30385320),'LG13':('NC_036792.1',2398101,2405661,23071814,23124402),'LG20':('NC_036799.1',19616371,19616393,33288065,33254389)}

linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
							  'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
							  'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
							  'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

fm_objs = {}

# Download genome versions
# 2 additional versions of Metriaclima zebra, 3 versions of Aulonocara sp. ‘chitande type north’ Nkhata Bay, 1 version of Rhamphochromis chilingali
# As outgroups, one version of Oreochromis niloticus and one version of Pundamilia nyererei
genome_versions = ['kocher_N_Met_zebra_Female','MZ4f_ptm','kocher_H_Aulon_yelhead_Female',
					'kocher_G_Aulon_yelhead_Male','YH7f_ptm','Rhamp_chilingali',
					'O_niloticus_UMD_NMBU','P_nyererei_v2']
for gv in genome_versions:
	fm_objs[gv] = FM(genome_version = gv)
	fm_objs[gv].downloadData(fm_objs[gv].localHybridScaffoldFile)

# Align each genome to the Mzebra_GT3 reference.
for gv in genome_versions:
	if gv == 'Mzebra_GT3':
		continue
	if gv in ['O_niloticus_UMD_NMBU','P_nyrerei_v2','Rhamp_chilingali']:
		align_genomes('Mzebra_GT3',gv)
	else:
		align_genomes('Mzebra_GT3',gv, hs_flag = True) # hs flag is hybrid scaffolded file that hasn't been anchored
	

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

