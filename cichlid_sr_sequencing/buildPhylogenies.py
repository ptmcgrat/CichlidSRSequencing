import pandas as pd
import seaborn as sns 
import subprocess, pdb, os, ete3
from helper_modules.file_manager import FileManager as FM
from Bio import AlignIO
import matplotlib.pyplot as plt
inversions = {'LG2':('NC_036781.1',19743639,20311160,43254805,43658853),'LG9':('NC_036788.1',14453796,15649299,32255605,33496468),'LG10':('NC_036789.1',11674905,11855817,29898615,29898615),
				'LG11':('NC_036790.1',8302039,8309764,30371888,30459686),'LG13':('NC_036792.1',2249541,2453698,23046928,23131968),'LG20a':('NC_036799.1',19614379,19689710,29716509,29673642),
				'LG20b':('NC_036799.1',29716509,29673642,32872827,33764042)}


def parallel_filter(input_vcf, output_vcf, samples):
	linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
							  'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
							  'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
							  'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

	base_path = os.path.dirname(output_vcf)
	processes = []
	vcf_files = []

	for lg in linkageGroups.keys():
		lg_vcf = base_path + '/' + lg + '.vcf.gz'
		vcf_files.append(lg_vcf)
		processes.append(subprocess.Popen(['bcftools','view','-m','2','-M','2','-V','indels,other','-s',','.join(samples), '--min-ac','1:minor', '-r', lg, '-o', lg_vcf, '-O', 'b', input_vcf]))

	for p in processes:
		p.communicate()

	for lg in linkageGroups.keys():
	   subprocess.run(['bcftools','index', base_path + '/' + lg + '.vcf.gz'])

	subprocess.run(['bcftools', 'concat'] + vcf_files + ['-o',output_vcf, '-O','z'])

	for lg in linkageGroups.keys():
		subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz'])
		subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz.csi'])

	subprocess.run(['bcftools','index', output_vcf])


# Object to keep track of filenames and download/upload data from lab Dropbox
fm_obj = FM(genome_version = 'Mzebra_GT3')

# Download vcf file for our phylogenies
main_vcf = fm_obj.localMasterDir + 'Outputs/MasterFile/GT3Cohort/PhylogenyFigure/phylogenyfigure_pass_variants_master_file.vcf.gz'
fm_obj.downloadData(main_vcf)
fm_obj.downloadData(main_vcf + '.tbi')

# Download and read in database for each sample (e.g. Ecogroup information)
fm_obj.downloadData(fm_obj.localSampleFile_v2)
s_dt = pd.read_excel(fm_obj.localSampleFile_v2, sheet_name = 'SampleLevel', index_col = 0)

# Filter vcf file based on allele counts, snv/indel, etc
malawi_samples = s_dt[(s_dt.PhylogenyFigure == 'Yes')].index.to_list()
malawi_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals.vcf.gz'
parallel_filter(main_vcf, malawi_vcf, malawi_samples)

#Create vcf files for each inversion from the malawi samples vcf file
for lg,(contig,temp,left,right,temp2) in inversions.items():
	region = contig + ':' + str(left) + '-' + str(right)
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	subprocess.run(['bcftools','view','-r', region, '-o', lg_vcf, '-O', 'z', malawi_vcf])
	subprocess.run(['bcftools','index', lg_vcf])

#Create relaxed phylip file from from vcf. Requires vcf2phylip to be installed in your PYTHONPATH
subprocess.run(['vcf2phylip.py','-i', malawi_vcf, '-m', '50','--output-folder',fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/'])
for lg,region in inversions.items():
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	subprocess.run(['vcf2phylip.py','-i',lg_vcf, '-m', '50','--output-folder', fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/'])

# Convert relaxed phylip file to stockholm using BioPYTHON
align = AlignIO.read(malawi_vcf.replace('.vcf.gz','.min50.phy'), 'phylip-relaxed')
AlignIO.write(align,malawi_vcf.replace('.vcf.gz','.min50.stk'),'stockholm')
for lg,region in inversions.items():
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	align = AlignIO.read(lg_vcf.replace('.vcf.gz','.min50.phy'), 'phylip-relaxed')
	AlignIO.write(align,lg_vcf.replace('.vcf.gz','.min50.stk'),'stockholm')

# Create phylogeny using iqtree or quicktree
treefiles = []
iqtree = True

for lg,region in inversions.items():
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	
	if iqtree:
		treefiles.append(lg_vcf.replace('.vcf.gz','.min50.treefile'))
		subprocess.run(['iqtree2', '-s', lg_vcf.replace('.vcf.gz','.min50.phy'),'-nt','24','-mem','96G','-v','--seqtype', 'DNA','-m','GTR+I+G','-B','1000','-o','SAMN15891801,SAMN15685499','--prefix', lg_vcf.replace('.vcf.gz','.min50')])
	else:
		treefiles.append(lg_vcf.replace('.vcf.gz','.min50.tree'))
		subprocess.run(['quicktree', lg_vcf.replace('.vcf.gz','.min50.stk')], stdout = open(lg_vcf.replace('.vcf.gz','.min50.tree'),'w'))
if iqtree:
	treefiles.append(malawi_vcf.replace('.vcf.gz','.min50.treefile'))
	subprocess.run(['iqtree2', '-s', malawi_vcf.replace('.vcf.gz','.min50.phy'),'-nt','24','-mem','96G','-v','--seqtype','DNA','-m','GTR+I+G','-B','1000','-o','SAMN15891801,SAMN15685499','--prefix', malawi_vcf.replace('.vcf.gz','.min50')])
else:
	treefiles.append(malawi_vcf.replace('.vcf.gz','.min50.tree'))
	subprocess.run(['quicktree', malawi_vcf.replace('.vcf.gz','.min50.stk')], stdout = open(malawi_vcf.replace('.vcf.gz','.min50.tree'),'w'))

# Replace the sample ID with names that are more interesting to display on the tree
for tf in treefiles:
	f = open(tf)
	data = f.read()
	for i,row in c_dt.iterrows():
		sample = row.SampleID
		species = row.Organism
		ecogroup = row.Ecogroup_PTM
		small_species = species[0] + '_' + species.split()[1]
		new_text = ecogroup + '_' + small_species + '_' + sample
		data = data.replace(sample,new_text)
	out = open(tf.replace('.tree','.modified.tree'), 'w')
	print(data, file = out, end = '')


# Plot phylogenies using ete3
nst = {}
nst['Shallow_Benthic'] = ete3.NodeStyle()
nst['Shallow_Benthic']["bgcolor"] = "LightSalmon"
nst['Deep_Benthic'] = ete3.NodeStyle()
nst['Deep_Benthic']["bgcolor"] = "LightSteelBlue"
nst['Utaka'] = ete3.NodeStyle()
nst['Utaka']["bgcolor"] = "SeaGreen"
nst['Mbuna'] = ete3.NodeStyle()
nst['Mbuna']["bgcolor"] = "MediumOrchid"
nst['AC'] = ete3.NodeStyle()
nst['AC']["bgcolor"] = "PaleGreen"
nst['Diplotaxodon'] = ete3.NodeStyle()
nst['Diplotaxodon']["bgcolor"] = "Gold"
nst['Rhamphochromis'] = ete3.NodeStyle()
nst['Rhamphochromis']["bgcolor"] = "Sienna"
nst['LakeVictoriaHap'] = ete3.NodeStyle()

treefiles = [malawi_vcf.replace('.vcf.gz','.min50.treefile')]
for lg,region in inversions.items():
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	treefiles.append(lg_vcf.replace('.vcf.gz','.min50.treefile'))

for tf in treefiles:
	tree = ete3.Tree(tf)
	ts = ete3.TreeStyle()
	ts.show_branch_support = True
	ts.scale =  1000 # 120 pixels per branch length unit

	for leaf in tree.get_leaves():
		eg = s_dt.loc[leaf.name,'Ecogroup_PTM']
		species = s_dt.loc[leaf.name,'Organism']
		leaf.name = species.split(' ')[0][0] + '. ' + species.split(' ')[1] 
		leaf.set_style(nst[eg])

	if 'LG' in tf:
		ts.title.add_face(ete3.TextFace(tf.split('/')[-1].split('_')[1].split('.')[0], fsize=20), column=0)
	else:
		ts.title.add_face(ete3.TextFace('WholeGenome', fsize=20), column=0)

	tree.render(tf.replace('.treefile','.pdf'), tree_style=ts, w=1600)

from pypdf import PdfMerger

pdfs = [x.replace('.treefile','.pdf') for x in treefiles]

merger = PdfMerger()

for pdf in pdfs:
    merger.append(pdf)

merger.write(fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_AllPhylogenies.pdf')
merger.close()
