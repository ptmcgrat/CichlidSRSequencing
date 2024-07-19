import pandas as pd
import seaborn as sns 
import subprocess, pdb, os
from helper_modules.file_manager_Replacement import FileManager as FM
from Bio import AlignIO
import matplotlib.pyplot as plt
inversions = {'LG2':('NC_036781.1',19743639,20311160,43254805,43658853),'LG9':('NC_036789.1',14453796,15649299,32255605,33496468),'LG10':('NC_036790.1',11674905,11855817,29898615,29898615),
				'LG11':('NC_036791.1',8302039,8309764,30371888,30459686),'LG13':('NC_036792.1',2249541,2453698,23046928,23131968),'LG20':('NC_036799.1',19614379,19689710,32872827,33764042)}


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
		processes.append(subprocess.Popen(['bcftools','view','-V','indels,other','-s',','.join(samples), '--min-ac','1:minor', '-r', lg, '-o', lg_vcf, '-O', 'b', input_vcf]))

	for p in processes:
		p.communicate()

	for lg in linkageGroups.keys():
	   subprocess.run(['bcftools','index', base_path + '/' + lg + '.vcf.gz'])

	subprocess.run(['bcftools', 'concat'] + vcf_files + ['-o',output_vcf, '-O','z'])

	for lg in linkageGroups.keys():
		subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz'])
		subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz.csi'])

	subprocess.run(['bcftools','index', output_vcf])


# Read and download essential data
fm_obj = FM(genome_version = 'Mzebra_GT3')
main_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/gt3_cohort_pass_variants.vcf.gz'
#fm_obj.downloadData(main_vcf)
fm_obj.downloadData(main_vcf + '.tbi')
fm_obj.downloadData(fm_obj.localSampleFile_v2)
s_dt = pd.read_excel(fm_obj.localSampleFile_v2, sheet_name = 'SampleLevel')

#############################################################
# 1. Create phylogenies for whole genome and each inversion #
#############################################################
# Create vcf file for lab sample data and for lake malawi phylogeny data (defined in Sample Database)
malawi_samples = s_dt[s_dt.CorePCA == 'Yes'].SampleID.to_list()
c_dt = s_dt[s_dt.CorePCA == 'Yes']
malawi_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals.vcf.gz'
"""
#parallel_filter(main_vcf, malawi_vcf, malawi_samples)
subprocess.run(['bcftools','index', malawi_vcf])


#Create lg subsets for each inversion from the malawi samples vcf file
for lg,(contig,temp,left,right,temp2) in inversions.items():
	region = contig + ':' + str(left) + '-' + str(right)
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	subprocess.run(['bcftools','view','-r', region, '-o', lg_vcf, '-O', 'z', malawi_vcf])
	subprocess.run(['bcftools','index', lg_vcf])

#Create relaxed phylip file from from vcf. Requires vcf2phylip to be installed in your PYTHONPATH
subprocess.run(['vcf2phylip.py','-i', malawi_vcf, '-m', '150','--output-folder',fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/'])

for lg,region in inversions.items():
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	subprocess.run(['vcf2phylip.py','-i',lg_vcf, '-m', '150','--output-folder', fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/'])

# Convert relaxed phylip file to stockholm using BioPYTHON
align = AlignIO.read(malawi_vcf.replace('.vcf.gz','.min150.phy'), 'phylip-relaxed')
AlignIO.write(align,malawi_vcf.replace('.vcf.gz','.min150.stk'),'stockholm')

for lg,region in inversions.items():
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	align = AlignIO.read(lg_vcf.replace('.vcf.gz','.min150.phy'), 'phylip-relaxed')
	AlignIO.write(align,lg_vcf.replace('.vcf.gz','.min150.stk'),'stockholm')

# Run quicktree to create phylogeny file
treefiles = []
treefiles.append(malawi_vcf.replace('.vcf.gz','.min150.tree'))
subprocess.run(['quicktree', malawi_vcf.replace('.vcf.gz','.min150.stk')], stdout = open(malawi_vcf.replace('.vcf.gz','.min150.tree'),'w'))
for lg,region in inversions.items():
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	subprocess.run(['quicktree', lg_vcf.replace('.vcf.gz','.min150.stk')], stdout = open(lg_vcf.replace('.vcf.gz','.min150.tree'),'w'))
	treefiles.append(lg_vcf.replace('.vcf.gz','.min150.tree'))
	


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
"""
#subprocess.run(['iqtree2', '-s', malawi_vcf.replace('.vcf.gz','.min150.phy'),'-nt','AUTO','-alrt','1000','-B','1000'])
subprocess.run(['iqtree2', '-s', malawi_vcf.replace('.vcf.gz','.min150.phy'),'-nt','24', '-v'])

###########################################################
# 2. Analyze Yellowhead individuals for pedigree analysis #
###########################################################
yellowhead_samples = list(set(s_dt[s_dt.Organism == 'Aulonocara Nhkata Yellow Head'].SampleID.to_list()))
yellowhead_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/YellowHeadIndividuals.vcf.gz'
#parallel_filter(main_vcf, yellowhead_vcf, yellowhead_samples)

sample_name = ['YH_1_m', 'YH_002_m_', 'YH_003_m', 'YH_004_m', 'YH_005_m', 'YH_006_f', 'YH_007_f', 'YH_008_m', 'YH_009_m', 'YH_010_f', 'YH_011_m', 'YH_013_f', 'YH_016', 'YH_017', 'YH_018', 'YH_019', 'YH_020', 'YH_021', 'YH_022', 'YH_023', 'YH_024', 'YH_025', 'YH_026', 'YH_027', 'YH_028', 'YH_029', 'YH_030', 'YH_031', 'YH_032', 'YH_033', 'YH_037', 'YH_038', 'YH_039', 'YH_040', 'YH_041', 'YH_042', 'YH_046', 'YH_047', 'YH_048', 'YH_049', 'YH_050', 'YH_051', 'YH_052', 'YH_053', 'YH_054', 'YH_055', 'YH_056', 'YH_057', 'YH_058', 'YH_059', 'YH_060', 'YH_061']
grouping = ['P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P13','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood2','Brood1','Brood1','Brood1','Brood1','Brood1','Brood1','Brood1','Brood1','Brood1','Brood1','Brood1','Brood1','Brood1','Brood1','Brood1','Brood1']
sample_to_brood = pd.DataFrame({'SampleName':sample_name,'Grouping':grouping})
#subprocess.call(['vcftools','--gzvcf',yellowhead_vcf,'--relatedness2','--out',yellowhead_vcf.replace('.vcf.gz','')])
#vcftools --vcf PATH/Example.vcf \ --relatedness2 \ --out PATH/Example
#subprocess.run(['plink2','--vcf', main_vcf,'--make-pgen','--double-id','--allow-extra-chr','--out',main_vcf.replace('.vcf.gz',''), '--vcf-min-gq', '5'])
#subprocess.run(['plink2', '-pfile', main_vcf.replace('.vcf.gz',''), '--allow-extra-chr', '--make-king-table', '--out', main_vcf.replace('.vcf.gz','')])
#subprocess.run(['plink2','-pfile',hf_vcf.replace('.vcf.gz',''),'--freq','cols=+pos','--allow-extra-chr','-out',fm_obj.localMasterDir + 'AllSamples_af'])
dt = pd.read_csv(main_vcf.replace('.vcf.gz','.kin0'), sep = '\t')
hm = dt.pivot(columns = 'IID1', index = 'IID2', values = 'KINSHIP')
for sample1 in hm.columns:
	hm.loc[sample1,sample1] = 0.5
	for sample2 in hm.columns:
		if hm.loc[sample1,sample2] != hm.loc[sample1,sample2]:
			hm.loc[sample1,sample2] = hm.loc[sample2,sample1]

yellowhead_data = hm[hm.index.str.contains('YH')][hm.columns[hm.columns.str.contains('YH')]]
#pdb.set_trace()
#plink --vcf PATH/Example.vcf \ --allow-extra-chr \ --make-bed \ --out PATH/Example
#plink2 --bfile PATH/Example \ --allow-extra-chr \ --make-king-table \ --out PATH/Example


