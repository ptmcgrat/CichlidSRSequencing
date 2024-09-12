import pandas as pd
import seaborn as sns 
import subprocess, pdb, os, allel
from helper_modules.file_manager_Replacement import FileManager as FM
from Bio import AlignIO
import matplotlib.pyplot as plt
import numpy as np

inversions = {'LG2':('NC_036781.1',19743639,20311160,43254805,43658853),'LG9':('NC_036788.1',14453796,15649299,32255605,33496468),'LG10':('NC_036789.1',11674905,11855817,29898615,29898615),
				'LG11a':('NC_036790.1',8302039,8309764,18425216,18425216),'LG11b':('NC_036790.1',18425216,18425216,30371888,30459686),'LG13':('NC_036792.1',2249541,2453698,23046928,23131968),
				'LG20a':('NC_036799.1',19614379,19689710,29716509,29673642),'LG20b':('NC_036799.1',29716509,29673642,32872827,33764042)}

linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
							  'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
							  'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
							  'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}


def inversion_location(row):
	out = 'Normal'
	for (inv,loc) in inversions.items():
		if row['Contig'] == linkageGroups[loc[0]] and row['Pos'] > loc[1] and row['Pos'] < loc[4]:
			out = inv
	return out

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
		processes.append(subprocess.Popen(['bcftools','view','-m','2','-M','2','-V','indels,other','-s',','.join(samples), '--min-ac','5:minor', '-r', lg, '-o', lg_vcf, '-O', 'b', input_vcf]))

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
#malawi_samples = s_dt[s_dt.CorePCA == 'Yes'].SampleID.to_list()

s_dt = s_dt[(s_dt.PCAFigure == 'Yes') & (s_dt.BionanoData == 'No')]
malawi_samples = s_dt.SampleID.to_list()
#c_dt = s_dt[s_dt.PCAFigure == 'Yes']
malawi_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals.vcf.gz'

print('Filtering main vcf')
#parallel_filter(main_vcf, malawi_vcf, malawi_samples)
#subprocess.run(['bcftools','index', '-f', malawi_vcf])
"""
for lg,(contig,temp,left,right,temp2) in inversions.items():
	region = contig + ':' + str(left) + '-' + str(right)
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	subprocess.run(['bcftools','view','-r', region, '-o', lg_vcf, '-O', 'z', malawi_vcf])
	subprocess.run(['bcftools','index', lg_vcf])
"""
vcf_obj = allel.read_vcf(malawi_vcf)
test = allel.GenotypeArray(vcf_obj['calldata/GT'])




rhampho = [vcf_obj['samples'].tolist().index(x) for x in s_dt[s_dt.Ecogroup_PTM == 'Rhamphochromis'].SampleID.to_list()]
diplo = [vcf_obj['samples'].tolist().index(x) for x in s_dt[s_dt.Ecogroup_PTM == 'Diplotaxodon'].SampleID.to_list()]
pelagic = rhampho + diplo
acs = [vcf_obj['samples'].tolist().index(x) for x in s_dt[s_dt.Ecogroup_PTM == 'AC'].SampleID.to_list()]

deep_ben = [vcf_obj['samples'].tolist().index(x) for x in s_dt[s_dt.LG10 == 'Inverted'].SampleID.to_list()]
shallow_ben = [vcf_obj['samples'].tolist().index(x) for x in s_dt[(s_dt.LG2 == 'Normal') & (s_dt.Ecogroup_PTM.isin(['Shallow_Benthic','Deep_Benthic','Utaka']))].SampleID.to_list()]

mbuna = [vcf_obj['samples'].tolist().index(x) for x in s_dt[s_dt.Ecogroup_PTM == 'Mbuna'].SampleID.to_list()]
benthics = [vcf_obj['samples'].tolist().index(x) for x in s_dt[s_dt.Ecogroup_PTM.isin(['Shallow_Benthic','Deep_Benthic','Utaka'])].SampleID.to_list()]



ac_r = test.count_alleles(subpop=rhampho)
ac_d =  test.count_alleles(subpop=diplo)
ac_p =  test.count_alleles(subpop=pelagic)
ac_ac =  test.count_alleles(subpop=acs)
ac_db =  test.count_alleles(subpop=deep_ben)
ac_sb =  test.count_alleles(subpop=shallow_ben)
ac_mb =  test.count_alleles(subpop=mbuna)
ac_b =  test.count_alleles(subpop=benthics)



num, den = allel.hudson_fst(ac_p, ac_ac)
fst_pelagic = num/den
num, den = allel.hudson_fst(ac_d, ac_r)
fst_diplo = num/den
num, den = allel.hudson_fst(ac_db, ac_sb)
fst_deep = num/den
num, den = allel.hudson_fst(ac_mb, ac_b)
fst_rs = num/den


dt = pd.DataFrame(columns = ['Contig','Pos'])
dt['Contig'] = [linkageGroups[x] for x in vcf_obj['variants/CHROM'][:]]
dt['Pos'] = vcf_obj['variants/POS'][:]
dt['Fst_pel'] = fst_pelagic
dt['Fst_diplo'] = fst_diplo
dt['Fst_deep'] = fst_deep
dt['Fst_rs'] = fst_rs
dt['Inversion'] = dt.apply(inversion_location, axis=1)

print(dt.groupby('Inversion')['Fst_pel'].mean())
print(dt.groupby('Inversion')['Fst_diplo'].mean())
print(dt.groupby('Inversion')['Fst_deep'].mean())
print(dt.groupby('Inversion')['Fst_rs'].mean())

sns.histplot(dt, x = 'Fst_pel', hue = 'Inversion', stat="density", common_norm = False, binwidth = .1, multiple="dodge", shrink = .8)
plt.show()
sns.histplot(dt, x = 'Fst_diplo', hue = 'Inversion', stat="density", common_norm = False, binwidth = .1, multiple="dodge", shrink = .8)
plt.show()
sns.histplot(dt, x = 'Fst_deep', hue = 'Inversion', stat="density", common_norm = False, binwidth = .1, multiple="dodge", shrink = .8)
plt.show()
sns.histplot(dt, x = 'Fst_rs', hue = 'Inversion', stat="density", common_norm = False, binwidth = .1, multiple="dodge", shrink = .8)
plt.show()

pdb.set_trace()
sns.lineplot(dt, x = 'Pos', y = 'Fst', hue = 'Contig', legend = False)
plt.show()



num, den = allel.hudson_fst(ac_p, ac_ac)
fst = num/den

vcf_obj = allel.read_vcf(fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + 'LG20a' + '.vcf.gz')
test = allel.GenotypeArray(vcf_obj['calldata/GT'])

ac_r = test.count_alleles(subpop=rhampho)
ac_d =  test.count_alleles(subpop=diplo)
ac_p =  test.count_alleles(subpop=pelagic)
ac_ac =  test.count_alleles(subpop=acs)

num, den = allel.hudson_fst(ac_p, ac_ac)
fst2 = num/den




"""

print('Creating vcfs for inversions')
#Create lg subsets for each inversion from the malawi samples vcf file



import scipy.spatial
out_dt = pd.DataFrame(columns = ['Ecogroup1','Ecogroup2','Sample1','Sample2','Inversion','Distance'])

print('Calculating genetic distance for each inversion')
for lg,(contig,temp,left,right,temp2) in inversions.items():
	print('lg = ' + lg)
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	
	size = right-left+1

	vcf_obj = allel.read_vcf(lg_vcf)
	

	test = allel.GenotypeArray(vcf_obj['calldata/GT'])
	out = allel.pairwise_distance(test.to_n_alt(), metric='cityblock')
	sq = scipy.spatial.distance.squareform(out)
	for i in range(len(sq)):
		for j in range(len(sq)):
			sam_i = vcf_obj['samples'][i]
			sam_j = vcf_obj['samples'][j]
			
			if sam_i not in malawi_samples or sam_j not in malawi_samples:
				continue
			if i == j:
				continue
			try:
				ecogroup_i = s_dt[s_dt.SampleID == sam_i]['Ecogroup_PTM'].values[0]
				ecogroup_j = s_dt[s_dt.SampleID == sam_j]['Ecogroup_PTM'].values[0]
			except:
				pdb.set_trace()
			out_dt.loc[len(out_dt)] = [ecogroup_i,ecogroup_j,sam_i,sam_j,lg,sq[i,j]/size]

print('Calculating genetic distance for total genome')
"""
vcf_obj = allel.read_vcf(malawi_vcf)
test = allel.GenotypeArray(vcf_obj['calldata/GT'])

window = 100

sex = [s_dt.loc[s_dt.SampleID == sample, 'Sex'].values[0] for sample in vcf_obj['samples']]
males = [i for i, x in enumerate(sex) if x=='m']
females = [i for i, x in enumerate(sex) if x=='f']
fst = allel.average_weir_cockerham_fst(test, [males,females], window)



for i in range(len(sq)):
	for j in range(len(sq)):
		sam_i = vcf_obj['samples'][i]
		sam_j = vcf_obj['samples'][j]
		if sam_i not in malawi_samples or sam_j not in malawi_samples:
			continue
		if i == j:
			continue
		try:
			ecogroup_i = s_dt[s_dt.SampleID == sam_i]['Ecogroup_PTM'].values[0]
			ecogroup_j = s_dt[s_dt.SampleID == sam_j]['Ecogroup_PTM'].values[0]
		except:
			pdb.set_trace()
		out_dt.loc[len(out_dt)] = [ecogroup_i,ecogroup_j,sam_i,sam_j,'WholeGenome',sq[i,j]/size]
"""
try:
	out_dt
except:
	out_dt = pd.read_csv('AllData.csv', index_col = 0)
	out_dt = out_dt[(out_dt.Sample1 != 'SAMEA4032048') & (out_dt.Sample2 != 'SAMEA4032048')]
print('Calculating genetic distance for each inversion')

LG11_inv = s_dt[(s_dt.LG11 == 'Inverted') & (s_dt.Ecogroup_PTM != 'Diplotaxodon')]['SampleID'].to_list()
LG11_noninv = s_dt[(s_dt.LG11 == 'Normal') & (s_dt.Ecogroup_PTM.isin(['Shallow_Benthic','Deep_Benthic','Utaka']))]['SampleID'].to_list()
LG11_CV = s_dt[(s_dt.LG11 == 'NormalCV')]['SampleID'].to_list()
LG11_exclude = s_dt[(s_dt.LG11 == 'Heterozygous')]['SampleID'].to_list()

LG9_inv = s_dt[(s_dt.LG9 == 'Inverted') & (s_dt.Ecogroup_PTM.isin(['Shallow_Benthic','Deep_Benthic','Utaka']))]['SampleID'].to_list()
LG9_het = s_dt[(s_dt.LG9 == 'Heterozygous') & (s_dt.Ecogroup_PTM.isin(['Shallow_Benthic','Deep_Benthic','Utaka']))]['SampleID'].to_list()

LG20_Rhamph = s_dt[s_dt.LG11 == 'Inverted']['SampleID'].to_list()

out_dt['LG11Type1'] = out_dt['Ecogroup1']
out_dt.loc[out_dt.Sample1.isin(LG11_inv),'LG11Type1'] = 'LG11Inv'
out_dt.loc[out_dt.Sample1.isin(LG11_noninv),'LG11Type1'] = 'LG11NonInv'
out_dt.loc[out_dt.Sample1.isin(LG11_CV),'LG11Type1'] = 'LG11CVInv'

out_dt['LG11Type2'] = out_dt['Ecogroup2']
out_dt.loc[out_dt.Sample2.isin(LG11_inv),'LG11Type2'] = 'LG11Inv'
out_dt.loc[out_dt.Sample2.isin(LG11_noninv),'LG11Type2'] = 'LG11NonInv'
out_dt.loc[out_dt.Sample2.isin(LG11_CV),'LG11Type2'] = 'LG11CVInv'

out_dt['LG9Type1'] = out_dt['Ecogroup1']
out_dt.loc[out_dt.Sample1.isin(LG9_inv),'LG9Type1'] = 'BenthicInv'
out_dt.loc[out_dt.Sample1.isin(LG9_het),'LG9Type1'] = 'BenthicHet'

out_dt['LG9Type2'] = out_dt['Ecogroup2']
out_dt.loc[out_dt.Sample2.isin(LG9_inv),'LG9Type2'] = 'BenthicInv'
out_dt.loc[out_dt.Sample2.isin(LG9_het),'LG9Type2'] = 'BenthicHet'

out_dt.to_csv('AllData.csv')

out_dt.loc[out_dt.Ecogroup2.isin(['Shallow_Benthic','Deep_Benthic','Utaka']),'Ecogroup2'] = 'Benthic/Utaka'

fig,axes = plt.subplots(3,1, figsize = (11,8.5))
fig.subplots_adjust(hspace=0.6)
fig.subplots_adjust(wspace=0.6)

sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Rhamphochromis')&(out_dt.Inversion == 'LG20a')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[0], binrange = (0,.0015), legend = True, bins = 100)
axes[0].set_title('Rhampho_Inv20a')
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Rhamphochromis')&(out_dt.Inversion == 'LG20b')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[1], binrange = (0,.0015), legend = False, bins = 100)
axes[1].set_title('Rhampho_Inv20b')
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Rhamphochromis')&(out_dt.Inversion == 'WholeGenome')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[2], binrange = (0,.0015), legend = False, bins = 100)
axes[2].set_title('Rhampho_WholeGenome')
plt.savefig('Rhampho_20.pdf')
plt.show()

fig,axes = plt.subplots(3,1, figsize = (11,8.5))
fig.subplots_adjust(hspace=0.6)
fig.subplots_adjust(wspace=0.6)

sns.histplot(out_dt[(out_dt.Ecogroup1 == 'AC')&(out_dt.Inversion == 'LG20a')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[0], binrange = (0,.0015), legend = True, bins = 100)
axes[0].set_title('AC_Inv20a')
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'AC')&(out_dt.Inversion == 'LG20b')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[1], binrange = (0,.0015), legend = False, bins = 100)
axes[1].set_title('AC_Inv20b')
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'AC')&(out_dt.Inversion == 'WholeGenome')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[2], binrange = (0,.0015), legend = False, bins = 100)
axes[2].set_title('AC_WholeGenome')
plt.savefig('AC_20.pdf')

plt.show()

fig,axes = plt.subplots(3,1, figsize = (11,8.5))
fig.subplots_adjust(hspace=0.6)
fig.subplots_adjust(wspace=0.6)

sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Diplotaxodon')&(out_dt.Inversion == 'LG20a')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[0], binrange = (0,.0015), legend = True, bins = 100)
axes[0].set_title('Diplo_Inv20a')
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Diplotaxodon')&(out_dt.Inversion == 'LG20b')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[1], binrange = (0,.0015), legend = False, bins = 100)
axes[1].set_title('Diplo_Inv20b')
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Diplotaxodon')&(out_dt.Inversion == 'WholeGenome')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[2], binrange = (0,.0015), legend = False, bins = 100)
axes[2].set_title('Diplo_WholeGenome')
plt.savefig('Diplo_20.pdf')
plt.show()

fig,axes = plt.subplots(3,1, figsize = (11,8.5))
fig.subplots_adjust(hspace=0.6)
fig.subplots_adjust(wspace=0.6)

sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Mbuna')&(out_dt.Inversion == 'LG20a')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[0], binrange = (0,.0015), legend = True, bins = 100)
axes[0].set_title('Mbuna_Inv20a')
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Mbuna')&(out_dt.Inversion == 'LG20b')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[1], binrange = (0,.0015), legend = False, bins = 100)
axes[1].set_title('Mbuna_Inv20b')
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Mbuna')&(out_dt.Inversion == 'WholeGenome')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[2], binrange = (0,.0015), legend = False, bins = 100)
axes[2].set_title('Mbuna_WholeGenome')
plt.savefig('Diplo_20.pdf')
plt.show()

fig,axes = plt.subplots(3,1, figsize = (11,8.5))
fig.subplots_adjust(hspace=0.6)
fig.subplots_adjust(wspace=0.6)

out_dt.loc[out_dt.Ecogroup1.isin(['Shallow_Benthic','Deep_Benthic','Utaka']),'Ecogroup1'] = 'Benthic/Utaka'

sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Benthic/Utaka')&(out_dt.Inversion == 'LG20a')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[0], binrange = (0,.0015), legend = True, bins = 100)
axes[0].set_title('Benthic_Inv20a')
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Benthic/Utaka')&(out_dt.Inversion == 'LG20b')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[1], binrange = (0,.0015), legend = False, bins = 100)
axes[1].set_title('Benthic_Inv20b')
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Benthic/Utaka')&(out_dt.Inversion == 'WholeGenome')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[2], binrange = (0,.0015), legend = False, bins = 100)
axes[2].set_title('Benthic_WholeGenome')
plt.savefig('Benthic_20.pdf')
plt.show()


fig,axes = plt.subplots(2,1, figsize = (8.5,11))
fig.subplots_adjust(hspace=0.6)
fig.subplots_adjust(wspace=0.6)


sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Benthic/Utaka')&(out_dt.Inversion == 'WholeGenome')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[0], binrange = (0,.002), legend = True, bins = 100)
axes[0].set_title('Benthic_WholeGenome')
sns.histplot(out_dt[(out_dt.LG9Type1 == 'BenthicInv')&(out_dt.Inversion == 'LG9')&(out_dt.LG9Type2 != 'BenthicHet')], x = 'Distance', hue = 'LG9Type2', stat = 'proportion', common_norm = False, ax = axes[1], binrange = (0,.002), legend = True, bins = 100)
axes[1].set_title('BenthicInverted_LG9')
plt.savefig('Inversion_9.pdf')
plt.show()


fig,axes = plt.subplots(3,1, figsize = (8.5,11))
fig.subplots_adjust(hspace=0.6)
fig.subplots_adjust(wspace=0.6)

out_dt = out_dt[(out_dt.LG11Type1 != 'Shallow_Benthic') & (out_dt.LG11Type2 != 'Shallow_Benthic')]
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Benthic/Utaka')&(out_dt.Inversion == 'WholeGenome')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[0], binrange = (0,.002), legend = True, bins = 100)
axes[0].set_title('Benthic_WholeGenome')
sns.histplot(out_dt[(out_dt.LG11Type1 == 'LG11Inv')&(out_dt.Inversion == 'LG11a')], x = 'Distance', hue = 'LG11Type2', stat = 'proportion', common_norm = False, ax = axes[1], binrange = (0,.002), legend = True, bins = 100)
axes[1].set_title('LG11Inverted_Inv11a')
sns.histplot(out_dt[(out_dt.LG11Type1 == 'LG11Inv')&(out_dt.Inversion == 'LG11b')], x = 'Distance', hue = 'LG11Type2', stat = 'proportion', common_norm = False, ax = axes[2], binrange = (0,.002), legend = False, bins = 100)
axes[2].set_title('LG11Inverted_Inv11b')

plt.savefig('Benthic_Inversion11_1st.pdf')

plt.show()

fig,axes = plt.subplots(4,1, figsize = (8.5,11))
fig.subplots_adjust(hspace=0.6)
fig.subplots_adjust(wspace=0.6)

sns.histplot(out_dt[(out_dt.LG11Type1 == 'LG11NonInv')&(out_dt.Inversion == 'LG11a')], x = 'Distance', hue = 'LG11Type2', stat = 'proportion', common_norm = False, ax = axes[0], binrange = (0,.002), legend = True, bins = 100)
axes[0].set_title('LG11Ancestral_Inv11a')
sns.histplot(out_dt[(out_dt.LG11Type1 == 'LG11NonInv')&(out_dt.Inversion == 'LG11b')], x = 'Distance', hue = 'LG11Type2', stat = 'proportion', common_norm = False, ax = axes[1], binrange = (0,.002), legend = True, bins = 100)
axes[1].set_title('LG11Ancestral_Inv11b')
sns.histplot(out_dt[(out_dt.LG11Type1 == 'LG11CVInv')&(out_dt.Inversion == 'LG11a')], x = 'Distance', hue = 'LG11Type2', stat = 'proportion', common_norm = False, ax = axes[2], binrange = (0,.002), legend = True, bins = 100)
axes[2].set_title('LG11CV_Inv11a')
sns.histplot(out_dt[(out_dt.LG11Type1 == 'LG11CVInv')&(out_dt.Inversion == 'LG11b')], x = 'Distance', hue = 'LG11Type2', stat = 'proportion', common_norm = False, ax = axes[3], binrange = (0,.002), legend = True, bins = 100)
axes[3].set_title('LG11CV_Inv11b')
plt.savefig('Benthic_Inversion11_2nd.pdf')
plt.show()

#sns.move_legend(axes[2,2], 'upper center', bbox_to_anchor=(-.5, -.5))

"""
		