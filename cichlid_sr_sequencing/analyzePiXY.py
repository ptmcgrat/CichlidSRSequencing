import pandas as pd
import seaborn as sns 
import subprocess, pdb, os, allel, sys
from helper_modules.file_manager_Replacement import FileManager as FM
from Bio import AlignIO
import matplotlib.pyplot as plt
import scipy.spatial

# Dictionary for the location of each inversion
inversions = {'LG2':('NC_036781.1',19743639,20311160,43254805,43658853),'LG9':('NC_036788.1',14453796,15649299,32255605,33496468),'LG10':('NC_036789.1',11674905,11855817,29898615,29898615),
				'LG11':('NC_036790.1',8302039,8309764,30371888,30459686),'LG13':('NC_036792.1',2249541,2453698,23046928,23131968),
				'LG20a':('NC_036799.1',19614379,19689710,29716509,29673642),'LG20b':('NC_036799.1',29716509,29673642,32872827,33764042)}


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

# Download master vcf file
main_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/gt3_cohort_pass_variants.vcf.gz'
fm_obj.downloadData(main_vcf)
fm_obj.downloadData(main_vcf + '.tbi')

# Download and read in database for each sample (e.g. Ecogroup information)
fm_obj.downloadData(fm_obj.localSampleFile_v2)
s_dt = pd.read_excel(fm_obj.localSampleFile_v2, sheet_name = 'SampleLevel')

# Filter vcf file for lab sample data and for lake malawi phylogeny data (defined in Sample Database)
malawi_samples = s_dt[s_dt.CorePCA == 'Yes'].SampleID.to_list()
s_dt = s_dt[(s_dt.PCAFigure == 'Yes') & (s_dt.BionanoData == 'No')]

malawi_samples = s_dt.SampleID.to_list()

malawi_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals.vcf.gz'

print('Filtering main vcf')
parallel_filter(main_vcf, malawi_vcf, malawi_samples)
subprocess.run(['bcftools','index', '-f', malawi_vcf])

print('Creating vcfs for inversions')

#Create vcf files for each each inversion from the malawi samples vcf file
for lg,(contig,temp,left,right,temp2) in inversions.items():
	region = contig + ':' + str(left) + '-' + str(right)
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	subprocess.run(['bcftools','view','-r', region, '-o', lg_vcf, '-O', 'z', malawi_vcf])
	subprocess.run(['bcftools','index', lg_vcf])


# Use scikit allel to calculate dxy comparing all samples. Store it in a pandas dataframe.
# Do this for the whole genome as well as the indepednent inversions
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
vcf_obj = allel.read_vcf(malawi_vcf)
test = allel.GenotypeArray(vcf_obj['calldata/GT'])
out = allel.pairwise_distance(test.to_n_alt(), metric='cityblock')
sq = scipy.spatial.distance.squareform(out)
size = 950000000
for i in range(len(sq)):
	for j in range(len(sq)):
		sam_i = vcf_obj['samples'][i]
		sam_j = vcf_obj['samples'][j]
		if sam_i not in malawi_samples or sam_j not in malawi_samples:
			continue
		if i == j:
			continue
		if s_dt[s_dt.SampleID == sam_i]['Organism'].values[0] == s_dt[s_dt.SampleID == sam_j]['Organism'].values[0]:
			continue
		try:
			ecogroup_i = s_dt[s_dt.SampleID == sam_i]['Ecogroup_PTM'].values[0]
			ecogroup_j = s_dt[s_dt.SampleID == sam_j]['Ecogroup_PTM'].values[0]
		except:
			pdb.set_trace()
		out_dt.loc[len(out_dt)] = [ecogroup_i,ecogroup_j,sam_i,sam_j,'WholeGenome',sq[i,j]/size]

# Create figures analyzing the sample comparisons based on inversion genotype and ecogroup
print('Calculating genetic distance for each inversion')

LG11_inv = s_dt[(s_dt.LG11 == 'Inverted') & (s_dt.Ecogroup_PTM != 'Diplotaxodon')]['SampleID'].to_list()
LG11_noninv = s_dt[(s_dt.LG11 == 'Normal') & (s_dt.Ecogroup_PTM.isin(['Shallow_Benthic','Deep_Benthic','Utaka']))]['SampleID'].to_list()
#LG11_CV = s_dt[(s_dt.LG11 == 'NormalCV')]['SampleID'].to_list()
LG11_exclude = s_dt[(s_dt.LG11 == 'Heterozygous')]['SampleID'].to_list()

LG9_inv = s_dt[(s_dt.LG9 == 'Inverted') & (s_dt.Ecogroup_PTM.isin(['Shallow_Benthic','Deep_Benthic','Utaka']))]['SampleID'].to_list()
LG9_het = s_dt[(s_dt.LG9 == 'Heterozygous') & (s_dt.Ecogroup_PTM.isin(['Shallow_Benthic','Deep_Benthic','Utaka']))]['SampleID'].to_list()

LG20_Rhamph = s_dt[s_dt.LG11 == 'Inverted']['SampleID'].to_list()

out_dt['LG11Type1'] = out_dt['Ecogroup1']
out_dt.loc[out_dt.Sample1.isin(LG11_inv),'LG11Type1'] = 'LG11Inv'
out_dt.loc[out_dt.Sample1.isin(LG11_noninv),'LG11Type1'] = 'LG11NonInv'
out_dt.loc[out_dt.Sample1.isin(LG11_exclude),'LG11Type1'] = 'LG11Het'

#out_dt.loc[out_dt.Sample1.isin(LG11_CV),'LG11Type1'] = 'LG11CVInv'

out_dt['LG11Type2'] = out_dt['Ecogroup2']
out_dt.loc[out_dt.Sample2.isin(LG11_inv),'LG11Type2'] = 'LG11Inv'
out_dt.loc[out_dt.Sample2.isin(LG11_noninv),'LG11Type2'] = 'LG11NonInv'
out_dt.loc[out_dt.Sample2.isin(LG11_exclude),'LG11Type2'] = 'LG11Het'

#out_dt.loc[out_dt.Sample2.isin(LG11_CV),'LG11Type2'] = 'LG11CVInv'

out_dt['LG9Type1'] = out_dt['Ecogroup1']
out_dt.loc[out_dt.Sample1.isin(LG9_inv),'LG9Type1'] = 'BenthicInv'
out_dt.loc[out_dt.Sample1.isin(LG9_het),'LG9Type1'] = 'BenthicHet'

out_dt['LG9Type2'] = out_dt['Ecogroup2']
out_dt.loc[out_dt.Sample2.isin(LG9_inv),'LG9Type2'] = 'BenthicInv'
out_dt.loc[out_dt.Sample2.isin(LG9_het),'LG9Type2'] = 'BenthicHet'

out_dt.to_csv('AllData.csv')

out_dt.loc[out_dt.Ecogroup1.isin(['Shallow_Benthic','Deep_Benthic','Utaka']),'Ecogroup1'] = 'Benthic/Utaka'
out_dt.loc[out_dt.Ecogroup2.isin(['Shallow_Benthic','Deep_Benthic','Utaka']),'Ecogroup2'] = 'Benthic/Utaka'

pdb.set_trace()

fig,axes = plt.subplots(5,1, figsize = (8.5,11))
fig.subplots_adjust(hspace=0.8)
fig.subplots_adjust(wspace=0.6)

malinksy_color_map = {'Mbuna': '#A020F0', 'AC': '#A2CD5A', 'LG11Inv': '#FF6347', 'BenthicInv': '#FF6347','LG11NonInv': 'grey','Benthic/Utaka': '#FF6347', 'Rhamphochromis': '#8B4513', 'Diplotaxodon': '#FFA54F'}

sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Benthic/Utaka')&(out_dt.Inversion == 'WholeGenome')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[0], binrange = (0,.002), legend = False, bins = 100, palette = malinksy_color_map)
axes[0].set_title('Benthic_WholeGenome')
sns.histplot(out_dt[(out_dt.LG9Type1 == 'BenthicInv')&(out_dt.Inversion == 'LG9')&(out_dt.LG9Type2 != 'BenthicHet') & (out_dt.LG9Type2 != 'Shallow_Benthic')], x = 'Distance', hue = 'LG9Type2', stat = 'proportion', common_norm = False, ax = axes[1], binrange = (0,.002), legend = False, bins = 100,palette = malinksy_color_map)
axes[1].set_title('BenthicInverted_LG9')
sns.histplot(out_dt[(out_dt.LG11Type1 == 'LG11Inv')&(out_dt.Inversion == 'LG11')&(out_dt.LG11Type2 != 'LG11Het')], x = 'Distance', hue = 'LG11Type2', stat = 'proportion', common_norm = False, ax = axes[2], binrange = (0,.002), legend = True, bins = 100,palette = malinksy_color_map)
axes[2].set_title('LG11Inverted_Inv11a')
#sns.histplot(out_dt[(out_dt.LG11Type1 == 'LG11Inv')&(out_dt.Inversion == 'LG11b')], x = 'Distance', hue = 'LG11Type2', stat = 'proportion', common_norm = False, ax = axes[2], binrange = (0,.002), legend = False, bins = 100)
#axes[2].set_title('LG11Inverted_Inv11b')
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Benthic/Utaka')&(out_dt.Inversion == 'LG20a')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[3], binrange = (0,.002), legend = False, bins = 100,palette = malinksy_color_map)
axes[3].set_title('Benthic_Inv20a')
sns.histplot(out_dt[(out_dt.Ecogroup1 == 'Benthic/Utaka')&(out_dt.Inversion == 'LG20b')], x = 'Distance', hue = 'Ecogroup2', stat = 'proportion', common_norm = False, ax = axes[4], binrange = (0,.002), legend = False, bins = 100,palette = malinksy_color_map)
axes[4].set_title('Benthic_Inv20b')
plt.savefig('ForPaper.pdf')

plt.show()

sys.exit()

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
#sns.histplot(out_dt[(out_dt.LG11Type1 == 'LG11CVInv')&(out_dt.Inversion == 'LG11a')], x = 'Distance', hue = 'LG11Type2', stat = 'proportion', common_norm = False, ax = axes[2], binrange = (0,.002), legend = True, bins = 100)
axes[2].set_title('LG11CV_Inv11a')
#sns.histplot(out_dt[(out_dt.LG11Type1 == 'LG11CVInv')&(out_dt.Inversion == 'LG11b')], x = 'Distance', hue = 'LG11Type2', stat = 'proportion', common_norm = False, ax = axes[3], binrange = (0,.002), legend = True, bins = 100)
axes[3].set_title('LG11CV_Inv11b')
plt.savefig('Benthic_Inversion11_2nd.pdf')
plt.show()

#sns.move_legend(axes[2,2], 'upper center', bbox_to_anchor=(-.5, -.5))


		