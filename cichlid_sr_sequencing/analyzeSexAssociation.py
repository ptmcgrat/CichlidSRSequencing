
import pandas as pd
import seaborn as sns 
import subprocess, pdb, os, allel
from helper_modules.file_manager import FileManager as FM
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d

def parallel_filter(input_vcf, output_vcf, samples):
	linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
							  'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
							  'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
							  'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

	base_path = os.path.dirname(output_vcf)
	processes = []
	vcf_files = []

	min_count = int(len(samples)*2/5)

	for lg in linkageGroups.keys():
		lg_vcf = base_path + '/' + lg + '.vcf.gz'
		vcf_files.append(lg_vcf)
		processes.append(subprocess.Popen(['bcftools','view','-m','2','-M','2','-V','indels,other','-s',','.join(samples), '--min-ac',str(min_count) + ':minor', '-r', lg, '-o', lg_vcf, '-O', 'b', input_vcf]))

	for p in processes:
		p.communicate()

	for lg in linkageGroups.keys():
	   subprocess.run(['bcftools','index', base_path + '/' + lg + '.vcf.gz'])

	subprocess.run(['bcftools', 'concat'] + vcf_files + ['-o',output_vcf, '-O','z'])

	for lg in linkageGroups.keys():
		subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz'])
		subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz.csi'])

	subprocess.run(['bcftools','index', output_vcf])


#Dictionary to convert NCBI names of chromosomes to easier to read names
linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
							  'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
							  'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
							  'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

# Object to keep track of filenames and download/upload data from lab Dropbox
fm_obj = FM(genome_version = 'Mzebra_GT3') 

# Download vcf file for our two YH pedigrees
yh_raw_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/YHPedigree/YHPedigree_pass_variants_master_file.vcf.gz' # Download vcf file of pedigree
fm_obj.downloadData(yh_raw_vcf)
fm_obj.downloadData(yh_raw_vcf + '.tbi')

# Download and read in database for each sample (e.g. Sex information)
fm_obj.downloadData(fm_obj.localSampleFile_v2)
s_dt = pd.read_excel(fm_obj.localSampleFile_v2, sheet_name = 'SampleLevel')

# Filter vcf file based on allele counts, snv/indel, etc
yh_samples = s_dt[(s_dt.YHPedigree == 'Yes') & (s_dt.Organism == 'Aulonocara Nhkata Yellow Head')].SampleID.to_list()
yh_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/YHPedigreeFiltered.vcf.gz'
parallel_filter(yh_raw_vcf, yh_vcf, yh_samples)

# VCF file is ready to analyze with scikit-allel 
window = 100
vcf_obj = allel.read_vcf(yh_vcf)
yh_genotype_array = allel.GenotypeArray(vcf_obj['calldata/GT'])

# Get sex information for each sample
sex = [s_dt.loc[s_dt.SampleID == sample, 'Sex'].values[0] for sample in vcf_obj['samples']]
males = [i for i, x in enumerate(sex) if x=='m']
females = [i for i, x in enumerate(sex) if x=='f']

# Calculate Fst for male vs. female samples
fst = allel.average_weir_cockerham_fst(yh_genotype_array, [males,females], window)

# Calculate heterozygosity in male vs. female samples
males_test = test.subset(sel1 = males)
females_test = test.subset(sel1 = females)
males_het = allel.heterozygosity_observed(males_test)
females_het = allel.heterozygosity_observed(females_test)
het_dif = males_het - females_het
het_dif_avg = uniform_filter1d(het_dif, window)


# Create dataframe for plotting
dt = pd.DataFrame(columns = ['Contig','Pos','Fst'])
dt['Contig'] = [linkageGroups[x] for x in vcf_obj['variants/CHROM'][::window][:-1]]
dt['Pos'] = vcf_obj['variants/POS'][::window][:-1]
dt['Fst'] = fst[2]
dt['Het'] = het_dif_avg[::window][:-1]
lengths = dt.groupby('Contig')['Pos'].max().to_dict()
lg_order = list(linkageGroups.values())
total_lengths = {}
for i,lg in enumerate(lg_order):
	total_length = 0
	for j in range(i):
		total_length += lengths[lg_order[j]]
	total_lengths[lg] = total_length
dt['Offset'] = dt.Pos + dt.Contig.map(total_lengths)
dt = dt.sort_values('Offset')
dt_grouped = dt.groupby(('Contig'), sort=False)

# Create figure
fig = plt.figure(figsize=(14, 8)) # Set the figure size
ax = fig.add_subplot(111)
colors = ['black','grey']
x_labels = []
x_labels_pos = []
for num, (name, group) in enumerate(dt_grouped):
	print(num)
	print(name)
	group.plot(kind='scatter', x='Offset', y='Fst',color=colors[num % len(colors)], s = 3, ax=ax)
	x_labels.append(name.replace('LG',''))
	x_labels_pos.append((group['Offset'].iloc[-1] - (group['Offset'].iloc[-1] - group['Offset'].iloc[0])/2))
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)

# set axis limits
ax.set_xlim([0, dt.Offset.max()])
#ax.set_ylim([0, 3.5])

# x axis label
ax.set_xlabel('Chromosome')

# show the graph
plt.savefig('Fst_YH_Manhattan.pdf')

plt.show()
fig = plt.figure(figsize=(14, 8)) # Set the figure size
ax = fig.add_subplot(111)
colors = ['black','grey']
x_labels = []
x_labels_pos = []
for num, (name, group) in enumerate(dt_grouped):
    group.plot(kind='scatter', x='Offset', y='Het',color=colors[num % len(colors)], s = 3, ax=ax)
    x_labels.append(name.replace('LG',''))
    x_labels_pos.append((group['Offset'].iloc[-1] - (group['Offset'].iloc[-1] - group['Offset'].iloc[0])/2))
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)

# set axis limits
ax.set_xlim([0, dt.Offset.max()])
#ax.set_ylim([0, 3.5])

# x axis label
ax.set_xlabel('Chromosome')

# show the graph
plt.savefig('Fst_YH_Het.pdf')

plt.show()

fig,axes = plt.subplots(2,2, figsize = (11,8.5))
sns.lineplot(dt, x = 'Pos', y = 'Fst', hue = 'Contig', ax = axes[0,0], legend = False)
axes[0,0].set_ylim([0,.5])
axes[0,0].set_title('YH_Fst')

sns.lineplot(dt, x = 'Pos', y = 'Het', hue = 'Contig', ax = axes[0,1], legend = False)
axes[0,1].set_ylim([-.5,.5])
axes[0,1].set_title('YH_Het')

plt.savefig('SexAnalysis.pdf')

plt.show()
pdb.set_trace()