import pandas as pd
import seaborn as sns 
import subprocess, os, allel
from helper_modules.file_manager_Replacement import FileManager as FM
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
import sys

linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
                 'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
                 'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
                 'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

fm_obj = FM(genome_version = 'Mzebra_GT3')

mc_raw_vcf = '/Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/vcf_concat_output/MC_fst/MCAnalysis_master_file.vcf.gz'
mc_vcf = '/Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/vcf_concat_output/MC_fst/MCAnalysisFiltered.vcf.gz'

mc_samples = [
    'Kocher_MC1m', 'Kocher_MC2f', 'MC-008-m', 'MC-1G2G-f', 'MC-1O7O-f', 
    'MC-1R6R-f', 'MC-2B4B-f', 'MC-2G8G-f', 'MC-2P4P-f', 'MC-3P6P-f', 
    'MC-3R6R-f', 'MC-4O5O-f', 'MC-5B6B-f', 'MC-5G11G-f', 'MC_010_f'
]

print("Filtering master VCF for the 15 MC samples and core linkage groups...")
lgs = ','.join(linkageGroups.keys())

# Added Error Handling and --mac 1 to drop sites that are monomorphic in the subset
try:
    subprocess.run([
        'bcftools', 'view',
        '-m', '2', '-M', '2',       
        '-V', 'indels,other',       
        '-s', ','.join(mc_samples), 
        '--min-ac', '1',            # Changed from --mac to --min-ac
        '-r', lgs,                  
        '-o', mc_vcf,
        '-O', 'z',                  
        mc_raw_vcf
    ], check=True, capture_output=True, text=True)
except subprocess.CalledProcessError as e:
    print(f"CRITICAL ERROR: bcftools failed. Details below:\n{e.stderr}")
    sys.exit(1)

print("Filtering complete. Indexing the new filtered VCF...")
subprocess.run(['bcftools', 'index', '-t', mc_vcf])

print("Beginning scikit-allel Analysis...")
vcf_obj = allel.read_vcf(mc_vcf)

# Safety check to ensure variants exist before proceeding
if vcf_obj is None or 'calldata/GT' not in vcf_obj:
    print("CRITICAL ERROR: The filtered VCF contains no variants. Check your bcftools logic.")
    sys.exit(1)

test = allel.GenotypeArray(vcf_obj['calldata/GT'])
vcf_samples = vcf_obj['samples']

window = 100

males = [i for i, sample in enumerate(vcf_samples) if str(sample).endswith('m')]
females = [i for i, sample in enumerate(vcf_samples) if str(sample).endswith('f')]

print(f"Identified {len(males)} males and {len(females)} females in the VCF.")
print("Calculating Fst and Heterozygosity. This may take a moment...")

fst = allel.average_weir_cockerham_fst(test, [males, females], window)
males_test = test.subset(sel1 = males)
females_test = test.subset(sel1 = females)
males_het = allel.heterozygosity_observed(males_test)
females_het = allel.heterozygosity_observed(females_test)

het_dif = males_het - females_het
het_dif_avg = uniform_filter1d(het_dif, window)

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

print("Generating FST Manhattan plot...")
fig = plt.figure(figsize=(14, 8))
ax = fig.add_subplot(111)
colors = ['black','grey']
x_labels = []
x_labels_pos = []
for num, (name, group) in enumerate(dt_grouped):
    group.plot(kind='scatter', x='Offset', y='Fst',color=colors[num % len(colors)], s = 3, ax=ax)
    x_labels.append(name.replace('LG',''))
    x_labels_pos.append((group['Offset'].iloc[-1] - (group['Offset'].iloc[-1] - group['Offset'].iloc[0])/2))
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, dt.Offset.max()])
ax.set_xlabel('Chromosome')
plt.savefig('Fst_MC_Manhattan.pdf')
plt.close()

print("Generating Heterozygosity Manhattan plot...")
fig = plt.figure(figsize=(14, 8))
ax = fig.add_subplot(111)
for num, (name, group) in enumerate(dt_grouped):
    group.plot(kind='scatter', x='Offset', y='Het',color=colors[num % len(colors)], s = 3, ax=ax)
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, dt.Offset.max()])
ax.set_xlabel('Chromosome')
plt.savefig('Het_MC_Manhattan.pdf')
plt.close()

print("Generating targeted LG10 line plots...")
fig,axes = plt.subplots(1, 2, figsize = (14, 6))
sns.lineplot(data=dt, x = 'Pos', y = 'Fst', hue = 'Contig', ax = axes[0], legend = False)
sns.lineplot(data=dt[dt.Contig == 'LG10'], x = 'Pos', y = 'Fst', linewidth = 2.5, color = 'black', ax = axes[0], legend = False)
axes[0].set_ylim([0,.5])
axes[0].set_title('MC Fst')

sns.lineplot(data=dt, x = 'Pos', y = 'Het', hue = 'Contig', ax = axes[1], legend = False)
sns.lineplot(data=dt[dt.Contig == 'LG10'], x = 'Pos', y = 'Het', linewidth = 2.5, color = 'black', ax = axes[1], legend = False)
axes[1].set_ylim([-.5,.5])
axes[1].set_title('MC Heterozygosity Difference')

plt.savefig('SexAnalysis_MC.pdf')
plt.close()

print("Analysis complete. Check your directory for the PDF outputs.")