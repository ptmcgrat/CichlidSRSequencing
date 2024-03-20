import pandas as pd 
import subprocess, pdb, os, vcf
from helper_modules.file_manager import FileManager as FM
from Bio import AlignIO
from io import StringIO
import pandas as pd 
import seaborn as sns
#mc_males = ['1B11','1B18','1B4','1B5','1C11','1C17','1C4','1C5','2B10','2B13','2B17','2B19','2C10','2C14','2C18',
#           '2C19','3B2','3B6','3B7','3B9','3C2','3C6','3C7','3C9','4B12','4B14','4B25','4C12','4C13','4C25',
#           '5B22','5B23','5B24','5B26','5C22','5C23','5C24','5C26','MC_1_m','MC-008-m']
mc_males = ['1B11','1B18','1B4','1B5','1C11','1C17','1C4','1C5','3B2','3B6','3B7','3B9','3C2','3C6','3C7','4B12','4B14','4B25','4C12','4C13','4C25',
            '5B22','5B23','5B24','5B26','5C22','5C23','5C24','5C26','MC_1_m','MC-008-m']

mc_females = ['MC_2_f','MC-002-f','MC-006-f','MC-007-f','MC-009-f','MC-1G2G-f','MC-1O7O-f','MC-1R6R-f','MC-2B4B-f',
            'MC-2G8G-f','MC-2P4P-f','MC-3P6P-f','MC-3R6R-f','MC-4O5O-f','MC-5B6B-f','MC-5G11G-f']

inversions = {'LG2':'NC_036781.1:14000000-28000000','LG9':'NC_036788.1:5300000-15500000', 'LG10':'NC_036789.1:10400000-28500000','LG11':'NC_036790.1:6000001-25200000',
                'LG13':'NC_036792.1:4000000-21700000','LG20':'NC_036799.1:13900000-23200000'}

lengths = {'NC_036780.1': 38676823, 'NC_036781.1': 32660920, 'NC_036782.1': 37314939, 'NC_036783.1': 30518969, 'NC_036784.1': 36170306, 'NC_036785.1': 39774503, 
            'NC_036786.1': 64916660, 'NC_036787.1': 23971387, 'NC_036788.1': 21023328, 'NC_036789.1': 32356376, 'NC_036790.1': 32446080, 'NC_036791.1': 34086040, 
            'NC_036792.1': 32072427, 'NC_036793.1': 37870038, 'NC_036794.1': 34548429, 'NC_036795.1': 34742138, 'NC_036796.1': 35776103, 'NC_036797.1': 29505383, 
            'NC_036798.1': 25963618, 'NC_036799.1': 29789237, 'NC_036800.1': 34725690, 'NC_036801.1': 42088218}


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

    #for lg in linkageGroups.keys():
    #   subprocess.run(['bcftools','index', base_path + '/' + lg + '.vcf.gz'])

    subprocess.run(['bcftools', 'concat'] + vcf_files + ['-o',output_vcf, '-O','z'])

    for lg in linkageGroups.keys():
        subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz'])

def parallel_filter_af(input_vcf, output_vcf, samples):
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
        processes.append(subprocess.Popen(['bcftools','view','-V','indels,other','-s',','.join(samples),'-r', lg, '-o', lg_vcf, '-O', 'b', input_vcf]))

    for p in processes:
        p.communicate()

    #for lg in linkageGroups.keys():
    #   subprocess.run(['bcftools','index', base_path + '/' + lg + '.vcf.gz'])

    subprocess.run(['bcftools', 'concat'] + vcf_files + ['-o',output_vcf, '-O','z'])

    for lg in linkageGroups.keys():
        subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz'])

def parallel_filter_population_af(input_vcf, output_vcf, samples):
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
        processes.append(subprocess.Popen(['bcftools','view','-V','indels,other','-m','2','-M','2','-r', lg,'-s',','.join(samples), '--min-ac','20:minor', '-o', lg_vcf, '-O', 'b', input_vcf]))

    for p in processes:
        p.communicate()

    #for lg in linkageGroups.keys():
    #   subprocess.run(['bcftools','index', base_path + '/' + lg + '.vcf.gz'])

    subprocess.run(['bcftools', 'concat'] + vcf_files + ['-o',output_vcf, '-O','z'])

    for lg in linkageGroups.keys():
        subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz'])


def make_inversion_genefile(fm_obj):
    annotation_file = fm_obj.localGenomeDir + 'GCF_000238955.4_M_zebra_UMD2a_genomic.gtf'
    ortholog_file = fm_obj.localGenomeDir + 'gene_list_mzebra_v2ref.csv'
    gene_info_file = fm_obj.localGenomeDir + 'MzUMD2a.geneinfo.txt'

    combined_gene_file = fm_obj.localGenomeDir + 'tree_fam_orthologs.xlsx'

    fm_obj.downloadData(annotation_file)
    fm_obj.downloadData(ortholog_file)
    fm_obj.downloadData(gene_info_file)

    o_dt = pd.read_csv(ortholog_file)
    o_dt = o_dt.groupby('mzebra').first()[['human','mzebra_description','human_description','ens']].reset_index()
    gi_dt = pd.read_csv(gene_info_file, sep = '\t')
    g_dt = pd.DataFrame(columns = ['Chr', 'Start','Stop','Strand','GeneName','Type'])
    

    with open(annotation_file) as f:
        for line in f:
            if line[0] == '#':
                continue
            chrom, temp, level, start, stop, temp, strand, temp, info = line.rstrip().split('\t')
            if level != 'gene':
                continue
            gene_name = info.split('"')[1]
            biotype = info.split('gene_biotype "')[1].split('"')[0]
            g_dt.loc[len(g_dt)] = [chrom,int(start), int(stop), strand, gene_name, biotype]
    
    all_dt = pd.merge(g_dt, o_dt, how = 'left', left_on = 'GeneName', right_on = 'mzebra')
    summary_data = pd.DataFrame(index = ['immunoglobin_Tcell','ncRNA','protein','pseudogene','Error','NoHomolog','Ortholog','Paralog','Error','OneToMany','OneToOne','OneToTwo'])
    data = gi_dt.groupby('PrimaryBiotype').count()['GeneID'].values.tolist() + gi_dt.groupby('HumanRelationship').count()['GeneID'].values.tolist() + gi_dt.groupby('HumanOrthologRelationship').count()['GeneID'].values.tolist()
    summary_data['All'] = data
    with pd.ExcelWriter(combined_gene_file) as writer:  
        gi_dt.to_excel(writer, sheet_name='AllGenes')
        for lg, location in inversions.items():
            contig = location.split(':')[0]
            start,stop = [int(x) for x in location.split(':')[1].split('-')]
            lg_dt = gi_dt[(gi_dt.Scaffold == contig) & (gi_dt.Start > start) & (gi_dt.Stop < stop)]
            data = lg_dt.groupby('PrimaryBiotype').count()['GeneID'].values.tolist() + lg_dt.groupby('HumanRelationship').count()['GeneID'].values.tolist() + lg_dt.groupby('HumanOrthologRelationship').count()['GeneID'].values.tolist()
            try:
                summary_data[lg] = data
            except:
                summary_data[lg] = data + [0]*(12-len(data))

            lg_dt.to_excel(writer, sheet_name = lg)
    pdb.set_trace()
    #all_dt.to_csv(combined_gene_file)
    fm_obj.uploadData(combined_gene_file)

# Create fm_obj and grab sample file
fm_obj = FM(genome_version = 'Mzebra_UMD2a')
sample_file = fm_obj.localReadsDir + 'SampleDatabase_v2.xlsx'

fm_obj.downloadData(sample_file)
dt = pd.read_excel(sample_file, sheet_name = 'SampleLevel')

#make_inversion_genefile(fm_obj)
#pdb.set_trace()
# Create vcf file for lab sample data and for lake malawi phylogeny data (defined in Sample Database)
lab_samples = dt[dt.LabStrain == 'Yes'].SampleID.to_list()
malawi_samples = dt[dt.CorePCA == 'Yes'].SampleID.to_list()
dt = dt[dt.CorePCA == 'Yes']

mbuna_samples = dt[dt.Ecogroup_PTM == 'Mbuna'].SampleID.to_list()
benthic_samples = dt[(dt.Ecogroup_PTM == 'Deep_Benthic') | (dt.Ecogroup_PTM == 'Shallow_Benthic') | (dt.Ecogroup_PTM == 'Utaka')].SampleID.to_list()
rhamp_samples = dt[dt.Ecogroup_PTM == 'Rhamphochromis'].SampleID.to_list()
diplo_samples = dt[dt.Ecogroup_PTM == 'Diplotaxodon'].SampleID.to_list()
ac_samples = dt[dt.Ecogroup_PTM == 'AC'].SampleID.to_list()
castle_samples = dt[dt.LG11 == 'Inverted'].SampleID.to_list()
pit_samples = dt[((dt.Ecogroup_PTM == 'Deep_Benthic') | (dt.Ecogroup_PTM == 'Shallow_Benthic') | (dt.Ecogroup_PTM == 'Utaka')) & (dt.LG11 == 'Normal_G1')].SampleID.to_list()
deep_samples = dt[dt.LG2 == 'Inverted'].SampleID.to_list()
inversion_samples = {'Mbuna':mbuna_samples,'Benthic':benthic_samples,'Rhamp':rhamp_samples,'Diplo':diplo_samples,'Castle':castle_samples,'Pit':pit_samples,'Deep':deep_samples}


inv2_samples = dt[dt.LG2 == 'Inverted'].SampleID.to_list()
inv9_samples = dt[(dt.LG9 == 'Inverted_Group1') | (dt.LG9 == 'Inverted_Group2')].SampleID.to_list()
inv10_samples = dt[dt.LG10 == 'Inverted'].SampleID.to_list()
inv11_samples = dt[dt.LG11 == 'Inverted'].SampleID.to_list()
inv11_non_samples = dt[(dt.LG11 == 'Normal_G1') & (dt.Ecogroup_PTM != 'Mbuna')].SampleID.to_list()
inv13_samples = dt[dt.LG13 == 'Inverted'].SampleID.to_list()
inv20_samples = dt[dt.LG20 == 'Inverted'].SampleID.to_list()
#inversion_samples = {'Inv2':inv2_samples,'Inv9':inv9_samples,'Inv10':inv10_samples,'Inv11':inv11_samples,'Inv13':inv13_samples,'Inv20':inv20_samples}
inversion_dts = []
# This is assumed to be downloaded in the right spot
input_vcf = fm_obj.localMasterDir + 'pass_variants_master_file.vcf.gz'
hf_vcf = fm_obj.localMasterDir + 'Malawi_AC_20.vcf.gz'

inv2_vcf = fm_obj.localMasterDir + 'inv2_samples.vcf.gz'
# Name of output files
lab_samples_vcf = fm_obj.localMasterDir + 'LabIndividuals.vcf.gz'
malawi_samples_vcf = fm_obj.localMasterDir + 'MalawiIndividuals.vcf.gz'

# This filters by sample and removes genotypes that aren't present in the samples and also indels and other SNP types
#parallel_filter_population_af(input_vcf, hf_vcf, malawi_samples)
#subprocess.run(['plink2','--vcf',hf_vcf,'--make-pgen','--double-id','--allow-extra-chr','--out',hf_vcf.replace('.vcf.gz',''), '--vcf-min-gq', '10'])
subprocess.run(['plink2','-pfile',hf_vcf.replace('.vcf.gz',''),'--freq','cols=+pos','--allow-extra-chr','-out',fm_obj.localMasterDir + 'AllSamples_af'])
inversion_dts.append(pd.read_csv(fm_obj.localMasterDir + 'AllSamples_af.afreq', sep = '\t'))
for inv,samples in inversion_samples.items():
    with open(fm_obj.localMasterDir + inv + '_samples.txt','w') as f:
        for sample in samples:
            print(sample + '\t' + sample, file = f)
    subprocess.run(['plink2','-pfile',hf_vcf.replace('.vcf.gz',''),'--keep',fm_obj.localMasterDir + inv + '_samples.txt', '--freq','cols=+pos','--allow-extra-chr','-out',fm_obj.localMasterDir + inv + '_af'])
    inversion_dts.append(pd.read_csv(fm_obj.localMasterDir + inv +'_af.afreq', sep = '\t', usecols = ['#CHROM','POS','ALT_FREQS','OBS_CT']))
    inversion_dts[-1] = inversion_dts[-1].rename(columns = {'ALT_FREQS':inv + '_af','OBS_CT':inv + '_ct'})
from functools import reduce
df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['#CHROM','POS']), inversion_dts)
ann_vcf = fm_obj.localMasterDir + 'eff.vcf.gz'
vcf_obj = vcf.VCFReader(filename = ann_vcf)

chr11 = df_merged[((df_merged.Castle_af > 0.95) & (df_merged.Pit_af < 0.2)) | ((df_merged.Pit_af > 0.8) & (df_merged.Castle_af < 0.05))]
chr11[['#CHROM','POS']].to_csv(fm_obj.localMasterDir + 'chr11_pos.tsv', sep = '\t', index = False)
subprocess.call(['bcftools','view','-R', fm_obj.localMasterDir + 'chr11_pos.tsv', '-o', fm_obj.localMasterDir + 'chr11_eff.vcf', ann_vcf])
pdb.set_trace()
i = 0
positions = chr11['POS'].to_list()
for rec in vcf_obj:
    if i % 10000 == 0:
        print(i)
    i+=1
    if rec.CHROM == 'NC_036790.1' and rec.POS in positions:
        pdb.set_trace()

pdb.set_trace()

parallel_filter_af(input_vcf, inv2_vcf, inv2_samples)
parallel_filter(input_vcf, lab_samples_vcf, lab_samples)
parallel_filter(input_vcf, malawi_samples_vcf, malawi_samples)

#Create lg subsets for each inversion from the malawi samples vcf file
for lg,region in inversions.items():
    lg_vcf = fm_obj.localMasterDir + 'MalawiIndividuals_' + lg + '.vcf.gz'
    subprocess.run(['bcftools','view','-r', region, '-o', lg_vcf, '-O', 'z', malawi_samples_vcf])

#Create relaxed phylip file from from vcf. Requires vcf2phylip to be installed in your PYTHONPATH
subprocess.run(['python3', '-m', 'vcf2phylip','-i',lab_samples_vcf, '-m', '40', '-o',fm_obj.localMasterDir])
subprocess.run(['python3', '-m', 'vcf2phylip','-i',malawi_samples_vcf, '-m', '150','-o',fm_obj.localMasterDir])

for lg,region in inversions.items():
    if lg not in ['LG13','LG20']:
        continue
    lg_vcf = fm_obj.localMasterDir + 'MalawiIndividuals_' + lg + '.vcf.gz'
    subprocess.run(['python3', '-m', 'vcf2phylip','-i',lg_vcf, '-m', '150','-o', fm_obj.localMasterDir])

# Convert relaxed phylip file to stockholm using BioPYTHON
align = AlignIO.read(lab_samples_vcf.replace('.vcf.gz','.min40.phy'), 'phylip-relaxed')
AlignIO.write(align,lab_samples_vcf.replace('.vcf.gz','.min40.stk'),'stockholm')

align = AlignIO.read(malawi_samples_vcf.replace('.vcf.gz','.min150.phy'), 'phylip-relaxed')
AlignIO.write(align,malawi_samples_vcf.replace('.vcf.gz','.min150.stk'),'stockholm')

for lg,region in inversions.items():
    lg_vcf = fm_obj.localMasterDir + 'MalawiIndividuals_' + lg + '.vcf.gz'
    align = AlignIO.read(lg_vcf.replace('.vcf.gz','.min150.phy'), 'phylip-relaxed')
    AlignIO.write(align,lg_vcf.replace('.vcf.gz','.min150.stk'),'stockholm')

# Run quicktree to create phylogeny file
subprocess.run(['quicktree', lab_samples_vcf.replace('.vcf.gz','.min40.stk')], stdout = open(lab_samples_vcf.replace('.vcf.gz','.min40.tree'),'w'))
subprocess.run(['quicktree', malawi_samples_vcf.replace('.vcf.gz','.min150.stk')], stdout = open(malawi_samples_vcf.replace('.vcf.gz','.min150.tree'),'w'))
for lg,region in inversions.items():
    lg_vcf = fm_obj.localMasterDir + 'MalawiIndividuals_' + lg + '.vcf.gz'
    subprocess.run(['quicktree', lg_vcf.replace('.vcf.gz','.min150.stk')], stdout = open(lg_vcf.replace('.vcf.gz','.min150.tree'),'w'))


# Replace the sample ID with names that are more interesting to display on the tree
treefiles = [lab_samples_vcf.replace('.vcf.gz','.min40.tree'),malawi_samples_vcf.replace('.vcf.gz','.min150.tree')] + [fm_obj.localMasterDir + 'MalawiIndividuals_' + x + '.min150.tree' for x in inversions.keys()]

for tf in treefiles:
    f = open(tf)
    data = f.read()
    for i,row in dt.iterrows():
        sample = row.SampleID
        species = row.Organism
        ecogroup = row.Ecogroup_PTM
        small_species = species[0] + '_' + species.split()[1]
        new_text = ecogroup + '_' + small_species + '_' + sample
        data = data.replace(sample,new_text)
    out = open(tf.replace('.tree','.modified.tree'), 'w')
    print(data, file = out, end = '')

# This is work to try to identify high heterozygosity regions but it doesn't really identify anything
output = subprocess.run(['vcftools','--gzvcf','LabIndividuals.vcf.gz','--het','-c'], capture_output = True)
pd.read_csv(StringIO(output.stdout.decode('utf-8')), sep = '\t')

lab_samples_haf_vcf = fm_obj.localMasterDir + 'LabIndividuals_AF.1.recode.vcf.gz'
for contig,length in lengths.items():
    start = 0
    stop = 200000
    while stop < length:
        position = contig + ':' + str(start) + '-' + str(stop)
        p1 = subprocess.Popen(['bcftools','view',lab_samples_haf_vcf,'-r',position], stdout = subprocess.PIPE)
        p2 = subprocess.Popen(['vcftools','--vcf', '-', '--het','-c'], stdin = p1.stdout, stdout = subprocess.PIPE)
        data = StringIO(p2.communicate()[0].decode('utf-8'))
        dt = pd.read_csv(data, sep='\t')
        dt[position] = (dt['N_SITES'] - dt['O(HOM)'])/dt['N_SITES']

        p_dt = pd.pivot_table(data = dt, columns = 'INDV', values = position)
        p_dt['N_Sites'] = dt.iloc[0]['N_SITES']
        try:
            out_data = pd.concat([out_data,p_dt])
        except NameError:
            out_data = p_dt

        start += 20000
        stop += 20000
