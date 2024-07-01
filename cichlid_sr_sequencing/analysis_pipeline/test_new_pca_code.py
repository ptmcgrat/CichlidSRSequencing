import subprocess, pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from helper_modules.nikesh_file_manager import FileManager as FM
import plotly.express as px
import plotly.graph_objs as go
import pathlib

"""
TODO:
issue where repeat sampleIDs end up in the filtered sample df
A temp fix has been implemented but a more thorough fix that will sample until 3 unique samples are chosen needs to be implemented in the first if block that checks if over 3 samples exist for a given organism


Dec. 12, 2023
The code is ready to be built into pca_maker.py
Here's what each part of pca_maker does and what pieces below will go where:

- __init__ initializes many variables:
    - The Regions list aka the linkage groups or the special exploratory regions for which a PCA is being generated 
    - the ecogroups used in the analysis which informs which samples will be used for generating the subset vcf file
        - NOTE: The subset vcf file when using the --sample_subset flag is incorrect. The vcf file generated only contains at max 3 samples per species and only maps those samples.
        - This must change so that teh vcf file is generated using the whole cohort as informed by the ecogroups and a subset file is created with the samples that will be used for the subsequent analysis.
        - May be a good idea to include the ecogroup, and species name to see how many per species are actually being used. The first column (SampleIDs) can be used as input into the subsequent code

- _create_sample_filter_file is used to filter samples based on the ecogroup. Lots of updates to do here:
    - This function will generate the filter file as informed by the eco groups. Note that the subset sample filename should be perhaps changed from "samples_to_keep" to "subset_samples"
    - After the ecogroups are defined, the SampleDatabase.csv file is downloaded and passed to pandas
    - The new filtering code will get incorporated into this section 
    - NOTE: The new code should be incorporated into the if args.sample_subset section (I already have a Pacbio fix in that code that I may want to use)
    - NOTE: The bcftools command that generates the subset vcf file will need to operate on a sample list informed by just filtering based on the ecogroup and not on the sample subset 
        - The subset sample file is used exclusively to inform the subset PCA on which samples will be used during the initial PCA tpo create the .acounts and .eigenvec.alles files
    - This command also takes care of generating the subset vcf file... Strongly consider making this its own function for ease of editing in the future

-_split_VCF_to_LG split large vcf to subregion vcfs:
    - This function takes in the regions list (linkage_group_list) and starts to split up the vcf file's variants by those regions
    - I think most of this could be left alone?? 

- _create_exploratory_region_eigen_files
    - Ignore for now

- _create_eigenfiles_per_LG is where plink2 is called:
    - For each regions VCF file, this runs the plink code to generate the eigenvector files that will be used to plot samples in the next function.
    - Should be straight forward command replacement and file path remapping to do here

- _create_interactive_pca plotly code for plotting the eigenvectors from the first 2 PCs per sample
    - should be some minor edits to do following this script's code
"""

# Load PCA results from the subset
fm = FM('Mzebra_UMD2a')
out_dir = '/Users/kmnike/Data/CichlidSequencingData/Outputs/pca_testing/'

"""
# old samples filtering pandas code
s_df = pd.read_csv(fm.localSampleFile)
no_pacbio_df = s_df[s_df['Platform'].isin(['ILLUMINA'])]
df_filtered = pd.DataFrame()

# generate a random sample of 3 individuals per species at max for vcf subsampling. Ensure no sampleID gets repeated as its being added:
for organism in no_pacbio_df['Organism'].unique():
    organism_rows = no_pacbio_df[no_pacbio_df['Organism'] == organism]
    if organism_rows.shape[0] > 3: # organism_rows.shape[0] gives the number of row. If greater than 3 do the following. If exactly 3 or less, execute the else code
        sampled_rows = organism_rows.sample(n=3, random_state=42) # randomly sample 3 samples from the species.
        df_filtered = pd.concat([df_filtered, sampled_rows], ignore_index=True).drop_duplicates('SampleID') # add them to df_filtered, and get rid of any duplicates 
    else:
        sampled_rows = organism_rows # if only 3 or less samples exist, the sampled_rows are the organism_rows
        df_filtered = pd.concat([df_filtered, sampled_rows], ignore_index=True).drop_duplicates('SampleID') # they are added to df_filtered and repeatr samples are dropped as above.

df_filtered.to_csv('subset_samples.csv', columns = ['SampleID'], index=False, header=False)
df_filtered.to_csv('subset_samples_with_metadata.csv', columns = ['SampleID', 'Ecogroup', 'Organism'], sep = '\t', index=False)
pdb.set_trace()
egs = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'Riverine', 'AC']
ecogroup_df = s_df[s_df.Ecogroup.isin(egs)]
ecogroup_df.drop_duplicates(subset='SampleID').to_csv('all_ecogroup_samples_csv', columns = ['SampleID'], header = False, index=False)
"""

# generate a subset vcf file from the ecogroup-specific vcf file. For testing, I am working with "All" ecogroups
# NOTE: there is a difference in the s_df that the test script and the pipeline is reading in casuign the pca_maker to get 609 samples vs 615 in this script.... 
# NOTE: the goal upon rwestarting this tomrorow is to replciate the pca analysis that this script can do, but using the pca_maker pipeline. This may help tease out things going wrong earlier in the script before we introduce LG specific complication


# Trying a new filtering of the SampleDatabase:
df_filtered = pd.DataFrame()
s_df = pd.read_csv(fm.localSampleFile)
s_df = s_df[s_df['Platform'].isin(['ILLUMINA'])].drop_duplicates(subset='SampleID') # get rid of PacBio Samples and drops duplicates in the SampleID column leaving 612 samples that we can now filter
# confirmed that samples in the above s_df exactly match the sampels returned from bcftools query -l 612_sample.cvf using s_df['SampleID'].isin(df_612_samples_list['SampleID']]).all() (returns True)
# Now we can filter for ecogroup samples, like excluding the 'Outlier' and 'Deep/Shallow Benthic' Samples:
egs = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'Riverine', 'AC']
e_df = s_df[s_df.Ecogroup.isin(egs)]
e_df.to_csv('ecogroup_samples.csv', columns = ['SampleID'], header=False, index=False)
for organism in e_df.Organism.unique():
    organism_rows = e_df[e_df['Organism'] == organism]
    if organism_rows.shape[0] > 3: # organism_rows.shape[0] gives the number of row. If greater than 3 do the following. If exactly 3 or less, execute the else code
        sampled_rows = organism_rows.sample(n=3, random_state=42) # randomly sample 3 samples from the species.
        df_filtered = pd.concat([df_filtered, sampled_rows], ignore_index=True)
    else:
        sampled_rows = organism_rows # if only 3 or less samples exist, the sampled_rows are the organism_rows
        df_filtered = pd.concat([df_filtered, sampled_rows], ignore_index=True)
df_filtered.to_csv('subset_samples.csv', columns = ['SampleID'], index=False, header=False)
df_filtered.to_csv('subset_samples_with_metadata.csv', columns = ['SampleID', 'Ecogroup', 'Organism'], sep = '\t', index=False)




"""
The subset VCF file is generated per every linkage group
This will allow a subset PCA to be run for each region (lg) 
The ecogroup specific VCF file is used as the input to generate the subset VCF file per LG

A second subset VCF file for ALL samples but region specific will need to be written into each LGs directory
Pfiles need to be written for both the subset vcf file for the region, as well as the all_ecogroup speciifc samples for that region
I was thinking that one set of files needs to be written for the ecogroup vcf file, but that's nt true since those pfiles will summarize variants across ALL regions which is incorerct  


Steps to take: 
write a subset vcf file in each region's dir where samples are pulled from the subset_samples.csv file
write a subset_ecogroup vcf file in each regions dir where ALL of the samples for that analysis are pulled from the all_ecogroyp_samples.csv file. 


"""

print('FILTERING THE WHOLE VCF FILE TO INCLUDE SAMPLES IN ECOGROUPS')
subprocess.run(['bcftools', 'view', fm.localOutputDir + 'vcf_concat_output/10k_variants_612_cohort_3_lg_subset.vcf.gz', '--samples-file', 'ecogroup_samples.csv', '--min-ac=1', '--no-update', '-o', fm.localOutputDir + 'pca_testing/ecogroup_samples.vcf.gz'])

print('GENERATING SUBSET WITH BCFTOOLS...')
subprocess.run(['bcftools', 'view', fm.localOutputDir + 'vcf_concat_output/10k_variants_612_cohort_3_lg_subset.vcf.gz', '--samples-file', 'subset_samples.csv', '-o', fm.localOutputDir + 'pca_testing/subset_samples.vcf.gz', '-O', 'z'])

# Once the subset and ecogroup files are generated, I need to recalculate the AF field in the subset VCF file:
# seems like I can't overwrite the old subset_samples.vcf.gz file
print('RECALCULATING AF IN THE SUBSET FILE')
subprocess.run(['bcftools', '+fill-tags', fm.localOutputDir + 'pca_testing/subset_samples.vcf.gz', '-Oz', '-o', fm.localOutputDir + 'pca_testing/AF_recalculated_subset_samples.vcf.gz', '--', '-t', 'AF'])

print('FILTERING OUT VARIANTS WITH AF < 0.01')
subprocess.run(['bcftools', 'view', '--exclude', 'AF<0.01', fm.localOutputDir + 'pca_testing/AF_recalculated_subset_samples.vcf.gz', '-o', fm.localOutputDir + 'pca_testing/af_filtered_subset_samples.vcf.gz', '-Oz'])

print('WRITING CHRMOSOME-POSITION FILE')
subprocess.run(['bcftools', 'query', '-f', '%CHROM\t%POS', fm.localOutputDir + 'pca_testing/af_filtered_subset_samples.vcf.gz', '-o', fm.localOutputDir + 'pca_testing/variants_to_keep.txt'])

print('PRUNING VARIANTS FROM ECGOGROUP VCF FILE USING THOSE PRESENT IN THE SUBSET')
subprocess.run(['vcftools', '--positions', fm.localOutputDir + 'pca_testing/variants_to_keep.txt', '--gzvcf', fm.localOutputDir + 'pca_testing/ecogroup_samples.vcf.gz', '--recode', '--out', fm.localOutputDir + 'pca_testing/variants_filtered_ecogroup_samples'])
# vcftools --positions test.txt --gzvcf ecogroup_samples.vcf.gz --recode --out variants_filtered_ecogroup_samples

# before the below code can work, I will now implement the allele frequency filtering above
# generate the pfiles from the subset vcf and whole vcf file 
# code will work again if the pfiles generated are for the original vcf file instead of the ecogroup vcf file. Uncomment the line with the path to the original vcf file and it will work 
print('GENERATING PFILES FOR WHOLE VCF FILE')
# subprocess.run(['plink2', '--vcf', fm.localOutputDir + 'vcf_concat_output/10k_variants_612_cohort_3_lg_subset.vcf.gz', '--out', fm.localOutputDir + 'pca_testing/whole', '--allow-extra-chr'])
subprocess.run(['plink2', '--vcf', fm.localOutputDir + 'pca_testing/variants_filtered_ecogroup_samples.recode.vcf', '--out', fm.localOutputDir + 'pca_testing/whole', '--allow-extra-chr'])

print('GENERATING PFILES FOR SUBSET VCF FILE')
subprocess.run(['plink2', '--vcf', fm.localOutputDir + 'pca_testing/af_filtered_subset_samples.vcf.gz', '--out', fm.localOutputDir + 'pca_testing/subset', '--allow-extra-chr'])

# add variant IDs to the whole dataset files to avoid crashing at pca steps
print('SETTING MISSING VARIANT IDS IN WHOLE PGEN FILES...')
subprocess.run(['plink2', '--pfile', fm.localOutputDir + 'pca_testing/whole', '--set-missing-var-ids', '@:#', '--make-pgen', '--out', fm.localOutputDir + 'pca_testing/whole_corrected', '--allow-extra-chr'])

print('PCA STEP 1 STARTING...')
subprocess.run(['plink2', '--pfile', fm.localOutputDir + 'pca_testing/subset', '--freq', 'counts', '--pca', 'allele-wts', '--out', fm.localOutputDir + 'pca_testing/subset_pca_test', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--max-alleles', '2'])

# cleaning .acounts file: to eliminate variants where ALT counts are zero likely due to no variants at that position from the random sampling 
# pdb.set_trace()
acounts_df = pd.read_csv(fm.localOutputDir + 'pca_testing/subset_pca_test.acount', sep='\t')
nonzero_acounts_df = acounts_df[acounts_df['ALT_CTS'] != 0]
nonzero_acounts_df.to_csv(fm.localOutputDir + 'pca_testing/subset_pca_test.acount', sep='\t', index=False)

print('PCA STEP 2 STARTING...')
subprocess.run(['plink2', '--pfile', fm.localOutputDir + 'pca_testing/whole_corrected', '--read-freq', fm.localOutputDir + 'pca_testing/subset_pca_test.acount', '--score', fm.localOutputDir + 'pca_testing/subset_pca_test.eigenvec.allele', '2', '5', 'header-read',  'no-mean-imputation', 'variance-standardize', '--score-col-nums', '6-15', '--out', fm.localOutputDir + 'pca_testing/new_projection', '--allow-extra-chr'])

# Plotting code which will need updating with new eigenvector file names 
pca_files = ('subset_pca_test.eigenvec', 'new_projection.sscore')
for file in pca_files:
    plotly_out = out_dir + 'interactive_PCA_outputs/'
    pathlib.Path(plotly_out).mkdir(parents=True, exist_ok=True)
    color_map = {'Mbuna': 'purple', 'AC': 'limegreen', 'Shallow_Benthic': 'red', 'Deep_Benthic': 'blue', 'Rhamphochromis': 'brown', 'Diplotaxodon': 'orange', 'Utaka': 'darkgreen', 'Riverine': 'pink'}
    shape_map = {'PRJEB15289': 'square', 'PRJEB1254': 'circle', 'RockSand_v1': 'diamond', 'ReferenceImprovement': 'x', 'BrainDiversity_s1': 'star', 'BigBrain': 'triangle-up', 'Multiome': 'cross'}
    eigen_df = pd.read_csv(out_dir + file, sep='\t')
    eigen_df = eigen_df.rename(columns = {'#IID':'SampleID'})
    metadata_df = pd.read_csv(fm.localSampleFile)
    metadata_df = metadata_df.drop_duplicates(subset='SampleID', keep='first')
    df_merged = pd.merge(eigen_df, metadata_df, on=['SampleID'])
    if file.endswith('.sscore'):
        fig = px.scatter(df_merged, x='PC1_AVG', y='PC2_AVG', color='Ecogroup', symbol='ProjectID', color_discrete_map=color_map, symbol_map=shape_map, title=file, hover_data=['SampleID', 'Ecogroup', 'Organism', 'ProjectID'])
    else:
        fig = px.scatter(df_merged, x='PC1', y='PC2', color='Ecogroup', symbol='ProjectID', color_discrete_map=color_map, symbol_map=shape_map, title=file, hover_data=['SampleID', 'Ecogroup', 'Organism', 'ProjectID'])
    larger_size = 10  # You can adjust this value
    fig.update_traces(marker=dict(size=larger_size), selector=dict(marker_symbol='cross'))
    fig.write_html(plotly_out + file + '_plotlyPCA.html')

print('PIPELINE RUN SUCCESS')



"""
NOTE:
The PCs in the resultant sscore file will be scaled a bit differently from subset_data.eigenvec; you need to multiply or divide the PCs by a multiple of sqrt(eigenvalue) to put them on the same scale.
Thus, I believe this means each value in PC1 column of the sscore file needs to get multiplied by the sqrt(first_eigenvalue) in the subset_data eigenval file

"""