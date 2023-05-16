import argparse, pdb, os, subprocess, pathlib
import pandas as pd
from cyvcf2 import VCF
import plotly.express as px
import plotly.graph_objs as go

class PCA_Maker:
    def __init__(self, input_vcffile, output_directory, sample_database, ecogroups, linkage_groups): # The PCA_Maker class will create an object (self) which will then take in (currently) 4 pieces of information (input file out dir, sample_database excel sheet, ecogroup names)
        # self.attr indicates that the object made with PCA_Maker will have the attribute named "attr"
        # The linkage_group_map attribute is a hard coded list of LG names
        self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1', 
                             'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1', 
                             'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1', 
                             'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1'}
        self.in_vcf = input_vcffile # The in_vcf attriubute is sequal to the input file name for the Object. 
        self.sample_database = sample_database # The sample_database attribute equals the sample_database name for the object
        self.ecogroups = ecogroups # This attribute is the list of ecgogroups used for filtering samples
        self.r_script = os.getcwd() + '/modules/pca.R' # use the dir from which the pipeline is called to get the file path to the R script. 

        self.vcf_obj = VCF(self.in_vcf) # The VCF object is made using the CVF class from cyvcf2. The VCF class takes an input vcf. For the object, this input vcf file is the "input_vcfcfile" which is defined under self.in_vcf
        self.linkage_groups = linkage_groups # the object will now take in linkage groups using args.regions. If default, it will default to the first 22 lgs
        # Create a filepath that includes name of input file and ecogroups on whcih the analysis is performed
        file_version = self.in_vcf.split('/')[-1].split('.')[0]
        analysis_ecogroups = '_'.join(str(eg).replace('_', '') for eg in self.ecogroups)
        pathlib.Path(output_directory + '/' + file_version + '/' + analysis_ecogroups).mkdir(parents=True, exist_ok=True) # This generates the outdir if it doesn't exist so later tools don't run into errors making files.
        self.out_dir = output_directory + '/' + file_version + '/' + analysis_ecogroups

        regions_list = []
        for region in self.linkage_groups:
            if region in self.linkage_group_map.keys():
                regions_list.append(self.linkage_group_map[region])
            elif region == "All":
                regions_list.extend(self.vcf_obj.seqnames[0:22])
            elif region == "Inversion":
                regions_list.extend(['pre_inversion', 'inversion1', 'inversion2', 'post_inversion'])
            elif region == "Whole":
                regions_list.append('Whole')
            else:
                raise Exception(region + ' is not a valid option')
        self.linkage_groups = regions_list
        duplicate_set_test = set(self.linkage_groups)
        if len(duplicate_set_test) != len(self.linkage_groups):
            raise Exception('A repeat region has been provided')
        # Ensure index file exists
        assert os.path.exists(self.in_vcf + '.tbi') # uses os.path.exists to see if the input file + 'tbi' extension exists. The object will be made using args.input_vcffile and args.input_vcffile will be passed to the script as an absolute file path so the path to the dir is taken care of

# Step 1: Identify unique species
unique_species = [...]  # List of unique species in your dataset

# Step 2: Randomly select up to three samples per species
selected_samples = []
for species in unique_species:
    species_samples = [...]  # List of samples belonging to the current species
    num_samples = min(3, len(species_samples))  # Select up to three samples
    random_samples = random.sample(species_samples, num_samples)
    selected_samples.extend(random_samples)

# Step 3: Create a new VCF file with selected samples
subprocess.run(['plink', '--vcf', self.samples_filtered_master_vcf, '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--keep', 'sample_list.txt', '--make-bed', '--out', self.out_dir + '/PCA/' + 'Whole/' + 'selected_samples'])

# Step 4: Perform PCA analysis on the new VCF file
subprocess.run(['plink', '--bfile', self.out_dir + '/PCA/' + 'Whole/' + 'selected_samples', '--pca', '--out', self.out_dir + '/PCA/' + 'Whole/' + 'selected_samples'])

# Step 5: Project eigenvalues of remaining samples onto the selected principal components
subprocess.run(['plink', '--bfile', self.out_dir + '/PCA/' + 'Whole/' + 'selected_samples', '--score', 'remaining_samples_list.txt', 'eigenvec', '--out', self.out_dir + '/PCA/' + 'Whole/' + 'projected_scores'])
