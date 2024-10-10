import argparse, pdb, os, subprocess, pathlib
import pandas as pd # sometimes pandas will need to be conda remove'd and pip uninstall'd to install python >=3.7 to get this script to work. pip install pandas again afterwards.
from helper_modules.nikesh_file_manager import FileManager as FM
from cyvcf2 import VCF
import plotly.express as px # Do not conda install plotly . Use this: pip install plotly==5.11.0
import plotly.graph_objs as go # Do not conda install plotly . Use this: pip install plotly==5.11.0
# import umap

parser = argparse.ArgumentParser(usage = "This pipeline is for running pca analysis on a filtered vcf file. Note that the script will assume tha name of the vcf file to analyze is pass_variants_master_file.vcf.gz")
parser.add_argument('genome', help = 'name of reference genome used in the creation of the VCF files')
parser.add_argument('output_dir', help = 'absolute filepath to an output directory')
parser.add_argument('--sample_subset', help = 'This flag will restruct the samples used for generating the PCA plot to a max of three three samples for any given organism.', action= 'store_true')
parser.add_argument('-e', '--ecogroups', help = 'one or multiple eco group names for filtering the data', choices = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'Riverine', 'AC', 'Non_Riverine', 'All', 'Lake_Malawi', 'Rock_Sand', 'Sand', 'Custom'], nargs = '*', default = ['All'])
parser.add_argument('-r', '--regions', help = 'list of linkage groups for which analyses will run', nargs = '*', default = ['All'])
parser.add_argument('-l', '--local_test', help = 'call this flag to predefine variables for testing on local machine', action='store_true')
parser.add_argument('-p', '--plink', help = 'use this flag to generate new eigenvalue/vector files using plink', action='store_true')
args = parser.parse_args()

"""
TODO:
2024.09.16
Patrick wants some specific things done with pca_maker for the bionanop_paper analyses
1. Run the 498 sample cohort normally
2. include an LG7 inversion region based on the coordinated I found in the Bionano data. I can just go ahead and run this region for all sub-ecogroups
3. Run the pipeline in a way that uses all "core_PCA" samples to calculate the PCs, then project only the MCxYH F1s on the axes to see if anything separates out in PC space for these small # of samples. 

time python pca_maker.py Mzebra_GT3 /Users/kmnike/Data/pca_testing -e Custom --local_test --sample_subset
time python pca_maker.py Mzebra_GT3 /Users/kmnike/Data/pca_testing -e Rock_Sand --local_test --sample_subset --plink -r Exploratory
time python pca_maker.py Mzebra_GT3 /Users/kmnike/Data/pca_testing -e Custom --local_test --sample_subset --plink
CVAnalysis
"""
# The class PCA_Maker will create objects that will take in a variety of inputs (generally determined by what input parametrs are being passed into the script).
# These objects will have many attributes which will serve to help build directory structure, define valid inputs, etc.
# There will be many fucntions defined within the class besides the __init__ function which give objects of the class their attributes.
class PCA_Maker:
    def __init__(self, genome, ecogroups, linkage_groups, output_directory): # The PCA_Maker class will create an object (self) which will then take in (currently) 4 pieces of information (input file out dir, sample_database excel sheet, ecogroup names)
        # self.attr indicates that the object made with PCA_Maker will have the attribute named "attr"
        # The linkage_group_map attribute is a hard coded list of LG names
        self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1', 
                             'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1', 
                             'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1', 
                             'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1', 'mito': 'NC_027944.1'}
        self.genome = genome
        self.fm_obj = FM(self.genome)
        self.in_vcf = self.fm_obj.localOutputDir + 'vcf_concat_output/498_cohort_pass_variants_master_file.vcf.gz' # The in_vcf attriubute is equal to the input file name for the Object.
        if args.local_test:
            self.in_vcf = self.fm_obj.localOutputDir + 'vcf_concat_output/local_testing_498_master_file.vcf.gz' # This file is a subset of the 612 cohort pass_variants_master_file. By default, the script will use this whole file as input for local testing 
        self.ecogroups = ecogroups # This attribute is the list of ecogroups used for filtering samples

        self.vcf_obj = VCF(self.in_vcf) # The VCF object is made using the CVF class from cyvcf2. The VCF class takes an input vcf. For the object, this input vcf file is the "input_vcfcfile" which is defined under self.in_vcf
        self.linkage_groups = linkage_groups # the object will now take in linkage groups using args.regions. If default, it will default to the first 22 lgs
        # Create a filepath that includes name of input file and ecogroups on which the analysis is performed. The below code is necessary to make sure the storage of PCAs are within a filepath unique to each analysis and unique for any VCF file passed in
        self.out_dir = output_directory
        if self.out_dir.endswith('/'): # ensure that the file path does not end with a forwardslash(/)
            self.out_dir = self.out_dir[:-1]

        regions_list = [] # regions list will add various linkage group names, or add the names of regions of interest like "Whole" or the various special regions of interest if passed "Exploratory"
        # predefine one exploratory_regions_list so the variable can be cointinually reused.
        self.exploratory_regions_list = ['lg2_Mbuna_inversion','lg2_Deep_Benthic_Inversion', 'lg2_non_inverted_region', 'lg7_Inversion', 'lg8_benthic_inversion', 'lg9_Whole_RockSand_Inversion', 'lg9a_RockSand_Inversion', 'lg9b_RockSand_Inversion', 'lg9c_RockSand_Inversion', 'lg9d_RockSand_Inversion', 'lg9e_RockSand_Inversion', 'lg9_non_inverted_region', 'lg10_Deep_Benthic_Inversion', 'lg10_non_inverted_region', 'lg11a_Inversion', 'lg11b_Inversion', 'lg11_non_inverted_region', 'lg12_Pelagic_inversion', 'lg13_Deep_Benthic_Inversion', 'lg13_non_inverted_region', 'lg20a_RockSand_Inversion', 'lg20b_RockSand_Inversion', 'lg20_non_inverted_region']
        for region in self.linkage_groups:
            if region in self.linkage_group_map.keys():
                regions_list.append(self.linkage_group_map[region])
            elif region == "All":
                regions_list.extend(self.vcf_obj.seqnames[0:22])
            elif region == "Exploratory":
                regions_list.extend(self.exploratory_regions_list)
            elif region == "Whole":
                regions_list.append('Whole')
            elif region in self.exploratory_regions_list:
                regions_list.append(region)
            else:
                raise Exception(region + ' is not a valid option')

        self.linkage_groups = regions_list # regions list is assigned to self.linkage_groups which will be used throughout the pipeline
        duplicate_set_test = set(self.linkage_groups) # test if any linkage groups are passed twice and throws error if True:
        if len(duplicate_set_test) != len(self.linkage_groups):
            raise Exception('A repeat region has been provided')
        # reset self.linkage_groups to include only the 3 lgs tested locally  if the --local_test flag is called
        if args.local_test:
            self.linkage_groups = ['NC_036788.1', 'NC_036781.1', 'NC_036782.1']
            if 'Whole' in args.regions:
                self.linkage_groups.extend(['Whole'])
            if 'Exploratory' in args.regions:
                self.linkage_groups.extend(self.exploratory_regions_list)

        # Ensure index file exists
        assert os.path.exists(self.in_vcf + '.tbi') # uses os.path.exists to see if the input file + 'tbi' extension exists. The object will be made using args.input_vcffile and args.input_vcffile will be passed to the script as an absolute file path so the path to the dir is taken care of

    def _create_sample_filter_files(self):
        self.fm_obj.downloadData(self.fm_obj.localSampleFile_v2) # download fresh SampleDatabase_v2.xlsx so we can read it in and get ecogroup information. Changed for now to get data from the grant specific file. - 2024.02.04
        self.df = pd.read_excel(self.fm_obj.localSampleFile_v2, sheet_name = 'SampleLevel') # generate df from SampleDatabase.csv
        self.df = self.df[self.df['Platform'].isin(['ILLUMINA'])].drop_duplicates(subset='SampleID') # get rid of PacBio Samples and drops duplicates in the SampleID column leaving 612 (or eventually more) samples that we can filter below
        ####### this code block is used to generate a sample file of ALL samples that are present in a given ecogroup OR from a given column from which you want to run in the analysis. NOTE: IF --sample_subet IS ON, THEN THE CorePCA COLUMN DICTATES WHICH SAMPLES ARE USED TO GENERATE THE 1ST AND 2ND PCs. THE SAMPLES TO BE PROJECT ON TO THOSE PCs ARE DICTATED BY THIS CODE BLOCK.
        analysis_ecogroups = '_'.join(str(eg).replace('_', '') for eg in self.ecogroups)
        file_version = self.in_vcf.split('/')[-1].split('.')[0] # if we run pca_maker on multiple vcf files, this allows us to ensure all plots from a particular file go to its own dir within the self.our_dir defined below

        if self.ecogroups == ['Custom']: # If samples in non-Ecogroup_PTM columns should be run, use the Custom sample flag
            while True:
                custom_column = input('Enter name of column you want to use for sample_subsetting from SampleDatabase_v2.xlsx: ')
                if custom_column not in self.df.columns.to_list():
                    print('Ivalid Column. Check spelling and re-enter a column name')
                    continue
                else:
                    break
            self.out_dir = self.out_dir + '/' + file_version + '/' + custom_column
            pathlib.Path(self.out_dir).mkdir(parents=True, exist_ok=True) # This generates the outdir if it doesn't exist so later tools don't run into errors making files.
            self.all_ecogroup_samples_csv = self.out_dir + '/all_samples_in_this_ecogroup.csv' # this file conatins all samples in the ecogroup(s) the pipeline is being run for. Samples are not subset yet in this file
            self.all_ecogroup_samples_metadata = self.out_dir + '/all_samples_in_this_ecogroup_with_metadata.csv'
            self.subset_samples_csv = self.out_dir + '/subset_samples.csv' # This will be the name of the output file containing names of samples that have the ecogroups specified.
            self.subset_samples_metadata = self.out_dir + '/subset_samples_with_metadata.csv'
            ecogroup_df = self.df[self.df[custom_column] == 'Yes']
            ecogroup_df.to_csv(self.all_ecogroup_samples_csv, columns = ['SampleID'], header=False, index=False)
            ecogroup_df.to_csv(self.all_ecogroup_samples_metadata, columns = ['SampleID', 'Ecogroup_PTM', 'Organism'], sep='\t', index=False)
        else:
            self.out_dir = self.out_dir + '/' + file_version + '/' + analysis_ecogroups
            pathlib.Path(self.out_dir).mkdir(parents=True, exist_ok=True) # This generates the outdir if it doesn't exist so later tools don't run into errors making files.
            self.all_ecogroup_samples_csv = self.out_dir + '/all_samples_in_this_ecogroup.csv' # this file conatins all samples in the ecogroup(s) the pipeline is being run for. Samples are not subset yet in this file
            self.all_ecogroup_samples_metadata = self.out_dir + '/all_samples_in_this_ecogroup_with_metadata.csv'
            self.subset_samples_csv = self.out_dir + '/subset_samples.csv' # This will be the name of the output file containing names of samples that have the ecogroups specified.
            self.subset_samples_metadata = self.out_dir + '/subset_samples_with_metadata.csv'
            # Code block to hard code in the eco group infomation and change the value of the "ecogroups" attribute depending on what the "args.ecogroups" value will be.
            if self.ecogroups == ['All']:
                self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'Riverine', 'AC']
            elif self.ecogroups == ['Non_Riverine']:
                self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon']
            elif self.ecogroups == ['Lake_Malawi']:
                self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'AC']
            elif self.ecogroups == ['Rock_Sand']:
                self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic']
            elif self.ecogroups == ['Sand']:
                self.ecogroups = ['Utaka', 'Shallow_Benthic', 'Deep_Benthic']

            ecogroup_df = self.df[self.df.Ecogroup_PTM.isin(self.ecogroups)]
            ecogroup_df.to_csv(self.all_ecogroup_samples_csv, columns = ['SampleID'], header=False, index=False)
            ecogroup_df.to_csv(self.all_ecogroup_samples_metadata, columns = ['SampleID', 'Ecogroup_PTM', 'Organism'], sep='\t', index=False)
        
        # This code block is to create the SUBSET sample file and an extra metadata file that may be useful for troubleshooting. This subset file is used by plink to generate the subset PC space that all the samples will be projected on to. Samples here are determined solely by the CorePCA column in SampleDatabase_v2
        if args.sample_subset:
            # We only want the samples in CorePCA to be used to generate the PCs. Thus, we can use the original self.df and just take the samples where CorePCA == True. These samples will be the same no matter what analysis is run, assuming the VCF file has variant data for all of these samples. 
            # NOTE: IF THE INPUT VCF FILE DOES NOT HAVE ALL OF THE SAMPLES IN CorePCA THEN DO NOT USE THE --sample_subset FLAG. PCA WILL STILL WORK, BUT THE PCs WILL NOT BE GENERATED USING THE 284 CorePCA SAMPLES. Even if this mistake is made, the erorr check at the end of this funciton will catch it.
            subset_df = self.df[self.df['CorePCA']=='Yes']
            subset_df.to_csv(self.subset_samples_csv, columns = ['SampleID'], header=False, index=False) # store the subset samples into a csv file in the directory for that ecogroup's analysis
            subset_df.to_csv(self.subset_samples_metadata, columns = ['SampleID', 'Ecogroup_PTM', 'Organism'], sep='\t', index=False) # write an extra csv file containing metadata information for samples included in the subsetting
        else: # we will not generate a subset vcf file that will be used in the analyses. Thus, self.subset_samples_csv & self.subset_samples_metadata will not be generated. 
            # To make it easier to keep the code revolving around the subset code that already exists, I can simply generate the same files but just make them equal to what exists in all_samples files.
            warning = input("WARNING: You have called on the script to run without subsetting samples.\nYour PCA results may be driven by overrepresented samples in the ecogroup you're analyzing. Are you sure you wish to continue? (Y/N): ")
            if warning.upper() == 'Y':
                print("Continuing... Note that files with 'subset' in the name will not actually be a subset of your samples. Filenames are preserved out of laziness...")
            else:
                exit()
            ecogroup_df.to_csv(self.subset_samples_csv, columns = ['SampleID'], header=False, index=False) # write the same samples that created the self.all_ecogroup_samples_csv and self.all_ecogroup_samples_metadata into the subset files so that I don't have to change code in other functions
            ecogroup_df.to_csv(self.subset_samples_metadata, columns = ['SampleID', 'Ecogroup_PTM', 'Organism'], sep='\t', index=False)

        # code that checks if the samples in self.subset_samples.csv are in self.in_vcf:
        in_vcf_samples = subprocess.check_output(['bcftools', 'query', '-l', self.in_vcf], encoding = 'utf-8').strip().split('\n')
        in_vcf_samples = {'SampleID':in_vcf_samples}
        in_vcf_samples_df = pd.DataFrame(in_vcf_samples)
        unmatched_samples = ecogroup_df[~ecogroup_df['SampleID'].isin(in_vcf_samples_df['SampleID'])]
        if unmatched_samples.shape[0] != 0: # if all samples in ecogroup_df match those in in_vcf, then the row length should be 0. If the value is different, raise an exception
            unmatched_samples.to_csv('unmatched_samples.csv', index=False)
            raise Exception(f"{unmatched_samples.shape[0]} samples in ecogroup_df are not in the input vcf file. See samples in unmatched_samples.csv.")
            

    def _create_ecogroup_specific_vcf(self):
        # path to the vcf file containing samples only in the ecogroups needed for the analysis
        self.ecogroup_specific_master_vcf = self.out_dir + '/all_ecogroup_samples.vcf.gz' # for each analysis, depending on the ecogroups chosen, we generate a vcf file with all samples of only that eco group. The path and name of that file is defined here. Old file name = samples_filtered_master.vcf.gz
        self.subset_master_vcf = self.out_dir + '/subset_samples.vcf.gz' # for each analysis, this subset vcf file will need to be generated to pull region specific variants from in the next step
        
        if pathlib.Path(self.ecogroup_specific_master_vcf).exists():
            # NOTE: DO SAMPLE ID CHECKS HERE 
            # The ecogroup_specific master file will always have the same number of samples as the all_ecogroup_samples.csv. However, depending on if the sample_subset flag is called, the subset_samples.csv and subset_samples.vcf.gz may have different numbers of samples. If they do, rerun the first level vcf file creation`.`
            if subprocess.run(f"bcftools query -l {self.subset_master_vcf}", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout == subprocess.run(f"cat {self.subset_samples_csv}", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout: # checks if the output from printing the sample names from samples_to_keep.csv and the column names from samples_filtered_master.vcf.gz are the sample
                print(f'\nThe file {self.subset_master_vcf} exists and samples within the file match those in subset_samples.csv. New all_ecogroup_samples.vcf.gz and subset_samples.vcf.gz files will not be built.\n')
            else:
                print('\nSamples in the existing subset_samples.vcf.gz file and the subset_samples.csv file are different. The exisiting all_ecogroup_samples.vcf.gz and subset_samples.vcf.gz files will be overwritten to match the samples present for this analysis\n')
                p1 = subprocess.Popen(['bcftools', 'view', self.in_vcf, '--samples-file', self.all_ecogroup_samples_csv, '--min-ac=1', '--no-update', '-o', self.ecogroup_specific_master_vcf, '-O', 'z'])
                p2 = subprocess.Popen(['bcftools', 'view', self.in_vcf, '--samples-file', self.subset_samples_csv, '-o', self.subset_master_vcf, '-O', 'z'])
                print('\nREGENERATING A VCF FILE CONTAINING ONLY SAMPLES OF THE ECOGROUPS SPECIFIED.')
                print('\nREGENERATING A VCF FILE WITH THE SUBSET SAMPLES FOR THE ECOGROUPS SPECIFIED.')
                p1.communicate()
                p2.communicate()

        else: # if no ecogroup_specific_master_vcf file exists yet, go ahead and use the samples to keep file to generate a cvf file containing only samples of the ecogroup being analyzed.
            p1 = subprocess.Popen(['bcftools', 'view', self.in_vcf, '--samples-file', self.all_ecogroup_samples_csv, '--min-ac=1', '--no-update', '-o', self.ecogroup_specific_master_vcf, '-O', 'z'])
            p2 = subprocess.Popen(['bcftools', 'view', self.in_vcf, '--samples-file', self.subset_samples_csv, '-o', self.subset_master_vcf, '-O', 'z'])
            print('\nGENERATING A VCF FILE CONTAINING ONLY SAMPLES OF THE ECOGROUPS SPECIFIED.')
            print('\nGENERATING A VCF FILE WITH THE SUBSET SAMPLES FOR THE ECOGROUPS SPECIFIED.')
            p1.communicate()
            p2.communicate()

    def _preprocess_vcfs(self): # this function prunes all variants from self.ecogroup_specific_master_vcf that are not present in self.subset_master_vcf following recalculation and filtering of variants with AF < 0.05. The 5% threshold was used in the Malinksy E&E papee according to Patrick Unverified by me
        # File check code to ensure these files don't already exist. If they do, the samples between existing files and the cohort's samples must match:
        if pathlib.Path(self.out_dir + '/af_recalculated_subset_samples.vcf.gz').exists():
            # if the ecogroup_specific_master_vcf (contains all variants per sample for the ecogroups specified) exists, then this checks that the samples match exactly. If not, a new file is built by filtering for samples in the self.good_samples_csv file.
            if subprocess.run(f"bcftools query -l {self.out_dir + '/af_recalculated_subset_samples.vcf.gz'}", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout == subprocess.run(f"cat {self.out_dir + '/subset_samples.csv'}", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout: # checks if the output from printing the sample names from samples_to_keep.csv and the column names from samples_filtered_master.vcf.gz are the sample
                print(f'\nThe file af_recalculated_subset_samples.vcf.gz exists in the ecogroup dir and samples within the file match those in subset_samples.csv. VCF preprocessing will be skipped\n')
            else: # if the samples differ between the existing af_recalculated_subset_samples.vcf.gz and the subset_samples.csv, redo all preprocessing steps. 
                print('\nSamples between the existing af_recalculated_subset_samples.vcf.gz and subset_samples.csv are different. Rerunning VCF preprocessing...\n')
                print('RESTARTING VCF FILE PREPROCESSING...\n')
                print('RECALCULATING AF IN THE SUBSET FILE\n')
                subprocess.run(['bcftools', '+fill-tags', self.subset_master_vcf, '-Oz', '-o', self.out_dir + '/af_recalculated_subset_samples.vcf.gz', '--', '-t', 'AF'])
                print('FILTERING OUT VARIANTS WITH AF < 0.05\n')
                subprocess.run(['bcftools', 'view', '--exclude', 'AF<0.05 || AF == 1', self.out_dir + '/af_recalculated_subset_samples.vcf.gz', '-o', self.out_dir + '/af_filtered_subset_samples.vcf.gz', '-Oz'])
                print('WRITING CHRMOSOME-POSITION FILE\n')
                subprocess.run(['bcftools', 'query', '-f', '%CHROM\t%POS\n', self.out_dir + '/af_filtered_subset_samples.vcf.gz', '-o', self.out_dir + '/variants_to_keep.txt'])
                print('PRUNING VARIANTS FROM ECGOGROUP VCF FILE USING THOSE PRESENT IN THE SUBSET\n')
                subprocess.run(['vcftools', '--positions', self.out_dir + '/variants_to_keep.txt', '--gzvcf', self.ecogroup_specific_master_vcf, '--recode', '--out', self.out_dir + '/variants_filtered_ecogroup_samples'])
                # the outputs that go into the plink code from here to generate the pfiles are self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf' & self.out_dir + '/af_filtered_subset_samples.vcf.gz'
                print('ZIPPING THE VARIANTS FILTERED ECOGROUP FILE\n')
                subprocess.run(['bgzip', '-f', self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf'])

                # This code needs to be incorporated at the vcf preprocessing stage since the outputs at the end of the file need to get indexed
                p3 = subprocess.Popen(['tabix', '-p', 'vcf', self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf.gz'])
                p4 = subprocess.Popen(['tabix', '-p', 'vcf', self.out_dir + '/af_filtered_subset_samples.vcf.gz'])
                print('FILTERED SAMPLES FILE GENERATED. INDEXING ALL SAMPLE VCF FILE...')
                print('INDEXING SUBSET SAMPLE VCF FILE...')
                p3.communicate()
                p4.communicate()
                print('INDEX CREATED...\nVCF PREPROCESSING SUCCESSS...\n')

        else: # if no ecogroup_specific_master_vcf file exists yet, go ahead and use the samples to keep file to generate a cvf file containing only samples of the ecogroup being analyzed.
            print('STARTING VCF FILE PREPROCESSING...\n')
            print('RECALCULATING AF IN THE SUBSET FILE\n')
            subprocess.run(['bcftools', '+fill-tags', self.subset_master_vcf, '-Oz', '-o', self.out_dir + '/af_recalculated_subset_samples.vcf.gz', '--', '-t', 'AF'])
            print('FILTERING OUT VARIANTS WITH AF < 0.05\n')
            subprocess.run(['bcftools', 'view', '--exclude', 'AF<0.05 || AF == 1', self.out_dir + '/af_recalculated_subset_samples.vcf.gz', '-o', self.out_dir + '/af_filtered_subset_samples.vcf.gz', '-Oz'])
            print('WRITING CHRMOSOME-POSITION FILE\n')
            subprocess.run(['bcftools', 'query', '-f', '%CHROM\t%POS\n', self.out_dir + '/af_filtered_subset_samples.vcf.gz', '-o', self.out_dir + '/variants_to_keep.txt'])
            print('PRUNING VARIANTS FROM ECGOGROUP VCF FILE USING THOSE PRESENT IN THE SUBSET\n')
            subprocess.run(['vcftools', '--positions', self.out_dir + '/variants_to_keep.txt', '--gzvcf', self.ecogroup_specific_master_vcf, '--recode', '--out', self.out_dir + '/variants_filtered_ecogroup_samples'])
            # the outputs that go into the plink code from here to generate the pfiles are self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf' & self.out_dir + '/af_filtered_subset_samples.vcf.gz'
            print('ZIPPING THE VARIANTS FILTERED ECOGROUP FILE\n')
            subprocess.run(['bgzip', '-f', self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf'])

            # This code needs to be incorporated at the vcf preprocessing stage since the outputs at the end of the file need to get indexed 
            p3 = subprocess.Popen(['tabix', '-p', 'vcf', self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf.gz'])
            p4 = subprocess.Popen(['tabix', '-p', 'vcf', self.out_dir + '/af_filtered_subset_samples.vcf.gz'])
            print('FILTERED SAMPLES FILE GENERATED. INDEXING ALL SAMPLE VCF FILE...')
            print('INDEXING SUBSET SAMPLE VCF FILE...')
            p3.communicate()
            p4.communicate()
            print('INDEX CREATED...\nVCF PREPROCESSING SUCCESSS...\n')

    def _split_VCF_to_LG(self, linkage_group_list):
        processes1 = []
        processes2 = []
        for lg in linkage_group_list:
            if lg == 'Whole':
                continue
            elif lg in self.exploratory_regions_list:
                self._create_exploratory_region_eigen_files()
            elif lg in self.linkage_group_map.values():
                pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True) # generate the file paths to he split LG dirs within a dir named "PCA"
                if pathlib.Path(self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset.vcf.gz').exists(): # For each linkage groups' dir in the PCA dir, a subset_sample.vcf.gz exists, then check that the samples in that file are the same as what's in the subset_samples.csv file 
                    if subprocess.run(f"bcftools query -l {self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset.vcf.gz'}", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout == subprocess.run(f"cat {self.subset_samples_csv}", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout: # checks if the subset vcf file that already exists has matching sample names as those in the samples_to_keep.csv file.
                        print('\nThe file ' + lg + '_sample_subset.vcf.gz exists and has the same samples as subset_sample.csv. New VCF files for ' + lg + ' will not be generated.\n')
                    else:
                        print('The sample_subset.vcf.gz file for ' + lg + ' has samples that does not match with those found in subset_samples.csv. New whole_ecogroup and sample_subset vcf files will be generated.\n')
                        print('Generating a VCF file containing the subset samples for ' + lg + '...')
                        p1 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf.gz', '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_ecogroup.vcf.gz', '-O', 'z']) # takes in variants filtered ecogroup_specific_master_vcf and filters out each LG and writes the file into appropriate PCA dir.
                        p2 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.out_dir + '/af_filtered_subset_samples.vcf.gz', '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset.vcf.gz', '-O', 'z']) # takes in the subset_samples file with the remaining variants with AF > 0.05 ecogroup_specific_master_vcf and filters out each LG and writes the file into appropriate PCA dir.
                        processes1.append(p1)
                        processes2.append(p2)
                        if len(processes1) == len(linkage_group_list): # parallelization code
                            for proc1 in processes1:
                                proc1.communicate()
                            processes1 = []
                        if len(processes2) == len(linkage_group_list):
                            for proc2 in processes2:
                                proc2.communicate()
                else:
                    print('Generating a VCF file containing the subset samples for ' + lg + '...')
                    p1 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf.gz', '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_ecogroup.vcf.gz', '-O', 'z']) # takes in variants filtered ecogroup_specific_master_vcf and filters out each LG and writes the file into appropriate PCA dir.
                    p2 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.out_dir + '/af_filtered_subset_samples.vcf.gz', '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset.vcf.gz', '-O', 'z']) # takes in the subset_samples file with the remaining variants with AF > 0.05 ecogroup_specific_master_vcf and filters out each LG and writes the file into appropriate PCA dir.
                    processes1.append(p1)
                    processes2.append(p2)
                    if len(processes1) == len(linkage_group_list): # parallelization code
                        for proc1 in processes1:
                            proc1.communicate()
                        processes1 = []
                    if len(processes2) == len(linkage_group_list):
                        for proc2 in processes2:
                            proc2.communicate()

    def _create_exploratory_region_eigen_files(self):
        # changes list to include a split lg11 and lg20 interval. The middle value for LG11 is 18425216 and the middle value for lg20 is 29716509. These values assume the Mzebra_GT3 genome.
        #  also expanded the list to equally split the lg9 inversion into 5 equivalent regions. This will be run alongside the whole lg9 inverted region and the whole of lg9
        
        inversion_regions = self.exploratory_regions_list
        exploratory_regions_list = []
        # if only one or two exploratory regions are passed, below code allows only those select few to be run instead of all of them every time. 
        for region in self.linkage_groups:
            if region in inversion_regions:
                exploratory_regions_list.append(region)
        processes = []

        for region in exploratory_regions_list:
            pathlib.Path(self.out_dir + '/PCA/' + region + '/').mkdir(parents=True, exist_ok=True) # make the filepath to the exploratory region output dirs
            if pathlib.Path(self.out_dir + '/PCA/' + region + '/' + region + '_sample_subset.vcf.gz').exists(): # if one of the files already exists, test if the number of sample is the same as in the subset_csv file in self.out_dir
                if subprocess.run(f"bcftools query -l {self.out_dir + '/PCA/' + region + '/' + region + '_sample_subset.vcf.gz'}", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout == subprocess.run(f"cat {self.subset_samples_csv}", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout: # checks if an existing exploratory region vcf file contains the same samples as samples_to_keep.csv. If not, make a new one
                    print('The file ' + region + '_sample_subset.vcf.gz exists and has the same samples as subset_sample.csv. New VCF files for ' + region + ' will not be generated.\n')
                else:
                    print('The sample_subset.vcf.gz file for ' + region + ' has samples that does not match with those found in subset_samples.csv. New whole_ecogroup and sample_subset vcf files will be generated.\n')
                    print('Regenerating subset and whole_ecogroup vcf files for the regions in ' + region + '...\n')
                    p1 = subprocess.Popen(['bcftools', 'view', self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf.gz', '-R', os.getcwd() + '/special_intervals/' + region + '.csv', '-o', self.out_dir + '/PCA/' + region + '/' + region + '_whole_ecogroup.vcf.gz', '-Oz'])
                    p2 = subprocess.Popen(['bcftools', 'view', self.out_dir + '/af_filtered_subset_samples.vcf.gz', '-R', os.getcwd() + '/special_intervals/' + region + '.csv', '-o', self.out_dir + '/PCA/' + region + '/' + region + '_sample_subset.vcf.gz', '-Oz'])
                    processes.append(p1)
                    processes.append(p2)
                    if len(processes) == 2*len(exploratory_regions_list): # parallelization code
                        for proc in processes:
                            proc.communicate()
                        processes = []
            else: # if the files don't exist in the first place, go alead and generate each region's sample_subset.vcf.gz and _whole_ecogroup.vcf.gz files. 
                # Generate 2 vcf files... a subset file one using the regions and a filtered ecogroup one using the regions. Do it all in parallel eventually
                print('Generating subset and whole_ecogroup vcf files for the regions in ' + region + '...\n')
                p1 = subprocess.Popen(['bcftools', 'view', self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf.gz', '-R', os.getcwd() + '/special_intervals/' + region + '.csv', '-o', self.out_dir + '/PCA/' + region + '/' + region + '_whole_ecogroup.vcf.gz', '-Oz'])
                p2 = subprocess.Popen(['bcftools', 'view', self.out_dir + '/af_filtered_subset_samples.vcf.gz', '-R', os.getcwd() + '/special_intervals/' + region + '.csv', '-o', self.out_dir + '/PCA/' + region + '/' + region + '_sample_subset.vcf.gz', '-Oz'])
                processes.append(p1)
                processes.append(p2)
                if len(processes) == 2*len(exploratory_regions_list): # parallelization code
                    for proc in processes:
                        proc.communicate()
                    processes = []

    def _create_eigenfiles_per_LG(self, linkage_group_list): # new hidden method that will create PCA plots for each LG in sample. It will define attributes for the object and also takes in a lingage grouup. Calling on this method in a for lopp should generate the eigenvalue/vector files needed per lg in self.contigs
        for lg in linkage_group_list:
            if lg == 'Whole':
                pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True)
                subprocess.run(['plink2', '--vcf', self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf.gz', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_generated_step_1', '--allow-extra-chr']) # whole vcf file pfile generation. Use the variants_filtered_ecogroup_samples.recode.vcf.gz file as the "whole" file.
                subprocess.run(['plink2', '--vcf', self.out_dir + '/af_filtered_subset_samples.vcf.gz', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_generated_step_1', '--allow-extra-chr']) # subset sample vcf file pfile generation. The subset file is the af_filtered_subset_samples.vcf.gz in self.out_dir
                subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_generated_step_1', '--set-missing-var-ids', '@:#', '--make-pgen', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected_generated_step_2', '--allow-extra-chr'])
                if int(subprocess.check_output(f"bcftools query -l {self.out_dir + '/af_filtered_subset_samples.vcf.gz'} | wc -l", shell=True, encoding='utf-8').strip()) >= 50: # NOTE: Since for "Whole", the subset samples are pulled from af_filtered_subset_samples.vcf.gz, this file must be queried to see if it contains less than 50 samples
                    # For running the whole genome, LD pruning is indeed needed
                    # Run linkage pruning only to generate the prune.in file using indep-pairwise 50 5 0.1
                    subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_generated_step_1', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_no_pca_generation_yes_ld_pruning_generated_step_3.1', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '5', '0.1'])
                    # use the variants in the pfiles from subset samples generated in step 1 and load in the prune.in varinats only and run pca. Generates the .account and .allele files needed for step4
                                    
                    subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_generated_step_1', '--freq', 'counts', '--pca', 'allele-wts', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--max-alleles', '2', '--extract', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_no_pca_generation_yes_ld_pruning_generated_step_3.1.prune.in'])
                    # modify the .acounts file to eliminate 0 count alleles:
                    acounts_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2.acount', sep='\t')
                    no_extremes_df = acounts_df[(acounts_df['ALT_CTS'] != 0) & (acounts_df['ALT_CTS'] != max(acounts_df['OBS_CT']))]
                    no_extremes_df.to_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2.acount', sep='\t', index=False)

                    # genenrate the .sscore file with eigenvectors for each samples after projection on to the PC space from the subset PCA
                    subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected_generated_step_2', '--read-freq', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2.acount', '--score', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2.eigenvec.allele', '2', '5', 'header-read', 'no-mean-imputation', 'variance-standardize', '--score-col-nums', '6-15', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection_yes_ld', '--allow-extra-chr'])
                # else: # code for mbuna and samples with <50 samples. Ignore for now.
                #     subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset', '--freq', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '5', '0.1', '--max-alleles', '2', '--bad-ld']) 
                #     subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected', '--pca', 'allele-wts', '--read-freq', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.afreq', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca', '--allow-extra-chr']) # this is where things are going wrong. Before, "Whole" was always running this faulty code that's generating .eigenvec files with a number of samples equal to the whole ecogroup and not the subset This needs to  be fixed, else Mbuna and other small dataset groups will not run right.
                #     subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected', '--read-freq', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.afreq', '--score', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.eigenvec.allele', '2', '5', 'header-read', 'no-mean-imputation', 'variance-standardize', '--score-col-nums', '6-15', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection', '--allow-extra-chr'])

            else: #lg in self.linkage_group_map.values():
                pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True)
                if not pathlib.Path(self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_ecogroup.vcf.gz').exists(): # NOTE: error checking...
                    print('ERROR: THE FILE ' + lg + '.VCF.GZ DOES NOT EXIST. MUST RUN _SPLIT_VCF_TO_LG TO CREATE IT...')
                    raise Exception
                else: # do all the plink pca magic here assuming you're going one linkage group at a time... figure out parallelization as I code this referencing the previous function. Parallelization may not be possible because plink likes to execute immediately and doesn't listen to the Popen constructor. It all executes before Popen.communicate() is called.
                    pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True) # make the path to a PCA directory in self.out_dir if it doesn't exist yet

                    # regardless of whether the region will be linkage pruned or not, the first 3 steps are the same.
                    subprocess.run(['plink2', '--vcf', self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf.gz', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_generated_step_1', '--allow-extra-chr']) # whole vcf file pfile generation for a given lg
                    subprocess.run(['plink2', '--vcf', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset.vcf.gz', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_generated_step_1', '--allow-extra-chr']) # subset sample vcf file pfile generation for a given lg
                    subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_generated_step_1', '--set-missing-var-ids', '@:#', '--make-pgen', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected_generated_step_2', '--allow-extra-chr'])
                    # If number of samples in the analysis is > 50: run code as normal
                    if int(subprocess.check_output(f"bcftools query -l {self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset.vcf.gz '} | wc -l", shell=True, encoding='utf-8').strip()) >= 50:
                        # Run the NO LD case, which is when regions in the exploratory regions list (that include ethe inversions) are run
                        if lg in self.exploratory_regions_list: # if lf is in an inverted region, we should not be doing linkage pruning
                            subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_generated_step_1', '--freq', 'counts', '--pca', 'allele-wts', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_generation_no_ld_pruning_generated_step_3', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--max-alleles', '2'])
                            # modify the .acounts file to eliminate 0 count alleles for the NO LD case:
                            acounts_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_generation_no_ld_pruning_generated_step_3.acount', sep='\t')
                            no_extremes_df = acounts_df[(acounts_df['ALT_CTS'] != 0) & (acounts_df['ALT_CTS'] != max(acounts_df['OBS_CT']))]
                            no_extremes_df.to_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_generation_no_ld_pruning_generated_step_3.acount', sep='\t', index=False)
                            # Generate the .sscore matrix by projecting the variants in the whole_corrected pfiles from step 2 on to the PCs generated by the subset samples 
                            subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected_generated_step_2', '--read-freq', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_generation_no_ld_pruning_generated_step_3.acount', '--score', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_generation_no_ld_pruning_generated_step_3.eigenvec.allele', '2', '5', 'header-read', 'no-mean-imputation', 'variance-standardize', '--score-col-nums', '6-15', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection_no_ld', '--allow-extra-chr'])
                        else: # if the lg is a standard chromosome, run the YES LD case
                            # Run linkage pruning only to generate the prune.in file using indep-pairwise 50 5 0.1 
                            subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_generated_step_1', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_no_pca_generation_yes_ld_pruning_generated_step_3.1', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '5', '0.1'])
                            # use the variants in the pfiles from subset samples generated in step 1 and load in the prune.in varinats only and run pca. Generates the .account and .allele files needed for step4
                            subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_generated_step_1', '--freq', 'counts', '--pca', 'allele-wts', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--max-alleles', '2', '--extract', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_no_pca_generation_yes_ld_pruning_generated_step_3.1.prune.in'])

                            # modify the .acounts file to eliminate 0 count alleles:
                            acounts_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2.acount', sep='\t')
                            no_extremes_df = acounts_df[(acounts_df['ALT_CTS'] != 0) & (acounts_df['ALT_CTS'] != max(acounts_df['OBS_CT']))]
                            no_extremes_df.to_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2.acount', sep='\t', index=False)

                            # Generate the .sscore matrix by projecting the variants in the whole_corrected pfiles from step 2 on to the PCs generated by the subset samples 
                            subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected_generated_step_2', '--read-freq', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2.acount', '--score', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2.eigenvec.allele', '2', '5', 'header-read', 'no-mean-imputation', 'variance-standardize', '--score-col-nums', '6-15', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection_yes_ld', '--allow-extra-chr'])
                    # else: # code for mbuna and samples with <50 samples. Ignore for now.
                    #     subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset', '--freq', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '5', '0.1', '--max-alleles', '2', '--bad-ld'])
                    #     subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected', '--pca', 'allele-wts', '--read-freq', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.afreq', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca', '--allow-extra-chr'])
                    #     subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected', '--read-freq', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.afreq', '--score', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.eigenvec.allele', '2', '5', 'header-read', 'no-mean-imputation', 'variance-standardize', '--score-col-nums', '6-15', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection', '--allow-extra-chr'])

    def _create_interactive_pca(self, linkage_group_list): # uses plotly to generate interactive PCA html outputs
        """
        TODO:
        - make it output some pdf versions too to ensure color isn't changing like crazy 
        - remove values on axes
        - remove default background color. make transparent
        """

        self.plotly_out = self.out_dir + '/interactive_PCA_outputs/' # define outdir for final plotly outputs
        self.plotly_pdf_out = self.out_dir + '/static_PCA_outputs/' # define PDF output dir
        pathlib.Path(self.plotly_out).mkdir(parents=True, exist_ok=True) # build the file paths with pathlib.Path
        pathlib.Path(self.plotly_pdf_out).mkdir(parents=True, exist_ok=True) # build pdf output filepath 
        # The below color hexcodes match what's in the Malinksy paper, extracted from Illustrator using the PDF of the publication 
        malinksy_color_map = {'Mbuna': '#A020F0', 'AC': '#A2CD5A', 'Shallow_Benthic': '#FF6347', 'Deep_Benthic': '#4876FF', 'Rhamphochromis': '#8B4513', 'Diplotaxodon': '#FFA54F', 'Utaka': '#006400'}
        bionano_shape_map = {'No': 'circle', 'Yes': 'x'}
        # try to get "x-thin" working while maintaining its line width and other attributes

        for lg in linkage_group_list:
            print('GENERATING PCA FOR ' + lg)
            # calculate percent variance explained by pc1 and 2. Round to 2 decimals
            variance_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.eigenval', header=None)
            pc1_variance = (variance_df.loc[0][0] / variance_df.sum())[0]*100
            pc2_variance = (variance_df.loc[1][0] / variance_df.sum())[0]*100
            pc1_variance = round(pc1_variance, 2)
            pc2_variance = round(pc2_variance, 2)
            eigen_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection.sscore', sep='\t')
            eigen_df = eigen_df.rename(columns = {'#IID':'SampleID'})
            df_merged = pd.merge(eigen_df, self.df, on=['SampleID'])
            df_merged['Color'] = df_merged['Ecogroup_PTM'].map(malinksy_color_map)  # Map Ecogroup_PTM to the malinsky_color_map
            df_merged.loc[df_merged['BionanoData'] == 'Yes', 'Color'] = 'black'  # Override to black for BionanoData 'Yes'
            df_merged['PDF_Size'] = df_merged['BionanoData'].apply(lambda x: 4 if x == 'Yes' else 3)
            df_merged['HTML_Size'] = df_merged['BionanoData'].apply(lambda x: 7 if x == 'Yes' else 5)
            df_merged['Opacity'] = df_merged['BionanoData'].apply(lambda x: 0.6 if x == 'Yes' else 0.9)
            if lg.startswith('NC'):
                plot_title = list(self.linkage_group_map.keys())[list(self.linkage_group_map.values()).index(lg)]
            else:
                plot_title = lg

            # Do the plotting magic
            # basically just pass in plotly.express.scatter a dataframe, tell it what the x and y axes are, and labels for the x and y axes. The plot title and hover data for the interactive plot can be overlayed after
            pdf_fig = px.scatter(df_merged, x='PC1_AVG', y='PC2_AVG', width = 330, height = 330,
                            labels={
                                'PC1_AVG': 'PC1 ' + str(pc1_variance) + '%',
                                'PC2_AVG': 'PC2 ' + str(pc2_variance) + '%'
                            },
                            # symbol_map=bionano_shape_map, # if you have a shape_map, uncomment and add that here
                            title=plot_title, hover_data=['SampleID', 'Ecogroup_PTM', 'Organism', 'ProjectID_PTM'])

            # Apply custom colors and marker shapes
            pdf_fig.update_traces(marker=dict(color=df_merged['Color'])) # maps color based on the Color column in df_merged
            pdf_fig.update_traces(marker=dict(symbol=df_merged['BionanoData'].map(bionano_shape_map))) # maps the X's for BionanoData samples
            pdf_fig.update_traces(marker=dict(line_width=0))
            pdf_fig.update_traces(marker=dict(size=df_merged['PDF_Size'])) # maps the size we want for the points based on the "Size" column.
            pdf_fig.update_traces(marker=dict(opacity=df_merged['Opacity']))
            pdf_fig.update_layout(margin=dict(l=20, r=20, t=30, b=20))
            
            # write the interacitve and static files (adjusted for publication)
            pdf_fig.write_image(self.plotly_pdf_out + lg + '_PCA.pdf')

            # Generate HTML figures with default sizes
            html_fig = px.scatter(df_merged, x='PC1_AVG', y='PC2_AVG',
                            labels={
                                'PC1_AVG': 'PC1 ' + str(pc1_variance) + '%',
                                'PC2_AVG': 'PC2 ' + str(pc2_variance) + '%'
                            },
                            # symbol_map=bionano_shape_map, # if you have a shape_map, uncomment and add that here
                            title=plot_title, hover_data=['SampleID', 'Ecogroup_PTM', 'Organism', 'ProjectID_PTM'])

            # Apply custom colors and marker shapes
            html_fig.update_traces(marker=dict(color=df_merged['Color'])) # maps color based on the Color column in df_merged
            html_fig.update_traces(marker=dict(symbol=df_merged['BionanoData'].map(bionano_shape_map))) # maps the X's for BionanoData samples
            html_fig.update_traces(marker=dict(line_width=0)) # make the white outline width 0
            html_fig.update_traces(marker=dict(size=df_merged['HTML_Size'])) # maps the size we want for the points based on the "Size" column.
            html_fig.update_traces(marker=dict(opacity=df_merged['Opacity']))

            # write the interacitve files (full_sized)
            html_fig.write_html(self.plotly_out + lg + '_PCA.html')


    def _create_umap(self, linkage_group_list):
        # code to generate and merge the sampledatabase_df and the eigen_df
        self.umap_out = self.out_dir + '/umap_outputs/'
        pathlib.Path(self.umap_out).mkdir(parents=True, exist_ok=True) # build a file path to the umap out dir
        color_map = {'Mbuna': 'purple', 'AC': 'limegreen', 'Shallow_Benthic': 'red', 'Deep_zBenthic': 'blue', 'Rhamphochromis': 'brown', 'Diplotaxodon': 'orange', 'Utaka': 'darkgreen', 'Hybrid': 'darkpink'}
        shape_map = {'MalinskyData': 'square', 'Streelman_McGrathData': 'diamond', 'BrainDiversity_s1': 'star', 'MC_males': 'circle', 'MC_females': 'circle-open'}
        for lg in linkage_group_list:
            print('GENERATING UMAP FOR ' + lg)
            eigen_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection.sscore', sep='\t')
            eigen_df = eigen_df.rename(columns = {'#IID':'SampleID'})
            df_merged = pd.merge(eigen_df, self.df, on=['SampleID'])
            ten_pcs = eigen_df.drop(['SampleID', 'ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM'], axis=1)
            ten_pcs_np_array = ten_pcs.to_numpy()
            # code to use the 10PCs to generate 2 columns of transformed data
            fit = umap.UMAP()
            umap_transformed_data = fit.fit_transform(ten_pcs_np_array)
            if lg.startswith('NC'):
                plot_title = list(self.linkage_group_map.keys())[list(self.linkage_group_map.values()).index(lg)]
            else:
                plot_title = lg
            fig = px.scatter(df_merged, x=umap_transformed_data[:,0], y=umap_transformed_data[:,1], color='Ecogroup_PTM', symbol='ProjectID_PTM', color_discrete_map=color_map, symbol_map=shape_map, title=plot_title, hover_data=['SampleID', 'Ecogroup_PTM', 'Organism', 'ProjectID_PTM'])
            larger_size = 9
            smaller_size = 3
            fig.update_traces(marker=dict(size=smaller_size), selector=dict(marker_symbol='square'))
            fig.update_traces(marker=dict(size=smaller_size), selector=dict(marker_symbol='diamond'))
            fig.update_traces(marker=dict(size=larger_size), selector=dict(marker_symbol='circle'))
            fig.update_traces(marker=dict(size=larger_size), selector=dict(marker_symbol='circle-open'))
            fig.update_traces(marker=dict(size=larger_size), selector=dict(marker_symbol='star'))
            fig.write_html(self.umap_out + lg + '_UMAP.html')

    def create_PCA(self):
        # Order of hidden methods to perform the analysis
        self._create_sample_filter_files()
        self._create_ecogroup_specific_vcf()
        self._preprocess_vcfs()
        self._split_VCF_to_LG(self.linkage_groups)
        if args.plink:
            try:
                self._create_eigenfiles_per_LG(self.linkage_groups)
            except:
                pass
        self._create_interactive_pca(self.linkage_groups)
        # if args.umap:
        #     self._create_umap(self.linkage_groups)

if __name__ == "__main__":
    pca_obj = PCA_Maker(args.genome, args.ecogroups, args.regions, args.output_dir)
    pca_obj.create_PCA()
    print('PIPELINE RUN SUCCESSFUL')

"""
Code used to test the new commands for PCA
Ok so generally here's how it works rn:
1. pfiles are generated for the subset and the whole cohort vcf files 
plink2 --vcf /Users/kmnike/Data/pca_testing/local_testing_498_master_file/PCAFigure/variants_filtered_ecogroup_samples.recode.vcf.gz --out test_whole_generated_step_1 --allow-extra-chr
plink2 --vcf lg_sample_subset.vcf.gz --out test_subset_generated_step_1 --allow-extra-chr

2. some missing IDs are corrected in the whole file  No read changes are made by the whole_corrected suite of pfiles is made
plink2 --pfile test_whole_generated_step_1 --set-missing-var-ids @:# --make-pgen --out test_whole_corrected_generated_step_2 --allow-extra-chr


3. The subset sampels are used to generate a PCA. This is where the issue is rn  I need to start testing here 
    - The indep-pairside only needs to be called and applied in the following scenarios:
        - The whole genome PCA is done
        - Independent linkeg groups 
    - The indep-pairwise pruned samples should not be called if we are running PCA on the variants within inverted regions
        - I also think that linkage pruning should not be done if the inversion seems to be fixed (lg13 lg20)

######## If No LD will be done (like for inverted regions) here's the one command used:
#### below is the original pca calling code. This was adjusted to get the below line of code that's run when no ld pruning is needed
plink2 --pfile _subset --freq counts --pca allele-wts --out _sample_subset_pca --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 5 0.1 --max-alleles 2
####
plink2 --pfile test_subset_generated_step_1 --freq counts --pca allele-wts --out test_subset_yes_pca_generation_no_ld_pruning_generated_step_3 --allow-extra-chr --set-missing-var-ids @:# --max-alleles 2

######## If LD pruning is needed, here's the 2 sommands used
plink2 --pfile test_subset_generated_step_1 --out test_subset_no_pca_generation_yes_ld_pruning_generated_step_3.1 --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 5 0.1
plink2 --pfile test_subset_generated_step_1 --freq counts --pca allele-wts --out test_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2 --allow-extra-chr --set-missing-var-ids @:# --max-alleles 2 --extract test_subset_no_pca_generation_yes_ld_pruning_generated_step_3.1.prune.in

4. The samples in the whole_corrected are projected on to the PCs created from the subset sampels
####### If running NOT LD pruned (inverted regions):
plink2 --pfile test_whole_corrected_generated_step_2 --read-freq test_subset_yes_pca_generation_no_ld_pruning_generated_step_3.acount --score test_subset_yes_pca_generation_no_ld_pruning_generated_step_3.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15 --out test_new_projection_no_ld --allow-extra-chr

####### If running LD pruned (whole genome, individual chromosomes):
plink2 --pfile test_whole_corrected_generated_step_2 --read-freq test_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2.acount --score test_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15 --out test_new_projection_yes_ld --allow-extra-chr



So here's the summary of what PCA analyses that need to be run:
For all inverted regions:
    - The LD pruning of variants will not be performed.
For all whole genome and individual LGs:
    - The LD pruning based on 50 5 0.1 will be performed and used.

Ok so here's what I need to test:
    - The whole PCA can be the start point. I need to run the plink commands on this file but use the prune.in variants only.

Sumamry of commands:
STEP1: BOTH SCENARIOS
plink2 --vcf /Users/kmnike/Data/pca_testing/local_testing_498_master_file/PCAFigure/variants_filtered_ecogroup_samples.recode.vcf.gz --out test_whole_generated_step_1 --allow-extra-chr
plink2 --vcf lg_sample_subset.vcf.gz --out test_subset_generated_step_1 --allow-extra-chr

STEP2: BOTH SCENARIOS
plink2 --pfile test_whole_generated_step_1 --set-missing-var-ids @:# --make-pgen --out test_whole_corrected_generated_step_2 --allow-extra-chr

STEP3: NO LD (EXPLORATORY REGIONS ONLY)
plink2 --pfile test_subset_generated_step_1 --freq counts --pca allele-wts --out test_subset_yes_pca_generation_no_ld_pruning_generated_step_3 --allow-extra-chr --set-missing-var-ids @:# --max-alleles 2
STEP3: YES LD (INDIVIDUAL CHROMOSOMES AND WHOLE GENOME)
plink2 --pfile test_subset_generated_step_1 --out test_subset_no_pca_generation_yes_ld_pruning_generated_step_3.1 --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 5 0.1
plink2 --pfile test_subset_generated_step_1 --freq counts --pca allele-wts --out test_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2 --allow-extra-chr --set-missing-var-ids @:# --max-alleles 2 --extract test_subset_no_pca_generation_yes_ld_pruning_generated_step_3.1.prune.in

STEP4: NO LD (EXPLORATORY REGIONS ONLY)
plink2 --pfile test_whole_corrected_generated_step_2 --read-freq test_subset_yes_pca_generation_no_ld_pruning_generated_step_3.acount --score test_subset_yes_pca_generation_no_ld_pruning_generated_step_3.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15 --out test_new_projection_no_ld --allow-extra-chr
STEP4: YES LD (INDIVIDUAL CHROMOSOMES AND WHOLE GENOME)
plink2 --pfile test_whole_corrected_generated_step_2 --read-freq test_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2.acount --score test_subset_yes_pca_geration_yes_ld_pruning_generated_step_3.2.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15 --out test_new_projection_yes_ld --allow-extra-chr



"""



"""
Legacy Plotly Code in case it's ever needed
Started 2024.09.20

color_map = {'Mbuna': 'purple', 'AC': 'limegreen', 'Shallow_Benthic': 'red', 'Deep_Benthic': 'blue', 'Rhamphochromis': 'brown', 'Diplotaxodon': 'orange', 'Utaka': 'darkgreen', 'Riverine': 'pink'}
project_ID_shape_map = {'MalinskyData': 'square', 'Streelman_McGrathData': 'diamond', 'BrainDiversity_s1': 'star', 'MC_males': 'circle', 'MC_females': 'circle-open'} # removed for now to exclude shapes when generating data for Patrick's grant.

# if you need to select a specific symbol and adjust an attribute do it like this:
fig.update_traces(marker=dict(width=1), selector=dict(symbol='x')) 
# the marker=dict() will allow you to edit settings of marker level attributes. The sleector=dict(symbol='') tells plotly which symbol you want the settings on the right to ghet applied to. 


Update code to output full HTMLs but keep small PDFs to maintain spacing.
"""
##############################################################################################################################################################################################################################################

"""
For local testing:
time python pca_maker.py Mzebra_GT3 /Users/kmnike/Data/pca_testing -e Custom --local_test --sample_subset --plink
PCAFigure

For running on Utaka:
NOTE: we are no longer running the "All" & "Non_Riverine" groups

time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -p -r All Whole Exploratory -e Lake_Malawi 2> pca_logs/error_lm_240713.txt 1> pca_logs/log_lm_24713.txt
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -p -r All Whole Exploratory -e Rock_Sand 2> pca_logs/error_rock_sand_240713.txt 1> pca_logs/log_rock_sand_240713.txt
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -p -r All Whole Exploratory -e Sand 2> pca_logs/error_sand_240713.txt 1> pca_logs/log_sand_240713.txt
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -p -r All Whole Exploratory -e Mbuna 2> pca_logs/error_mbuna_240713.txt 1> pca_logs/log_mbuna_240713.txt

time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -p -r All Whole Exploratory -e Core_and_SD 2> pca_logs/error_Core_and_SD_240717.txt 1> pca_logs/log_Core_and_SD_240717.txt
Local_Testing_Code:
python pca_maker.py Mzebra_GT3 /Users/kmnike/Data/pca_testing --sample_subset -p -r All Whole Exploratory -e Core_and_SD --local_test



time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r All Whole Exploratory -e Lake_Malawi 2> pca_logs/error_lm_240918.txt 1> pca_logs/log_lm_240918.txt
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r All Whole Exploratory -e Rock_Sand 2> pca_logs/error_rock_sand_240918.txt 1> pca_logs/log_rock_sand_240918.txt
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r All Whole Exploratory -e Sand 2> pca_logs/error_sand_240918.txt 1> pca_logs/log_sand_240918.txt
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r All Whole Exploratory -e Mbuna 2> pca_logs/error_mbuna_240918.txt 1> pca_logs/log_mbuna_240918.txt

time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r All Whole Exploratory -e Custom 2> pca_logs/error_hybrid_analysis_240918.txt 1> pca_logs/log_hybrid_analysis_240918.txt
hybrid_analysis

time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r All Whole Exploratory -e Custom 2> pca_logs/error_cvanalysis_240918.txt 1> pca_logs/log_cvanalysis_240918.txt
CVAnalysis

time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r All Whole Exploratory -e Custom 2> pca_logs/error_yhpedigree_240918.txt 1> pca_logs/log_yhpedigree_240918.txt
YHPedigree

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# to get the exploratory regions to work 240918
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r Exploratory -e Lake_Malawi
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r Exploratory -e Rock_Sand
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r Exploratory -e Sand
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r Exploratory -e Mbuna

time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r Exploratory -e Custom
hybrid_analysis

time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r Exploratory -e Custom
CVAnalysis

time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r Exploratory -e Custom
YHPedigree



time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r All Whole Exploratory -e Rhamphochromis Diplotaxodon 2> pca_logs/error_pelagic_240918.txt 1> pca_logs/log_pelagic_240918.txt

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Rerunning analyses with edits from Patrick

# Not working b/c Nyererei are in this analysis and are absent in the 498 cohort vcf file
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r All Whole Exploratory -e Custom 2> pca_logs/error_phylogenyfigure_240918.txt 1> pca_logs/log_phylogenyfigure_240918.txt
PhylogenyFigure

# Running
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r All Whole Exploratory -e Custom 2> pca_logs/error_pcafigure_240918.txt 1> pca_logs/log_pcafigure_240918.txt
PCAFigure

# Running with having exluded the mom MC-5G11G-f and adding in kocherMC-female and MC-010-f as proxies. 
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r All Whole Exploratory -e Custom 2> pca_logs/error_yhpedigree_240918.txt 1> pca_logs/log_yhpedigree_240918.txt
YHPedigree

# Running 
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset --plink -r All Whole Exploratory -e Custom 2> pca_logs/error_cvanalysis_240918.txt 1> pca_logs/log_cvanalysis_240918.txt
CVAnalysis

# reran with the updated figure parameters
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -r All Whole Exploratory -e Mbuna

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Regeneerating figures with updated colors and sizes:
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -r All Whole Exploratory -e Lake_Malawi 
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -r All Whole Exploratory -e Rock_Sand 
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -r All Whole Exploratory -e Sand
time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -r All Whole Exploratory -e Mbuna

time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -r All Whole Exploratory -e Custom
hybrid_analysis

time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -r All Whole Exploratory -e Custom
CVAnalysis

time python pca_maker.py Mzebra_GT3 /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -r All Whole Exploratory -e Custom
YHPedigree

"""


