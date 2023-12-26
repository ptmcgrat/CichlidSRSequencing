import argparse, pdb, os, subprocess, pathlib
import pandas as pd
from helper_modules.nikesh_file_manager import FileManager as FM
from cyvcf2 import VCF
import plotly.express as px
import plotly.graph_objs as go

parser = argparse.ArgumentParser(usage = "This pipeline is for running pca analysis on a filtered vcf file. Note that the script will assume tha name of the vcf file to analyze is pass_variants_master_file.vcf.gz")
parser.add_argument('genome', help = 'name of reference genome used in the creation of the VCF files')
parser.add_argument('output_dir', help = 'absolute filepath to an output directory')
parser.add_argument('--sample_subset', help = 'This flag will restruct the samples used for generating the PCA plot to a max of three three samples for any given organism.', action= 'store_true')
parser.add_argument('-e', '--ecogroups', help = 'one or multiple eco group names for filtering the data', choices = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'Riverine', 'AC', 'Non_Riverine', 'All', 'Lake_Malawi', 'Rock_Sand', 'Sand'], nargs = '*', default = ['All'])
parser.add_argument('-r', '--regions', help = 'list of linkage groups for which analyses will run', nargs = '*', default = ['All'])
parser.add_argument('-l', '--local_test', help = 'call this flag to predefine variables for testing on local machine', action='store_true')
parser.add_argument('-p', '--plink', help = 'use this flag to generate new eigenvalue/vector files using plink', action='store_true')
# parser.add_argument('input_vcffile', help = 'absolute filepath to the filtered, gzipped input file')
# parser.add_argument('sample_database', help = 'sample database that lists ecotype for each sample')
args = parser.parse_args()

"""
To Do:
- The location of the intervals for LG 11 are also stored in the 'lg11_intervals' directory which the script will assume is stored in the dir from which pca_maker.py is called. If this is untrue there will be issues... 
- The script needs to be able to take in a "Whole" option for regions and generate a PCA for the whole genome's variants taken together, instead of being split by EG.
- The script also needs to be able to taek in multiple regions per LG and run analyses sequentially for all regiosn provided. The regions taht will often be called are "All," "Inversion," & "Whole"
- All and Inversion regions take the whole file, then split the genome or LG11 into subfiles in parallel, then runs plink in parallel for each and outputs eigenvec files which are used for PCA creation
    - For "Whole," the splitting step needs to be skipped, and the whole file must be used as input into PLINK. 
    - The splitting methods are _split_VCF_to_LG for "all" or lg specific analyses, and "_create_inversion_files" for "Inversion"
    - These must be skipped... The above methods take in the self.samples_filtered_master_file as inpout and output to the PCA dir. Instead of outputting a split file to a "Whole" dir in .../PCA I can nmake the "Whole" dir, then pipe outputs from plink there. 

- Add code to check if the number of columns for each vcf file in the PCA dir is the same as the number of sampels in samples_to_keep.csv. If not, then make a new file. If so, keep the original
    - this will ensure overwrites of the files in the event of a --sample_subset call

- Add code to make the Multiome projectID samples larger on the PCA

NOTE:
- update the vis code to account for new projectIDs and resizing stuff
- add code to not rewrite preprocessed files if they exist. Maybe create a list of all files that get created and verify they exist in the dir. If they do and the number of samples 

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
        self.in_vcf = self.fm_obj.localOutputDir + 'vcf_concat_output/pass_variants_master_file.vcf.gz' # The in_vcf attriubute is equal to the input file name for the Object.
        if args.local_test:
            self.in_vcf = self.fm_obj.localOutputDir + 'vcf_concat_output/612_cohort_3_lg_subset.vcf.gz' # This file is a subset of the 612 cohort pass_variants_master_file. By default, the script will use this whole file as input for local testing 
        self.ecogroups = ecogroups # This attribute is the list of ecogroups used for filtering samples

        self.vcf_obj = VCF(self.in_vcf) # The VCF object is made using the CVF class from cyvcf2. The VCF class takes an input vcf. For the object, this input vcf file is the "input_vcfcfile" which is defined under self.in_vcf
        self.linkage_groups = linkage_groups # the object will now take in linkage groups using args.regions. If default, it will default to the first 22 lgs
        # Create a filepath that includes name of input file and ecogroups on whcih the analysis is performed. The below code is necessary to make sure the storage of PCAs are within a filepath unique to each analysis and unique for any VCF file passed in
        file_version = self.in_vcf.split('/')[-1].split('.')[0] # if we run pca_maker on multiple vcf files, this allows us to ensure all plots from a particular file go to its own dir within the self.our_dir defined below
        analysis_ecogroups = '_'.join(str(eg).replace('_', '') for eg in self.ecogroups)
        self.out_dir = output_directory 
        if self.out_dir.endswith('/'): # ensure that the file path does not end with a forwardslash(/)
            self.out_dir = self.out_dir[:-1]
        pathlib.Path(self.out_dir + '/' + file_version + '/' + analysis_ecogroups).mkdir(parents=True, exist_ok=True) # This generates the outdir if it doesn't exist so later tools don't run into errors making files.
        self.out_dir = self.out_dir + '/' + file_version + '/' + analysis_ecogroups
        regions_list = [] # regions list will add various linkage group names, or add the names of regions of interest like "Whole" or the various special regions of interest if passed "Exploratory"
        for region in self.linkage_groups:
            if region in self.linkage_group_map.keys():
                regions_list.append(self.linkage_group_map[region])
            elif region == "All":
                regions_list.extend(self.vcf_obj.seqnames[0:22])
            elif region == "Exploratory":
                regions_list.extend(['lg2_YH_Inversion', 'lg4_YH_CV_Inversion', 'lg5_OB', 'lg9_RockSand_Inversion', 'lg10_MC_Insertion', 'lg10_YH_Inversion', 'lg11_Inversion', 'lg13_YH_Inversion', 'lg20_RockSand_Inversion'])
            elif region == "Whole":
                regions_list.append('Whole')
            else:
                raise Exception(region + ' is not a valid option')

        self.linkage_groups = regions_list # regions list is assigned to self.linkage_groups which will be used throughout the pipeline
        duplicate_set_test = set(self.linkage_groups) # test if any linkage groups are passed twice and throws error if True:
        if len(duplicate_set_test) != len(self.linkage_groups):
            raise Exception('A repeat region has been provided')
        
        # reset self.linkage_groups to include only the 3 lgs tested locally  if the --local_test flag is called
        if args.local_test:
            self.linkage_groups = ['NC_036787.1', 'NC_036788.1', 'NC_036798.1']
            if 'Whole' in args.regions:
                self.linkage_groups.extend(['Whole'])
            if 'Exploratory' in args.regions:
                self.linkage_groups.extend(['lg2_YH_Inversion', 'lg4_YH_CV_Inversion', 'lg5_OB', 'lg9_RockSand_Inversion', 'lg10_MC_Insertion', 'lg10_YH_Inversion', 'lg11_Inversion', 'lg13_YH_Inversion', 'lg20_RockSand_Inversion'])

        # Ensure index file exists
        assert os.path.exists(self.in_vcf + '.tbi') # uses os.path.exists to see if the input file + 'tbi' extension exists. The object will be made using args.input_vcffile and args.input_vcffile will be passed to the script as an absolute file path so the path to the dir is taken care of

    def _create_sample_filter_files(self):
        self.all_ecogroup_samples_csv = self.out_dir + '/all_samples_in_this_ecogroup.csv'
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
        self.fm_obj.downloadData(self.fm_obj.localSampleFile) # download fresh SampleDatabase.csv so we can read it in and get ecogroup information
        self.df = pd.read_csv(self.fm_obj.localSampleFile) # generate df from SampleDatabase.csv
        self.df = self.df[self.df['Platform'].isin(['ILLUMINA'])].drop_duplicates(subset='SampleID') # get rid of PacBio Samples and drops duplicates in the SampleID column leaving 612 (or eventually more) samples that we can filter below

        # this code block is used to generate a sample file of ALL samples that are present in a given ecogroup that's being run in the analysis
        ecogroup_df = self.df[self.df.Ecogroup.isin(self.ecogroups)]
        ecogroup_df.to_csv(self.all_ecogroup_samples_csv, columns = ['SampleID'], header=False, index=False)
        ecogroup_df.to_csv(self.all_ecogroup_samples_metadata, columns = ['SampleID', 'Ecogroup', 'Organism'], sep='\t', index=False)

        # This code block is to create the subset sample file and an extra metadata file that may be useful for troubleshooting. This subset file is used by plink to generate the subset PC space that all the samples will be projected on to
        if args.sample_subset:
            subset_df = pd.DataFrame()
            ecogroup_df = self.df[self.df.Ecogroup.isin(self.ecogroups)]
            for organism in ecogroup_df['Organism'].unique():
                organism_rows = ecogroup_df[ecogroup_df['Organism'] == organism]
                if organism_rows.shape[0] > 3: # organism_rows.shape[0] gives the number of row. If greater than 3 do the following. If exactly 3 or less, execute the else code
                    sampled_rows = organism_rows.sample(n=3, random_state=42) # randomly sample 3 samples from the species.
                    subset_df = pd.concat([subset_df, sampled_rows], ignore_index=True) # add them to df_filtered
                else: # if only 3 or less samples exist, the sampled_rows are the organism_rows so just concat the organism rows into subset_df
                    subset_df = pd.concat([subset_df, organism_rows], ignore_index=True) # add them to df_filtered
            subset_df.to_csv(self.subset_samples_csv, columns = ['SampleID'], header=False, index=False) # store the subset samples into a csv file in the directory for that ecogroup's analysis
            subset_df.to_csv(self.subset_samples_metadata, columns = ['SampleID', 'Ecogroup', 'Organism'], sep='\t', index=False) # write an extra csv file containing metadata information for samples included in the subsetting

        else: # we will not generate a subset vcf file that will be used in the analyses. Thus, self.subset_samples_csv & self.subset_samples_metadata will not be generated. 
            # To make it easier to keep the code revolving around the subset code that already exists, I can simply generate the same files but just make them equal to what exists in all_samples files. 
            warning = input("WARNING: You have called on the script to run without subsetting samples.\nYour PCA results may be driven by overrepresented samples in the ecogroup you're analyzing. Are you sure you wish to continue? (Y/N): ")
            if warning.upper() == 'Y':
                print("Continuing... Note that files with 'subset' in the name will not actually be a subset of your samples. Filenames are preserved out of laziness...")
            else:
                exit()
            ecogroup_df.to_csv(self.subset_samples_csv, columns = ['SampleID'], header=False, index=False) # write the same samples that created the self.all_ecogroup_samples_csv and self.all_ecogroup_samples_metadata into the subset files so that I don't have to change code in other functions
            ecogroup_df.to_csv(self.subset_samples_metadata, columns = ['SampleID', 'Ecogroup', 'Organism'], sep='\t', index=False)

    def _create_ecogroup_specific_vcf(self):
        # path to the vcf file containing samples only in the ecogroups needed for the analysis
        self.ecogroup_specific_master_vcf = self.out_dir + '/all_ecogroup_samples.vcf.gz' # for each analysis, depending on the ecogroups chosen, we generate a vcf file with all samples of only that eco group. The path and name of that file is defined here. Old file name = samples_filtered_master.vcf.gz
        self.subset_master_vcf = self.out_dir + '/subset_samples.vcf.gz' # for each analysis, thi subset vcf file will need to be generated to pull region specific variants from in the next step

        if pathlib.Path(self.ecogroup_specific_master_vcf).exists():
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

    def _preprocess_vcfs(self): # this function prunes all variants from self.ecogroup_specific_master_vcf that are not present in self.subset_master_vcf following recalculation and filtering of variants with AF < 0.01.
        # File check code to ensure these files don't already exist. If they do, the samples between existing files and the cohort's samples must match:
        # rm -r PCA interactive_PCA_outputs
        # rm af_filtered_subset_samples.vcf.gz af_filtered_subset_samples.vcf.gz.tbi af_recalculated_subset_samples.vcf.gz variants_filtered_ecogroup_samples.recode.vcf.gz variants_filtered_ecogroup_samples.recode.vcf.gz.tbi variants_to_keep.txt
        if pathlib.Path(self.out_dir + '/af_recalculated_subset_samples.vcf.gz').exists():
            # if the ecogroup_specific_master_vcf (contains all variants per sample for the ecogroups specified) exists, then this checks that the samples match exactly. If not, a new file is built by filtering for samples in the self.good_samples_csv file.
            if subprocess.run(f"bcftools query -l {self.out_dir + '/af_recalculated_subset_samples.vcf.gz'}", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout == subprocess.run(f"cat {self.out_dir + '/subset_samples.csv'}", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout: # checks if the output from printing the sample names from samples_to_keep.csv and the column names from samples_filtered_master.vcf.gz are the sample
                print(f'\nThe file af_recalculated_subset_samples.vcf.gz exists in the ecogroup dir and samples within the file match those in subset_samples.csv. VCF preprocessing will be skipped\n')
            else: # if the samples differ between the existing af_recalculated_subset_samples.vcf.gz and the subset_samples.csv, redo all preprocessing steps. 
                print('\nSamples between the existing af_recalculated_subset_samples.vcf.gz and subset_samples.csv are different. Rerunning VCF preprocessing...\n')
                print('RESTARTING VCF FILE PREPROCESSING...\n')
                print('RECALCULATING AF IN THE SUBSET FILE\n')
                subprocess.run(['bcftools', '+fill-tags', self.subset_master_vcf, '-Oz', '-o', self.out_dir + '/af_recalculated_subset_samples.vcf.gz', '--', '-t', 'AF'])
                print('FILTERING OUT VARIANTS WITH AF < 0.01\n')
                subprocess.run(['bcftools', 'view', '--exclude', 'AF<0.01 || AF == 1', self.out_dir + '/af_recalculated_subset_samples.vcf.gz', '-o', self.out_dir + '/af_filtered_subset_samples.vcf.gz', '-Oz'])
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
            print('FILTERING OUT VARIANTS WITH AF < 0.01\n')
            subprocess.run(['bcftools', 'view', '--exclude', 'AF<0.01 || AF == 1', self.out_dir + '/af_recalculated_subset_samples.vcf.gz', '-o', self.out_dir + '/af_filtered_subset_samples.vcf.gz', '-Oz'])
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
            elif lg == "lg2_YH_Inversion":
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
                        p2 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.out_dir + '/af_filtered_subset_samples.vcf.gz', '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset.vcf.gz', '-O', 'z']) # takes in the subset_samples file with the remaining variants with AF > 0.01 ecogroup_specific_master_vcf and filters out each LG and writes the file into appropriate PCA dir.
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
                    p2 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.out_dir + '/af_filtered_subset_samples.vcf.gz', '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset.vcf.gz', '-O', 'z']) # takes in the subset_samples file with the remaining variants with AF > 0.01 ecogroup_specific_master_vcf and filters out each LG and writes the file into appropriate PCA dir.
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
        exploratory_regions_list = ['lg2_YH_Inversion', 'lg4_YH_CV_Inversion', 'lg5_OB', 'lg9_RockSand_Inversion', 'lg10_MC_Insertion', 'lg10_YH_Inversion', 'lg11_Inversion', 'lg13_YH_Inversion', 'lg20_RockSand_Inversion']
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
                subprocess.run(['plink2', '--vcf', self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf.gz', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole', '--allow-extra-chr']) # whole vcf file pfile generation. Use the variants_filtered_ecogroup_samples.recode.vcf.gz file as the "whole" file.
                subprocess.run(['plink2', '--vcf', self.out_dir + '/af_filtered_subset_samples.vcf.gz', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset', '--allow-extra-chr']) # subset sample vcf file pfile generation. The subset file is the af_filtered_subset_samples.vcf.gz in self.out_dir
                subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole', '--set-missing-var-ids', '@:#', '--make-pgen', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected', '--allow-extra-chr'])
                subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset', '--freq', 'counts', '--pca', 'allele-wts', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--max-alleles', '2'])

                # modify the .acounts file to eliminate 0 count alleles:
                acounts_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.acount', sep='\t')
                no_extremes_df = acounts_df[(acounts_df['ALT_CTS'] != 0) & (acounts_df['ALT_CTS'] != max(acounts_df['OBS_CT']))]
                no_extremes_df.to_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.acount', sep='\t', index=False)

                # genenrate the .sscore file with eigenvectors for each samples after projection on to the PC space from the subset PCA
                subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected', '--read-freq', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.acount', '--score', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.eigenvec.allele', '2', '5', 'header-read', 'no-mean-imputation', 'variance-standardize', '--score-col-nums', '6-15', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection', '--allow-extra-chr'])


            else: #lg in self.linkage_group_map.values():
                pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True)
                if not pathlib.Path(self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_ecogroup.vcf.gz').exists(): # NOTE: error checking... revisit later 
                    print('ERROR: THE FILE ' + lg + '.VCF.GZ DOES NOT EXIST. MUST RUN _SPLIT_VCF_TO_LG TO CREATE IT...')
                    raise Exception
                else: # do all the plink pca magic here assuming you're going one linkage group at a time... figure out parallelization as I code this referencing the previous function. Parallelization may not be possible because plink likes to execute immediately and doesn't listen to the Popen constructor. It all executres before Popen.communicate() is called.
                    pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True)

                    subprocess.run(['plink2', '--vcf', self.out_dir + '/variants_filtered_ecogroup_samples.recode.vcf.gz', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole', '--allow-extra-chr']) # whole vcf file pfile generation for a given lg
                    subprocess.run(['plink2', '--vcf', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset.vcf.gz', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset', '--allow-extra-chr']) # subset sample vcf file pfile generation for a given lg
                    subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole', '--set-missing-var-ids', '@:#', '--make-pgen', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected', '--allow-extra-chr'])
                    subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset', '--freq', 'counts', '--pca', 'allele-wts', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--max-alleles', '2'])

                    # modify the .acounts file to eliminate 0 count alleles:
                    acounts_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.acount', sep='\t')
                    no_extremes_df = acounts_df[(acounts_df['ALT_CTS'] != 0) & (acounts_df['ALT_CTS'] != max(acounts_df['OBS_CT']))]
                    no_extremes_df.to_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.acount', sep='\t', index=False)

                    # genenrate the .sscore file with eigenvectors for each samples after projection on to the PC space from the subset PCA
                    subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected', '--read-freq', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.acount', '--score', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.eigenvec.allele', '2', '5', 'header-read', 'no-mean-imputation', 'variance-standardize', '--score-col-nums', '6-15', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection', '--allow-extra-chr'])


    def _create_interactive_pca(self, linkage_group_list): # uses plotly to generate interactive PCA html outputs
        self.plotly_out = self.out_dir + '/interactive_PCA_outputs/' # define outdir
        pathlib.Path(self.plotly_out).mkdir(parents=True, exist_ok=True) # build the file path with pathlib.Path
        color_map = {'Mbuna': 'purple', 'AC': 'limegreen', 'Shallow_Benthic': 'red', 'Deep_Benthic': 'blue', 'Rhamphochromis': 'brown', 'Diplotaxodon': 'orange', 'Utaka': 'darkgreen', 'Riverine': 'pink'}
        shape_map = {'MalinskyData': 'square', 'Streelman_McGrathData': 'diamond', 'BrainDiversity_s1': 'star', 'MC_males': 'circle', 'MC_females': 'circle-open'}

        for lg in linkage_group_list:
            eigen_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection.sscore', sep='\t')
            eigen_df = eigen_df.rename(columns = {'#IID':'SampleID'})
            df_merged = pd.merge(eigen_df, self.df, on=['SampleID'])
            if lg.startswith('NC'):
                plot_title = list(self.linkage_group_map.keys())[list(self.linkage_group_map.values()).index(lg)]
            else:
                plot_title = lg
            fig = px.scatter(df_merged, x='PC1_AVG', y='PC2_AVG', color='Ecogroup', symbol='ProjectID_2', color_discrete_map=color_map, symbol_map=shape_map, title=plot_title, hover_data=['SampleID', 'Ecogroup', 'Organism', 'ProjectID_2'])
            # else: # NOTE: if we want to generate a subset PCA uncomment and incorporate below code into the function
            #     fig = px.scatter(df_merged, x='PC1', y='PC2', color='Ecogroup', symbol='ProjectID', color_discrete_map=color_map, symbol_map=shape_map, title=lg, hover_data=['SampleID', 'Ecogroup', 'Organism', 'ProjectID'])
            larger_size = 9
            smaller_size = 3
            fig.update_traces(marker=dict(size=smaller_size), selector=dict(marker_symbol='square'))
            fig.update_traces(marker=dict(size=smaller_size), selector=dict(marker_symbol='diamond'))
            fig.update_traces(marker=dict(size=larger_size), selector=dict(marker_symbol='circle'))
            fig.update_traces(marker=dict(size=larger_size), selector=dict(marker_symbol='circle-open'))
            fig.update_traces(marker=dict(size=larger_size), selector=dict(marker_symbol='star'))
            fig.write_html(self.plotly_out + lg + '_plotlyPCA.html')

    def create_PCA(self):
        # Order of hidden methods to perform the analysis
        self._create_sample_filter_files()
        self._create_ecogroup_specific_vcf()
        self._preprocess_vcfs()
        self._split_VCF_to_LG(self.linkage_groups)
        if args.plink:
            try:
                self._create_eigenfiles_per_LG(self.linkage_groups) # This line is used to test the _create_PCA_linakge hidden method using only LG1.
            except:
                pass
        self._create_interactive_pca(self.linkage_groups)

if __name__ == "__main__":
    pca_obj = PCA_Maker(args.genome, args.ecogroups, args.regions, args.output_dir)
    pca_obj.create_PCA()
    print('PIPELINE RUN SUCCESSFUL')

"""
For running on Utaka:
time python pca_maker.py Mzebra_UMD2a /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -p -r All Whole Exploratory -e All 2> pca_logs/error_all_.txt 1> pca_logs/log_all.txt
time python pca_maker.py Mzebra_UMD2a /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -p -r All Whole Exploratory -e Lake_Malawi 2> pca_logs/error_lm_.txt 1> pca_logs/log_lm.txt
time python pca_maker.py Mzebra_UMD2a /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -p -r All Whole Exploratory -e Non_Riverine 2> pca_logs/error_non_riverine_.txt 1> pca_logs/log_non_riverine.txt
time python pca_maker.py Mzebra_UMD2a /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -p -r All Whole Exploratory -e Rock_Sand 2> pca_logs/error_rock_sand_.txt 1> pca_logs/log_rock_sand.txt
time python pca_maker.py Mzebra_UMD2a /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs --sample_subset -p -r All Whole Exploratory -e Sand 2> pca_logs/error_sand_.txt 1> pca_logs/log_sand.txt

"""
