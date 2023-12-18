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

NOTE: The sample filtering and randomization seems to be working, but I count that there should be 613 unique samples excluding outliers and Shallow/Deep Benthic sample and the 2 other samples that cannot run through Bamfile generation... Yet I get 610. Need to address this


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
            self.in_vcf = self.fm_obj.localOutputDir + 'vcf_concat_output/10k_variants_612_cohort_3_lg_subset.vcf.gz' # This file is a subset of the 612 cohort pass_variants_master_file. By default, the script will use this whole file as input for local testing 
        self.ecogroups = ecogroups # This attribute is the list of ecgogroups used for filtering samples

        self.vcf_obj = VCF(self.in_vcf) # The VCF object is made using the CVF class from cyvcf2. The VCF class takes an input vcf. For the object, this input vcf file is the "input_vcfcfile" which is defined under self.in_vcf
        self.linkage_groups = linkage_groups # the object will now take in linkage groups using args.regions. If default, it will default to the first 22 lgs

        # Create a filepath that includes name of input file and ecogroups on whcih the analysis is performed. The below code is necessary to make sure the storage of PCAs are within a filepath unique to each analysis and unique for any VCF file passed in
        file_version = self.in_vcf.split('/')[-1].split('.')[0] # if we run pca_maker on multiple vcf files, this allows us to ensure all plots from a particular file go to its own dir within the self.our_dir defined below
        analysis_ecogroups = '_'.join(str(eg).replace('_', '') for eg in self.ecogroups)
        self.out_dir = output_directory
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
            if 'Whole' in args.regions:
                self.linkage_groups == ['Whole']
            else:
                self.linkage_groups = ['NC_036787.1', 'NC_036788.1', 'NC_036798.1']
            

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
        # else: # NOTE: not sure what this is supposed to take care of.... if sample_subset is not called... then don't we not generate a subset and we just do a normal PCA? Commenting it all out for now 
        #     subset_df = self.df[self.df.Ecogroup.isin(self.ecogroups)]
        #     subset_df.drop_duplicates(subset='SampleID').to_csv(self.subset_samples_csv, columns = ['SampleID'], header = False, index = False) # store all unique subset samples into a csv file in the directory for that ecogroup's analysis
        #     subset_df.drop_duplicates(subset='SampleID').to_csv(self.subset_samples_metadata, columns = ['SampleID', 'Ecogroup', 'Organism'], sep = '\t', index=False) # write an extra csv file containing metadata information for the unique subset samples included in the above file



        """
        OLD CODE FOR SUBSETTING SAMPLES:
        # self.s_dt[self.s_dt.Ecogroup.isin(self.ecogroups)] is returning all rows where the self.ecogroups values are contained in the "Ecogroup" column of the excel sheet
        # The .SampleID.unique() is filtering these rows to include only unique values in the SampleIDs column. However, this creates a numpy array so the pd.DataFrame wrapper around the whole thing converts this numpy array to a pandas dataframe. 
        # The to_csv(self.good_samples_txt, header = False, index = False) converts into a csv file that bcftools will be able to take as input. The index/header = False eliminate any index/header information, leaving only filtered samples. 
        # Note that the name of the output file is self.good_samples.csv which is where the file pathing for the output file is taken care of
        self.df = pd.read_csv(self.sample_database) # generate df from SampleDatabase.csv
        if args.sample_subset: # if sample_subset flag is called then self.df_filtered is initialized as a blank DataFrame. 
            self.df_filtered = pd.DataFrame()
            ecogroup_df = self.df[self.df.Ecogroup.isin(self.ecogroups)] # gets samples matching eco groups provided to the pipeline
            for organism in ecogroup_df['Organism'].unique(): # for each unique organism in the ecogroups provided... 
                organism_rows = ecogroup_df[ecogroup_df['Organism'] == organism] # gets all rows that match that organism name
                if organism_rows.shape[0] >= 3: # if the number of rows (shape[0]) > 3, then randomly sample using seed 100 and add to self.df_filtered using pd.concat (note that df.append is no longer supported in pandas as of May 2023)
                    sampled_rows = organism_rows.sample(n=3, random_state=100)
                    self.df_filtered = pd.concat([self.df_filtered, sampled_rows], ignore_index=True)
                else:
                    self.df_filtered = pd.concat([self.df_filtered, organism_rows], ignore_index=True) # if total rows for organism are less than 3, add all rows to self.df_filtered
            pd.DataFrame(self.df_filtered[self.df_filtered.Platform != 'PACBIO'].SampleID.unique()).to_csv(self.good_samples_csv, header = False, index = False) # get rid of Pacbio samples and leave only Illumina ones
            # pd.DataFrame(self.df_filtered.SampleID.unique()).to_csv(self.good_samples_csv, header = False, index = False) # using the SampleID column, get all unique IDs and write them to self.good_samples_csv which is used to pull samples for VCF filtering in the pipeline 
        elif args.local_test: # if using a local file, just take in all of the samples in the small file and set those all to be the samples for testing
            self.df_filtered = self.df[self.df.Ecogroup.isin(self.ecogroups)]
            local_samples = subprocess.check_output(['bcftools', 'query', '-l', self.in_vcf], encoding='utf-8').split('\n')
            del local_samples[-1]
            with open(self.good_samples_csv, 'w') as fh:
                for sample in local_samples:
                    fh.write(sample + '\n')
        else:
            self.df_filtered = self.df[self.df.Ecogroup.isin(self.ecogroups)] # if sample_subset flag is not called, just get all samples for the given ecogroups and proceed
            pd.DataFrame(self.df_filtered[self.df_filtered.Platform != 'PACBIO'].SampleID.unique()).to_csv(self.good_samples_csv, header = False, index = False) # get rid of Pacbio samples and leave only Illumina ones

        """

    def _create_ecogroup_specific_vcf(self):
        # path to the vcf file containing samples only in the ecogroups needed for the analysis
        self.ecogroup_specific_master_vcf = self.out_dir + '/all_ecogroup_samples.vcf.gz' # for each analysis, depending on the ecogroups chosen, we generate a vcf file with all samples of only that eco group. The path and name of that file is defined here. Old file name = samples_filtered_master.vcf.gz
        self.subset_master_vcf = self.out_dir + '/subset_samples.vcf.gz' # for each analysis, thi subset vcf file will need to be generated to pull region specific variants from in the next step

        if pathlib.Path(self.ecogroup_specific_master_vcf).exists(): # if the ecogroup_specific_master_vcf (contains all variants per sample for the ecogroups specified) exists, then this checks that the samples match exactly. If not, a new file is built by filtering for samples in the self.good_samples_csv file.
            if subprocess.run(f"bcftools query -l {self.ecogroup_specific_master_vcf}", shell=True, stdout=subprocess.DEVNULL, encoding='utf-8').stdout == subprocess.run(f"cat {self.all_ecogroup_samples_csv}", shell=True, stdout=subprocess.DEVNULL, encoding='utf-8').stdout: # checks if the output from printing the sample names from samples_to_keep.csv and the column names from samples_filtered_master.vcf.gz are the sample
                print(f'\nThe file {self.ecogroup_specific_master_vcf} exists and samples within the file match those in self.good_samples_csv. New samples_filtered_master_vcf file will not be built.')
            else:
                p1 = subprocess.Popen(['bcftools', 'view', self.in_vcf, '--samples-file', self.all_ecogroup_samples_csv, '--min-ac=1', '--no-update', '-o', self.ecogroup_specific_master_vcf, '-O', 'z'])
                p2 = subprocess.Popen(['bcftools', 'view', self.in_vcf, '--samples-file', self.subset_samples_csv, '-o', self.subset_master_vcf, '-O', 'z'])
                print('\nGENERATING A VCF FILE CONTAINING ONLY SAMPLES OF THE ECOGROUPS SPECIFIED.')
                print('\nGENERATING A VCF FILE WITH THE SUBSET SAMPLES FOR THE ECOGROUPS SPECIFIED.')
                p1.communicate()
                p2.communicate()

                p3 = subprocess.Popen(['tabix', '-p', 'vcf', self.ecogroup_specific_master_vcf])
                p4 = subprocess.Popen(['tabix', '-p', 'vcf', self.subset_master_vcf])
                print('FILTERED SAMPLES FILE GENERATED. INDEXING ALL SAMPLE VCF FILE...')
                print('INDEXING SUBSET SAMPLE VCF FILE...')
                p3.communicate()
                p4.communicate()
                print('INDEX CREATED...')
        else: # if no ecogroup_specific_master_vcf file exists yet, go ahead and use the samples to keep file to generate a cvf file containing only samples of the ecogroup being analyzed.
            p1 = subprocess.Popen(['bcftools', 'view', self.in_vcf, '--samples-file', self.all_ecogroup_samples_csv, '--min-ac=1', '--no-update', '-o', self.ecogroup_specific_master_vcf, '-O', 'z'])
            p2 = subprocess.Popen(['bcftools', 'view', self.in_vcf, '--samples-file', self.subset_samples_csv, '-o', self.subset_master_vcf, '-O', 'z'])
            print('\nGENERATING A VCF FILE CONTAINING ONLY SAMPLES OF THE ECOGROUPS SPECIFIED.')
            print('\nGENERATING A VCF FILE WITH THE SUBSET SAMPLES FOR THE ECOGROUPS SPECIFIED.')
            p1.communicate()
            p2.communicate()

            p3 = subprocess.Popen(['tabix', '-p', 'vcf', self.ecogroup_specific_master_vcf])
            p4 = subprocess.Popen(['tabix', '-p', 'vcf', self.subset_master_vcf])
            print('FILTERED SAMPLES FILE GENERATED. INDEXING ALL SAMPLE VCF FILE...')
            print('INDEXING SUBSET SAMPLE VCF FILE...')
            p3.communicate()
            p4.communicate()
            print('INDEX CREATED...')

        """
        print('\nGENERATING A VCF FILE CONTAINING ONLY SAMPLES OF THE ECOGROUPS SPECIFIED.')
        subprocess.run(['bcftools', 'view', self.in_vcf, '--samples-file', self.all_ecogroup_samples_csv, '-o', self.ecogroup_specific_master_vcf, '-O', 'z']) # code to generate a master_vcf file of filtered samples
        print('FILTERED SAMPLES FILE GENERATED. INDEXING FILE...')
        subprocess.run(['tabix', '-p', 'vcf', self.ecogroup_specific_master_vcf]) # code to generate an index for this file using bcftools at the location of plink_master_vcf
        print('INDEX CREATED...')
        """

        """
        # Below block first checks if any existing samples_filtered_master_vcf file is present. If so, check if it contains the same sampleIDs as the samples we want to include in the analysis we're currently running  If it does don't regenerate it. If it doesn't, make a new samples_filtered_master_vcf file 
        if pathlib.Path(self.samples_filtered_master_vcf).exists(): # if the samples_filtered_master_vcf (contains all variants per sample for the ecogroups specified) exists, then this checks that the samples match exactly. If not, a new file is built by filtering for samples in the self.good_samples_csv file.
            if subprocess.run(f"bcftools query -l {self.samples_filtered_master_vcf}", shell=True, stdout=subprocess.DEVNULL, encoding='utf-8').stdout == subprocess.run(f"cat {self.good_samples_csv}", shell=True, stdout=subprocess.DEVNULL, encoding='utf-8').stdout: # checks if the output from printing the sample names from samples_to_keep.csv and the column names from samples_filtered_master.vcf.gz are the sample
                print(f'\nThe file {self.samples_filtered_master_vcf} exists and samples within the file match those in self.good_samples_csv. New samples_filtered_master_vcf file will not be built.')
                pass
            else:
                print('\nGENERATING A VCF FILE CONTAINING ONLY SAMPLES OF THE ECOGROUPS SPECIFIED.')
                subprocess.run(['bcftools', 'view', self.in_vcf, '--samples-file', self.good_samples_csv, '-o', self.samples_filtered_master_vcf, '-O', 'z']) # code to generate a master_vcf file of filtered samples
                print('FILTERED SAMPLES FILE GENERATED. INDEXING FILE...')
                subprocess.run(['tabix', '-p', 'vcf', self.samples_filtered_master_vcf]) # code to generate an index for this file using bcftools at the location of plink_master_vcf
                print('INDEX CREATED...')
        else: # if no samples_filtered_master_vcf file exists, go ahead and use the samples to keep file to generate a cvf file containing only samples of the ecogroup being analyzed. 
            print('\nGENERATING A VCF FILE CONTAINING ONLY SAMPLES OF THE ECOGROUPS SPECIFIED.')
            subprocess.run(['bcftools', 'view', self.in_vcf, '--samples-file', self.good_samples_csv, '-o', self.samples_filtered_master_vcf, '-O', 'z']) # code to generate a master_vcf file of filtered samples
            print('FILTERED SAMPLES FILE GENERATED. INDEXING FILE...')
            subprocess.run(['tabix', '-p', 'vcf', self.samples_filtered_master_vcf]) # code to generate an index for this file using bcftools at the location of plink_master_vcf
            print('INDEX CREATED...')
        """

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
                # TODO: come back to the below if block
                if pathlib.Path(self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz').exists(): # For each linkage groups' dir in the PCA dir, if the LG's vcf file exists, then skip the generation of that file from the samples_filtered_master.vcf.gz file.
                    if subprocess.run(f"bcftools query -l {self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz'}", shell=True, stdout=subprocess.DEVNULL, encoding='utf-8').stdout == subprocess.run(f"cat {self.all_ecogroup_samples_csv}", shell=True, stdout=subprocess.DEVNULL, encoding='utf-8').stdout: # checks if the subset vcf file that already exists has matching sample names as those in the samples_to_keep.csv file.
                        print('The file ' + lg + '.vcf.gz exists and has the same samples as self.good_samples_csv. A file for ' + lg + ' will not be generated.')
                    else:
                        print('Generating a VCF file containing the subset samples for ' + lg + '...')
                        p1 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.ecogroup_specific_master_vcf, '-o', self.out_dir + '/PCA' + lg + '/' + lg + '_whole_ecogroup.vcf.gz', '-O', 'z']) # takes in ecogroup_specific_master_vcf and filters out each LG and writes the file into appropriate PCA dir.
                        p2 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.subset_master_vcf, '-o', self.out_dir + '/PCA' + lg + '/' + lg + '_sample_subset.vcf.gz', '-O', 'z']) # takes in ecogroup_specific_master_vcf and filters out each LG and writes the file into appropriate PCA dir.])
                        # p1 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.samples_filtered_master_vcf, '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '-O', 'z']) 
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
                    p1 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.ecogroup_specific_master_vcf, '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_ecogroup.vcf.gz', '-O', 'z']) # takes in ecogroup_specific_master_vcf and filters out each LG and writes the file into appropriate PCA dir.
                    p2 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.subset_master_vcf, '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset.vcf.gz', '-O', 'z']) # takes in ecogroup_specific_master_vcf and filters out each LG and writes the file into appropriate PCA dir.])
                    processes1.append(p1)
                    processes2.append(p2)
                    if len(processes1) == len(linkage_group_list): # parallelization code
                        for proc1 in processes1:
                            proc1.communicate()
                        processes1 = []
                    if len(processes2) == len(linkage_group_list):
                        for proc2 in processes2:
                            proc2.communicate()

    """
    def _split_VCF_to_LG(self, linkage_group_list):
        processes = []
        for lg in linkage_group_list:
            if lg == 'Whole':
                continue
            elif lg == "lg2_YH_Inversion":
                self._create_exploratory_region_eigen_files()
            elif lg in self.linkage_group_map.values():
                pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True) # generate the file paths to he split LG dirs within a dir named "PCA"
                if pathlib.Path(self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz').exists(): # For each linkage groups' dir in the PCA dir, if the LG's vcf file exists, then skip the generation of that file from the samples_filtered_master.vcf.gz file.
                    if subprocess.run(f"bcftools query -l {self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz'}", shell=True, stdout=subprocess.DEVNULL, encoding='utf-8').stdout == subprocess.run(f"cat {self.all_ecogroup_samples_csv}", shell=True, stdout=subprocess.DEVNULL, encoding='utf-8').stdout: # checks if the subset vcf file that already exists has matching sample names as those in the samples_to_keep.csv file.
                        print('The file ' + lg + '.vcf.gz exists and has the same samples as self.good_samples_csv. A file for ' + lg + ' will not be generated.')
                    else:
                        print('Generating a subset VCF file for ' + lg + '...')
                        p1 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.samples_filtered_master_vcf, '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '-O', 'z']) # takes in samples_filtered_master.vcf and filters out each LG and writes the file into appropriate PCA dir.
                        processes.append(p1)
                        if len(processes) == len(linkage_group_list): # parallelization code
                            for proc in processes:
                                proc.communicate()
                            processes = []
                else:
                    print('Generating a subset VCF file for ' + lg + '...')
                    p1 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.samples_filtered_master_vcf, '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '-O', 'z']) # takes in samples_filtered_master.vcf and filters out each LG and writes the file into appropriate PCA dir.
                    processes.append(p1)
                    if len(processes) == len(linkage_group_list): # parallelization code
                        for proc in processes:
                            proc.communicate()
                        processes = []
    """

    def _create_exploratory_region_eigen_files(self):
        processes1 = []
        processes2 = []
        remake_index = False
        exploratory_regions_list = ['lg2_YH_Inversion', 'lg4_YH_CV_Inversion', 'lg5_OB', 'lg9_RockSand_Inversion', 'lg10_MC_Insertion', 'lg10_YH_Inversion', 'lg11_Inversion', 'lg13_YH_Inversion', 'lg20_RockSand_Inversion']
        for region in exploratory_regions_list:
            if pathlib.Path(self.out_dir + '/PCA/' + region + '/' + region + '.vcf.gz').exists():
                if subprocess.run(f"bcftools query -l {self.out_dir + '/PCA/' + region + '/' + region + '.vcf.gz'}", shell=True, stdout=subprocess.DEVNULL, encoding='utf-8').stdout == subprocess.run(f"cat {self.good_samples_csv}", shell=True, stdout=subprocess.DEVNULL, encoding='utf-8').stdout: # checks if an existing exploratory region vcf file contains the same samples as samples_to_keep.csv. If not, make a new one
                    print(f'VCF file for {region} exists and number of samples in the vcf file matches the number of samples in self.good_samples_csv. Skippping and jumping to next region')
                else:
                    remake_index = True
                    self.special_interval_file = os.getcwd() + '/special_intervals/' + region + '.interval_list' # gets the interval file for the region specified by searching in the "special_intervals" dir in the dir from where pca_maker.py was called from
                    pathlib.Path(self.out_dir + '/PCA/' + region + '/').mkdir(parents=True, exist_ok=True) # builds filepath, if needed for the region of interest
                    p1 = subprocess.Popen(['gatk', 'SelectVariants', '-V', self.samples_filtered_master_vcf, '-L', self.special_interval_file, '-O', self.out_dir + '/PCA/' + region + '/' + region + '.vcf', '--verbosity', 'WARNING']) # gatk selectvariants gets the variants for the samples in self.samples_filtered_master_vcf in the regions dictated by self.special_interval_file and outputs them to the appropriate directory in PCA
                    processes1.append(p1)
                    if len(processes1) == len(exploratory_regions_list):
                        for proc in processes1:
                            proc.communicate()
                        processes1 = []
            else: # if the vcf file for the region doesn't exist in the first place, then run the gatk SelectVariants code followed by zipping and indexing 
                self.special_interval_file = os.getcwd() + '/special_intervals/' + region + '.interval_list' # gets the interval file for the region specified by searching in the "special_intervals" dir in the dir from where pca_maker.py was called from
                pathlib.Path(self.out_dir + '/PCA/' + region + '/').mkdir(parents=True, exist_ok=True) # builds filepath, if needed for the region of interest
                p1 = subprocess.Popen(['gatk', 'SelectVariants', '-V', self.samples_filtered_master_vcf, '-L', self.special_interval_file, '-O', self.out_dir + '/PCA/' + region + '/' + region + '.vcf', '--verbosity', 'WARNING']) # gatk selectvariants gets the variants for the samples in self.samples_filtered_master_vcf in the regions dictated by self.special_interval_file and outputs them to the appropriate directory in PCA
                processes1.append(p1)
                if len(processes1) == len(exploratory_regions_list):
                    for proc in processes1:
                        proc.communicate()
                    processes1 = []

        for region in exploratory_regions_list: # note that currently, this will not generate a new index if the file already exists, which isn't good...
            if pathlib.Path(self.out_dir + '/PCA/' + region + '/' + region + '.vcf.gz.idx').exists() & remake_index:
                p2 = subprocess.Popen(['bgzip', '-f', self.out_dir + '/PCA/' + region + '/' + region + '.vcf'])
                processes2.append(p2)
                if len(processes2) == len(exploratory_regions_list):
                    for proc in processes2:
                        proc.communicate()
                    processes2 = []
            elif pathlib.Path(self.out_dir + '/PCA/' + region + '/' + region + '.vcf.gz.idx').exists():
                print(f'index file for {region} exists  skippping and jumping to next region')
            else:
                p2 = subprocess.Popen(['bgzip', '-f', self.out_dir + '/PCA/' + region + '/' + region + '.vcf'])
                processes2.append(p2)
                if len(processes2) == len(exploratory_regions_list):
                    for proc in processes2:
                        proc.communicate()
                    processes2 = []

    def _create_eigenfiles_per_LG(self, linkage_group_list): # new hidden method that will create PCA plots for each LG in sample. It will define attributes for the object and also takes in a lingage grouup. Calling on this method in a for lopp should generate the eigenvalue/vector files needed per lg in self.contigs
        """
        1. generate pfiles for whole and subset vcf files in each LG dir
        2. fix missing IDs in each linkage groups' whole vcf's pfile and rename that file as "whole_corrected"
        3. run step1 of pca
        4. remove any 0 allele count variants from the acounts file using pandas
        5. run pca step 2

        """

        for lg in linkage_group_list:
            if lg == 'Whole': # NOTE: ignore for now 
                pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True)
                subprocess.run(['plink2', '--vcf', self.ecogroup_specific_master_vcf, '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole', '--allow-extra-chr']) # whole vcf file pfile generation for a given lg
                subprocess.run(['plink2', '--vcf', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset.vcf.gz', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset', '--allow-extra-chr']) # subset sample vcf file pfile generation for a given lg
                subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole', '--set-missing-var-ids', '@:#', '--make-pgen', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected', '--allow-extra-chr'])
                subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset', '--freq', 'counts', '--pca', 'allele-wts', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--max-alleles', '2'])
                
                
                # modify the .acounts file to eliminate 0 count alleles:
                acounts_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.acount', sep='\t')
                nonzero_acounts_df = acounts_df[acounts_df['ALT_CTS'] != 0]
                nonzero_acounts_df.to_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.acount', sep='\t', index=False)

                # genenrate the .sscore file with eigenvectors for each samples after projection on to the PC space from the subset PCA
                subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected', '--read-freq', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.acount', '--score', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.eigenvec.allele', '2', '5', 'header-read', 'no-mean-imputation', 'variance-standardize', '--score-col-nums', '6-15', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection', '--allow-extra-chr'])

                # pathlib.Path(self.out_dir + '/PCA/' + 'Whole/').mkdir(parents=True, exist_ok=True)
                # print("RUNNING PLINK TO TRANSFORM THE WHOLE FILTERED VCF FILE'S DATA TO A PLINK OBJECT...")
                # subprocess.run(['plink', '--vcf', self.samples_filtered_master_vcf, '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--out', self.out_dir + '/PCA/' + 'Whole/' + 'test' ]) 
                # print('GENERATING EIGENVALUE AND EIGENVECTOR FILES...')
                # subprocess.run(['plink', '--vcf', self.samples_filtered_master_vcf, '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--extract',  self.out_dir + '/PCA/' + 'Whole/' + 'test.prune.in', '--make-bed', '--pca', '--out',  self.out_dir + '/PCA/' + 'Whole/' + 'test'])
            else: #lg in self.linkage_group_map.values():
                pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True)
                if not pathlib.Path(self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_ecogroup.vcf.gz').exists(): # NOTE: error checking... revisit later 
                    print('ERROR: THE FILE ' + lg + '.VCF.GZ DOES NOT EXIST. MUST RUN _SPLIT_VCF_TO_LG TO CREATE IT...')
                    raise Exception
                else: # do all the plink pca magic here assuming you're going one linkage group at a time... figure out parallelization as I code this referencing the previous function. Parallelization may not be possible because plink likes to execute immediately and doesn't listen to the Popen constructor. It all executres before Popen.communicate() is called.
                    pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True)
                    # pdb.set_trace()
                    subprocess.run(['plink2', '--vcf', self.ecogroup_specific_master_vcf, '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole', '--allow-extra-chr']) # whole vcf file pfile generation for a given lg
                    subprocess.run(['plink2', '--vcf', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset.vcf.gz', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset', '--allow-extra-chr']) # subset sample vcf file pfile generation for a given lg
                    subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole', '--set-missing-var-ids', '@:#', '--make-pgen', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected', '--allow-extra-chr'])
                    # below command will not work for cohorts less than 50 samples 
                    # pdb.set_trace()
                    # if int(subprocess.check_output(['wc', '-l', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset.psam'], encoding='utf-8').split(' ')[6]) < 50:
                    #     # subprocess.run([])
                    #     subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset', '--freq', 'counts', '--pca', 'allele-wts', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--max-alleles', '2', '--bad-ld'])
                    # else:
                    subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_subset', '--freq', 'counts', '--pca', 'allele-wts', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--max-alleles', '2'])

                    # modify the .acounts file to eliminate 0 count alleles:
                    # pdb.set_trace()
                    acounts_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.acount', sep='\t')
                    no_extremes_df = acounts_df[(acounts_df['ALT_CTS'] != 0) & (acounts_df['ALT_CTS'] != max(acounts_df['OBS_CT']))]
                    # nonzero_acounts_df = acounts_df[acounts_df['ALT_CTS'] != 0]
                    no_extremes_df.to_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.acount', sep='\t', index=False)

                    # genenrate the .sscore file with eigenvectors for each samples after projection on to the PC space from the subset PCA
                    subprocess.run(['plink2', '--pfile', self.out_dir + '/PCA/' + lg + '/' + lg + '_whole_corrected', '--read-freq', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.acount', '--score', self.out_dir + '/PCA/' + lg + '/' + lg + '_sample_subset_pca.eigenvec.allele', '2', '5', 'header-read', 'no-mean-imputation', 'variance-standardize', '--score-col-nums', '6-15', '--out', self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection', '--allow-extra-chr'])

        """
        if lg == 'Whole':
            pathlib.Path(self.out_dir + '/PCA/' + 'Whole/').mkdir(parents=True, exist_ok=True)
            print("RUNNING PLINK TO TRANSFORM THE WHOLE FILTERED VCF FILE'S DATA TO A PLINK OBJECT...")
            subprocess.run(['plink', '--vcf', self.samples_filtered_master_vcf, '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--out', self.out_dir + '/PCA/' + 'Whole/' + 'test' ]) 
            print('GENERATING EIGENVALUE AND EIGENVECTOR FILES...')
            subprocess.run(['plink', '--vcf', self.samples_filtered_master_vcf, '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--extract',  self.out_dir + '/PCA/' + 'Whole/' + 'test.prune.in', '--make-bed', '--pca', '--out',  self.out_dir + '/PCA/' + 'Whole/' + 'test'])
        else: #lg in self.linkage_group_map.values():
            pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True)
            if not pathlib.Path(self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz').exists():
                print('ERROR: THE FILE ' + lg + '.VCF.GZ DOES NOT EXIST. MUST RUN _SPLIT_VCF_TO_LG TO CREATE IT...')
                raise Exception
            else:
                print('RUNNING PLINK TO TRANSFORM THE VCF DATA TO A PLINK OBJECT...')
                subprocess.run(['plink', '--vcf', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--out', self.out_dir + '/PCA/' + lg + '/' + 'test' ]) 
                print('GENERATING EIGENVALUE AND EIGENVECTOR FILES...')
                subprocess.run(['plink', '--vcf', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--extract', self.out_dir + '/PCA/' + lg + '/' + 'test.prune.in', '--make-bed', '--pca', '--out', self.out_dir + '/PCA/' + lg + '/' + 'test'])
        """
    def _create_interactive_pca(self, linkage_group_list): # uses plotly to generate interactive PCA html outputs
        self.plotly_out = self.out_dir + '/interactive_PCA_outputs/' # define outdir
        pathlib.Path(self.plotly_out).mkdir(parents=True, exist_ok=True) # build the file path with pathlib.Path
        color_map = {'Mbuna': 'purple', 'AC': 'limegreen', 'Shallow_Benthic': 'red', 'Deep_Benthic': 'blue', 'Rhamphochromis': 'brown', 'Diplotaxodon': 'orange', 'Utaka': 'darkgreen', 'Riverine': 'pink'}
        shape_map = {'PRJEB15289': 'square', 'PRJEB1254': 'circle', 'RockSand_v1': 'diamond', 'ReferenceImprovement': 'x', 'BrainDiversity_s1': 'star', 'BigBrain': 'triangle-up', 'Multiome': 'cross'}
        
        for lg in linkage_group_list:
            eigen_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/' + lg + '_new_projection.sscore', sep='\t')
            eigen_df = eigen_df.rename(columns = {'#IID':'SampleID'})
            df_merged = pd.merge(eigen_df, self.df, on=['SampleID'])
            fig = px.scatter(df_merged, x='PC1_AVG', y='PC2_AVG', color='Ecogroup', symbol='ProjectID', color_discrete_map=color_map, symbol_map=shape_map, title=lg, hover_data=['SampleID', 'Ecogroup', 'Organism', 'ProjectID'])
            # else: # NOTE: if we want to generate a subset PCA uncomment and incorporate below code into the function
            #     fig = px.scatter(df_merged, x='PC1', y='PC2', color='Ecogroup', symbol='ProjectID', color_discrete_map=color_map, symbol_map=shape_map, title=lg, hover_data=['SampleID', 'Ecogroup', 'Organism', 'ProjectID'])
            larger_size = 10  # You can adjust this value
            fig.update_traces(marker=dict(size=larger_size), selector=dict(marker_symbol='cross'))
            fig.write_html(self.plotly_out + lg + '_plotlyPCA.html')

        
        """
        # This section will be in the test.eigenvec file per LG and generate an interactive PCA plot as an HTML file.
        # inputs: test.eigenvec per LG, SampleDatabase.xlsx file
        # Outputs: HTML file labeled per LG in in a new interactive_PCA directory
        self.plotly_out = self.out_dir + '/interactive_PCA_outputs/' # define outdir
        pathlib.Path(self.plotly_out).mkdir(parents=True, exist_ok=True) # build the file path with pathlib.Path
        color_map = {'Mbuna': 'purple', 'AC': 'limegreen', 'Shallow_Benthic': 'red', 'Deep_Benthic': 'blue', 'Rhamphochromis': 'brown', 'Diplotaxodon': 'orange', 'Utaka': 'darkgreen', 'Riverine': 'pink'}
        shape_map = {'PRJEB15289': 'square', 'PRJEB1254': 'circle', 'RockSand_v1': 'diamond', 'ReferenceImprovement': 'x', 'BrainDiversity_s1': 'star', 'BigBrain': 'triangle-up', 'Multiome': 'cross'}
        if self.df_filtered.shape[0] < 20:
            header = ['SampleID'] + ['PC{}'.format(i) for i in range(1, self.df_filtered.shape[0] +1)] # set header to 'SampleID' followed by PC1-20
        else:
            header = ['SampleID'] + ['PC{}'.format(i) for i in range(1, 21)] # set header to 'SampleID' followed by PC1-20
        for lg in linkage_group_list:
            self.eigen_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/test.eigenvec', sep=' ', header=None, index_col=0) # read in the lg's eigenvector file as a pandas dataframe
            self.eigen_df.columns = header # set the header for the eigen_df as SampleID followed by PC1-20
            self.metadata_df = pd.read_csv(self.sample_database) # read in SampleDatabase.csv
            self.metadata_df = self.metadata_df.drop_duplicates(subset='SampleID', keep='first') # remove the duplicate SampleIDs in the file and keep only the first instance
            self.df_merged = pd.merge(self.eigen_df, self.metadata_df, on=['SampleID']) # merge the dataframes on SampleID to get rid of samples not in the eigenvector file (which contains a filtered subset of samples based on eco groups provided to the script)
            if lg.startswith('NC'):
                plot_title = list(self.linkage_group_map.keys())[list(self.linkage_group_map.values()).index(lg)]
            else:
                plot_title = lg
            fig = px.scatter(self.df_merged, x='PC1', y='PC2', color='Ecogroup', symbol='ProjectID', color_discrete_map=color_map, symbol_map=shape_map, title=plot_title, hover_data=['SampleID', 'Ecogroup', 'Organism', 'ProjectID'])
            larger_size = 10  # You can adjust this value
            fig.update_traces(marker=dict(size=larger_size), selector=dict(marker_symbol='cross'))
            fig.write_html(self.plotly_out + lg + '_plotlyPCA.html')
        """
    def create_PCA(self):
        # code block of the hidden methods used to generate the PCA analysis for the PCA_Maker object
        self._create_sample_filter_files() # I think that when an object is initialized, the hidden method _create_sample_filter_file() is run automatically. This is needed so that when creating the object, a samples_filtered file will be created for use in the create_PCA method.
        self._create_ecogroup_specific_vcf()
        self._split_VCF_to_LG(self.linkage_groups)
        if args.plink:
            self._create_eigenfiles_per_LG(self.linkage_groups) # This line is used to test the _create_PCA_linakge hidden method using only LG1.
        # self._create_plots(self.linkage_groups) #commented out for now since interactive PCA plots are preferred
        self._create_interactive_pca(self.linkage_groups)


if __name__ == "__main__":
    pca_obj = PCA_Maker(args.genome, args.ecogroups, args.regions, args.output_dir)
    pca_obj.create_PCA()
    print('PIPELINE RUN SUCCESSFUL')



"""
For local testing after rehauling the pipeline 
python pca_maker.py Mzebra_UMD2a /Users/kmnike/Data/pca_testing -e Rock_Sand --local_test


CODE FOR RERUNNING ON UTAKA TO GET ALL OUTPUTS PER SAMPLE:

python3 pca_maker.py /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/original_data/pass_variants_master_file.vcf.gz /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs /home/mcgrath-lab/nkumar317/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Lake_Malawi -r Whole Exploratory All --sample_subset
python3 pca_maker.py /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/original_data/pass_variants_master_file.vcf.gz /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs /home/mcgrath-lab/nkumar317/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Non_Riverine -r Whole Exploratory All --sample_subset
python3 pca_maker.py /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/original_data/pass_variants_master_file.vcf.gz /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs /home/mcgrath-lab/nkumar317/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Mbuna Utaka Shallow_Benthic Deep_Benthic -r Whole Exploratory All --sample_subset
python3 pca_maker.py /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/original_data/pass_variants_master_file.vcf.gz /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs /home/mcgrath-lab/nkumar317/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e All -r Whole Exploratory All --sample_subset
python3 pca_maker.py /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/original_data/pass_variants_master_file.vcf.gz /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs /home/mcgrath-lab/nkumar317/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Utaka Shallow_Benthic Deep_Benthic -r Whole Exploratory All --sample_subset


local testing:
/Users/kmnike/miniforge3/envs/pipeline/bin/python3 pca_maker.py /Users/kmnike/Data/CichlidSequencingData/Outputs/small_out_files/small10percent_variants.vcf.gz /Users/kmnike/CichlidSRSequencing/576_pca_test /Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Utaka -r Whole Inversion All
/Users/kmnike/miniforge3/envs/pipeline/bin/python3 pca_maker.py /Users/kmnike/Data/CichlidSequencingData/Outputs/small_out_files/small10percent_variants.vcf.gz /Users/kmnike/CichlidSRSequencing/576_pca_test /Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Lake_Malawi -r Whole Exploratory All --sample_subset



# Old Code for generating ploits using R 
        self.r_script = os.getcwd() + '/modules/pca.R' # use the dir from which the pipeline is called to get the file path to the R script. Deprecated function since switching to plotly
        self.metadata_csv = self.out_dir + '/metadata_for_R.csv' # deprecated function. not needed with
        
        #### Below 2 lines are legacy code used to generate a PCA plot using R code  Since the pipeline has shifted to using Plotly instead, this metadata file is no longer needed and will not be gnerated anymore. If needed, add these into the _create_sample_filter_file() function
        # self.df_filtered['metadata_id'] = self.df_filtered['SampleID'] + "_" + self.df_filtered['Ecogroup']
        # self.df_filtered[['SampleID', 'metadata_id']].to_csv(self.metadata_csv, index = False)
    def _create_plots(self, linkage_group_list): # the method  has been commented out since interative PCAs are preferred. If they are ever needed, we can uncomment. I can also include a flag to determine if an interactive, static, or both plots are desired and then run the correspoding hidden methods. 
        self.pca_out = self.out_dir + '/PCA_outputs/' # define a path to an output dir where each PCA plot will go
        pathlib.Path(self.pca_out).mkdir(parents=True, exist_ok=True) # generate the filepath to said PCA_output directory
        for lg in linkage_group_list: # for each lg in the list, generate a PCA plot
            wd = self.out_dir + '/PCA/' + lg + '/' # for each iteration, this changes the working directory to the LG's eigenvalue/vector files so i dont have to name them uniquely when generating them. I can use the same prefix and just change the dir I'm working in.
            os.chdir(wd)
            # uses conda to run a script. -n specifies the env you need and the following are commands to run in that env.
            # Rscript is a command synonymous to "python3" and essentially invokes R to run the rscript. I set a path to the r_script so I don't have to hard code the filepath. I pass in the metadata file, output dir, and linkage group so I can write specific output fiel names
            subprocess.run(f"conda run -n R Rscript {self.r_script} {self.metadata_csv} {self.pca_out} {lg}", shell=True)

time bcftools view bad_lg7_pass_variants_master_file.vcf.gz --regions NC_036780.1,NC_036781.1,NC_036782.1,NC_036783.1,NC_036784.1,NC_036785.1,NC_036787.1,NC_036788.1,NC_036789.1,NC_036790.1,NC_036791.1,NC_036792.1,NC_036793.1,NC_036794.1,NC_036795.1,NC_036796.1,NC_036797.1,NC_036798.1,NC_036799.1,NC_036800.1,NC_036801.1,NC_027944.1 > pass_variants_master_file.vcf
time bgzip -c pass_variants_master_file.vcf > pass_variants_master_file.vcf.gz
time tabix -p vcf pass_variants_master_file.vcf.gz

time python pca_maker.py Mzebra_UMD2a /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs/ --sample_subset -r LG1 LG2 LG3 LG4 LG5 LG6 LG8 LG9 LG10 LG11 LG12 LG13 LG14 LG15 LG16 LG17 LG18 LG19 LG20 LG21 LG22 -e All
time python pca_maker.py Mzebra_UMD2a /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs/ --sample_subset -r LG1 LG2 LG3 LG4 LG5 LG6 LG8 LG9 LG10 LG11 LG12 LG13 LG14 LG15 LG16 LG17 LG18 LG19 LG20 LG21 LG22 -e Lake_Malawi
time python pca_maker.py Mzebra_UMD2a /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs/ --sample_subset -r LG1 LG2 LG3 LG4 LG5 LG6 LG8 LG9 LG10 LG11 LG12 LG13 LG14 LG15 LG16 LG17 LG18 LG19 LG20 LG21 LG22 -e Non_Riverine
time python pca_maker.py Mzebra_UMD2a /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs/ --sample_subset -r LG1 LG2 LG3 LG4 LG5 LG6 LG8 LG9 LG10 LG11 LG12 LG13 LG14 LG15 LG16 LG17 LG18 LG19 LG20 LG21 LG22 -e Rock_Sand
time python pca_maker.py Mzebra_UMD2a /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs/ --sample_subset -r LG1 LG2 LG3 LG4 LG5 LG6 LG8 LG9 LG10 LG11 LG12 LG13 LG14 LG15 LG16 LG17 LG18 LG19 LG20 LG21 LG22 -e Sand


python pca_maker.py Mzebra_UMD2a /Users/kmnike/Data/pca_testing --local_test -r LG8 LG9 LG19

"""
