import argparse, pdb, os, subprocess, pathlib
import pandas as pd
from cyvcf2 import VCF
import plotly.express as px
import plotly.graph_objs as go

parser = argparse.ArgumentParser(usage = "This pipeline is for running pca analysis on a filtered vcf file")
parser.add_argument('input_vcffile', help = 'absolute filepath to the filtered, gzipped input file')
parser.add_argument('output_dir', help = 'absolute filepath to an output directory')
parser.add_argument('sample_database', help = 'sample database that lists ecotype for each sample')
parser.add_argument('-e', '--ecogroups', help = 'one or multiple eco group names for filtering the data', choices = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'Riverine', 'AC', 'Non_Riverine', 'All', 'Lake_Malawi'], nargs = '*', default = ['All'])
parser.add_argument('-r', '--regions', help = 'list of linkage groups for which analyses will run', nargs = '*', default = ['All'])
args = parser.parse_args()

"""
To Do:
- The location of the pca.R script is hard coded in and assumes the pipeline will be called from the directory conatining pca_maker.py and that this directory contains the modules/pca.R script. See if this can be changed.
- The location of the intervals for LG 11 are also stored in the 'lg11_intervals' directory which the script will assume is stored in the dir from which pca_maker.py is called. If this is untrue there will be issues... 
- The script needs to be able to take in a "Whole" option for regions and generate a PCA for the whole genome's variants taken together, instead of being split by EG.
- The script also needs to be able to taek in multiple regions per LG and run analyses sequentially for all regiosn provided. The regions taht will often be called are "All," "Inversion," & "Whole"
- All and Inversion regions take the whole file, then split the genome or LG11 into subfiles in parallel, then runs plink in parallel for each and outputs eigenvec files which are used for PCA creation
    - For "Whole," the splitting step needs to be skipped, and the whole file must be used as input into PLINK. 
    - The splitting methods are _split_VCF_to_LG for "all" or lg specific analyses, and "_create_inversion_files" for "Inversion"
    - These must be skipped... The above methods take in the self.samples_filtered_master_file as inpout and output to the PCA dir. Instead of outputting a split file to a "Whole" dir in .../PCA I can nmake the "Whole" dir, then pipe outputs from plink there. 




"""

# The class PCA_Maker will create objects that will take in a variety of inputs (generally determined by what input parametrs are being passed into the script).
# These objects will have many attributes which will serve to help build directory structure, define valid inputs, etc.
# There will be many fucntions defined within the class besides the __init__ function which give objects of the class their attributes.
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

    def _create_sample_filter_file(self): # Here's a hidden function which will carry out a bunch of code in the background using the attributes defined in the __init__ block. 
        self.samples_filtered_master_vcf = self.out_dir + '/samples_filtered_master.vcf.gz'  # plink_master_vcf is an attribute that gives a filepath to an output file in the out dir. Edit to make the filepath more of what you want it to be & use pathlib to generate parent structure if it doesn't exist
        self.good_samples_csv = self.out_dir + '/samples_to_keep.csv' # This will be the name of the output file containing names of samples that have the ecogroups specified. 
        self.metadata_csv = self.out_dir + '/metadata_for_R.csv'
        # Code block to hard code in the eco group infomation and change the value of the "ecogroups" attribute depending on what the "args.ecogroups" value will be. 
        if self.ecogroups == ['All']:
            self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'Riverine', 'AC']
        elif self.ecogroups == ['Non_Riverine']:
            self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon']
        elif self.ecogroups == ['Lake_Malawi']:
            self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'AC']

        # self.s_dt[self.s_dt.Ecogroup.isin(self.ecogroups)] is returning all rows where the self.ecogroups values are contained in the "Ecogroup" column of the excel sheet
        # The .SampleID.unique() is filtering these rows to include only unique values in the SampleIDs column. However, this creates a numpy array so the pd.DataFrame wrapper around the whole thing converts this numpy array to a pandas dataframe. 
        # The to_csv(self.good_samples_txt, header = False, index = False) converts into a csv file that bcftools will be able to take as input. The index/header = False eliminate any index/header information, leaving only filtered samples. 
        # Note that the name of the output file is self.good_samples.txt which is where the file pathing for the output file is taken care of
        self.df = pd.read_excel(self.sample_database, sheet_name = 'vcf_samples') # This line generates a pandas dataframe using the "sample_database" attribute so we can use it below:
        self.df_filtered = self.df[self.df.Ecogroup.isin(self.ecogroups)]
        pd.DataFrame(self.df_filtered.SampleID.unique()).to_csv(self.good_samples_csv, header = False, index = False)
        self.df_filtered['metadata_id'] = self.df_filtered['SampleID'] + "_" + self.df_filtered['Ecogroup']
        self.df_filtered[['SampleID', 'metadata_id']].to_csv(self.metadata_csv, index = False)
        if pathlib.Path(self.samples_filtered_master_vcf).exists():
            print(f'\nThe file {self.samples_filtered_master_vcf} exists. New file will not be built.')
            pass
        else:
            print('\nGenerating a vcf file containing only samples of the EcoGroups Specified.')
            subprocess.run(['bcftools', 'view', self.in_vcf, '--samples-file', self.good_samples_csv, '-o', self.samples_filtered_master_vcf, '-O', 'z']) # code to generate a master_vcf file of filtered samples
            print('Filtered samples file generated. Indexing file...')
            subprocess.run(['tabix', '-p', 'vcf', self.samples_filtered_master_vcf]) # code to generate an index for this file using bcftools at the location of plink_master_vcf
            print('Index created...')

    def _split_VCF_to_LG(self, linkage_group_list):
        processes = []
        for lg in linkage_group_list:
            if lg == 'Whole':
                continue
            elif lg == "pre_inversion":
                self._create_inversion_files()
            elif lg in self.linkage_group_map.values():
                pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True) # generate the file paths to he split LG dirs within a dir named "PCA"
                if pathlib.Path(self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz').exists(): # For each linkage groups' dir in the PCA dir, if the LG's vcf file exists, then skip the generation of that file from the samples_filtered_master.vcf.gz file.
                    print('The file ' + lg + '.vcf.gz exists. A file for ' + lg + ' will not be generated.')
                else:
                    print('Generating a subset VCF file for ' + lg + '...')
                    p1 = subprocess.Popen(['bcftools', 'filter', '-r', lg, self.samples_filtered_master_vcf, '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '-O', 'z']) # takes in samples_filtered_master.vcf and filters out each LG and writes the file into appropriate PCA dir.
                    processes.append(p1)
                    if len(processes) == len(linkage_group_list): # parallelization code
                        for proc in processes:
                            proc.communicate()
                        processes = []

    def _create_inversion_files(self):
        processes1 = []
        processes2 = []
        inversion_regions_list = ['pre_inversion', 'inversion1', 'inversion2', 'post_inversion']
        for region in inversion_regions_list:
            self.inversion_interval_file = os.getcwd() + '/intervals_lg11/' + region + '.interval_list'
            pathlib.Path(self.out_dir + '/PCA/' + region + '/').mkdir(parents=True, exist_ok=True)
            p1 = subprocess.Popen(['gatk', 'SelectVariants', '-V', self.samples_filtered_master_vcf, '-L', self.inversion_interval_file, '-O', self.out_dir + '/PCA/' + region + '/' + region + '.vcf'])
            processes1.append(p1)
            if len(processes1) == 4:
                for proc in processes1:
                    proc.communicate()
                processes1 = []
                
        for region in inversion_regions_list:
            p2 = subprocess.Popen(['bgzip', '-f', self.out_dir + '/PCA/' + region + '/' + region + '.vcf'])
            processes2.append(p2)
            if len(processes2) == 4:
                for proc in processes2:
                    proc.communicate()
                processes2 = []

    def _create_eigenfiles_per_LG(self, linkage_group_list): # new hidden method that will create PCA plots for each LG in sample. It will define attributes for the object and also takes in a lingage grouup. Calling on this method in a for lopp should generate the eigenvalue/vector files needed per lg in self.contigs
        for lg in linkage_group_list:
            if lg == 'Whole':
                pathlib.Path(self.out_dir + '/PCA/' + 'Whole/').mkdir(parents=True, exist_ok=True)
                print("Running plink to transform the whole filtered VCF file's data to a plink object...")
                subprocess.run(['plink', '--vcf', self.samples_filtered_master_vcf, '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--out', self.out_dir + '/PCA/' + 'Whole/' + 'test' ]) 
                print('Generating eigenvalue and eigenvector files...')
                subprocess.run(['plink', '--vcf', self.samples_filtered_master_vcf, '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--extract',  self.out_dir + '/PCA/' + 'Whole/' + 'test.prune.in', '--make-bed', '--pca', '--out',  self.out_dir + '/PCA/' + 'Whole/' + 'test'])
            else:
                pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True)
                if not pathlib.Path(self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz').exists():
                    print('The file ' + lg + '.vcf.gz does not exist. Must run _split_VCF_to_LG to create it.')
                    raise Exception
                else:
                    print('Running plink to transform the VCF data to a plink object...')
                    subprocess.run(['plink', '--vcf', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--out', self.out_dir + '/PCA/' + lg + '/' + 'test' ]) 
                    print('Generating eigenvalue and eigenvector files...')
                    subprocess.run(['plink', '--vcf', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--extract', self.out_dir + '/PCA/' + lg + '/' + 'test.prune.in', '--make-bed', '--pca', '--out', self.out_dir + '/PCA/' + lg + '/' + 'test'])


    def _create_plots(self, linkage_group_list): # the method  has been commented out since interative PCAs are preferred. If they are ever needed, we can uncomment. I can also include a flag to determine if an interactive, static, or both plots are desired and then run the correspoding hidden methods. 
        self.pca_out = self.out_dir + '/PCA_outputs/' # define a path to an output dir where each PCA plot will go
        pathlib.Path(self.pca_out).mkdir(parents=True, exist_ok=True) # generate the filepath to said PCA_output directory
        for lg in linkage_group_list: # for each lg in the list, generate a PCA plot
            wd = self.out_dir + '/PCA/' + lg + '/' # for each iteration, this changes the working directory to the LG's eigenvalue/vector files so i dont have to name them uniquely when generating them. I can use the same prefix and just change the dir I'm working in.
            os.chdir(wd)
            # uses conda to run a script. -n specifies the env you need and the following are commands to run in that env.
            # Rscript is a command synonymous to "python3" and essentially invokes R to run the rscript. I set a path to the r_script so I don't have to hard code the filepath. I pass in the metadata file, output dir, and linkage group so I can write specific output fiel names
            subprocess.run(f"conda run -n R Rscript {self.r_script} {self.metadata_csv} {self.pca_out} {lg}", shell=True)

    def _create_interactive_pca(self, linkage_group_list): # uses plotly to generate interactive PCA html outputs
        # This section will ke in the test.eigenvec file per LG and generate an interactive PCA plot as an HTML file.
        # inputs: test.eigenvec per LG, SampleDatabase.xlsx file
        # Outputs: HTML file labeled per LG in in a new interactive_PCA directory
        self.plotly_out = self.out_dir + '/interactive_PCA_outputs/' # define outdir 
        pathlib.Path(self.plotly_out).mkdir(parents=True, exist_ok=True) # build the file path with pathlib.Path
        color_map = {'Mbuna': 'purple', 'AC': 'limegreen', 'Shallow_Benthic': 'red', 'Deep_Benthic': 'blue', 'Rhamphochromis': 'brown', 'Diplotaxodon': 'orange', 'Utaka': 'darkgreen', 'Riverine': 'pink'}
        shape_map = {'PRJEB15289': 'square', 'PRJEB1254': 'circle', 'RockSand_v1': 'diamond', 'ReferenceImprovement': 'x', 'BrainDiversity_s1': 'star', 'BigBrain': 'triangle-up'}
        header = ['SampleID'] + ['PC{}'.format(i) for i in range(1, 21)] # set header to 'SampleID' followed by PC1-20
        for lg in linkage_group_list:
            self.eigen_df = pd.read_csv(self.out_dir + '/PCA/' + lg + '/test.eigenvec', sep=' ', header=None, index_col=0) # read in the lg's eigenvector file as a pandas dataframe
            self.eigen_df.columns = header # set the header for the eigen_df as SampleID followed by PC1-20
            self.metadata_df = pd.read_excel(self.sample_database, sheet_name='vcf_samples') # read in SampleDatabase.xlsx 
            self.metadata_df = self.metadata_df.drop_duplicates(subset='SampleID', keep='first') # remove the duplicate SampleIDs in the file and keep only the first instance
            self.df_merged = pd.merge(self.eigen_df, self.metadata_df, on=['SampleID']) # merge the dataframes on SampleID to get rid of samples not in the eigenvector file (which contains a filtered subset of samples based on eco groups provided to the script)
            if lg.startswith('NC'):
                plot_title = list(self.linkage_group_map.keys())[list(self.linkage_group_map.values()).index(lg)]
            else:
                plot_title = lg
            fig = px.scatter(self.df_merged, x='PC1', y='PC2', color='Ecogroup', symbol='ProjectID', color_discrete_map=color_map, symbol_map=shape_map, title=plot_title, hover_data=['SampleID', 'Ecogroup', 'Organism', 'ProjectID'])
            fig.write_html(self.plotly_out + lg + '_plotlyPCA.html')

    def create_PCA(self):
        # code block of the hidden methods used to generate the PCA analysis for the PCA_Maker object
        self._create_sample_filter_file() # I think that when an object is initialized, the hidden method _create_sample_filter_file() is run automatically. This is needed so that when creating the object, a samples_filtered file will be created for use in the create_PCA method.
        self._split_VCF_to_LG(self.linkage_groups)
        self._create_eigenfiles_per_LG(self.linkage_groups) # This line is used to test the _create_PCA_linakge hidden method using only LG1.
        # self._create_plots(self.linkage_groups) #commented out for now since interactvie PCA plots are preferred
        self._create_interactive_pca(self.linkage_groups)

pca_obj = PCA_Maker(args.input_vcffile, args.output_dir, args.sample_database, args.ecogroups, args.regions)
pca_obj.create_PCA()
print('PIPELINE RUN SUCCESSFUL')

"""
CODE FOR RERUNNING ON UTAKA TO GET ALL OUTPUTS PER SAMPLE:

python3 pca_maker.py /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pass_variants_master_file.vcf.gz /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs /home/mcgrath-lab/nkumar317/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Lake_Malawi -r Whole All Inversion
python3 pca_maker.py /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pass_variants_master_file.vcf.gz /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs /home/mcgrath-lab/nkumar317/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Non_Riverine -r Whole All Inversion
python3 pca_maker.py /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pass_variants_master_file.vcf.gz /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs /home/mcgrath-lab/nkumar317/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Mbuna Utaka Shallow_Benthic Deep_Benthic -r Whole All Inversion
python3 pca_maker.py /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pass_variants_master_file.vcf.gz /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs /home/mcgrath-lab/nkumar317/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e All -r Whole All Inversion
python3 pca_maker.py /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pass_variants_master_file.vcf.gz /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/pca_outputs /home/mcgrath-lab/nkumar317/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Utaka Shallow_Benthic Deep_Benthic -r Whole All Inversion




lcoal testing:
/Users/kmnike/anaconda3/envs/pipeline/bin/python3 pca_maker.py /Users/kmnike/Data/CichlidSequencingData/Outputs/small_out_files/small_.1_variants.vcf.gz /Users/kmnike/CichlidSRSequencing/576_pca_test/ /Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Mbuna -r Whole Inversion All


The 'symbol' property is an enumeration that may be specified as:
      - One of the following enumeration values:
            [0, '0', 'circle', 100, '100', 'circle-open', 200, '200',
            'circle-dot', 300, '300', 'circle-open-dot', 1, '1',
            'square', 101, '101', 'square-open', 201, '201',
            'square-dot', 301, '301', 'square-open-dot', 2, '2',
            'diamond', 102, '102', 'diamond-open', 202, '202',
            'diamond-dot', 302, '302', 'diamond-open-dot', 3, '3',
            'cross', 103, '103', 'cross-open', 203, '203',
            'cross-dot', 303, '303', 'cross-open-dot', 4, '4', 'x',
            104, '104', 'x-open', 204, '204', 'x-dot', 304, '304',
            'x-open-dot', 5, '5', 'triangle-up', 105, '105',
            'triangle-up-open', 205, '205', 'triangle-up-dot', 305,
            '305', 'triangle-up-open-dot', 6, '6', 'triangle-down',
            106, '106', 'triangle-down-open', 206, '206',
            'triangle-down-dot', 306, '306', 'triangle-down-open-dot',
            7, '7', 'triangle-left', 107, '107', 'triangle-left-open',
            207, '207', 'triangle-left-dot', 307, '307',
            'triangle-left-open-dot', 8, '8', 'triangle-right', 108,
            '108', 'triangle-right-open', 208, '208',
            'triangle-right-dot', 308, '308',
            'triangle-right-open-dot', 9, '9', 'triangle-ne', 109,
            '109', 'triangle-ne-open', 209, '209', 'triangle-ne-dot',
            309, '309', 'triangle-ne-open-dot', 10, '10',
            'triangle-se', 110, '110', 'triangle-se-open', 210, '210',
            'triangle-se-dot', 310, '310', 'triangle-se-open-dot', 11,
            '11', 'triangle-sw', 111, '111', 'triangle-sw-open', 211,
            '211', 'triangle-sw-dot', 311, '311',
            'triangle-sw-open-dot', 12, '12', 'triangle-nw', 112,
            '112', 'triangle-nw-open', 212, '212', 'triangle-nw-dot',
            312, '312', 'triangle-nw-open-dot', 13, '13', 'pentagon',
            113, '113', 'pentagon-open', 213, '213', 'pentagon-dot',
            313, '313', 'pentagon-open-dot', 14, '14', 'hexagon', 114,
            '114', 'hexagon-open', 214, '214', 'hexagon-dot', 314,
            '314', 'hexagon-open-dot', 15, '15', 'hexagon2', 115,
            '115', 'hexagon2-open', 215, '215', 'hexagon2-dot', 315,
            '315', 'hexagon2-open-dot', 16, '16', 'octagon', 116,
            '116', 'octagon-open', 216, '216', 'octagon-dot', 316,
            '316', 'octagon-open-dot', 17, '17', 'star', 117, '117',
            'star-open', 217, '217', 'star-dot', 317, '317',
            'star-open-dot', 18, '18', 'hexagram', 118, '118',
            'hexagram-open', 218, '218', 'hexagram-dot', 318, '318',
            'hexagram-open-dot', 19, '19', 'star-triangle-up', 119,
            '119', 'star-triangle-up-open', 219, '219',
            'star-triangle-up-dot', 319, '319',
            'star-triangle-up-open-dot', 20, '20',
            'star-triangle-down', 120, '120',
            'star-triangle-down-open', 220, '220',
            'star-triangle-down-dot', 320, '320',
            'star-triangle-down-open-dot', 21, '21', 'star-square',
            121, '121', 'star-square-open', 221, '221',
            'star-square-dot', 321, '321', 'star-square-open-dot', 22,
            '22', 'star-diamond', 122, '122', 'star-diamond-open',
            222, '222', 'star-diamond-dot', 322, '322',
            'star-diamond-open-dot', 23, '23', 'diamond-tall', 123,
            '123', 'diamond-tall-open', 223, '223',
            'diamond-tall-dot', 323, '323', 'diamond-tall-open-dot',
            24, '24', 'diamond-wide', 124, '124', 'diamond-wide-open',
            224, '224', 'diamond-wide-dot', 324, '324',
            'diamond-wide-open-dot', 25, '25', 'hourglass', 125,
            '125', 'hourglass-open', 26, '26', 'bowtie', 126, '126',
            'bowtie-open', 27, '27', 'circle-cross', 127, '127',
            'circle-cross-open', 28, '28', 'circle-x', 128, '128',
            'circle-x-open', 29, '29', 'square-cross', 129, '129',
            'square-cross-open', 30, '30', 'square-x', 130, '130',
            'square-x-open', 31, '31', 'diamond-cross', 131, '131',
            'diamond-cross-open', 32, '32', 'diamond-x', 132, '132',
            'diamond-x-open', 33, '33', 'cross-thin', 133, '133',
            'cross-thin-open', 34, '34', 'x-thin', 134, '134',
            'x-thin-open', 35, '35', 'asterisk', 135, '135',
            'asterisk-open', 36, '36', 'hash', 136, '136',
            'hash-open', 236, '236', 'hash-dot', 336, '336',
            'hash-open-dot', 37, '37', 'y-up', 137, '137',
            'y-up-open', 38, '38', 'y-down', 138, '138',
            'y-down-open', 39, '39', 'y-left', 139, '139',
            'y-left-open', 40, '40', 'y-right', 140, '140',
            'y-right-open', 41, '41', 'line-ew', 141, '141',
            'line-ew-open', 42, '42', 'line-ns', 142, '142',
            'line-ns-open', 43, '43', 'line-ne', 143, '143',
            'line-ne-open', 44, '44', 'line-nw', 144, '144',
            'line-nw-open', 45, '45', 'arrow-up', 145, '145',
            'arrow-up-open', 46, '46', 'arrow-down', 146, '146',
            'arrow-down-open', 47, '47', 'arrow-left', 147, '147',
            'arrow-left-open', 48, '48', 'arrow-right', 148, '148',
            'arrow-right-open', 49, '49', 'arrow-bar-up', 149, '149',
            'arrow-bar-up-open', 50, '50', 'arrow-bar-down', 150,
            '150', 'arrow-bar-down-open', 51, '51', 'arrow-bar-left',
            151, '151', 'arrow-bar-left-open', 52, '52',
            'arrow-bar-right', 152, '152', 'arrow-bar-right-open', 53,
            '53', 'arrow', 153, '153', 'arrow-open', 54, '54',
            'arrow-wide', 154, '154', 'arrow-wide-open']












"""

#### I htink that things are getting messed up in the "_split_VCF_to_LG" method for the inversion regions since no values are being written to the VCF files
