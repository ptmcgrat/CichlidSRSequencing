import argparse, pdb, os, subprocess, pathlib
import pandas as pd
from cyvcf2 import VCF

parser = argparse.ArgumentParser(usage = "This pipeline is for running pca analysis on a filtered vcf file")
parser.add_argument('input_vcffile', help = 'absolute filepath to the filtered, gzipped input file')
parser.add_argument('output_dir', help = 'absolute filepath to an output directory')
parser.add_argument('sample_database', help = 'sample database that lists ecotype for each sample')
parser.add_argument('-e', '--ecogroups', help = 'one or multiple eco group names for filtering the data', choices = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhampochromis', 'Diplotaxodon', 'Riverine', 'AC', 'Non_Riverine', 'All'], nargs = '*', default = ['All'])
parser.add_argument('-r', '--regions', help = 'list of linkage groups for which analyses will run', nargs = '*', default = ['All'])
parser.add_argument('--PCA', help = 'generate a PCA analysis for the specified linkage groups', default = "All") # this shoudln't be needed since the script is designed to generate the PCA plot
parser.add_argument('-f', '--filters', help = 'list of tunable parametrs for filtering the raw input vcf file.', default  = ['DP > 11000', 'DP < 8000', 'InbreedingCoeff < -0.6', 'FS > 40.0', 'QD < 2.0', 'NCC> 125', 'MQ < 50', 'AF < 0.000958']) #remove this an implement anotehr script that can be used to filter the raw input file at a later time
args = parser.parse_args()

# The class PCA_Maker will create objects that will take in a variety of inputs (generally determined by what input parametrs are being passed into the script). 
# These objects will have many attributes which will serve to help build directory structure, define valid inputs, etc.
# There will be many fucntions defined within the class besides the __init__ function whichgive objects of the class their attributes. Additional hidden functions will serve to 
class PCA_Maker:
    def __init__(self, input_vcffile, output_directory, sample_database, ecogroups, linkage_groups): # The PCA_Maker class will create an object (self) which will then take in (currently) 4 pieces of information (input file out dir, sample_database excel sheet, ecogroup names)
        # self.attr indicates that the object made with PCA_Maker will have the attribute named "attr"
        # The linkage_group_map attribute is a hard coded list of LG names
        self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1', 
                             'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1', 
                             'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1', 
                             'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1'}
        self.in_vcf = input_vcffile # The in_vcf attriubute is sequal to the input file name for the Object. 
        self.out_dir = output_directory# The out_dir attribute equals the output_directory name for the object
        self.sample_database = sample_database # The sample_database attribute equals the sample_database name for the object
        self.ecogroups = ecogroups # This attribute is the list of ecgogroups used for filtering samples
        self.r_script = os.getcwd() + '/modules/pca.R' # use the dir from which the pipeline is called to get the file path to the R script. 
        
        self.vcf_obj = VCF(self.in_vcf) # The VCF object is made using the CVF class from cyvcf2. The VCF class takes an input vcf. For the object, this input vcf file is the "input_vcfcfile" which is defined under self.in_vcf 
        self.linkage_groups = linkage_groups # the object will now take in linkage groups using args.regions. If default, it will default to the first 22 lgs
        if self.linkage_groups == ['All']:
            self.linkage_groups = self.vcf_obj.seqnames[0:22]
        else: # else statement takes in LG names as just "LG1, LG2, etc. then converts them to the actual contig names and assigns these to self.linkage_groups"
            remapped_linkage_groups = []
            for lg,contig_name in self.linkage_group_map.items():
                if lg in self.linkage_groups:
                    remapped_linkage_groups.append(contig_name)
            self.linkage_groups = remapped_linkage_groups
        # Makes sure LGs are in vcf file provided
        for lg in self.linkage_groups: # Code to ensure each translated lg value is in the hard coded dictionary
            assert lg in self.linkage_group_map.values() # assert checks if the lg names in self.linkage_groups are in the dictionary.If anything doesn't match, this assertion fails and an error is thrown

        # Ensure index file exists
        assert os.path.exists(self.in_vcf + '.tbi') # uses os.path.exists to see if the input file + 'tbi' extension exists. The object will be made using args.input_vcffile and args.input_vcffile will be passed to the script as an absolute file path so the path to the dir is taken care of 

        pathlib.Path(output_directory).mkdir(parents=True, exist_ok=True) # This generates the outdir if it doesn't exist so later tools don't run into errors making files.
        self._create_sample_filter_file() # I think that when an object is initialized, the hidden method _create_sample_filter_file() is run automatically. This is needed so that when creating the object, a samples_filtered file will be created for use in the create_PCA method.
        self._create_PCA_per_LG(self.linkage_groups) # This line is used to test the _create_PCA_linakge magic method using only LG1. 
        self._create_plots(self.linkage_groups)

    def _create_sample_filter_file(self): # Here's a hidden function which will carry out a bunch of code in the background using the attributes defined in the __init__ block. 
        self.samples_filtered_master_vcf = self.out_dir + '/samples_filtered_master.vcf.gz'  # plink_master_vcf is an attribute that gives a filepath to an output file in the out dir. Edit to make the filepath more of what you want it to be & use pathlib to generate parent structure if it doesn't exist
        self.good_samples_csv = self.out_dir + '/samples_to_keep.csv' # This will be the name of the output file containing names of samples that have the ecogroups specified. 
        self.metadata_csv = self.out_dir + '/filtered_samples.csv'
        # Code block to hard code in the eco group infomation and change the value of the "ecogroups" attribute depending on what the "args.ecogroups" value will be. 
        if self.ecogroups == ['All']:
            self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhampochromis', 'Diplotaxodon', 'Riverine', 'AC']
        elif self.ecogroups == ['Non_Riverine']:
            self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhampochromis', 'Diplotaxodon']

        self.df = pd.read_excel(self.sample_database, sheet_name = 'vcf_samples') # This line generates a pandas dataframe using the "sample_database" attribute so we can use it below:
        # lot going on with the below lines.
        # self.s_dt[self.s_dt.Ecogroup.isin(self.ecogroups)] is returning all rows where the self.ecogroups values are contained in the "Ecogroup" column of the excel sheet
        # The .SampleID.unique() is filtering these rows to include only unique values in the SampleIDs column. However, this creates a numpy array so the pd.DataFrame wrapper around the whole thing converts this numpy array to a pandas dataframe. 
        # The to_csv(self.good_samples_txt, header = False, index = False) converts into a csv file that bcftools will be able to take as input. The index/header = False eliminate any index/header information, leaving only filtered samples. 
        # Note that the name of the output file is self.good_samples.txt which is where the file pathing for the output file is taken care of
        self.df_filtered = self.df[self.df.Ecogroup.isin(self.ecogroups)]
        pd.DataFrame(self.df_filtered.SampleID.unique()).to_csv(self.good_samples_csv, header = False, index = False)
        self.df_filtered['metadata_id'] = self.df_filtered['SampleID'] + "_" + self.df_filtered['Ecogroup']
        self.df_filtered[['SampleID', 'metadata_id']].to_csv(self.metadata_csv, index = False)
        subprocess.run(['bcftools', 'view', self.in_vcf, '--samples-file', self.good_samples_csv, '-o', self.samples_filtered_master_vcf, '-O', 'z']) # code to generate a master_vcf file of filtered samples
        subprocess.run(['tabix', '-p', 'vcf', self.samples_filtered_master_vcf]) # code to generate an index for this file using bcftools at the location of plink_master_vcf

    def _create_PCA_per_LG(self, linkage_group_list): # new magic method that will create PCA plots for each LG in sample. It will define attributes for the object and also takes in a lingage grouup. Calling on htis method in a for lopp should generate the eigenvalue/vector files needed per lg in self.contigs
        for lg in linkage_group_list:
            pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True)
            subprocess.run(['bcftools', 'filter', '-r', lg, self.samples_filtered_master_vcf, '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '-O', 'z'])
            subprocess.run(['plink', '--vcf', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--out', self.out_dir + '/PCA/' + lg + '/' + 'test' ])
            subprocess.run(['plink', '--vcf', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--extract', self.out_dir + '/PCA/' + lg + '/' + 'test.prune.in', '--make-bed', '--pca', '--out', self.out_dir + '/PCA/' + lg + '/' + 'test'])

    def _create_plots(self, linkage_group_list):
        self.pca_out = self.out_dir + '/PCA_outputs/' # define a path to an output dir where each PCA plot will go
        pathlib.Path(self.pca_out).mkdir(parents=True, exist_ok=True) # generate the filepath to said PCA_output directory
        for lg in linkage_group_list: # for each lg in the list, generate a PCA plot
            wd = self.out_dir + '/PCA/' + lg + '/' # for each iteration, this changes the working directory to the LG's eigenvalue/vector files so i dont have to name them uniquely when generating them. I can use the same prefix and just change the dir I'm working in.
            os.chdir(wd)
            # uses conda to run a script. -n specifies the env you need and the following are commands to run in that env.
            # Rscript is a command synonymous to "python3" and essentially invokes R to run the rscript. I set a path to the r_script so I don't have to hard code the filepath. I pass in the metadata file, output dir, and linkage group so I can write specific output fiel names
            subprocess.run(f"conda run -n R Rscript {self.r_script} {self.metadata_csv} {self.pca_out} {lg}", shell=True) 

pca_obj = PCA_Maker(args.input_vcffile, args.output_dir, args.sample_database, args.ecogroups, args.regions)


"""
TEST THE CODE ON SERVER:
python3 pca_maker.py /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/small_test_files/small_lg1-22_master_file.vcf.gz ~/Test ~/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Non_Riverine -r LG11 LG15

RUN CODE FOR FILTERED VCF FILE:
python3 pca_maker.py /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/original_data/filtered_variants_v1.vcf.gz /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/pipeline_outpus /home/ad.gatech.edu/bio-mcgrath-dropbox/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Non_Riverine

TEST THE CODE LOCALLY
~/anaconda3/envs/mcgrath/bin/python3 pca_maker.py ~/Data/CichlidSequencingData/Pipeline/raw_data/small_lg1-22_master_file.vcf.gz ~/Test ~/CichlidSRSequencing/cichlid_sr_sequencing/SampleDatabase.xlsx -e Non_Riverine
"""