import subprocess as sp, argparse, shlex, os
from helper_modules.nikesh_file_manager import FileManager as FM
import pdb

"""
USE THE process_vcf.py SCRIPT INSTEAD. IT INCLUDES ALL BELOW CODE ALONG WITH VCF STATS AND FILTERING CODE  
The script will assume that the individual VCF files per chromosome are stored on the Utaka server itself in the fm_obj.localOutputDir.
Within this directory, I need to find files that match '*.vcf' and concatenate these only  

Functionality has been added to allow inidividual LG VCF files stored in the storage dir in Utaka to be read in and used for concatenation.
I will add functionality that will automatically zip and index the VCF file following creation of the master file. 
"""

parser = argparse.ArgumentParser(usage = 'This script will concatenate VCF ouputs from the CallSmallSNVs.py pipeline into a master VCF file.')
parser.add_argument('reference_genome', type = str, help = 'Version of the genome use')
parser.add_argument('-u', '--utaka', help = 'use this flag when concatenating on the Utaka server and the individual files are stored in the /Output directory', action = 'store_true')
parser.add_argument('-z', '--zip', help = 'call on this flag to zip the concatenated masterfile', action = 'store_true')
parser.add_argument('-i', '--index', help = 'call on this flag to index the zipped masterfile', action = 'store_true')
parser.add_argument('--local_test', help = 'call this flag to predefine variables for testing on local machine', action='store_true')
args = parser.parse_args()

class VCF_Concatenator:
    def __init__(self, genome):
        self.fm_obj = FM(args.reference_genome)
        self.genome = self.fm_obj.localGenomeFile
        self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1', 
                             'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1', 
                             'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1', 
                             'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1', 'mito': 'NC_027944.1'}
        self.linkage_groups = list(self.linkage_group_map.values())
        if args.local_test:
            self.linkage_groups = ['NC_036787.1', 'NC_036788.1', 'NC_036798.1']

    def _make_file_structure(self):
        try:
            os.makedirs(self.fm_obj.localOutputDir + 'vcf_concat_output/') # make a vcf_concat_output directory 
        except FileExistsError:
            pass
        self.output_dir = self.fm_obj.localOutputDir + 'vcf_concat_output/' # set the output_dir to the above created value
        
        #### Quick check to see if a master vcf file has already been created. 
        self.file_num = 0 # initialize a file_num variable. This will ensure a new VCF concatenation file will be generated each time to avoid accidental overwrite when concatenating multiple files or data sources 
        for file in os.listdir(self.output_dir):
            self.file_num += 1

    def _concat_vcfs(self):
        #### if no master file exists, run through the code normally and name the file "master_file.vcf"
        if self.file_num == 0:
            with open(f"{self.output_dir}master_file.vcf", 'w+') as f1: # open the master list file if the vcf_concat_output dir is empty
                for file in self.linkage_groups:
                    print('Reading in and concatenating ' + file + '...')
                    if args.utaka:
                        file_to_write = self.fm_obj.StorageOutputDir + file + '_output.vcf'
                    else:
                        file_to_write = self.fm_obj.localOutputDir + file + '_output.vcf'
                    if os.path.getsize(f"{self.output_dir}master_file.vcf") == 0: # if the master file is empty, take the file and write the whole contents. This will write the header
                        f1.write(sp.check_output(shlex.split(f"cat {file_to_write}"), encoding='utf-8')) # the newline is needed at the end
                    else: # if file has contents, open the file in read mode then go through each line. If it doesnt start with a "#" then write it to the master file.
                        with open(file_to_write, 'r') as f2:
                            for line in f2:
                                if not line.startswith('#'): # avoids writing header information from all subsequent files 
                                    f1.write(line)
        else: # if the master file exists, create a new file so that the data in previous iterations will be preserved to prevent overwrite of an existing file 
            with open(f"{self.output_dir}master_file{self.file_num}.vcf", 'w+') as f1: # open a new master_file.vcf with a new file_num if one already exists
                for file in self.linkage_groups:
                    print('Reading in and concatenating ' + file + '...')
                    if args.utaka:
                        file_to_write = self.fm_obj.StorageOutputDir + file + '_output.vcf'
                    else:
                        file_to_write = file_to_write = self.fm_obj.localOutputDir + file + '_output.vcf'
                    if os.path.getsize(f"{self.output_dir}master_file{self.file_num}.vcf") == 0: # if the master file is empty, take the file and write the whole contents. This will write the header
                        f1.write(sp.check_output(shlex.split(f"cat {file_to_write}"), encoding='utf-8')) # the newline is needed at the end
                    else: # if file has contents, open the file in read mode then go through each line. If it doesnt start with a "#" then write it to the master file.
                        with open(file_to_write, 'r') as f2:
                            for line in f2:
                                if not line.startswith('#'): # avoids writing header information from all subsequent files 
                                    f1.write(line)
        print('MASTER_FILE GENERATED')

    def zip_masterfile(self):
        print('Starting master_file compression...')
        if self.file_num == 0:
            with open(f"{self.output_dir}master_file.vcf.gz", 'w+') as f1:
                # I've confirmed that the below command gives a zipped file that can be unzipped and has the same data as the original using the diff command and using ls -l and confirming that the bytes are the same
                sp.call(shlex.split(f"bgzip -c {self.output_dir}master_file.vcf"), stdout=f1)
        else:
            with open(f"{self.output_dir}master_file{self.file_num}.vcf.gz", 'w') as f2:
                sp.call(shlex.split(f"bgzip -c {self.output_dir}master_file.vcf"), stdout=f2)
        print('MASTER_FILE COMPRESSION COMPLETE')

    def index_masterfile(self):
        print('Generating index using zipped master_file...')
        if self.file_num == 0:
            sp.run(shlex.split(f"tabix -p vcf {self.output_dir}master_file.vcf.gz"))
        else:
            sp.run(shlex.split(f"tabix -p vcf {self.output_dir}master_file{self.file_num}.vcf.gz"))
        print('MASTER_FILE INDEXING COMPLETE')
    def run_methods(self):
        self._make_file_structure()
        self._concat_vcfs()
        if args.zip:
            self.zip_masterfile()
        if args.index:
            self.index_masterfile()

concat_obj = VCF_Concatenator(args.reference_genome)
concat_obj.run_methods()
print("Pipeline run complete")

"""
Below is the old code KEEP THE CODE TO PRESVERVE CODE THAT ALLOWS FOR DOWNLOAD AND REMOVE CONCATENATION IF FILES ARE STORED ONLY ON DROPBOX 


#### use rclone lsf to create a list of existing vcf files
outputs = sp.run(shlex.split("rclone lsf ptm_dropbox:BioSci-McGrath/Apps/CichlidSequencingData/Outputs/Mzebra_UMD2a/Outputs576Cohort/ --include '*.vcf'"), stdout=sp.PIPE, encoding='utf-8').stdout.splitlines()
path_to_outputs = 'ptm_dropbox:BioSci-McGrath/Apps/CichlidSequencingData/Outputs/Mzebra_UMD2a/Outputs576Cohort/'

#### check if an output dir exists. If not, make it.
try:
    os.makedirs(os.getcwd() + '/vcf_concat_output/')
except FileExistsError:
    pass
output_dir = os.getcwd() + '/vcf_concat_output/'

#### make a Temp dir to temporarily store the VCF files as they're concatenated
try:
    os.makedirs(os.getcwd() + '/Temp/')
except FileExistsError:
    pass
temp_dir = os.getcwd() + '/Temp/'

#### Quick check to see if a master vcf file has already been created. 
file_num = 0
for file in os.listdir(output_dir):
    file_num += 1

#### if no master file exists, run through the code normally and name the file "master_file.vcf"
if file_num == 0:
    with open(f"{output_dir}master_file.vcf", 'w+') as f1: # open the master list file if the vcf_concat_output dir is empty
        for file in outputs: # the output file names were defined above. Iterate through each file
            sp.run(shlex.split(f"rclone copy {path_to_outputs + file} {temp_dir} -P")) # download the file to a temp dir
            file_to_write = temp_dir + file 
            if os.path.getsize(f"{output_dir}master_file.vcf") == 0: # if the master file is empty, take the file and write the whole contents. This will write the header
                f1.write(sp.check_output(shlex.split(f"cat {file_to_write}"), encoding='utf-8') + '\n') #the newline is needed at the end
            else: # if file has contents, open the file in read mode then go through each line. If it doesnt start with a "#" then write it to the master file.
                with open(temp_dir + file, 'r') as f2:
                    for line in f2:
                        if not line.startswith('#'):
                            f1.write(line)
                    f1.write('\n') # last newline needed to not merge lines between files together. 
            sp.run(shlex.split(f"rm {file_to_write}"))
            sp.run(shlex.split(f"rm  -r {temp_dir}"))
else: # if the master file exists, create a new file so that the data in previous iterations will be preserved. 
    with open(f"{output_dir}master_file{file_num}.vcf", 'w+') as f1: #open the master list file if the vcf_concat_output dir is empty
            for file in outputs: # the output file names were defined above. Iterate through each file
                sp.run(shlex.split(f"rclone copy {path_to_outputs + file} {temp_dir} -P")) # download the file to a temp dir
                file_to_write = temp_dir + file 
                if os.path.getsize(f"{output_dir}master_file{file_num}.vcf") == 0: # if the master file is empty, take the file and write the whole contents. This will write the header
                    f1.write(sp.check_output(shlex.split(f"cat {file_to_write}"), encoding='utf-8') + '\n') #the newline is needed at the end
                else: # if file has contents, open the file in read mode then go through each line. If it doesnt start with a "#" then write it to the master file. 
                    with open(temp_dir + file, 'r') as f2:
                        for line in f2:
                            if not line.startswith('#'):
                                f1.write(line)
                        f1.write('\n') # last newline needed to not merge lines between files together.
                sp.run(shlex.split(f"rm {file_to_write}"))
                sp.run(shlex.split(f"rm  -r {temp_dir}"))
"""