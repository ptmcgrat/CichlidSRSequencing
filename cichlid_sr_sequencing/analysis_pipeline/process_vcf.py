import argparse, pdb, glob, shlex, os, datetime, time
from helper_modules.nikesh_file_manager import FileManager as FM
import subprocess
from multiprocessing import Process

parser = argparse.ArgumentParser(usage = 'Once VCF files per chromopsome have been generated, this script will combine the VCFs per LG into one file then optionally perfform compression, filtering, etc.')
parser.add_argument('genome', help = 'name of reference genome used in the creation of the VCF files')
parser.add_argument('-m', '--merge', help = 'this flag will merge the VCF outputs into a master file', action = 'store_true')
parser.add_argument('-c', '--compress_and_index', help = 'call this flag to compress and index the output vcf file using bgzip, then tabix', action = 'store_true')
parser.add_argument('-v', '--vcf_stats', help = 'this flag will run vcftools to generate vcf stats for a zipped vcf file', action = 'store_true')
parser.add_argument('-p', '--prefix', help = 'define a prefix for filtering stats files', default = ['cohort_stats'], nargs = 1, type = str)
parser.add_argument('-u', '--utaka', help = 'use this flag when concatenating on the Utaka server and the individual files are stored in the /Output directory', action = 'store_true')
parser.add_argument('-f', '--filter', help = 'this flag will invoke filtering of the merged vcf file', action = 'store_true')
parser.add_argument('--local_test', help = 'call this flag to predefine variables for testing on local machine', action='store_true')
args = parser.parse_args()

"""
TODO: 
This script will ideally take the output from the callVariants.py pipeline and perform concatenation, compression, indexing, stats generation, and filtering.
Ideally it will incorporate the rclone uploading as well.

- implement gatk VariantFiltration
- implement rclone uploads at the end of every piece of the pipeline, but use a flag to dictate whether the upload will happen
    - be sure to implement chunker uplaods for massive output files.
- I believe a final pass_variant file indexing step is needed. 

time python process_vcf.py -m -c -v
"""

class VCFProcessor:
    def __init__(self, genome):
        self.genome = genome
        self.fm_obj = FM(self.genome)
        self.master_file = self.fm_obj.localOutputDir + 'vcf_concat_output/master_file.vcf'
        self.zipped_master_file = self.master_file + '.gz'
        self.stats_dir = self.fm_obj.localOutputDir + 'filtering_stats/'
        # getting rid of linkage_group code since it doesn;t apply to the new interval based processing 2024.07.02
        # self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1', 
        #                      'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1', 
        #                      'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1', 
        #                      'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1', 'mito': 'NC_027944.1'}
        # self.linkage_groups = list(self.linkage_group_map.values())
        intervals = list(range(1,97))
        self.intervals = list(map(str, intervals))
        self.current_time = datetime.datetime.now()

        if args.local_test:
            self.linkage_groups = ['NC_036787.1', 'NC_036788.1', 'NC_036798.1'] # old code based on LG based file processing. Now, even for local_test, intervals will get processed. 

    def _make_file_structure(self):
        # check if vcf_concat_output dir exists and make if it doesn't exist 
        try:
            os.makedirs(self.fm_obj.localOutputDir + 'vcf_concat_output/') # make a vcf_concat_output directory 
        except FileExistsError:
            pass
        # check if filtering_stats dir exists and make if it doesnt't exist
        try:
            os.makedirs(self.fm_obj.localOutputDir + 'filtering_stats/')
        except FileExistsError:
            pass
        self.vcfConcatDir = self.fm_obj.localOutputDir + 'vcf_concat_output/'
        self.file_num = 0
        for file in os.listdir(self.vcfConcatDir):
            if file.endswith('vcf'): # added so any index or zipped files in self.vcfConcatDir do no increment the counter
                self.file_num += 1

    def merge_vcfs(self):
        print(f"VCF Master File concatenation started at {self.current_time}")
    #### if no master file exists, run through the code normally and name the file "master_file.vcf"
        if self.file_num == 0:
            with open(f"{self.vcfConcatDir}master_file.vcf", 'w+') as f1: # open the master list file if the vcf_concat_output dir is empty
                for interval in self.intervals:
                    print(f"Concatenation of interval {interval} started at {self.current_time}")
                    if args.utaka:
                        file_to_write = self.fm_obj.StorageOutputDir + interval + '_output.vcf'
                    else:
                        file_to_write = self.fm_obj.localOutputDir + interval + '_output.vcf'
                    if os.path.getsize(f"{self.vcfConcatDir}master_file.vcf") == 0: # if the master file is empty, take the file and write the whole contents. This will write the header
                        f1.write(subprocess.check_output(shlex.split(f"cat {file_to_write}"), encoding='utf-8')) # the newline is needed at the end
                    else: # if file has contents, open the file in read mode then go through each line. If it doesnt start with a "#" then write it to the master file.
                        with open(file_to_write, 'r') as f2:
                            for line in f2:
                                if not line.startswith('#'): # avoids writing header information from all subsequent files 
                                    f1.write(line)
                    print(f"Concatenation of interval {interval} cpmplete at {self.current_time}")
        else: # if the master file exists, create a new file so that the data in previous iterations will be preserved to prevent overwrite of an existing file 
            with open(f"{self.vcfConcatDir}master_file{self.file_num}.vcf", 'w+') as f1: # open a new master_file.vcf with a new file_num if one already exists
                for interval in self.intervals:
                    print(f"Concatenation of interval {interval} started at {self.current_time}")
                    if args.utaka:
                        file_to_write = self.fm_obj.StorageOutputDir + interval + '_output.vcf'
                    else:
                        file_to_write = file_to_write = self.fm_obj.localOutputDir + interval + '_output.vcf'
                    if os.path.getsize(f"{self.vcfConcatDir}master_file{self.file_num}.vcf") == 0: # if the master file is empty, take the file and write the whole contents. This will write the header
                        f1.write(subprocess.check_output(shlex.split(f"cat {file_to_write}"), encoding='utf-8')) # the newline is needed at the end
                    else: # if file has contents, open the file in read mode then go through each line. If it doesnt start with a "#" then write it to the master file.
                        with open(file_to_write, 'r') as f2:
                            for line in f2:
                                if not line.startswith('#'): # avoids writing header information from all subsequent files 
                                    f1.write(line)
                    print(f"Concatenation of interval {interval} complete at {self.current_time}")
        print(f"VCF Master File concatenation complete at {self.current_time}")

    def compress_and_index_vcf(self):
        print(f"master_file.vcf compression staretd at {self.current_time}")
        if self.file_num == 0:
            with open(f"{self.vcfConcatDir}master_file.vcf.gz", 'w+') as f1:
                # I've confirmed that the below command gives a zipped file that can be unzipped and has the same data as the original using the diff command and using ls -l and confirming that the bytes are the same
                subprocess.call(shlex.split(f"bgzip -c {self.vcfConcatDir}master_file.vcf"), stdout=f1)
        else:
            with open(f"{self.vcfConcatDir}master_file{self.file_num}.vcf.gz", 'w') as f2:
                subprocess.call(shlex.split(f"bgzip -c {self.vcfConcatDir}master_file.vcf"), stdout=f2)
        print(f"master_file.vcf compression complete at {self.current_time}")

        #Below code can easily just be written into its own function if ever needed
        print(f"master_file.vcf indexing staretd at {self.current_time}")
        if self.file_num == 0:
            subprocess.run(shlex.split(f"tabix -p vcf {self.vcfConcatDir}master_file.vcf.gz"))
        else:
            subprocess.run(shlex.split(f"tabix -p vcf {self.vcfConcatDir}master_file{self.file_num}.vcf.gz"))
        print(f"master_file.vcf indexing complete at {self.current_time}")

    def generate_stats(self):
        print('Generating stats for master_file.vcf...')
        command1 = shlex.split(f'vcftools --gzvcf {self.zipped_master_file} --freq2 --max-alleles 2 --out {self.stats_dir}/{args.prefix[0]} ') # allele freq
        command2 = shlex.split(f'vcftools --gzvcf {self.zipped_master_file} --depth --out {self.stats_dir}/{args.prefix[0]}') # depth per sample 
        command3 = shlex.split(f'vcftools --gzvcf {self.zipped_master_file} --site-mean-depth --out {self.stats_dir}/{args.prefix[0]}') # mean depth per variant site
        command4 = shlex.split(f'vcftools --gzvcf {self.zipped_master_file} --site-depth --out {self.stats_dir}/{args.prefix[0]}') # depth at every variant site
        command5 = shlex.split(f'vcftools --gzvcf {self.zipped_master_file} --site-quality --out {self.stats_dir}/{args.prefix[0]}') # quality per variant site
        command6 = shlex.split(f'vcftools --gzvcf {self.zipped_master_file} --missing-indv --out {self.stats_dir}/{args.prefix[0]}') # fraction of missing data by sample
        command7 = shlex.split(f'vcftools --gzvcf {self.zipped_master_file} --missing-site --out {self.stats_dir}/{args.prefix[0]}') # fraction of missing data per variant site
        command8 = shlex.split(f'vcftools --gzvcf {self.zipped_master_file} --het --out {self.stats_dir}/{args.prefix[0]}') # het & inbreeding coefficient
        commands = [command1, command2, command3, command4, command5, command6, command7, command8]

        procs = [subprocess.Popen(i) for i in commands]
        for p in procs:
            print(f"starting process at {self.current_time}")
            p.wait()
        print(f"Stat calculation for vcf file complete at {self.current_time}")

    def filter_variants(self):
        print('RUNNING VARIANT FILTERING')
        self.filtered_file = self.fm_obj.localOutputDir + 'vcf_concat_output/filtered_master_file.vcf.gz'
        self.pass_file = self.fm_obj.localOutputDir + 'vcf_concat_output/pass_variants_master_file.vcf.gz'
        subprocess.run(shlex.split(f"gatk VariantFiltration \
                                    -R {self.fm_obj.localGenomeFile} \
                                    -V {self.zipped_master_file} \
                                    -O {self.filtered_file} \
                                    --filter-name 'allele_freq' \
                                    --filter-expression 'AF < 0.05' \
                                    --filter-name 'inbreeding_test' \
                                    --filter-expression 'InbreedingCoeff < -0.6' \
                                    --filter-name 'depth_Qual' \
                                    --filter-expression 'QD < 2.0' \
                                    --filter-name 'max_DP' \
                                    --filter-expression 'DP > 10000' \
                                    --filter-name 'min_DP' \
                                    --filter-expression 'DP < 7500' \
                                    --filter-name 'strand_bias' \
                                    --filter-expression 'FS > 40.0' \
                                    --filter-name 'mapping_quality' \
                                    --filter-expression 'MQ < 50.0' \
                                    --filter-name 'no_calls' \
                                    --filter-expression 'NCC > 113' \
                                    --verbosity ERROR"))
        print('FILTERING COMPLETE... EXTRACTING ONLY PASS VARIANTS')
        subprocess.run(['gatk', 'SelectVariants', '-V', self.filtered_file, '--exclude-filtered', '-O', self.pass_file])
        print('INDEXING PASS VARIANTS FILE...')
        subprocess.run(shlex.split(f"tabix -p vcf self.pass_file"))
        print('PASS VARIANT FILE GENERATED')

    def run_methods(self):
        self._make_file_structure()
        if args.merge:
            self.merge_vcfs()
        if args.compress_and_index:
            self.compress_and_index_vcf()
        if args.vcf_stats:
            self.generate_stats()
        if args.filter:
            self.filter_variants()

if __name__ == "__main__":
    vcf_processor_obj = VCFProcessor(args.genome)
    vcf_processor_obj.run_methods()
    print('Pipeline Run Successful')

    


"""
NOTE:
Find filtering stats notes in the "VCF Filtering" notebook on Dropbox

time python process_vcf.py Mzebra_GT3 --filter 2> error_filter_gt3_master_file_240712.txt 1> log_filter_gt3_master_file_240712.txt

time python process_vcf.py Mzebra_UMD2a --merge --compress_and_index --vcf_stats --prefix 419_cohort 2> vcf_pipeline_logs/error_419cohort_concat_24.01.19.txt 1> vcf_pipeline_logs/log_419cohort_concat_24.01.19.txt

time python process_vcf.py Mzebra_GT3 --merge --compress_and_index


time python process_vcf.py Mzebra_GT3 --filter 2> error_filter_gt3_master_file_240723.txt 1> log_filter_gt3_master_file_240723.txt


"""