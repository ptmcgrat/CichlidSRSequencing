import argparse, pdb, glob, shlex, os
from helper_modules.nikesh_file_manager import FileManager as FM
import subprocess

parser = argparse.ArgumentParser(usage = 'Once VCF files per chromopsome have been generated, this script will combine the VCFs per LG into one file then optionally perfform compression, filtering, etc.')
parser.add_argument('genome', help = 'name of reference genome used in the creation of the VCF files', choices = ['Mzebra_UMD2a', 'Mzebra_GT1'], nargs = 1)
parser.add_argument('-m', '--merge', help = 'this flag will merge the VCF outputs into a master file', action = 'store_true')
parser.add_argument('-f', '--filter', help = 'this flag will invoke filtering of the merged vcf file', action = 'store_true')
parser.add_argument('-o', '--output', help = 'name of output VCF file. This will be used as the prefix for all other outputs from the pipeline', default = 'master_file.vcf')
parser.add_argument('-c', '--compress', help = 'callt his flag to compress and index the output vcf file using bgzip, then tabix', action = 'store_true')
args = parser.parse_args()

"""
TODO: 
- implement gatk MergeVcfs
- implement gatk VarinatFiltration
- implement bgzip and tabx 
- implement rclone uploads at the end of every piece of the pipeline, but use a flag to dictate whether the upload will happen
    - be sure to implement chunker uplaods for massive output files. 

"""

class VCFProcessor:
    def __init__(self, output_file_name):
        self.fm_obj = FM(args.genome[0])
        self.output = output_file_name
        self.master_file = self.fm_obj.localOutputDir + self.output + '.gz'

    def _generate_merge_file_list(self):
        file_list = glob.glob(self.fm_obj.localOutputDir + '*_output.vcf')
        print(file_list)
        with open(self.fm_obj.localOutputDir + 'files_to_merge.txt', 'w') as f:
            for file in file_list:
                f.write(file + '\n')

    def merge_vcfs(self):
        file_list = glob.glob(self.fm_obj.localOutputDir + '*_output.vcf')
        with open('files_to_merge.txt', 'w') as f:
            for file in file_list:
                f.write(file + '\n')
        self.files_to_merge = os.getcwd() + '/files_to_merge.txt'
        if os.path.isfile(self.fm_obj.localOutputDir + self.output):
            raise Exception(self.output + ' already exists. Please chose a different output name to avoid restarting a merge which may take hours to complete.')
        else:
            subprocess.run(['gatk', 'MergeVcfs', '-I', self.files_to_merge, '-O', self.master_file])
            subprocess.run(['rm', self.files_to_merge])

    def filter_variants(self):
        print('RUNNING VARIANT FILTERING')
        self.filtered_file = self.fm_obj.localOutputDir + 'filtered_'+ self.output
        self.pass_file = self.fm_obj.localOutputDir + 'pass_variants_' + self.output
        subprocess.run(shlex.split(f"gatk VariantFiltration \
                                    -R {self.fm_obj.localGenomeFile} \
                                    -V {self.master_file} \
                                    -O {self.filtered_file} \
                                    --filter-name 'allele_freq' \
                                    --filter-expression 'AF < 0.000958' \
                                    --filter-name 'inbreeding_test' \
                                    --filter-expression 'InbreedingCoeff < -0.6' \
                                    --filter-name 'depth_Qual' \
                                    --filter-expression 'QD < 2.0' \
                                    --filter-name 'max_DP' \
                                    --filter-expression 'DP > 12000' \
                                    --filter-name 'min_DP' \
                                    --filter-expression 'DP < 9000' \
                                    --filter-name 'strand_bias' \
                                    --filter-expression 'FS > 40.0' \
                                    --filter-name 'mapping_quality' \
                                    --filter-expression 'MQ < 50.0' \
                                    --filter-name 'no_calls' \
                                    --filter-expression 'NCC > 125.0' \
                                    --verbosity ERROR"))
        print('FILTERING COMPLETE... EXTRACTING ONLY PASS VARIANTS')
        subprocess.run(['gatk', 'SelectVariants', '-V', self.filtered_file, '--exclude-filtered', '-O', self.pass_file])
        print('PASS VARIANT FILE GENERATED')

    def compress_vcf(self):
        print('COMPRESSING MASTER VCF FILE...')
        subprocess.run(f'bgzip -c {self.master_file} > {self.master_file}.gz', shell=True)
        print('MASTER VCF FILE COMPRESSION COMPLETE')
        
        print('COMPRESSING FILTERED VCF FILE...')
        subprocess.run(f'bgzip -c {self.pass_file} > {self.pass_file}.gz', shell=True)
        print('FILTERED VCF FILE COMPRESSION COMPLETE')

        print('INDEXING MASTER VCF FILE')
        subprocess.run(f'tabix -p vcf {self.master_file}.gz', shell=True)
        print('MASTER VCF FILE INDEXING COMPLETE')

        print('INDEXING FILTERED VCF FILE')
        subprocess.run(f'tabix -p vcf {self.pass_file}.gz', shell=True)
        print('FILTERED VCF FILE INDEXING COMPLETE')

    def run_methods(self):
        self._generate_merge_file_list()
        if args.merge:
            self.merge_vcfs()
        if args.filter:
            self.filter_variants()
        if args.compress:
            self.compress_vcf()

vcf_processor_obj = VCFProcessor(args.output)
vcf_processor_obj.run_methods()

"""
LOCAL TESTING COMMAND:
/Users/kmnike/anaconda3/envs/pipeline/bin/python3 process_vcf.py Mzebra_UMD2a --filter


"""
