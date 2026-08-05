import argparse, pdb, os, subprocess, datetime, time, multiprocessing, pathlib, random
import pandas as pd
from helper_modules.nikesh_file_manager import FileManager as FM

parser= argparse.ArgumentParser(usage = "This script will pipeup variants from high coverage parents and assign genotypes for low coverage samples")
parser.add_argument('reference_genome', help = 'full file path to the reference genome that reads will be aligned to for varaint calling')
parser.add_argument('bed_file', help = 'name of BED file used for the analysis. Please enter just the name of the file. File paths are not permitted.', type=str)
parser.add_argument('-p', '--pileup', help = 'This flag will pileup variants for all samples in the analysis cohort', action='store_true')
parser.add_argument('-s', '--sampleIDs', help = 'Use this flag to customize the sampleIDs that the pipeline will run for', default = ['All'], choices = ['custom', 'v2_column'], nargs = 1)
parser.add_argument('-b', '--download_bams', help = 'Download the BAM files from the cloud on which to run this pipeline', action = 'store_true')
parser.add_argument('--concurrent_processes', help = 'specify the number of processes to start concurrently', type = int, default = 96)
parser.add_argument('--local_test', help = 'when this flag is called, variables will be preset to test the code locally', action = 'store_true')
args = parser.parse_args()

"""
bcftools mpileup MCYHF1-002/MCYHF1-002.all.bam YH_024/YH_024.all.bam -f Mconophoros_GT1/anchored_kocher_E_Mchenga_conof_Male_contigs_hs_with_kocher_MC_female_molecules_mito_corrected.fasta -R YH_MC_F1s.bed -a AD -Q 0 | more
-R is the bed file to genotype on
-a tells it to output allele depth
-Q is to remove the base quality filter
general command:
bcftools mpileup <bam1> <bam2> etc

How it works:
    1. script takes in a path to a reference genome, *.all.bam files for all samples that you want to include in the analysis, and a bed file that says where we want to genotype 
        a. The -a and -Q flags are also passed to output allele depth, and to remove the base quality filter, respetively 
    2. The principle here is that GATK performs our initial SNP discovery using the parents. However, there seems to be a 0/0 bias where genotypes per variant are called as 0/0 if not enough depth exists for that sample to genotype it at that locus 
    3. bcftools mpileup uses a different algo to perform variant calling. It also allows us to use bam data for low coverage samples and genotype at the locations we deem as informative from the parent genotyping.
        a. This is done using a bed file

How to code:
    1. Figure out a way to encode how samples for joint genotyping are chosen 
        a. Recycle code from callVariants.py 
    2. Be sure different BED fiels can be added.
        a. Maybe add a directory in Dropbox for genotyping BED files that can be chosen like RunInfoFiles for the downloadShortReadData.py script
    3. Direct Outputs somewhere efficiently
    4. Include code to check for and download paths to BAM files 
"""


"""
For Running the code:
time python PileupVariants.py Mconophoros_GT1 YH_MC_F1s.bed --local_test --pileup
real    14m39.037s
user    14m20.219s
sys     0m11.744s


time python PileupVariants.py Mconophoros_GT1 YH_MC_F1s.bed --pileup -s v2_column 2> error_pileup.txt 1> log_pileup.txt
"""

class pileupVariants():
    def __init__(self, genome, bed_file, sampleIDs, processes):
        self.genome = genome
        self.bed_file = bed_file
        self.fm_obj = FM(self.genome)
        self.sampleIDs = sampleIDs
        self.concurrent_processes = processes
        self.current_time = datetime.datetime.now()
        
        # grab a prefix for the vcf file that matches the bed file.
        self.localBedFile = self.fm_obj.localBEDDir + args.bed_file
        self.out_prefix = self.localBedFile.split('/')[-1].split('.')[0]
        self.out_file_path = self.fm_obj.localOutputDir + 'vcf_concat_output/'+ self.out_prefix + '_whole_cohort_variant_pileup.vcf.gz'

        # block for defining sampleIDs if you want specific samples from a SampleDatabase column, using the alignmnetdatabase, or anything else
        if self.sampleIDs == ['custom']:
            self.fm_obj.downloadData(self.fm_obj.localSampleFile_v2) # downloads the most up-to-date SampleDatabase_v2.xlsx file
            s_df = pd.read_excel(self.fm_obj.localSampleFile_v2, sheet_name='SampleLevel') # reads in the SampleLevel sheet from SampleDatabase_v2.xlsx
            custom_samples = input('Enter the custom samples you want to run as a comma separated list without quotation marks or spaces: ')
            self.sampleIDs = custom_samples.split(',')
            # error checking for sampleIDs in custom file
            for sample in self.sampleIDs:
                if not s_df['SampleID'].eq(sample).any():
                    raise Exception(f"{sample} not found in the Sample Database")
        elif self.sampleIDs == ['v2_column']:
            self.fm_obj.downloadData(self.fm_obj.localSampleFile_v2)
            s_df = pd.read_excel(self.fm_obj.localSampleFile_v2, sheet_name='SampleLevel')
            while True:
                column = input('Enter name of column you want to use from SampleDatabase_v2.xlsx: ')
                if column not in s_df.columns.to_list():
                    print('Ivalid Column. Check spelling and re-enter a column name')
                    continue
                else:
                    break
            self.sampleIDs = s_df[s_df[column] == 'Yes'].SampleID.to_list() # basically allows filtering based on columns that use a Yes or No classification per sample. 
        # pre-defining samples for local testing. Pass in the first 3 LGs only since the interval file has been created for only these.
        if args.local_test:
            self.sampleIDs = ['YH_008_m', 'MC-5G11G', 'MCYHBC1-615-1', 'MCYHBC1-677-1']

        print(f"Number of samples for this pipeline run is {len(self.sampleIDs)}")
    
    def _download_BED_file(self):
        print(f"STARTING BED file download at {self.current_time}")
        self.fm_obj.downloadData(self.localBedFile)
        print(f"BED file download COMPLETED at {self.current_time}")
    
    def _write_Sample_File(self):
        self.samples_file = 'pileup_samples.txt'
        print(f"Samples inlcuded in this analysis are {(' ').join(self.sampleIDs)}")
        print(f"Total number of samples is {len(self.sampleIDs)}")
        with open(self.samples_file, 'w') as f:
            for sample in self.sampleIDs:
                self.fm_obj.createSampleFiles(sample)
                f.write(self.fm_obj.localBamFile)
                f.write('\n')

    def DownloadBams(self, sampleID):
            self.fm_obj.createSampleFiles(sampleID)
            if  not pathlib.Path(self.fm_obj.localBamFile).exists():
                print(f"No Bam File found for {sampleID}... Starting Download...")
                print(f"STARTING BAM file and index download for {sampleID} at {self.current_time}")
                self.fm_obj.downloadData(self.fm_obj.localBamFile)
                self.fm_obj.downloadData(self.fm_obj.localBamIndex)
                print(f"BAM file and index for {sampleID} COMPLETED at {self.current_time}")
            else:
                print(f"All samples for this analysis have BAM files already downloaded.")

    def PileupVariants(self):
        # command to perform variant pileup of LC fish alongside high coverage parents
        print(f"STARTING variant pileup for {len(self.sampleIDs)} samples at {self.current_time}")
        subprocess.run(['bcftools', 'mpileup', '-f', self.fm_obj.localGenomeFile, '-R', self.localBedFile, '-b', self.samples_file, '-a', 'AD', '-Q', '0', '-o', self.out_file_path, '-Oz', '--threads', str(self.concurrent_processes)])
        print(f"Variant pileup for {len(self.sampleIDs)} samples COMPLETED at {self.current_time}")

        # Remove the sample file needed to pass samples to bcftools 
        subprocess.run(['rm', 'pileup_sample.txt'])

        """
        command line testing 
        bcftools mpileup -f /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mconophoros_GT1/anchored_kocher_E_Mchenga_conof_Male_contigs_hs_with_kocher_MC_female_molecules_mito_corrected.fasta -R /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/BEDFiles/YH_MC_F1s.bed -b pileup_samples.txt -a AD -Q 0 -o test_pileup.vcf.gz -Oz --threads 16
        """

    def multiprocess(self, function, sample_type):
            # Author: Lauren Sabo; edits made by NK
            if sample_type == 'sampleID':
                inputs = self.sampleIDs
            concurrent_processes = min(self.concurrent_processes, len(inputs))
            try:
                with multiprocessing.Pool(processes=concurrent_processes) as pool:
                    pool.map(function, inputs)
            except Exception as e:
                print(f"Error occurred during multiprocessing: {e}")

    def run_methods(self):
        self._download_BED_file()
        self._write_Sample_File()
        if args.download_bams:
            self.multiprocess(self.DownloadBams, 'sampleID')
        if args.pileup:
            self.PileupVariants()
        

if __name__ == "__main__":
    variant_caller_obj = pileupVariants(args.reference_genome, args.bed_file, args.sampleIDs, args.concurrent_processes)
    variant_caller_obj.run_methods()
    print('PIPELINE RUN COMPLETE')