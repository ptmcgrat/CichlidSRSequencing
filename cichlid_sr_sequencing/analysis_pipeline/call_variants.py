import argparse, pdb, os, subprocess
import pandas as pd
from helper_modules.nikesh_file_manager import FileManager as FM
from cyvcf2 import VCF

parser = argparse.ArgumentParser(usage='This pipeline will take in a set of unaligned bam files generated from Illumina sequencing reads, and call variants using GATK HaplotypeCaller')
parser.add_argument('reference_genome', help = 'full file path to the reference genome that reads will be aligned to for varaint calling')
parser.add_argument('-p', '--projectIDs', help = 'list of projectIDs on which to run the pipeline', nargs = '*', default = ['All'])
parser.add_argument('-i', '--import_databses', help = 'Call this flag to run GenomicsDBImport on the samples in the database', action = 'store_true')
parser.add_argument('-g', '--genotype', help = 'Call this flag to run GenotypeGVCFs on the samples in the database', action = 'store_true')
parser.add_argument('-d', '--download_data', help = 'Use this flag if you need to Download Data from Dropbox to include in the analysis', action = 'store_true')
parser.add_argument('-r', '--regions', help = 'list of linkage groups for which analyses will run', nargs = '*', default = ['All'])
parser.add_argument('-b', '--bam-download', help = 'Download the BAM files from the cloud on which to call HaplotypeCaller', action = 'store_true')
parser.add_argument('-h', '--halplotype-caller', help = 'run the gatk HaplotypeCaller algorithm to re-generate GVCF files on which to call the pipeline', action = 'store_true')
parser.add_argument('--local_test', help = 'when this flag is called, variables will be preset to test the code locally', action = 'store_true')
args = parser.parse_args()

"""
in order to run variant calling, the script will run 2 main commands
gatk GenomicsDBImport
# gatk GenomicsDBImport --genomicsdb-workspace-path {'/Data/mcgrath-lab/Data/CichlidSequencingData/TestingDatabases/' + lg + '_database'} --intervals small_contig.interval_list --sample-name-map sample_map_utaka.txt --max-num-intervals-to-import-in-parallel 4"])
gatk GenotypeGVCFs
# "gatk GenotypeGVCFs -R /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {'gendb://../../../../../Data/mcgrath-lab/Data/CichlidSequencingData/TestingDatabases/' + lg + '_database/'}  -O {'/Data/mcgrath-lab/Data/CichlidSequencingData/TestingOutputs/' + lg + '_output.vcf'} --heterozygosity 0.0012"))
"""

class VariantCaller:
    def __init__(self, genome, project_ids, linkage_groups):
        self.genome = genome
        self.fm_obj = FM("Mzebra_UMD2a")
        self.projectIDs = project_ids

        # Code block to define the set of SampleIDs based on the passed ProjectIDs
        if self.projectIDs == ['All']:
            self.projectIDs = ['ReferenceImprovement', 'RockSand_v1', 'PRJEB15289', 'PRJEB1254','BrainDiversity_s1', 'BigBrain']
        self.alignment_df = pd.read_csv(self.fm_obj.localAlignmentFile)
        filtered_df = self.alignment_df[self.alignment_df['ProjectID'].isin(self.projectIDs)]
        self.sampleIDs = filtered_df['SampleID'].tolist()

        # # Get contig names which will be passed into the GATK commands
        # self.fasta_obj = pysam.FastaFile(self.fm_obj.localGenomeFile)
        # self.contigs = self.fasta_obj.references[0:22] # defines the LG names for the first 22 LGs in the genome.

        # Code block for determining which linkage groups will be processed by the script:
        self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1', 
                             'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1', 
                             'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1', 
                             'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1'}
        self.linkage_groups = linkage_groups
        if self.linkage_groups == ['All']:
            self.linkage_groups = self.linkage_group_map.values()
        else:
            regions_list = []
            for region in self.linkage_groups:
                if region in self.linkage_group_map.keys():
                    regions_list.append(self.linkage_group_map[region])
                else:
                    raise Exception(region + ' is not a valid option')
            self.linkage_groups = regions_list
            duplicate_set_test = set(self.linkage_groups)
            if len(duplicate_set_test) != len(self.linkage_groups):
                raise Exception('A repeat region has been provided')
            
        # pre-defining samples for local testing. Pass in the first 3 LGs only since the interval file has been created for only these.
        if args.local_test:
            self.sampleIDs = ['MC_1_m', 'SAMEA2661294', 'SAMEA2661322', 'SAMEA4032100', 'SAMEA4033261']

    def _generate_sample_map(self):
        if args.local_test:
            pass
        else:
            sampleIDs = self.sampleIDs
            with open('sample_map.txt', 'w') as fh:
                for sampleID in sampleIDs:
                    self.fm_obj.createSampleFiles(sampleID)
                    fh.write(sampleID + '\t' + self.fm_obj.localGVCFFile + '\n')

    def data_downloader(self):
        print('Downloading new Alignment File')
        self.fm_obj.downloadData(self.fm_obj.localAlignmentFile)
        print('Downloading most recent Genome Dir')
        self.fm_obj.downloadData(self.fm_obj.localGenomeDir)

        # Download the GVCF and GVCF.idx files for each sample in the AlignmentDatabase
        for sampleID in self.sampleIDs:
            if sampleID in []:
                continue
            print('Downloading ' + sampleID + '...')
            self.fm_obj.createSampleFiles(sampleID)
            self.fm_obj.downloadData(self.fm_obj.localGVCFFile)
            self.fm_obj.downloadData(self.fm_obj.localGVCFFile + '.tbi')
            print('Done Downloading ' + sampleID)

    def RunGenomicsDBImport(self):
        processes = []
        for lg in self.linkage_groups:
            
            if args.local_test:
                p = subprocess.Popen(['gatk', '--java-options', '-Xmx450G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + lg + '_database', '--intervals', os.getcwd() + '/all_lg_intervals/test_intervals/' + lg + '.interval_list', '--sample-name-map', os.getcwd() + '/local_test_sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4', '--overwrite-existing-genomicsdb-workspace'])
                processes.append(p)
            else:
                p = subprocess.Popen(['gatk', '--java-options', '-Xmx450G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + lg + '_database', '--intervals', os.getcwd() + '/all_lg_intervals/' + lg + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4', '--overwrite-existing-genomicsdb-workspace'])
                processes.append(p)
            if args.local_test:
                if len(processes) == 3:
                    for proc in processes:
                        proc.communicate()
                    processes = []
            else:
                if len(processes) == 22:
                    for proc in processes:
                        proc.communicate()
                    processes = []

    def RunGenotypeGVCFs(self):
        processes = []
        for lg in self.linkage_groups:
            if args.local_test:
                p = subprocess.Popen(['gatk', '--java-options', '-Xmx450G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../' + self.fm_obj.localDatabasesDir + lg + '_database/', '-O', self.fm_obj.localOutputDir + lg + 'output.vcf', '--heterozygosity', '0.0012'])
                processes.append(p)
            else:
                p = subprocess.Popen(['gatk', '--java-options', '-Xmx450G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../../' + self.fm_obj.localDatabasesDir + lg + '_database/', '-O', self.fm_obj.localOutputDir + lg + 'output.vcf', '--heterozygosity', '0.0012'])
                processes.append(p)

            if args.local_test:
                if len(processes) == 3:
                    for proc in processes:
                        proc.communicate()
                    processes = []
            else:
                if len(processes) == 22:
                    for proc in processes:
                        proc.communicate()
                    processes = []

    def run_methods(self):
        self._generate_sample_map()
        if args.download_data:
            print('download data will run')
            self.data_downloader()
        if args.import_databses:
            self.RunGenomicsDBImport()
        if args.genotype:
            self.RunGenotypeGVCFs()

variant_caller_obj = VariantCaller(args.reference_genome, args.projectIDs, args.regions)
variant_caller_obj.run_methods()

"""
LOCAL TESTING COMMAND TO DOWNLOAD DATA
/Users/kmnike/anaconda3/envs/variant/bin/python3 call_variants.py /Users/kmnike/Data_backup/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz -d

LOCAL TESTING COMMAND TO RUN PIPELINE ON BIGBRAIN AND BRAINDIVERSITY SAMPLES
/Users/kmnike/anaconda3/envs/variant/bin/python3 call_variants.py /Users/kmnike/Data_backup/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz --local_test --import_databses --genotype --regions LG1 LG2 LG3

RUNNING WHOLE PIPELINE ON UTAKA SERVER, DOWNLOADING ALL NEEDED DATA, AND RUNNING EACH GATK COMMAND IN PARALLEL:
python3 call_variants.py /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz -p BrainDiversity_s1 BigBrain --download_data --import_databses --genotype

"""