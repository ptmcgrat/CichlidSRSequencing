import argparse, pdb, os, subprocess, pathlib, pysam, shlex
import pandas as pd
from cyvcf2 import VCF
from helper_modules.nikesh_file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage='This pipeline will take in a set of unaligned bam files generated from Illumina sequencing reads, and call variants using GATK HaplotypeCaller')
parser.add_argument('reference_genome', help = 'full file path to the reference genome that reads will be aligned to for varaint calling')
parser.add_argument('-p', '--projectIDs', help = 'list of projectIDs on which to run the pipeline', nargs = '*', default = ['All'])
parser.add_argument('-i', '--import_databses', help = 'Call this flag to run GenomicsDBImport on the samples in the database', action = 'store_true')
parser.add_argument('-g', '--genotype', help = 'Call this flag to run GenotypeGVCFs on the samples in the database', action = 'store_true')
parser.add_argument('-d', '--download_data', help = 'Use this flag if you need to Download Data from Dropbox to include in the analysis', action = 'store_true')
args = parser.parse_args()

"""
in order to run variant calling, the script will run 2 main commands
gatk GenomicsDBImport
gatk --java-options '-Xmx450G' GenomicsDBImport --genomicsdb-workspace-path {path_to_database_dir + contig_database} --intervals {list_of_4_intervals_for_that_contig} --sample-name-map {path_to_sample_map} --interval-merging-rule OVERLAPPING_ONLY --reader-threads {not_sure_if_this_applies_due_to_next_flag} --max-num-intervals-to-import-in-parallel 4 [optional_during_testing:--overwrite-existing-genomicsdb-workspace]
gatk GenotypeGVCFs
gatk --java-options '-Xmx450G' GenotypeGVCFs -R {reference_genome_path} -V {'gendb://../../../../../Data/mcgrath-lab/Data/CichlidSequencingData/Databases/' + lg7 + '_database/'}  -O {'/Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/' + lg7 + '_output.vcf'} --heterozygosity 0.0012
"""

class VariantCaller:
    def __init__(self, genome, project_ids):
        self.genome = genome
        self.fm_obj = FM("Mzebra_UMD2a")
        self.projectIDs = project_ids

        # Code block to define the set of SampleIDs based on the passed ProjectIDs
        if self.projectIDs == ['All']:
            self.projectIDs = ['ReferenceImprovement', 'RockSand_v1', 'PRJEB15289', 'PRJEB1254','BrainDiversity_s1', 'BigBrain']
        self.alignment_df = pd.read_csv(self.fm_obj.localAlignmentFile)
        filtered_df = self.alignment_df[self.alignment_df['ProjectID'].isin(self.projectIDs)]
        self.sampleIDs = filtered_df['SampleID'].tolist()

        # Get contig names which will be passed into the GATK commands
        self.fasta_obj = pysam.FastaFile(self.fm_obj.localGenomeFile)
        self.contigs = self.fasta_obj.references[0:22] # defines the LG names for the first 22 LGs in the genome.

    def _generate_sample_map(self):
        print('hi from _generate_sample_map')
        sampleIDs = self.sampleIDs
        with open('sample_map.txt', 'w') as fh:
            for sampleID in sampleIDs:
                self.fm_obj.createSampleFiles(sampleID)
                fh.write(sampleID + '\t' + self.fm_obj.localGVCFFile + '\n')

    def data_downloader(self):
        self.fm_obj.downloadData(self.fm_obj.localAlignmentFile)
        # due to memory constructions locally, instead of downloading the full reference genome dict, I only downloaded the localGenomeFile (GCF_000238955.4_M_zebra_UMD2a_genomic.fna). Replace with the following code for server:
        # self.fm_obj.downloadData(self.fm_obj.localGenomeDir)
        self.fm_obj.downloadData(self.fm_obj.localGenomeFile)

        # Download the GVCF and GVCF.idx files for each sample in the AlignmentDatabase
        pdb.set_trace()
        for sampleID in self.sampleIDs:
            if sampleID in []:
                continue
            print(sampleID)
            self.fm_obj.createSampleFiles(sampleID)
            self.fm_obj.downloadData(self.fm_obj.localGVCFFile)
            self.fm_obj.downloadData(self.fm_obj.localGVCFFile + '.tbi')

    def RunGenomicsDBImport(self):
        processes = []
        pdb.set_trace()
        #### First gatk command takes in chromosome names and a tab delimited cohort of samples for which to generate a genomicsdb workspace. The location of the workspace, per chromosome, must be specified using an absolute filepath.
        #### The loop is parallelized to run each chromosome in parallel on 4 cores.
        #### Below is the new code to use after splitting the contigs into intervals and running by importing 4 intervals at a time. 100kbp took about 2.41 mins
        for contig in self.contigs:
            # gatk GenomicsDBImport --genomicsdb-workspace-path {'/Data/mcgrath-lab/Data/CichlidSequencingData/TestingDatabases/' + contig + '_database'} --intervals small_contig.interval_list --sample-name-map sample_map_utaka.txt --max-num-intervals-to-import-in-parallel 4"])
            p = subprocess.Popen(['gatk', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + contig + '_database', '--intervals', os.getcwd() + '/all_lg_intervals/' + contig + '.interval_list', '--sample-map', os.getcwd() + '/sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4'])
            processes.append(p)

        if len(processes) == 22:
            for proc in processes:
                proc.communicate()
            processes = []

    def RunGenotypeGVCFs(self):
        processes = []
        pdb.set_trace()
        #### First gatk command takes in chromosome names and a tab delimited cohort of samples for which to generate a genomicsdb workspace. The location of the workspace, per chromosome, must be specified using an absolute filepath.
        #### The loop is parallelized to run each chromosome in parallel on 4 cores.
        #### Below is the new code to use after splitting the contigs into intervals and running by importing 4 intervals at a time. 100kbp took about 2.41 mins
        for contig in self.contigs:
            # "gatk GenotypeGVCFs -R /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {'gendb://../../../../../Data/mcgrath-lab/Data/CichlidSequencingData/TestingDatabases/' + contig + '_database/'}  -O {'/Data/mcgrath-lab/Data/CichlidSequencingData/TestingOutputs/' + contig + '_output.vcf'} --heterozygosity 0.0012"))
            p = subprocess.Popen(['gatk', 'GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../..' + self.fm_obj.localDatabasesDir + contig + '_database/', '-O', self.fm_obj.localOutputDir + contig + 'output.vcf', '--heterozygosity', '0.0012'])
            processes.append(p)

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

variant_caller_obj = VariantCaller(args.reference_genome, args.projectIDs)
variant_caller_obj.run_methods()

"""
LOCAL TESTING COMMAND TO DOWNLOAD DATA
/Users/kmnike/anaconda3/envs/variant/bin/python3 call_variants.py /Users/kmnike/Data_backup/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz -d

LOCAL TESTING COMMAND TO RUN PIPELINE ON BIGBRAIN AND BRAINDIVERSITY SAMPLES
/Users/kmnike/anaconda3/envs/variant/bin/python3 call_variants.py /Users/kmnike/Data_backup/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz -p BrainDiversity_s1 BigBrain

"""