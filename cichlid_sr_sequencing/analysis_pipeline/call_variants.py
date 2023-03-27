import argparse, pdb, os, subprocess, pathlib, pysam, shlex
import pandas as pd
from helper_modules.nikesh_file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage='This pipeline will take in a set of unaligned bam files generated from Illumina sequencing reads, and call variants using GATK HaplotypeCaller')
parser.add_argument('reference_genome', help = 'full file path to the reference genome that reads will be aligned to for varaint calling')
parser.add_argument('-p', '--projectIDs', help = 'list of projectIDs on which to run the pipeline', nargs = '*', default = ['All'])
parser.add_argument('-i', '--import_databses', help = 'Call this flag to run GenomicsDBImport on the samples in the database', action = 'store_true')
parser.add_argument('-g', '--genotype', help = 'Call this flag to run GenotypeGVCFs on the samples in the database', action = 'store_true')
parser.add_argument('-d', '--download_data', help = 'Use this flag if you need to Download Data from Dropbox to include in the analysis', action = 'store_true')
parser.add_argument('--local_test', help = 'when this flag is called, variables will be preset to test the code locally', action = 'store_true')
args = parser.parse_args()

"""
in order to run variant calling, the script will run 2 main commands
gatk GenomicsDBImport
# gatk GenomicsDBImport --genomicsdb-workspace-path {'/Data/mcgrath-lab/Data/CichlidSequencingData/TestingDatabases/' + contig + '_database'} --intervals small_contig.interval_list --sample-name-map sample_map_utaka.txt --max-num-intervals-to-import-in-parallel 4"])
gatk GenotypeGVCFs
# "gatk GenotypeGVCFs -R /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {'gendb://../../../../../Data/mcgrath-lab/Data/CichlidSequencingData/TestingDatabases/' + contig + '_database/'}  -O {'/Data/mcgrath-lab/Data/CichlidSequencingData/TestingOutputs/' + contig + '_output.vcf'} --heterozygosity 0.0012"))
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

        # redefining samples and contigs for local testing
        if args.local_test:
            self.sampleIDs = ['MC_1_m', 'SAMEA2661294', 'SAMEA2661322', 'SAMEA4032100', 'SAMEA4033261']
            self.contigs = ['NC_036780.1', 'NC_036781.1', 'NC_036782.1']

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
        for contig in self.contigs:
            
            if args.local_test:
                p = subprocess.Popen(['gatk', '--java-options', '-Xmx450G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + contig + '_database', '--intervals', os.getcwd() + '/all_lg_intervals/test_intervals/' + contig + '.interval_list', '--sample-name-map', os.getcwd() + '/local_test_sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4', '--overwrite-existing-genomicsdb-workspace'])
                processes.append(p)
            else:
                p = subprocess.Popen(['gatk', '--java-options', '-Xmx450G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + contig + '_database', '--intervals', os.getcwd() + '/all_lg_intervals/' + contig + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4', '--overwrite-existing-genomicsdb-workspace'])
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
        for contig in self.contigs:
            if args.local_test:
                p = subprocess.Popen(['gatk', '--java-options', '-Xmx450G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../' + self.fm_obj.localDatabasesDir + contig + '_database/', '-O', self.fm_obj.localOutputDir + contig + 'output.vcf', '--heterozygosity', '0.0012'])
                processes.append(p)
            else:
                p = subprocess.Popen(['gatk', '--java-options', '-Xmx450G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../../' + self.fm_obj.localDatabasesDir + contig + '_database/', '-O', self.fm_obj.localOutputDir + contig + 'output.vcf', '--heterozygosity', '0.0012'])
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

variant_caller_obj = VariantCaller(args.reference_genome, args.projectIDs)
variant_caller_obj.run_methods()

"""
LOCAL TESTING COMMAND TO DOWNLOAD DATA
/Users/kmnike/anaconda3/envs/variant/bin/python3 call_variants.py /Users/kmnike/Data_backup/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz -d

LOCAL TESTING COMMAND TO RUN PIPELINE ON BIGBRAIN AND BRAINDIVERSITY SAMPLES
/Users/kmnike/anaconda3/envs/variant/bin/python3 call_variants.py /Users/kmnike/Data_backup/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz --local_test --import_databses --genotype

RUNNING WHOLE PIPELINE ON UTAKA SERVER, DOWNLOADING ALL NEEDED DATA, AND RUNNING EACH GATK COMMAND IN PARALLEL:
python3 call_variants.py /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz -p BrainDiversity_s1 BigBrain --download_data --import_databses --genotype

"""