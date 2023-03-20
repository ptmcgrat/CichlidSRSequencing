import argparse, pdb, os, subprocess, pathlib
from cyvcf2 import VCF
from helper_modules.nikesh_file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage='This pipeline will take in a set of unaligned bam files generated from Illumina sequencing reads, and call variants using GATK HaplotypeCaller')
parser.add_argument('reference_genome', help = 'full file path to the reference genome that reads will be aligned to for varaint calling')
parser.add_argument('analysis_ID', help = 'Name of Analysis ID for the samples to be processed', type = str, nargs = 1)
parser.add_argument('-o', '--output_dir', help = 'Full file path to an output dir where analyses will be carried out.')
parser.add_argument('-d', '--download_data', help = 'Use this flag if you need to Download Data from Dropbox to include in the analysis', action = 'store_true')
args = parser.parse_args()

"""
in order to run variant calling, the script will run 2 main commands
gatk GenomicsDBImport
gatk --java-options '-Xmx450G' GenomicsDBImport --genomicsdb-workspace-path {path_to_database_dir + contig_database} --intervals {list_of_4_inteervals_for_that_contig} --sample-name-map {path_to_sample_map} --interval-merging-rule OVERLAPPING_ONLY --reader-threads {not_sure_if_this_applies_due_to_next_flag} --max-num-intervals-to-import-in-parallel 4 [optional_during_testing:--overwrite-existing-genomicsdb-workspace]
gatk GenotypeGVCFs
gatk --java-options '-Xmx450G' GenotypeGVCFs -R {reference_genome_path} -V {'gendb://../../../../../Data/mcgrath-lab/Data/CichlidSequencingData/Databases/' + lg7 + '_database/'}  -O {'/Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/' + lg7 + '_output.vcf'} --heterozygosity 0.0012

Additonal files the script needs to create/take in as input:
1. directory of intervals per linkage group (already manually curated)
2. sample_map.tsv file containing the file paths to each sample we want to run for analysis. Here's the format per sample:
    MZ_1_m	/Data/mcgrath-lab/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/MZ_1_m/MZ_1_m.g.vcf.gz

I should incorporate the Download Data into this script... The downloadData method is called but first defining an object of the FileManager class (like fm_obj), then calling fm_obj.downloadData(fm_obj.localAnalysisFile)
The fm_obj.localanalysisfile
since the data generation will follow the struicture set by DownloadData, I should make some assumptions about the sample file storage location since I'll adopt the file manager for this script... 

"""

class VariantCaller:
    def __init__(self, genome, analysis_ID):
        self.genome = genome
        self.analysisID = analysis_ID[0]
        # pdb.set_trace()

    def data_downloader(self, analysis_ID):
        analysis_ID = self.analysisID
        fm_obj = FM(self.genome)
        fm_obj.createAnalysisIDFiles(analysis_ID)
        fm_obj.downloadData(fm_obj.localAnalysisFile)

    def _generate_sample_map(self, sample_dir):
        print('hi from _generate_sample_map')
    
    def run_methods(self):
        if args.download_data:
            print('download data will run')
            self.data_downloader(self.analysisID)


variant_caller_obj = VariantCaller(args.reference_genome, args.analysis_ID)
variant_caller_obj.run_methods()



"""
LOCAL TESTING COMMAND:
/Users/kmnike/anaconda3/envs/variant/bin/python3 call_variants.py /Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz BrainDiversity_s1 -d

"""