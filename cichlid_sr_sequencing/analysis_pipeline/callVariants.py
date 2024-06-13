import argparse, pdb, os, subprocess
import pandas as pd
from multiprocessing import Process
from helper_modules.nikesh_file_manager import FileManager as FM # type: ignore

parser = argparse.ArgumentParser(usage='This pipeline will take in a set of unaligned bam files generated from Illumina sequencing reads, and call variants using GATK HaplotypeCaller')
parser.add_argument('reference_genome', help = 'full file path to the reference genome that reads will be aligned to for varaint calling')
parser.add_argument('-e', '--ecogroups', help = 'one or multiple eco group names that will filter the samples on which the pipeline is run', choices = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'Riverine', 'AC', 'Non_Riverine', 'All', 'Lake_Malawi', 'Rock_Sand', 'Sand'], nargs = '*', default = ['All'])
parser.add_argument('-p', '--projectIDs', help = 'list of projectIDs on which to run the pipeline. Custom IDs can be provided as a csv file containing one sample per line. Save the file as "custom_samples.csv" in the directory containing this script.', nargs = '*', default = ['All'])
parser.add_argument('-i', '--import_databases', help = 'Call this flag to run GenomicsDBImport on the samples in the database', action = 'store_true')
parser.add_argument('-g', '--genotype', help = 'Call this flag to run GenotypeGVCFs on the samples in the database', action = 'store_true')
parser.add_argument('-d', '--download_GVCF_data', help = 'Use this flag if you need to download the GVCF files from Dropbox to include in the VCF analysis', action = 'store_true')
parser.add_argument('-r', '--regions', help = 'list of linkage groups for which analyses will run', nargs = '*', default = ['All'])
parser.add_argument('-b', '--download_bams', help = 'Download the BAM files from the cloud on which to call HaplotypeCaller', action = 'store_true')
parser.add_argument('-H', '--efficient_haplotypecaller', help = 'use this flag to download BAM files and run HaplotypeCaller on samples', action = 'store_true')
parser.add_argument('-l', '--local_test', help = 'when this flag is called, variables will be preset to test the code locally', action = 'store_true')
parser.add_argument('-m', '--memory', help = 'How much memory, in GB, to allocate to each child process', default = [4], nargs = 1)
parser.add_argument('-u', '--unmapped', help = 'Use this flag to run -i and -g on the unmapped contigs in the genome', action = 'store_true')
parser.add_argument('--concurrent_processes', help = 'specify the number of processes to start concurrently', type = int, default = 96)
args = parser.parse_args()

"""
NOTE: 
as of 2024 June 5, I can conda install -c bioconda gatk4 into a fresh conda env. The version installed is gatk4-4.0.5.1-0
The GATK 4.5.0.0 binary I have seems to require a new version of Java... maybe Java runtime 17 according to one post I read here: https://gatk.broadinstitute.org/hc/en-us/community/posts/17535513525915-Annotation-for-variant-calling
I have gatk working in an env called 'gatk' for now.
NOTE:
2024 May 30
We have a GT3 genome! Patrick has been running alignment to the new genome recently.
The new genome has 22 LGs + 40 unmapped contigs + mito + one unplaced scaffold 
Total conitgs = 64

The first 22 contigs are 933Mbp 
The unmapped content is 28Mbp. 
    - The largest of these is 9Mbp
    - 2nd largest is 4Mbp and it drops off from there

Ther goal is to parallelize GenotypeGVCFs
Issue is that GenotypeGVCFs operates on a database
This means each contigs needs to be separated into smaller intervals
This will give smaller VCF files that will then need to be end-to-end concatenated


How GenomicsDBImport parallelization works:
currently for each LG, I have 4 intervals allowing each Lg to be read in 4 sections at a time and written into the same database
To use up as many cores as possible, I can leave one core per unmapped contig? 
Or I can run all the unmapped ones after the fact... lets see....

Ok so to maximize efficiency for GenotypeGVCFs, I need to break the 933Mbp into 96 equivalent sections. Issues may pop up when bridging contigs 
Look at notes in the GT3_intervals/make_GT3_intervals.py script

2024 June 04 
The 96 intervals have been generated to maximize the number of cores that will be used when running everything 
I'll downlaod 3 new BAM files and test this on my mac usign a max of 10 cores (or 8 cores to be safe?)
Note that the --local_test flag will be updated to work with the GT3 genome and maximize cores to generate 10 databases then 10 VCFs which will be concatenated into one VCF file 

TODO:
1. sequential data download, GVCF generation, file upload and data deletion
    - I'll upload the partial BAM files for testing and run it 1 at a time or 2 at a time.
    - DONE
2. work on HaplotypeCaller changes 
    - BAMS should be downloaded I think 96 at a time and processed through haplotypecaller 
    - The GVCFs should be uploaded and the BAMs removed
    - DONE

3. Parallelize GenomicsDBImport 
    - the 96 intervals should be used 
    - paths to the databases will need to be updated in the sample file or something like that for GenotypeGVCFs

4. Parallelize GenotypeGVCFs (should be easy)
"""

class VariantCaller:
    def __init__(self, genome, project_ids, linkage_groups, memory, ecogroups, processes):
        self.genome = genome
        self.fm_obj = FM(self.genome)
        self.projectIDs = project_ids
        self.memory = memory
        self.ecogroups = ecogroups
        self.concurrent_processes = processes
        # TODO: need to be able to read a list of the valid ecogroups in the Ecogroup_PTM column and let all samples be captured if 'All' is called as teh ecogroup. For now, only the 'LakeMalawi' ecogroups are recognized. 
        if self.ecogroups == ['All']:
            self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Shallow_Benthic2', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'Riverine', 'AC'] # Note that Patrick added a 'Shallow_Benthic2' ecogroup for some reason aorund May 2024... Unsure what this is. need to talk to him about it 
        elif self.ecogroups == ['Non_Riverine']:
            self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Shallow_Benthic2', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon']
        elif self.ecogroups == ['Lake_Malawi']:
            self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Shallow_Benthic2', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'AC']
        elif self.ecogroups == ['Rock_Sand']:
            self.ecogroups = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Shallow_Benthic2', 'Deep_Benthic']
        elif self.ecogroups == ['Sand']:
            self.ecogroups = ['Utaka', 'Shallow_Benthic', 'Shallow_Benthic2', 'Deep_Benthic']

        # Code block to set the ecogroups
        self.fm_obj.downloadData(self.fm_obj.localSampleFile_v2) # downloads the most up-to-date SampleDatabase_v2.xlsx file
        s_df = pd.read_excel(self.fm_obj.localSampleFile_v2, sheet_name='SampleLevel') # reads in the SampleLevel sheet from SampleDatabase_v2.xlsx
        eg_filtered_df = s_df[s_df['Ecogroup_PTM'].isin(self.ecogroups)] # filter samples that match the ecogroups in the analysis
        self.sampleIDs = eg_filtered_df['SampleID'].to_list()

        for sample in self.sampleIDs:
            if not s_df['SampleID'].eq(sample).any():
                raise Exception(f"{sample} not found in the Sample Database")

        # Code block for determining which linkage groups will be processed by the script:
        self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1', 
                             'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1', 
                             'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1', 
                             'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1', 'mito': 'NC_027944.1'}
        self.linkage_groups = linkage_groups
        if self.linkage_groups == ['All']:
            self.linkage_groups = list(self.linkage_group_map.values())
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
            self.sampleIDs = ['CJ_2204_m', 'CV-006-m', 'LA_3006_m', 'MC-008-m', 'OC-001-m']
            self.memory = [3]
            self.linkage_groups = ['NC_036780.1', 'NC_036781.1', 'NC_036782.1']
            self.concurrent_processes = 1


    def _generate_sample_map(self):
        sampleIDs = self.sampleIDs
        with open('sample_map.txt', 'w') as fh:
            for sampleID in sampleIDs:
                self.fm_obj.createSampleFiles(sampleID)
                # commenting out below line since BAM and GVCF files are stored at /Output/Bamfiles/Mzebra_UMD2a
                # fh.write(sampleID + '\t' + self.fm_obj.localGVCFFile + '\n')
                if not args.local_test:
                    fh.write(sampleID + '\t' + self.fm_obj.StorageGVCFFile + '\n')
                else:
                    fh.write(sampleID + '\t' + self.fm_obj.localTestGVCFFile + '\n')
        pdb.set_trace()
    def GVCF_downloader(self):
        print('Downloading new Alignment File')
        self.fm_obj.downloadData(self.fm_obj.localAlignmentFile)
        print('Downloading most recent Genome Dir')
        self.fm_obj.downloadData(self.fm_obj.localGenomeDir)

        # Download the GVCF and GVCF.idx files for each sample in the AlignmentDatabase
        for sampleID in self.sampleIDs:
            if sampleID in []: # confused what this is doing... rememdy it...
                continue
            print('Downloading ' + sampleID + '...')
            self.fm_obj.createSampleFiles(sampleID)
            self.fm_obj.downloadData(self.fm_obj.localGVCFFile)
            self.fm_obj.downloadData(self.fm_obj.localGVCFFile + '.tbi')
            print('Done Downloading ' + sampleID)

    def EfficientHaplotypeCaller(self, sample):
        # to use multiprocess, this function will sequenctially need to perform BAM/BAI download, GVCF generation, and then deleting the BAM/BAI files
        if args.local_test:
            self.fm_obj.createSampleFiles(sample)
            print('Downloading BAM and BAM index for ' + sample + '...')
            self.fm_obj.downloadData(self.fm_obj.localTestBamFile)
            self.fm_obj.downloadData(self.fm_obj.localTestBamIndex)

            print('Generating GCVF file for ' + sample + '...')
            subprocess.run(['gatk', 'HaplotypeCaller', '--emit-ref-confidence', 'GVCF', '-R', self.fm_obj.localGenomeFile, '-I', self.fm_obj.localTestBamFile, '-O', self.fm_obj.localTestGVCFFile])
            
            print('Removing BAM file for ' + sample + '...')
            os.remove(self.fm_obj.localTestBamFile)
            os.remove(self.fm_obj.localTestBamIndex)
        else:
            self.fm_obj.createSampleFiles(sample)
            print('Downloading BAM and BAM index for ' + sample + '...')
            self.fm_obj.downloadData(self.fm_obj.localBamFile)
            self.fm_obj.downloadData(self.fm_obj.localBamIndex)
            
            print('Generating GCVF file for ' + sample + '...')
            subprocess.run(['gatk', 'HaplotypeCaller', '--emit-ref-confidence', 'GVCF', '-R', self.fm_obj.localGenomeFile, '-I', self.fm_obj.localBamFile, '-O', self.fm_obj.localGVCFFile])

            print('Removing BAM file for ' + sample + '...')
            os.remove(self.fm_obj.localBamFile)
            os.remove(self.fm_obj.localBamIndex)

    def RunGenomicsDBImport(self, interval):
        # this funciton is what exists in the multiprocvessGATK.py file. It may have some edits not in the other script so u can copy this over to there if u wanna code on the other script instead of this one which has the whole pipeline. 
        print('starting processing ' + interval + ' for all samples in cohort')
        if args.local_test:
            subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + interval + '_database', '--intervals', os.getcwd() + '/GT3_intervals/' + interval + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt','--max-num-intervals-to-import-in-parallel', '1', '--overwrite-existing-genomicsdb-workspace'])
        elif not args.unmapped:
            subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + interval + '_database', '--intervals', os.getcwd() + '/all_lg_intervals/' + interval + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4',])
        else:
            subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + interval + '_database', '--intervals', os.getcwd() + '/unmapped_contig_intervals/' + interval + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4'])
        print('GENOMICSDBIMPORT RUN SUCCESSFULLY FOR interval ' + interval)

    def RunGenotypeGVCFs(self, interval):
        # You can ignore this function for now, but I implemented the parallelization code into it
        # Still need to add code to run the unmapped contigs 
        print('starting processing ' + interval + ' for all samples in cohort')
        if args.local_test: # update the location of the GenDB if "local testing" on Utaka vs on the mac. As of Oct 13, 2023, gatk4 does not work via a conda install on my M2 mac on Ventura 13.4.1 and it also doesn't work on a fresh conda env on Utaka (gatk 4.0.-.- gets installed when I need at least 4.3.0.0). gatk 4.3.0.0 still works on Utaka in the 'genomics' env
            # path is to the gendb located on the Utaka server 
            local_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../../' + self.fm_obj.localDatabasesDir + interval + '_database/', '-O', self.fm_obj.localOutputDir + interval + '_output.vcf', '--heterozygosity', '0.00175'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
            local_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
            local_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
            local_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
            local_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
            subprocess.run(local_command)
            print('GENOTYPEGVCFS RUN SUCCESSFULLY FOR ' + interval)
        elif not args.unmapped: # make sure this will work with the same parameters
            genotypegvcfs_unmapped_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../../' + self.fm_obj.localDatabasesDir + interval + '_database/', '-O', self.fm_obj.localOutputDir + interval + '_output.vcf', '--heterozygosity', '0.00175'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
            genotypegvcfs_unmapped_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
            genotypegvcfs_unmapped_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
            genotypegvcfs_unmapped_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
            genotypegvcfs_unmapped_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
            subprocess.run(genotypegvcfs_unmapped_command)
        else:
            genotypegvcfs_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../../' + self.fm_obj.localDatabasesDir + interval + '_database/', '-O', self.fm_obj.localOutputDir + interval + '_output.vcf', '--heterozygosity', '0.00175'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
            genotypegvcfs_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
            genotypegvcfs_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
            genotypegvcfs_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
            genotypegvcfs_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
            subprocess.run(genotypegvcfs_command)
            print('GENOTYPEGVCFS RUN SUCCESSFULLY FOR ' + interval)

    def multiprocess(self, function, sample_type):
        # Author: Lauren Sabo; edits made by NK
        concurrent_processes = self.concurrent_processes
        # the below code will allow multiprocess to be run on a function based on SampleIDs or by LG. Similar code can be added to breakup the run by either projectID, intervals, etc. 
        if sample_type == 'lg':
            blocks_to_process = [self.linkage_groups[i:int(i+concurrent_processes)] for i in range(0, len(self.linkage_groups), int(concurrent_processes))] # this generates the list of lists to process together. For instance, if processing 3 LGs, 2 at a time,  concurrent_processes is [['NC_036780.1', 'NC_036781.1'], ['NC_036782.1']]
        elif sample_type == 'sampleID':
            blocks_to_process = [self.sampleIDs[i:int(i+concurrent_processes)] for i in range(0, len(self.sampleIDs), int(concurrent_processes))]
        elif sample_type == 'interval':
            intervals = list(range(1,97))
            str_intervals = list(map(str, intervals))
            blocks_to_process = [str_intervals[i:int(i+concurrent_processes)] for i in range(0, len(str_intervals), int(concurrent_processes))]

        jobs = []
        for block in blocks_to_process: # for each sublist of processes to start in the larger list of processes:
            for contig in block:
                j = Process(target = function, args = (contig,)) # define the processes
                jobs.append(j) # append the processes we want to start to the list of jobs 
                j.start()
                print(contig, " - started")
        
            i = 0
            for j in jobs:
                j.join()
                print(i, " - finished")
                i += 1

            del jobs[:]
        pdb.set_trace()

    def run_methods(self):
        self._generate_sample_map()
        if args.download_GVCF_data:
            print('download data will run and GVCF files per sample will be downloaded')
            self.GVCF_downloader()
        if args.efficient_haplotypecaller:
            self.multiprocess(self.EfficientHaplotypeCaller, 'sampleID')
        if args.import_databases:
            self.multiprocess(self.RunGenomicsDBImport, 'interval')
        if args.genotype:
            self.multiprocess(self.RunGenotypeGVCFs, 'interval')


if __name__ == "__main__":
    variant_caller_obj = VariantCaller(args.reference_genome, args.projectIDs, args.regions, args.memory, args.ecogroups, args.concurrent_processes)
    variant_caller_obj.run_methods()
    print('PIPELINE RUN SUCCESSFUL')

"""
LOCAL TESTING COMMAND SKELETON
/Users/kmnike/anaconda3/envs/variant/bin/python3 call_variants.py Mzebra_UMD2a --local_test --regions LG1 LG2 LG3

time python callVariants.py Mzebra_UMD2a -e Lake_Malawi --import_databases --genotype --memory 40 --concurrent_processes 23 2> vcf_pipeline_logs/error_419cohort_24.01.17.txt 1> vcf_pipeline_logs/log_419cohort_24.01.17.txt



RUNNING WHOLE PIPELINE ON UTAKA SERVER, DOWNLOADING ALL NEEDED DATA, AND RUNNING EACH GATK COMMAND IN PARALLEL:
python3 call_variants.py Mzebra_UMD2a -p BrainDiversity_s1 BigBrain --import_databases --genotype -m 21

The pipeline is now ready to run on Utaka server since all of the HaplotypeCaller reruns for the previously failed LG7 samples is complete and since the new BigBrain & BrainDiversity_s1 datasets have been rerun, and teh files have been renamed.
Note that I need to process samples using 4 cores each. Here is the command to run:

python3 call_variants.py Mzebra_UMD2a -p All --import_databases --genotype --regions All --memory 2

For rerunning LG14:
python3 call_variants.py Mzebra_UMD2a -p All --import_databases --genotype --regions LG14 --memory 100


Old Code:
    # parser.add_argument('-f', '--failed_samples', help = 'Run HaplotypeCaller for 74 samples that fail VCF generation for LG7', action = 'store_true')
    # def RunFailedSamplesHaplotypeCaller(self):
    #     processes = []
    #     error_files = ['SAMEA4033321', 'SAMEA2661241', 'SAMEA4032067', 'SAMEA4033318', 'SAMEA4032108', 'SAMEA4033341', 'SAMEA4033314', 'SAMEA4032129', 'SAMEA4033252', 'SAMEA4032105', 'SAMEA3388855', 'SAMEA4032131', 'SAMEA4033248', 'SAMEA2661359', 'SAMEA1877414', 'SAMN06685761', 'SAMEA4032096', 'SAMEA2661294', 'SAMEA4032053', 'SAMEA1920091', 'SAMEA4032046', 'SAMEA4033254', 'SAMEA2661396', 'SAMEA4033283', 'SAMEA4032136', 'SAMEA2661367', 'SAMEA3388862', 'SAMEA4033301', 'SAMEA2661292', 'SAMEA4033271', 'SAMEA4032127', 'SAMEA2661414', 'SAMEA4033317', 'SAMEA2661221', 'SAMEA2661310', 'SAMEA1904322', 'SAMEA4032051', 'SAMEA4032033', 'SAMEA4033277', 'SAMEA4033307', 'SAMEA4032038', 'SAMEA4032042', 'SAMEA2661239', 'SAMN08051112', 'SAMEA4032125', 'SAMEA2661258', 'SAMEA2661306', 'SAMEA4033322', 'SAMEA1920095', 'SAMEA4033284', 'SAMEA4033304', 'SAMEA4033305', 'SAMEA4032137', 'SAMEA1920096', 'SAMEA3388871', 'SAMEA4033278', 'SAMEA4033286', 'SAMEA4032049', 'SAMEA1920092', 'SAMEA4033289', 'SAMEA4033324', 'SAMEA1904832', 'SAMEA2661339', 'SAMEA1877499', 'SAMEA2661277', 'SAMEA2661248', 'SAMEA2661280', 'SAMEA2661381', 'SAMEA1904328', 'SAMEA4033261', 'SAMEA4033287', 'SAMEA4032088', 'SAMEA4032070', 'SAMEA4032069']
    #     for file in error_files:
    #         print('Generating new GVCF file for ' + file)
    #         self.fm_obj.createSampleFiles(file)
    #         p = subprocess.Popen(['gatk', 'HaplotypeCaller', '--emit-ref-confidence', 'GVCF', '-R', self.fm_obj.localGenomeFile, '-I', self.fm_obj.localBamFile, '-O', self.fm_obj.localGVCFFile])
    #         processes.append(p)

    #         if len(processes) == len(error_files):
    #             for proc in processes:
    #                 proc.communicate()
    #             processes = []
    # if args.failed_samples:
    #     self.RunFailedSamplesHaplotypeCaller()


    # def old_RunGenomicsDBImport(self):
    #     #IMPLEMENT BATCH SIZE CODE TO BE ABLE TO STILL RUN ALL CHROMOSOMES AT ONCE
    #     processes = []
    #     for lg in self.linkage_groups:
    #         if args.local_test:
    #             p = subprocess.Popen(['gatk', '--java-options', '-Xmx' + str(self.memory) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + lg + '_database', '--intervals', os.getcwd() + '/all_lg_intervals/test_intervals/' + lg + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4', '--overwrite-existing-genomicsdb-workspace'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    #             processes.append(p)
    #         else:
    #             p = subprocess.Popen(['gatk', '--java-options', '-Xmx' + str(self.memory) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + lg + '_database', '--intervals', os.getcwd() + '/all_lg_intervals/' + lg + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4'])
    #             processes.append(p)

    #         if len(processes) == len(self.linkage_groups):
    #             for proc in processes:
    #                 proc.communicate()
    #             processes = []
    #     print('GENOMICSDBIMPORT RUN SUCCESSFULLY FOR ALL SAMPLES')

       # def old_RunGenotypeGVCFs(self):
    #     processes = []
    #     # update GenotypeGVCFs command with below annotations
    #     for lg in self.linkage_groups:
    #         if args.local_test:
    #             local_command = ['gatk', '--java-options', '-Xmx' + str(self.memory) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../' + self.fm_obj.localDatabasesDir + lg + '_database/', '-O', self.fm_obj.localOutputDir + lg + '_output.vcf', '--heterozygosity', '0.0012']
    #             local_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
    #             local_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
    #             local_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
    #             local_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
    #             p = subprocess.Popen(local_command)
    #             processes.append(p)
    #         else:
    #             genotypegvcfs_command = ['gatk', '--java-options', '-Xmx' + str(self.memory) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../../' + self.fm_obj.localDatabasesDir + lg + '_database/', '-O', self.fm_obj.localOutputDir + lg + '_output.vcf', '--heterozygosity', '0.0012']
    #             genotypegvcfs_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
    #             genotypegvcfs_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
    #             genotypegvcfs_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
    #             genotypegvcfs_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
    #             p = subprocess.Popen(genotypegvcfs_command)
    #             processes.append(p)

    #         if len(processes) == len(self.linkage_groups):
    #             for proc in processes:
    #                 proc.communicate()
    #             processes = []
    #     print('GENOTYPEGVCFS RUN SUCCESSFULLY FOR ALL SAMPLES')
        
        # commenting out the projectID code and replacing with the Ecogroup code instead. Revisit if I need to filter by projectIDs and not by Ecogroups in the future - NK, 2024.01.17
        # code block to set the ProjectIDs
        self.fm_obj.downloadData(self.fm_obj.localSampleFile_v2) # downloads the most up-to-date SampleDatabase.csv file
        s_df = pd.read_excel(self.fm_obj.localSampleFile_v2, sheet_name='SampleLevel')
        valid_project_IDs = s_df['ProjectID'].unique().tolist() # reads in SampleDatabase.csv and converts unique ProjectIDs to a list
        valid_project_IDs.extend(['All', 'custom']) # will allow 'All' which is the defualt projectID, and 'custom' which will allow a file containing vaid sample names to be valid IDs and not trip the below Exception
        for projectID in args.projectIDs: # for loop that will check each ProjectID passed to the script
            if projectID not in valid_project_IDs:
                raise Exception(projectID + ' is not a ProjectID in SampleDatabase.csv; valid Project IDs are:' + '\n' + str(valid_project_IDs))
            elif self.projectIDs == ['All'] or self.projectIDs == 'All':
                self.projectIDs = valid_project_IDs[:-2] # the -2 ensures the "All" and "custom" projectIDs are not incldued in the IDs
        # Code block to define the set of SampleIDs based on the passed ProjectIDs
        if self.projectIDs != ['custom']:
            self.fm_obj.downloadData(self.fm_obj.localAlignmentFile) # downloads the most up-to-date AlignmentDatabase.csv file from Dropbox. This file contains a list of all samples that have BAM files generated from the previous steps in the pipeline 
            self.alignment_df = pd.read_csv(self.fm_obj.localAlignmentFile) # read in the AlignmentDatabase.csv file 
            filtered_df = self.alignment_df[self.alignment_df['ProjectID'].isin(self.projectIDs)] # filter rows to only include those that are a part of the projectIDs we want to analyze 
            filtered_df = filtered_df[filtered_df['GenomeVersion'] == self.genome] # filter by genome version to eliminate duplicate samples that have been aligned to multiple genome versions
            self.sampleIDs = filtered_df['SampleID'].tolist() # set a sampleID list equal to the samples filtered by ProjectIDs. These will be processed in subsequent methods
        else:
            custom_samples_df = pd.read_csv('custom_samples.csv', header=None) # read in a custom_file.csv file (assumes no header)
            self.sampleIDs = custom_samples_df[custom_samples_df.columns[0]].values.tolist() # converts the column of samples to a list 





    def download_BAMs(self):
        for sampleID in self.sampleIDs:
            if sampleID in []:
                continue
            print('Downloading Bam file for ' + sampleID + '...')
            self.fm_obj.createSampleFiles(sampleID)
            self.fm_obj.downloadData(self.fm_obj.localBamFile)
            self.fm_obj.downloadData(self.fm_obj.localBamIndex)

    def RunHaplotypeCaller(self):
        # Find notes in the CallSmallSNVs Pipeline Notebook on Benchling. Date of entry: Wednesday March 29th, 2023
        processes = []
        for sampleID in self.sampleIDs:
            print('Generating new GVCF file for ' + sampleID)
            self.fm_obj.createSampleFiles(sampleID) # creates a sample specific directory structure so outputs are stored on the server and cloud correctly
            if args.local_test:
                self.fm_obj.localBamFile = self.fm_obj.localSampleBamDir + sampleID + '_small.all.bam'
            # try running with the basic haplotypecaller command and see if something changes
            p = subprocess.Popen(['gatk', 'HaplotypeCaller', '--emit-ref-confidence', 'GVCF', '-R', self.fm_obj.localGenomeFile, '-I', self.fm_obj.localBamFile, '-O', self.fm_obj.localGVCFFile])
            processes.append(p)

            if len(processes) == len(self.sampleIDs):
                for proc in processes:
                    proc.communicate()
                processes = []

        # if args.download_bams:
            # self.download_BAMs()
        # if args.haplotypecaller:
            # self.RunHaplotypeCaller()
"""

