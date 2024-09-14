import argparse, pdb, os, subprocess, datetime, time, multiprocessing, pathlib, random
import pandas as pd
from pyfaidx import Fasta
from helper_modules.nikesh_file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage='This pipeline will take in a set of unaligned bam files generated from Illumina sequencing reads, and call variants using GATK HaplotypeCaller')
parser.add_argument('reference_genome', help = 'full file path to the reference genome that reads will be aligned to for varaint calling')
parser.add_argument('-e', '--ecogroups', help = 'one or multiple eco group names that will filter the samples on which the pipeline is run', choices = ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic','Rhamphochromis', 'Diplotaxodon', 'Riverine', 'AC', 'Non_Riverine', 'All', 'Lake_Malawi', 'Rock_Sand', 'Sand'], nargs = '*', default = ['All'])
parser.add_argument('-i', '--import_databases', help = 'Call this flag to run GenomicsDBImport on the samples in the database', action = 'store_true')
parser.add_argument('-g', '--genotype', help = 'Call this flag to run GenotypeGVCFs on the samples in the database', action = 'store_true')
parser.add_argument('-d', '--download_GVCF_data', help = 'Use this flag if you need to download the GVCF files from Dropbox to include in the VCF analysis', action = 'store_true')
parser.add_argument('-r', '--regions', help = 'list of linkage groups for which analyses will run', nargs = '*', default = ['All'])
parser.add_argument('-b', '--download_bams', help = 'Download the BAM files from the cloud on which to call HaplotypeCaller', action = 'store_true')
parser.add_argument('-H', '--efficient_haplotypecaller', help = 'use this flag to download BAM files and run HaplotypeCaller on samples', action = 'store_true')
parser.add_argument('-m', '--memory', help = 'How much memory, in GB, to allocate to each child process', default = [4], nargs = 1)
parser.add_argument('-u', '--unmapped', help = 'Use this flag to run -i and -g on the unmapped contigs in the genome', action = 'store_true')
parser.add_argument('-s', '--sampleIDs', help = 'Use this flag to customize the sampleIDs that the pipeline will run for', default = ['All'], choices = ['All', 'alignment_file', 'custom', 'v2_column'], nargs = 1)
parser.add_argument('-c', '--concat_and_index', help = 'Use this flag to concatenate the vcf files output by GenotypeGVCFs into a master_file.vcf.', action = 'store_true')
parser.add_argument('--Output', help = 'Use this flag to specify that the files for use in the pipleine or that need to be downloaded are located at /Output in Utaka', action = 'store_true')
parser.add_argument('--upload', help = 'Use this flag to upload GVCF files from the server to Dropbox', action =  'store_true')
parser.add_argument('--all_sites', help = 'use this flag if you want to gun GenotypeGVCFs in --include-non-varinat-sites mode', action = 'store_true')
parser.add_argument('--concurrent_processes', help = 'specify the number of processes to start concurrently', type = int, default = 4)
parser.add_argument('--local_test', help = 'when this flag is called, variables will be preset to test the code locally', action = 'store_true')
args = parser.parse_args()

"""
time python callVariants.py Mzebra_GT3 -i -g -m 10 -s v2_column --Output --concurrent_processes 96 2> error_phylogenyfigure_240812.txt 1> log_phylogenyfigure_240812.txt
For running on the 498 sample Cohort:
time python callVariants.py Mzebra_GT3 -d --concurrent_processes 26 -s alignment_file 2> error_download_bionano_paper_data_240914.txt 1> log_download_bionano_paper_data_240914.txt
"""

"""
NOTE: 
as of 2024 June 5, I can conda install -c bioconda gatk4 into a fresh conda env. The version installed is gatk4-4.0.5.1-0
The GATK 4.5.0.0 binary I have seems to require a new version of Java... maybe Java runtime 17 according to one post I read here: https://gatk.broadinstitute.org/hc/en-us/community/posts/17535513525915-Annotation-for-variant-calling
I have gatk working in an env called 'gatk' for now.
TODO:
1. Paralellize gatk Haplotypecaller to run much more efficiently
    Be sure to include code that will not overwrite files if they exist 
2. CODE DOES NOT EXIST ANYMORE TO DOWNLOAD BAM OR GVCF FILES OR ITS BROKEN. NEED TO FIX THIS!!
    GVCF file downloader has been fixed but I do not know how to download it straight to /Output b/c the directory structure is not mirrored on Dropbox and /Output. i think i need to add an option in FIleManager that allows me to specify that the master dir is /Output instead of /Data 
"""

class VariantCaller:
    def __init__(self, genome, sampleIDs, linkage_groups, memory, ecogroups, processes):
        self.genome = genome
        self.fm_obj = FM(self.genome)
        self.sampleIDs = sampleIDs
        self.memory = memory
        self.ecogroups = ecogroups
        self.concurrent_processes = processes
        self.current_time = datetime.datetime.now()
        # TODO: need to be able to read a list of the valid ecogroups in the Ecogroup_PTM column and let all samples be captured if 'All' is called as the ecogroup. For now, only the 'LakeMalawi' ecogroups are recognized.
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

        # block for defining sampleIDs if you want specific samples from a SampleDatabase column, using the alignmnetdatabase, or anythign else
        if self.sampleIDs == ['All']:
            self.fm_obj.downloadData(self.fm_obj.localSampleFile_v2) # downloads the most up-to-date SampleDatabase_v2.xlsx file
            s_df = pd.read_excel(self.fm_obj.localSampleFile_v2, sheet_name='SampleLevel') # reads in the SampleLevel sheet from SampleDatabase_v2.xlsx
            eg_filtered_df = s_df[s_df['Ecogroup_PTM'].isin(self.ecogroups)] # filter samples that match the ecogroups in the analysis
            self.sampleIDs = eg_filtered_df['SampleID'].to_list()
        elif self.sampleIDs == ['alignment_file']:
            self.fm_obj.downloadData(self.fm_obj.localAlignmentFile) # download the AlignmentDatabase.csv file 
            s_df = pd.read_csv(self.fm_obj.localAlignmentFile)
            self.sampleIDs = s_df[s_df['GenomeVersion'] == self.genome].SampleID.to_list() # get sampleIDs by filtering on
        elif self.sampleIDs == ['custom']:
            # TODO: allow a file with custom sample names to be used as input for the pipleine. 
            print('ERROR: CUSTOM SAMPLE FILE NOT YET IMPLEMENTED')
            exit
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
            self.sampleIDs = ['CJ_2204_m', 'CV-006-m', 'LA_3006_m', 'MC-008-m', 'OC-001-m', 'small_file1', 'small_file2', 'small_file3', 'small_file4']
            self.memory = [1]
            self.linkage_groups = ['NC_036780.1', 'NC_036781.1', 'NC_036782.1']
            self.concurrent_processes = 10
        print(f"Number of samples for this pipeline run is {len(self.sampleIDs)}")

    def _generate_sample_map(self):
        sampleIDs = self.sampleIDs
        with open('sample_map.txt', 'w') as fh:
            for sampleID in sampleIDs:
                self.fm_obj.createSampleFiles(sampleID)
                if args.local_test:
                    fh.write(sampleID + '\t' + self.fm_obj.localTestGVCFFile + '\n')
                elif args.Output:
                    fh.write(sampleID + '\t' + self.fm_obj.StorageGVCFFile + '\n')
                else:
                    fh.write(sampleID + '\t' + self.fm_obj.localGVCFFile + '\n')

    def GVCF_downloader(self, sampleID):
        # Download the GVCF and GVCF.idx files for each sample in self.sampleIDs
        self.fm_obj.createSampleFiles(sampleID)
        if args.local_test: # for downloading files locally and testing if the script works 
            if pathlib.Path(self.fm_obj.localTestGVCFFile).exists():
                print('GCVF file for ' + sampleID + ' exists. Skipping sample...')
            else:
                print(f"Starting GVCF and GVCF index download for sample {sampleID} at {self.current_time}")
                self.fm_obj.downloadData(self.fm_obj.localTestGVCFFile)
                self.fm_obj.downloadData(self.fm_obj.localTestGVCFIndex)
                print(f"Download of GVCF and GVCF index for sample {sampleID} complete at {self.current_time}")
        elif args.Output: # NOTE: DOES NOT WORK. THE DROPBOX DTRUCTURE MUST BE MIRRORED. SO IF I WANT TO DOWNLOAD INTO /OUTPUTS, IT MUST EXIST IN /OUTPUTS IN DROPBOX. 
            if pathlib.Path(self.fm_obj.StorageGVCFFile).exists():
                print(print('GCVF file for ' + sampleID + ' exists at /Output. Skipping sample...'))
            else:
                print(f"Starting GVCF and GVCF index download for sample {sampleID} at {self.current_time}. Sample is being stored at /Outputs")
                self.fm_obj.downloadData(self.fm_obj.StorageGVCFFile)
                self.fm_obj.downloadData(self.fm_obj.StorageGVCFIndex)
                print(f"Download of GVCF and GVCF index for sample {sampleID} complete at {self.current_time}")
        else:
            if pathlib.Path(self.fm_obj.localGVCFFile).exists():
                print('GCVF file for ' + sampleID + ' exists. Skipping sample...')
            else:
                print(f"Starting GVCF and GVCF index download for sample {sampleID} at {self.current_time}")
                self.fm_obj.downloadData(self.fm_obj.localGVCFFile)
                self.fm_obj.downloadData(self.fm_obj.localGVCFIndex)
                print(f"Download of GVCF and GVCF index for sample {sampleID} complete at {self.current_time}")

    def uploadGVCFs(self, sampleID):
        # upload GVCF file and index for sampleID in self.sampleIDs if it doesn't already exist on Dropbox. 
        self.fm_obj.createSampleFiles(sampleID)
        if args.local_test: # code for testing using small local files. 
            # check to see if the file has already been uploaded to Dropbox. If it does. Skip it
            subprocess.check_output(['rclone', 'lsf', self.fm_obj.rcloneRemote + self.fm_obj.localBamFile], encoding='utf-8')
            print(f"Checking if {sampleID}'s GVCF file and Index exists on Dropbox...")
            cloud_gvcf_path = self.fm_obj.localTestGVCFFile.replace(self.fm_obj.localMasterDir, self.fm_obj.cloudMasterDir)
            cloud_gvcf_index_path = self.fm_obj.localTestGVCFIndex.replace(self.fm_obj.localMasterDir, self.fm_obj.cloudMasterDir)

            gvcf_code = subprocess.run(['rclone', 'lsf', cloud_gvcf_path]).returncode
            if gvcf_code == 0:
                print(f"Rclone found the GVCF file for {sampleID} on Dropbox. Skipping")
            else:
                print(f"GVCF file for {sampleID} is not on dropbox. Starting upload for GVCF file at {self.current_time}")
                self.fm_obj.uploadData(self.fm_obj.localTestGVCFFile)
                print(f"GVCF file upload for {sampleID} complete at {self.current_time}")
                
            index_code = subprocess.run(['rclone', 'lsf', cloud_gvcf_index_path]).returncode
            if index_code == 0:
                print(f"Rclone found the GVCF Index for {sampleID} on Dropbox. Skipping...")
            else:
                print(f"Index file for {sampleID} is not on dropbox. Starting upload for GVCF Index at {self.current_time}")
                self.fm_obj.uploadData(self.fm_obj.localTestGVCFIndex)
                print(f"GVCF Inxex upload for {sampleID} complete at {self.current_time}")
        else: # code for uploading server files
            print(f"Checking if {sampleID}'s GVCF file and Index exists on Dropbox...")
            cloud_gvcf_path = self.fm_obj.localGVCFFile.replace(self.fm_obj.localMasterDir, self.fm_obj.cloudMasterDir)
            cloud_gvcf_index_path = self.fm_obj.localGVCFIndex.replace(self.fm_obj.localMasterDir, self.fm_obj.cloudMasterDir)

            gvcf_code = subprocess.run(['rclone', 'lsf', cloud_gvcf_path]).returncode
            if gvcf_code == 0:
                print(f"Rclone found the GVCF file for {sampleID} on Dropbox. Skipping")
            else:
                print(f"GVCF file for {sampleID} is not on dropbox. Starting upload for GVCF file at {self.current_time}")
                self.fm_obj.uploadData(self.fm_obj.localGVCFFile)
                print(f"GVCF file upload for {sampleID} complete at {self.current_time}")
                
            index_code = subprocess.run(['rclone', 'lsf', cloud_gvcf_index_path]).returncode
            if index_code == 0:
                print(f"Rclone found the GVCF Index for {sampleID} on Dropbox. Skipping...")
            else:
                print(f"Index file for {sampleID} is not on dropbox. Starting upload for GVCF Index at {self.current_time}")
                self.fm_obj.uploadData(self.fm_obj.localGVCFIndex)
                print(f"GVCF Inxex upload for {sampleID} complete at {self.current_time}")

    def EfficientHaplotypeCaller(self, sample):
        """
        NOTE: IMPLEMENT A FILE CHECKER TO MAKE SURE EXISTING FILES DO NOT GET DELETED. INCLUDE CODE FOR UPLOADING GVCF FILES AND INDICIES AS WELL 
        """
        print(f"Processing for sample {sample} started at {self.current_time}")
        # to use multiprocess, this function will sequenctially need to perform BAM/BAI download, GVCF generation, and then deleting the BAM/BAI files
        if args.local_test:
            self.fm_obj.createSampleFiles(sample)
            if pathlib.Path(self.fm_obj.localTestBamFile).exists():
                print(f"BAM file for {sample} exists. Skipping download for this sample...")
            else:
                print(f"Downloading BAM and BAM index for {sample}...")
                print(f"BAM Download for sample {sample} started at {self.current_time}")
                self.fm_obj.downloadData(self.fm_obj.localTestBamFile)
                self.fm_obj.downloadData(self.fm_obj.localTestBamIndex)
                print(f"BAM Download for sample {sample} ended at {self.current_time}...")

            if pathlib.Path(self.fm_obj.localTestGVCFFile).exists():
                print('GCVF file for ' + sample + ' exists. Skipping sample... ')
            else:
                print('Generating GCVF file for ' + sample + '...')
                subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'HaplotypeCaller', '--emit-ref-confidence', 'GVCF', '-R', self.fm_obj.localGenomeFile, '-I', self.fm_obj.localTestBamFile, '-O', self.fm_obj.localTestGVCFFile])

            # commenting out bam removal code while GCVF code isn't yet working
            # print('Removing BAM file for ' + saple + '...')
            # os.remove(self.fm_obj.localTestBamFile)
            # os.remove(self.fm_obj.localTestBamIndex)
        else:
            self.fm_obj.createSampleFiles(sample) # create the file structure for the sample
            if pathlib.Path(self.fm_obj.localGVCFIndex).exists(): # If the index for that file exists, then the GVCF file is fine and Haplotypecaller does not need to be rerun 
                print(f"GVCF and Index created for sample {sample}. Skipping HaplotypeCaller run for this sample")
            else: # if index doesn't exist, process the sample
                print(f"No GVCF Index exists for sample {sample}. Starting processing")
                if pathlib.Path(self.fm_obj.localGVCFFile).exists(): # if a GVCF file is found for the sample without the index, that means it will be a partial one. Thus, it needs to be removed. 
                    print(f"Partial GVCF with no Index Detected. Deleting the GVCF file for sample {sample} and regenerating it.")
                    os.remove(self.fm_obj.localGVCFFile)
                new_gcvf_created = False
                if pathlib.Path(self.fm_obj.localBamFile).exists(): # If the BAM file exists for the sample, do not download it.
                    print('BAM file for ' + sample + ' exists. Skipping download for this sample... ')
                else: # download BAM since it does not exist on the server.
                    print('Downloading BAM and BAM index for ' + sample + '...')
                    self.fm_obj.downloadData(self.fm_obj.localBamFile)
                    self.fm_obj.downloadData(self.fm_obj.localBamIndex)
                
                # if the GVCFindex did not exist and a partial file was found, at this stage, the partial file would have been deleted and a new GVCF would need to be generated with the below code:
                print(f"Generating GCVF file for sample {sample}...")
                subprocess.run(['gatk', '--java-options', '-Xmx' +  str(self.memory[0]) + 'G', 'HaplotypeCaller', '--emit-ref-confidence', 'GVCF', '-R', self.fm_obj.localGenomeFile, '-I', self.fm_obj.localBamFile, '-O', self.fm_obj.localGVCFFile])
                new_gcvf_created = True

                if new_gcvf_created:
                    print('Removing BAM file for ' + sample + '...')
                    os.remove(self.fm_obj.localBamFile)
                    os.remove(self.fm_obj.localBamIndex)

        print(f"Sample {sample} finished processing at {self.current_time}")

    def RunGenomicsDBImport(self, interval):
        print(f"Processing for interval {interval} started at {self.current_time}")
        if args.local_test: # for testing code locally. Keep or remove --overwrite-existing-genomicsdb-workspace as needed. Default is that it's true
            subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + interval + '_database', '--intervals', os.getcwd() + '/GT3_intervals/' + interval + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt','--overwrite-existing-genomicsdb-workspace'])
        elif not args.unmapped: # for running on normal 96 intervals on the servers. Note that --overwrite-existing-genomicsdb-workspace is false, so nothing accidentally gets overwritten on the server. 
            subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + interval + '_database', '--intervals', os.getcwd() + '/GT3_intervals/' + interval + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt'])
          # subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + lg + '_database', '--intervals', os.getcwd() + '/unmapped_contig_intervals/' + lg + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4']) # this is code used if lgs are passed instaed of intervals. This looks to the all_lg_intervals dir and runs based on the old parallelization methods.
        else: # If running on unmapped contigs, use intervals described in the unmapped intervals dir within analysis_pipeline/
            subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + interval + '_database', '--intervals', os.getcwd() + '/unmapped_contig_intervals/' + interval + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt'])
          
        print(f"Interval {interval} finished processing at {self.current_time}")

    def RunGenotypeGVCFs(self, interval):
        # You can ignore this function for now, but I implemented the parallelization code into it
        # check if a vcf index file already exists for the intervalk you're working on  if it does, skip it and move to another interval
        if pathlib.Path(self.fm_obj.localOutputDir + interval + '_output.vcf.gz.idx').exists():
            print(f"VCF file's index for interval {interval} found. Skipping GenotypeGVCFs for this interval")
        else:
            print(f"Processing for interval {interval} started at {self.current_time}")
            if args.local_test: # section is for local_testing code
                if args.all_sites: # This section is for running locally with --include-non-variant-sites true. IMPORTANT: This option is true only in gatk v4.3.0.0 - 24.07.23 NK. This flag is --all-sites on older versions!! Be careful with this option depending on what gatk version you're using!!
                    local_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../' + self.fm_obj.localDatabasesDir + interval + '_database/', '-O', self.fm_obj.localOutputDir + interval + '_output.vcf.gz', '--heterozygosity', '0.00175', '--include-non-variant-sites', 'true', '--intervals', os.getcwd() + '/GT3_intervals/' + interval + '.interval_list'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
                    local_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
                    local_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
                    local_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
                    local_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
                    subprocess.run(local_command)
                else: # this section is for running locally in normal mode, i.e. WITHOUT --include-non-variant-sites
                    local_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../' + self.fm_obj.localDatabasesDir + interval + '_database/', '-O', self.fm_obj.localOutputDir + interval + '_output.vcf.gz', '--heterozygosity', '0.00175'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
                    local_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
                    local_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
                    local_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
                    local_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
                    subprocess.run(local_command)
            elif not args.unmapped: # For running the standard 96 intervals (not unmappped and not local_test). gendb path is to the gendb located on the Utaka server. Will need to change if running on Mzebra 
                if args.all_sites: # This section is for running the 96 intervals with --include-non-variant-sites true. IMPORTANT: This option is true only in gatk v4.3.0.0 - 24.07.23 NK. This flag is --all-sites on older versions!! Be careful with this option depending on what gatk version you're using!!
                    genotypegvcfs_unmapped_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../../' + self.fm_obj.localDatabasesDir + interval + '_database/', '-O', self.fm_obj.localOutputDir + interval + '_output.vcf.gz', '--heterozygosity', '0.00175', '--include-non-variant-sites', 'true', '--intervals', os.getcwd() + '/GT3_intervals/' + interval + '.interval_list'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
                    genotypegvcfs_unmapped_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
                    genotypegvcfs_unmapped_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
                    genotypegvcfs_unmapped_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
                    genotypegvcfs_unmapped_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
                    subprocess.run(genotypegvcfs_unmapped_command)
                else: # this section is for running the 96 intervals in normal mode, i.e. WITHOUT --include-non-variant-sites
                    genotypegvcfs_unmapped_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../../' + self.fm_obj.localDatabasesDir + interval + '_database/', '-O', self.fm_obj.localOutputDir + interval + '_output.vcf.gz', '--heterozygosity', '0.00175'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
                    genotypegvcfs_unmapped_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
                    genotypegvcfs_unmapped_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
                    genotypegvcfs_unmapped_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
                    genotypegvcfs_unmapped_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
                    subprocess.run(genotypegvcfs_unmapped_command)
            else: # This section is for running the unmapped contigs.
                if args.all_sites: # This section is for running unmapped contigs with --include-non-variant-sites true. IMPORTANT: This option is true only in gatk v4.3.0.0 - 24.07.23 NK. This flag is --all-sites on older versions!! Be careful with this option depending on what gatk version you're using!!
                    genotypegvcfs_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../../' + self.fm_obj.localDatabasesDir + interval + '_database/', '-O', self.fm_obj.localOutputDir + interval + '_output.vcf.gz', '--heterozygosity', '0.00175', '--include-non-variant-sites', 'true', '--intervals', os.getcwd() + '/unmapped_contig_intervals/' + interval + '.interval_list'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
                    genotypegvcfs_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
                    genotypegvcfs_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
                    genotypegvcfs_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
                    genotypegvcfs_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
                    subprocess.run(genotypegvcfs_command)
                else: # this section is for running unmapped contigs in normal mode, i.e. WITHOUT --include-non-variant-sites
                    genotypegvcfs_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', 'gendb://../../../../../../' + self.fm_obj.localDatabasesDir + interval + '_database/', '-O', self.fm_obj.localOutputDir + interval + '_output.vcf.gz', '--heterozygosity', '0.00175'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
                    genotypegvcfs_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
                    genotypegvcfs_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
                    genotypegvcfs_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
                    genotypegvcfs_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
                    subprocess.run(genotypegvcfs_command)
            print(f"Interval {interval} finished processing at {self.current_time}")

    def _merge_vcfs(self):
        vcfConcatDir = self.fm_obj.localOutputDir + 'vcf_concat_output/'
        file_num = 1
        for file in os.listdir(vcfConcatDir):
            if file.startswith('master_file') and file.endswith('.vcf.gz'): # If any other compressed master_files (but not index files) exist, then they do not get overwritten
                file_num += 1

        if args.local_test: # local testing interval zipping
            intervals = list(range(1,97))
            with open('compression_file_list.txt', 'w') as fh:
                for interval in intervals:
                    fh.write(f"{self.fm_obj.localOutputDir}{str(interval)}_output.vcf.gz\n")
            if file_num == 1:
                print(f"Starting compression of individual VCF files at {self.current_time}")
                subprocess.run(['bcftools', 'concat', '-f', os.getcwd() + '/compression_file_list.txt', '-Wtbi', '-O', 'z', '-o', self.fm_obj.localOutputDir + 'vcf_concat_output/master_file.vcf.gz', '--threads', str(self.concurrent_processes)])
                print(f"master_file.vcf.gz finished concatenating and inxing at {self.current_time}")
            else:
                print(f"Another master_file has been found in the vcf_concat_output directory. Writing master_file{file_num}.vcf.gz at {self.current_time}")
                subprocess.run(['bcftools', 'concat', '-f', os.getcwd() + '/compression_file_list.txt', '-Wtbi', '-O', 'z', '-o', self.fm_obj.localOutputDir + 'vcf_concat_output/master_file'+ str(file_num) + '.vcf.gz', '--threads', str(self.concurrent_processes)])

        elif not args.unmapped: # normal 96 interval zipping
            intervals = list(range(1,97))
            with open('compression_file_list.txt', 'w') as fh:
                for interval in intervals:
                    fh.write(f"{self.fm_obj.localOutputDir}{str(interval)}_output.vcf.gz\n")
            if file_num == 1:
                print(f"Starting compression of individual VCF files at {self.current_time}")
                subprocess.run(['bcftools', 'concat', '-f', os.getcwd() + '/compression_file_list.txt', '-Wtbi', '-O', 'z', '-o', self.fm_obj.localOutputDir + 'vcf_concat_output/master_file.vcf.gz', '--threads', str(self.concurrent_processes)])
                print(f"master_file.vcf.gz finished concatenating and inxing at {self.current_time}")
            else:
                print(f"Another master_file has been found in the vcf_concat_output directory. Writing master_file{file_num}.vcf.gz at {self.current_time}")
                subprocess.run(['bcftools', 'concat', '-f', os.getcwd() + '/compression_file_list.txt', '-Wtbi', '-O', 'z', '-o', self.fm_obj.localOutputDir + 'vcf_concat_output/master_file'+ str(file_num) + '.vcf.gz', '--threads', str(self.concurrent_processes)])

        else: # unmapped contig vcf file output concatenation and compression
            self.fasta = Fasta(self.fm_obj.localGenomeFile)
            self.unmapped_contigs = [contig for contig in self.fasta.keys() if contig.startswith('NW') or contig == 'NC_027944.1' or contig == 'ptg000146l_obj_unaligned']
            with open('compression_file_list.txt', 'w') as fh:
                for contig in self.unmapped_contig:
                    fh.write(f"{self.fm_obj.localOutputDir}{contig}\n")
            if file_num == 1:
                print(f"Starting compression of individual VCF files at {self.current_time}")
                subprocess.run(['bcftools', 'concat', '-f', os.getcwd() + '/compression_file_list.txt', '-Wtbi', '-O', 'z', '-o', self.fm_obj.localOutputDir + 'vcf_concat_output/master_file.vcf.gz', '--threads', str(self.concurrent_processes)])
                print(f"master_file.vcf.gz finished concatenating and inxing at {self.current_time}")
            else:
                print(f"Another master_file has been found in the vcf_concat_output directory. Writing master_file{file_num}.vcf.gz at {self.current_time}")
                subprocess.run(['bcftools', 'concat', '-f', os.getcwd() + '/compression_file_list.txt', '-Wtbi', '-O', 'z', '-o', self.fm_obj.localOutputDir + 'vcf_concat_output/master_file'+ file_num + '.vcf.gz', '--threads', str(self.concurrent_processes)])

    def multiprocess(self, function, sample_type):
        # TODO: find a way to load processes that would take the lonegst time to go first
        # Author: Lauren Sabo; edits made by NK
        if sample_type == 'lg':
            inputs = self.linkage_groups
        elif sample_type == 'sampleID':
            inputs = self.sampleIDs
        elif sample_type == 'interval':
            intervals = list(range(1,97))
            inputs = list(map(str, intervals))
        elif sample_type == 'unmapped':
            # define the unmapped contigs
            self.fasta = Fasta(self.fm_obj.localGenomeFile)
            self.unmapped_contigs = [contig for contig in self.fasta.keys() if contig.startswith('NW') or contig == 'NC_027944.1' or contig == 'ptg000146l_obj_unaligned']
            inputs = self.unmapped_contigs
        if args.all_sites and not args.local_test: # if running GenotypeGVCFs in all-sites mode, 96 concurrent processes will deplete all available memory, so change the max concurrent processes to be 48 instead. Makes sure this only applies on server.
            print('All sites run detected. Limiting the number of concurrent processes to 48 to ensure each process has enough RAM to complete.')
            self.concurrent_processes = 48
        elif args.all_sites and args.local_test:
            self.concurrent_processes = 5
        concurrent_processes = min(self.concurrent_processes, len(inputs))
        try:
            with multiprocessing.Pool(processes=concurrent_processes) as pool:
                pool.map(function, inputs)
        except Exception as e:
            print(f"Error occurred during multiprocessing: {e}")

    def run_methods(self):
        self._generate_sample_map()
        if args.download_GVCF_data:
            self.multiprocess(self.GVCF_downloader, 'sampleID')
        if args.upload:
            self.multiprocess(self.uploadGVCFs, 'sampleID')
        if args.efficient_haplotypecaller:
            self.multiprocess(self.EfficientHaplotypeCaller, 'sampleID')
        if args.import_databases and not args.unmapped:
            self.multiprocess(self.RunGenomicsDBImport, 'interval')
        if args.import_databases and args.unmapped:
            self.multiprocess(self.RunGenomicsDBImport, 'unmapped')
        if args.genotype and not args.unmapped:
            self.multiprocess(self.RunGenotypeGVCFs, 'interval')
        if args.genotype and args.unmapped:
            self.multiprocess(self.RunGenomicsDBImport, 'unmapped')
        if args.concat_and_index:
            self._merge_vcfs()

if __name__ == "__main__":
    variant_caller_obj = VariantCaller(args.reference_genome, args.sampleIDs, args.regions, args.memory, args.ecogroups, args.concurrent_processes)
    variant_caller_obj.run_methods()
    print('PIPELINE RUN COMPLETE')


"""
time python callVariants.py Mzebra_GT3 -b  -H -a --concurrent_processes 24 -m 40 2> error_sd_rerun_240721.txt 1> log_sd_rerun240721.txt
time python callVariants.py Mzebra_GT3 -b -a --concurrent_processes 24 -m 40
time python callVariants.py Mzebra_GT3 -g -a --concurrent_processes 96 -m 10
time python callVariants.py Mzebra_GT3 --concurrent_processes 96 -m 10 --temp_zip 2> error_zip_allsites_vcfs_240802.txt 1> log_zip_allsites_vcfs_240802.txt

time python callVariants.py Mzebra_GT3 --concurrent_processes 96 -m 10 --concat_and_index 2> error_concat_all_sites_240805.txt  1> log_concat_all_sites_240805.txt


time python callVariants.py Mzebra_GT3 --local_test --concat_and_index --memory 1 --concurrent_processes 10

"""

"""
Old Code & Unused functions:
    def mp_test_function(self, interval):
        print(f"Task {interval} started at {self.current_time}")
        # Simulate some work with sleep
        time.sleep(interval)
        print(f"Task {interval} finished at {self.current_time}")
    def new_test(self, sample_name):
        print(f"Task {sample_name} started at {self.current_time}")
        # Simulate some work with sleep
        print(sample_name)
        time.sleep(random.randint(1,7))
        print(f"Task {sample_name} finished at {self.current_time}")
            below code was used to test the new_test() and mp_test_function() functions in the multiprocess function.
            if not args.local_test:
                intervals = list(range(1,97))
                inputs = list(map(str, intervals))
            else:
                inputs = [5, 3, 8, 2, 6, 1, 7, 4]
        Below code goes in the run_methods() function
        if args.local_test:
            self.multiprocess(self.mp_test_function, 'interval')
        if args.local_test:
            self.multiprocess(self.new_test, 'sampleID')
"""