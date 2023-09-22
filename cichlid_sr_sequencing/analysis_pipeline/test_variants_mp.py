import subprocess as sp
from multiprocessing import Process, freeze_support
import argparse, pdb, os
import pandas as pd
from pyfaidx import Fasta
from helper_modules.nikesh_file_manager import FileManager as FM


parser = argparse.ArgumentParser(usage='This pipeline will take in a set of unaligned bam files generated from Illumina sequencing reads, and call variants using GATK HaplotypeCaller')
parser.add_argument('reference_genome', help = 'full file path to the reference genome that reads will be aligned to for varaint calling')
parser.add_argument('-p', '--projectIDs', help = 'list of projectIDs on which to run the pipeline. Custom IDs can be provided as a csv file containing one sample per line. Save the file as "custom_samples.csv" in the directory containing this script.', nargs = '*', default = ['All'])
parser.add_argument('-i', '--import_databases', help = 'Call this flag to run GenomicsDBImport on the samples in the database', action = 'store_true')
parser.add_argument('-g', '--genotype', help = 'Call this flag to run GenotypeGVCFs on the samples in the database', action = 'store_true')
parser.add_argument('-d', '--download_data', help = 'Use this flag if you need to dwonload the GVCF files from Dropbox to include in the VCF analysis', action = 'store_true')
parser.add_argument('-r', '--regions', help = 'list of linkage groups for which analyses will run', nargs = '*', default = ['All'])
parser.add_argument('-b', '--download_bams', help = 'Download the BAM files from the cloud on which to call HaplotypeCaller', action = 'store_true')
parser.add_argument('-H', '--haplotypecaller', help = 'run the gatk HaplotypeCaller algorithm to re-generate GVCF files on which to call the pipeline', action = 'store_true')
parser.add_argument('-l', '--local_test', help = 'when this flag is called, variables will be preset to test the code locally', action = 'store_true')
parser.add_argument('-m', '--memory', help = 'How much memory, in GB, to allocate to each child process', default = 4, nargs = 1)
parser.add_argument('-u', '--unmapped', help = 'Use this flag to run -i and -g on the unmapped contigs in the genome', action = 'store_true')
args = parser.parse_args()

"""
TODO:
Lauren's mp script allows us to pass in any number of concureent processes which will run together using the mp module

Currently it works by just creating blank files, each with the name of an unmapped contig, and as mnay unmappedtontigs can be passed in at once 
    the main subprocess call is the touch command 


The code need to be edited to take in gatk subprocess calls instead
Lauren's code used pyfaidx to get the contig names. 
Since callVariants uses a hard coded contig dict (self.linkage_group_map), pyfaidx may be good to use ot gte the rest of the unmapped contig names 
Laren's code already filters out the main LGs and the mito one so we can use that here



Let's start by asking the user if they want ot run the unmapped contigs
    - create a flag and set self.linakeg_groups to the unmapped contigs

Pass in these contigs to the import_databases funnction and process with Lauren's code

May need to edit or figure out how the unmapped contig splits will work with the --intervals flags in genomicsDBImport and GenotypeGVCFs

"""

class VariantCaller:

    def __init__(self, genome, project_ids, linkage_groups, memory):
        self.genome = genome
        self.fm_obj = FM(self.genome)
        self.projectIDs = project_ids
        self.memory = memory

        # code block to set the ProjectIDs
        self.fm_obj.downloadData(self.fm_obj.localSampleFile) # downloads the most up-to-date SampleDatabase.csv file
        s_df = pd.read_csv(self.fm_obj.localSampleFile)
        valid_project_IDs = s_df['ProjectID'].unique().tolist() # reads in SampleDatabase.csv and converts unique ProjectIDs to a list
        valid_project_IDs.extend(['All', 'custom']) # will allow 'All' which is teh defualt projectID, and 'custom' which will allow a file containing vaid sample names to be valid IDs and not trip the below Exception
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
        if args.unmapped:
            self.linkage_groups = []
            contigs = Fasta(self.fm_obj.localGenomeFile).keys()
            for contig in contigs:
                if contig.startswith('NC'):
                    continue
                self.linkage_groups.append(contig)

        # pre-defining samples for local testing. Pass in the first 3 LGs only since the interval file has been created for only these.
        if args.local_test:
            self.sampleIDs = ['MC_1_m', 'SAMEA2661294', 'SAMEA2661322', 'SAMEA4032100', 'SAMEA4033261']
            self.memory = 5

    def _generate_sample_map(self):
        sampleIDs = self.sampleIDs
        with open('sample_map.txt', 'w') as fh:
            for sampleID in sampleIDs:
                self.fm_obj.createSampleFiles(sampleID)
                fh.write(sampleID + '\t' + self.fm_obj.localGVCFFile + '\n')

    def multiprocess(self, function):
        concurrent_processes = 22 # make concurrent processes an argument later and add it as a required argument later
        contigs_to_process = [self.linkage_groups[i:int(i+concurrent_processes)] for i in range(0, len(self.linkage_groups), int(concurrent_processes))]

        jobs = []
        for parallel_processes in contigs_to_process:
            j = Process(target = function, args = (parallel_processes,))
            jobs.append(j)
        for job in jobs:
            job.start()
        jobs.join()


            # ### SPLIT JOBS FUNCTION
            #     # 1. The function starts off by using the SPLICE Function and saving the lists of lists into a variable called "scopes".
            #     #       Remember: Each of the lists' lengths within "scopes" are equal to the run size you inputted
            #     # 2. A directory is then made with your desired name in your desired location.
            #     # 3. Now... (here comes the multiprocessing magic)
            #     #       For each of the lists within the "scopes" list, we are going to tell the computer to run the contigs within 
            #     #           each list, simultaneously. For example, if each list within "scopes" contains a list of length 5, then 5 
            #     #           contigs will run at once, and then the next group will run together, and so on.
            #     #       A. To do this, we must create a new list (AKA "jobs") with our scopes's lists + each of the list's contigs and their 
            #     #           commands. For example, if our scope is currently [A,B,C], then the altered list (AKA "jobs") 
            #     #           will be [do(A), do(B), do(C)]. We have to make a new list of "scopes" and not alter the current one. It's 
            #     #           simpler.
            #     #       B. Once we have successfully copied over the "jobs" list with all of the inner lists' contigs + the contigs' commands, 
            #     #           now we run it. Since it is a for-loop, we're going to run each of the scopes sequentially, and the n-number
            #     #           of items within each scope will run together.
            #     #       
            # def splitJobs(size, fastaLoc, contigKeyFileLoc, fileName):
            #     # 1
            #     scopes = splice(fastaLoc, size)

            #     # 2
            #     sp.run(["mkdir", contigKeyFileLoc + "/" + fileName])

            #     # 3A
            #     jobs = []
            #     for scope in scopes:
            #         j = Process(target=makeFile, args=(scope, contigKeyFileLoc + "/" + fileName))
            #         jobs.append(j)

            #     # 3B
            #     for j in jobs:
            #         j.start()


            # if __name__ == "__main__":
            #     splitJobs(400, '/Users/laurengsabo/Documents/mcgrath/2023/GCF_000238955.4_M_zebra_UMD2a_genomic.fna', "/Users/laurengsabo/Documents/mcgrath/2023/VSCode", "contigTrash")

    def RunGenomicsDBImport(self, regions):
        for lg in regions:
            if not args.unmapped:
                sp.run(['gatk', '--java-options', '-Xmx' + str(self.memory) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + lg + '_database', '--intervals', os.getcwd() + '/all_lg_intervals/test_intervals/' + lg + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4', '--overwrite-existing-genomicsdb-workspace'])
            else:
                sp.run(['gatk', '--java-options', '-Xmx' + str(self.memory) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + lg + '_database', '--intervals', os.getcwd() + '/unmapped_contig_intervals/' + lg + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4', '--overwrite-existing-genomicsdb-workspace'])
        print('GENOMICSDBIMPORT RUN SUCCESSFULLY FOR ALL SAMPLES')






        # #  below is the original code
        # # IMPLEMENT BATCH SIZE CODE TO BE ABLE TO STILL RUN ALL CHROMOSOMES AT ONCE
        # processes = []
        # for lg in self.linkage_groups:
        #     if args.local_test:
        #         p = sp.run(['gatk', '--java-options', '-Xmx' + str(self.memory) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + lg + '_database', '--intervals', os.getcwd() + '/all_lg_intervals/test_intervals/' + lg + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4', '--overwrite-existing-genomicsdb-workspace'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        #         processes.append(p)
        #     else:
        #         p = sp.run(['gatk', '--java-options', '-Xmx' + str(self.memory) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + lg + '_database', '--intervals', os.getcwd() + '/all_lg_intervals/' + lg + '.interval_list', '--sample-name-map', os.getcwd() + '/sample_map.txt', '--max-num-intervals-to-import-in-parallel', '4'])
        #         processes.append(p)

        #     if len(processes) == len(self.linkage_groups):
        #         for proc in processes:
        #             proc.communicate()
        #         processes = []
        

    def run_methods(self):
        self._generate_sample_map()
        if args.download_data:
            print('download data will run and GVCF files per sample will be downloaded')
            self.data_downloader()
        if args.import_databases:
            self.multiprocess(self.RunGenomicsDBImport)
        if args.genotype:
            self.RunGenotypeGVCFs()
        if args.download_bams:
            self.download_BAMs()
        if args.haplotypecaller:
            self.RunHaplotypeCaller()

if __name__ == "__main__":
    variant_caller_obj = VariantCaller(args.reference_genome, args.projectIDs, args.regions, args.memory)
    variant_caller_obj.run_methods()
    print('PIPELINE RUN SUCCESSFUL')