import argparse, pdb, subprocess
from helper_modules.nikesh_file_manager import FileManager as FM
import pandas as pd

parser = argparse.ArgumentParser(usage = 'This script will take in a project ID and a reference genome and use BWA or minimap2 to align Illumina paired-end reads or PacBio reads respectively.')
parser.add_argument('reference', help = 'name of reference genome to which the reads will be aligned', choices  =['Mzebra_UMD2a', 'Mzebra_GT1'])
parser.add_argument('projectID', help = 'list of projectIDs on which to run the pipeline', nargs = '*', default = ['All'])
parser.add_argument('platform', help = 'name of the platform which generate dthe reads. For now, the script will support running reads for one platform at a time', choices = ['illumina', 'pacbio'])
parser.add_argument('-d', '--download_data', help = 'Use this flag if you need to Download Data from Dropbox to include in the analysis', action = 'store_true')
args = parser.parse_args()

"""
TO DO:
- For Illumina reads, simply download teh ubam file from the Dropbox, and pipe, in parallel, into bwa using the mem alrogithm found in the alignment_worker.py. 
- Downloaded reads will need to be aligned using BWA if the reads are Illumina & minimap2 if from PacBio
    - utilize the alignFastQ.py script & alignmentworker.py scripts to see how Patrick did this for Illumina reads (especially since there's 4 files per sample; L1 + L2 pe reads).

"""

class AlignReads:
    def __init__(self, genome, project_ids, platform):
        self.genome = genome
        self.platform = platform
        self.fm_obj = FM(self.genome)

        self.projectIDs = project_ids
        if self.projectIDs == ['All']:
            self.projectIDs = ['ReferenceImprovement', 'RockSand_v1', 'PRJEB15289', 'PRJEB1254','BrainDiversity_s1', 'BigBrain']
        self.alignment_df = pd.read_csv(self.fm_obj.localAlignmentFile)
        filtered_df = self.alignment_df[self.alignment_df['ProjectID'].isin(self.projectIDs)]
        self.sampleIDs = filtered_df['SampleID'].tolist()
        

    def data_downloader(self):
        print('Downloading Unmapped Alignment File')
        self.fm_obj.downloadData(self.fm_obj.localAlignmentFile)
        print('Downloading most recent Genome Dir')
        self.fm_obj.downloadData(self.fm_obj.localGenomeDir)
        # Download the GVCF and GVCF.idx files for each sample in the AlignmentDatabase
        for sampleID in self.sampleIDs:
            if sampleID in []:
                continue
            print('Downloading ' + sampleID + '...')
            self.fm_obj.createSampleFiles(sampleID)
            self.fm_obj.downloadData(self.fm_obj.localUnmappedBamFile)
            print('Done Downloading unmapped BAM for ' + sampleID)

    def uBamtoBam(self):        
        for sampleID in self.sampleIDs:
            self.fm_obj.createSampleFiles(sampleID)
            command1 = ['gatk', 'SamToFastq', '-I', self.fm_obj.localUnmappedBamFile, '--FASTQ', '/dev/stdout', '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2']
            command1 += ['--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', self.fm_obj.localTempDir]

            # Second command aligns fastq data to reference
            command2 = ['bwa', 'mem', '-t', '12', '-M', '-p', self.fm_obj.localGenomeFile, '/dev/stdin']

            # Final command reads read group information to aligned bam file and sorts it
            # Figure out how to keep hard clipping
            command3 = ['gatk', 'MergeBamAlignment', '-R', self.fm_obj.localGenomeFile, '--UNMAPPED_BAM', self.fm_obj.localUnmappedBamFile, '--ALIGNED_BAM', '/dev/stdin']
            command3 += ['-O', t_bam, '--ADD_MATE_CIGAR', 'true', '--CLIP_ADAPTERS', 'false', '--CLIP_OVERLAPPING_READS', 'true']
            command3 += ['--INCLUDE_SECONDARY_ALIGNMENTS', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1', '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant']
            command3 += ['--ATTRIBUTES_TO_RETAIN', 'XS', '--TMP_DIR', self.fm_obj.localTempDir]
            t_bam = self.fm_obj.localTempDir + sampleID + '.' + str(i) + '.sorted.bam'
            p1 = subprocess.Popen(command1, stdout=subprocess.PIPE, stderr = subprocess.DEVNULL)
            p2 = subprocess.Popen(command2, stdin = p1.stdout, stdout = subprocess.PIPE, stderr = subprocess.DEVNULL)
            p1.stdout.close()
            p3 = subprocess.Popen(command3, stdin = p2.stdout, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
            p2.stdout.close()
            output = p3.communicate()

    def run_methods(self):
        if args.download_data:
            self.data_downloader()
        self.uBamtoBam()

align_obj = AlignReads(args.reference, args.projectID, args.platform)
align_obj.run_methods()
print('PIPELINE RUN SUCCESSFUL')
"""
COMMAND FOR LOCAL TESTING:
/Users/kmnike/anaconda3/envs/pipeline/bin/python3 alignPB_Illumina.py Mzebra_GT1 ReferenceImprovement illumina -d

FOR TESTING ON UTAKA SERVER:
python3 alignPB_Illumina.py Mzebra_GT1 ReferenceImprovement illumina -d

gatk SamToFastq -I CV_1_m.unmapped.bam --FASTQ /dev/stdout --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true --NON_PF true --TMP_DIR ./Temp

"""