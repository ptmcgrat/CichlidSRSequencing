import argparse, pdb
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
        pdb.set_trace()
        # Download the GVCF and GVCF.idx files for each sample in the AlignmentDatabase
        for sampleID in self.sampleIDs:
            if sampleID in []:
                continue
            print('Downloading ' + sampleID + '...')
            self.fm_obj.createSampleFiles(sampleID)
            self.fm_obj.downloadData(self.fm_obj.localGVCFFile)
            self.fm_obj.downloadData(self.fm_obj.localGVCFFile + '.tbi')
            print('Done Downloading ' + sampleID)

    def run_methods(self):
        if args.download_data:
            self.data_downloader()
    
    
    print("Pipeline Run Success")

align_obj = AlignReads(args.reference, args.projectID, args.platform)
align_obj.run_methods()
"""
COMMAND FOR LOCAL TESTING:
/Users/kmnike/anaconda3/envs/pipeline/bin/python3 alignPB_Illumina.py Mzebra_GT1 ReferenceImprovement illumina -d

FOR TESTING ON UTAKA SERVER:
python3 alignPB_Illumina.py Mzebra_GT1 ReferenceImprovement illumina -d

"""