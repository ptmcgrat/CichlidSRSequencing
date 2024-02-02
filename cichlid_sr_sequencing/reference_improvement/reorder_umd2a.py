import argparse, pdb, re, pathlib
import pandas as  pd
from pyfaidx import Fasta
from helper_modules.nikesh_file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage = 'This script will reorder the nucleotides in a reference genome as informed by OGM Data')
parser.add_argument('--genome', help = 'name of the reference to reorganize', default='Mzebra_UMD2a')
parser.add_argument('-o', '--output_dir', help = 'absolute filepath to an output directory', nargs=1)
parser.add_argument('-r', '--regions', help = 'list of linakge groups where script will operate. All linkage groups will be output into the same directory', nargs = '*', default = ['All'])
parser.add_argument('-f', '--file_name', help = 'name of the output reorganized genoem file.', default='Mzebra_GT1.fna')
parser.add_argument('-w', '--write_unmapped_contigs', help = 'flag that will write in the unmapped contigs if called', action='store_true', default=False)
args = parser.parse_args()

"""
Link top isntructions for how to get pandas to read the google sheets document
https://medium.com/geekculture/2-easy-ways-to-read-google-sheets-data-using-python-9e7ef366c775#c6bb

URL backup:
url = f'https://docs.google.com/spreadsheets/d/e/2PACX-1vRLiFIHXSXrljCLzO58qCwqcvdU3sw-lNAA8P_riJcCS14mSQRcNjZMRnsW3eKrAJcZENUbpEJoiakH/pub?output=xlsx'
The sheet is called "UMD2a_remap_v1" in my google sheets.


TODO:
The reordering of the main linkage groups is messing up. Originally, the script was popping the final contig out. I think this was a remnant from when I used to add in NNNs as their own elements in the self.remapped_genome list
I removed this, yet the expected length of the LG still doesn't match with, and is less than the original length of the genome. 

I think that the script is fixed at this point, but now , I need to re-check whether each LG actually has all of the contigs actually present in the remapping instriuctions by going tthrough the google sheet in detail, per LG. 
There is still 

There is a minor issue where LGs 1, 6, 8, & 15 have an extra N at the end of their sequences. I think this may stem from the fact that newlines are inserted every 80 bases. Sometimes, the last 101 bases have 2 newlines, so one takes up where an N would be removed.
Resolve by adding newlines after splitting, or remove Ns before adding newlines.
"""

class ReorderGenome:
    def __init__(self, genome, output_directory, linkage_groups, output_filename):
        # The linkage_group_map attribute is a hard coded list of LG names
        self.fm_obj = FM(genome)
        self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1', 
                             'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1', 
                             'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1', 
                             'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1'}
        self.genome = self.fm_obj.localGenomeFile
        if not pathlib.Path(output_directory).exists():
            pathlib.Path(output_directory).mkdir(parents=True, exist_ok=True)
            self.out_dir = output_directory
        else:
            self.out_dir = output_directory

        # code to take in simplified region inputs like "LG1", "LG7", etc, as well as check for invalid inputs & duplciate LG names.
        self.linkage_groups = linkage_groups
        if self.linkage_groups == ['All']:
            self.linkage_groups = [*self.linkage_group_map.keys()]
        for region in self.linkage_groups:
            if region not in self.linkage_group_map.keys():
                raise Exception(region + ' is not a valid option')
            duplicate_set_test = set(self.linkage_groups)
            if len(duplicate_set_test) != len(self.linkage_groups):
                raise Exception('A repeat region has been provided')
        self.filename = output_filename

    def make_remap_instructions(self, linkage_group):
        url = f'https://docs.google.com/spreadsheets/d/e/2PACX-1vRLiFIHXSXrljCLzO58qCwqcvdU3sw-lNAA8P_riJcCS14mSQRcNjZMRnsW3eKrAJcZENUbpEJoiakH/pub?output=xlsx'
        self.df = pd.read_excel(url, sheet_name = linkage_group)
        remapping_instructions = []
        for index, rows in self.df.iterrows():
            remapping_instructions.append([rows.contig_name, rows.Start, rows.Stop, rows.Direction])
        return(remapping_instructions)

    def _rearrange_genome(self):
        self.pyfaidx_obj = Fasta(self.genome) # load in the reference genome as a Fasta class object from pyfaidx
        allowed_region_classes = ['Normal', 'Inversion', 'Gap']
        with open(self.out_dir + '/' + self.filename, 'w') as fh:
            self.contigs_used_list = []
            self.contigs_unused_list = []
            for lg in self.linkage_groups:
                print('Getting Instructions for how to rebuild ' + lg + '...')

                linkage_group_instructions = self.make_remap_instructions(lg)
                fh.write('>' + self.linkage_group_map[lg] + '\n')
                self.remapped_genome = []
                for region in linkage_group_instructions:
                    if region[3] not in allowed_region_classes:
                        raise Exception(f"the region '{region[3]}' is not a valid region for genome reordering. Fix the mistakes in the excel sheet for {lg} for the {region[1]} to {region[2]} interval.")
                    if region[0] == 'Gap':
                        print('Gap detected in instructions... Adding 200 Ns to the reference to indicate presence of a gap')
                        self.remapped_genome.append('N'*200)
                    else:
                        self.bases_to_write = self.pyfaidx_obj[region[0]][int(region[1]):int(region[2])]
                        if region[0] not in self.contigs_used_list:
                            self.contigs_used_list.append(region[0])
                        if region[3] == "Normal":
                            print('Normal contig found. Adding bases in the same orientation.')
                            self.normal_transformation(self.bases_to_write)
                        elif region[3] == 'Inversion': # apply Inverted Method to restructure self.bases_to_write
                            print('Inverted contig found. Adding bases in reverse orientation')
                            self.inverted_transformation(self.bases_to_write)
                self.remapped_genome = "".join([str(contig) for contig in self.remapped_genome]) # self.remapped_genome is a list of all of the pieces of each LG that have been oriented or added into the correct orientation/order. This combines them into one string to write into the genome file
                self.remapped_genome = re.sub("(.{80})", "\\1\n", self.remapped_genome, 0, re.DOTALL)[:-101] # the -102 removes the last 102 characgers which inclues 100 Ns and 2 newlines
                fh.write(self.remapped_genome)
                fh.write('\n')
            print('MAIN LGs REORGANIZED AND WRITTEN INTO THE NEW GENOME FILE...')

    def normal_transformation(self, basepair_list):
        self.bases = basepair_list
        self.bases = self.bases.__str__() + 'N'*100
        self.remapped_genome.append(self.bases)

    def inverted_transformation(self, basepair_list):
        self.bases = basepair_list
        self.bases = self.bases.reverse.complement
        self.bases = self.bases.__str__() + 'N'*100
        self.remapped_genome.append(self.bases)

    def _write_unmapped_contigs(self):
        print('WRITING UNUSED CONTIGS TO THE END OF THE FILE...')
        with open(self.out_dir + '/' + self.filename, 'a') as fh: # unmapped regions are appended to the end of the file instead of written 
            # Code that will get all unused contigs and write them to the end of the genome file
            for contig in self.pyfaidx_obj.keys(): # self.pyfaidx_obj.keys() gives an iterable list of all contig names in UMD2a
                if contig == 'NC_027944.1': # to account for the mt DNA contig
                    self.contigs_unused_list.append(contig)
                elif contig.startswith('NC'): # ensures that any unused main linkage groups will not be written to the file
                    continue
                elif contig not in self.contigs_used_list: # if the contig in the iteration is not in the list of contigs that has been used in the reordering, then add it to the contigs_unused_list
                    self.contigs_unused_list.append(contig)

            for unused_contig in self.contigs_unused_list:
                print(f'WRITING {unused_contig} TO THE FILE')
                self.bases_to_write = self.pyfaidx_obj[unused_contig][:].__str__() # This line uses the pyfaidx_obj and fetches the contig name and all bases in that contig
                fh.write('>' + unused_contig + '\n') # the contig name header is writen
                fh.write(re.sub("(.{80})", "\\1\n", self.bases_to_write, 0, re.DOTALL)) # the bases are split to wrap every 80 basepairs into a newline
                fh.write('\n') # necessary to put each subsequent unused contig into a newline
            # Important Note: The number of used contigs is 152. This means that 130 unmapped contigs were added into the genome. 
            # This leaves 1538 contigs unused. Since the 130 unmapped ones were placed into the existing 22 LGs, the new total number of contigs in Mzebra_GT1 will be 22 + 1538 remaining unused contigs = 1560 total contigs. 
            # This was verified using pyfaidx


    def run_methods(self):
        self._rearrange_genome()
        if args.write_unmapped_contigs:
            self._write_unmapped_contigs()

remap_obj = ReorderGenome(args.genome, args.output_dir[0], args.regions, args.file_name)
remap_obj.run_methods()
print('GENOME REMAPPED')

"""
COMMAND FOR LOCAL TESTING:
/Users/kmnike/anaconda3/envs/pipeline/bin/python3 reorder_umd2a.py Mzebra_UMD2a /Users/kmnike/CichlidSRSequencing/reorder_genome_testing -r All

COMMAND FOR LOCAL TESTING ON LESS LGs:
/Users/kmnike/anaconda3/envs/pipeline/bin/python3 reorder_umd2a.py Mzebra_UMD2a /Users/kmnike/CichlidSRSequencing/reorder_genome_testing -r LG1 LG2


/Users/kmnike/miniforge3/envs/pipeline_x86/bin/python3 reorder_umd2a.py --genome Mzebra_UMD2a -o /Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_GT1 -r All --write_unmapped_contigs -f Mzebra_GT1_v2.1.fna

LG21 Testing:
/Users/kmnike/miniforge3/envs/pipeline_x86/bin/python3 reorder_umd2a.py --genome Mzebra_UMD2a -o /Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/reference_improvement/ -r LG21 -f lg21_test.fna

Code for testing the number of contigs in the new GT1 file  
/Users/kmnike/miniforge3/envs/pipeline_x86/bin/python3
from pyfaidx import Fasta
gt1 = Fasta('/Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_GT1/Mzebra_GT1_v2.fna')
umd2a = Fasta('/Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna')

Index the genome with bwa index:
bwa index <GT.fna>
This one takes a little while.
Upload the genome and all files to Dropbox. Be sure to include a non-zipped version too. 

build an index using samtools faidx <GT.fa>

build dictionary using gatk CreateSequenceDictionary -R <GT.fa>
Be sure to remove the genome files from the git directory before pushing anything or I'll crash the repo!
"""