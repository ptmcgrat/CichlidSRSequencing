import argparse, pdb, re
import pandas as  pd
from pyfaidx import Fasta
from helper_modules.nikesh_file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage = 'This script will reorder the nucleotides in a reference genome as informed by OGM Data')
parser.add_argument('genome', help = 'absolute path to UMD2a reference genome')
parser.add_argument('output_dir', help = 'absolute filepath to an output directory')
parser.add_argument('-r', '--regions', help = 'list of linakge groups where script will operate. All linkage groups will be output into the same directory', nargs = '*', default = ['All'])
args = parser.parse_args()

"""
Link to mapping google doc:
https://docs.google.com/spreadsheets/d/e/2PACX-1vRLiFIHXSXrljCLzO58qCwqcvdU3sw-lNAA8P_riJcCS14mSQRcNjZMRnsW3eKrAJcZENUbpEJoiakH/pub?output=csv

TODO:
Nothing for now...
"""

class ReorderGenome:
    def __init__(self, genome, output_directory, linkage_groups):
        # The linkage_group_map attribute is a hard coded list of LG names
        self.fm_obj = FM(genome)
        self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1', 
                             'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1', 
                             'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1', 
                             'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1'}
        self.genome = self.fm_obj.localGenomeFile
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

    def make_remap_instructions(self, linkage_group):
        url = f'https://docs.google.com/spreadsheets/d/e/2PACX-1vRLiFIHXSXrljCLzO58qCwqcvdU3sw-lNAA8P_riJcCS14mSQRcNjZMRnsW3eKrAJcZENUbpEJoiakH/pub?output=xlsx'
        self.df = pd.read_excel(url, sheet_name = linkage_group)
        remapping_instructions = []
        for index, rows in self.df.iterrows():
            remapping_instructions.append([rows.contig_name, rows.Start, rows.Stop, rows.Direction])
        return(remapping_instructions)

    def _rearrange_genome(self):
        self.pyfaidx_obj = Fasta(self.genome) # load in the reference genome as a Fasta class object from pyfaidx
        with open(self.out_dir + '/Mzebra_GT1.fna', 'w') as fh:
            for lg in self.linkage_groups:
                print('Getting Instructions for how to rebuild ' + lg + '...')
                linkage_group_instructions = self.make_remap_instructions(lg)
                fh.write('>' + self.linkage_group_map[lg] + '\n')
                self.remapped_genome = []
                for region in linkage_group_instructions:
                    if region[0] == 'Gap':
                        print('Gap detected in instructions... Adding 200 Ns to the reference to indicate presence of a gap')
                        self.remapped_genome.append('N'*200)
                    else:
                        self.bases_to_write = self.pyfaidx_obj[region[0]][region[1]:region[2]]
                        if region[3] == "Normal":
                            print('Normal contig found. Adding bases in the same orientation.')
                            self.normal_transformation(self.bases_to_write)
                        elif region[3] == 'Inversion': # apply Inverted Method to restructure self.bases_to_write
                            print('Inverted contig found. Adding bases in reverse orientation')
                            self.inverted_transformation(self.bases_to_write)
                self.remapped_genome.pop() # This line removes the extra Ns written after the final contig in the linkage group before writing the bases into the file.
                self.remapped_genome = "".join([str(contig) for contig in self.remapped_genome])
                self.remapped_genome = re.sub("(.{80})", "\\1\n", self.remapped_genome, 0, re.DOTALL)[:-102] # the -102 removes the last 102 characgers which inclues 100 Ns and 2 newlines
                fh.write(self.remapped_genome)
                fh.write('\n')

    def normal_transformation(self, basepair_list):
        self.bases = basepair_list
        self.bases = self.bases.__str__() + 'N'*100
        self.remapped_genome.append(self.bases)

    def inverted_transformation(self, basepair_list):
        self.bases = basepair_list
        self.bases = self.bases.reverse
        self.bases = self.bases.__str__() + 'N'*100
        self.remapped_genome.append(self.bases)

    def run_methods(self):
        self._rearrange_genome()

remap_obj = ReorderGenome(args.genome, args.output_dir, args.regions)
remap_obj.run_methods()
print('Genome Remapped')

"""
COMMAND FOR LOCAL TESTING:
/Users/kmnike/anaconda3/envs/pipeline/bin/python3 reorder_umd2a.py Mzebra_UMD2a /Users/kmnike/CichlidSRSequencing/reorder_genome_testing -r All

COMMAND FOR LOCAL TESTING ON LESS LGs:
/Users/kmnike/anaconda3/envs/pipeline/bin/python3 reorder_umd2a.py Mzebra_UMD2a /Users/kmnike/CichlidSRSequencing/reorder_genome_testing -r LG1 LG2

"""