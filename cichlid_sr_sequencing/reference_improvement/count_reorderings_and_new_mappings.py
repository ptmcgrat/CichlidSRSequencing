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
Script was used on 24.02.12 to calculate how many contigs of UMD2a were reordered and how many unmapped contigs were placed into the the GT1, then GT2 genomes  
"""

class CountStuff:
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
    

    def _make_stats(self):
        with open('stats.txt', 'w') as fh:
            reorder_running_total = 0
            new_contig_running_total = 0
            for lg in self.linkage_groups:
                print('working on ' + lg)
                fh.write(lg + '\n')
                reorder_counter = 0
                new_contig_placed = 0
                linkage_group_instructions = self.make_remap_instructions(lg)
                for region in linkage_group_instructions:
                    if region[3] != 'Normal' and region[3] != 'Gap':
                        reorder_counter += 1
                        reorder_running_total += 1
                    if region[0] != self.linkage_group_map[lg] and region[0] != 'Gap':
                        new_contig_placed += 1
                        new_contig_running_total += 1
                fh.write('reorder count: ' + str(reorder_counter) + '\n')
                fh.write('new contigs placed: ' + str(new_contig_placed) + '\n')
            fh.write('Total Reordered Contigs: ' + str(reorder_running_total) + '\n')
            fh.write('Total Added  Unmapped Contigs: ' + str(new_contig_running_total) + '\n')

    def run_methods(self):
        self._make_stats()

remap_obj = CountStuff(args.genome, args.output_dir[0], args.regions, args.file_name)
remap_obj.run_methods()
print('STATS CALCULATED')
