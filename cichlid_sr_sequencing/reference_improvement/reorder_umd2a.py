import argparse, os, pdb, pathlib
import pandas as  pd
from pyfaidx import Fasta
import gzip

parser = argparse.ArgumentParser(usage = 'This script will reorder the nucleotides in a reference genome as informed by OGM Data')
parser.add_argument('genome', help = 'absolute path to referecne genome')
parser.add_argument('output_dir', help = 'absolute filepath to an output directory')
parser.add_argument('reorder_map', help = 'file containing the reordering instructions for the reference genome')
parser.add_argument('-r', '--regions', help = 'list of linakge groups where script will operate. All linkage groups will be output into the same directory', nargs = '*')
args = parser.parse_args()

"""
To do:
1. read in the map file and extract the information needed for reordering the genome. 
    - Feilds Needed:
        a. start coordinate
        b. end coordinate
        c. attribute
        d. Ns needed for padding between missing regions
    - How to encode the reordering data:
        a. list of tuples should work fine. They accept both int & string values.
        b. value locations in tuples will be hard coded in, but this should be ok... 

2. read in genome file's linkage group based on the region
    - based on attribute, write genomic coordinates into a new file
        - start this file with a new header reflecting that this is a reorganized LG (i.e don't use same name as the UMD2a genome's headers)
    - each attribute can have it's own hidden method for how to write the genome reads
"""

class ReorderGenome:
    def __init__(self, genome, output_directory, map_file, linkage_groups):
        # The linkage_group_map attribute is a hard coded list of LG names
        self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1', 
                             'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1', 
                             'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1', 
                             'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1'}
        self.genome = genome
        self.map = map_file
        self.out_dir = output_directory

        # code to take in simplified region inputs like "LG1", "LG7", etc, and convert them to the actual names of contigs, as well as check for invalid inputs & duplciate LG names.
        self.linkage_groups = linkage_groups
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

    def _read_map(self):
        # This method will take in the map file and create a pandas dataframe which we can use to build either a dictionary or a list of lists that will encode how each region will be remapped.
        for lg in self.linkage_groups:
            self.df = pd.read_excel(self.map, sheet_name=lg)
            self.remapping_dict = {}
            for index, rows in self.df.iterrows():
                row_map_data = [rows.Start, rows.Stop, rows.sample_based_region_order, rows.region_attributes, rows.approx_missing_len, rows.anchor]
                self.remapping_dict[row_map_data[2]] = row_map_data[:2] + row_map_data[3:]

    def _rearrange_genome(self):
        self.pyfaidx_obj = Fasta(self.genome) # load in the reference genome as a Fasta class object from pyfaidx
        self.remapped_genome = [] # define the list that will be used to write the genome.

        """
        Testing methods on genome subset
        # self.small_data = {3: [0, 613372, 'Normal', 0, 1], 0: [613472, 1343998, 'Inverted', 0, 1], 2: [1344098, 2170558, 'Inverted', 0, 1], 1: [2170658, 4254576, 'Normal', 0, 1], 4: [4254676, 4318228, 'Missing', 700000, 0], 5: [17493698, 17534073, 'Translocation', 0, 362]}
        # for lg in self.linkage_groups:
        #     for region in range(0,len(self.small_data)):
        #         self.bases_to_write = self.pyfaidx_obj[lg][self.small_data[region][0]:self.small_data[region][1]]
        #         # apply transformations to the bases_to_write:
        #         if self.small_data[region][2] == 'Normal': # apply Normal Method to restructure self.bases_to_write
        #             self.normal_transformation(self.bases_to_write)
        #         elif self.small_data[region][2] == 'Inverted': # apply Inverted Method to restructure self.bases_to_write
        #             self.inverted_transformation(self.bases_to_write)
        #         elif self.small_data[region][2] == 'Missing': # apply Missing Method to restructure self.bases_to_write
        #             self.missing_transformation(self.bases_to_write, region, self.small_data[region][3])
        #         elif self.small_data[region][2] == 'Translocation': # apply Translocation Method to restructure self.bases_to_write; for now region will be appended until I figure out how to better integrate these specifically into the ref. 
        #             self.translocation_transformation(self.bases_to_write, self.small_data[region][4])
        #         elif self.small_data[region][2] == 'Irregular': # apply Irregular Method to restructure self.bases_to_write
        #             self.irregular_transformation(self.bases_to_write)
        """

        for lg in self.linkage_groups:
            for region in range(0,len(self.remapping_dict)):
                self.bases_to_write = self.pyfaidx_obj[lg][self.remapping_dict[region][0]:self.remapping_dict[region][1]]
                # apply transformations to the bases_to_write:
                if self.remapping_dict[region][2] == 'Normal': # apply Normal Method to restructure self.bases_to_write
                    self.normal_transformation(self.bases_to_write)
                elif self.remapping_dict[region][2] == 'Inverted': # apply Inverted Method to restructure self.bases_to_write
                    self.inverted_transformation(self.bases_to_write)
                elif self.remapping_dict[region][2] == 'Missing': # apply Missing Method to restructure self.bases_to_write
                    self.missing_transformation(self.bases_to_write, region, self.remapping_dict[region][3])
                elif self.remapping_dict[region][2] == 'Translocation': # apply Translocation Method to restructure self.bases_to_write
                    self.translocation_transformation(self.bases_to_write, self.remapping_dict[region][4])
                elif self.remapping_dict[region][2] == 'Irregular': # apply Irregular Method to restructure self.bases_to_write
                    self.irregular_transformation(self.bases_to_write)

            with open(self.out_dir + '/remapped_genome.fna', 'w') as fh:
                fh.write('>' + lg + '\n')
                for bases in self.remapped_genome:
                    fh.write(bases)

    def normal_transformation(self, basepair_list):
        self.bases = basepair_list
        self.remapped_genome.append(self.bases.__str__())

    def inverted_transformation(self, basepair_list):
        self.bases = basepair_list
        self.bases = self.bases.reverse
        self.remapped_genome.append(self.bases.__str__())

    def missing_transformation(self, basepair_list, region, Ns_to_add):
        self.bases = basepair_list
        self.region = region
        self.Ns = Ns_to_add
        N_list = ['N'] * int(self.Ns)
        self.remapped_genome.append(self.bases.__str__() + ''.join(N_list))

    def translocation_transformation(self, basepair_list, translocated_region):
        self.bases = basepair_list
        self.extra_contig = translocated_region
        self.remapped_genome.append(self.bases.__str__() + self.pyfaidx_obj[self.extra_contig][::].__str__())

    def irregular_transformation(self, basepair_list):
        self.bases = basepair_list
        self.remapped_genome.append(self.bases.__str__())

    def run_methods(self):
        self._read_map()
        self._rearrange_genome()


remap_obj = ReorderGenome(args.genome, args.output_dir, args.reorder_map, args.regions)
remap_obj.run_methods()
print('Genome Remapped')

"""
COMMAND FOR LOCAL TESTING:
/Users/kmnike/anaconda3/envs/pipeline/bin/python3 reorder_umd2a.py /Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz /Users/kmnike/CichlidSRSequencing/reorder_genome_testing /Users/kmnike/Desktop/UMD2a_remap_v1.xlsx -r LG1


"""