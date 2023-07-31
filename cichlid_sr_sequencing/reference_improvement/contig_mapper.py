import pdb, argparse, pysam, pathlib, re
import subprocess as sp
import pandas as pd

parser = argparse.ArgumentParser(usage = 'Script will take in a BAM file containing contigs aligned to a MZebra Reference and map the contigs, in order, according to Linkage Groups')
parser.add_argument('-i', '--input', help = 'absolute filepath to the input BAM file')
parser.add_argument('-o', '--output', help= 'absolute filepath to the output file', default = 'mapped_assembly.fasta')
args = parser.parse_args()

"""
TO DO:

"""

class ContigMapper:
    def __init__(self, input_bam, output_fasta):
        self.bam = input_bam
        self.out = output_fasta
        # Check if bam index exists. If not, make it sort the bam file then generate the index. Note that this does edit the bam fiel in place durign the sorting.
        if not pathlib.Path(self.bam + '.bai').exists():
            sp.run(['samtools', 'sort', self.bam, '-o', self.bam])
            sp.run(['samtools', 'index', self.bam])
        self.pysam_obj = pysam.AlignmentFile(self.bam, 'rb')
        self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1', 
                             'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1', 
                             'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1', 
                             'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1'}
        self.reference_names = list(self.linkage_group_map.values())

    def _map_contigs_to_LGs(self):
        bases_to_write = []
        with open(self.out, 'w') as fh: # open the outfile to write
            for lg in self.linkage_group_map.keys(): # start by iterating through each LG
                print('MAPPING CONTIGS TO', lg)
                fh.write('>' + lg + '\n') # write LG name and a newline
                for alignment in self.pysam_obj.fetch(): # per alignment (line) in the bam file,  if the alignment reference starts with the lg name,
                    if alignment.flag in [0,16] and alignment.reference_name.startswith(self.linkage_group_map[lg]):
                        bases_to_write.append(alignment.seq + 'N'*100) # add the sequence for that line into bases_to_write and do it for all lines in the file that start with the lg name
                    else:
                        continue
                bases_to_write = "".join([str(contig) for contig in bases_to_write])[:-len('N' * 100)] # bases_to_write is a list of all of the pieces of each LG that have been combined, in order, into a list This combines them into one string to write into the genome file, then removes the final 100 extra Ns
                bases_to_write = re.sub("(.{80})", "\\1\n", bases_to_write, 0, re.DOTALL) # this wraps every 80 bp into a newline 
                fh.write(bases_to_write) # write the wrapped bases for the right formatting
                fh.write('\n') # insert a newline to signifiy the start of a new LG
                bases_to_write = [] # reset the bases_to_write variable and got through it again for the next LG
                print(lg, " MAPPED")


    def _write_unmapped_contigs(self):
        bases_to_write = []
        print('WRITING UNMAPPED CONTIGS TO THE END OF THE FILE...')
        with open(self.out, 'a') as fh: # unmapped regions are appended to the end of the file instead of written 
            for alignment in self.pysam_obj.fetch():
                if alignment.reference_name not in self.reference_names: # self.reference_names already contains the names of the main linkage groups at this point, so only the non-LGs will be mapped here.
                    if alignment.flag in [0,16]:
                        print('WRITING UNMAPPED CONTIG ', alignment.reference_name)
                        fh.write('>' + alignment.reference_name + '\n') # write the name of the unmapped contig
                        self.reference_names.append(alignment.reference_name) # by adding the reference name to the list of stored references, we will keep looping through the "bases_to_write.append(alignment.seq + 'N'*100)" code...
                    else:
                        continue
                    if len(bases_to_write) > 0: # This code block formats and adds the bases_to_write to the file, only if bases_to_write is not empty. This prevents an extra newline from being added when processing the first unmapped contig 
                        bases_to_write = "".join([str(contig) for contig in bases_to_write])[:-len('N' * 100)] # This combines the list of bases found in bases_to_write into one string to write into the genome file and removes the final 100Ns from the last unmapped contig 
                        bases_to_write = re.sub("(.{80})", "\\1\n", bases_to_write, 0, re.DOTALL) # parses every 80 characters with a newline and removes the final 100bp of NNN nts 
                        fh.write(bases_to_write) # write the formatted bases
                        fh.write('\n') # add a newline before the next linkage group's header 
                        bases_to_write = [] # reset the variable for the next unmapped contig 
                    bases_to_write.append(alignment.seq + 'N'*100) # If the reference the code is processing was not in self.reference_names, then add the bases to a blank bases_to_write list. Once the first bases got added, then the below else statement applies until the code runs into a new reference name. Then the if statement above applies and the code formats and writes the bases
                else:
                    bases_to_write.append(alignment.seq + 'N'*100) # if the reference name already exists in self.reference_names, then we keep adding those bases to bases_to_add

    def run_methods(self):
        self._map_contigs_to_LGs()
        self._write_unmapped_contigs()

contig_mapper_obj = ContigMapper(args.input, args.output)
contig_mapper_obj.run_methods()

print('CONTIGS SUCCESSFULLY MAPPED ')
"""
Command for local testing on subset bam file:
python contig_mapper.py -i /Users/kmnike/Data/CichlidSequencingData/Outputs/alignment/pbmm2/lja_assembly_to_gt1_v2_subset.bam -o subset_mapped_assembly.fasta

Command for local genome building on all lja assembly contigs
python contig_mapper.py -i /Users/kmnike/Data/CichlidSequencingData/Outputs/alignment/pbmm2/sorted_lja_assembly_to_gt1_v2.bam

Command for mapping lja assembly contigs to LGs using GT1_v2 as the reference on the Mzebra Server:
python contig_mapper.py -i /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/alignment/pbmm2/sorted_lja_assembly_to_gt1_v2.bam -o /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/alignment/mapped_fastas/lja_contigs_mapped_to_lgs.fasta
"""