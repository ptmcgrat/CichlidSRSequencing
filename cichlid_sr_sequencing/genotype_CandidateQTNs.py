from helper_modules.file_manager import FileManager as FM
from helper_modules.CandidateGenotyper import CandidateGenotyper as CG
import argparse, pdb

fm_obj = FM()

parser = argparse.ArgumentParser(usage = 'This script will download fastq data the McGrath lab dropbox and align it to the Genome version of choice. It will also create gvcf files')
parser.add_argument('Genome', type = str, choices = fm_obj.returnOptions('Genomes'), help = 'Version of the genome to align to')
parser.add_argument('Candidate_TSV', type = str, help = 'TSV file of candidate QTNs')
parser.add_argument('-e', '--Ecogroups', nargs = '+', metavar = '', choices = fm_obj.returnOptions('Ecogroups'), help = 'Restrict analysis to a specific Ecogroup: ' + ','.join(fm_obj.returnOptions('Ecogroups')))
args = parser.parse_args()

cg_obj = CG(args.Genome, args.Candidate_TSV, args.Ecogroups)
cb_obj.genotypeSamples()

