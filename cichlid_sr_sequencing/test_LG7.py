# for testing on the mzebra server only
import subprocess as sp, argparse, shlex, pysam, os
from helper_modules.file_manager import FileManager as FM


parser = argparse.ArgumentParser(usage = 'Testing gatk GenomicsDBImport parallelization optimization')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Make sure genome version is a valid option
if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))

# fm_obj.downloadData(fm_obj.localAlignmentFile)
#### Download the reference genome files and build the proper directry structure for them locally.
# fm_obj.downloadData(fm_obj.localGenomeDir)

lg7 = 'NC_036786.1'
sp.run(shlex.split(f"gatk GenomicsDBImport --genomicsdb-workspace-path {'/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Databases/' + lg7 + '_database'} --intervals 4_per_LG.interval_list --sample-name-map lg7_sample_map.txt --interval-merging-rule OVERLAPPING_ONLY --max-num-intervals-to-import-in-parallel 4 --overwrite-existing-genomicsdb-workspace"))
# sp.run(shlex.split(f"gatk GenotypeGVCFs -R /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {'gendb://../../Data/CichlidSequencingData/Databases/' + lg7 + '_database/'}  -O {'/home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/' + lg7 + '_output.vcf'} --heterozygosity 0.0012"))

print('Pipeline Completed')