import argparse, pdb, subprocess as sp, pysam
from helper_modules.file_manager import FileManager as FM
from helper_modules.Timer import Timer
import shlex
import pandas as pd

parser = argparse.ArgumentParser(usage = 'This script is for troubleshooting why LGs 3, 7, 15, and 20 failed to run ')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Make sure genome version is a valid option
if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))

fm_obj.downloadData(fm_obj.localAlignmentFile)

a_dt = pd.read_csv(fm_obj.localAlignmentFile)

a_dt = a_dt[a_dt.GenomeVersion == args.Genome]

sampleIDs = set(a_dt.SampleID)

gvcffiles = []

# Download the reference genome files and build the proper directry structure for them locally.
# fm_obj.downloadData(fm_obj.localGenomeDir)
fasta_obj = pysam.FastaFile(fm_obj.localGenomeFile)

# Download the GVCF and GVCF.idx files for each sample in the AlignmentDatabase
# for sampleID in sampleIDs:
# 	if sampleID in []:
# 		continue
# 	print(sampleID)
# 	fm_obj.createSampleFiles(sampleID)
# 	fm_obj.downloadData(fm_obj.localGVCFFile)
# 	fm_obj.downloadData(fm_obj.localGVCFFile + '.tbi')
# 	gvcffiles.append(fm_obj.localGVCFFile)
	# pdb.set_trace() # used for troubleshooting

processes = []
processes2 = []

#### Since the intervals list will generate one database for all of the unmapped contigs, I can run GenomcisDBImport and GenotypeGVCFs back to back without the need to parallelize further (since it's not possible to parallelize GenotypeGVCFs this way anyways)
sp.run(shlex.split(f"gatk GenomicsDBImport --genomicsdb-workspace-path {'/Data/mcgrath-lab/Data/CichlidSequencingData/Databases/' + 'all_unmapped_contigs_database'} --intervals unmapped_contigs.interval_list --sample-name-map sample_map_utaka.txt --interval-merging-rule OVERLAPPING_ONLY --max-num-intervals-to-import-in-parallel 4 --overwrite-existing-genomicsdb-workspace"))
sp.run(shlex.split(f"gatk GenotypeGVCFs -R /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {'gendb://../../../../../Data/mcgrath-lab/Data/CichlidSequencingData/Databases/all_unmapped_contigs_database'}  -O {'/Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/' + 'all_unmapped_contigs_output.vcf'} --heterozygosity 0.0012"))

print('Pipeline Completed')