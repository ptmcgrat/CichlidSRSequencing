import argparse, pdb, subprocess as sp, pysam
from helper_modules.file_manager import FileManager as FM
from helper_modules.Timer import Timer
import shlex
import pandas as pd

parser = argparse.ArgumentParser(usage = 'This script is for running unmapped contigs in Mzebra_UMD2a')
# parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
genome = 'Mzebra_UMD2a'
fm_obj = FM(genome)

# Make sure genome version is a valid option

a_dt = pd.read_csv(fm_obj.localAlignmentFile)
a_dt = a_dt[a_dt.GenomeVersion == genome]
sampleIDs = set(a_dt.SampleID)

gvcffiles = []

# Download the reference genome files and build the proper directry structure for them locally.
# fm_obj.downloadData(fm_obj.localGenomeDir)
fasta_obj = pysam.FastaFile(fm_obj.localGenomeFile)
contigs = fasta_obj.references[22:]

processes1 = []
processes2 = []

#### Process unmapped contigs 22 at a time through GenomicsDBImport. I'm giving each process 20Gb --> 440Gb memory being used per 22 processes
for contig in contigs:
	p1 = sp.Popen(shlex.split(f"gatk --java-options '-Xmx20G' GenomicsDBImport --genomicsdb-workspace-path {'/Data/mcgrath-lab/Data/CichlidSequencingData/Databases/' + contig + '_database'} --intervals {'/home/mcgrath-lab/nkumar317/CichlidSRSequencing/cichlid_sr_sequencing/intervals_unmapped_contigs/' + contig + '.interval_list'} --sample-name-map sample_map_utaka.txt --interval-merging-rule OVERLAPPING_ONLY --max-num-intervals-to-import-in-parallel 4 --overwrite-existing-genomicsdb-workspace"))
	processes1.append(p1)
	for proc in processes1:
		if len(processes1) == 22:
			proc.communicate()
	processes1 = []


for contig in contigs:
	p2 = sp.Popen(shlex.split(f"gatk --java-options '-Xmx20G' GenotypeGVCFs -R /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {'gendb://../../../../../Data/mcgrath-lab/Data/CichlidSequencingData/Databases/' + contig + '_database'}  -O {'/Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/' + contig + '_output'} --heterozygosity 0.0012"))
	processes2.append(p2)
	for proc in processes2:
		if len(processes2) == 22:
			proc.communicate()
	processes2 = []
	
print('Pipeline Completed')
