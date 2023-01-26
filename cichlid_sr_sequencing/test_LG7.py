# for testing on the mzebra server only
import subprocess as sp, argparse, shlex, pandas as pd, pysam
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

# a_dt = pd.read_csv(fm_obj.localAlignmentFile)
# a_dt = a_dt[a_dt.GenomeVersion == args.Genome]
# sampleIDs = set(a_dt.SampleID)

data = pd.read_csv('sample_IDs.csv')
sampleIDs = data['sample_IDs'].tolist()

#### Confirm absence of SAMEA4033252 and that only SAMEA4033252_redo will be run:
# for i in sampleIDs:
# 	if "redo" in i:
# 		print(i)
# 	if i == 'SAMEA4033252':
# 		print(i)

lg7 = 'NC_036786.1'
#### After a discussion with Patrick, the most efficient and quick way to figure out which samples are not running for LG7 is to simply run each sample one by 1 and store its data into a separate Database. If a database's creation is not working, then on exception, we can record the sample name into a file and pass that one
# processes = []
# for sample in sampleIDs:
# 	try:
# 		p = sp.run(shlex.split(f"gatk --java-options '-Xmx450G' GenomicsDBImport --genomicsdb-workspace-path {'/Data/mcgrath-lab/Data/CichlidSequencingData/TestingDatabases/' + sample + '_database'} --intervals lg7.interval_list -V {'/Data/mcgrath-lab/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/' + sample + '/' + sample + '.g.vcf.gz'} --interval-merging-rule OVERLAPPING_ONLY --max-num-intervals-to-import-in-parallel 4 --overwrite-existing-genomicsdb-workspace"))
# 		processes.append(p)
# 		if len(processes) == 22:
# 			for p in processes:
# 				p.communicate()
# 			processes = []
# 	except:
# 		with open('failed_samples_230126.txt', 'a+') as fh:
# 			fh.write(f"sample {sample} failed GenomicsDBImport\n")
# 		pass

#### The above code was used to ID the failed sample (SAMEA4033252). The below code is to create an LG7 database for the samples including the SAMEA4033252 and create a VCF file once the GenomicsDBImport is complete
sp.run(shlex.split(f"gatk --java-options '-Xmx450G' GenomicsDBImport --genomicsdb-workspace-path {'/Data/mcgrath-lab/Data/CichlidSequencingData/Databases/' + lg7 + '_database'} --intervals lg7.interval_list --sample-name-map sample_map_utaka.txt --interval-merging-rule OVERLAPPING_ONLY --max-num-intervals-to-import-in-parallel 4 --overwrite-existing-genomicsdb-workspace"))
sp.run(shlex.split(f"gatk --java-options '-Xmx450G' GenotypeGVCFs -R /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {'gendb://../../../../../Data/mcgrath-lab/Data/CichlidSequencingData/Databases/' + lg7 + '_database/'}  -O {'/Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/' + lg7 + '_output.vcf'} --heterozygosity 0.0012"))

print('Pipeline Completed')