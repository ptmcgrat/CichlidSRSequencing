import subprocess as sp, argparse, shlex, pandas as pd, pysam, os
from helper_modules.file_manager import FileManager as FM


parser = argparse.ArgumentParser(usage = 'Testing gatk GenomicsDBImport parallelization optimization')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Make sure genome version is a valid option
if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))

# Get the error files
df = pd.read_csv('returncodes.txt', sep='\t', header=None, names=['col1','col2'])
error_files=df['col1'].tolist()
print(error_files)
#### Just to get this running on the utaka server:
for file in error_files:
	sp.run(shlex.split(f"{'rclone copy ptm_dropbox:BioSci-McGrath/Apps/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/'+file+'/'+file+'all.bam'} {'/Data/mcgrath-lab/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/'+file+'/'}"))




# lg7 = 'NC_036786.1'
# processes = []

# local_files = ['MC_1_m', 'SAMEA2661294', 'SAMEA2661322', 'SAMEA4032100', 'SAMEA4033261']
# for file in local_files:
# 	if os.path.isfile(f"/Users/kmnike/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/{file}.all.bam"):
# 		os.remove(f"/Users/kmnike/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/{file}.all.bam")

# for file in error_files:
# 	if os.path.isfile(f"/Data/mcgrath-lab/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/{file}.all.bam"):
# 		os.remove(f"/Data/mcgrath-lab/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/{file}.all.bam")


### parallel download test on mzebra
# for file in error_files:lab
# 	p = sp.Popen(shlex.split(f"rclone copy ptm_dropbox:BioSci-McGrath/Apps/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/{file}.all.bam /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/ -P"))
# 	processes.append(p)
# 	if len(processes) == 74:


"""
for file in error_files:
	# sp.run(shlex.split(f"rclone copy ptm_dropbox:BioSci-McGrath/Apps/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/{file}.all.bam /Data/mcgrath-lab/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/ -P"))
	sp.run(shlex.split(f"rclone copy ptm_dropbox:BioSci-McGrath/Apps/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/{file}.all.bai /Data/mcgrath-lab/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/ -P"))
"""
"""
for file in error_files:
	command = f"gatk --java-options '-Xmx16G' HaplotypeCaller --emit-ref-confidence GVCF  -L {lg7} -R /Data/mcgrath-lab/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -I {'/Data/mcgrath-lab/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/' + file + '.all.bam'} -O {'/Data/mcgrath-lab/Data/' + file + '/' + file + '_rerun_230203.g.vcf.gz'}"
	p = sp.Popen(shlex.split(command))
	processes.append(p)
	if len(processes) == 22:
		for proc in processes:
			proc.communicate()
		proc = []
"""


#### for local testing:


# fm_obj.downloadData(fm_obj.localAlignmentFile)
#### Download the reference genome files and build the proper directry structure for them locally.
# fm_obj.downloadData(fm_obj.localGenomeDir)

# a_dt = pd.read_csv(fm_obj.localAlignmentFile)
# a_dt = a_dt[a_dt.GenomeVersion == args.Genome]
# sampleIDs = set(a_dt.SampleID)
