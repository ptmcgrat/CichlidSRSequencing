import subprocess as sp, argparse, shlex
from helper_modules.file_manager import FileManager as FM


parser = argparse.ArgumentParser(usage = 'Testing gatk GenomicsDBImport parallelization optimization')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Make sure genome version is a valid option
if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))

# Get the error files; for some reason, the .tolist() funciton is not outputting a list that rclone likes, so I just printed the list and hard coded it in for now. 
error_files = ['SAMEA4033321', 'SAMEA2661241', 'SAMEA4032067', 'SAMEA4033318', 'SAMEA4032108', 'SAMEA4033341', 'SAMEA4033314', 'SAMEA4032129', 'SAMEA4033252', 'SAMEA4032105', 'SAMEA3388855', 'SAMEA4032131', 'SAMEA4033248', 'SAMEA2661359', 'SAMEA1877414', 'SAMN06685761', 'SAMEA4032096', 'SAMEA2661294', 'SAMEA4032053', 'SAMEA1920091', 'SAMEA4032046', 'SAMEA4033254', 'SAMEA2661396', 'SAMEA4033283', 'SAMEA4032136', 'SAMEA2661367', 'SAMEA3388862', 'SAMEA4033301', 'SAMEA2661292', 'SAMEA4033271', 'SAMEA4032127', 'SAMEA2661414', 'SAMEA4033317', 'SAMEA2661221', 'SAMEA2661310', 'SAMEA1904322', 'SAMEA4032051', 'SAMEA4032033', 'SAMEA4033277', 'SAMEA4033307', 'SAMEA4032038', 'SAMEA4032042', 'SAMEA2661239', 'SAMN08051112', 'SAMEA4032125', 'SAMEA2661258', 'SAMEA2661306', 'SAMEA4033322', 'SAMEA1920095', 'SAMEA4033284', 'SAMEA4033304', 'SAMEA4033305', 'SAMEA4032137', 'SAMEA1920096', 'SAMEA3388871', 'SAMEA4033278', 'SAMEA4033286', 'SAMEA4032049', 'SAMEA1920092', 'SAMEA4033289', 'SAMEA4033324', 'SAMEA1904832', 'SAMEA2661339', 'SAMEA1877499', 'SAMEA2661277', 'SAMEA2661248', 'SAMEA2661280', 'SAMEA2661381', 'SAMEA1904328', 'SAMEA4033261', 'SAMEA4033287', 'SAMEA4032088', 'SAMEA4032070', 'SAMEA4032069']
#### Just to get this running on the utaka server:
# for file in error_files:
	# sp.run(shlex.split(f"rclone copy ptm_dropbox:BioSci-McGrath/Apps/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/{file}.all.bam /Data/mcgrath-lab/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/ -P"))
	# sp.run(shlex.split(f"rclone copy ptm_dropbox:BioSci-McGrath/Apps/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/{file}.all.bai /Data/mcgrath-lab/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/ -P"))


lg7 = 'NC_036786.1'
processes = []


#### parallel download test on mzebra
for file in error_files:
	p = sp.Popen(shlex.split(f"rclone copy ptm_dropbox:BioSci-McGrath/Apps/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/{file}.all.bam /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/ -P"))
	processes.append(p)
	if len(processes) == 2:
		for proc in processes:
			proc.communicate()
		processes = []


"""
for file in error_files:
	# sp.run(shlex.split(f"rclone copy ptm_dropbox:BioSci-McGrath/Apps/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/{file}.all.bam /Data/mcgrath-lab/Data/CichlidSequencingData/Bamfiles/Mzebra_UMD2a/{file}/ -P"))
	
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
