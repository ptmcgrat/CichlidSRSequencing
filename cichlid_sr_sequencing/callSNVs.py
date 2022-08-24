import argparse, pdb, timer
from helper_modules.file_manager import FileManager as FM
import pandas as pd



# Need to make SampleIDs and ProjectIDs mutually exclusive
parser = argparse.ArgumentParser(usage = 'This script will download fastq data the McGrath lab dropbox and align it to the Genome version of choice')
parser.add_argument('Genome', type = str, help = 'Version of the genome to align to')
parser.add_argument('-AnalysisID', type = str, help = 'Restrict analysis to a specific ProjectID. Project are akin to BioProject and contain multiple Samples')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj = FM(args.Genome)

# Make sure genome version is a valid option
if args.Genome not in fm_obj.returnGenomeVersions():
	raise argparse.ArgumentTypeError('Genome version does not exist. Options are: ' + ','.join(fm_obj.returnGenomeVersions()))

# Download master sample database and read it in (this contains info on valid sampleIDs, projectIDs, and file locations)
fm_obj.downloadData(fm_obj.localAlignmentFile)
a_dt = pd.read_csv(fm_obj.localAlignmentFile)

pdb.set_trace()

# If running on projectID, make sure it is valid and subset sample database to those with the right projectID
if args.ProjectID is not None:
	if args.ProjectID not in set(a_dt.ProjectID):
		raise argparse.ArgumentTypeError('ProjectID ' + args.ProjectID + ' does not exist. Options are: ' + ','.join(set(a_dt.ProjectID)))
	a_dt = a_dt[a_dt.ProjectID == args.ProjectID]
	good_samples = set(a_dt.SampleID)

# If running on sampleIDs, make sure they are valid
if args.SampleIDs is not None:
	bad_samples = []
	for sample in args.SampleIDs:
		if sample not in list(a_dt.SampleID):
			bad_samples.append(sample)

	if len(bad_samples) > 0:
		raise argparse.ArgumentTypeError('The following samples were not found in sample database: ' + ','.join(bad_samples))
	good_samples = set(args.SampleIDs)

# Make directories necessary for analysis
os.makedirs(fm_obj.localMasterDir, exist_ok = True)
os.makedirs(fm_obj.localTempDir, exist_ok = True)

# Download genome data necessary for analysis
timer.start('Downloading genome')		
fm_obj.downloadData(fm_obj.localGenomeDir)
timer.stop()

for sample in good_samples:
	fm_obj.createBamFiles(sample)
	fm_obj.createPileupFiles(sample)



	command = ['gatk', 'HaplotypeCaller', 'I', fm_obj.localBamFile, '-O', fm_obj.localGVCFFile, '-R', fm_obj.localGenomeFile, '--emit-ref-confidence', 'GVCF']

java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R REFERENCE.fa --
emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter
128000 –I SAMPLEn.bam –o GATK_SAMPLEn.g.vcf
3
Haplotype caller per-sample files were combined using the GATK GenotypeGVCFs tool with
the --includeNonVariantSites option so that every basepair of the assembly was
represented in the multisample VCF file.
Hard filters

Hard filters were applied to the following GATK annotations:
Minimal inbreeding coefficient: 'InbreedingCoeff < -0.6'
Minimum overall read depth: 'DP < 1000'
Maximum overall read depth: 'DP > 3000' (except for mtDNA: scaffolds
747,2036)
Max phred-scaled p-value from Fisher's exact test to detect strand bias: 'FS
> 40.0' (except for mtDNA: scaffolds 747,2036)
QualityByDepth: 'QD < 2.0’
Excess Missingness: 'NCC > 32' (>16 individuals with missing data)
More than 10% of reads have mapping quality zero: '(MQ0/(1.0*DP)) > 0.10'
Low mapping quality: 'MQ < 40.0'

samtools calling (multisample):
samtools mpileup -t DP,DPR,INFO/DPR -C50 -pm2 -F0.2 –ugf REFERENCE.fa
SAMPLE1.bam SAMPLE2.bam ... | bcftools call -vmO z -f GQ -o
samtools_VARIANTS.vcf.gz
The consensus GATK and samtools call set was obtained using the GATK:
java -Xmx10000m -jar GenomeAnalysisTK.jar -T SelectVariants -R REFERENCE.fa -
-variant onlyVariants_filtered.vcf.gz -o onlyVariants_filtered_concord.vcf.gz
--concordance samtools_unfiltered.vcf.gz
"""