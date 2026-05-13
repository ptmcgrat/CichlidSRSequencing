import argparse,os,sys,pysam,subprocess
sys.path.append("..") # Adds higher directory to python modules path.
from helper_modules.file_manager import FileManager as FM
from helper_modules.smallVariantGenotyper import genotype_at_sites
from helper_modules.InsertionGenotyper import IndelReadClassifier, ClassifierCall, write_classifier_vcf, concat_sample_vcfs
import pandas as pd

# Need to make SampleIDs and ProjectIDs mutually exclusive
parser = argparse.ArgumentParser(usage = 'This script will split bamfiles ')
parser.add_argument('SV_VCF', type = str, help = 'VCF of small variants to genotype')
parser.add_argument('LV_CSV', type = str, help = 'CSV of large variants to genotype')
parser.add_argument('OUT_VCF', type = str, help = 'CSV of large variants to genotype')
parser.add_argument('genome_version', type = str, help = 'Version of the genome to align to')
parser.add_argument('SampleID', type = str, help = 'Version of the genome to align to')
args = parser.parse_args()

fm_obj = FM(genome_version = args.genome_version)
fm_obj.createSampleFiles(args.SampleID, reads = False)
refObj = pysam.FastaFile(fm_obj.localGenomeFile)

fm_obj.downloadData(fm_obj.localSampleBamDir)

sv_temp_vcf = fm_obj.localSampleTempDir + args.SampleID + '.sv.vcf.gz'
genotype_at_sites(fm_obj.localBamFile, args.SV_VCF, fm_obj.localGenomeFile, sv_temp_vcf)

dt = pd.read_csv(args.LV_CSV)

lv_temp_vcf = fm_obj.localSampleTempDir + args.SampleID + '.lv.vcf.gz'

calls = []
for i,row in dt.iterrows():
    classifier = IndelReadClassifier.from_vcf_record(
        ref_genome=refObj, chrom=row.Chromosome, position=row.Position,
        ref=row.Reference.upper(), alt=row.Alt.upper(), flanking=250,
    )

    # Pull reads from your aligned BAM (and optionally a discordant/unmapped BAM)
    reads = classifier.fetch_reads_near_insertion(
        bam=fm_obj.localBamFile,
        discordant_bam=fm_obj.localDiscordantBamFile,   # optional
    )

    pair_results = classifier.classify_read_pairs(reads)
    call = classifier.call_genotype(pair_results)
    calls.append(ClassifierCall(
        chrom=row.Chromosome, pos=row.Position, ref=row.Reference, alt=row.Alt,
            gt=call.genotype, ad_ref=call.n_ref, ad_alt=call.n_alt,
    ))
 
write_classifier_vcf(calls, args.SampleID, lv_temp_vcf, template_vcf=sv_temp_vcf)
concat_sample_vcfs(lv_temp_vcf, lv_temp_vcf, args.OUT_VCF)

fm_obj.uploadData(args.OUT_VCF)
subprocess.run(['rm','-rf',fm_obj.localSampleBamDir])
