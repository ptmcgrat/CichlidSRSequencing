#!/Users/kmnike/anaconda3/envs/mcgrath/bin/python3
#import argparse, pysam, shlex, subprocess as sp

import argparse, pdb, os, subprocess, pathlib
import pandas as pd
from cyvcf2 import VCF

parser = argparse.ArgumentParser(usage = "This pipeline is for running pca analysis on a filtered vcf file")
parser.add_argument('input_vcffile', help = 'absolute filepath to the filtered, gzipped input file')
parser.add_argument('output_dir', help = 'absolute filepath to an output directory')
parser.add_argument('sample_database', help = 'sample database that lists ecotype for each sample')
parser.add_argument('-e', '--ecogroups', help = 'one or multiple eco group names for filtering the data', choices = ['Mbuna', 'Utaka', 'Shallow Benthic', 'Deep Benthic','Rhampochromis', 'Diplotaxodon', 'Riverine', 'AC', 'Non_Riverine', 'All'], nargs = '*', default = ['All'])
parser.add_argument('-r', '--regions', help = 'list of linkage groups for which analyses will run', nargs = '*', default = 'All')
parser.add_argument('--PCA', help = 'generate a PCA analysis for the specified linkage groups', default = "All")
parser.add_argument('-f', '--filters', help = 'list of tunable parametrs for filtering the raw input vcf file.', default  = ['DP > 11000', 'DP < 8000', 'InbreedingCoeff < -0.6', 'FS > 40.0', 'QD < 2.0', 'NCC> 125', 'MQ < 50', 'AF < 0.000958'])
args = parser.parse_args()

class PCA_Maker:
    def __init__(self, input_vcffile, output_directory, sample_database, ecogroups):
        self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1', 
                             'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1', 
                             'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1', 
                             'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1'}
        self.in_vcf = input_vcffile
        self.out_dir = output_directory
        self.sample_database = sample_database
        self.ecogroups = ecogroups

        self._create_sample_filter_file()

        self.vcf_obj = VCF(self.in_vcf)
        self.contigs = self.vcf_obj.seqnames
        # Makes sure LGs are in vcf file provided
        for lg,contig in self.linkage_group_map.items():
            assert contig in self.contigs
        # Ensure index file exists
        assert os.path.exists(self.in_vcf + '.tbi')
        self._create_PCA_linkage('LG1')

    def _create_sample_filter_file(self):
        self.plink_master_vcf = self.out_dir + '/vcf_filtered_samples.vcf.gz'
        self.good_samples_txt = self.out_dir + '/samples_to_keep.txt'

        if self.ecogroups == ['All']:
            self.ecogroups = ['Mbuna', 'Utaka', 'Shallow Benthic', 'Deep Benthic','Rhampochromis', 'Diplotaxodon', 'Riverine', 'AC']
        elif self.ecogroups == ['Non_Riverine']:
            self.ecogroups = ['Mbuna', 'Utaka', 'Shallow Benthic', 'Deep Benthic','Rhampochromis', 'Diplotaxodon']

        self.s_dt = pd.read_excel(self.sample_database)
        pd.DataFrame(self.s_dt[self.s_dt.Ecogroup.isin(self.ecogroups)].SampleID.unique()).to_csv(self.good_samples_txt, header = False, index = False)
        subprocess.run(['bcftools', 'view', self.in_vcf, '--samples-file', self.good_samples_txt, '-o', self.plink_master_vcf, '-O', 'z'])
        subprocess.run(['bcftools', 'index', self.plink_master_vcf])
   
    def _create_PCA_linkage(self, lg):
        contig_name = self.linkage_group_map[lg]
        pathlib.Path(self.out_dir + '/PCA/' + lg + '/').mkdir(parents=True, exist_ok=True) 

        subprocess.run(['bcftools', 'filter', '-r', contig_name, self.plink_master_vcf, '-o', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '-O', 'z'])
        subprocess.run(['plink', '--vcf', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--indep-pairwise', '50', '10', '0.1', '--out', 'test'])
        subprocess.run(['plink', '--vcf', self.out_dir + '/PCA/' + lg + '/' + lg + '.vcf.gz', '--double-id', '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--extract', 'test.prune.in', '--make-bed', '--pca', '--out', 'test'])

        # pdb.set_trace()

pca_obj = PCA_Maker(args.input_vcffile, args.output_dir, args.sample_database, args.ecogroups)

"""
PCA Steps:
1. ensure raw file is indexed and gzipped
2. filter variants by lg (bcftools filter -r <region> <args.file> > <args.output_dir/lg>)
3. Implement filtering to get rid of unwanted samples (vcftools view <prev_output> --samples_file <file> > <filtered_output; similar pathing as above>)
4. define filepath to file for input to plink... can probably define this earlier? 
5. run plik code and output to a plink dir within each LG
6. implement R code
    a. import libraries
    b. set path to files in the LG you'e analyzing
    c. run through each line of code and output to out_dir



gzip the whole file with LG7 added at the end while I try to resolve issues with bcftools.
bgzip -c master_file.vcf > lg1-22_master_file.vcf.gz --threads 88
tabix -p vcf lg1-22_master_file.vcf.gz
filter out variants for lg1 and lg11 only
bcftools filter -r NC_036780.1 lg1-22_master_file.vcf.gz > /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/PCA/lg1/lg1.vcf
bcftools filter -r NC_036790.1 lg1-22_master_file.vcf.gz > /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/PCA/lg11/lg11.vcf
BEFORE THIS STEP, RESOLVE THE ISSUES WITH THE MISSING FILES USING PANDAS: 
Filter out riverine, AC, and missing Samples:
Create a file with only the sample names & no extra metadata
bcftools view /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/PCA/lg1/lg1.vcf --samples-file only_sample_names.csv > /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/PCA/lg1/filtered_lg1.vcf
bcftools view /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/PCA/lg11/lg11.vcf --samples-file only_sample_names.csv > /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/PCA/lg11/filtered_lg11.vcf
MAY NOT BE REQUIRED: Remove any potential duplicate entries:
bcftools norm --rm-dup all /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/PCA/lg1/filtered_lg1.vcf > /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/PCA/lg1/no_dup_filtered_lg1.vcf
bcftools norm --rm-dup all /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/PCA/lg11/filtered_lg11.vcf > /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/PCA/lg11/no_dup_filtered_lg11.vcf

Define path to VCF file:
From here on out, work with one VCF file at a time and change into the corresponding dir for ease
LG1: VCF=$PWD/filtered_lg1.vcf
Linkage Pruning:
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out cichlids
Perform PCA Calculations
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract cichlids.prune.in \
--make-bed --pca --out cichlids
R is not cooperating on the server, but the output files are manageable in size, so I'll run everything on Rstudio
Load up R packages and R environment
Import in the filtered sample file with the needed metadata_id column
conda activate R
library(gdsfmt)
library(tidyverse)
setwd(</path/to/above/files/>)
Run all the R Commands:
pca <- read_table("./cichlids.eigenvec", col_names = FALSE)
eigenval <- scan("./cichlids.eigenval")
pca <- pca[,-1]
# set names
names(pca)[1] <- "sample"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
new_names <- read.csv('filtered_samples.csv')
pca$sample=new_names$metadata_id
eco_group <- rep(NA, length(pca$sample))
eco_group[grep("Mbuna", pca$sample)] <- "Mbuna"
eco_group[grep("Utaka", pca$sample)] <- "Utaka"
eco_group[grep("Deep_Benthic", pca$sample)] <- "Deep_Benthic"
eco_group[grep("Shallow_Benthic", pca$sample)] <- "Shallow_Benthic"
eco_group[grep("Diplotaxodon", pca$sample)] <- "Diplotaxodon"
eco_group[grep("Rhampochromis", pca$sample)] <- "Rhampochromis"
pca <- as_tibble(data.frame(pca, eco_group))
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
b <- ggplot(pca, aes(PC1, PC2, col = eco_group)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue", "green", "brown", "purple", "pink"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
Save the output:
ggsave('filename', plot=b, device='png')
"""