import pandas as pd
import seaborn as sns 
import subprocess, pdb, os, allel
from helper_modules.file_manager_Replacement import FileManager as FM
from Bio import AlignIO
import matplotlib.pyplot as plt

inversions = {'LG2':('NC_036781.1',19743639,20311160,43254805,43658853),'LG9':('NC_036788.1',14453796,15649299,32255605,33496468),'LG10':('NC_036789.1',11674905,11855817,29898615,29898615),
				'LG11':('NC_036790.1',8302039,8309764,30371888,30459686),'LG13':('NC_036792.1',2249541,2453698,23046928,23131968),'LG20a':('NC_036799.1',19614379,19689710,29716509,29673642),
				'LG20b':('NC_036799.1',29716509,29673642,32872827,33764042)}


def parallel_filter(input_vcf, output_vcf, samples):
	linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
							  'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
							  'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
							  'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

	base_path = os.path.dirname(output_vcf)
	processes = []
	vcf_files = []

	for lg in linkageGroups.keys():
		lg_vcf = base_path + '/' + lg + '.vcf.gz'
		vcf_files.append(lg_vcf)
		processes.append(subprocess.Popen(['bcftools','view','-m','2','-M','2','-V','indels,other','-s',','.join(samples), '--min-ac','1:minor', '-r', lg, '-o', lg_vcf, '-O', 'b', input_vcf]))

	for p in processes:
		p.communicate()

	for lg in linkageGroups.keys():
	   subprocess.run(['bcftools','index', base_path + '/' + lg + '.vcf.gz'])

	subprocess.run(['bcftools', 'concat'] + vcf_files + ['-o',output_vcf, '-O','z'])

	for lg in linkageGroups.keys():
		subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz'])
		subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz.csi'])

	subprocess.run(['bcftools','index', output_vcf])


# Read and download essential data
fm_obj = FM(genome_version = 'Mzebra_GT3')
main_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/gt3_cohort_pass_variants.vcf.gz'
#fm_obj.downloadData(main_vcf)
fm_obj.downloadData(main_vcf + '.tbi')
fm_obj.downloadData(fm_obj.localSampleFile_v2)
s_dt = pd.read_excel(fm_obj.localSampleFile_v2, sheet_name = 'SampleLevel')


#############################################################
# 1. Create phylogenies for whole genome and each inversion #
#############################################################
# Create vcf file for lab sample data and for lake malawi phylogeny data (defined in Sample Database)
malawi_samples = s_dt[s_dt.CorePCA == 'Yes'].SampleID.to_list()
c_dt = s_dt[s_dt.CorePCA == 'Yes']
malawi_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals.vcf.gz'
#parallel_filter(main_vcf, malawi_vcf, malawi_samples)
subprocess.run(['bcftools','index', '-f', malawi_vcf])

#Create lg subsets for each inversion from the malawi samples vcf file
for lg,(contig,temp,left,right,temp2) in inversions.items():
	region = contig + ':' + str(left) + '-' + str(right)
	lg_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/MalawiIndividuals_' + lg + '.vcf.gz'
	subprocess.run(['bcftools','view','-r', region, '-o', lg_vcf, '-O', 'z', malawi_vcf])
	subprocess.run(['bcftools','index', lg_vcf])


vcf_obj = allel.read_vcf(malawi_vcf)
test = allel.GenotypeArray(vcf_obj['calldata/GT'])
out = allel.pairwise_distance(test.to_n_alt(), metric='cityblock')
import scipy.spatial
sq = scipy.spatial.distance.squareform(out)
pdb.set_trace()

