import allel, pdb, os, subprocess

from helper_modules.file_manager_Replacement import FileManager as FM

def parallel_filter(input_vcf, output_vcf, samples):
	linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
							  'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
							  'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
							  'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

	base_path = os.path.dirname(output_vcf)
	processes = []
	vcf_files = []

	min_count = int(len(samples)*2/5)

	for lg in linkageGroups.keys():
		lg_vcf = base_path + '/' + lg + '.vcf.gz'
		vcf_files.append(lg_vcf)
		processes.append(subprocess.Popen(['bcftools','view','-m','2','-M','2','-V','indels,other','-s',','.join(samples), '--min-ac',str(min_count) + ':minor', '-r', lg, '-o', lg_vcf, '-O', 'b', input_vcf]))

	for p in processes:
		p.communicate()

	for lg in linkageGroups.keys():
	   subprocess.run(['bcftools','index', base_path + '/' + lg + '.vcf.gz'])

	subprocess.run(['bcftools', 'concat'] + vcf_files + ['-o',output_vcf, '-O','z'])

	for lg in linkageGroups.keys():
		subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz'])
		subprocess.run(['rm', base_path + '/' + lg + '.vcf.gz.csi'])

	subprocess.run(['bcftools','index', output_vcf])

fm_obj = FM(genome_version = 'Mzebra_GT3')
parents1 = ['YH_006_f','YH_011_m']
progeny1 = ['YH_016','YH_017','YH_018','YH_019','YH_020','YH_021','YH_022','YH_023','YH_024','YH_025',
			'YH_026','YH_027','YH_028','YH_029','YH_030','YH_031','YH_032','YH_033','YH_037',
			'YH_038','YH_039','YH_040','YH_041','YH_042']
parents2 = ['YH_005_m','YH_010_f']
progeny2 = ['YH_046','YH_047','YH_048','YH_049','YH_050','YH_051','YH_052','YH_053','YH_054','YH_055',
			'YH_056','YH_057','YH_058','YH_059','YH_060','YH_061','YH_062','YH_063','YH_064',
			'YH_065','YH_066','YH_068','YH_069','YH_070']
yh_brood1 = parents1 + progeny1
yh_brood2 = parents2 + progeny2

yh_raw_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/YHPedigree/YHPedigree_pass_variants_master_file.vcf.gz'
fm_obj.downloadData(yh_raw_vcf)
fm_obj.downloadData(yh_raw_vcf + '.tbi')

yh_brood1_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/YHPedigree/YHPedigree_Brood1.vcf.gz'
yh_brood2_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/YHPedigree/YHPedigree_Brood2.vcf.gz'

#parallel_filter(yh_raw_vcf, yh_brood1_vcf, yh_brood1)
#parallel_filter(yh_raw_vcf, yh_brood2_vcf, yh_brood2)

vcf_obj1 = allel.read_vcf(yh_brood1_vcf)

master = allel.GenotypeArray(vcf_obj1['calldata/GT'])
master = master[(master[:,0,0] != -1)&(master[:,1,0] != -1)]
for i in range(26,2,-1):
	temp = master[:,:i,:]
	g = allel.phase_by_transmission(temp, 25)
	print(g.is_phased.sum(axis = 0)[0:2])
pdb.set_trace()

#parents = ['MC-5B6B-f','YH_005_m']
#progeny = ['MCYHF1_010_m','MCYHF1_011_m','MCYHF1_012_m','MCYHF1_013_m','MCYHF1_014_f','MCYHF1_015_f','MCYHF1_016_f','MCYHF1_017_f']


pdb.set_trace()