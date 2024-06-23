from helper_modules.file_manager_Replacement import FileManager as FM
from Bio import AlignIO
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import subprocess, pdb

#https://dgenies.toulouse.inra.fr/

inversions = {'LG2':('NC_036781.1',19705000,19748000,43254805,43658853),'LG9':('NC_036789.1',14453796,15649299,32255605,33496468),'LG10':('NC_036790.1',11674905,11855817,29898615,29898615),
				'LG11':('NC_036791.1',8302039,8309764,30371888,30459686),'LG13':('NC_036792.1',2249541,2453698,23046928,23131968),'LG20':('NC_036799.1',19614379,19689710,32872827,33764042)}


fm_obj_mz = FM(genome_version = 'Mzebra_GT3')
fm_obj_yh = FM(genome_version = 'kocher_YH_female')

#fm_obj_mz.downloadData(fm_obj_mz.localGenomeFile)
#fm_obj_yh.downloadData(fm_obj_yh.localGenomeFile)

#subprocess.run(['GSAlign','-dp','-i',fm_obj_mz.localGenomeFile,'-q',fm_obj_yh.localGenomeFile, '-o', fm_obj_mz.localGenomesDir + 'MZ_YH_Alignment'])

processes = []
for inv_contig in ['LG2','LG9','LG10','LG11','LG13','LG20']:
	subprocess.run(['faidx', fm_obj_mz.localGenomeFile, inversions[inv_contig][0], '-o', fm_obj_mz.localGenomesDir + inversions[inv_contig][0] + '_MZ.fa'])
	subprocess.run(['faidx', fm_obj_yh.localGenomeFile, inversions[inv_contig][0], '-o', fm_obj_mz.localGenomesDir + inversions[inv_contig][0] + '_YH.fa'])

	command = ['lastz', fm_obj_mz.localGenomesDir + inversions[inv_contig][0] + '_MZ.fa', fm_obj_mz.localGenomesDir + inversions[inv_contig][0] + '_YH.fa']
	command += ['--strand=both','--nochain','--gap=600,150', '--gappedthresh=3000', '--masking=254', '--hspthresh=4500']
	command += ['--allocate:traceback=1.99G', '--output=' + fm_obj_mz.localGenomesDir + inv_contig + '.maf', '--format=maf']

	processes.append(subprocess.Popen(command))
for p in processes:
	p.communicate()

pdb.set_trace()
aln = AlignIO.parse(fm_obj_mz.localGenomesDir + 'MZ_YH_Alignment.maf', "maf")

dt = pd.read_csv('map_A_spYH_GT1_to_Mzebra_GT3.paf', sep = '\t', names = ['QueryContig','QueryContigSize','QueryStart','QueryStop','QueryStrand','RefContig','RefContigSize','RefStart','RefStop','ResiduesMatch','AlignmentBlockLength','MappingQuality','b','c','d','e','f','g'])
dt = dt[['RefContig','RefStart','RefStop','QueryContig','QueryStart','QueryStop','QueryStrand','ResiduesMatch','AlignmentBlockLength','MappingQuality']]
LG = dt[(dt.QueryContig == 'NC_036781.1') & (dt.RefContig == 'NC_036781.1')].sort_values('RefStart')
pdb.set_trace()




dt = dt[['RefContig','RefStart','RefStop','QueryContig','QueryStart','QueryStop','QueryStrand']]
mz = []
yh = []
dt = pd.DataFrame(columns = ['RefContig','RefStart','RefStop','RefSize','QueryContig','QueryStart','QueryStop','QuerySize'])
for ma in aln:
	ref_contig = ma[0].id 
	ref_start = ma[0].annotations['start']
	ref_length = ma[0].annotations['size']
	ref_stop = ref_start + ma[0].annotations['strand'] * ref_lengthLG
	qry_contig = ma[1].id 
	qry_start = ma[1].annotations['start']
	qry_length = ma[1].annotations['size']
	qry_stop = qry_start + ma[1].annotations['strand'] * qry_length

	dt.loc[len(dt)] = [ref_contig,ref_start,ref_stop,ref_length,qry_contig,qry_start,qry_stop,qry_length]

LG = dt[(dt.QueryContig == 'NC_036781.1') & (dt.TargetContig == 'NC_036781.1')].sort_values('QueryStart')
pdb.set_trace()

plt.scatter(mz,yh)
plt.show()
