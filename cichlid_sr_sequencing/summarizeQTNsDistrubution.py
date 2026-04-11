import pdb, pysam
import pandas as pd
import numpy as np
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from helper_modules.file_manager import FileManager as FM

dt = pd.read_csv('candidateQTNs_LakeMalawiGenotypes.csv', index_col = 0)
dt2 = pd.read_csv('candidateQTNs_LakeMalawiGenotypes_Blumer.csv', index_col = 0)

fm_obj = FM('Mzebra_GT3')
fm_obj.downloadData(fm_obj.localGenomeFile)
fasta_old = pysam.FastaFile(fm_obj.localGenomeFile)
fm_obj = FM('Mzebra_GT3_NCBI')
fm_obj.readSampleDatabase()
fm_obj.downloadData(fm_obj.localGenomeFile)
fasta_new = pysam.FastaFile(fm_obj.localGenomeFile)
mapping = {k:v for k,v in zip(fasta_old.references[:-2],fasta_new.references[:-1])}
pdb.set_trace()

#order = ['YH_fam1_males', 'YH_fam2_males','YH_parental_males','DeepBenthicHets','YH_fam1_females','YH_fam2_females','YH_parental_females','DeepBenthicInverted','MC_males', 'MC_females','LabCVs_male','LabCVs_female','DeepBenthicNormal','ShallowBenthics','Utaka','Mbuna','ACs','Diplotaxodon','Rhamphochromis']
special = [24869679,24870603,24918628,24921178,24921394] #[i,d,i,i,i]

for lv in special:
    category_counts = defaultdict(lambda: np.zeros(2))
    row = dt2.loc[lv]
    for sample,geno in row.items():
        inv10_geno = fm_obj.sample_dt.loc[fm_obj.sample_dt.SampleID == sample,'Inversion10'].values[0]
        species = fm_obj.sample_dt.loc[fm_obj.sample_dt.SampleID == sample,'Species'].values[0]
        ac = np.array([int(x) for x in row[sample].split(',')])
        if ac.sum() == 0:
            category_counts[(species,inv10_geno)] += (1,0)
        else:
            category_counts[(species,inv10_geno)] += (0,1) 
    pdb.set_trace()

for position, row in dt2.iterrows():
    category_counts = defaultdict(lambda: np.zeros(2))
    for sample in dt2.columns:
        inv10_geno = fm_obj.sample_dt.loc[fm_obj.sample_dt.SampleID == sample,'Inversion10'].values[0]
        species = fm_obj.sample_dt.loc[fm_obj.sample_dt.SampleID == sample,'Species'].values[0]
        try:
            ac = np.array([int(x) for x in row[sample].split(',')])
            if ac.sum() > 4:
                if ac[0]/ac.sum() > .9: # Reference
                    category_counts[(species,inv10_geno)] += np.array([2,0])
                elif ac[0]/ac.sum() < .1: # Mutant
                    category_counts[(species,inv10_geno)] += np.array([0,2])
                else:
                    category_counts[(species,inv10_geno)] += np.array([1,1])
        except AttributeError:
            pdb.set_trace()
            continue
    try:
        out_dt.loc[position] = [0 if category_counts[x].sum() == 0 else category_counts[x][1]/category_counts[x].sum() for x in order]
    except:
        order = category_counts.keys()
        out_dt = pd.DataFrame(columns = order)
        out_dt.loc[position] = [0 if category_counts[x].sum() == 0 else category_counts[x][1]/category_counts[x].sum() for x in order]
        pdb.set_trace()
pdb.set_trace()

out_dt = out_dt[(out_dt.YH_fam1_males > 0.3) & (out_dt.YH_fam2_males > 0.3) & (out_dt.YH_parental_males > 0.3) & (out_dt.DeepBenthicHets > 0.3)]

row_linkage = linkage(out_dt[['YH_fam1_males', 'YH_fam2_males','YH_parental_males','DeepBenthicHets','YH_fam1_females','YH_fam2_females','YH_parental_females','DeepBenthicInverted','MC_males', 'MC_females','LabCVs_male','LabCVs_female']], method='average')
sns.clustermap(out_dt, row_linkage=row_linkage, col_cluster = False)
plt.show()
sns.clustermap(out_dt[['YH_fam1_males', 'YH_fam2_males','YH_parental_males','DeepBenthicHets','DeepBenthicNormal']], col_cluster=False)
plt.show()
pdb.set_trace()