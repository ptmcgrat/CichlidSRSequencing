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

fm_obj = FM('Mzebra_GT3_NCBI')
fm_obj.readSampleDatabase()

#order = ['YH_fam1_males', 'YH_fam2_males','YH_parental_males','DeepBenthicHets','YH_fam1_females','YH_fam2_females','YH_parental_females','DeepBenthicInverted','MC_males', 'MC_females','LabCVs_male','LabCVs_female','DeepBenthicNormal','ShallowBenthics','Utaka','Mbuna','ACs','Diplotaxodon','Rhamphochromis']
special = [24869679,24870603,24918628,24921178,24921394] #[i,d,i,i,i]

out_dt = pd.DataFrame(columns = ['Species','Inv10'])
for lv in special:
    out_dt[lv] = ''
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
    for k,v in category_counts.items():
        if len(out_dt[(out_dt.Species == k[0]) & (out_dt.Inv10 == k[1])]) == 0:
            out_dt.loc[len(out_dt)] = [k[0],k[1],','.join([str(int(x)) for x in v])]
        else:
            out_dt.loc[(out_dt.Species == k[0]) & (out_dt.Inv10 == k[1]),lv] = ','.join([str(int(x)) for x in v])


out_dt2 = pd.DataFrame(columns = [0,1,2])

for position, row in dt2.iterrows():
    if position in special:
        continue
    category_counts = {0: np.zeros(3),1: np.zeros(3),2: np.zeros(3)}
    for sample in dt2.columns:
        inv10_geno = int(fm_obj.sample_dt.loc[fm_obj.sample_dt.SampleID == sample,'Inversion10'].values[0])
        species = fm_obj.sample_dt.loc[fm_obj.sample_dt.SampleID == sample,'Species'].values[0]
        try:
            ac = np.array([int(x) for x in row[sample].split(',')])
            if ac.sum() > 4:
                if ac[0]/ac.sum() > .9: # Reference
                    category_counts[inv10_geno] += np.array([1,0,0])
                elif ac[0]/ac.sum() < .1: # Mutant
                    category_counts[inv10_geno] += np.array([0,0,1])
                else:
                    category_counts[inv10_geno] += np.array([0,1,0])
        except AttributeError:
            continue
    out_dt2.loc[position] = [category_counts[0],category_counts[1],category_counts[2]]
pdb.set_trace()
"""
try:
    out_dt.loc[position] = [0 if category_counts[x].sum() == 0 else category_counts[x][1]/category_counts[x].sum() for x in order]
except:
    order = category_counts.keys()
    out_dt = pd.DataFrame(columns = order)
    out_dt.loc[position] = [0 if category_counts[x].sum() == 0 else category_counts[x][1]/category_counts[x].sum() for x in order]
    pdb.set_trace()
"""
pdb.set_trace()

out_dt = out_dt[(out_dt.YH_fam1_males > 0.3) & (out_dt.YH_fam2_males > 0.3) & (out_dt.YH_parental_males > 0.3) & (out_dt.DeepBenthicHets > 0.3)]

row_linkage = linkage(out_dt[['YH_fam1_males', 'YH_fam2_males','YH_parental_males','DeepBenthicHets','YH_fam1_females','YH_fam2_females','YH_parental_females','DeepBenthicInverted','MC_males', 'MC_females','LabCVs_male','LabCVs_female']], method='average')
sns.clustermap(out_dt, row_linkage=row_linkage, col_cluster = False)
plt.show()
sns.clustermap(out_dt[['YH_fam1_males', 'YH_fam2_males','YH_parental_males','DeepBenthicHets','DeepBenthicNormal']], col_cluster=False)
plt.show()
pdb.set_trace()