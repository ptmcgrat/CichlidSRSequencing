from helper_modules.file_manager import FileManager as FM
import os, subprocess, pdb, allel
from multiprocessing import cpu_count
from types import SimpleNamespace
import pandas as pd
import numpy as np
from hmmlearn.hmm import CategoricalHMM
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
 
class GenotypeBC:
    def __init__(self):
        self.fm_obj = FM()
        self.fm_obj.setGenome('Mconophoros_GT1')
        self.master_directory = self.fm_obj.localMasterDir + '2026_Nikesh/'
        self.broods = ['MCYH-BC1-24-01-10-1','MCYH-BC1-24-03-17-1','MCYH-BC1-24-03-11-1','MCYH-BC1-24-04-22-1','MCYH-BC1-24-03-27-1']
        self.yh_sire = 'YH_008_m'
        self.mc_dam = 'MC-1R6R-f'
        self.mothers = ['MC-3R6R-f','MC-1G2G-f','MC-5G11G-f']
        self.f1_dads = ['MCYHF1-002','MCYHF1-003']
        self.parents_unfiltered_vcf = self.master_directory + 'AllBroodParents_Unfiltered.vcf.gz'
        self.parents_filtered_vcf = self.master_directory + 'AllBroodParents_Filtered.vcf.gz'
        self.all_pileup = self.master_directory + 'BC_Pileup.vcf.gz'
        self.all_vcf = self.master_directory + 'BC_Variants.vcf.gz'
        self.all_samples = self.master_directory + 'BC_Pileup_Samples.txt'
        self.linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
                              'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
                              'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
                              'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

    def readPhenoDatabase(self):
        dt = pd.read_csv("https://docs.google.com/spreadsheets/d/1BovaQm-FaOzchci9By71xTh3MzKCqMSpK3DCWB9sGiE/export?gid=1566739159&format=csv")
        pheno_dt = dt[(dt.species.isin(['MCYHBC1'])) & (dt.tagtype == 'PIT')][['samplename','sex','sex_100','implant_date','DOF','standard_length_cm','body_mass_g','Dead']]           
        pheno_dt['implant_date'] = pd.to_datetime(pheno_dt.implant_date)
        pheno_dt['DOF'] = pd.to_datetime(pheno_dt.DOF)
        pheno_dt['AgeAtInjection'] = (pheno_dt.implant_date - pheno_dt.DOF).dt.days
        pheno_dt['StandardLength'] = pheno_dt.standard_length_cm.astype('float')
        pheno_dt['BodyMass'] = pheno_dt.body_mass_g.astype('float')
        #pheno_dt = dt[(dt.species.isin(['MCYHBC1'])) & (dt.tagtype == 'PIT') & (dt.samplename.isin(geno_dt.index))][['samplename','sex','sex_100']]           
        pheno_dt.loc[pheno_dt.sex == 'PF','sex'] = 'F'
        #pheno_dt.loc[pheno_dt.sex == 'M','sex'] = '1'
        #pheno_dt.loc[pheno_dt.sex == 'F','sex'] = '0'

        pheno_dt.loc[pheno_dt.sex_100.isna(),'sex_100'] = '-'
        #pheno_dt.loc[pheno_dt.sex_100 == 'M','sex_100'] = '1'
        #pheno_dt.loc[pheno_dt.sex_100 == 'F','sex_100'] = '0'
        self.pheno_dt = pheno_dt

    def createSampleList(self):
        self.fm_obj.readSampleDatabase()
        with open(self.all_samples, 'w') as f:
            all_samples = [self.yh_sire] + [self.mc_dam] + self.mothers + self.f1_dads
            all_samples += self.fm_obj.sample_dt[self.fm_obj.sample_dt.Species == 'YHxMC x MC BC'].SampleID.tolist()
            for sample in all_samples:
                print('Downloading sample: ' + sample)
                self.fm_obj.createSampleFiles(sample)
                self.fm_obj.downloadData(self.fm_obj.localSampleBamDir)
                if not os.path.exists(self.fm_obj.localGVCFFile + '.tbi'):
                    output = subprocess.run(['gatk','IndexFeatureFile','-I', self.fm_obj.localGVCFFile], capture_output = True)
                    self.fm_obj.uploadData(self.fm_obj.localGVCFFile + '.tbi')
                print(self.fm_obj.localBamFile, file = f)

    def createMasterVCFs(self):
        # This function assumes that all parents and backcross samples have been 
        # sequenced and aligned to the Mconophoros_GT3 reference

        fm_obj = self.fm_obj
        print('Downloading reference data')
        fm_obj.downloadData(fm_obj.localGenomeDir)

        #print('Downloading YH sire data')
        #fm_obj.createSampleFiles(self.yh_sire) # Start with the YH father
        #fm_obj.downloadData(fm_obj.localSampleBamDir)
        
        print('Identifying small variants from YH sire and MC mothers')
        command = ['gatk','CombineGVCFs','-R',fm_obj.localGenomeFile]
        fm_obj.createSampleFiles(self.yh_sire) # Start with the YH father
        command += ['--variant',fm_obj.localGVCFFile]
        fm_obj.createSampleFiles(self.mc_dam) # Add the MC mother
        command += ['--variant',fm_obj.localGVCFFile]
        for mother in self.mothers:
            fm_obj.createSampleFiles(mother) # Add the MC mother
            command += ['--variant',fm_obj.localGVCFFile]

        command += ['-O', self.parents_unfiltered_vcf + '.temp_cohort.g.vcf']
        #output = subprocess.run(command, capture_output = True)
        command = ['gatk','GenotypeGVCFs','-R',fm_obj.localGenomeFile,'-V',self.parents_unfiltered_vcf + '.temp_cohort.g.vcf', '-O',self.parents_unfiltered_vcf]
        #output = subprocess.run(command, capture_output = True)
        #subprocess.run(['rm','-f', self.parents_unfiltered_vcf + '.temp_cohort.g.vcf'])

        command = ['bcftools','view','-v','snps','-i','GT[4]="AA" && GT[1]="RR" & FORMAT/DP[4] >12 && FORMAT/DP[1] >12','-o', self.parents_filtered_vcf, '-O','z', self.parents_unfiltered_vcf]
        #output = subprocess.run(command, capture_output = True)
        #pdb.set_trace()
        print('Piling up variants from BCs')
        subprocess.run(['ulimit','-n','4096'])
        output = subprocess.run(['bcftools', 'mpileup', '-f', self.fm_obj.localGenomeFile, '-R', self.parents_filtered_vcf, '-b', self.all_samples, '-a', 'AD', '-Q', '0', '-o', self.all_pileup, '-Oz', '--threads', str(cpu_count())], capture_output = True)
        output = subprocess.run(['bcftools', 'call', '-c','-o',self.all_vcf,'-O','z',self.all_pileup], capture_output = True)

    def createHaplotypesBCs(self):
        all_data = []
        bc_pileup = self._loadVCF(self.all_pileup)
        bc_vcf = self._loadVCF(self.all_vcf)
        brood_masks = self.createBroodMasks(bc_vcf)
        os.makedirs(self.master_directory + 'BC_HMM_Maps/', exist_ok = True)
        for brood_idx in range(5):
            mother, bcs = self.retBroodInfo(brood_idx)
            brood_mask = brood_masks[mother]
            available_bcs = [x for x in bcs if x in bc_pileup.samples]
            for offspring in available_bcs:
                offspring_idx = int((bc_pileup.samples == offspring).argmax())
                off_ad = bc_pileup.allelic_depth[brood_mask, offspring_idx]
                off_chromosomes = bc_pileup.chromosomes[brood_mask]
                off_positions = bc_pileup.positions[brood_mask]
                str_ad, str_chr, str_pos = self._stretchADVector(off_ad, off_chromosomes, off_positions)
                data = self._hmmVector(offspring, str_ad, str_chr, str_pos)
                all_data.append(data)
        merged_dt = reduce(lambda left, right: pd.merge(left, right, on=['Chromosome','Position'], how='inner'), all_data)
        merged_dt.to_csv(self.master_directory + 'OutputGenotypes.csv')
        self.plotRecombinationByChromosome()
        self.fm_obj.uploadData(self.master_directory)

    def createBroodMasks(self, data):
        
        yh_sire_idx = int((data.samples == self.yh_sire).argmax())
        mc_dam_idx = int((data.samples == self.mc_dam).argmax())

        yh_sire_mask = data.converted_genotype[:,yh_sire_idx].is_hom_alt()
        mc_dam_mask = data.converted_genotype[:,mc_dam_idx].is_hom_ref()
        brood_masks = {}
        for mother in self.mothers:
            mother_idx = int((data.samples == mother).argmax())
            mother_mask = data.converted_genotype[:,mother_idx].is_hom_ref()

            brood_masks[mother] = yh_sire_mask & mc_dam_mask & mother_mask
        return brood_masks

    def retBroodInfo(self, index):
        #dt = pd.read_csv("https://docs.google.com/spreadsheets/d/1BovaQm-FaOzchci9By71xTh3MzKCqMSpK3DCWB9sGiE/export?gid=1566739159&format=csv")
        #f_dt = dt[(dt.species.isin(['MCYHBC1'])) & (dt.tagtype == 'PIT')][['brood_ID','samplename']]
        self.fm_obj.readSampleDatabase()
        bc_dt = self.fm_obj.sample_dt[self.fm_obj.sample_dt.BroodID == self.broods[index]]
        mother = bc_dt.Parents.iloc[0].split(',')[1]
        return mother, bc_dt.SampleID.tolist()
        
    def fineMapHelper(self, vcf_file, samples, chrom, start, stop):
        self._loadVCF(vcf_file)
        fig, axs = plt.subplots(len(samples),1, figsize = (5,10))
        brood_masks = self.createBroodMasks2(self.mothers)
        location_mask = (self.chromosomes == chrom) & (self.positions > start) & (self.positions < stop) 
        tot_idx = 0
        for brood_idx in range(5):
            mother, bcs = self.retBroodInfo(brood_idx)
            brood_mask = brood_masks[mother]
            available_bcs = [x for x in bcs if x in samples]
            for offspring in available_bcs:
                offspring_idx = int((self.samples == offspring).argmax())
                off_ad = self.allelic_depth[brood_mask & location_mask, offspring_idx]
                off_chromosomes = self.chromosomes[brood_mask & location_mask]
                off_positions = self.positions[brood_mask & location_mask]
                axs[tot_idx].plot(off_positions, off_ad)
                axs[tot_idx].set_title(offspring)
                tot_idx += 1
        plt.tight_layout()
        plt.show()

    def fineMapRegion(self):
        self.readPhenoDatabase()
        self.fm_obj.downloadData(self.master_directory + 'OutputGenotypes.csv')
        merged_dt = pd.read_csv(self.master_directory + 'OutputGenotypes.csv', index_col = 0)
        chrom_dt = merged_dt[merged_dt.Chromosome == 'NC_036789.1']
        mapped_region = (chrom_dt.Position >= 26000000) & (chrom_dt.Position <= 28200000)
        
        chrom_dt = chrom_dt.drop(columns=['Chromosome','Position'])
            
        
        temp_dt = chrom_dt.copy()
        temp_dt['QTL'] = mapped_region
        small_genos = temp_dt[temp_dt.QTL].mean()
        informative_bc = small_genos[(small_genos>0)&(small_genos<1)]
        self.informative_bcs = informative_bc.index.to_list()
        sex_dt = self.pheno_dt[self.pheno_dt.samplename.isin(self.informative_bcs)][['samplename','sex','sex_100']]
        #pd.merge(temp_dt)
        plot_dt = pd.merge(temp_dt[temp_dt.QTL].T, sex_dt, left_index = True, right_on = 'samplename').set_index('samplename')
        order = ['MCYHBC1-709-1','MCYHBC1-825-1','MCYHBC1-737-1','MCYHBC1-765-1','MCYHBC1-800-1','MCYHBC1-638-1','MCYHBC1-788-1','MCYHBC1-814-1','MCYHBC1-619-1','MCYHBC1-790-1','MCYHBC1-656-1','MCYHBC1-645-1','MCYHBC1-580-1','MCYHBC1-668-1','MCYHBC1-919-1','MCYHBC1-642-1','MCYHBC1-828-1','MCYHBC1-755-1','MCYHBC1-701-1','MCYHBC1-795-1','MCYHBC1-927-1']
        plot_dt = plot_dt.reindex(order)
        plot_dt = plot_dt[plot_dt.sex_100 != '-']

        plot_color = plot_dt['sex_100']

        plot_dt = plot_dt.drop(columns = ['sex','sex_100'])
        new_names = [merged_dt.loc[x].Position for x in plot_dt.columns.tolist()]
        plot_dt.columns = new_names
        ax = sns.clustermap(plot_dt.astype(float), col_cluster=False, row_cluster = False, row_colors = plot_color.map({'M':'blue','F':'pink'}), cmap = 'cividis')
        for label in ax.ax_heatmap.get_yticklabels():
            labels_text = label.get_text()
            if self.pheno_dt[self.pheno_dt.samplename == labels_text].Dead.values[0] == 'TRUE':
                label.set_color("red")
            else:
                label.set_color("black")    
        plt.show()

    def plotRecombinationByChromosome(self):

        layout = """
        ABCDEFGH
        IJKLMNOP
        QRSTUVWW
        """
        fig, axs = plt.subplot_mosaic(layout, layout="constrained", figsize = (18,9))
        for chrom,letter in zip(self.linkageGroups.keys(),'ABCDEFGHIJKLMNOPQRSTUV'):
            chrom_dt = merged_dt[merged_dt.Chromosome == chrom]
            chrom_dt = chrom_dt.drop(columns=['Chromosome','Position'])
            sort_data = pd.concat([chrom_dt.iloc[0], chrom_dt.mean()], axis=1, keys = ['First','Mean'])
            sorted_cols = sort_data.sort_values(['First','Mean']).index
            chrom_dt = chrom_dt[sorted_cols]
            if chrom == 'NC_036789.1':
                temp_dt = chrom_dt.copy()
                temp_dt['QTL'] = mapped_region
                small_genos = temp_dt[temp_dt.QTL].mean()
                informative_bc = small_genos[(small_genos>0)&(small_genos<1)]
                self.informative_bcs = informative_bc.index.to_list()
                pdb.set_trace()
            g = sns.heatmap(chrom_dt.T, xticklabels = False, yticklabels=False, cbar = False, ax = axs[letter])

            axs[letter].set_title(self.linkageGroups[chrom]) 
        plt.savefig(self.master_directory + 'RecombinantsByChromosome.pdf', bbox_inches='tight', format='pdf')
        plt.close()

    def _loadVCF(self, vcf_file):

        vcf_dict = allel.read_vcf(vcf_file, fields=['variants/CHROM','variants/POS','calldata/GT','calldata/AD','calldata/DP', 'samples'])
        return SimpleNamespace(samples = vcf_dict['samples'], chromosomes = vcf_dict['variants/CHROM'],
            positions = vcf_dict['variants/POS'], genotype = vcf_dict['calldata/GT'],
            depth = vcf_dict['calldata/DP'], allelic_depth = vcf_dict['calldata/AD'][:,:,0:2],
            converted_genotype = allel.GenotypeArray(vcf_dict['calldata/GT']))

    def _hmmVector(self, offspring, ads, chromosomes, positions, window_size = 200, trans_prob = .0000000000000000000001):
        model = CategoricalHMM(n_components = 2, n_features = 2)
        model.startprob_ = np.array([0.5, 0.5])
        model.transmat_ = np.array([[1-trans_prob, trans_prob], [trans_prob, 1-trans_prob]])
        model.emissionprob_ = np.array([[.97, .03],[.6, .4]])
        kernel = np.ones(window_size) / window_size
        layout = """
        ABCDEFGH
        IJKLMNOP
        QRSTUVWW
        """
        fig, axs = plt.subplot_mosaic(layout, layout="constrained", figsize = (18,9))
        for chrom, letter in zip(dict.fromkeys(chromosomes),'ABCDEFGHIJKLMNOPQRSTUV'):
            
            try:
                model_data = ads[chromosomes == chrom]
                pos_chr = positions[chromosomes == chrom]
                chr_chr = chromosomes[chromosomes == chrom]
            except:
                pdb.set_trace()
            y1=model.predict(model_data.reshape(-1,1))
            y1b = np.convolve(model_data, kernel, mode='same')
            y1b[0:100] = y1b[100] # correct boundaries
            y1b[-100:] = y1b[-100] # correct boundaries
            handle = axs[letter].plot(range(len(y1b)),y1b,range(len(y1)),y1)
            axs[letter].set_ylim([-.1,1.1])
            y1_interp = np.interp(range(0,max(pos_chr),100000), pos_chr, y1, left=None, right=None, period=None)
            out_pos = [x for x in range(0,max(pos_chr),100000)]
            out_chr = [chrom] * len(out_pos)
            try:
                out_dt = pd.concat([out_dt, pd.DataFrame({'Chromosome':out_chr, 'Position':out_pos,offspring:y1_interp})])
            except NameError:
                out_dt = pd.DataFrame({'Chromosome':out_chr, 'Position':out_pos,offspring:y1_interp})
        plt.savefig(self.master_directory + 'BC_HMM_Maps/' + offspring + '.pdf', bbox_inches='tight', format='pdf')
        plt.close()
        return out_dt

    def _stretchADVector(self, ad, chromosomes, positions):
        output_shape = ad.sum()
        out_obs = np.zeros(shape = (output_shape), dtype = 'uint8')
        out_chromosomes = np.empty(shape=(output_shape,), dtype=object)
        out_positions = np.zeros(shape = (output_shape), dtype = 'int')
        idx = 0
        for ads, chrom, pos in zip(ad, chromosomes, positions):
            temp_obs = np.zeros(shape = (ads.sum()), dtype = 'uint8')
            for read in range(ads[0]):
                temp_obs[read] = 0
            for read in range(ads[1]):
                temp_obs[ads[0]+read] = 1
            np.random.shuffle(temp_obs)
            for obs in temp_obs:
                out_obs[idx] = obs
                out_chromosomes[idx] = chrom
                out_positions[idx] = pos
                idx+=1
        return out_obs, out_chromosomes, out_positions

gt_bc = GenotypeBC()
#gt_bc.createSampleList()
#gt_bc.createMasterVCFs()
#gt_bc.createHaplotypesBCs()
gt_bc.fineMapRegion()

#gt_bc.identifyHaplotypesInBCs('All_BCs.vcf')
#gt_bc.fineMapHelper('YH_MC_F1s_variant_pileup.vcf', ['MCYHBC1-642-1','MCYHBC1-668-1','MCYHBC1-656-1', 'MCYHBC1-755-1'], 'NC_036789.1', 26500000, 27800000)
linkageGroups = {'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
                              'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 'NC_036789.1':'LG10', 'NC_036790.1':'LG11',
                              'NC_036791.1':'LG12', 'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 'NC_036795.1':'LG16', 'NC_036796.1':'LG17',
                              'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'}

merged_dt = pd.read_csv(gt_bc.fm_obj.localNikeshDir + 'OutputGenotypes.csv', index_col = 0)
#gt_bc.plotRecombinationByChromosome()

gmap_dt = pd.DataFrame()
gmap_dt['chr'] = merged_dt.Chromosome.map(linkageGroups)
gmap_dt['pos'] = merged_dt.Position
gmap_dt.index.name = 'marker'
gmap_dt.to_csv('~/Desktop/QTL_Data/nikesh_gmap.csv')

geno_dt = merged_dt.drop(columns=['Chromosome','Position']).T
geno_dt = merged_dt.drop(columns=['Chromosome','Position']).T.round(0).astype(int)
geno_map = {0:'MC',1:'YH'}
for x in geno_dt.columns:
    geno_dt[x] = geno_dt[x].map(geno_map)
geno_dt.to_csv('~/Desktop/QTL_Data/nikesh_geno.csv')
dt = pd.read_csv("https://docs.google.com/spreadsheets/d/1BovaQm-FaOzchci9By71xTh3MzKCqMSpK3DCWB9sGiE/export?gid=1566739159&format=csv")
pheno_dt = dt[(dt.species.isin(['MCYHBC1'])) & (dt.tagtype == 'PIT') & (dt.samplename.isin(geno_dt.index))][['samplename','sex','sex_100','implant_date','DOF','standard_length_cm','body_mass_g']]           
pheno_dt['implant_date'] = pd.to_datetime(pheno_dt.implant_date)
pheno_dt['DOF'] = pd.to_datetime(pheno_dt.DOF)
pheno_dt['AgeAtInjection'] = (pheno_dt.implant_date - pheno_dt.DOF).dt.days
pheno_dt['StandardLength'] = pheno_dt.standard_length_cm.astype('float')
pheno_dt['BodyMass'] = pheno_dt.body_mass_g.astype('float')
#pheno_dt = dt[(dt.species.isin(['MCYHBC1'])) & (dt.tagtype == 'PIT') & (dt.samplename.isin(geno_dt.index))][['samplename','sex','sex_100']]           
pheno_dt.loc[pheno_dt.sex == 'PF','sex'] = 'F'
pheno_dt.loc[pheno_dt.sex == 'M','sex'] = '1'
pheno_dt.loc[pheno_dt.sex == 'F','sex'] = '0'

pheno_dt.loc[pheno_dt.sex_100.isna(),'sex_100'] = '-'
pheno_dt.loc[pheno_dt.sex_100 == 'M','sex_100'] = '1'
pheno_dt.loc[pheno_dt.sex_100 == 'F','sex_100'] = '0'

#dt2 = pd.read_csv("https://docs.google.com/spreadsheets/d/1YdN1R2na-J8AGbFYluA4yPZp2DQOaQ2qQ7vWNfrUstk/export?gid=0&format=csv")[['samplename','BowerShape']]
#pheno_dt = pd.merge(pheno_dt, dt2, on='samplename', how = 'left')
#pheno_dt['BowerShape'] = pheno_dt['BowerShape'].fillna('-')


pheno_dt = pheno_dt.set_index('samplename')
pheno_dt.index.name = 'id'
pheno_dt[['sex','sex_100','StandardLength','BodyMass']].to_csv('~/Desktop/QTL_Data/nikesh_pheno.csv')
pheno_dt[['sex','sex_100']].to_csv('~/Desktop/QTL_Data/sex_pheno.csv')
pdb.set_trace()
informative_bcs = ['MCYHBC1-709-1', 'MCYHBC1-825-1', 'MCYHBC1-737-1', 'MCYHBC1-800-1', 'MCYHBC1-645-1', 'MCYHBC1-788-1', 'MCYHBC1-814-1', 'MCYHBC1-619-1', 'MCYHBC1-790-1', 'MCYHBC1-656-1', 'MCYHBC1-668-1', 'MCYHBC1-642-1', 'MCYHBC1-701-1', 'MCYHBC1-927-1', 'MCYHBC1-765-1', 'MCYHBC1-638-1', 'MCYHBC1-580-1', 'MCYHBC1-919-1', 'MCYHBC1-828-1', 'MCYHBC1-755-1', 'MCYHBC1-795-1']

sex_dt = merged_dt[(merged_dt.Chromosome == 'NC_036789.1') & (merged_dt.Position >= 26100000) & (merged_dt.Position <= 28100000)]
sex_dt = sex_dt.set_index('Position')
sex_dt = sex_dt.drop(columns=['Chromosome']).T

sex_dt = sex_dt[(sex_dt.mean(axis = 1) > 0) & (sex_dt.mean(axis = 1) < 1)]

plot_dt = pd.merge(sex_dt, pheno_dt[~pheno_dt.sex_100.isna()][['sex_100']], left_index=True, right_index=True)

plot_dt = plot_dt.reindex(['MCYHBC1-642-1', 'MCYHBC1-668-1','MCYHBC1-919-1', 'MCYHBC1-755-1', 'MCYHBC1-656-1', 'MCYHBC1-619-1', 'MCYHBC1-788-1', 'MCYHBC1-814-1', 'MCYHBC1-765-1', 'MCYHBC1-638-1'])
plot_color = plot_dt['sex_100']

plot_dt = plot_dt.drop(columns = ['sex_100'])
sns.clustermap(plot_dt, col_cluster=False, row_cluster = False, row_colors = plot_color.map({'M':'blue','F':'pink'}))
plt.show()

