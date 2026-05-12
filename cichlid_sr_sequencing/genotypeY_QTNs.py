from helper_modules.file_manager import FileManager as FM
from helper_modules.smallVariantGenotyper import normalize_sites
import pdb, os, subprocess
import pandas as pd
from types import SimpleNamespace

def print_sv(outfile, dt):
    with open(outfile,'w') as fp:
        print('##fileformat=VCFv4.2', file = fp)
        print('##source=MSAUniqueVariantCaller', file = fp)
        print('##sample=Y_reg', file = fp)
        print('##contig=<ID=NC_135176.1>', file = fp)
        print('##INFO=<ID=TYPE,Number=1,Type=String,Description="SNV, INS, or DEL">', file = fp)
        print('##INFO=<ID=ALN_COL,Number=1,Type=Integer,Description="1-based alignment column where the variant starts">', file = fp)
        print('##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length of inserted or deleted sequence">', file = fp)
        print('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']), file = fp)
        
        for i,row in dt.iterrows():
            print('\t'.join([row.Chromosome,str(row.Position),row.Name,row.Reference.upper(), row.Alt.upper(), str(row.Q), 'PASS', row.Info]), file = fp)

fm_obj = FM(genome_version = 'Mzebra_GT3_NCBI')
print(fm_obj.localNikeshDir + 'QTG_Candidates/')
os.makedirs(fm_obj.localNikeshDir + 'QTG_Candidates/', exist_ok = True)
fm_obj.downloadData(fm_obj.localGenomeFile)

fm_obj.readSampleDatabase()
fm_obj.readAlignmentDatabase()
all_samples = fm_obj.sample_dt[fm_obj.sample_dt.Ecogroup.isin(['Deep_Benthic','Shallow_Benthic'])].SampleID.to_list()
aligned_samples = fm_obj.alignment_dt[fm_obj.alignment_dt.SampleID.isin(all_samples)].SampleID.to_list()

dt = pd.read_csv('candidateQTNs_all.tsv', sep = '\t')
threshold = 10
num_parallel = 48
sv_vcf_file = fm_obj.localNikeshDir + 'QTG_Candidates/' + 'candidateQTNs_sv.vcf'
sv_norm_vcf_file = fm_obj.localNikeshDir + 'QTG_Candidates/' + 'candidateQTNs_sv.norm.vcf.gz'
lv_csv_file = fm_obj.localNikeshDir + 'QTG_Candidates/' + 'candidateQTNs_lv.csv'
print_sv(sv_vcf_file, dt[(dt.Alt.str.len() <= threshold) & (dt.Reference.str.len() <= threshold)])
normalize_sites(sv_vcf_file, fm_obj.localGenomeFile, sv_norm_vcf_file)
dt[(dt.Alt.str.len() > threshold) | (dt.Reference.str.len() > threshold)].to_csv(lv_csv_file)

commands = []
for i,sampleID in enumerate(aligned_samples):
    out_vcf = fm_obj.localNikeshDir + 'QTG_Candidates/' + sampleID + '_candidate_QTNs.vcf.gz'
    command = ['python','-m', 'unit_scripts.genotypeCandidates',sv_norm_vcf_file,lv_csv_file,out_vcf,'Mzebra_GT3_NCBI',sampleID]
    if i < num_parallel:
        commands.append(SimpleNamespace(sampleID=sampleID, command = command, process = subprocess.Popen(data.command, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)))
    else:
        commands.append(SimpleNamespace(sampleID=sampleID, command = command, process = None))
    
while commands:
    finished_processes = [x for x in commands if x.process is not None and x.process.poll() is not None]
    for data in finished_processes:
        print(f'..{data.sampleID} complete..', end = '')
        if data.process.returncode != 0:
            bad_samples.append(data.sampleID)
        command_list.remove(data)  # Remove finished process from monitoring list
        next_command = next((x for x in command_list if x.process is None),None)
        if next_command is not None:
            next_command.process = subprocess.Popen(next_command.command, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

print(bad_samples)
