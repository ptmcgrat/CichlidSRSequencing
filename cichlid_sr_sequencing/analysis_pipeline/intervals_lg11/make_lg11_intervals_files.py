import subprocess as sp
import shlex

regions  = ['pre_inversion', 'inversion1', 'inter_inversion', 'inversion2', 'post_inversion']

for region in regions:
    # sp.run(shlex.split(f"gatk SelectVariants -V /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/pipeline_outputs/samples_filtered_master.vcf.gz -L /home/ad.gatech.edu/bio-mcgrath-dropbox/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/intervals_lg11/{region}.interval_list -O {region}.vcf"))
    sp.run(shlex.split(f"gatk SelectVariants -V /Users/kmnike/CichlidSRSequencing/Test/samples_filtered_master.vcf.gz -L /Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/intervals_lg11/{region}.interval_list -O {region}.vcf"))
    sp.run(shlex.split(f"bgzip {region}.vcf"))
    