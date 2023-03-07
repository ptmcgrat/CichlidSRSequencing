import subprocess as sp
import shlex, pdb

regions  = ['pre_inversion', 'inversion1', 'inversion2', 'post_inversion']
processes1 = []
processes2 = []

for region in regions:
    # # Mzebra Server Code; use genomics env
    p1 = sp.Popen(shlex.split(f"gatk SelectVariants -V /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/pipeline_outputs/samples_filtered_master.vcf.gz -L /home/ad.gatech.edu/bio-mcgrath-dropbox/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/intervals_lg11/{region}.interval_list -O {region}.vcf"))
    p2 = sp.Popen(shlex.split(f"bgzip {region}.vcf"))
    processes1.append(p1)
    processes2.append(p2)
    if len(processes1) == 4 & len(processes2) == 4:
        for proc1 in processes1:
            proc1.communicate()
        for proc2 in processes2:
            proc2.communicate()
        processes1 = []
        processes2 = []
        

    # Local Code; use gatk_test env: ~/anaconda3/envs/gatk_test/bin/python3 make_lg11_intervals_files.py 
    # p1 = sp.Popen(shlex.split(f"gatk SelectVariants -V /Users/kmnike/CichlidSRSequencing/Test/samples_filtered_master.vcf.gz -L /Users/kmnike/CichlidSRSequencing/cichlid_sr_sequencing/analysis_pipeline/intervals_lg11/{region}.interval_list -O {region}.vcf"))
    # p2 = sp.Popen(shlex.split(f"bgzip {region}.vcf"))
    # processes1.append(p1)
    # processes2.append(p2)
    # if len(processes1) == 4 and len(processes2) == 4:
    #     for proc1 in processes1:
    #         proc1.communicate()
    #     for proc2 in processes2:
    #         proc2.communicate()
    #     processes1 = []
    #     processes2 = []