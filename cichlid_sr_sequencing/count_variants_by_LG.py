import subprocess as sp, pysam
from helper_modules.file_manager import FileManager as FM

fm_obj = FM('Mzebra_UMD2a')

fasta_obj = pysam.FastaFile(fm_obj.localGenomeFile)
# below line defines the LG names for the first 22 LGs in the genome.
contigs = fasta_obj.references[0:22]

variant_list = {}
for i in contigs:
    print(f"counting variants for {i}")
    # cmd = f"zgrep {i} /Users/kmnike/Data/CichlidSequencingData/Outputs/small_lg1-22_master_file.vcf.gz | wc -l"
    cmd = f"zgrep {i} /home/ad.gatech.edu/bio-mcgrath-dropbox/Data/CichlidSequencingData/Outputs/vcf_concat_output/original_data/lg1-22_master_file.vcf.gz | wc -l"
    ps = sp.Popen(cmd, shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
    output = int(ps.communicate()[0].decode('utf-8').strip())
    variant_list[i] = output
print(variant_list)

for key in variant_list:
    print(f"{key}\t{variant_list[key]}")
