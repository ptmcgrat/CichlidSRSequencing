import subprocess as sp
import os
import shlex

files = os.listdir()
processes = []
rc = []
for file in files:
    if file.startswith('SAM'):
        command = shlex.split(f"/Users/kmnike/bin/gatk-4.2.6.1/gatk ValidateVariants -R ~/Data/CichlidSequencingData/Genome/GCF_000238955.4_M_zebra_UMD2a_genomic.fna -V {file + '/' + file + '.g.vcf.gz'}")
        p = sp.Popen(command, stdout=sp.DEVNULL, stderr=sp.STDOUT)
        processes.append(p)
    if len(processes) == 3:
        for p in processes:
            p.communicate()
            rc.append(p.returncode)
        processes = []
print(rc)