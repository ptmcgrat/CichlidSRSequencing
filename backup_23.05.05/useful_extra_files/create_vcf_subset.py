import subprocess, shlex

subprocess.run(shlex.split('bcftools view -R regions_file.txt -O z -o small.vcf.gz'))