import subprocess as sp, argparse, shlex, pysam, os
# from helper_modules.file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage = 'This script will concatenate VCF ouputs from the CallSmallSNVs.py pipeline into a master VCF file')
parser.add_argument('Genome', type = str, help = 'Version of the genome use')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
# fm_obj = FM(args.Genome)

#### use rclone lsf to create a list of existing vcf files 
outputs = sp.run(shlex.split("rclone lsf ptm_dropbox:BioSci-McGrath/Apps/CichlidSequencingData/Outputs/Mzebra_UMD2a/Outputs576Cohort/ --include '*.vcf'"), stdout=sp.PIPE, encoding='utf-8').stdout.splitlines()
path_to_outputs = 'ptm_dropbox:BioSci-McGrath/Apps/CichlidSequencingData/Outputs/Mzebra_UMD2a/Outputs576Cohort/'

#### check if an output dir exists. If not, make it.
try:
    os.makedirs(os.getcwd() + '/vcf_concat_output/')
except FileExistsError:
    pass
output_dir = os.getcwd() + '/vcf_concat_output/'

#### make a Temp dir to temporrily store the VCF files as they're concatenated
try:
    os.makedirs(os.getcwd() + '/Temp/')
except FileExistsError:
    pass
temp_dir = os.getcwd() + '/Temp/'

#### Quick check to see if a master vcf file has already been created. 
file_num = 0
for file in os.listdir(output_dir):
    file_num += 1

#### if no master file exists, run through the code normally and name the file "master_file.vcf"
if file_num == 0:
    with open(f"{output_dir}master_file.vcf", 'w+') as f1: #open the master list file if the vcf_concat_output dir is empty
        for file in outputs: # the output file names were defined above. Iterate through each file
            sp.run(shlex.split(f"rclone copy {path_to_outputs + file} {temp_dir} -P")) # download the file to a temp dir
            file_to_write = temp_dir + file 
            if os.path.getsize(f"{output_dir}master_file.vcf") == 0: # if the master file is empty, take the file and write the whole contents. This will write the header
                f1.write(sp.check_output(shlex.split(f"cat {file_to_write}"), encoding='utf-8') + '\n') #the newline is needed at the end
            else: # if file has contents, open the file in read mode then go through each line. If it doesnt start with a "#" then write it to the master file.
                with open(temp_dir + file, 'r') as f2:
                    for line in f2:
                        if not line.startswith('#'):
                            f1.write(line)
                    f1.write('\n') # last newline needed to not merge lines between files together. 
            sp.run(shlex.split(f"rm {file_to_write}"))
            sp.run(shlex.split(f"rm  -r {temp_dir}"))
else: # if the master file exists, create a new file so that the data in previous iterations will be preserved. 
    with open(f"{output_dir}master_file{file_num}.vcf", 'w+') as f1: #open the master list file if the vcf_concat_output dir is empty
            for file in outputs: # the output file names were defined above. Iterate through each file
                sp.run(shlex.split(f"rclone copy {path_to_outputs + file} {temp_dir} -P")) # download the file to a temp dir
                file_to_write = temp_dir + file 
                if os.path.getsize(f"{output_dir}master_file{file_num}.vcf") == 0: # if the master file is empty, take the file and write the whole contents. This will write the header
                    f1.write(sp.check_output(shlex.split(f"cat {file_to_write}"), encoding='utf-8') + '\n') #the newline is needed at the end
                else: # if file has contents, open the file in read mode then go through each line. If it doesnt start with a "#" then write it to the master file. 
                    with open(temp_dir + file, 'r') as f2:
                        for line in f2:
                            if not line.startswith('#'):
                                f1.write(line)
                        f1.write('\n') # last newline needed to not merge lines between files together.
                sp.run(shlex.split(f"rm {file_to_write}"))
                sp.run(shlex.split(f"rm  -r {temp_dir}"))


print("CONCATENATION COMPLETE")

