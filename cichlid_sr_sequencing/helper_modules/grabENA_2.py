import subprocess, argparse, datetime, os, pdb,sys, pysam
from file_manager import FileManager as FM
import pandas as pd

parser = argparse.ArgumentParser(usage = 'This script grabs the ENA data for a run and uploads it to dropbox')
parser.add_argument('RunInfoFile', type = str, help = 'File containing info on the runs you would like to download')
parser.add_argument('RunID', type = str, help = 'RunID to download')
parser.add_argument('OutputBam', type = str, help = 'Location of output bam files')
parser.add_argument('Temp_directory', type = str, help = 'Directory to temporarily hold datafiles')
args = parser.parse_args()

fm_obj = FM()

target_directory = args.Temp_directory
dt = pd.read_csv(args.RunInfoFile)

for index, row in dt.iterrows():
	if row.RunID == args.RunID:
		asp_fq1 = row.FastqAspera.split(';')[0]
		asp_fq2 = row.FastqAspera.split(';')[1]
		local_fq1 = args.Temp_directory + asp_fq1.split('/')[-1]
		local_fq2 = args.Temp_directory + asp_fq2.split('/')[-1]
		break

temp_bam_file = args.Temp_directory + args.RunID + '_temp.bam'
print('  Processing files for ' + args.RunID + ', Time:' + str(datetime.datetime.now()))	
for i in range(3):
	output = subprocess.run(['ascp', '-QT', '-k','3', '-l', '1000m', '-P', '33001', '-i', os.getenv('HOME') + '/miniconda3/envs/phase/etc/asperaweb_id_dsa.openssh','era-fasp@' + asp_fq1,args.Temp_directory], capture_output = True)
	if output.returncode == 0:
		print(['md5sum',local_fq1])
		md5_output = subprocess.run(['md5sum',local_fq1], capture_output = True)
		
		if md5_output.returncode != 0:
			pdb.set_trace()
		if md5_output.stdout.decode().split()[0] == row.FastqMd5.split(';')[0]:
			print(output.stdout.decode())
			break
		else:
			print('Redownloading')
			continue

	if i == 2:
		pdb.set_trace()

for i in range(3):
	output = subprocess.run(['ascp', '-QT', '-k','3','-l', '1000m', '-P', '33001', '-i', os.getenv('HOME') + '/miniconda3/envs/phase/etc/asperaweb_id_dsa.openssh','era-fasp@' + asp_fq2,args.Temp_directory], capture_output = True)
	if output.returncode == 0:
		print(['md5sum',local_fq2])
		md5_output = subprocess.run(['md5sum',local_fq2], capture_output = True)
		if md5_output.returncode != 0:
			pdb.set_trace()

		if md5_output.stdout.decode().split()[0] == row.FastqMd5.split(';')[1]:
			print(output.stdout.decode())
			break
		else:
			print('Redownloading')
			continue
	if i == 2:
		pdb.set_trace()


# Quality control fastq files
f1 = pysam.FastqFile(local_fq1)
f2 = pysam.FastqFile(local_fq2)

fixed_fq1 = local_fq1.replace(local_fq1.split('/')[-1],'fixed_' + local_fq1.split('/')[-1]).replace('.gz','')
fixed_fq2 = local_fq2.replace(local_fq2.split('/')[-1],'fixed_' + local_fq2.split('/')[-1]).replace('.gz','')


with open(fixed_fq1, 'w') as outfq1, open(fixed_fq2, 'w') as outfq2:
	for r1,r2 in zip(f1,f2):
		if r1.name == '' or r2.name == '' or r1.sequence == '' or r2.sequence == '' or r1.quality == '' or r2.quality == '':
			continue
		elif r1.name is None or r2.name is None or r1.sequence is None or r2.sequence is None or r1.quality is None or r2.quality is None:
			continue
		else:
			outfq1.write('@' + r1.name + ' 1:N:0:2\n' + r1.sequence + '\n+\n' + r1.quality + '\n')
			outfq2.write('@' + r2.name + ' 2:N:0:2\n' + r2.sequence + '\n+\n' + r2.quality + '\n')


command = ['gatk', 'FastqToSam', '--FASTQ', fixed_fq1, '--FASTQ2', fixed_fq2, '--READ_GROUP_NAME', args.RunID, '--TMP_DIR', args.Temp_directory]
command += ['--OUTPUT', temp_bam_file, '--SAMPLE_NAME', row.SampleID, '--LIBRARY_NAME', row.LibraryID, '--PLATFORM', row.Platform]

output1 = subprocess.run(command, capture_output = True)
if output1.returncode != 0:
	with open(args.OutputBam + '.FastQToSamErrors.txt', 'w') as f:
		print(output1.stderr.decode('utf-8'), file = f)
	#fm_obj.uploadData(args.OutputBam + '.FastQToSamErrors.txt')
	print('  Failed converting to Sam ' + args.RunID + ', Time:' + str(datetime.datetime.now()))	
	sys.exit()

# Mark illumina adapters
if row.Platform == 'ILLUMINA' or 'ELEMENT' in row.Platform:
	#print('  Marking Illumina adapters')
	command = ['gatk', 'MarkIlluminaAdapters', '-I', temp_bam_file, '-O', args.OutputBam, '-M', args.OutputBam + '.metrics.txt', '--TMP_DIR', args.Temp_directory]
	output2 = subprocess.run(command, capture_output = True)
	if output2.returncode != 0:
		with open(args.OutputBam + '.MarkIlluminaErrors.txt', 'w') as f:
			print(output2.stderr.decode('utf-8'), file = f)
		#fm_obj.uploadData(args.OutputBam + '.MarkIlluminaErrors.txt')
		print('  Failed converting to Sam ' + args.RunID + ', Time:' + str(datetime.datetime.now()))	
		sys.exit()

# Upload data to dropbox
#print('  Uploading uBam files for ' + args.RunID + ', Time:' + str(datetime.datetime.now()))
fm_obj.uploadData(args.OutputBam)
fm_obj.uploadData(args.OutputBam + '.metrics.txt')
print('  Successfully processed ' + args.RunID + ', Time:' + str(datetime.datetime.now()))

# Remove files that were created
if args.LibraryLayout == 'PAIRED':
	subprocess.run(['rm', local_fq1, local_fq2, fixed_fq1, fixed_fq2, temp_bam_file, args.OutputBam, args.OutputBam + '.metrics.txt'])
else:
	subprocess.run(['rm', local_fq1, fixed_fq1, temp_bam_file, args.OutputBam, args.OutputBam + '.metrics.txt'])
