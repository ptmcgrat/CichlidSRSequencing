import subprocess, argparse, datetime, os, pdb,sys, pysam
from helper_modules.file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage = 'This script grabs the ENA data for a run and uploads it to dropbox')
parser.add_argument('RunID', type = str, help = 'File containing information on each run')
parser.add_argument('fq1', type = str, help = 'File containing information on each run')
parser.add_argument('fq2', type = str, help = 'File containing information on each run')
parser.add_argument('OutputBam', type = str, help = 'Location of output bam files')
parser.add_argument('Temp_directory', type = str, help = 'Directory to temporarily hold datafiles')
parser.add_argument('SampleName', type = str, help = 'Info for creating read group tag for bam file')
parser.add_argument('LibraryName', type = str, help = 'Info for creating read group tag for bam file')
parser.add_argument('Platform', type = str, help = 'Info for creating read group tag for bam file')
parser.add_argument('LibraryLayout', type = str, help = 'Paired end or single reads')
parser.add_argument('-t', '--TestData', action = 'store_true', help = 'Use this flag if you want to create a small test file (1000 reads) instead of the entire read set')
parser.add_argument('-l', '--Local', action = 'store_true', help = 'Use this flag if the read data is on Dropbox')

args = parser.parse_args()

fm_obj = FM()

target_directory = args.Temp_directory

if args.Local:
	local_fq1 = args.fq1
	local_fq2 = args.fq2
else:
	local_fq1 = args.Temp_directory + args.RunID + '_1.fastq.gz'
	local_fq2 = args.Temp_directory + args.RunID + '_2.fastq.gz'
temp_bam_file = args.Temp_directory + args.RunID + '_temp.bam'

if args.Local:
	fm_obj.downloadData(args.fq1)
	if args.LibraryLayout == 'PAIRED':
		fm_obj.downloadData(args.fq2)
else:

	#target_directory = args.Local_fq1.replace(args.Local_fq1.split('/')[-1],'')
	print('  Fastq files acsping for ' + args.RunID + ', Time:' + str(datetime.datetime.now()))
	for i in range(3):
		output = subprocess.run(['ascp', '-QT', '-l', '1000m', '-P', '33001', '-i', os.getenv('HOME') + '/anaconda3/envs/CichlidSRSequencing/etc/asperaweb_id_dsa.openssh', args.ENA_fq1.replace('ftp.sra.ebi.ac.uk/','era-fasp@fasp.sra.ebi.ac.uk:'),target_directory], capture_output = True)
		if output.returncode == 0:
			break
		elif i == 2:
			sys.exit()
		print('Redownloading ' + args.RunID + ' try ' + str(i + 1))

	for i in range(3):
		output = subprocess.run(['ascp', '-QT', '-l', '1000m', '-P', '33001', '-i', os.getenv('HOME') + '/anaconda3/envs/CichlidSRSequencing/etc/asperaweb_id_dsa.openssh', args.ENA_fq2.replace('ftp.sra.ebi.ac.uk/','era-fasp@fasp.sra.ebi.ac.uk:'),target_directory], capture_output = True)
		if output.returncode == 0:
			break
		elif i == 2:
			sys.exit()
		print('Redownloading ' + args.RunID + ' try ' + str(i + 1))

# Convert fastq files to unmapped bam
print('  Converting fastq files to uBam file')


# Quality control fastq files
f1 = pysam.FastqFile(local_fq1)
if args.LibraryLayout == 'PAIRED':
	f2 = pysam.FastqFile(local_fq2)

fixed_fq1 = local_fq1.replace(local_fq1.split('/')[-1],'fixed_' + local_fq1.split('/')[-1]).replace('.gz','')
if args.LibraryLayout == 'PAIRED':
	fixed_fq2 = local_fq2.replace(local_fq2.split('/')[-1],'fixed_' + local_fq2.split('/')[-1]).replace('.gz','')

if args.LibraryLayout == 'PAIRED':

	with open(fixed_fq1, 'w') as outfq1, open(fixed_fq2, 'w') as outfq2:
		for r1,r2 in zip(f1,f2):
			if r1.sequence == '' or r2.sequence == '':
				continue
			else:
				outfq1.write('@' + r1.name + ' 1:N:0:2\n' + r1.sequence + '\n+\n' + r1.quality + '\n')
				outfq2.write('@' + r2.name + ' 2:N:0:2\n' + r2.sequence + '\n+\n' + r2.quality + '\n')

else:
	with open(fixed_fq1, 'w') as outfq1:
		for r1 in f1,f:
			if r1.sequence == '':
				continue
			else:
				outfq1.write('@' + r1.name + ' 1:N:0:2\n' + r1.sequence + '\n+\n' + r1.quality + '\n')

if args.LibraryLayout == 'PAIRED':
	command = ['gatk', 'FastqToSam', '--FASTQ', fixed_fq1, '--FASTQ2', fixed_fq2, '--READ_GROUP_NAME', args.RunID, '--TMP_DIR', args.Temp_directory]
	command += ['--OUTPUT', temp_bam_file, '--SAMPLE_NAME', args.SampleName, '--LIBRARY_NAME', args.LibraryName, '--PLATFORM', args.Platform]
else:
	command = ['gatk', 'FastqToSam', '--FASTQ', fixed_fq1, '--READ_GROUP_NAME', args.RunID, '--TMP_DIR', args.Temp_directory]
	command += ['--OUTPUT', temp_bam_file, '--SAMPLE_NAME', args.SampleName, '--LIBRARY_NAME', args.LibraryName, '--PLATFORM', args.Platform]

output1 = subprocess.run(command, capture_output = True)
if output1.returncode != 0:
	with open(args.OutputBam + '.FastQToSamErrors.txt', 'w') as f:
		print(output1.stderr.decode('utf-8'), file = f)
	fm_obj.uploadData(args.OutputBam + '.FastQToSamErrors.txt')
	sys.exit()

# Mark illumina adapters
if args.Platform == 'ILLUMINA':
	print('  Marking Illumina adapters')
	command = ['gatk', 'MarkIlluminaAdapters', '-I', temp_bam_file, '-O', args.OutputBam, '-M', args.OutputBam + '.metrics.txt', '--TMP_DIR', args.Temp_directory]
	output2 = subprocess.run(command, capture_output = True)
	if output2.returncode != 0:
		with open(args.OutputBam + '.MarkIlluminaErrors.txt', 'w') as f:
			print(output2.stderr.decode('utf-8'), file = f)
		fm_obj.uploadData(args.OutputBam + '.MarkIlluminaErrors.txt')
		sys.exit()

# Upload data to dropbox
print('  Uploading uBam files for ' + args.RunID + ', Time:' + str(datetime.datetime.now()))
fm_obj.uploadData(args.OutputBam)
fm_obj.uploadData(args.OutputBam + '.metrics.txt')
print('  Finished for ' + args.RunID + ', Time:' + str(datetime.datetime.now()))

# Remove files that were created
if args.LibraryLayout == 'PAIRED':
	subprocess.run(['rm', local_fq1, local_fq2, fixed_fq1, fixed_fq2, temp_bam_file, args.OutputBam, args.OutputBam + '.metrics.txt'])
else:
	subprocess.run(['rm', local_fq1, fixed_fq1, temp_bam_file, args.OutputBam, args.OutputBam + '.metrics.txt'])
