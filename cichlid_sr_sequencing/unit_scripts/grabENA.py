import subprocess, argparse, datetime, os, pdb,sys
from helper_modules.file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage = 'This script grabs the ENA data for a run and uploads it to dropbox')
parser.add_argument('RunID', type = str, help = 'File containing information on each run')
parser.add_argument('ENA_fq1', type = str, help = 'File containing information on each run')
parser.add_argument('ENA_fq2', type = str, help = 'File containing information on each run')
parser.add_argument('OutputBam', type = str, help = 'Location of output bam files')
parser.add_argument('Temp_directory', type = str, help = 'Directory to temporarily hold datafiles')
parser.add_argument('SampleName', type = str, help = 'Info for creating read group tag for bam file')
parser.add_argument('LibraryName', type = str, help = 'Info for creating read group tag for bam file')
parser.add_argument('Platform', type = str, help = 'Info for creating read group tag for bam file')
parser.add_argument('-t', '--TestData', action = 'store_true', help = 'Use this flag if you want to create a small test file (1000 reads) instead of the entire read set')

args = parser.parse_args()

fm_obj = FM()

target_directory = args.Temp_directory
local_fq1 = args.Temp_directory + args.RunID + '_1.fastq.gz'
local_fq2 = args.Temp_directory + args.RunID + '_2.fastq.gz'
temp_bam_file = args.Temp_directory + args.RunID + '_temp.bam'

#target_directory = args.Local_fq1.replace(args.Local_fq1.split('/')[-1],'')
print('  Fastq files acsping for ' + args.RunID + ', Time:' + str(datetime.datetime.now()))
for i in range(3):
	pdb.set_trace()
	output = subprocess.run(['ascp', '-QT', '-l', '1000m', '-P', '33001', '-i', os.getenv('HOME') + '/anaconda3/envs/CichlidSRSequencing/etc/asperaweb_id_dsa.openssh', args.ENA_fq1.replace('ftp.sra.ebi.ac.uk/','era-fasp@fasp.sra.ebi.ac.uk:'),target_directory], capture_output = True)
	if output.returncode == 0:
		print(output.stderr.decode('utf-8'))
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
command = ['gatk', 'FastqToSam', '--FASTQ', local_fq1, '--FASTQ2', local_fq2, '--READ_GROUP_NAME', args.RunID]
command += ['--OUTPUT', temp_bam_file, '--SAMPLE_NAME', args.SampleName, '--LIBRARY_NAME', args.LibraryName, '--PLATFORM', args.Platform]
output1 = subprocess.run(command, capture_output = True)
if output1.returncode != 0:
	with open(args.outputBam + '.FastQToSamErrors.txt', 'w') as f:
		print(output1.stderr.decode('utf-8'), file = f)
	fm_obj.uploadData(args.OutputBam + '.FastQToSamErrors.txt')
	sys.exit()

# Mark illumina adapters
print('  Marking Illumina adapters')
command = ['gatk', 'MarkIlluminaAdapters', '-I', temp_bam_file, '-O', args.OutputBam, '-M', args.OutputBam + '.metrics.txt', '--TMP_DIR', args.Temp_directory]
output2 = subprocess.run(command, capture_output = True)
if output2.returncode != 0:
	with open(args.outputBam + '.MarkIlluminaErrors.txt', 'w') as f:
		print(output2.stderr.decode('utf-8'), file = f)
	fm_obj.uploadData(args.OutputBam + '.MarkIlluminaErrors.txt')
	sys.exit()

# Upload data to dropbox
print('  Uploading uBam files for ' + args.RunID + ', Time:' + str(datetime.datetime.now()))
fm_obj.uploadData(args.OutputBam)
fm_obj.uploadData(args.OutputBam + '.metrics.txt')
print('  Finished for ' + args.RunID + ', Time:' + str(datetime.datetime.now()))

# Remove files that were created
subprocess.run(['rm', local_fq1, local_fq2, temp_bam_file, args.OutputBam, args.OutputBam + '.metrics.txt'])
