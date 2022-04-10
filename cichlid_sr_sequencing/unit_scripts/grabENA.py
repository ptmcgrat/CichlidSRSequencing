import subprocess, argparse, datetime, os, pdb,sys
from helper_modules.file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage = 'This script grabs the ENA data for a run and uploads it to dropbox')
parser.add_argument('RunID', type = str, help = 'File containing information on each run')
parser.add_argument('ENA_fq1', type = str, help = 'File containing information on each run')
parser.add_argument('ENA_fq2', type = str, help = 'File containing information on each run')
parser.add_argument('Local_fq1', type = str, help = 'File containing information on each run')
parser.add_argument('Local_fq2', type = str, help = 'File containing information on each run')

args = parser.parse_args()

fm_obj = FM()

target_directory = args.Local_fq1.replace(args.Local_fq1.split('/')[-1],'')
print('  Fastq files acsping for ' + args.RunID + ', Time:' + str(datetime.datetime.now()))

for i in range(3):
	output = subprocess.run(['ascp', '-QT', '-l', '1000m', '-P', '33001', '-i', os.getenv('HOME') + '/anaconda3/envs/CichlidSRSequencing/etc/asperaweb_id_dsa.openssh', args.ENA_fq1.replace('ftp.sra.ebi.ac.uk/','era-fasp@fasp.sra.ebi.ac.uk:'),target_directory], capture_output = True)
	if output.returncode == 0:
		break
	elif i == 2:
		sys.exit()

for i in range(3):
	output = subprocess.run(['ascp', '-QT', '-l', '1000m', '-P', '33001', '-i', os.getenv('HOME') + '/anaconda3/envs/CichlidSRSequencing/etc/asperaweb_id_dsa.openssh', args.ENA_fq2.replace('ftp.sra.ebi.ac.uk/','era-fasp@fasp.sra.ebi.ac.uk:'),target_directory], capture_output = True)
	if output.returncode == 0:
		break
	elif i == 2:
		sys.exit()
	print('Redownloading ' + args.RunID + ' try ' + str(i + 1))

print('  Rcloning files for ' + args.RunID + ', Time:' + str(datetime.datetime.now()))
fm_obj.uploadData(args.Local_fq1)
fm_obj.uploadData(args.Local_fq2)
print('  Finished for' + args.RunID + ', Time:' + str(datetime.datetime.now()))

subprocess.run(['rm', args.Local_fq1])
subprocess.run(['rm', args.Local_fq2])
