import subprocess, argparse, contextlib, urllib, shutil
from helper_modules.file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage = 'This script grabs the ENA data for a run and uploads it to dropbox')
parser.add_argument('RunID', type = str, help = 'File containing information on each run')
parser.add_argument('ENA_fq1', type = str, help = 'File containing information on each run')
parser.add_argument('ENA_fq2', type = str, help = 'File containing information on each run')
parser.add_argument('Local_fq1', type = str, help = 'File containing information on each run')
parser.add_argument('Local_fq2', type = str, help = 'File containing information on each run')

args = parser.parse_args()

fm_obj = FM()

print('Fastq files ftping for ' + args.RunID + ', Time:' + str(datetime.datetime.now()))
with contextlib.closing(urllib.request.urlopen(args.ENA_fq1)) as r:
	with open(args.Local_fq1, 'wb') as f:
		shutil.copyfileobj(r, f)

with contextlib.closing(urllib.request.urlopen(args.ENA_fq2)) as r:
	with open(args.Local_fq2, 'wb') as f:
		shutil.copyfileobj(r, f)

print('Rcloning files for ' + args.RunID + ', Time:' + str(datetime.datetime.now()))
fm_obj.uploadData(fargs.Local_fq1)
fm_obj.uploadData(fargs.Local_fq2)
print('Finished for' + args.RunID + ', Time:' + str(datetime.datetime.now()))

subprocess.run(['rm', args.Local_fq1])
subprocess.run(['rm', args.Local_fq2])
