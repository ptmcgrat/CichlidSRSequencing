from helper_modules.file_manager import FileManager as FM
import subprocess

fm_obj = FM()
fm_obj.downloadData(fm_obj.localSampleFile)
s_dt = pd.read_csv(fm_obj.localSampleFile)

fq1 = fm_obj.localReadsDir + 'PRJEB15289/ERR3634107_1.fastq.gz'
fq2 = fm_obj.localReadsDir + 'PRJEB15289/ERR3634107_2.fastq.gz'

fm_obj.downloadData(fq1)
fm_obj.downloadData(fq2)

subprocess.call(['gunzip', fq1])
subprocess.call(['gunzip', fq2])

subprocess.call(['head', '1600', fq1.replace('.gz','')], stdout = open(fm_obj.localReadsDir + 'PRJEB15289/Test_1.fastq', 'w'))
subprocess.call(['head', '1600', fq2.replace('.gz','')], stdout = open(fm_obj.localReadsDir + 'PRJEB15289/Test_2.fastq', 'w'))

pdb.set_trace()

subprocess.call(['gunzip', fm_obj.localReadsDir + 'PRJEB15289/Test_1.fastq'])
subprocess.call(['gunzip', fm_obj.localReadsDir + 'PRJEB15289/Test_2.fastq'])

