import subprocess as sp
import shlex

commands = [shlex.split('python3 timer.py 3'), shlex.split('python3 timer.py 6'), shlex.split('python3 timer.py 2')]
processes = []
for command in commands:
	p = sp.Popen(command)#, stdout=sp.PIPE, stderr=sp.PIPE)
	processes.append(p)
	# print('processes list:', processes)
	
	if len(processes) == 3:
		for p in processes:
			p.communicate()
		processes = []



"""
processes = []
for contig in fasta_obj.references:
	#timer.start('Calling SNVs for ' + str(len(bamfiles)) + ' bamfiles.')

	p1 = subprocess.Popen(['bcftools', 'mpileup', '-r', contig, '-C', '50', '-pm2', '-F', '0.2', '-f', fm_obj.localGenomeFile] + bamfiles, stdout = subprocess.PIPE)
	processes.append(subprocess.Popen(['bcftools', 'call', '-vmO', 'v', '-f', 'GQ', '-o', contig + '.vcf'], stdin = p1.stdout))

	if len(processes) > 23:
		for p in processes:	
			p.communicate()
		processes = []
	#timer.stop()



f_list = glob.glob('./*bz2')
cmds_list = [['./bunzip2_file.py', file_name] for file_name in f_list]
procs_list = [Popen(cmd, stdout=PIPE, stderr=PIPE) for cmd in cmds_list]
for proc in procs_list:
	proc.wait()

"""