







bwa_command = ['bwa', 'mem', '-t', str(cpu_count()), '-R', run.RG_info, '-M', ref_file, run.fqfile1, run.fqfile2]
p1 = Popen(['samtools', 'view', '-bhS', '-@', str(cpu_count()), tfile1], stdout=PIPE)
p2 = Popen(['samtools', 'sort','-o', tfile2, '-@', str(cpu_count()), '-'], stdin = p1.stdout, stderr = FNULL)
p2.communicate()
