from helper_modules.file_manager import FileManager as FM
import datetime, pdb, subprocess
import pandas as pd
from helper_modules.Timer import Timer
from multiprocessing import cpu_count


fm_obj = FM('Mzebra_GT3')
pdb.set_trace()
#fm_obj.downloadData(fm_obj.localGenomeDir)
#subprocess.run(['bwa-mem2','index', fm_obj.localGenomeFile])
fm_obj.createSampleFiles('SAMEA117806073')
timer = Timer()
for uBam_file in fm_obj.localRawBamFiles:
	fm_obj.downloadData(uBam_file)
	"""
	timer.start('Create fastq files')
	command1 = ['gatk', 'SamToFastq', '-I', uBam_file, '--FASTQ', uBam_file.replace('.bam','.fastq'), '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2']
	command1 += ['--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', fm_obj.localSampleTempDir]
	output = subprocess.run(command1, capture_output = True)
	timer.stop()
	"""
	timer.start('Do it the new way')
	command1 = ['gatk', 'SamToFastq', '-I', uBam_file, '--FASTQ', '/dev/stdout', '--CLIPPING_ATTRIBUTE', 'XT', '--CLIPPING_ACTION', '2']
	command1 += ['--INTERLEAVE', 'true', '--NON_PF', 'true', '--TMP_DIR', fm_obj.localSampleTempDir]
	command2 = ['minimap2', '-x','sr', '-a', '-t', str(cpu_count()), fm_obj.localGenomeFile, '-']
	command3 = ['gatk', 'MergeBamAlignment', '-R', fm_obj.localGenomeFile, '--UNMAPPED_BAM', uBam_file, '--ALIGNED_BAM', '/dev/stdin']
	command3 += ['-O', 'newway.bam', '--ADD_MATE_CIGAR', 'true', '--CLIP_ADAPTERS', 'false', '--CLIP_OVERLAPPING_READS', 'true']
	command3 += ['--INCLUDE_SECONDARY_ALIGNMENTS', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1', '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant']
	command3 += ['--ATTRIBUTES_TO_RETAIN', 'XS', '--TMP_DIR', fm_obj.localSampleTempDir]

	error_file = open(fm_obj.localSampleTempDir + 'Alignment_errors.txt', 'w')
	p1 = subprocess.Popen(command1, stdout=subprocess.PIPE, stderr = error_file)
	p2 = subprocess.Popen(command2, stdin = p1.stdout, stdout = subprocess.PIPE, stderr = error_file)
	p1.stdout.close()
	p3 = subprocess.Popen(command3, stdin = p2.stdout, stderr = error_file, stdout = subprocess.DEVNULL)
	p2.stdout.close()
	output = p3.communicate()
	timer.stop()
	timer.start('Just mem')
	command1 = ['bwa', 'mem', '-t', str(cpu_count()), '-M', '-p', fm_obj.localGenomeFile, uBam_file.replace('.bam','.fastq'), '-o','mem_fq.sam']
	output1 = subprocess.run(command1, capture_output = True)
	timer.stop()
	timer.start('Just mem2')
	command1 = ['bwa-mem2', 'mem', '-t', str(cpu_count()), '-M', '-p', fm_obj.localGenomeFile, uBam_file.replace('.bam','.fastq'), '-o','mem2_fq.sam']
	output2 = subprocess.run(command1, capture_output = True)
	timer.stop()
	
	timer.start('Just minimap2')
	command1 = ['minimap2', '-x','sr', '-a', '-t', str(cpu_count()), fm_obj.localGenomeFile, uBam_file.replace('.bam','.fastq'), '-o','mmap2_fq.sam']
	output3 = subprocess.run(command1, capture_output = True)
	timer.stop()
	pdb.set_trace()


