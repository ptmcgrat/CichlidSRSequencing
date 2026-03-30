import pysam, pdb, math, subprocess

from collections import defaultdict
import pandas as pd

from helper_modules.file_manager import FileManager as FM
from helper_modules.ChimericData import ChimericRead

def genotype_normal(bam_obj, position, mtype):
	for puc in bam_obj.pileup('NC_036789.1', position - 3, position + 3):
	    # Check if the current column is the exact target position (important when fetching larger regions)
	    if puc.pos == position - 1:
	        
	        # Iterate over all reads in this specific pileup column
	        bases = ''
	        for pileupread in puc.pileups:
	            # Skip reads that are deleted or reference skips at this position
	            if pileupread.indel>0:
	            	bases += '-'
	            elif pileupread.is_del or pileupread.is_refskip:
	            	bases += '*'
	            else:
	                # The query_position is the 0-based index of the base in the read's sequence
	                base = pileupread.alignment.query_sequence[pileupread.query_position]
	                bases += base
	       	refs = sum([row.Reference.upper() == x for x in bases])
	       	if mtype == 'SUBSTITUTE':
	       		alts = sum([row.Alt.upper() == x for x in bases])
	       	elif mtype == 'DELETION':
	       		refs = sum([row.Reference[0].upper() == x for x in bases])
	       		alts = sum(['*' == x for x in bases])
	       	elif mtype == 'INSERTION':
	       		alts = sum(['-' == x for x in bases])

	       	return str(refs) + ',' + str(alts)

def genotype_longdeletion(bam_obj):
	deletion = (24870603,(9, 24870601, 'r', 9, 24872359, 'l', '', 'del'))
	deletion_dict = defaultdict(int)
	long_templates = 0
	for aln in bam_obj.fetch('NC_036789.1',deletion[0] - 2500, deletion[0] + 2500):
		try:
			SA = aln.get_tag('SA')
			out = [aln.reference_name, aln.reference_start, aln.cigarstring, aln.template_length, SA]
			#print(out)
			if SA.split(',')[0]!='NC_036789.1':
				continue

			aln.set_tag('SA',SA.replace(SA.split(',')[0], '9'))

			try:
				newRead = ChimericRead(aln)
			except TypeError:
				continue
			if newRead.data[7] == 'del':
				deletion_dict[newRead.data] += 1
		except KeyError:
			if aln.template_length > 1000 and aln.reference_start < deletion[0] and aln.next_reference_start > deletion[0]:
				long_templates += 1
	return '0,' + str(long_templates+ deletion_dict[deletion[1]])

deletion = (24870603,(9, 24870601, 'r', 9, 24872359, 'l', '', 'del'))
fm_obj = FM(genome_version = 'Mzebra_GT3')
samples = fm_obj.a_dt[fm_obj.a_dt.GenomeVersion == 'Mzebra_GT3'].SampleID.to_list()

dt = pd.read_csv('candidateQTNs_all.tsv', sep = '\t')
dt['Type'] = dt.Info.str.split('TYPE=').str[1].str.split(',').str[0]
b_dt = pd.read_csv('LargeInsertionsPairedHits.csv', index_col = 0)
out_dt = pd.DataFrame(index = dt.Position)

for sample in samples:
	try:
		organism = fm_obj.sample_dt[fm_obj.sample_dt.SampleID == sample].Species.values[0]
	except IndexError:
		organism = 'UnknownSpecies'
	print(sample + '\t' + organism)
	fm_obj.createSampleFiles(sample)
	fm_obj.downloadData(fm_obj.localSampleBamDir)
	bam_obj = pysam.AlignmentFile(fm_obj.localBamFile)
	bam_obj_dis = pysam.AlignmentFile(fm_obj.localDiscordantBamFile)

	out_data = []
	for i,row in dt.iterrows():
		position = row['Position']
		mtype = row['Type']
		if position not in b_dt.InsertionLocation.values and position != deletion[0]:
			out_data.append(genotype_normal(bam_obj, position, mtype))

		elif position == deletion[0]:
			out_data.append(genotype_longdeletion(bam_obj))

		elif position in b_dt.InsertionLocation.values:
			hits = 0
			for aln in bam_obj_dis.fetch('NC_036789.1',position - 2000,position+2000):
				if not aln.is_proper_pair and position - aln.reference_start <= 700 and aln.reference_start - position <= 600:
					paired_pos = 10000*int(aln.next_reference_start/10000)
					paired_chr = aln.next_reference_name
					if ((b_dt.PairChrom == paired_chr) & (b_dt.BinnedPair == paired_pos)).sum() >= 1:
						hits+=1
			out_data.append('0,' + str(hits))
		else:
			print('Error')
	out_dt[sample] = out_data

	subprocess.run(['rm','-rf',fm_obj.localSampleBamDir])
	out_dt.to_csv('candidateQTNs_LakeMalawiGenotypes.csv')
