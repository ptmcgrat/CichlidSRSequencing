import pysam, pdb, math, subprocess, os

from collections import defaultdict
import pandas as pd
from helper_modules.file_manager import FileManager as FM
from helper_modules.ChimericData import ChimericRead

def genotype_insertion(chrom, position, ref, alt, bamfiles, flanking = 10):
	ref_seq = refObj[chrom][position - flanking:position - 1] + ref + refObj[chrom][position:position + flanking]
	alt_seq = refObj[chrom][position - flanking:position - 1] + alt + refObj[chrom][position:position + flanking]
	if len(alt) > 100:
		third_seq = row.Alt

    ref_SSW = StripedSmithWaterman(ref_seq)
    alt_SSW = StripedSmithWaterman(alt_seq)
   	
   	for read in bam_obj.fetch('NC_135176.1', position - 200, position + 200):

    for read in unique_reads.values():
        if not read.is_secondary and read.mapq > 10:
            ref_score = max(ref_SSW(read.seq).optimal_alignment_score, ref_SSW(reverse_complement(read.seq)).optimal_alignment_score)
            alt_score = max(alt_SSW(read.seq).optimal_alignment_score, alt_SSW(reverse_complement(read.seq)).optimal_alignment_score)
            if ref_score - alt_score > 10:
                geno[0]+=1
            elif ref_score - alt_score < -10:
                geno[1]+=1
            else:
                geno[2]+=1

def genotype_normal(bam_obj, position, mtype):
	for puc in bam_obj.pileup('NC_135176.1', position - 3, position + 3):
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
	for aln in bam_obj.fetch('NC_135176.1',deletion[0] - 2500, deletion[0] + 2500):
		try:
			SA = aln.get_tag('SA')
			out = [aln.reference_name, aln.reference_start, aln.cigarstring, aln.template_length, SA]
			#print(out)
			if SA.split(',')[0]!='NC_135176.1':
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
fm_obj = FM(genome_version = 'Mzebra_GT3_NCBI')
fm_obj.downloadData(fm_obj.localGenomeFile)
refObj = pysam.FastaFile(fm_obj.localGenomeFile)
#fm_obj.readSampleDatabase()
#samples = os.listdir(fm_obj.localBamRefDir)
#samples = fm_obj.a_dt[fm_obj.a_dt.GenomeVersion == 'Mzebra_GT3'].SampleID.to_list()

flanking = 10

dt = pd.read_csv('candidateQTNs_all.tsv', sep = '\t')
for i,row in dt.iterrows():
	#print(row.Reference + ' ' + refObj[row.Chromosome][row.Position-1])
	if len(row.Alt) > 8:
		#Insertion
		ref_allele = refObj[row.Chromosome][row.Position - flanking:row.Position - 1] + row.Reference + refObj[row.Chromosome][row.Position:row.Position + flanking]
		alt_allele = refObj[row.Chromosome][row.Position - flanking:row.Position - 1] + row.Alt + refObj[row.Chromosome][row.Position:row.Position + flanking]
		if len(row.Alt) > 100:
			third_allele = row.Alt
		#print(refObj[row.Chromosome][row.Position - flanking:row.Position + flanking])
		#print(row.Reference + ' ' + row.Alt)
	elif len(row.Reference) > 8 and '-' not in row.Reference and row.Alt == '-':
		ref_allele = refObj[row.Chromosome][row.Position - flanking:row.Position + flanking + + len(row.Reference)]
		alt_allele = refObj[row.Chromosome][row.Position - flanking:row.Position - 1] + refObj[row.Chromosome][row.Position + len(row.Reference) - 1:row.Position + len(row.Reference) + flanking]
		print(ref_allele)
		print(alt_allele)
		print(row.Reference)
		pdb.set_trace()
		#Deletion
		continue
		print(row.Position)
pdb.set_trace()
dt['Type'] = dt.Info.str.split('TYPE=').str[1].str.split(',').str[0]
b_dt = pd.read_csv('LargeInsertionsPairedHits.csv', index_col = 0)
mapping = {'NC_036780.1': 'NC_135167.1', 'NC_036781.1': 'NC_135168.1', 'NC_036782.1': 'NC_135169.1', 'NC_036783.1': 'NC_135170.1', 'NC_036784.1': 'NC_135171.1', 'NC_036785.1': 'NC_135172.1', 'NC_036786.1': 'NC_135173.1', 'NC_036787.1': 'NC_135174.1', 'NC_036788.1': 'NC_135175.1', 'NC_036789.1': 'NC_135176.1', 'NC_036790.1': 'NC_135177.1', 'NC_036791.1': 'NC_135178.1', 'NC_036792.1': 'NC_135179.1', 'NC_036793.1': 'NC_135180.1', 'NC_036794.1': 'NC_135181.1', 'NC_036795.1': 'NC_135182.1', 'NC_036796.1': 'NC_135183.1', 'NC_036797.1': 'NC_135184.1', 'NC_036798.1': 'NC_135185.1', 'NC_036799.1': 'NC_135186.1', 'NC_036800.1': 'NC_135187.1', 'NC_036801.1': 'NC_135188.1', 'NW_020192332.1': 'NW_027490031.1', 'NW_020192379.1': 'NW_027490032.1', 'NW_020192400.1': 'NW_027490033.1', 'NW_020192408.1': 'NW_027490034.1', 'NW_020192422.1': 'NW_027490035.1', 'NW_020192424.1': 'NW_027490036.1', 'NW_020192441.1': 'NW_027490037.1', 'NW_020192451.1': 'NW_027490038.1', 'NW_020192466.1': 'NW_027490039.1', 'NW_020192561.1': 'NW_027490040.1', 'NW_020192592.1': 'NW_027490041.1', 'NW_020192594.1': 'NW_027490042.1', 'NW_020192618.1': 'NW_027490043.1', 'NW_020192670.1': 'NW_027490044.1', 'NW_020192695.1': 'NW_027490045.1', 'NW_020192700.1': 'NW_027490046.1', 'NW_020192775.1': 'NW_027490047.1', 'NW_020192777.1': 'NW_027490048.1', 'NW_020192802.1': 'NW_027490049.1', 'NW_020192833.1': 'NW_027490050.1', 'NW_020192874.1': 'NW_027490051.1', 'NW_020192879.1': 'NW_027490052.1', 'NW_020192943.1': 'NW_027490053.1', 'NW_020192953.1': 'NW_027490054.1', 'NW_020192955.1': 'NW_027490055.1', 'NW_020193025.1': 'NW_027490056.1', 'NW_020193076.1': 'NW_027490057.1', 'NW_020193110.1': 'NW_027490058.1', 'NW_020193115.1': 'NW_027490059.1', 'NW_020193238.1': 'NW_027490060.1', 'NW_020193277.1': 'NW_027490061.1', 'NW_020193282.1': 'NW_027490062.1', 'NW_020193366.1': 'NW_027490063.1', 'NW_020193479.1': 'NW_027490064.1', 'NW_020193668.1': 'NW_027490065.1', 'NW_020193671.1': 'NW_027490066.1', 'NW_020193692.1': 'NW_027490067.1', 'NW_020193834.1': 'NW_027490068.1', 'NW_020193902.1': 'NW_027490069.1', 'NW_020193974.1': 'NW_027490070.1'}
b_dt['PairChrom'] = b_dt.PairChrom.map(mapping)

out_dt = pd.DataFrame(index = dt.Position)

for sample in samples:
	try:
		organism = fm_obj.sample_dt[fm_obj.sample_dt.SampleID == sample].Species.values[0]
	except IndexError:
		organism = 'UnknownSpecies'
	print(sample + '\t' + organism)
	fm_obj.createSampleFiles(sample)
	#fm_obj.downloadData(fm_obj.localSampleBamDir)
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
			for aln in bam_obj_dis.fetch('NC_135176.1',position - 2000,position+2000):
				if not aln.is_proper_pair and position - aln.reference_start <= 700 and aln.reference_start - position <= 600:
					paired_pos = 10000*int(aln.next_reference_start/10000)
					paired_chr = aln.next_reference_name
					if ((b_dt.PairChrom == paired_chr) & (b_dt.BinnedPair == paired_pos)).sum() >= 1:
						hits+=1
			out_data.append('0,' + str(hits))
		else:
			print('Error')
	out_dt[sample] = out_data
	out_dt = out_dt.copy()
	#subprocess.run(['rm','-rf',fm_obj.localSampleBamDir])
	out_dt.to_csv('candidateQTNs_LakeMalawiGenotypes_Blumer.csv')
