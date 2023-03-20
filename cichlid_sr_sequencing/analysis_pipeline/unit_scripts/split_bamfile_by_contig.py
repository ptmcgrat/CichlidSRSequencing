import pysam, argparse, pdb

# Need to make SampleIDs and ProjectIDs mutually exclusive
parser = argparse.ArgumentParser(usage = 'This script will split bamfiles ')
parser.add_argument('Bamfile', type = str, help = 'Version of the genome to align to')
parser.add_argument('Contig', type = str, help = 'Version of the genome to align to')
args = parser.parse_args()

align_file = pysam.AlignmentFile(args.Bamfile) 

unmapped = pysam.AlignmentFile(args.Bamfile.replace('bam', args.Contig + '.unmapped.bam'), mode = 'wb', template = align_file)
discordant = pysam.AlignmentFile(args.Bamfile.replace('bam', args.Contig + '.discordant.bam'), mode = 'wb', template = align_file)
inversion = pysam.AlignmentFile(args.Bamfile.replace('bam', args.Contig + '.inversion.bam'), mode = 'wb', template = align_file)
duplication = pysam.AlignmentFile(args.Bamfile.replace('bam', args.Contig + '.duplication.bam'), mode = 'wb', template = align_file)
clipped = pysam.AlignmentFile(args.Bamfile.replace('bam', args.Contig + '.clipped.bam'), mode = 'wb', template = align_file)
chimeric = pysam.AlignmentFile(args.Bamfile.replace('bam', args.Contig + '.chimeric.bam'), mode = 'wb', template = align_file)

for read in align_file.fetch(contig = args.Contig):
	total_read_quality = sum([ord(x) - 33 for x in read.qual])/len(read.qual)
	if read.is_duplicate:
		continue
	elif read.is_unmapped:
		if total_read_quality > 25:
			unmapped.write(read)
		continue
	else:

		# Check if read is chimeric
		if read.has_tag('SA'):
			chimeric.write(read)

		# Check if read is soft clipped
		elif read.cigarstring.count('S') == 1:
			# Ensure read is not clipped due to adapter sequence, isn't secondary, and maps uniquely
			if not read.has_tag('XM') and not read.is_secondary and read.mapq != 0: #
				if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 5: #Soft clipping first and longer than 5 bp
					clipped_pos = read.cigartuples[0][1]
					clipped_read_quality = sum([ord(x) - 33 for x in read.qual[0:clipped_pos]])/clipped_pos
					if clipped_read_quality > 25: # phred score < 25
						clipped.write(read)
				elif read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 5:
					clipped_pos = read.cigartuples[-1][1]
					clipped_read_quality = sum([ord(x) - 33 for x in read.qual[-1*clipped_pos:]])/clipped_pos
					if clipped_read_quality > 25: # phred score < 25
						clipped.write(read)
	
	if not read.mate_is_unmapped:
     	# Chromosome fusion
		if read.reference_id!=read.next_reference_id:
			discordant.write(read)
		# Inversion
		elif read.is_reverse == read.mate_is_reverse:
			inversion.write(read)
		# Duplication
		elif ((read.pos < read.mpos and read.is_reverse) or (read.pos > read.mpos and read.mate_is_reverse)) and abs(read.isize) > 102:
			duplication.write(read)

align_file.close()
unmapped.close()
discordant.close()
inversion.close()
duplication.close()
clipped.close()
chimeric.close()

