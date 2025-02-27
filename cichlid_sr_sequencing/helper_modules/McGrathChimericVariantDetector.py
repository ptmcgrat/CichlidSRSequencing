import pysam, time, operator, sys, scipy, collections, shutil, os, random, pdb
import numpy as np
from functools import reduce
from statistics import median
from collections import defaultdict
from skbio.alignment import StripedSmithWaterman
from subprocess import call, Popen, PIPE
from operator import attrgetter, itemgetter
from Bio import pairwise2
        

strand_converter = {True: -1,
                    False: 1,
                    '+': 1,
                    '-': -1}

COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}                
    
def reverse_complement(sequence):
    return ''.join(COMPLEMENT[b] for b in sequence[::-1])
   
def genotype_bam(bam_obj, chrom, pos):
    pu = bam_obj.pileup(chrom, pos, pos+1)
    for puc in pu:
        if puc.pos == pos:
            return genotype_pucpileup(puc.pileups)
    return []

def genotype_pucpileup(puc_pileup):
    base_g = []
    for pur in puc_pileup:
        if pur.is_del or pur.is_refskip:
            continue
        elif pur.indel < 0:
            qp = pur.query_position
            base_g.append(pur.alignment.seq[qp] + str(pur.indel) )
        else:
            qp = pur.query_position
            if pur.indel > 0:
                base_g.append(pur.alignment.seq[qp:qp+1+pur.indel])
            else:
                base_g.append(pur.alignment.query_sequence[qp])
    return base_g

def overlap_bam(bam_obj, chrom, pos):
    reads = bam_obj.fetch(chrom, pos, pos+1)
    out_dist = []
    clipped = 0
    for read in reads:
        if not read.is_unmapped and not read.is_secondary and read.mapq > 0:
            if 'S' in read.cigarstring:
                clipped += 1
            start = read.get_blocks()[0][0]
            stop = read.get_blocks()[-1][1]-1
            out_dist.append(min([pos - start, stop - pos]))
    return out_dist, clipped

class ChimericRead():
    def __init__(self, read):
        # read is a pysam read object
        # Assumes pos attribute of read is 0 based but SA tag position is 1 based

        # Make sure the read is chimeric
        try:
            SA_tag = read.get_tag('SA').split(',')
        except KeyError:
            raise TypeError('ChimericRead Error: read must contain the SA tag')

        # Make sure the primary alignment is simple
        if len(read.cigar) != 2:
            raise TypeError('ChimericRead Warning: read cigar too complicated ' + str(read.cigar))

        # Make sure that the secondary alignment is simple
        chi_cigarstring = SA_tag[3]
        if chi_cigarstring.count('S') != 1 or chi_cigarstring.count('M') != 1 or chi_cigarstring.count('I') != 0:
            raise TypeError('ChimericRead Warning: read chimeric cigar too complicated ' + str(chi_cigarstring) + ' ' + read.qual)

        first = {} # Used to keep track of info for the left alignment
        secondary = {} # Used to keep track of the info for the right alignment

        # Matched part of read first e.g. 35M65S
        if read.cigar[0][0] == 0:
            first['MB'] = read.cigar[0][1] # Mapped bases
            first['Chr'] = read.rname
            first['Pos'] = read.pos + first['MB'] - 1 # Already 0-coordinate; Position of the junction point
            first['Dir'] = 'right' # Direction of chimeric junction
            first['Strand'] = '-' if read.is_reverse else '+' 
            secondary['Chr'] = int(SA_tag[0])
            secondary['Strand'] = SA_tag[2]
            if first['Strand'] == secondary['Strand']:
                # Matched part of read second e.g. Secondary hit: 34S67M
                try:
                    secondary['MB'] = int(chi_cigarstring.split('S')[1].split('M')[0])
                except ValueError:
                    pdb.set_trace()
                secondary['Pos'] = int(SA_tag[1]) - 1 # Convert 1-coordinate to 0-coordinate
                secondary['Dir'] = 'left'
            else:
                # Matched part of read first: e.g. Secondary hit: 67M34S
                secondary['MB'] = int(chi_cigarstring.split('M')[0])
                secondary['Pos'] = int(SA_tag[1]) + secondary['MB'] - 2
                secondary['Dir'] = 'right'
                
            overlap = len(read.seq) - secondary['MB'] - first['MB']
            if overlap > 0:
                ins_bases = read.seq[first['MB']: len(read.seq) - secondary['MB']]

        if read.cigar[0][0] == 4:
            # Matched part of read second e.g. 34S67M first hit is Soft clipping
            first['MB'] = read.cigar[1][1]
            first['Chr'] = read.rname # Chromosome/contig id
            first['Pos'] = read.pos # Already 0-coordinate system
            first['Dir'] = 'left'
            first['Strand'] = '-' if read.is_reverse else '+' # True for positive strand, False for negative strand
            secondary['Chr'] = str(SA_tag[0])
            secondary['Strand'] = SA_tag[2]
            if first['Strand'] == secondary['Strand']:
                # Matched part of read first: e.g. Secondary hit: 67M34S
                secondary['MB'] = int(chi_cigarstring.split('M')[0])
                secondary['Pos'] = int(SA_tag[1]) + secondary['MB'] - 2
                secondary['Dir'] = 'right'
            else:
                # Matched part of read second e.g. Secondary hit: 34S67M
                secondary['MB'] =  int(chi_cigarstring.split('S')[1].split('M')[0])
                secondary['Pos'] = int(SA_tag[1]) - 1 # Convert 1-coordinate to 0-coordinate
                secondary['Dir'] = 'left'
 
            overlap = len(read.seq) - secondary['MB'] - first['MB']
            if overlap > 0: # There are extra bases in the read that we need to keep track of
                ins_bases = read.seq[secondary['MB']: len(read.seq) - first['MB']]
            
        # Determine order of mappings
        flip_order = False
        if first['Chr'] == secondary['Chr']: #Both alignments are on the same chromosome, use location info
            if first['Pos'] > secondary['Pos']:
                flip_order = True
        elif first['Chr'] > secondary['Chr']:
            flip_order = True

        # Must deal with parts of read mapped with both alignments. Subtract from the less than alignment
        if overlap <=0: 
            ins_bases = ''
            if not flip_order:
                if first['Dir'] == 'left':
                    first['Pos'] -= overlap
                elif first['Dir'] == 'right':
                    first['Pos'] += overlap
                else:
                    raise(NameError('first_direction unknown value'))
            else:
                if secondary['Dir'] == 'left':
                    secondary['Pos'] -= overlap
                elif secondary['Dir'] == 'right':
                    secondary['Pos'] += overlap
                else:
                    raise(NameError('first_direction unknown value'))

        # Make sure order is consistent
        if not flip_order:
            self.data = (first['Chr'], first['Pos'], first['Dir'][0], secondary['Chr'], secondary['Pos'], secondary['Dir'][0], ins_bases)
        else:
            self.data = (secondary['Chr'], secondary['Pos'], secondary['Dir'][0], first['Chr'], first['Pos'], first['Dir'][0], ins_bases)

        if self.data[0] == self.data[3]:
            if self.data[2] == 'r' and self.data[5] == 'l':
                self.d_type = 'del'
            elif self.data[2] == 'l' and self.data[5] == 'r':
                self.d_type = 'dup'
            else:
                self.d_type = 'inv'
        else:
            self.d_type = 'cf'

        self.data = self.data + (self.d_type,)
        
class ChimericCaller():
    def __init__(self, reffile, outvcffile):
                    
        self.refObj = pysam.FastaFile(reffile)
        self.contigNames = self.refObj.references
        self.outvcffile = outvcffile

        #self.outLocationsFile = outLocationsFile

        self.t_polys = {} #Stores all polys

    def identifyChimericVariants(self, discoveryBams, minChimericReads = 5, maxLength = 100000, minMapQ = 10):

        discoveryBamObjs = []

        for bamfile in discoveryBams:
            discoveryBamObjs.append(pysam.AlignmentFile(bamfile))

        f = open(self.outvcffile.replace('.vcf', '.tsv'), 'w')

        for contig in self.contigNames:
            discoveryChimeras = defaultdict(int)
            discoveryChimerasL = defaultdict(list)

            for i,bamObj in enumerate(discoveryBamObjs):
                for read in bamObj.fetch(contig):
                    if not read.is_secondary and read.mapq > minMapQ:
                        try:
                            SA = read.get_tag('SA')
                        except KeyError:
                            continue
                        read.set_tag('SA',SA.replace(SA.split(',')[0], str(self.contigNames.index(SA.split(',')[0]))))
                        try:
                            newRead = ChimericRead(read)
                        except TypeError:
                            continue
                        discoveryChimeras[newRead.data] += 1
                        discoveryChimerasL[newRead.data].append(i)

            print(contig + ': ' + str(len(discoveryChimeras)) + ' potential chimeric sites', file = sys.stderr)
            samples = [x.split('/')[-1].split('.')[0] for x in discoveryBams]
            print('\t'.join(['Location'] + samples), file = f)
            for loc, counts in discoveryChimeras.items():
                if counts < minChimericReads:
                    continue
                
                start = loc[1]
                stop = loc[4]
                d_type = loc[7]

                if abs(stop - start) < maxLength and d_type != 'cf':
                    out_counts = []
                    for i,bam_file in enumerate(discoveryBams):
                        out_counts.append(discoveryChimerasL[loc].count(i))
                    print('\t'.join(str(x) for x in [loc] + out_counts ), file = f)
                    continue
                    
                if d_type == 'del':
                    refbase = self.refObj.fetch(contig,start,stop).upper()
                    altbase = self.refObj.fetch(contig,start,start+1).upper() + loc[6]
                    tpoly = poly(contig, start, refbase, altbase, self.refObj)

                elif d_type == 'dup':
                    refbase = self.refObj.fetch(contig,stop,stop+1).upper()
                    altbase = refbase + self.refObj.fetch(contig,start,stop+1).upper() + loc[6]
                    tpoly = poly(contig, stop, refbase, altbase, self.refObj, is_dup = True)

                elif d_type == 'inv':
                    if start > stop:
                        continue
                    refbase = ''
                    altbase = ''
                    tpoly = poly(contig, start, stop, altbase, self.refObj, is_inv=True)

                elif d_type == 'cf':
                    continue
                    
                else:
                    print('Unknown d_type: ' + str(d_type), file = sys.stderr)
                    sys.exit()

                self.t_polys[(contig, start, refbase, altbase)] = tpoly
                
            print('\t' + str(len(self.t_polys)) + ' total candidate chimeric events')

        print('Creating Poly object')
        tPolys = Polys(self.t_polys, self.refObj)
        print('Genotyping Poly object')
        f.close()
        #tPolys.genotypePolys(genotypeBamObjs)
        #tPolys.create_VCFfile(self.outvcffile)

    def identifyLargeInsertions(self, discordant_bamfiles, binsize = 500, overlap = 250, min_reads=1):
        # Identify large insertions
        
        for i,bam_file in enumerate(discordant_bamfiles):
            discordant = pysam.AlignmentFile(bam_file)
        
            for i in range(0,discordant.nreferences):
                contig = discordant.references[i]
                length = discordant.lengths[i]
                current_hit = False
                hit_bin = False
                for j in range(0, length-bin_size, overlap):
                    reads = self.coverage_calculator(discordant, contig, j, j+bin_size)
                    #Does bin have necessary number of forward and reverse discordant reads?
                    if reads[0] > min_reads and reads[0] > 5 * reads[1]:
                        if j + 4*bin_size < length:
                            for k in range(1,4):
                                reads2 = self.coverage_calculator(discordant, contig, j + k*bin_size, j + (k+1)*bin_size)
                                if reads2[1] > min_reads and reads2[1] > 5*reads2[0]:
                                    hit_bin = True
                                    
                                    # Have we reached the end of the hit?
                    if current_hit == True and hit_bin == False:
                        mid_pos = int(j + self.bin_size/2)
                        #Try to identify sequence of 
                        g = DBGraph(50, contig + ':' + str(mid_pos) + '_' + str(mid_pos), self.refObj)
                        all_reads = self.collect_reads(contig + ':' + str(mid_pos) + '_' + str(mid_pos), discordant)
                        for read in all_reads:
                            g.add_read(read)
                        a = g.identify_polymorphism()
                        print(a.chrom + '\t' + str(a.pos) + '\t' + a.ref + '\t' + a.alt)

                        #if a is not None:
                        #    mut_geno = poly.genotype_bam(mut_all, sample)
                        #    if mut_geno[0] == '1/1':
                        #        self.tempPolys[(a.chrom, a.start, a.stop, a.alt)] = a

                        current_hit = False
                        
                    if hit_bin == True:
                        current_hit = True
                    hit_bin = False
                
    def _coverage(self,bam_obj):
        stats = pysam.idxstats(bam_obj.filename).rstrip().split('\n')
        tot_reads = sum([int(x.split('\t')[2]) for x in stats])

        tot_bases = sum(bam_obj.lengths)
        
        return tot_reads/tot_bases*self._avg_read_len(bam_obj)

    def _avg_coverage(self,bam_obj, chrom, pos1, pos2):
        n_reads = 0
        if pos2 - pos1 > 10000:
            return bam_obj.count(chrom, pos1, pos2)*avg_read_len(bam_obj)/(pos2-pos1+1)
        else:
            for puc in bam_obj.pileup(chrom, pos1, pos2):
                if puc.pos >= pos1 and puc.pos <= pos2:
                    for pur in puc.pileups:
                        if pur.is_del or pur.is_refskip:
                            pass
                        else:
                            n_reads += 1
        return n_reads/(pos2 - pos1 + 1)

    def avg_read_len(self, bam_obj):
        count = 0
        read_len = []
        for read in bam_obj.fetch():
            read_len.append(len(read.seq))
            count += 1
            if count > 100000:
                break
        return median(read_len)
    
    def ret_long_clip(self, position):
        lengths = []
        consensus = []
        for cl_read in self.single_reads[position]:
            lengths.append(len(cl_read.clipped_seq))
        if len(lengths) == 0:
            return ''
        max_len = max(lengths)
        for cl_read in self.single_reads[position]:
            if len(cl_read.clipped_seq) == max_len:
                return cl_read
                   
                        
    def ret_consensus(self, position):
                   
        lengths = []
        consensus = []
        if position[0] is not None:
            index = 1
        if position[1] is not None:
            index = 0
        for cl_read in self.single_reads[position]:
            lengths.append(len(cl_read.clipped_seq[index]))
        if len(lengths) == 0:
            return ''
        for i in range(0,max(lengths)):
            temp_dict = defaultdict(int)
            total = 0
            for cl_read, length in zip(self.single_reads[position], lengths):
                if i < length:
                    temp_dict[cl_read.clipped_seq[index][i*(2*index-1)]] += 1
                    total += 1
                    cons_char = max(iter(temp_dict.items()), key=operator.itemgetter(1))[0]
            if temp_dict[cons_char]/total >= .74 and total > 1:
                consensus.append(cons_char)
            elif total == 1:
                consensus.append('')
            else: 
                consensus.append('N')
        consensus = ''.join(consensus)
        return consensus

    def ret_nearby_cov(self, bam_obj, position, delta):
        if position[0] == None:
            chrom = bam_obj.references[position[1][0]]
            stop = position[1][1]
            return avg_coverage(bam_obj, chrom, stop-delta, stop)
        else:
            chrom = bam_obj.references[position[0][0]]
            start = position[0][1]
            return avg_coverage(bam_obj, chrom, start, start+delta)
                   
    def ret_position(self, position):
        if position[0] is not None:
            return position[0][0]*25000000 + position[0][1]
        else:
            return position[1][0]*25000000 + position[1][1]


def collect_reads(transposon, q_bam_file, or_bam_file):
    all_reads = []
    
    reads = q_bam_file.fetch(transposon[0], transposon[1] - 250, transposon[1])
    for read in reads:
        if not read.is_reverse:
            all_reads.append(read.seq)
            read_name = read.query_name
            mate_contig = read.next_reference_name
            mate_pos = read.mpos
            o_reads = or_bam_file.fetch(mate_contig, max(0,mate_pos - 1), mate_pos + 1)
            count = 0
            for o_read in o_reads:
                if o_read.query_name == read_name:
                    count += 1
                    if o_read.is_reverse:
                        all_reads.append(o_read.seq)
                    else:
                        all_reads.append(reverse_complement(o_read.seq))
                
    reads = q_bam_file.fetch(transposon[0], transposon[1], transposon[1]+250)
    for read in reads:
        if read.is_reverse:
            all_reads.append(read.seq)
            read_name = read.query_name
            mate_contig = read.next_reference_name
            mate_pos = read.mpos
            o_reads = or_bam_file.fetch(mate_contig, max(0,mate_pos - 1), mate_pos + 1)
            count = 0
            for o_read in o_reads:
                if o_read.query_name == read_name:
                    count += 1
                    if not o_read.is_reverse:
                        all_reads.append(o_read.seq)
                    else:
                        all_reads.append(reverse_complement(o_read.seq))
def coverage_calculator(bam_file, contig, start, stop):
    f_reads = 0
    r_reads = 0
    reads = bam_file.fetch(contig, start, stop)
    for read in reads:
        if 'S' not in read.cigarstring and 'H' not in read.cigarstring:
            if read.is_reverse:
                r_reads += 1
            else:
                f_reads += 1
    return (f_reads, r_reads)

        
class Polys():

    def __init__(self, poly_list, ref_obj):
        if type(poly_list) is dict:
            self.polys = list(poly_list.values())
        else:
            self.polys = poly_list
        self.ref_obj = ref_obj
        self.geno_strains = set()
        #self.out_vcf = SO.get_vcf_file(strains[0], self.genome_version)
        
    def genotypePolys(self, genotypeBams):
        for bamobj in genotypeBams:
            strain = bamobj.header['RG'][0]['SM']
            print('Genotyping ' + strain + '...')
            self.geno_strains.add(strain)
            for poly in self.polys:
                poly.genotype_bam(bamobj, strain, add = True)
                    
    def vcfLoader(self, infile, add_geno = False):
        with open(infile) as f:
            sample_order = []
            for line in f:
                if line[:2] == '##':
                    continue
                else:
                    tokens = line.rstrip().split('\t')
                    if len(tokens) < 2:
                        continue
                    if tokens[0] == '#CHROM':
                        if add_geno:
                            for strain in tokens[9:]:
                                self.geno_strains.add(strain.replace('Sample',''))
                                sample_order.append(strain.replace('Sample', ''))
                            continue
                        else:
                            continue
                    chrom = tokens[0]
                    pos = int(tokens[1]) - 1
                    ref = tokens[3]
                    alt = tokens[4]
                    if alt == '<DEL>':
                        stop = int(line.split('END=')[1].split(';')[0])
                        ref = self.ref_obj.fetch(chrom, pos, stop)
                        alt = self.ref_obj.fetch(chrom, pos, pos+1)

                    snp_type = tokens[7].split('SVTYPE=')[1].split(';')
                    if snp_type == 'DUP':
                        tpoly = poly(chrom, pos, ref, alt, self.ref_obj, is_dup = True)
                    elif snp_type == 'INV':
                        tpoly = poly(chrom, pos, ref, alt, self.ref_obj, is_inv = True)
                    else:
                        tpoly = poly(chrom, pos, ref, alt, self.ref_obj)

                    if add_geno:
                        for i,tgeno in enumerate(tokens[9:]):
                            genos = tgeno.split(':')
                            tpoly.geno[sample_order[i]] = (genos[0],(int(genos[1]),int(genos[2]),int(genos[3])))
                            
                    self.polys.append(tpoly)


    def create_VCFfile(self, filename, geno = True):
        self.polys = sorted(self.polys, key = attrgetter('chrom', 'pos'))
        strain_order = sorted(self.geno_strains)
        with open(filename, 'w') as f:
            print(self.create_VCFheader(geno), file = f)
            prev_poly = None
            count = 0
            for poly in self.polys:
                if poly.empty:
                    continue
                if not self.polys_overlap(poly, prev_poly):
                    print(poly.ret_VCFrecord(strain_order, geno), file = f)
                    count += 0
                    prev_poly = poly
                else:
                    print('Overlapping polys (2nd poly not added):', file = sys.stderr)
                    print(prev_poly.ret_VCFrecord(strain_order, geno), file = sys.stderr)
                    print(poly.ret_VCFrecord(strain_order, geno), file = sys.stderr)
        print(str(count) + ' total polys printed to a vcf file', file = sys.stderr)

    def create_VCFheader(self, geno):
        out_string = '##fileformat=VCFv4.2\n'
        out_string += '##fileDate=' + time.strftime('%Y%m%d') + '\n'
        out_string += '##source=BamTools.py\n'
        out_string += '##reference=' + self.ref_obj.filename.decode('utf-8') + '\n'
        for i in range(0, self.ref_obj.nreferences):
            name = self.ref_obj.references[i]
            length = self.ref_obj.lengths[i]
            out_string += '##contig=<ID=' + name + ',length=' +str(length) + '>\n'
        out_string += '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n'
        out_string += '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of variant">\n'
        out_string += '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of Variant">\n'
        out_string += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        out_string += '##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Reference Reads">\n'
        out_string += '##FORMAT=<ID=AR,Number=1,Type=Float,Description="Alt Reads">\n'
        out_string += '##FORMAT=<ID=BR,Number=1,Type=Float,Description="Ambiguous Reads">\n'
        if geno:
            out_string += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(sorted(self.geno_strains))
        else:
            out_string += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'

        return out_string

    def annotate_SNPs(self, outfile, snpeff_database):
        print ('Annotating')
        call(['snpeff', 'eff', snpeff_database, self.out_vcf], stdout = open(outfile, 'w'))
        #call(['python3', '/Users/pmcgrath7/Dropbox/Projects/Sequencing/Scripts/modify_snpEff.py', tempfile, self.genome_version], stdout = open(SO.get_vcf_file(self.mut_strain, self.genome_version, annotation = True), 'w'))
        #call(['rm','-f',tempfile, self.out_vcf])

    def polys_overlap(self, poly1, poly2):
        if poly2 is None:
            return False
        if poly1.chrom != poly2.chrom:
            return False
        if poly1.is_inv:
            if not poly2.is_inv:
                return False
        if abs(poly1.start - poly2.start) < 6:
            return True
        if abs(poly1.stop - poly2.stop) < 6:
            return True
        return False
            
    def plot_cluster(self):
        num_polys = len(self.polys)
        num_samples = len(self.geno_strains)
        geno_data = numpy.empty(shape = (num_polys, num_samples))
        for i, poly in enumerate(self.polys):
            geno_data[i] = poly.ret_geno_vector()
        
        
class poly():
    def __init__(self, chrom, pos, ref, alt, ref_o, is_inv = False, is_dup = False, is_ins = False):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.geno = {}
        self.ref_o = ref_o
        self.start = self.pos
        self.empty = True
        if is_inv:
            self.stop = ref
            self.alt = 'INV'
        else:
            self.stop = self.pos + len(self.ref)
        self.is_inv = is_inv
        self.is_dup = is_dup
        self.is_ins = is_ins
        self.max_len = 20000

    def __cmp__(self, other):
        if self.chrom != other.chrom:
            return self.chrom.__cmp__(other.chrom)
        else:
            return self.pos.__cmp__(other.pos)
        
    def ret_ref_allele(self, flanking):
        if self.is_inv:
            seq = self.ref_o.fetch(self.chrom, max(0,self.start-flanking), self.stop + flanking).upper()
            if len(seq) > 3*flanking:
                return seq[:2*flanking] + seq[-2*flanking:]
            else:
                return seq
        else:
            return self.ref_o.fetch(self.chrom, max(0,self.start - flanking), self.stop + flanking).upper()

    def ret_alt_allele(self, flanking):
        if self.is_inv:
            seq = self.ref_o.fetch(self.chrom, max(0,self.start - flanking), self.start) + reverse_complement(self.ref_o.fetch(self.chrom, self.start, self.stop+1).upper()) + self.ref_o.fetch(self.chrom, self.stop+2, self.stop + flanking)
            if len(seq) > 3*flanking:
                return seq[:2*flanking] + seq[-2*flanking:]
            else:
                return seq
        else:
            return self.ref_o.fetch(self.chrom, max(0,self.start - flanking), self.start).upper() + self.alt  + self.ref_o.fetch(self.chrom, self.stop, self.stop + flanking).upper() 
        
    def genotype_bam(self, bam_file, strain, read_length = 100, add = False, troubleshooting = False):
        unique_reads = {}
        ref_seq = self.ret_ref_allele(read_length)
        alt_seq = self.ret_alt_allele(read_length)
        if troubleshooting != False:
            print(ref_seq, file = troubleshooting)
            print(alt_seq, file = troubleshooting)
        reads1 = bam_file.fetch(self.chrom, max(0,self.start-10), self.start + 10)
        for read in reads1:
            if read.seq is None:
                continue
            if read.pos - 5 < self.start and read.pos + len(read.seq) > self.start:
                unique_reads[read.qname] = read
        reads2 = bam_file.fetch(self.chrom, max(0,self.stop-10),  self.stop + 10)
        for read in reads2:
            if read.seq is None:
                continue
            if read.pos -5 < self.stop and read.pos + len(read.seq) > self.stop:
                unique_reads[read.qname] = read
        geno = [0,0,0]
        ref_SSW = StripedSmithWaterman(ref_seq)
        alt_SSW = StripedSmithWaterman(alt_seq)
       
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
                if troubleshooting != False:
                    print(read.seq + '\t' + str(ref_SSW(read.seq).optimal_alignment_score) + '\t' + str(alt_SSW(read.seq).optimal_alignment_score), file = troubleshooting)
                    print(reverse_complement(read.seq) + '\t' + str(ref_SSW(reverse_complement(read.seq)).optimal_alignment_score) + '\t' + str(alt_SSW(reverse_complement(read.seq)).optimal_alignment_score), file = troubleshooting)
                    
        if len(unique_reads) == 0:
            geno_o = ('./.',tuple(geno))
        if geno[0] == 0 and geno[1] == 0:
            geno_o = ('./.',tuple(geno))
        elif geno[0] > 4*geno[1]:
            geno_o = ('0/0',tuple(geno))
        elif geno[0] < 1/4*geno[1]:
            geno_o = ('1/1',tuple(geno))
            self.empty = False
        else:
            geno_o = ('0/1',tuple(geno))
            self.empty = False
        if add == True:
            self.geno[strain] = geno_o
        return geno_o
    
    def ret_VCFrecord(self, strain_order, geno):

        out_string = self.chrom + '\t' + str(self.pos+1) + '\t.\t'
        if self.is_inv:
            if self.stop - self.start < self.max_len and len(self.alt) < self.max_len:
                out_string += str(self.ref_o.fetch(self.chrom, self.start-1, self.stop).upper()) + '\t' + str(reverse_complement(self.ref_o.fetch(self.chrom, self.start-1, self.stop).upper())) 
            else:
                out_string += self.ref_o.fetch(self.chrom, self.start-1, self.start).upper() + '\t' + '<INV>'
        else:
            if self.stop - self.start < self.max_len:
                out_string += str(self.ref) + '\t' + self.alt
            elif self.is_dup:
                out_string += self.ref_o.fetch(self.chrom, self.start-1, self.start).upper() + '\t' + '<DUP>'
            else:
                out_string += self.ref_o.fetch(self.chrom, self.start-1, self.start).upper() + '\t' + '<DEL>'
                
        out_string += '\t30\tPASS\tNS='  + str(len(self.geno)) + ';END=' + str(self.stop) + ';SVTYPE='
        if self.is_inv:
            out_string += 'INV'
        elif self.is_dup:
            out_string += 'DUP'
        elif self.is_ins:
            out_string += 'INS'
        else:
            out_string += 'DEL'
            
        out_string += ';SVLEN=' + str(self.stop - self.start + 1) + '\tGT:RR:AR:BR'

        if geno:
            for strain in strain_order:
                out_string += '\t' + self.geno[strain][0] + ':' + str(self.geno[strain][1][0]) + ':' + str(self.geno[strain][1][1]) + ':' + str(self.geno[strain][1][2])
        return out_string

    def ret_geno_vector(self):
        out_genos = [0]*len(self.geno)
        for i, geno in enumerate(self.geno):
            if geno[0] == '0/0':
                out_genos[i] = 0
            if geno[0] == '0/1':
                out_genos[i] = 1
            if geno[0] == '1/1':
                out_genos[i] = 2
        return out_genos
                
    
class DBGraph:
    def __init__(self, k_size, name, ref_o):
        self.k = k_size
        self.ref_o = ref_o
        self.nodes = {}
        self.name = name
                
    def add_read(self, read):
        for i in range(0,len(read) - self.k):
            seq1 = read[i:self.k+i]
            seq2 = read[i+1:self.k+1+i]
            if seq1 not in self.nodes:
                self.nodes[seq1] = db_node(seq1)
            else:
                self.nodes[seq1].num += 1
            if seq2 not in self.nodes:
                self.nodes[seq2] = db_node(seq2)
            self.nodes[seq1].f_edges[seq2] += 1
            self.nodes[seq2].r_edges[seq1] += 1

    def clean_graph(self):
        del_keys = []
        for node in self.nodes.values():
            if sum(node.f_edges.values()) <= 2 and sum(node.r_edges.values()) <= 2:
                for key, value in node.f_edges.items():
                    self.nodes[key].r_edges.pop(node.sequence)
                for key, value in node.r_edges.items():
                    self.nodes[key].f_edges.pop(node.sequence)
                del_keys.append(node.sequence)
        for key in del_keys:
            self.nodes.pop(key)
            
    def follow_node(self, seq):
        out_seq = seq
        out_count = 0
        count = 0
        while(True and count < 2000):
            try:
                new_seq = max(self.nodes[seq].f_edges.items(), key=operator.itemgetter(1))[0]
                out_seq += new_seq[-1]
                seq = new_seq
                count += 1
            except ValueError:
                break
        if count != 2000:
            return out_seq
        else:
            return ''
            
    def identify_polymorphism(self):
        self.clean_graph()
        long_seqs = []
        for seq, node in self.nodes.items():
            if sum(node.r_edges.values()) == 0:
                created_seq = self.follow_node(node.sequence)
                if len(created_seq) > 300:
                    long_seqs.append(created_seq)
        good_seqs = [True]*len(long_seqs)
        for i in range(0,len(long_seqs)):
            for j in range(0,len(long_seqs)):
                if i != j:
                    if long_seqs[i][-50:] == long_seqs[j][-50:]:
                        if len(long_seqs[i]) > len(long_seqs[j]):
                            good_seqs[j] = False
                        else:
                            good_seqs[i] = False
        out_seqs = []
        for i in range(0,len(long_seqs)):
            if good_seqs[i]:
                out_seqs.append(long_seqs[i])
        #print(out_seqs)
            
        if len(out_seqs) in [1,2]:
            contig = self.name.split(':')[0]
            pos1 = int(self.name.split(':')[1].split('_')[0])
            pos2 = int(self.name.split(':')[1].split('_')[1])            

            ref_seq = str(self.ref_o.fetch(contig, pos1-250, pos2+250)).replace('N','')
            self.out_seqs = out_seqs
            self.ref_seq = ref_seq
            if len(out_seqs) == 1:                
                c = out_seqs[0]
                out_aln = pairwise2.align.globalms(c, ref_seq, 2, -1, -2, -.1)
            else:
                a = out_seqs[0]
                b = out_seqs[1]
                a1 = pairwise2.align.globalms(a+'NNNNN' + b, ref_seq, 2, -1, -2, -.1)
                a2 = pairwise2.align.globalms(b+'NNNNN' + a, ref_seq, 2, -1, -2, -.1)
                if a1[0][2] > a2[0][2]:
                    out_aln = a1
                else:
                    out_aln = a2

            out_str = ''
            out_strs = []
            for c in out_aln[0][1]:
                if c is not '-':
                    out_str += c
                else:
                    if len(out_str) > 50:
                        out_strs.append(out_str)
                    out_str = ''
            if len(out_str) > 50:
                out_strs.append(out_str)
                
            if len(out_strs) == 2:
                i1 = out_aln[0][1].find(out_strs[0]) + len(out_strs[0])
                i2 = out_aln[0][1].find(out_strs[1])
                count = 0
                for i,c in enumerate(out_aln[0][1]):
                    if c is not '-':
                        count += 1
                    if i == i1:
                        break
                inserted_bases = out_aln[0][0][i1-1:i2].replace('-','')
                deleted_bases = out_aln[0][1][i1:i2].replace('-','')
                position = pos1-250+count
                reference_bases = str(self.ref_o.fetch(contig, position-1, position + len(out_aln[0][1][i1:i2].replace('-',''))))
                vcf_line = '\t'.join([contig, str(position), '.', reference_bases, inserted_bases, '30', 'PASS', 'NS=0', 'GT:GQ:DP'])
                print('InsertionCandidate: ' + inserted_bases + '\t' + str(out_aln[0][2]))
                if len(inserted_bases) > 50 and out_aln[0][2] > 850:
                    return self.poly(contig, position, reference_bases, inserted_bases, self.ref_o, is_ins = True)

class db_node:
    def __init__(self, sequence):
        self.sequence = sequence
        self.num = 1
        self.f_edges = defaultdict(int)
        self.r_edges = defaultdict(int)
                
