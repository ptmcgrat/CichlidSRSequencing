import pysam

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

        if read.rname != int(SA_tag[0]):
            raise TypeError('ChimericRead Error: both alignments must be on same chromosome')

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
            secondary['Chr'] = int(SA_tag[0])
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