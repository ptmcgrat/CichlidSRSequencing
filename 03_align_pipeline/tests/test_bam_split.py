"""route() must make the SAME decision as the old unit_scripts/split_bamfile_by_contig.py for every read it could ever handle. pysam is not needed: reads are mocks that offer
both the old pysam API (qual string, is_duplicate, pos, mpos, isize, ...) and the new one (flag, query_qualities, reference_start, ...)."""
import os, random, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from align_pipeline import bam_split
from align_pipeline.bam_split import route, is_chimeric

CIG = 'MIDNSHP=X'


class Read:
    def __init__(self, flag, quals, cigar, mapq, ref, nref, pos, mpos, isize, tags=()):
        self.flag, self.query_qualities, self.cigartuples, self.mapping_quality = flag, quals, cigar, mapq
        self.reference_id, self.next_reference_id, self.reference_start, self.next_reference_start, self.template_length = ref, nref, pos, mpos, isize
        self.tags = set(tags)

    def has_tag(self, t):
        return t in self.tags

    # ---- the OLD pysam-style API used by the original script
    qual = property(lambda s: None if s.query_qualities is None else ''.join(chr(q + 33) for q in s.query_qualities))
    is_duplicate = property(lambda s: bool(s.flag & 0x400))
    is_unmapped = property(lambda s: bool(s.flag & 0x4))
    is_secondary = property(lambda s: bool(s.flag & 0x100))
    is_reverse = property(lambda s: bool(s.flag & 0x10))
    mate_is_reverse = property(lambda s: bool(s.flag & 0x20))
    mate_is_unmapped = property(lambda s: bool(s.flag & 0x8))
    mapq = property(lambda s: s.mapping_quality)
    pos = property(lambda s: s.reference_start)
    mpos = property(lambda s: s.next_reference_start)
    isize = property(lambda s: s.template_length)
    cigarstring = property(lambda s: ''.join(f'{n}{CIG[op]}' for op, n in s.cigartuples) if s.cigartuples else '')


def old_destinations(read):
    """A transcription of the decision logic of the old script (write() calls replaced by names)."""
    total_read_quality = sum([ord(x) - 33 for x in read.qual]) / len(read.qual)
    if read.is_duplicate:
        return ()
    elif read.is_unmapped:
        return ('unmapped',) if total_read_quality > 25 else ()
    out = []
    if read.has_tag('SA'):
        out.append('chimeric')
    elif read.cigarstring.count('S') == 1:
        if not read.has_tag('XM') and not read.is_secondary and read.mapq != 0:
            if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 5:
                clipped_pos = read.cigartuples[0][1]
                if sum([ord(x) - 33 for x in read.qual[0:clipped_pos]]) / clipped_pos > 25:
                    out.append('clipped')
            elif read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 5:
                clipped_pos = read.cigartuples[-1][1]
                if sum([ord(x) - 33 for x in read.qual[-1 * clipped_pos:]]) / clipped_pos > 25:
                    out.append('clipped')
    if not read.mate_is_unmapped:
        if read.reference_id != read.next_reference_id:
            out.append('discordant')
        elif read.is_reverse == read.mate_is_reverse:
            out.append('inversion')
        elif ((read.pos < read.mpos and read.is_reverse) or (read.pos > read.mpos and read.mate_is_reverse)) and abs(read.isize) > 102:
            out.append('duplication')
    return tuple(out)


def random_read(rng):
    flag = 0x1
    for bit, p in ((0x4, 0.12), (0x8, 0.12), (0x10, 0.5), (0x20, 0.5), (0x100, 0.05), (0x400, 0.15), (0x800, 0.03)):
        if rng.random() < p:
            flag |= bit
    n = rng.choice([100, 125, 150, 151])
    hi = rng.random() < 0.5                                                   # high- or low-quality read, so means straddle 25
    quals = [rng.randint(15, 41) if hi else rng.randint(2, 38) for _ in range(n)]
    cigar = None
    if not flag & 0x4:
        k, m = rng.randint(1, 30), rng.randint(1, 30)
        cigar = rng.choice([[(0, n)], [(4, k), (0, n - k)], [(0, n - k), (4, k)], [(4, k), (0, n - k - m), (4, m)], [(0, 40), (1, 2), (0, n - 42)], [(4, 5), (0, n - 5)], [(4, 6), (0, n - 6)],
                            [(0, 50), (2, 3), (0, n - 50 - 5), (4, 5)], [(4, 10), (0, 60), (3, 100), (0, n - 70)]])
    tags = [t for t, p in (('SA', 0.1), ('XM', 0.05)) if rng.random() < p]
    ref = rng.choice([0, 0, 0, 1, 2])
    nref = ref if rng.random() < 0.8 else rng.choice([-1, 0, 1, 2, 3])
    pos = rng.randint(0, 10_000)
    mpos = pos + rng.choice([-5000, -300, -101, -50, 0, 50, 101, 103, 300, 5000])
    return Read(flag, quals, cigar, rng.choice([0, 0, 1, 20, 60, 60]), ref, nref, pos, mpos, rng.choice([-600, -103, -102, -101, 0, 101, 102, 103, 600]), tags)


class TestRoute(unittest.TestCase):
    def test_same_decisions_as_the_old_script_on_60000_random_reads(self):
        rng = random.Random(20261105)
        seen = {}
        for _ in range(60_000):
            r = random_read(rng)
            new, old = route(r), old_destinations(r)
            self.assertEqual(new, old, (r.flag, r.cigarstring, r.mapq, r.tags, r.reference_id, r.next_reference_id, r.pos, r.mpos, r.isize))
            seen[old] = seen.get(old, 0) + 1
        for needed in (('unmapped',), ('chimeric',), ('clipped',), ('discordant',), ('inversion',), ('duplication',), ('chimeric', 'discordant'), ('clipped', 'inversion'), ()):
            self.assertGreater(seen.get(needed, 0), 20, f'the random reads never produced {needed}: the comparison would prove little')

    def mk(self, flag=0x1 | 0x20, quals=None, cigar=((0, 150),), mapq=60, ref=0, nref=0, pos=100, mpos=300, isize=300, tags=()):
        return Read(flag, quals if quals is not None else [37] * 150, [tuple(c) for c in cigar] if cigar else None, mapq, ref, nref, pos, mpos, isize, tags)

    def test_boundaries(self):
        q = [37] * 5 + [20] * 145
        self.assertEqual(route(self.mk(cigar=((4, 5), (0, 145)), quals=q, flag=0x1 | 0x10 | 0x8)), (), 'a soft clip of exactly 5 bases is not enough')
        self.assertEqual(route(self.mk(cigar=((4, 6), (0, 144)), flag=0x1 | 0x10 | 0x8)), ('clipped',))
        self.assertEqual(route(self.mk(flag=0x1 | 0x4, quals=[25] * 150, cigar=None)), (), 'a mean quality of exactly 25 is not above 25')
        self.assertEqual(route(self.mk(flag=0x1 | 0x4, quals=[25] * 149 + [26], cigar=None)), ('unmapped',))

    def test_duplication_needs_an_outward_pair_with_an_insert_above_102(self):
        outward = dict(flag=0x1 | 0x10, pos=100, mpos=500)                    # reverse read upstream of its (forward) mate
        self.assertEqual(route(self.mk(isize=103, **outward)), ('duplication',))
        self.assertEqual(route(self.mk(isize=-103, **outward)), ('duplication',))
        self.assertEqual(route(self.mk(isize=102, **outward)), ())
        self.assertEqual(route(self.mk(isize=300, flag=0x1 | 0x20, pos=100, mpos=500)), (), 'a normal forward/reverse (inward) pair is neither duplication nor inversion evidence')

    def test_mate_on_another_chromosome_is_discordant_same_strand_is_inversion(self):
        self.assertEqual(route(self.mk(flag=0x1, ref=0, nref=2)), ('discordant',))
        self.assertEqual(route(self.mk(flag=0x1 | 0x10 | 0x20, ref=1, nref=1)), ('inversion',))
        self.assertEqual(route(self.mk(flag=0x1 | 0x8, ref=0, nref=2)), (), 'no mate evidence when the mate is unmapped')

    def test_duplicates_and_reads_without_position_are_never_routed(self):
        self.assertEqual(route(self.mk(flag=0x1 | 0x400, ref=0, nref=2)), ())
        self.assertEqual(route(self.mk(flag=0x1 | 0x400 | 0x4, cigar=None)), ())

    def test_read_without_qualities_does_not_crash(self):
        """The old script raised for these (and killed the whole contig job); the new one skips the quality-dependent rules."""
        self.assertEqual(route(Read(0x1 | 0x4, None, None, 0, 0, 0, 0, 0, 0)), ())
        self.assertEqual(route(Read(0x1 | 0x8, None, [(4, 20), (0, 130)], 60, 0, 0, 100, 100, 0)), ())


class FakePysam:
    """Just enough of pysam for split_bam: a reader that yields mock reads and writers that remember what they were given and leave a file behind."""
    reads = []
    written = {}

    class AlignmentFile:
        def __init__(self, path, mode='r', template=None, threads=1):
            self.path, self.mode = path, mode
            if mode == 'wb':
                FakePysam.written[path] = []

        def __iter__(self):
            return iter(FakePysam.reads)

        def write(self, read):
            FakePysam.written[self.path].append(read)

        def close(self):
            if self.mode == 'wb':
                open(self.path, 'w').write('bam')


class TestChimericOnly(unittest.TestCase):
    """Decision 41: only the chimeric BAM is written. The fast path must agree with route(), and the writer logic is checked with a fake pysam (pysam itself is not needed)."""
    def setUp(self):
        import tempfile
        self.dir = tempfile.mkdtemp() + '/'
        self.saved = sys.modules.get('pysam')
        FakePysam.written = {}
        fake = type(sys)('pysam'); fake.AlignmentFile = FakePysam.AlignmentFile
        sys.modules['pysam'] = fake

    def tearDown(self):
        if self.saved is None:
            sys.modules.pop('pysam', None)
        else:
            sys.modules['pysam'] = self.saved

    def test_is_chimeric_agrees_with_route_on_40000_random_reads(self):
        rng = random.Random(20261108)
        n = 0
        for _ in range(40_000):
            r = random_read(rng)
            self.assertEqual(is_chimeric(r), 'chimeric' in route(r), (r.flag, r.tags))
            n += is_chimeric(r)
        self.assertGreater(n, 500, 'the random reads hardly contained chimeric reads: the comparison would prove little')

    def mk(self, flag=0x1, tags=(), ref=0, quals=True):
        return Read(flag, [37] * 100 if quals else None, [(0, 100)], 60, ref, ref, 100, 300, 300, tags)

    def test_only_the_chimeric_file_is_written_and_counted(self):
        FakePysam.reads = [self.mk(tags=('SA',)), self.mk(tags=('SA',)), self.mk(flag=0x1 | 0x400, tags=('SA',)), self.mk(flag=0x1 | 0x4, tags=('SA',), ref=-1), self.mk(flag=0x1 | 0x10 | 0x20),
                           self.mk(flag=0x1 | 0x4, ref=0), self.mk(tags=('SA',), quals=False)]
        counts = bam_split.split_bam(self.dir + 'S1.all.bam', self.dir + 'S1')
        self.assertEqual(sorted(os.listdir(self.dir)), ['S1.chimeric.bam'], 'no other split file, and no .partial file left behind')
        self.assertEqual(len(FakePysam.written[self.dir + 'S1.chimeric.bam.partial']), 3)       # the two SA reads and the SA read without qualities; not the duplicate, not the unplaced read
        self.assertEqual((counts['reads_seen'], counts['chimeric'], counts['duplicates_skipped'], counts['no_position_skipped'], counts['no_quality']), (7, 3, 1, 1, 1))
        self.assertNotIn('discordant', counts)

    def test_the_other_types_can_still_be_requested_and_follow_route(self):
        INWARD = 0x1 | 0x20                                                  # a forward read whose mate is reverse: a normal pair, no rule applies to it
        FakePysam.reads = [self.mk(flag=INWARD, tags=('SA',)), self.mk(flag=0x1 | 0x10 | 0x20), self.mk(flag=INWARD)]
        counts = bam_split.split_bam(self.dir + 'S2.all.bam', self.dir + 'S2', outputs=('chimeric', 'inversion'))
        self.assertEqual(sorted(os.listdir(self.dir)), ['S2.chimeric.bam', 'S2.inversion.bam'])
        self.assertEqual((counts['chimeric'], counts['inversion']), (1, 1))

    def test_main_runs_by_file_path_and_writes_the_summary(self):
        import json
        FakePysam.reads = [self.mk(tags=('SA',))]
        self.assertEqual(bam_split.main([self.dir + 'S3.all.bam', self.dir + 'S3', '--threads', '1']), 0)
        summary = json.load(open(self.dir + 'S3.split_summary.json'))
        self.assertEqual((summary['reads_seen'], summary['chimeric']), (1, 1))
        self.assertEqual(bam_split.OUTPUTS, ('chimeric',))


if __name__ == '__main__':
    unittest.main()
