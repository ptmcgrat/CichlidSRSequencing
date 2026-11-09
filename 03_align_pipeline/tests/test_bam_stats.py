"""The statistics the database gets. Numbers are chosen so every derived value can be checked by hand."""
import math, os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from align_pipeline import bam_stats as bs

STATS = '\n'.join([
    '# This file was produced by samtools stats (1.17+htslib-1.17)', '# The command line was:  stats -@ 4 x.all.bam', '# CHK, Checksum\t[2]Read Names\t[3]Sequences\t[4]Qualities',
    'SN\traw total sequences:\t100000000\t# excluding supplementary and secondary reads by default', 'SN\tfiltered sequences:\t0', 'SN\tsequences:\t100000000',
    'SN\treads mapped:\t97000000', 'SN\treads mapped and paired:\t95000000\t# paired-end technology bit set + both mates mapped', 'SN\treads unmapped:\t3000000',
    'SN\treads properly paired:\t93000000\t# proper-pair bit set', 'SN\treads duplicated:\t8730000\t# PCR or optical duplicate bit set', 'SN\treads MQ0:\t4000000\t# mapped and MQ=0',
    'SN\tsupplementary alignments:\t500000', 'SN\tbases mapped (cigar):\t14300000000\t# more accurate', 'SN\tbases duplicated:\t1309500000',
    'SN\terror rate:\t3.496503e-03\t# mismatches / bases mapped (cigar)', 'SN\taverage length:\t150', 'SN\taverage quality:\t36.5',
    'SN\tinsert size average:\t380.5', 'SN\tinsert size standard deviation:\t72.1', 'FFQ\t1\t2\t3']) + '\n'

IDXSTATS = 'NC_135167.1\t45606000\t5000000\t10000\nNC_135168.1\t48807247\t5500000\t12000\nNW_027490041.1\t900000\t1000\t5\n*\t0\t0\t1200000\n'

# the layout of a GATK/Picard MarkDuplicates metrics file: header comments, the DuplicationMetrics table, a blank line, then the histogram
PICARD = '\n'.join([
    '## htsjdk.samtools.metrics.StringHeader', '# MarkDuplicates INPUT=[x.sorted.bam] OUTPUT=x.all.bam METRICS_FILE=x.metrics.txt', '## htsjdk.samtools.metrics.StringHeader', '# Started on: Thu Nov 05 12:00:00 EST 2026', '',
    '## METRICS CLASS\tpicard.sam.DuplicationMetrics',
    'LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE',
    'libA\t1000000\t48000000\t500000\t1200000\t80000\t4325000\t300000\t0.09\t250000000', '',
    '## HISTOGRAM\tjava.lang.Double', 'BIN\tVALUE', '1.0\t1.0']) + '\n'


class TestParsing(unittest.TestCase):
    def test_samtools_stats(self):
        sn = bs.parse_samtools_stats(STATS)
        self.assertEqual(sn['raw total sequences'], 100_000_000)
        self.assertEqual(sn['bases mapped (cigar)'], 14_300_000_000)
        self.assertAlmostEqual(sn['error rate'], 0.003496503)
        self.assertEqual(sn['average length'], 150)
        self.assertNotIn('FFQ', sn)

    def test_idxstats(self):
        rows = bs.parse_idxstats(IDXSTATS)
        t = bs.idxstats_totals(rows)
        self.assertEqual(t['unplaced_unmapped'], 1_200_000, 'the * row = reads without a position (both mates unmapped)')
        self.assertEqual(t['records'], 5_010_000 + 5_512_000 + 1_005 + 1_200_000)
        self.assertEqual(t['mapped'], 5_000_000 + 5_500_000 + 1_000)

    def test_picard_metrics(self):
        libs = bs.parse_picard_metrics(PICARD)
        self.assertEqual(len(libs), 1)
        self.assertEqual(libs[0]['READ_PAIRS_EXAMINED'], 48_000_000)
        self.assertAlmostEqual(libs[0]['PERCENT_DUPLICATION'], 0.09)

    def test_picard_percent_duplication_is_reproduced_by_the_formula(self):
        """(UNPAIRED_READ_DUPLICATES + 2 * READ_PAIR_DUPLICATES) / (UNPAIRED_READS_EXAMINED + 2 * READ_PAIRS_EXAMINED) must equal Picard's own PERCENT_DUPLICATION."""
        lib = bs.parse_picard_metrics(PICARD)[0]
        s = bs.duplication_summary([lib])
        self.assertEqual(s['examined_reads'], 97_000_000)
        self.assertEqual(s['duplicate_reads'], 8_730_000)
        self.assertAlmostEqual(s['fraction'], lib['PERCENT_DUPLICATION'], places=9)
        self.assertEqual(s['optical_reads'], 600_000)

    def test_several_libraries_are_pooled_by_reads_examined_not_averaged(self):
        a = dict(UNPAIRED_READS_EXAMINED=0, READ_PAIRS_EXAMINED=1_000_000, UNPAIRED_READ_DUPLICATES=0, READ_PAIR_DUPLICATES=500_000, READ_PAIR_OPTICAL_DUPLICATES=0, UNMAPPED_READS=0)   # 50 %
        b = dict(UNPAIRED_READS_EXAMINED=0, READ_PAIRS_EXAMINED=9_000_000, UNPAIRED_READ_DUPLICATES=0, READ_PAIR_DUPLICATES=900_000, READ_PAIR_OPTICAL_DUPLICATES=0, UNMAPPED_READS=0)     # 10 %
        s = bs.duplication_summary([a, b])
        self.assertAlmostEqual(s['fraction'], 1_400_000 / 10_000_000)
        self.assertNotAlmostEqual(s['fraction'], (0.5 + 0.1) / 2)


class TestDerived(unittest.TestCase):
    def setUp(self):
        self.sn = bs.parse_samtools_stats(STATS)
        self.idx = bs.idxstats_totals(bs.parse_idxstats(IDXSTATS))
        self.libs = bs.parse_picard_metrics(PICARD)

    def test_new_columns(self):
        c = bs.derive_columns(self.sn, self.idx, self.libs, reference_size=962_000_000, ubam_reads=100_000_000)
        self.assertEqual(c['PrimaryReads'], 100_000_000)
        self.assertAlmostEqual(c['MappedFraction'], 0.97)
        self.assertAlmostEqual(c['UnmappedFraction'], 0.03)
        self.assertAlmostEqual(c['UnmappedNoPositionFraction'], 0.012, msg='both mates unmapped: no position in the BAM')
        self.assertAlmostEqual(c['UnmappedOrphanFraction'], 0.018, msg='unmapped but placed at the mate: what unmapped.bam draws from')
        self.assertAlmostEqual(c['UnmappedNoPositionFraction'] + c['UnmappedOrphanFraction'], c['UnmappedFraction'])
        self.assertAlmostEqual(c['DuplicateFraction'], 0.09)
        self.assertAlmostEqual(c['DuplicateFractionFromBam'], 8_730_000 / 97_000_000)
        self.assertAlmostEqual(c['OpticalDuplicateFraction'], 600_000 / 97_000_000)
        self.assertAlmostEqual(c['CoverageMapped'], 14_300_000_000 / 962_000_000)
        self.assertAlmostEqual(c['CoverageNonDuplicate'], (14_300_000_000 - 1_309_500_000) / 962_000_000)
        self.assertEqual(c['QCWarnings'], '', c['QCWarnings'])
        self.assertEqual(c['PipelineVersion'], 'align_pipeline')
        self.assertEqual(bs.derive_columns(self.sn, self.idx, self.libs, 962_000_000, pipeline_version='align_pipeline 3.0 abc1234')['PipelineVersion'], 'align_pipeline 3.0 abc1234')

    def test_independent_measurements_that_disagree_raise_warnings(self):
        sn = dict(self.sn, **{'reads duplicated': 11_640_000})                                       # flags say 12 %, Picard says 9 %
        w = bs.derive_columns(sn, self.idx, self.libs, 962_000_000)['QCWarnings']
        self.assertIn('duplicate rate differs', w)
        w = bs.derive_columns(self.sn, self.idx, self.libs, 962_000_000, ubam_reads=99_000_000)['QCWarnings']
        self.assertIn('unmapped BAM', w)
        w = bs.derive_columns(self.sn, dict(self.idx, unplaced_unmapped=2_400_000), self.libs, 962_000_000)['QCWarnings']
        self.assertIn('fully unmapped reads differ', w)

    def test_without_picard_metrics_the_flag_based_rate_is_used_and_flagged(self):
        c = bs.derive_columns(self.sn, self.idx, [], 962_000_000)
        self.assertAlmostEqual(c['DuplicateFraction'], c['DuplicateFractionFromBam'])
        self.assertIn('no Picard duplication metrics', c['QCWarnings'])
        self.assertTrue(math.isnan(c['OpticalDuplicateFraction']))

    def test_the_two_kept_legacy_columns_keep_their_original_meaning(self):
        """TotalReads = every record in all.bam; Coverage = TotalReads * (ReadLength / 2) / genome size. The six split fractions are gone (decisions 37 and 41)."""
        c = bs.legacy_columns(dict(all=101_000_000), read_length=300, reference_size=962_000_000)
        self.assertEqual(set(c), {'TotalReads', 'Coverage'})
        self.assertEqual(c['TotalReads'], 101_000_000)
        self.assertAlmostEqual(c['Coverage'], 101_000_000 * 150 / 962_000_000)
        self.assertTrue(math.isnan(bs.legacy_columns(dict(all=5), None, 962_000_000)['Coverage']), 'no read length: no coverage, not a wrong one')
        self.assertEqual((bs.LEGACY_BAM_TYPES, bs.LEGACY_COLUMNS), (('all',), ('TotalReads', 'Coverage')))


if __name__ == '__main__':
    unittest.main()
