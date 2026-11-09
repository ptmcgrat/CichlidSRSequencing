"""The slice table: END means the whole chromosome, errors are loud, and the table that ships (decision 42) is LG2, LG10 and LG13 as whole chromosomes."""
import os, sys, tempfile, unittest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from align_pipeline import bam_slices as sl

SHIPPED = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'slice_regions', 'Mzebra_GT3_NCBI.tsv')
LENGTHS = {'NC_135168.1': 45_000_000, 'NC_135176.1': 34_949_860, 'NC_135179.1': 30_000_000, 'NC_135171.1': 40_000_000}


def table(text):
    p = tempfile.mkdtemp() + '/slices.tsv'
    open(p, 'w').write(text)
    return p


def fai(lengths=LENGTHS):
    p = tempfile.mkdtemp() + '/ref.fa.fai'
    open(p, 'w').write(''.join(f'{k}\t{v}\t0\t60\t61\n' for k, v in lengths.items()))
    return p


class TestSlices(unittest.TestCase):
    def test_end_is_the_length_of_the_chromosome_from_the_index(self):
        s = sl.load_slices(table('A\tNC_135176.1\t1\tEND\tx\nB\tNC_135168.1\t100\t200\n'), fai=fai())
        self.assertEqual([(x.name, x.start, x.end, x.region) for x in s], [('A', 1, 34949860, 'NC_135176.1:1-34949860'), ('B', 100, 200, 'NC_135168.1:100-200')])
        self.assertEqual(sl.load_slices(table('A\tNC_135176.1\t1\tend\n'), fai=LENGTHS)[0].end, 34949860, 'END is case-insensitive and a dict of lengths works as well as a file')

    def test_end_without_an_index_or_with_an_unknown_sequence_is_an_error(self):
        with self.assertRaisesRegex(ValueError, 'END needs the reference index'):
            sl.load_slices(table('A\tNC_135176.1\t1\tEND\n'))
        with self.assertRaisesRegex(ValueError, 'not in the reference index'):
            sl.load_slices(table('A\tNC_999.1\t1\tEND\n'), fai=fai())

    def test_other_errors_are_still_loud(self):
        for text, why in (('A\tc\t1\n', 'need name'), ('A b\tc\t1\t5\n', 'may only contain'), ('A\tc\t5\t5\n', 'need 1 <= start < end'), ('A\tc\t1\t5\nA\tc\t1\t6\n', 'duplicate')):
            with self.assertRaisesRegex(ValueError, why):
                sl.load_slices(table(text))

    def test_check_against_fai_reports_unknown_sequences_and_ends_past_the_chromosome(self):
        s = sl.load_slices(table('A\tNC_135176.1\t1\t99999999\nB\tNC_X\t1\t5\nC\tNC_135176.1\t1\t10\n'))
        problems = sl.check_against_fai(s, fai())
        self.assertEqual(len(problems), 2); self.assertTrue(any('past the end' in p for p in problems) and any('not in the reference' in p for p in problems))

    def test_the_shipped_table_is_lg2_lg10_and_lg13_as_whole_chromosomes(self):
        s = sl.load_slices(SHIPPED, fai=fai())
        self.assertEqual([x.name for x in s], ['LG2', 'LG10', 'LG13'])
        self.assertEqual([x.chrom for x in s], ['NC_135168.1', 'NC_135176.1', 'NC_135179.1'])
        self.assertEqual([(x.start, x.end) for x in s], [(1, LENGTHS['NC_135168.1']), (1, 34949860), (1, LENGTHS['NC_135179.1'])])
        self.assertEqual(sl.check_against_fai(s, fai()), [])
        self.assertFalse({'LG9', 'LG11', 'LG20'} & {x.name for x in s}, 'decision 42: LG9, LG11 and LG20 are not sliced')

    def test_commands(self):
        s = sl.load_slices(table('A\tNC_135176.1\t1\tEND\n'), fai=fai())[0]
        self.assertEqual(sl.view_command('/d/S.all.bam', '/d/out.bam', s), ['samtools', 'view', '-b', '-@', '2', '-o', '/d/out.bam', '/d/S.all.bam', 'NC_135176.1:1-34949860'])
        self.assertEqual(sl.slice_path('/d/slices/', 'S1', s), '/d/slices/S1.slice.A.bam')


if __name__ == '__main__':
    unittest.main()
