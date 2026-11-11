"""The published inversion table (SPEC decision 70): Kumar et al. 2025, Supplementary File 4, for LG2, 9, 10, 11a/11b, 13 and 20a/20b; LG5 is this pipeline's own call until its numbers are in."""
import os, shutil, tempfile, unittest
import pandas as pd
from landmarks import inversions_table as T

# Kumar et al. 2025, eLife, Supplementary File 4, exactly as given (chromosome, breakpoint 1 interval start and end, breakpoint 2 interval start and end)
KUMAR = {'2': (19751534, 19845509, 43556655, 43559056), '9': (14912159, 14966448, 32492000, 32996230), '10': (11459086, 11818098, 29881276, 29902559),
         '11a': (8281988, 8303032, 18425177, 18425215), '11b': (18425177, 18425215, 30372599, 30387602), '13': (2398108, 2412755, 23074431, 23124401),
         '20a': (19616377, 19623494, 29682528, 29698815), '20b': (29682528, 29698815, 33170004, 33254382)}
LG5 = (19657762, 19657848, 26726676, 26726685)      # this pipeline's call, run 20261110-183552_aa6580ae, high-confidence consensus row (LatrZW:A, MezeOBm:A): left_min, left_max, right_min, right_max
CHROM = {'LG2': 'NC_135168.1', 'LG5': 'NC_135171.1', 'LG9': 'NC_135175.1', 'LG10': 'NC_135176.1', 'LG11': 'NC_135177.1', 'LG13': 'NC_135179.1', 'LG20': 'NC_135186.1'}


def raw():
    return pd.read_csv(T.table_path(), sep='\t', comment='#', dtype=str, keep_default_na=False)


def with_edit(fn):
    """The shipped table written to a temp file after `fn(DataFrame)` edited it."""
    d = raw(); fn(d)
    p = os.path.join(tempfile.mkdtemp(), 't.tsv'); d.to_csv(p, sep='\t', index=False)
    return p


class TestShippedTable(unittest.TestCase):
    def test_every_kumar_interval_is_exactly_as_published(self):
        d = raw().set_index('part')
        for part, want in KUMAR.items():
            got = tuple(int(d.at[part, c]) for c in T.COORDS)
            self.assertEqual(got, want, f'{part} differs from Supplementary File 4')
            self.assertEqual((d.at[part, 'status'], d.at[part, 'source']), ('published', 'Kumar2025 Supplementary File 4'))
            self.assertEqual(d.at[part, 'reference_state'], 'ancestral')

    def test_the_table_has_these_inversions_and_no_others(self):
        self.assertEqual(sorted(raw().part, key=lambda s: (int(''.join(ch for ch in s if ch.isdigit())), s)), ['2', '5', '9', '10', '11a', '11b', '13', '20a', '20b'])

    def test_each_row_is_on_its_chromosome(self):
        for r in raw().itertuples():
            self.assertEqual(r.chrom, CHROM[r.lg], r.part)

    def test_the_file_is_sound(self):
        self.assertEqual(T.validate(raw()), [])

    def test_the_tandem_pairs_share_their_breakpoint_exactly(self):
        d = raw().set_index('part')
        for a, b in (('11a', '11b'), ('20a', '20b')):
            self.assertEqual((d.at[a, 'bp2_start'], d.at[a, 'bp2_end']), (d.at[b, 'bp1_start'], d.at[b, 'bp1_end']))

    def test_lg5_is_the_pipelines_call_as_recorded(self):
        r = raw().set_index('part').loc['5']
        self.assertEqual((r.source, r.status, r.reference_state), ('inversion_breakpoints (LatrZW, MezeOBm)', 'pipeline', 'unset'))
        self.assertEqual(tuple(int(r[c]) for c in T.COORDS), LG5)
        self.assertTrue(0 < LG5[0] <= LG5[1] <= LG5[2] <= LG5[3] <= 42_954_235, 'inside NC_135171.1 (42,954,235 bp)')

    def test_the_comparison_record_is_in_the_repository(self):
        c = pd.read_csv(os.path.join(T.DATA, 'comparison_kumar_vs_pipeline_20261110.tsv'), sep='\t')
        self.assertEqual(len(c), 12); self.assertEqual(set(c.inversion), {'9', '10', '13', '11a/11b', '20a/20b'})


class TestLoader(unittest.TestCase):
    def test_the_shipped_table_loads_with_no_pending_row(self):
        d = T.load_inversions()
        self.assertEqual(len(d), 9); self.assertFalse((d.status == 'pending').any())
        l5 = d.set_index('part').loc['5']
        self.assertEqual((int(l5.length_min), int(l5.length_max)), (26726676 - 19657848, 26726685 - 19657762))

    def test_a_pending_row_is_refused_by_default_and_is_returned_when_asked_for(self):
        def pend(d):
            i = d.index[d.part == '5'][0]
            for c in T.COORDS:
                d.at[i, c] = ''
            d.at[i, 'status'] = 'pending'
        p = with_edit(pend)
        with self.assertRaises(ValueError) as cm:
            T.load_inversions(path=p)
        self.assertIn('no breakpoints yet (5)', str(cm.exception))
        d = T.load_inversions(path=p, allow_pending=True)
        self.assertTrue(d.loc[d.part == '5', 'bp1_start'].isna().all()); self.assertEqual(len(d), 9)

    def test_integer_coordinates_and_the_derived_columns(self):
        d = T.load_inversions(allow_pending=True).set_index('part')
        self.assertEqual(int(d.at['13', 'bp1_start']), 2398108); self.assertEqual(int(d.at['13', 'bp1_mid']), (2398108 + 2412755) // 2)
        self.assertEqual((int(d.at['20b', 'length_min']), int(d.at['20b', 'length_max'])), (3471189, 3571854), "20b: 3.47 to 3.57 Mb from the supplement's own intervals, not the 4.1 Mb of the text")
        self.assertEqual((int(d.at['20a', 'length_min']), int(d.at['20a', 'length_max'])), (10059034, 10082438))
        for p in ('2', '9', '10', '11a', '11b', '13', '20a', '20b'):
            self.assertLess(d.at[p, 'length_min'], d.at[p, 'length_max'])
            self.assertGreater(d.at[p, 'length_min'], 3_000_000); self.assertLess(d.at[p, 'length_max'], 25_000_000)

    def test_the_two_parts_of_a_tandem_inversion_add_up_to_the_whole(self):
        d = T.load_inversions(allow_pending=True).set_index('part')
        for a, b in (('11a', '11b'), ('20a', '20b')):
            self.assertEqual(d.at[a, 'bp2_mid'], d.at[b, 'bp1_mid'])

    def test_an_unsound_table_is_refused_with_the_reason(self):
        cases = {
            'a changed tandem breakpoint': (lambda d: d.__setitem__('bp1_start', d.bp1_start.where(d.part != '11b', '18425100')), 'tandem pair'),
            'an interval that starts after it ends': (lambda d: d.__setitem__('bp1_end', d.bp1_end.where(d.part != '13', '2398000')), 'starts after it ends'),
            'overlapping breakpoints': (lambda d: d.__setitem__('bp2_start', d.bp2_start.where(d.part != '9', '14912000')), 'ends after the right one starts'),
            'the wrong chromosome': (lambda d: d.__setitem__('chrom', d.chrom.where(d.part != '10', 'NC_135175.1')), 'is not the accession of LG10'),
            'a duplicated part': (lambda d: d.__setitem__('part', d.part.where(d.part != '9', '2')), 'duplicate part names'),
            'a bad status': (lambda d: d.__setitem__('status', d.status.where(d.part != '9', 'maybe')), 'status'),
            'no source': (lambda d: d.__setitem__('source', d.source.where(d.part != '9', '')), 'no source'),
            'a word in a coordinate': (lambda d: d.__setitem__('bp2_end', d.bp2_end.where(d.part != '9', '33M')), 'whole numbers'),
            'a pending row with its numbers still in': (lambda d: d.__setitem__('status', d.status.where(d.part != '5', 'pending')), 'pending row must have empty'),
        }
        for name, (fn, why) in cases.items():
            with self.assertRaises(ValueError, msg=name) as cm:
                T.load_inversions(path=with_edit(fn), allow_pending=True)
            self.assertIn(why, str(cm.exception), name)


class TestWindows(unittest.TestCase):
    def setUp(self):
        self.inv = T.load_inversions()

    def test_one_window_per_breakpoint_and_the_tandem_breakpoint_once(self):
        w = T.breakpoint_windows(self.inv)
        self.assertEqual(len(w), 16, '2, 5, 9, 10 and 13 have two each; 11 and 20 have three (left, internal, right)')
        self.assertEqual(sorted(n for n in w.name if 'internal' in n), ['inv11_internal', 'inv20_internal'])
        self.assertEqual(len(set(w.name)), len(w))
        for n in ('inv2_bp1', 'inv2_bp2', 'inv5_bp1', 'inv5_bp2', 'inv11a_bp1', 'inv11b_bp2', 'inv20a_bp1', 'inv20b_bp2', 'inv13_bp1'):
            self.assertIn(n, set(w.name))

    def test_windows_are_the_published_intervals(self):
        w = T.breakpoint_windows(self.inv).set_index('name')
        self.assertEqual((w.at['inv13_bp1', 'chrom'], w.at['inv13_bp1', 'start'], w.at['inv13_bp1', 'end']), ('NC_135179.1', 2398108, 2412755))
        self.assertEqual((w.at['inv11_internal', 'start'], w.at['inv11_internal', 'end']), (18425177, 18425215))
        self.assertEqual((w.at['inv20_internal', 'start'], w.at['inv20_internal', 'end']), (29682528, 29698815))

    def test_pad_widens_every_window_and_never_goes_below_zero(self):
        w0 = T.breakpoint_windows(self.inv).set_index('name'); w = T.breakpoint_windows(self.inv, pad=10_000).set_index('name')
        self.assertTrue(((w.start == w0.start - 10_000) | (w.start == 0)).all()); self.assertTrue((w.end == w0.end + 10_000).all())
        tiny = self.inv.copy(); tiny.loc[tiny.part == '13', 'bp1_start'] = 5
        self.assertEqual(int(T.breakpoint_windows(tiny, pad=1_000).set_index('name').at['inv13_bp1', 'start']), 0)

    def test_a_pending_row_has_no_windows(self):
        d = self.inv.copy(); i = d.index[d.part == '5'][0]
        d.loc[i, T.COORDS] = pd.NA; d.loc[i, 'status'] = 'pending'
        w = T.breakpoint_windows(d)
        self.assertEqual(len(w), 14); self.assertFalse(any(n.startswith('inv5_') for n in w.name))


if __name__ == '__main__':
    unittest.main()
