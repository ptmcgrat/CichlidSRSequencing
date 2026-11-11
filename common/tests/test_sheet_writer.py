"""cichlid_common.sheet_writer against a fake worksheet: the plan, what is and is not written, concurrent edits, verification, new columns, the journal."""
import os, tempfile, unittest
from cichlid_common import sheet_writer as SW
from cichlid_common.fake_sheet import FakeWorksheet

HEADER = ['SampleID', 'GenomeVersion', 'Coverage', 'UsableCoverage', 'Notes']
ROWS = [['A', 'g', '9.9', '6.8', 'hand note'], ['B', 'g', '19.9', '', ''], ['C', 'g', '14.6', '', ''], ['D', 'g', '10', '', '']]


def updater(ws, **kw):
    return SW.ColumnUpdater(lambda: ws, log=lambda m: None, **kw)


class TestHelpers(unittest.TestCase):
    def test_column_letters(self):
        self.assertEqual([SW.col_letter(n) for n in (1, 26, 27, 52, 53, 702, 703)], ['A', 'Z', 'AA', 'AZ', 'BA', 'ZZ', 'AAA'])
        self.assertEqual(SW.a1(17, 13), 'M17')
        with self.assertRaises(ValueError):
            SW.col_letter(0)

    def test_a_number_matches_the_rounding_its_display_shows(self):
        self.assertTrue(SW.same_value('4.0123', 4.0123)); self.assertTrue(SW.same_value('4.01', 4.0123)); self.assertTrue(SW.same_value('4', 4.0123 - 0.4))
        self.assertFalse(SW.same_value('4.02', 4.0123)); self.assertFalse(SW.same_value('', 4.0)); self.assertFalse(SW.same_value('abc', 4.0))
        self.assertTrue(SW.same_value('72.6%', 0.726)); self.assertTrue(SW.same_value('1,234.5', 1234.5)); self.assertTrue(SW.same_value(' text ', 'text')); self.assertFalse(SW.same_value('text', 'other'))
        self.assertFalse(SW.same_value('1', True), 'a bool is not a number here')


class TestPlan(unittest.TestCase):
    def test_every_status(self):
        ws = FakeWorksheet(HEADER, ROWS + [['D', 'g', '1', '', '']])                       # D twice
        ws.grid[2][3] = '5.0'                                                              # B's UsableCoverage already holds something else
        ws.grid[3][3] = '3.8'                                                              # C's already holds the value we would write
        u = {'A': {'UsableCoverage': 3.8}, 'B': {'UsableCoverage': 4.2}, 'C': {'UsableCoverage': 3.8}, 'D': {'UsableCoverage': 1.0}, 'Z': {'UsableCoverage': 1.0}, 'E?': {'UsableCoverage': 1}}
        p = updater(ws).plan(u)
        st = {c['sample']: c['status'] for c in p.changes}
        self.assertEqual(st, {'A': 'conflict', 'B': 'conflict', 'C': 'same', 'D': 'duplicate_row', 'Z': 'no_row', 'E?': 'no_row'})
        self.assertEqual(updater(ws).plan({'A': {'UsableCoverage': 3.8}}, overwrite=True).changes[0]['status'], 'overwrite')
        self.assertEqual(updater(ws).plan({'B': {'Notes': 'x'}}).changes[0]['status'], 'write', 'an empty cell is written')

    def test_a_plan_writes_nothing(self):
        ws = FakeWorksheet(HEADER, ROWS)
        p = updater(ws).plan({'B': {'EffectiveCoverage': 4.0, 'Notes': 'n'}, 'C': {'EffectiveCoverage': 5.0}})
        self.assertEqual(ws.writes(), []); self.assertEqual(p.header_missing, ['EffectiveCoverage'])
        self.assertEqual([c['status'] for c in p.changes], ['write', 'write', 'write'])
        self.assertIn('columns to ADD at the end of the header: EffectiveCoverage', p.format()); self.assertIn('first cells to write', p.format())

    def test_a_missing_column_is_planned_at_the_end_of_the_header(self):
        ws = FakeWorksheet(HEADER, ROWS)
        p = updater(ws).plan({'B': {'M1': 1.0, 'M2': 2.0}})
        self.assertEqual([(c['column'], c['col']) for c in p.changes], [('M1', 6), ('M2', 7)])

    def test_the_header_must_be_unambiguous(self):
        ws = FakeWorksheet(HEADER + ['Coverage'], [r + [''] for r in ROWS])
        with self.assertRaises(SW.SheetWriteError):
            updater(ws).plan({'A': {'Coverage': 1.0}})
        with self.assertRaises(SW.SheetWriteError):
            updater(FakeWorksheet(['Sample', 'x'], [['A', '1']])).plan({'A': {'x': 1}})
        with self.assertRaises(SW.SheetWriteError):
            updater(FakeWorksheet([], [])).plan({'A': {'x': 1}})

    def test_ids_are_matched_after_stripping(self):
        ws = FakeWorksheet(HEADER, [[' A ', 'g', '1', '', '']])
        self.assertEqual(updater(ws).plan({'A ': {'UsableCoverage': 1.0}}).changes[0]['status'], 'write')


class TestApply(unittest.TestCase):
    def test_only_the_named_cells_change(self):
        ws = FakeWorksheet(HEADER, ROWS, extra_columns=2)
        before = [list(r) for r in ws.grid]
        rep = updater(ws).apply({'B': {'Notes': 'measured', 'UsableCoverage': 4.5}, 'C': {'UsableCoverage': 5.5}})
        self.assertTrue(rep.ok); self.assertEqual(rep.written, 3); self.assertEqual(rep.added_columns, [])
        changed = [(r, c) for r in range(len(ws.grid)) for c in range(len(ws.grid[0])) if ws.grid[r][c] != before[r][c]]
        self.assertEqual(changed, [(2, 3), (2, 4), (3, 3)], 'B Usable, B Notes, C Usable and nothing else')
        self.assertEqual(ws.grid[2][3], 4.5); self.assertEqual(ws.grid[1][4], 'hand note', "A's hand-written note is untouched")
        self.assertEqual([c for c in ws.calls if c[0] == 'batch_update'], [('batch_update', 3, 'RAW')])

    def test_a_conflict_is_left_alone_unless_overwrite(self):
        ws = FakeWorksheet(HEADER, ROWS)
        rep = updater(ws).apply({'A': {'UsableCoverage': 3.8}})
        self.assertEqual((rep.written, ws.grid[1][3]), (0, '6.8')); self.assertEqual(rep.plan.count('conflict'), 1); self.assertEqual(ws.writes(), [], 'not even a call when there is nothing to write')
        rep = updater(ws).apply({'A': {'UsableCoverage': 3.8}}, overwrite=True)
        self.assertEqual((rep.written, ws.grid[1][3]), (1, 3.8))

    def test_new_columns_are_added_at_the_end_then_filled(self):
        ws = FakeWorksheet(HEADER, ROWS)
        rep = updater(ws).apply({'B': {'EffectiveCoverage': 4.0, 'MapqPassFraction': 0.8}, 'C': {'EffectiveCoverage': 5.0}})
        self.assertEqual(rep.added_columns, ['EffectiveCoverage', 'MapqPassFraction']); self.assertEqual(ws.grid[0], HEADER + ['EffectiveCoverage', 'MapqPassFraction'])
        self.assertEqual(ws.column('EffectiveCoverage'), ['', 4.0, 5.0, '']); self.assertEqual(ws.column('MapqPassFraction'), ['', 0.8, '', ''])
        self.assertTrue(rep.ok)
        self.assertEqual(ws.column('SampleID'), ['A', 'B', 'C', 'D']); self.assertEqual(ws.column('Coverage'), ['9.9', '19.9', '14.6', '10'])

    def test_the_grid_is_widened_when_the_header_fills_it(self):
        ws = FakeWorksheet(HEADER, ROWS)
        self.assertEqual(ws.col_count, 5)
        updater(ws).apply({'B': {'NewCol': 1.0}})
        self.assertIn(('add_cols', 1), ws.calls); self.assertEqual(ws.column('NewCol'), ['', 1.0, '', ''])
        ws2 = FakeWorksheet(HEADER, ROWS, extra_columns=3)
        updater(ws2).apply({'B': {'NewCol': 1.0}}); self.assertFalse([c for c in ws2.calls if c[0] == 'add_cols'], 'room was already there')

    def test_a_second_run_is_a_no_op(self):
        ws = FakeWorksheet(HEADER, ROWS)
        u = {'B': {'EffectiveCoverage': 4.0123}, 'C': {'EffectiveCoverage': 5.0}}
        updater(ws).apply(u); n = len(ws.calls)
        rep = updater(ws).apply(u)
        self.assertEqual((rep.written, rep.plan.count('same')), (0, 2)); self.assertEqual(len(ws.calls), n, 'nothing was written, not even a header')

    def test_a_value_shown_rounded_still_counts_as_written(self):
        ws = FakeWorksheet(HEADER + ['EffectiveCoverage'], [r + [''] for r in ROWS], show_decimals={'EffectiveCoverage': 2})
        u = {'B': {'EffectiveCoverage': 4.0123}}
        rep = updater(ws).apply(u)
        self.assertTrue(rep.ok); self.assertEqual(ws.get_all_values()[2][5], '4.01'); self.assertEqual(updater(ws).apply(u).plan.count('same'), 1)

    def test_limit_keeps_the_first_samples_alphabetically(self):
        ws = FakeWorksheet(HEADER, ROWS)
        rep = updater(ws).apply({'D': {'Notes': 'd'}, 'B': {'Notes': 'b'}, 'C': {'Notes': 'c'}}, limit=2)
        self.assertEqual(rep.written, 2); self.assertEqual(ws.column('Notes'), ['hand note', 'b', 'c', ''])

    def test_batches_are_chunked(self):
        rows = [[f'S{i}', 'g', '1', '', ''] for i in range(25)]
        ws = FakeWorksheet(HEADER, rows)
        rep = updater(ws, chunk=10).apply({f'S{i}': {'Notes': i} for i in range(25)})
        self.assertEqual([c[1] for c in ws.calls if c[0] == 'batch_update'], [10, 10, 5]); self.assertEqual(rep.written, 25)


class TestWhileOthersWork(unittest.TestCase):
    def test_a_row_inserted_since_the_plan_does_not_misplace_a_write(self):
        ws = FakeWorksheet(HEADER, ROWS)
        u = {'B': {'UsableCoverage': 4.5}}
        up = updater(ws)
        up.plan(u)                                                                   # the dry run
        ws.insert_row_at_top(['NEW', 'g', '1', '', ''])                            # somebody inserts a row afterwards
        rep = up.apply(u)
        self.assertTrue(rep.ok); self.assertEqual(dict(zip(ws.column('SampleID'), ws.column('UsableCoverage')))['B'], 4.5)

    def test_a_cell_filled_by_someone_since_the_plan_is_not_overwritten(self):
        ws = FakeWorksheet(HEADER, ROWS)
        u = {'B': {'UsableCoverage': 4.5}, 'C': {'UsableCoverage': 5.5}}
        up = updater(ws); first = up.plan(u); self.assertEqual(first.count('write'), 2)
        ws.grid[2][3] = 'typed by hand'
        rep = up.apply(u)
        self.assertEqual((rep.written, rep.plan.count('conflict')), (1, 1)); self.assertEqual(ws.grid[2][3], 'typed by hand')

    def test_a_shift_between_the_write_and_the_check_is_reported_loudly(self):
        ws = FakeWorksheet(HEADER, ROWS)
        ws.after_write = lambda w: w.insert_row_at_top(['NEW', 'g', '1', '', ''])        # rows shift after the batch was sent
        rep = updater(ws).apply({'B': {'UsableCoverage': 4.5}})
        self.assertFalse(rep.ok); self.assertEqual(len(rep.failures), 1); self.assertIn('the sheet changed while writing', rep.failures[0]['problem'])

    def test_a_row_inserted_between_the_final_read_and_the_write_is_caught_by_the_check(self):
        ws = FakeWorksheet(HEADER, ROWS)
        ws.before_write = lambda w: w.insert_row_at_top(['NEW', 'g', '1', '', ''])                 # the unavoidable race: the write lands one row off
        rep = updater(ws).apply({'B': {'UsableCoverage': 4.5}})
        self.assertFalse(rep.ok); self.assertIn('the sheet changed while writing', rep.failures[0]['problem'])

    def test_the_write_uses_the_worksheet_that_was_just_read(self):
        opened = []
        ws = FakeWorksheet(HEADER, ROWS)
        def fn():
            opened.append(1); return ws
        SW.ColumnUpdater(fn, log=lambda m: None).apply({'B': {'UsableCoverage': 4.5}})
        self.assertEqual(len(opened), 3, 'the first plan, the second plan, the check: the write itself opens nothing')

    def test_a_write_that_did_not_take_is_reported(self):
        ws = FakeWorksheet(HEADER, ROWS)
        ws.after_write = lambda w: w.grid[2].__setitem__(3, '')                          # the value vanished
        rep = updater(ws).apply({'B': {'UsableCoverage': 4.5}})
        self.assertFalse(rep.ok); self.assertIn('reads back', rep.failures[0]['problem'])


class TestJournal(unittest.TestCase):
    def test_planned_then_verified_with_the_old_values(self):
        ws = FakeWorksheet(HEADER, ROWS)
        j = os.path.join(tempfile.mkdtemp(), 'sub', 'j.tsv')
        updater(ws).apply({'B': {'UsableCoverage': 4.5}, 'C': {'Notes': 'x'}}, journal=j)
        lines = [l.split('\t') for l in open(j).read().splitlines()]
        self.assertEqual(lines[0][1:], ['phase', 'sample', 'column', 'cell', 'old', 'new', 'note'])
        body = [l[1:] for l in lines[1:]]
        self.assertEqual(body, [['planned', 'B', 'UsableCoverage', 'D3', '', '4.5', 'write'], ['planned', 'C', 'Notes', 'E4', '', 'x', 'write'],
                                ['verified', 'B', 'UsableCoverage', 'D3', '', '4.5', ''], ['verified', 'C', 'Notes', 'E4', '', 'x', '']])

    def test_a_failure_is_in_the_journal(self):
        ws = FakeWorksheet(HEADER, ROWS); ws.after_write = lambda w: w.insert_row_at_top(['NEW', 'g', '1', '', ''])
        j = os.path.join(tempfile.mkdtemp(), 'j.tsv')
        updater(ws).apply({'B': {'UsableCoverage': 4.5}}, journal=j)
        self.assertIn('FAILED', open(j).read()); self.assertIn('the sheet changed while writing', open(j).read())


if __name__ == '__main__':
    unittest.main()
