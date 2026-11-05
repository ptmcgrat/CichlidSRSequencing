"""AppendingDatabaseWriter against a fake worksheet that behaves like the Google Sheets API: append is atomic, whole-sheet rewrites are not."""
import os, sys, threading, time, unittest
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from helper_modules.alignment_db import AppendingDatabaseWriter, LEGACY_KEYS


class FakeSheet:
    """Rows are lists of cell values; row 1 is the header. Several views (= several machines) can share the same rows."""
    def __init__(self, header, rows=(), col_count=None, shared=None):
        self.shared = shared or dict(cells=[list(header)] + [list(r) for r in rows], lock=threading.Lock(), appended=0)
        self.col_count = col_count or len(header)
        self.fail_appends, self.silent_loss, self.inject_before_append = 0, 0, None

    def view(self):
        return FakeSheet(None, shared=self.shared, col_count=self.col_count)

    def row_values(self, n):
        with self.shared['lock']:
            return list(self.shared['cells'][n - 1]) if n <= len(self.shared['cells']) else []

    def col_values(self, c):
        with self.shared['lock']:
            return [(r[c - 1] if c - 1 < len(r) else '') for r in self.shared['cells']]

    def append_rows(self, values, value_input_option=None, table_range=None, insert_data_option=None):
        assert value_input_option == 'RAW' and table_range == 'A1'
        if self.inject_before_append:
            self.inject_before_append(); self.inject_before_append = None
        if self.fail_appends:
            self.fail_appends -= 1
            raise RuntimeError('APIError: 429 quota exceeded')
        if self.silent_loss:
            self.silent_loss -= 1
            return
        with self.shared['lock']:
            self.shared['cells'].extend([list(v) for v in values]); self.shared['appended'] += len(values)

    def add_cols(self, n):
        self.col_count += n

    def update_cell(self, r, c, v):
        with self.shared['lock']:
            row = self.shared['cells'][r - 1]
            row.extend([''] * (c - len(row))); row[c - 1] = v


def writer(sheet, **kw):
    fm = type('FM', (), {'genome_version': 'Mzebra_GT3_NCBI'})()
    kw.setdefault('sleep', lambda s: None); kw.setdefault('backoff', 0); kw.setdefault('log', lambda m: None)
    return AppendingDatabaseWriter(fm, worksheet_fn=lambda: sheet, **kw)


def row(sample, **extra):
    return dict({k: 1 for k in LEGACY_KEYS}, SampleID=sample, **extra)


HEADER = ['SampleID', 'GenomeVersion', 'Notes', 'RunIDs', 'Coverage', 'TotalReads', 'UnmappedReads', 'DiscordantReads', 'InversionReads', 'DuplicationReads', 'ClippedReads',
          'ChimericReads', 'minimap2_version', 'gatk_version', 'pysam_version', 'BamSize', 'PrimaryReads', 'DuplicateFraction']


class TestAppend(unittest.TestCase):
    def test_the_row_follows_the_sheets_column_order_blank_where_unknown_and_python_types(self):
        sh = FakeSheet(HEADER, [['OLD', 'G', 'a note', 'r'] + [''] * 14])
        out = writer(sh).add_row(row('NEW', PrimaryReads=np.int64(7), DuplicateFraction=np.float64('nan')))
        self.assertEqual(out, 'added')
        new = sh.shared['cells'][-1]
        self.assertEqual(new[0], 'NEW'); self.assertEqual(new[HEADER.index('Notes')], '', 'a sheet column we know nothing about stays blank')
        self.assertEqual(new[HEADER.index('PrimaryReads')], 7); self.assertIsInstance(new[HEADER.index('PrimaryReads')], int)
        self.assertEqual(new[HEADER.index('DuplicateFraction')], '', 'NaN is written as an empty cell')
        self.assertEqual(len(sh.shared['cells']), 3)

    def test_an_existing_sample_is_left_alone(self):
        sh = FakeSheet(HEADER, [['DONE'] + [''] * 17])
        self.assertEqual(writer(sh).add_row(row('DONE')), 'exists')
        self.assertEqual(sh.shared['appended'], 0)

    def test_a_missing_column_is_refused_with_the_setup_instruction_and_nothing_is_written(self):
        sh = FakeSheet(HEADER[:-2])
        with self.assertRaises(RuntimeError) as cm:
            writer(sh).add_row(row('NEW', PrimaryReads=1, DuplicateFraction=0.1))
        self.assertIn('--setup-database', str(cm.exception)); self.assertEqual(sh.shared['appended'], 0)

    def test_setup_adds_the_new_columns_once_and_resizes_the_sheet(self):
        sh = FakeSheet(HEADER[:16], col_count=16)
        w = writer(sh)
        added = w.setup_columns()
        self.assertIn('PrimaryReads', added); self.assertIn('QCWarnings', added)
        header = sh.row_values(1)
        self.assertTrue(all(c in header for c in AppendingDatabaseWriter.needed_columns()))
        self.assertGreaterEqual(sh.col_count, len(header), 'the worksheet must be wide enough for the new header')
        self.assertEqual(w.setup_columns(), [], 'a second run adds nothing')
        self.assertEqual(header[:16], HEADER[:16], 'existing columns are never moved')

    def test_check_columns_reports_what_the_header_lacks(self):
        sh = FakeSheet(HEADER)
        missing = writer(sh).check_columns()
        self.assertIn('QCWarnings', missing); self.assertNotIn('PrimaryReads', missing)
        writer(sh).setup_columns()
        self.assertEqual(writer(sh).check_columns(), [])

    def test_two_machines_appending_at_the_same_time_lose_nothing(self):
        sh = FakeSheet(HEADER)
        machines = [writer(sh.view()), writer(sh.view())]
        errors = []

        def work(w, prefix):
            try:
                for i in range(60):
                    w.add_row(row(f'{prefix}{i:03d}', PrimaryReads=i))
            except Exception as e:
                errors.append(e)
        ts = [threading.Thread(target=work, args=(machines[0], 'A')), threading.Thread(target=work, args=(machines[1], 'B'))]
        [t.start() for t in ts]; [t.join() for t in ts]
        ids = [r[0] for r in sh.shared['cells'][1:]]
        self.assertEqual(errors, [])
        self.assertEqual((len(ids), len(set(ids))), (120, 120), 'every row from both machines is there exactly once')

    def test_the_old_whole_sheet_rewrite_does_lose_rows_in_the_same_race(self):
        """What the append replaces: read the sheet, add a row, write the WHOLE sheet back. Two machines whose read-modify-write overlap: the later write erases the earlier row."""
        sh = FakeSheet(HEADER)
        gate = threading.Barrier(2)

        def old_writer(sample):
            with sh.shared['lock']:
                snapshot = [list(r) for r in sh.shared['cells']]                      # read
            gate.wait()                                                               # both machines have read before either writes
            snapshot.append([sample] + [''] * 17)
            with sh.shared['lock']:
                sh.shared['cells'][:] = snapshot                                      # write the whole sheet back
        ts = [threading.Thread(target=old_writer, args=(s,)) for s in ('M1', 'M2')]
        [t.start() for t in ts]; [t.join() for t in ts]
        self.assertEqual(len(sh.shared['cells']) - 1, 1, 'one of the two rows is gone: the lost-row race the append avoids')

    def test_transient_api_errors_are_retried_and_a_silently_lost_append_is_detected(self):
        sh = FakeSheet(HEADER); sh.fail_appends = 2
        self.assertEqual(writer(sh).add_row(row('S1')), 'added'); self.assertEqual([r[0] for r in sh.shared['cells'][1:]], ['S1'])
        sh2 = FakeSheet(HEADER); sh2.silent_loss = 1                                  # the call "succeeds" but the row never appears
        self.assertEqual(writer(sh2).add_row(row('S2')), 'added'); self.assertEqual([r[0] for r in sh2.shared['cells'][1:]], ['S2'])
        sh3 = FakeSheet(HEADER); sh3.fail_appends = 10
        with self.assertRaises(RuntimeError) as cm:
            writer(sh3, attempts=3).add_row(row('S3'))
        self.assertIn('after 3 attempts', str(cm.exception))

    def test_a_row_added_by_another_machine_at_the_same_moment_is_reported_not_hidden(self):
        sh = FakeSheet(HEADER)
        other = sh.view()
        sh.inject_before_append = lambda: other.append_rows([['DUP'] + [''] * 17], 'RAW', 'A1')       # the other machine appends the same sample just before ours
        logs = []
        self.assertEqual(writer(sh, log=logs.append).add_row(row('DUP')), 'added')
        self.assertTrue(any('appears 2 times' in m for m in logs), logs)


if __name__ == '__main__':
    unittest.main()
