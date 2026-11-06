"""Tests of helper_modules/download_manifest.py. The Google part is a fake client that has the gspread methods the backend calls; it checks OUR use of them (ranges, chunking, A1 notation),
not Google's behaviour."""
import os, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
from helper_modules import download_manifest as M


def runinfo(runs, **kw):
    return pd.DataFrame([{c: '' for c in M.RUNINFO_COLUMNS} | {'Run': r, 'BioSample': 'SAMEA' + r[-3:], 'BioProject': 'PRJEB1', 'Bases': '1000', 'Organism': 'Fish', 'Platform': 'ILLUMINA'} | kw for r in runs])


class Clock:
    def __init__(self, t=1_800_000_000.0): self.t = t
    def __call__(self): return self.t
    def sleep(self, dt): self.t += dt


def manifest(backend=None, who='utaka-1', clock=None, tab='Blumer', **kw):
    clock = clock or Clock()
    return M.Manifest(backend or M.MemoryBackend(), tab, who, now=clock, sleep=clock.sleep, log=lambda *a: None, **kw), clock


def cellv(b, tab, run, col):
    v = b.data[tab]; hdr = v[0]
    row = next(r for r in v[1:] if r[hdr.index('Run')] == run)
    return row[hdr.index(col)]


class TestHelpers(unittest.TestCase):
    def test_a1_and_cells(self):
        self.assertEqual([M.col_letters(n) for n in (1, 26, 27, 28, 52, 53)], ['A', 'Z', 'AA', 'AB', 'AZ', 'BA']); self.assertEqual(M.a1(3, 28), 'AB3')
        self.assertEqual((M.cell(True), M.cell(False), M.cell(None), M.cell(7), M.cell(12.3456), M.cell(1234567.0)), ('TRUE', 'FALSE', '', '7', '12.35', '1234567'))
        self.assertEqual(M.parse_iso(M.utc_iso(1_800_000_000)), 1_800_000_000); self.assertIsNone(M.parse_iso('')); self.assertIsNone(M.parse_iso('garbage'))


class TestStructureAndImport(unittest.TestCase):
    def test_setup_creates_the_header_and_adds_missing_columns_without_reordering(self):
        m, _ = manifest(); m.setup()
        self.assertEqual(m.b.data['Blumer'][0], M.RUNINFO_COLUMNS + M.TRACK_COLUMNS)
        m.b.data['Blumer'][0] = M.RUNINFO_COLUMNS + ['DownloadStatus', 'MyOwnColumn']          # an older / hand-edited sheet
        added = m.setup()
        self.assertEqual(m.b.data['Blumer'][0][:len(M.RUNINFO_COLUMNS) + 2], M.RUNINFO_COLUMNS + ['DownloadStatus', 'MyOwnColumn'])
        self.assertIn('InDatabase', added); self.assertNotIn('DownloadStatus', added); self.assertIn('Skip', m.b.data['Blumer'][0])

    def test_import_sets_status_from_the_database_and_the_cloud_and_is_idempotent(self):
        m, _ = manifest()
        r = m.import_runinfo(runinfo(['ERR001', 'ERR002', 'ERR003']), in_db={'ERR001'}, on_cloud={'ERR002'})
        self.assertEqual(r, dict(added=3, existing=0))
        b = m.b
        self.assertEqual([cellv(b, 'Blumer', x, 'DownloadStatus') for x in ('ERR001', 'ERR002', 'ERR003')], ['in_database', 'uploaded', 'pending'])
        self.assertEqual([cellv(b, 'Blumer', x, 'InDatabase') for x in ('ERR001', 'ERR002', 'ERR003')], ['TRUE', 'FALSE', 'FALSE'])
        self.assertEqual([cellv(b, 'Blumer', x, 'Source') for x in ('ERR001', 'ERR002', 'ERR003')], ['legacy', 'legacy', ''])
        hdr = b.data['Blumer'][0]; b.data['Blumer'][3][hdr.index('Notes')] = 'check with Rob'; b.data['Blumer'][3][hdr.index('DownloadStatus')] = 'failed'     # a person edited ERR003
        r2 = m.import_runinfo(runinfo(['ERR003', 'ERR004']), in_db=set(), on_cloud=set())
        self.assertEqual(r2, dict(added=1, existing=1)); self.assertEqual(len(b.data['Blumer']), 5)
        self.assertEqual((cellv(b, 'Blumer', 'ERR003', 'Notes'), cellv(b, 'Blumer', 'ERR003', 'DownloadStatus')), ('check with Rob', 'failed'))     # existing rows untouched


class TestUpdates(unittest.TestCase):
    def setUp(self):
        self.m, self.clock = manifest(); self.m.import_runinfo(runinfo(['ERR001', 'ERR002', 'ERR003']))

    def test_flush_writes_only_the_changed_cells_in_one_batch(self):
        b = self.m.b; b.calls.clear()
        self.m.update('ERR002', DownloadStatus='downloading', Machine='utaka', Attempts=2)
        self.assertTrue(self.m.flush(force=True))
        self.assertEqual(b.calls.count('write_cells'), 1)
        self.assertEqual([cellv(b, 'Blumer', 'ERR002', c) for c in ('DownloadStatus', 'Machine', 'Attempts')], ['downloading', 'utaka', '2'])
        self.assertEqual(cellv(b, 'Blumer', 'ERR001', 'DownloadStatus'), 'pending')                       # others untouched
        self.assertTrue(cellv(b, 'Blumer', 'ERR002', 'UpdatedAt').endswith('Z'))

    def test_rows_are_found_by_run_id_even_after_the_sheet_is_sorted(self):
        b = self.m.b; b.data['Blumer'][1:] = list(reversed(b.data['Blumer'][1:]))
        self.m.update('ERR001', DownloadStatus='uploaded'); self.m.flush(force=True)
        self.assertEqual(cellv(b, 'Blumer', 'ERR001', 'DownloadStatus'), 'uploaded'); self.assertEqual(cellv(b, 'Blumer', 'ERR003', 'DownloadStatus'), 'pending')

    def test_flush_is_rate_limited_unless_forced(self):
        self.m.update('ERR001', Notes='a'); self.assertTrue(self.m.flush(force=True)); n = self.m.b.calls.count('write_cells')
        self.m.update('ERR001', Notes='b'); self.m.flush()
        self.assertEqual(self.m.b.calls.count('write_cells'), n)
        self.clock.sleep(20); self.m.flush(); self.assertEqual(self.m.b.calls.count('write_cells'), n + 1); self.assertEqual(cellv(self.m.b, 'Blumer', 'ERR001', 'Notes'), 'b')

    def test_an_unreachable_sheet_journals_the_updates_and_the_next_flush_replays_them(self):
        j = tempfile.mkdtemp() + '/j.jsonl'; m, clock = manifest(journal=j); m.import_runinfo(runinfo(['ERR001', 'ERR002']))
        def down(name):
            if name in ('header', 'column', 'write_cells'): raise ConnectionError('google unreachable')
        m.b.fail = down
        m.update('ERR001', DownloadStatus='downloading'); self.assertFalse(m.flush(force=True))
        self.assertTrue(os.path.exists(j)); self.assertEqual(cellv(m.b, 'Blumer', 'ERR001', 'DownloadStatus'), 'pending')
        m.b.fail = None; m.update('ERR002', DownloadStatus='failed'); m.update('ERR001', DownloadStatus='uploaded')      # newer than the journal
        self.assertTrue(m.flush(force=True)); self.assertFalse(os.path.exists(j))
        self.assertEqual((cellv(m.b, 'Blumer', 'ERR001', 'DownloadStatus'), cellv(m.b, 'Blumer', 'ERR002', 'DownloadStatus')), ('uploaded', 'failed'))   # pending beats the older journal

    def test_without_a_journal_file_the_updates_stay_queued(self):
        m, _ = manifest(); m.import_runinfo(runinfo(['ERR001']))
        m.b.fail = lambda n: (_ for _ in ()).throw(ConnectionError('x')) if n == 'write_cells' else None
        m.update('ERR001', Notes='x'); self.assertFalse(m.flush(force=True)); self.assertIn('ERR001', m.pending)
        m.b.fail = None; m.clock = None; self.assertTrue(m.flush(force=True)); self.assertEqual(cellv(m.b, 'Blumer', 'ERR001', 'Notes'), 'x')

    def test_a_run_that_is_not_in_the_sheet_is_reported_not_fatal(self):
        logs = []; self.m.log = logs.append
        self.m.update('ERR999', Notes='x'); self.assertTrue(self.m.flush(force=True)); self.assertTrue(any('ERR999' in str(x) for x in logs))

    def test_a_new_column_is_added_on_demand(self):
        self.m.update('ERR001', BrandNew='v'); self.m.flush(force=True)
        self.assertIn('BrandNew', self.m.b.data['Blumer'][0]); self.assertEqual(cellv(self.m.b, 'Blumer', 'ERR001', 'BrandNew'), 'v')


class TestClaims(unittest.TestCase):
    def setUp(self):
        self.clock = Clock(); self.b = M.MemoryBackend()
        self.a, _ = manifest(self.b, 'utaka-1', self.clock); self.a.import_runinfo(runinfo(['ERR001', 'ERR002', 'ERR003', 'ERR004']))
        self.o, _ = manifest(self.b, 'workstation-9', self.clock)

    def test_claim_succeeds_once_and_blocks_the_other_machine_until_it_expires(self):
        self.assertTrue(self.a.claim('ERR001')); self.assertEqual(cellv(self.b, 'Blumer', 'ERR001', 'ClaimedBy'), 'utaka-1'); self.assertEqual(cellv(self.b, 'Blumer', 'ERR001', 'DownloadStatus'), 'claimed')
        self.assertFalse(self.o.claim('ERR001'))
        self.clock.sleep(7 * 3600); self.assertTrue(self.o.claim('ERR001')); self.assertEqual(cellv(self.b, 'Blumer', 'ERR001', 'ClaimedBy'), 'workstation-9')
        self.assertTrue(self.o.claim('ERR001'))                                                          # our own claim can be renewed

    def test_final_skipped_and_unknown_runs_cannot_be_claimed(self):
        hdr = self.b.data['Blumer'][0]
        for run, col, val in (('ERR002', 'DownloadStatus', 'in_database'), ('ERR003', 'Skip', 'TRUE')):
            next(r for r in self.b.data['Blumer'][1:] if r[0] == run)[hdr.index(col)] = val
        self.assertFalse(self.a.claim('ERR002')); self.assertFalse(self.a.claim('ERR003')); self.assertFalse(self.a.claim('ERR999'))

    def test_a_claim_overwritten_by_a_simultaneous_claimer_is_lost(self):
        def rival_writes(dt):                                                                              # the other machine writes between our write and our re-read
            hdr = self.b.data['Blumer'][0]; row = next(r for r in self.b.data['Blumer'][1:] if r[0] == 'ERR004'); row[hdr.index('ClaimedBy')] = 'workstation-9'
        self.a.sleep = rival_writes
        self.assertFalse(self.a.claim('ERR004'))

    def test_runs_to_download(self):
        hdr = self.b.data['Blumer'][0]
        def setrow(run, **kw):
            row = next(r for r in self.b.data['Blumer'][1:] if r[0] == run)
            for k, v in kw.items(): row[hdr.index(k)] = v
        setrow('ERR001', DownloadStatus='failed', FailureKind='permanent'); setrow('ERR002', DownloadStatus='failed', FailureKind='throttle')
        setrow('ERR003', DownloadStatus='downloading', ClaimedBy='workstation-9', ClaimedAt=M.utc_iso(self.clock()))
        self.assertEqual(sorted(self.a.runs_to_download().Run), ['ERR002', 'ERR004'])                    # permanent failure and a fresh foreign claim are left out
        self.assertEqual(sorted(self.a.runs_to_download(retry_failed='all').Run), ['ERR001', 'ERR002', 'ERR004'])
        self.assertEqual(sorted(self.a.runs_to_download(retry_failed=False).Run), ['ERR004'])
        self.clock.sleep(7 * 3600); self.assertIn('ERR003', list(self.a.runs_to_download().Run))             # abandoned
        setrow('ERR004', Skip='TRUE'); self.assertNotIn('ERR004', list(self.a.runs_to_download().Run))

    def test_reconcile_marks_runs_that_reached_the_database_or_the_cloud(self):
        n = self.a.reconcile(in_db={'ERR001'}, on_cloud={'ERR002'}); self.a.flush(force=True)
        self.assertEqual(n, 2)
        self.assertEqual([cellv(self.b, 'Blumer', r, 'DownloadStatus') for r in ('ERR001', 'ERR002', 'ERR003')], ['in_database', 'uploaded', 'pending']); self.assertEqual(cellv(self.b, 'Blumer', 'ERR001', 'InDatabase'), 'TRUE')


class FakeWorksheet:
    def __init__(self, title, rows, cols):
        self.title, self.rows, self.col_count, self.batches = title, [[''] * cols for _ in range(rows)], cols, []
    def _fit(self, r, c):
        while len(self.rows) < r: self.rows.append([''] * self.col_count)
        for row in self.rows:
            while len(row) < c: row.append('')
        self.col_count = max(self.col_count, c)
    def get_all_values(self):
        vals = [[x for x in r] for r in self.rows]
        while vals and not any(vals[-1]): vals.pop()                          # Google trims empty rows
        return [list(r) for r in vals]
    def col_values(self, c): return [r[c - 1] if len(r) >= c else '' for r in self.get_all_values()]
    def row_values(self, r):
        v = self.rows[r - 1]
        while v and v[-1] == '': v = v[:-1]
        return list(v)
    def add_cols(self, n): self.col_count += n; [r.extend([''] * n) for r in self.rows]
    def freeze(self, rows=0): pass
    def batch_update(self, data, value_input_option='USER_ENTERED'):
        assert value_input_option == 'RAW'; self.batches.append(data)
        import re
        for d in data:
            m = re.fullmatch(r'([A-Z]+)(\d+):([A-Z]+)(\d+)', d['range']); assert m, d['range']
            def num(s):
                n = 0
                for ch in s: n = n * 26 + ord(ch) - 64
                return n
            c0, r0, c1, r1 = num(m[1]), int(m[2]), num(m[3]), int(m[4]); assert r0 == r1 and len(d['values']) == 1 and len(d['values'][0]) == c1 - c0 + 1
            self._fit(r0, c1)
            for k, v in enumerate(d['values'][0]): self.rows[r0 - 1][c0 - 1 + k] = v
    def append_rows(self, rows, value_input_option='USER_ENTERED'):
        assert value_input_option == 'RAW'
        n = len(self.get_all_values())
        for k, r in enumerate(rows):
            self._fit(n + k + 1, len(r)); self.rows[n + k][:len(r)] = r


class FakeSpreadsheet:
    def __init__(self): self.ws = {}
    def worksheet(self, t): return self.ws[t]
    def worksheets(self): return list(self.ws.values())
    def add_worksheet(self, title, rows, cols): self.ws[title] = FakeWorksheet(title, rows, cols); return self.ws[title]


class FakeGC:
    def __init__(self): self.sh, self.opened = FakeSpreadsheet(), []
    def open_by_key(self, k): self.opened.append(k); return self.sh


class TestGspreadBackend(unittest.TestCase):
    def test_the_same_flows_on_the_gspread_backend(self):
        gc = FakeGC(); clock = Clock(); b = M.GspreadBackend(gc, '1gNFAd_B8tWxQZRJeEmuhVl1qsqJf1KJtKYpehoEMfrQ')
        m = M.Manifest(b, 'Malinsky', 'utaka-1', now=clock, sleep=clock.sleep, log=lambda *a: None)
        m.import_runinfo(runinfo([f'ERR{i:03d}' for i in range(1, 31)]), in_db={'ERR001'})
        self.assertEqual(gc.opened, ['1gNFAd_B8tWxQZRJeEmuhVl1qsqJf1KJtKYpehoEMfrQ'] * len(gc.opened)); self.assertEqual(len(gc.sh.ws['Malinsky'].get_all_values()), 31)
        self.assertTrue(m.claim('ERR005')); m.update('ERR005', DownloadStatus='uploaded', Source='ENA-aspera', DownloadSeconds=61.2, MBPerSec=12.5); self.assertTrue(m.flush(force=True))
        df = m.load(); r = df[df.Run == 'ERR005'].iloc[0]
        self.assertEqual((r.DownloadStatus, r.Source, r.DownloadSeconds, r.MBPerSec, r.ClaimedBy), ('uploaded', 'ENA-aspera', '61.2', '12.5', 'utaka-1'))
        self.assertEqual(df[df.Run == 'ERR001'].iloc[0].DownloadStatus, 'in_database')
        self.assertEqual(sorted(m.runs_to_download().Run)[:2], ['ERR002', 'ERR003'])

    def test_adjacent_cells_become_one_range_and_big_writes_are_chunked(self):
        gc = FakeGC(); b = M.GspreadBackend(gc, 'x'); m = M.Manifest(b, 'T', 'w', log=lambda *a: None); m.import_runinfo(runinfo(['ERR001']))
        ws = gc.sh.ws['T']; ws.batches.clear()
        b.write_cells('T', [(2, 15, 'a'), (2, 16, 'b'), (2, 17, 'c'), (2, 20, 'd')])
        self.assertEqual(ws.batches[0], [{'range': 'O2:Q2', 'values': [['a', 'b', 'c']]}, {'range': 'T2:T2', 'values': [['d']]}])
        ws.batches.clear(); b.write_cells('T', [(r, 15, 'x') for r in range(2, 502)])
        self.assertEqual([len(x) for x in ws.batches], [200, 200, 100])


if __name__ == '__main__':
    unittest.main()
