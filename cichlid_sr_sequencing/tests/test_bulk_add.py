"""--add-pending --bulk: the same rows as addDNAReadRow run by run, with one read and one write of the two sheets, a backup, and a check afterwards. The sheets are a fake that counts calls and
enforces the grid limit like Google (set_with_dataframe does not enlarge the grid)."""
import json, os, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
import downloadShortReadData_v2 as D

S_COLS = ['SampleID', 'Sex', 'Species', 'DoB', 'BroodID', 'Parents', 'Ecogroup', 'LabReared', 'Inversion10', 'gsdf number']
R_COLS = ['SampleID', 'ProjectID', 'RunID', 'ReadLength', 'TotalBases', 'Instrument', 'LibraryID', 'LibraryLayout', 'LibrarySource', 'Platform', 'FileLocations', 'FileSize']


def old_samples(n):
    return pd.DataFrame([dict(SampleID=f'OLD{i}', Sex='f', Species='MC', DoB='', BroodID='', Parents='', Ecogroup='x', LabReared='', Inversion10='0/0', **{'gsdf number': 2}) for i in range(n)], columns=S_COLS)


def old_reads(n):
    return pd.DataFrame([dict(SampleID=f'OLD{i}', ProjectID='PRJ0', RunID=f'OLDRUN{i}', ReadLength=250, TotalBases=1, Instrument='HiSeq', LibraryID='l', LibraryLayout='PAIRED', LibrarySource='GENOMIC', Platform='ILLUMINA',
                         FileLocations=f'PRJ0/OLDRUN{i}.unmapped_marked_adapters.bam', FileSize=5) for i in range(n)], columns=R_COLS)


def line(run, sample, project='PRJ', size=''):
    rd = dict(SampleID=sample, ProjectID=project, RunID=run, ReadLength=250, TotalBases=1000, Instrument='HiSeq', LibraryID=sample + '_lib', LibraryLayout='PAIRED', LibrarySource='GENOMIC', Platform='ILLUMINA',
              FileLocations=f'{project}/{run}.unmapped_marked_adapters.bam', FileSize=size)
    sd = dict(SampleID=sample, Sex='', Species='x', DoB='', BroodID='', Parents='', Ecogroup='', LabReared='', Inversion10='')
    return json.dumps(dict(read_data=rd, sample_data=sd))


class WS:
    def __init__(self, fm, name): self.fm, self.name = fm, name
    row_count = property(lambda self: self.fm.grid[self.name])
    def resize(self, rows): self.fm.grid[self.name] = rows


class SheetFM:
    def __init__(self, n_old=5, grid=2000, fail_write=False, corrupt_write=False):
        self.localTempDir = tempfile.mkdtemp() + '/'; self.localReadsDir = self.localTempDir + 'Reads/'
        self.remote = dict(SampleDatabase=old_samples(n_old), DNAReadsDatabase=old_reads(n_old)); self.grid = dict(SampleDatabase=grid, DNAReadsDatabase=grid)
        self.cloud, self.s_ID, self.gc = {}, 'S', self
        self.reads_n = self.writes_n = self.size_calls = self.cloud_calls = 0; self.fail_write, self.corrupt_write = fail_write, corrupt_write
    def open_by_key(self, k): return self
    def worksheet(self, name): return WS(self, name)
    def readSampleDatabase(self):
        self.reads_n += 1; self.sample_dt, self.reads_dt = self.remote['SampleDatabase'].copy(), self.remote['DNAReadsDatabase'].copy()
    def setSampleDatabase(self):
        if self.fail_write:
            return False
        for name, df in (('SampleDatabase', self.sample_dt), ('DNAReadsDatabase', self.reads_dt)):
            if len(df) + 1 > self.grid[name]:
                return False                                                          # Google: range exceeds grid limits (the real method prints it three times and returns False)
        self.remote['SampleDatabase'] = self.sample_dt.copy()
        self.remote['DNAReadsDatabase'] = self.reads_dt.iloc[:-1].copy() if self.corrupt_write else self.reads_dt.copy()
        self.writes_n += 1; return True
    def returnCloudFiles(self, d): self.cloud_calls += 1; return sorted(self.cloud.get(d, set()))
    def returnFileSize(self, p): self.size_calls += 1; return 999
    def addDNAReadRow(self, dna_dict, sample_dict):                                   # the real method, line by line
        self.readSampleDatabase()
        assert set(dna_dict.keys()) == set(R_COLS) and dna_dict['RunID'] not in self.reads_dt.RunID.to_list()
        self.reads_dt = pd.concat([self.reads_dt, pd.DataFrame([dna_dict])], ignore_index=True)
        if sample_dict['SampleID'] not in self.sample_dt.SampleID.to_list():
            self.sample_dt = pd.concat([self.sample_dt, pd.DataFrame([sample_dict])], ignore_index=True)
        self.setSampleDatabase()
    def put(self, *runs, project='PRJ'):
        self.cloud.setdefault(self.localReadsDir + project + '/', set()).update(f'{r}.unmapped_marked_adapters.bam' for r in runs)


def pending_file(lines):
    p = tempfile.mkdtemp() + '/pending.jsonl'; open(p, 'w').write('\n'.join(lines) + '\n'); return p


LINES = [line(f'RUN{i}', f'S{i // 2}') for i in range(12)]                          # 12 runs of 6 samples (two runs each)


class T(unittest.TestCase):
    def test_bulk_gives_the_same_sheets_as_run_by_run_with_two_reads_and_one_write(self):
        path = pending_file(LINES); logs = []
        a = SheetFM(); a.put(*[f'RUN{i}' for i in range(12)]); added, skipped, missing = D.add_pending_bulk(a, path, log=logs.append)
        b = SheetFM(); b.put(*[f'RUN{i}' for i in range(12)]); b.readSampleDatabase(); D.add_pending(b, path)
        for tab in ('SampleDatabase', 'DNAReadsDatabase'):
            pd.testing.assert_frame_equal(a.remote[tab], b.remote[tab])
        self.assertEqual((added, skipped, missing, len(a.remote['DNAReadsDatabase']), len(a.remote['SampleDatabase'])), (12, 0, [], 17, 11))      # 6 new samples, not 12
        self.assertEqual((a.reads_n, a.writes_n), (2, 1)); self.assertEqual((b.reads_n, b.writes_n), (13, 12))                                 # the point of the bulk mode
        self.assertTrue(any('verified: DNAReadsDatabase 5 -> 17' in str(l) for l in logs))

    def test_runs_already_in_the_database_duplicates_and_missing_ubams_are_not_added(self):
        fm = SheetFM(); fm.put('RUN0', 'RUN1', 'RUN2', 'OLDRUN0')                                                       # RUN3 has no uBAM on the cloud
        path = pending_file([line('OLDRUN0', 'OLD0'), line('RUN0', 'S0'), line('RUN0', 'S0'), line('RUN1', 'S0'), line('RUN2', 'S1'), line('RUN3', 'S1')])
        added, skipped, missing = D.add_pending_bulk(fm, path, log=lambda *a: None)
        self.assertEqual((added, skipped, missing), (3, 2, ['RUN3']))
        self.assertEqual(sorted(fm.remote['DNAReadsDatabase'].RunID.tolist()[5:]), ['RUN0', 'RUN1', 'RUN2']); self.assertEqual(fm.remote['SampleDatabase'].SampleID.tolist()[5:], ['S0', 'S1'])

    def test_a_sample_that_is_already_in_the_sheet_is_not_added_again_and_blank_sizes_are_filled(self):
        fm = SheetFM(); fm.put('RUN0'); D.add_pending_bulk(fm, pending_file([line('RUN0', 'OLD1')]), log=lambda *a: None)
        self.assertEqual(len(fm.remote['SampleDatabase']), 5); self.assertEqual(fm.remote['DNAReadsDatabase'].iloc[-1].SampleID, 'OLD1')
        self.assertEqual((fm.remote['DNAReadsDatabase'].iloc[-1].FileSize, fm.size_calls), (999, 1))

    def test_both_sheets_are_backed_up_first_and_old_rows_are_untouched(self):
        fm = SheetFM(); fm.put(*[f'RUN{i}' for i in range(12)]); before = fm.remote['DNAReadsDatabase'].copy()
        D.add_pending_bulk(fm, pending_file(LINES), log=lambda *a: None)
        d = [x for x in os.listdir(fm.localTempDir) if x.startswith('database_backup_')]; self.assertEqual(len(d), 1)
        self.assertEqual(len(pd.read_csv(fm.localTempDir + d[0] + '/DNAReadsDatabase.csv')), 5); self.assertEqual(len(pd.read_csv(fm.localTempDir + d[0] + '/SampleDatabase.csv')), 5)
        pd.testing.assert_frame_equal(fm.remote['DNAReadsDatabase'].iloc[:5].reset_index(drop=True), before.reset_index(drop=True))

    def test_a_dry_run_reads_but_never_writes_or_backs_up(self):
        fm = SheetFM(); fm.put(*[f'RUN{i}' for i in range(12)])
        added, _, _ = D.add_pending_bulk(fm, pending_file(LINES), log=lambda *a: None, dry_run=True)
        self.assertEqual((added, fm.reads_n, fm.writes_n), (12, 1, 0)); self.assertFalse([x for x in os.listdir(fm.localTempDir) if x.startswith('database_backup_')])
        self.assertEqual(len(fm.remote['DNAReadsDatabase']), 5)

    def test_a_sheet_without_enough_rows_is_enlarged_before_the_write(self):
        fm = SheetFM(grid=8); fm.put(*[f'RUN{i}' for i in range(12)]); logs = []
        D.add_pending_bulk(fm, pending_file(LINES), log=logs.append)
        self.assertEqual(len(fm.remote['DNAReadsDatabase']), 17); self.assertTrue(any('grid enlarged from 8 to' in str(l) for l in logs)); self.assertGreaterEqual(fm.grid['DNAReadsDatabase'], 18)

    def test_a_failed_write_and_a_wrong_read_back_stop_with_the_backup_named(self):
        fm = SheetFM(fail_write=True); fm.put(*[f'RUN{i}' for i in range(12)])
        with self.assertRaises(SystemExit) as cm:
            D.add_pending_bulk(fm, pending_file(LINES), log=lambda *a: None)
        self.assertIn('restore from', str(cm.exception)); self.assertIn('database_backup_', str(cm.exception))
        fm = SheetFM(corrupt_write=True); fm.put(*[f'RUN{i}' for i in range(12)])
        with self.assertRaises(SystemExit) as cm:
            D.add_pending_bulk(fm, pending_file(LINES), log=lambda *a: None)
        self.assertIn('do not hold what was written', str(cm.exception))

    def test_a_malformed_line_stops_everything_before_any_write(self):
        fm = SheetFM(); fm.put('RUN0', 'RUN1')
        bad = json.loads(line('RUN1', 'S1')); del bad['read_data']['Platform']
        with self.assertRaises(SystemExit) as cm:
            D.add_pending_bulk(fm, pending_file([line('RUN0', 'S0'), json.dumps(bad), 'not json']), log=lambda *a: None)
        self.assertIn('nothing was written', str(cm.exception)); self.assertIn('RUN1', str(cm.exception))
        self.assertEqual((fm.writes_n, [x for x in os.listdir(fm.localTempDir) if x.startswith('database_backup_')]), (0, []))

    def test_one_listing_per_project_folder_gives_existence_and_sizes_without_any_per_run_call(self):
        fm = SheetFM(); calls = []
        def listing(folder):
            calls.append(folder); return {f'RUN{i}.unmapped_marked_adapters.bam': 1000 + i for i in range(12) if f'/{"PRJA" if i < 6 else "PRJB"}/' in folder}
        lines = [line(f'RUN{i}', f'S{i // 2}', project='PRJA' if i < 6 else 'PRJB') for i in range(12)] + [line('RUN99', 'S99', project='PRJA')]
        added, skipped, missing = D.add_pending_bulk(fm, pending_file(lines), log=lambda *a: None, listing=listing)
        self.assertEqual((added, missing), (12, ['RUN99'])); self.assertEqual(len(calls), 2)                                  # 13 runs, 2 folders, 2 listings
        self.assertEqual((fm.size_calls, fm.cloud_calls), (0, 0))                                                            # no rclone size / lsf per run
        r = fm.remote['DNAReadsDatabase'].set_index('RunID'); self.assertEqual((r.loc['RUN0', 'FileSize'], r.loc['RUN7', 'FileSize']), (1000, 1007))

    def test_a_dry_run_asks_for_no_sizes(self):
        fm = SheetFM(); fm.put(*[f'RUN{i}' for i in range(12)])
        D.add_pending_bulk(fm, pending_file(LINES), log=lambda *a: None, dry_run=True)
        self.assertEqual(fm.size_calls, 0)

    def test_the_rclone_listing_maps_the_path_and_distinguishes_missing_from_failed(self):
        class R:
            def __init__(self, rc, out='', err=''): self.returncode, self.stdout, self.stderr = rc, out, err
        class FM:
            localMasterDir, cloudMasterDir = '/home/u/Temp/CSD/', 'ptm_dropbox:/COS/CSD/'
        calls = []
        def run(cmd, capture_output=True, text=True):
            calls.append(cmd); return R(0, json.dumps([dict(Name='a.bam', Size=5, Path='a.bam'), dict(Name='b.bam', Size=7, Path='b.bam')]))
        self.assertEqual(D._cloud_listing(FM, '/home/u/Temp/CSD/Reads/PRJ1/', run=run), {'a.bam': 5, 'b.bam': 7})
        self.assertEqual(calls[0], ['rclone', 'lsjson', '--files-only', 'ptm_dropbox:/COS/CSD/Reads/PRJ1/'])
        self.assertEqual(D._cloud_listing(FM, '/home/u/Temp/CSD/Reads/X/', run=lambda *a, **k: R(3)), {})                 # directory not found = empty
        with self.assertRaises(RuntimeError) as cm:
            D._cloud_listing(FM, '/home/u/Temp/CSD/Reads/X/', run=lambda *a, **k: R(1, '', 'network unreachable'))
        self.assertIn('network unreachable', str(cm.exception))                                                           # a failed listing is never "no files"

    def test_the_default_listing_is_rclone_when_the_file_manager_has_cloud_paths(self):
        fm = SheetFM(); self.assertEqual(D._default_listing(fm)('/x/'), {})                                               # the fake has no cloud paths: returnCloudFiles, sizes later
        fm.cloudMasterDir, fm.localMasterDir = 'remote:/m/', fm.localTempDir
        import unittest.mock as um
        with um.patch.object(D.subprocess, 'run', lambda cmd, **k: type('R', (), dict(returncode=0, stdout=json.dumps([dict(Name='z.bam', Size=3)]), stderr=''))()):
            self.assertEqual(D._default_listing(fm)(fm.localTempDir + 'Reads/P/'), {'z.bam': 3})

    def test_the_command_line_options(self):
        fm = SheetFM(); fm.put(*[f'RUN{i}' for i in range(12)]); path = pending_file(LINES)
        D.run(D.parse_args(['--add-pending', path, '--dry-run']), fm, log=lambda *a: None)                                  # --dry-run never writes, even without --bulk
        self.assertEqual((fm.writes_n, len(fm.remote['DNAReadsDatabase'])), (0, 5))
        D.run(D.parse_args(['--add-pending', path, '--bulk']), fm, log=lambda *a: None)
        self.assertEqual((fm.writes_n, len(fm.remote['DNAReadsDatabase'])), (1, 17))


if __name__ == '__main__':
    unittest.main()
