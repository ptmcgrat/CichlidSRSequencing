"""The downloader in manifest mode (runs from Google-Sheets tabs, status written back), the management tool, and the file_manager patcher. A memory backend stands in for the sheet; the workers,
the scheduler, the orchestrator and the manifest are the real code."""
import copy, json, os, sys, tempfile, types, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
from test_ena_v2 import FakeFM, fake_gatk, fake_download, rec, csv_to_ena, HERE
from download_reads.sources import ena as E
from download_reads import worker as W, manifest as MF
from download_reads import run as D
from download_reads import manifest_tool as MT
from download_reads.tools import patch_file_manager as PFM


class ManifestFM(FakeFM):
    def __init__(self, tmp, in_db=()):
        super().__init__(tmp, in_db); self.backend = MF.MemoryBackend()
    def addDNAReadRow(self, read_data, sample_data):
        super().addDNAReadRow(read_data, sample_data); self.reads_dt = pd.concat([self.reads_dt, pd.DataFrame([dict(RunID=read_data['RunID'])])], ignore_index=True)
    def downloadManifest(self, tab, who=None, flush_seconds=15.0):
        return MF.Manifest(self.backend, tab, who or 'testhost-1', sleep=lambda s: None, journal=self.localTempDir + f'dm_{tab}.jsonl', flush_s=flush_seconds, log=lambda *a: None)


def sheet(fm, tab='Blumer'):
    return fm.downloadManifest(tab).load().set_index('Run')


class Harness(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp() + '/'
        self.csv = pd.read_csv(os.path.join(HERE, 'MalinskyPaperData.csv'), dtype=str, keep_default_na=False).head(14)
        self.runs = list(self.csv.Run); self.path = self.tmp + 'file.csv'; self.csv.to_csv(self.path, index=False)
        self.fm = ManifestFM(self.tmp, in_db=self.runs[0:2])
        proj = dict(zip(self.csv.Run, self.csv.BioProject)); self.fm.cloud[self.fm.localReadsDir + proj[self.runs[2]] + '/'] = {self.runs[2] + '.unmapped_marked_adapters.bam'}
        MT.cmd_import(self.fm, 'Blumer', self.path, log=lambda *a: None)
        self.ena = csv_to_ena(self.csv); r = self.runs
        no = self.ena.run_accession == r[3]
        self.ena.loc[no, ['fastq_ftp', 'fastq_aspera', 'fastq_md5', 'fastq_bytes']] = ''; self.ena.loc[no, 'submitted_ftp'] = 'ftp.sra.ebi.ac.uk/x/a.cram'; self.ena.loc[no, 'submitted_format'] = 'CRAM'
        self.throttled = {'done': False}; self.logs = []
        self.fetch = lambda acc: self.ena[(self.ena.study_accession == acc) | (self.ena.run_accession == acc)].reset_index(drop=True)

    def download(self, sel, dest, **kw):
        if sel['r1']['name'].startswith(self.runs[5]) and not self.throttled['done']:
            self.throttled['done'] = True; raise E.EnaError('Session Stop: Too many connections', 'throttle')
        return fake_download()(sel, dest)

    def launch(self, kind, spec_path, logf):
        spec = json.load(open(spec_path))
        st = W.fetch_stage(spec, download=self.download, log=lambda *a: None) if kind == 'fetch' else W.convert_stage(spec, self.fm, run=fake_gatk(fail_for=[self.runs[4]]), log=lambda *a: None)
        return types.SimpleNamespace(poll=lambda: 0 if st['state'] in ('done', 'fetched') else W.EXIT[st['kind']])

    def go(self, *extra, tab='Blumer'):
        clock = [0.0]
        args = D.parse_args(['--manifest', tab, '--stagger', '0', '--backoff', '0', '--workers', '2', '--dl-workers', '2', '--ascp', '/no/ascp', '--source', 'ena'] + list(extra))
        rc = D.run(args, self.fm, fetch=self.fetch, launch=self.launch, sleep=lambda s: clock.__setitem__(0, clock[0] + s), now=lambda: clock[0], free_fn=lambda p: 1e12, log=self.logs.append)
        return rc


class TestImportTool(Harness):
    def test_import_statuses_and_reimport(self):
        s = sheet(self.fm)
        self.assertEqual(len(s), 14); self.assertEqual(list(s.loc[self.runs[:3], 'DownloadStatus']), ['in_database', 'in_database', 'uploaded']); self.assertTrue((s.loc[self.runs[3:], 'DownloadStatus'] == 'pending').all())
        self.assertEqual(MT.cmd_import(self.fm, 'Blumer', self.path, log=lambda *a: None), dict(added=0, existing=14))
        self.assertEqual(MT.cmd_import(self.fm, 'Other', self.path, log=lambda *a: None)['added'], 14)                    # a second tab with the same runs is its own tab

    def test_status_skips_tabs_that_are_not_manifest_tabs(self):                      # the empty default 'Sheet1' of a new spreadsheet crashed `status`
        self.fm.backend.data['Sheet1'] = [[]]
        out = []; r = MT.cmd_status(self.fm, None, log=out.append)
        self.assertEqual(list(r), ['Blumer']); self.assertTrue(any('Sheet1: not a manifest tab' in x for x in out))
        self.assertTrue(self.fm.downloadManifest('Sheet1').load().empty)                                         # loading it is harmless too
        self.fm.backend.data['Blank'] = []; self.assertTrue(self.fm.downloadManifest('Blank').load().empty)

    def test_status_and_release(self):
        mf = self.fm.downloadManifest('Blumer'); mf.update(self.runs[5], DownloadStatus='downloading', ClaimedBy='old-1', ClaimedAt=MF.utc_iso(1_000_000_000)); mf.flush(force=True)
        out = []; r = MT.cmd_status(self.fm, ['Blumer'], log=out.append)
        self.assertEqual(r['Blumer']['stale'], [self.runs[5]]); self.assertTrue(any('STALE' in x for x in out)); self.assertTrue(any('in_database 2' in x for x in out))
        self.assertEqual(MT.cmd_release(self.fm, 'Blumer', log=lambda *a: None), 1); self.assertEqual(sheet(self.fm).loc[self.runs[5], 'DownloadStatus'], 'pending')


class TestManifestMode(Harness):
    def test_full_cycle_statuses_sources_attempts_and_claims(self):
        r = self.runs
        other = self.fm.downloadManifest('Blumer', 'other-1'); self.assertTrue(other.claim(r[6]))                         # another machine already has run 6
        mf = self.fm.downloadManifest('Blumer'); mf.update(r[7], Skip=True); mf.flush(force=True)                          # a person marked run 7 to skip
        rc = self.go()
        s = sheet(self.fm)
        self.assertEqual(rc, 1)
        self.assertEqual(list(s.loc[r[:3], 'DownloadStatus']), ['in_database', 'in_database', 'uploaded'])               # untouched
        self.assertEqual((s.loc[r[3], 'DownloadStatus'], s.loc[r[3], 'FailureKind']), ('skipped', 'permanent')); self.assertTrue(s.loc[r[3], 'LastError'].startswith('no_fastq'))
        self.assertEqual((s.loc[r[4], 'DownloadStatus'], s.loc[r[4], 'FailureKind'], s.loc[r[4], 'ClaimedBy']), ('failed', 'permanent', '')); self.assertIn('FastqToSam', s.loc[r[4], 'LastError'])
        row = s.loc[r[5]]
        self.assertEqual((row.DownloadStatus, row.InDatabase, row.Attempts, row.Source, row.Verified, row.UBamBytes), ('in_database', 'TRUE', '2', 'ENA-aspera', 'TRUE', '123'))
        self.assertTrue(float(row.ConvertSeconds) >= 0 and float(row.UploadSeconds) >= 0)                                   # stage timings reach the sheet
        self.assertTrue(row.DownloadedOn.endswith('Z') and row.UBamPath.endswith(f'{r[5]}.unmapped_marked_adapters.bam') and row.Machine and row.ClaimedBy == '' and float(row.MBPerSec) > 0)
        self.assertEqual((s.loc[r[6], 'DownloadStatus'], s.loc[r[6], 'ClaimedBy']), ('claimed', 'other-1'))                # not touched: someone else's
        self.assertEqual((s.loc[r[7], 'DownloadStatus'], s.loc[r[7], 'Skip']), ('pending', 'TRUE'))
        self.assertTrue((s.loc[r[8:], 'DownloadStatus'] == 'in_database').all())
        added = {x[0]['RunID'] for x in self.fm.rows}
        self.assertEqual(added, set(r[5:6] + r[8:]))                                                                       # database rows for exactly the successes

    def test_a_second_invocation_does_nothing_unless_failures_are_retried(self):
        self.go(); n = len(self.fm.rows); self.assertEqual(sheet(self.fm).loc[self.runs[4], 'Attempts'], '1'); self.logs.clear()
        self.assertEqual(self.go(), 0)
        self.assertTrue(any('nothing to download' in str(m) for m in self.logs)); self.assertEqual(len(self.fm.rows), n)
        self.go('--retry-failed', 'all')                                                                                    # the permanent failure is tried again only on request
        self.assertEqual(sheet(self.fm).loc[self.runs[4], 'Attempts'], '2'); self.assertEqual(sheet(self.fm).loc[self.runs[4], 'DownloadStatus'], 'failed')      # one attempt per invocation, added up

    def test_add_pending_with_the_manifest_turns_uploaded_rows_into_in_database(self):
        self.go('--no-db'); r = self.runs
        self.assertEqual(sheet(self.fm).loc[r[8], 'DownloadStatus'], 'uploaded')
        pend = self.fm.localTempDir + 'pending_db_rows.jsonl'
        clock = [0.0]
        D.run(D.parse_args(['--manifest', 'Blumer', '--add-pending', pend]), self.fm, fetch=self.fetch, launch=self.launch, log=lambda *a: None)
        s = sheet(self.fm)
        self.assertEqual((s.loc[r[8], 'DownloadStatus'], s.loc[r[8], 'InDatabase']), ('in_database', 'TRUE')); self.assertIn(r[8], {x[0]['RunID'] for x in self.fm.rows})

    def test_reconcile_command(self):
        self.fm.reads_dt = pd.concat([self.fm.reads_dt, pd.DataFrame([dict(RunID=self.runs[2])])], ignore_index=True)           # the uBAM on Dropbox got its database row meanwhile
        self.assertEqual(sheet(self.fm).loc[self.runs[2], 'DownloadStatus'], 'uploaded')
        r = MT.cmd_reconcile(self.fm, ['Blumer'], log=lambda *a: None); self.assertEqual(r['Blumer'], 1)
        self.assertEqual(sheet(self.fm).loc[self.runs[2], 'DownloadStatus'], 'in_database'); self.assertEqual(MT.cmd_reconcile(self.fm, ['Blumer'], log=lambda *a: None)['Blumer'], 0)

    def test_dry_run_reads_but_never_writes_the_sheet(self):
        before = copy.deepcopy(self.fm.backend.data); self.go('--dry-run')
        self.assertEqual(self.fm.backend.data, before)

    def test_no_db_marks_runs_uploaded_not_in_database(self):
        self.go('--no-db'); s = sheet(self.fm)
        self.assertEqual((s.loc[self.runs[8], 'DownloadStatus'], s.loc[self.runs[8], 'InDatabase']), ('uploaded', 'FALSE')); self.assertEqual(self.fm.rows, [])

    def test_an_unreachable_sheet_does_not_stop_the_downloads_and_the_journal_catches_up(self):
        calls = {'n': 0}
        def flaky(name):
            if name == 'write_cells':
                calls['n'] += 1
                if calls['n'] <= 6: raise ConnectionError('google unreachable')
        self.fm.backend.fail = flaky
        self.go(); self.fm.backend.fail = None
        self.assertEqual({x[0]['RunID'] for x in self.fm.rows}, set(self.runs[5:3:-1][:1] + self.runs[6:]))                      # every download happened although the sheet was down
        mf = self.fm.downloadManifest('Blumer'); mf.flush(force=True)                                                         # a later flush replays whatever was journalled
        s = sheet(self.fm)
        self.assertEqual(s.loc[self.runs[8], 'DownloadStatus'], 'in_database'); self.assertEqual(s.loc[self.runs[4], 'DownloadStatus'], 'failed')


class TestFileManagerPatch(unittest.TestCase):
    SRC = ("import platform, os\n\nclass FileManager():\n\tdef __init__(self):\n\t\tself.s_ID = 'S'\n\t\tself.a_ID = 'A' # Alignment Database\n\n\t\tself._createMasterDirs()\n\n"
           "\tdef _authenticateGS(self):\n\t\tpass\n")

    def test_two_insertions_idempotent_and_compilable(self):
        new, done = PFM.patch(self.SRC); compile(new, 'fm', 'exec')
        self.assertEqual(done, ['dl_ID', 'methods']); self.assertIn("self.dl_ID = '1gNFAd_B8tWxQZRJeEmuhVl1qsqJf1KJtKYpehoEMfrQ'", new)
        self.assertTrue(new.index('self.a_ID') < new.index('self.dl_ID') < new.index('self._createMasterDirs()'))
        self.assertTrue(new.index('def downloadManifest') < new.index('def _authenticateGS')); self.assertEqual(PFM.patch(new), (new, []))
        ns = {}; exec(new.replace('class FileManager', 'class FileManager'), ns); self.assertTrue(all(hasattr(ns['FileManager'], m) for m in ('downloadManifest', 'downloadTabs', 'readDownloadTab')))

    def test_a_missing_anchor_is_an_error_not_a_silent_skip(self):
        with self.assertRaises(ValueError):
            PFM.patch('class FileManager():\n\tdef x(self):\n\t\tpass\n')

    OLD = ("import platform, os\n\nclass FileManager():\n\tdef __init__(self):\n\t\tself.s_ID = 'S'\n\t\tself.a_ID = 'A' # Alignment Database\n\t\tself.dl_ID = 'X'\n\n"
           "\t# ---- Download database\n\tdef downloadManifest(self, tab, who = None, flush_seconds = 15.0):\n\t\ttry:\n\t\t\tfrom helper_modules import download_manifest as dm\n\t\texcept ImportError:\n\t\t\timport download_manifest as dm\n"
           "\t\twho = who or 'x'\n\t\treturn None\n\n\tdef _authenticateGS(self):\n\t\tpass\n")

    def test_a_file_patched_by_the_earlier_version_is_upgraded_to_the_new_import_chain(self):
        new, done = PFM.patch(self.OLD); compile(new, 'fm', 'exec')
        self.assertEqual(done, ['methods upgraded']); self.assertIn('from download_reads import manifest as dm', new); self.assertIn('from helper_modules import download_manifest as dm', new)      # new first, the old place as the fallback
        self.assertEqual(PFM.patch(new), (new, []))

    def test_the_methods_delegate_to_the_manifest_with_a_journal_and_the_right_spreadsheet(self):
        new, _ = PFM.patch(self.SRC.replace("pass\n", "pass\n"))
        ns = {'__file__': '/tmp/repo/cichlid_sr_sequencing/helper_modules/file_manager.py'}; exec(new, ns); fmc = ns['FileManager']
        o = fmc.__new__(fmc); o.dl_ID = 'SHEET'; o.gc = object(); o.localTempDir = '/tmp/x/'
        m = o.downloadManifest('My Tab', who='w-1')
        self.assertEqual((m.tab, m.who, m.b.id, m.journal), ('My Tab', 'w-1', 'SHEET', '/tmp/x/download_manifest_My_Tab.journal.jsonl'))


if __name__ == '__main__':
    unittest.main()
