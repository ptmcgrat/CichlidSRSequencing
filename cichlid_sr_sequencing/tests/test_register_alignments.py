"""unit_scripts/register_alignments.py: the row must be what alignment_worker builds, written once, only for samples without a row, with the BAMs of this machine."""
import csv, os, subprocess, sys, tempfile, types, unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'unit_scripts'))
import pandas as pd
import register_alignments as R

COLS = ['SampleID', 'GenomeVersion', 'RunIDs', 'Coverage', 'TotalReads', 'UnmappedReads', 'DiscordantReads', 'InversionReads', 'DuplicationReads', 'ClippedReads', 'ChimericReads', 'minimap2_version', 'gatk_version',
        'pysam_version', 'BamSize']
VERS = {'minimap2_version': '2.28', 'gatk_version': '4.5', 'pysam_version': '0.22'}
REF = 1_000_000


class WS:
    def __init__(self, fm): self.fm = fm
    row_count = property(lambda self: self.fm.grid)
    def resize(self, rows): self.fm.grid = rows


class FM:
    def __init__(self, registered=(), grid=2000, fail_write=False, corrupt=False):
        tmp = tempfile.mkdtemp() + '/'; self.localTempDir = tmp + 'Temp/'; self.localBamRefDir = tmp + 'Bam/'; os.makedirs(self.localTempDir); os.makedirs(self.localBamRefDir)
        self.a_ID, self.gc, self.grid, self.fail_write, self.corrupt = 'A', self, grid, fail_write, corrupt
        self.reads_n = self.writes_n = 0; self.remote = pd.DataFrame([self.old_row(s) for s in registered], columns=COLS); self.reads = {}
        for s in registered: self.make(s); self.add_reads(s, ['R_' + s])
    def old_row(self, s):
        return dict(SampleID=s, GenomeVersion='G1', RunIDs='R_' + s, Coverage=1000 * 125 / REF, TotalReads=1000, UnmappedReads=0.1, DiscordantReads=0.05, InversionReads=0.0, DuplicationReads=0.0, ClippedReads=0.2, ChimericReads=0.01,
                    BamSize=5, **VERS)
    def make(self, s, **counts):
        c = dict(all=1000, unmapped=100, discordant=50, inversion=0, duplication=0, clipped=200, chimeric=10); c.update(counts)
        os.makedirs(self.localBamRefDir + s, exist_ok=True)
        for k, e in zip(R.KEYS, R.ENDINGS):
            open(f'{self.localBamRefDir}{s}/{s}{e}', 'w').write(str(c[k]))
    def add_reads(self, s, runs, length=250, listed=True): self.reads[s] = (runs, length); self.unlisted = getattr(self, 'unlisted', set()) | (set() if listed else {s})
    def open_by_key(self, k): return self
    def worksheet(self, name): return WS(self)
    def setGenome(self, g): self.genome_version = g
    def readAlignmentDatabase(self): self.reads_n += 1; self.alignment_dt = self.remote.copy()
    def readSampleDatabase(self):
        rows = [(s, r, l) for s, (runs, l) in self.reads.items() for r in runs]
        self.reads_dt = pd.DataFrame(rows, columns=['SampleID', 'RunID', 'ReadLength']); self.merged_dt = self.reads_dt[~self.reads_dt.SampleID.isin(getattr(self, 'unlisted', set()))][['SampleID', 'ReadLength']].copy()
    def setAlignmentDatabase(self):
        if self.fail_write or len(self.alignment_dt) + 1 > self.grid:
            return False
        self.remote = self.alignment_dt.iloc[:-1].copy() if self.corrupt else self.alignment_dt.copy(); self.writes_n += 1; return True


def counter(path):
    return int(open(path).read())


def run(fm, argv, **kw):
    out = []
    res = R.main(['--genome', 'G1'] + argv, fm=fm, counter=counter, versions=VERS, reference_size=REF, out=out.append, **kw)
    return res, out


def audit_file(rows):
    p = tempfile.mkdtemp() + '/audit.tsv'
    with open(p, 'w') as f:
        f.write('SampleID\tclass\tdetails\tEcogroup\tfiles_dropbox\tfiles_local\n')
        for s, c in rows: f.write(f'{s}\t{c}\t\tMbuna\t17\t17\n')
    return p


class T(unittest.TestCase):
    def test_the_row_is_what_alignment_worker_builds(self):
        fm = FM(); fm.make('S1', all=2000, unmapped=200, discordant=100, inversion=20, duplication=40, clipped=400, chimeric=10); fm.add_reads('S1', ['RUN1', 'RUN2'], length=250)
        res, out = run(fm, ['--samples', 'S1'])
        r = res[0]
        self.assertEqual(set(r), set(COLS))
        self.assertEqual((r['SampleID'], r['GenomeVersion'], r['RunIDs'], r['TotalReads']), ('S1', 'G1', 'RUN1,,RUN2', 2000))                                       # ',,' like the original
        self.assertAlmostEqual(r['Coverage'], 2000 * 125 / REF)                                                                                                 # ReadLength / 2 = 125
        self.assertEqual((r['UnmappedReads'], r['DiscordantReads'], r['InversionReads'], r['DuplicationReads'], r['ClippedReads'], r['ChimericReads']), (0.1, 0.05, 0.01, 0.02, 0.2, 0.005))
        self.assertEqual((r['BamSize'], r['minimap2_version']), (len(str(2000)), '2.28'))

    def test_only_samples_without_a_row_and_with_local_bams_are_registered_with_two_reads_and_one_write(self):
        fm = FM(registered=['DONE1']); fm.add_reads('DONE1', ['R_DONE1'])
        for s in ('N1', 'N2', 'N3'): fm.make(s); fm.add_reads(s, ['R_' + s])
        fm.add_reads('NOLOCAL', ['R_NOLOCAL']); fm.add_reads('NOREADS_ROW', [])                                                                              # NOLOCAL: no BAMs here
        fm.make('NOREADS')                                                                                                                                     # BAMs but not in the DNAReadsDatabase
        path = audit_file([('N1', 'needs_row'), ('N2', 'needs_row'), ('N3', 'needs_row'), ('DONE1', 'needs_row'), ('NOLOCAL', 'needs_row'), ('NOREADS', 'needs_row'), ('X', 'not_started')])
        res, out = run(fm, ['--from-audit', path]); text = '\n'.join(out)
        self.assertEqual([r['SampleID'] for r in res], ['N1', 'N2', 'N3']); self.assertEqual(fm.remote.SampleID.tolist(), ['DONE1', 'N1', 'N2', 'N3'])
        self.assertEqual((fm.reads_n, fm.writes_n), (2, 1))                                                                                                     # the read before and the read-back
        self.assertIn('5 samples to register (1 already have a row', text); self.assertIn('NOLOCAL: missing locally', text); self.assertIn('NOREADS is not in the DNAReadsDatabase', text)
        self.assertIn('verified: AlignmentDatabase G1 1 -> 4 rows', text)

    def test_samples_that_are_not_in_the_sample_database_are_skipped_like_alignment_worker_unless_allowed(self):
        fm = FM(); fm.make('UNL'); fm.add_reads('UNL', ['RUNX'], length=300, listed=False)
        res, out = run(fm, ['--samples', 'UNL', '--dry-run'])
        self.assertEqual(res, []); self.assertIn('no row in the SampleDatabase', '\n'.join(out)); self.assertIn('--allow-unlisted', '\n'.join(out))
        res, out = run(fm, ['--samples', 'UNL', '--allow-unlisted'])
        self.assertEqual((res[0]['SampleID'], res[0]['RunIDs']), ('UNL', 'RUNX')); self.assertAlmostEqual(res[0]['Coverage'], 1000 * 150 / REF)                      # ReadLength 300 / 2
        self.assertEqual(fm.remote.SampleID.tolist(), ['UNL'])

    def test_a_dry_run_and_the_limit(self):
        fm = FM()
        for s in ('N1', 'N2', 'N3'): fm.make(s); fm.add_reads(s, ['R_' + s])
        res, out = run(fm, ['--samples', 'N1', 'N2', 'N3', '--limit', '2', '--dry-run'])
        self.assertEqual(([r['SampleID'] for r in res], fm.writes_n, len(fm.remote)), (['N1', 'N2'], 0, 0)); self.assertFalse([x for x in os.listdir(fm.localTempDir) if x.startswith('database_backup')])

    def test_check_existing_validates_the_counting_and_the_coverage_against_the_sheet(self):
        fm = FM(registered=['A1', 'A2']); res, out = run(fm, ['--check-existing', '2']); text = '\n'.join(out)
        self.assertTrue(res); self.assertIn('the counting MATCHES the sheet', text); self.assertIn('9 values compared with the sheet', text)                        # 7 counts-derived columns, coverage and RunIDs
        fm = FM(registered=['A1']); fm.make('A1', all=1001); res, out = run(fm, ['--check-existing', '1'])                                                       # the BAM now holds another count than the sheet
        self.assertFalse(res); self.assertIn('DIFFERENT', '\n'.join(out))
        fm = FM(registered=['A1']); fm.remote.loc[0, 'Coverage'] = 9.9; res, out = run(fm, ['--check-existing', '1'])                                              # same counts, but the coverage in the sheet is not what the formula gives
        self.assertFalse(res); self.assertIn('Coverage: sheet 9.9, computed 0.125', '\n'.join(out))
        fm = FM(registered=['A1']); fm.add_reads('A1', ['OTHER_RUN']); res, out = run(fm, ['--check-existing', '1'])
        self.assertFalse(res); self.assertIn('RunIDs: sheet R_A1, computed OTHER_RUN', '\n'.join(out))
        fm = FM(registered=['A1']); res, out = run(fm, ['--check-existing', '3']); self.assertTrue(res)                                                          # asks for 3, only 1 has BAMs here

    def test_backup_grid_and_failures(self):
        fm = FM(registered=['D1'], grid=3); fm.make('N1'); fm.add_reads('N1', ['R1']); fm.make('N2'); fm.add_reads('N2', ['R2'])
        res, out = run(fm, ['--samples', 'N1', 'N2']); text = '\n'.join(out)
        self.assertIn('grid enlarged from 3 to', text); self.assertEqual(len(fm.remote), 3); self.assertTrue([x for x in os.listdir(fm.localTempDir) if x.startswith('database_backup')])
        fm = FM(fail_write=True); fm.make('N1'); fm.add_reads('N1', ['R1'])
        with self.assertRaises(SystemExit) as cm:
            run(fm, ['--samples', 'N1'])
        self.assertIn('restore from', str(cm.exception))
        fm = FM(corrupt=True); fm.make('N1'); fm.add_reads('N1', ['R1'])
        with self.assertRaises(SystemExit) as cm:
            run(fm, ['--samples', 'N1'])
        self.assertIn('does not hold what was written', str(cm.exception))

    def test_a_sheet_without_a_column_for_the_row_is_refused_before_any_write(self):
        fm = FM(); fm.remote = fm.remote.drop(columns=['BamSize']); fm.make('N1'); fm.add_reads('N1', ['R1'])
        with self.assertRaises(SystemExit) as cm:
            run(fm, ['--samples', 'N1'])
        self.assertIn('BamSize', str(cm.exception)); self.assertEqual(fm.writes_n, 0)

    def test_one_failing_sample_does_not_stop_the_others(self):
        fm = FM(); fm.make('GOOD'); fm.add_reads('GOOD', ['R1']); fm.make('BAD', all=0); fm.add_reads('BAD', ['R2'])                                              # no reads in the all.bam
        res, out = run(fm, ['--samples', 'GOOD', 'BAD'])
        self.assertEqual([r['SampleID'] for r in res], ['GOOD']); self.assertEqual(fm.remote.SampleID.tolist(), ['GOOD'])
        self.assertIn('not registered: BAD: no reads in the all.bam', '\n'.join(out))

    def test_the_gatk_and_samtools_counters_parse_the_output(self):
        class P:
            def __init__(self, out, rc=0, err=''): self.stdout, self.returncode, self.stderr = out, rc, err
        orig = subprocess.run
        try:
            subprocess.run = lambda cmd, **k: P('Using GATK jar /x.jar\nTool returned:\n123456\n') if cmd[0] == 'gatk' else P('987654\n')
            self.assertEqual((R.count_gatk('f.bam'), R.count_samtools('f.bam')), (123456, 987654))
            subprocess.run = lambda cmd, **k: P('ERROR something\n', 1, 'boom')
            for fn in (R.count_gatk, R.count_samtools):
                with self.assertRaises(RuntimeError):
                    fn('f.bam')
        finally:
            subprocess.run = orig

    def test_progress_lines_and_per_file_tasks(self):
        fm = FM()
        for s in ('N1', 'N2', 'N3'): fm.make(s); fm.add_reads(s, ['R_' + s])
        seen = []
        def counting(path):
            seen.append(os.path.basename(path)); return counter(path)
        out = []; rows = R.main(['--genome', 'G1', '--samples', 'N1', 'N2', 'N3', '--jobs', '4', '--dry-run'], fm=fm, counter=counting, versions=VERS, reference_size=REF, out=out.append)
        self.assertEqual(len(seen), 21)                                                                                    # 3 samples x 7 BAMs, each its own task
        prog = [l for l in out if 'samples counted' in l]
        self.assertEqual(len(prog), 3); self.assertIn('3/3 samples counted', prog[-1]); self.assertEqual(sorted(r['SampleID'] for r in rows), ['N1', 'N2', 'N3'])

    def test_one_bam_that_cannot_be_counted_fails_only_its_sample(self):
        fm = FM()
        for s in ('GOOD', 'BAD'): fm.make(s); fm.add_reads(s, ['R_' + s])
        def flaky(path):
            if 'BAD.discordant' in path: raise RuntimeError('gatk crashed')
            return counter(path)
        out = []; rows = R.main(['--genome', 'G1', '--samples', 'GOOD', 'BAD', '--dry-run'], fm=fm, counter=flaky, versions=VERS, reference_size=REF, out=out.append)
        self.assertEqual([r['SampleID'] for r in rows], ['GOOD']); self.assertIn('not registered: BAD: gatk crashed', '\n'.join(out))

    def test_the_pysam_counter_uses_pysams_own_samtools(self):
        import types as _t
        calls = []
        fake = _t.ModuleType('pysam'); fake.view = lambda *a: (calls.append(a), '4242\n')[1]
        saved = sys.modules.get('pysam'); sys.modules['pysam'] = fake
        try:
            self.assertEqual(R.count_pysam('x.bam', threads=6), 4242); self.assertEqual(calls[0], ('-c', '-@', '6', 'x.bam'))
            fake.view = lambda *a: 'error text'
            with self.assertRaises(RuntimeError):
                R.count_pysam('x.bam')
        finally:
            if saved is None: sys.modules.pop('pysam', None)
            else: sys.modules['pysam'] = saved
        self.assertEqual(sorted(R.COUNTERS), ['gatk', 'pysam', 'samtools'])

    def test_the_audit_file_selects_by_class(self):
        p = audit_file([('A', 'needs_row'), ('B', 'ok'), ('C', 'partial_upload'), ('D', 'needs_row')])
        self.assertEqual(R.samples_from_audit(p, ['needs_row']), ['A', 'D']); self.assertEqual(R.samples_from_audit(p, ['needs_row', 'partial_upload']), ['A', 'C', 'D'])


if __name__ == '__main__':
    unittest.main()
