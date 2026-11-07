"""unit_scripts/verify_ubam.py: the comparison, the reading of the first reads, and the sheet write-back, with fakes for rclone, samtools, ENA and the sheet."""
import os, sys, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'unit_scripts'))
import pandas as pd
from helper_modules import download_manifest as MF
import verify_ubam as V

SEQ = 'ACGTACGTAC'


def reads(n, tweak=None):
    out = [(SEQ[i % 3:] + SEQ[:i % 3], 'BBFF/<7#B1'[:10]) for i in range(n)]
    if tweak is not None:
        s, q = out[tweak]; out[tweak] = (s, q[:-1] + '#')
    return out


def sam(rs):
    return [f'r{i}\t77\t*\t0\t0\t*\t*\t0\t0\t{s}\t{q}\tXT:i:5\n' for i, (s, q) in enumerate(rs)]


class FakeFM:
    cloudMasterDir = 'ptm_dropbox:/m/'
    def __init__(self):
        self.backend = MF.MemoryBackend(); mf = self.downloadManifest('Blumer')
        mf.import_runinfo(pd.DataFrame([{c: '' for c in MF.RUNINFO_COLUMNS} | {'Run': r, 'BioSample': 's'} for r in ('ERR1', 'ERR2', 'ERR3', 'ERR4')]))
        for r, src, qc in (('ERR1', 'SRA-aws', 'skipped (ENA unreachable)'), ('ERR2', 'SRA-aws', 'identical-to-ENA'), ('ERR3', 'SRA-aws', 'skipped (ENA unreachable)'), ('ERR4', 'ENA-https', 'md5 of the ENA files')):
            mf.update(r, Source=src, QualityCheck=qc, UBamPath=f'PRJ/{r}.unmapped_marked_adapters.bam')
        mf.flush(force=True)
    def downloadManifest(self, tab, who=None, flush_seconds=15.0):
        return MF.Manifest(self.backend, tab, who or 'h-1', log=lambda *a: None, sleep=lambda s: None)
    def readDownloadTab(self, tab):
        return self.downloadManifest(tab).load() if tab == 'Blumer' else pd.DataFrame()


class T(unittest.TestCase):
    def test_compare(self):
        self.assertEqual(V.compare(reads(300), reads(300))[0], 'identical-to-ENA')
        v, d = V.compare(reads(300, tweak=150), reads(300)); self.assertEqual(v, 'MISMATCH'); self.assertIn('1 quality', d)
        shifted = reads(300)[1:]
        v, d = V.compare(reads(300), shifted); self.assertEqual(v, 'MISMATCH')                                                              # a different order (as in ERR14734399) differs
        self.assertTrue(V.compare(reads(300), [])[0].startswith('skipped')); self.assertTrue(V.compare(reads(50), reads(300))[0].startswith('skipped'))

    def test_ubam_reads_takes_sequence_and_quality_and_stops_early(self):
        got = V.ubam_reads('remote', 100, lines=lambda r: iter(sam(reads(500))))
        self.assertEqual((len(got), got[0]), (100, reads(1)[0]))
        self.assertEqual(V.ubam_reads('remote', 100, lines=lambda r: iter(['@bad line\n', 'x\ty\n'])), [])

    def run_main(self, argv, mine, theirs, fm=None):
        out, fm = [], fm or FakeFM()
        res = V.main(argv, fm=fm, ubam_reader=lambda remote, n: (out.append(remote) or mine), ena_reader=lambda run, n: theirs, out=out.append)
        return res, out, fm

    def test_skipped_runs_are_found_in_the_sheet_checked_and_written_back_only_with_write(self):
        res, out, fm = self.run_main(['--skipped'], reads(300), reads(300))
        self.assertEqual([r[1] for r in res], ['ERR1', 'ERR3']); self.assertEqual({r[2] for r in res}, {'identical-to-ENA'})                  # ERR2 was verified, ERR4 came from ENA
        self.assertIn('ptm_dropbox:/m/Reads/PRJ/ERR1.unmapped_marked_adapters.bam', out)
        self.assertTrue(fm.readDownloadTab('Blumer').set_index('Run').loc['ERR1', 'QualityCheck'].startswith('skipped'))                       # nothing written without --write
        res, out, fm = self.run_main(['--skipped', '--write'], reads(300), reads(300))
        d = fm.readDownloadTab('Blumer').set_index('Run'); self.assertEqual((d.loc['ERR1', 'QualityCheck'], d.loc['ERR3', 'QualityCheck']), ('identical-to-ENA (rechecked)',) * 2)

    def test_a_mismatch_is_flagged_and_does_not_change_the_status(self):
        res, out, fm = self.run_main(['ERR1', '--tab', 'Blumer', '--write'], reads(300, tweak=7), reads(300))
        self.assertTrue(any('ACTION' in l for l in out)); d = fm.readDownloadTab('Blumer').set_index('Run')
        self.assertEqual(d.loc['ERR1', 'QualityCheck'], 'MISMATCH (rechecked)'); self.assertIn('rechecked from the uBAM', d.loc['ERR1', 'Notes']); self.assertEqual(d.loc['ERR1', 'DownloadStatus'], 'pending')

    def test_when_ena_still_does_not_answer_nothing_is_written(self):
        res, out, fm = self.run_main(['--skipped', '--write'], reads(300), [])
        self.assertTrue(all(r[2].startswith('skipped') for r in res)); self.assertTrue(fm.readDownloadTab('Blumer').set_index('Run').loc['ERR1', 'QualityCheck'].startswith('skipped (ENA'))


if __name__ == '__main__':
    unittest.main()
