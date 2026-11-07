"""unit_scripts/audit_alignments.py with a fake file manager, fake listings and a fake rclone."""
import json, os, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'unit_scripts'))
import pandas as pd
import audit_alignments as A

FILES = ['.all.bam', '.all.bam.bai', '.unmapped.bam', '.discordant.bam']


def files(s, which=FILES, size=100):
    return {s + w: size for w in which}


def scenario():
    cloud, local, rows = {}, {}, set()
    for i in range(10):                                                    # ten finished samples: row and every file on Dropbox
        s = f'OK{i}'; cloud[s] = files(s); rows.add(s)
    cloud['NEEDROW'] = files('NEEDROW')                                    # uploaded, row never added
    cloud['PARTIAL'] = files('PARTIAL', FILES[:2])                          # upload stopped half way
    cloud['EMPTY'] = files('EMPTY'); cloud['EMPTY']['EMPTY.all.bam'] = 0     # an empty BAM on Dropbox
    cloud['DIFFERS'] = files('DIFFERS'); local['DIFFERS'] = files('DIFFERS', size=999)      # Dropbox and local copies disagree
    local['NEEDUP'] = files('NEEDUP')                                       # only here, complete
    local['LOCINC'] = files('LOCINC', FILES[:1])                            # only here, incomplete
    rows.add('ROWONLY'); cloud['ROWONLY'] = files('ROWONLY', FILES[:1])     # a row, but the files are not all on Dropbox
    reads = set(rows) | {'NEEDROW', 'PARTIAL', 'EMPTY', 'DIFFERS', 'NEEDUP', 'LOCINC', 'NEVER1', 'NEVER2'}
    return rows, cloud, local, reads


class FM:
    localMasterDir, cloudMasterDir = '/Data/u/Temp/CSD/', 'ptm_dropbox:/COS/CSD/'
    localBamRefDir = '/Data/u/Temp/CSD/Bamfiles/G1/'

    def __init__(self, rows, reads):
        self.localTempDir = tempfile.mkdtemp() + '/'; self.rows, self.reads = rows, reads
    def setGenome(self, g): self.g = g
    def readAlignmentDatabase(self): self.alignment_dt = pd.DataFrame(dict(SampleID=sorted(self.rows) + [None]))
    def readSampleDatabase(self):
        self.reads_dt = pd.DataFrame(dict(SampleID=sorted(self.reads))); self.sample_dt = pd.DataFrame(dict(SampleID=sorted(self.reads), Ecogroup=['Mbuna'] * len(self.reads)))


class T(unittest.TestCase):
    def test_suffixes_and_the_usual_set(self):
        self.assertEqual(A.suffix_of('S1', 'S1.all.bam'), '.all.bam'); self.assertEqual(A.suffix_of('S1', 'other.txt'), 'other.txt')
        rows, cloud, local, reads = scenario()
        self.assertEqual(A.typical_suffixes(rows, cloud), set(FILES))                                                    # the 10 finished samples decide, ROWONLY (1 of 11) does not
        cloud2 = dict(cloud); cloud2['OK0'] = dict(files('OK0'), **{'OK0.extra.vcf': 5})
        self.assertNotIn('.extra.vcf', A.typical_suffixes(rows, cloud2))                                                 # an ending only one sample has is not "usual"

    def test_every_class(self):
        rows, cloud, local, reads = scenario()
        res = A.classify(rows, cloud, local, reads, set(FILES)); cls = {s: c for s, (c, d) in res.items()}
        self.assertEqual({cls[f'OK{i}'] for i in range(10)}, {'ok'})
        self.assertEqual((cls['NEEDROW'], cls['PARTIAL'], cls['EMPTY'], cls['DIFFERS']), ('needs_row', 'partial_upload', 'partial_upload', 'partial_upload'))
        self.assertEqual((cls['NEEDUP'], cls['LOCINC'], cls['ROWONLY'], cls['NEVER1']), ('needs_upload', 'local_incomplete', 'row_only', 'not_started'))
        self.assertIn('missing on Dropbox: .discordant.bam,.unmapped.bam', res['PARTIAL'][1]); self.assertIn('empty on Dropbox: .all.bam', res['EMPTY'][1]); self.assertIn('size differs from the local copy', res['DIFFERS'][1])

    def test_the_report_and_the_tsv(self):
        rows, cloud, local, reads = scenario(); fm = FM(rows, reads); out = []
        res = A.main(['--genome', 'G1'], fm=fm, cloud_fn=lambda d: cloud, local_fn=lambda d: local, out=out.append)
        text = '\n'.join(out)
        self.assertIn('11 rows in the AlignmentDatabase', text); self.assertIn('usual files: .all.bam, .all.bam.bai, .discordant.bam, .unmapped.bam', text)
        import re
        self.assertRegex(text, r'needs_row\s+1\s+NEEDROW'); self.assertRegex(text, r'not_started\s+2\s+NEVER1 NEVER2'); self.assertRegex(text, r'row_only\s+1\s+ROWONLY')
        self.assertIn('not_started by Ecogroup: Mbuna 2', text)
        tsv = open(fm.localTempDir + 'audit_alignments_G1.tsv').read().splitlines()
        self.assertEqual(tsv[0].split('\t')[:3], ['SampleID', 'class', 'details']); self.assertEqual(len(tsv), 1 + len(res))
        self.assertIn('NEEDROW\tneeds_row\t\tMbuna\t4\t0', tsv)

    def test_the_cloud_listing_is_one_call_and_distinguishes_missing_from_failed(self):
        class R:
            def __init__(self, rc, out='', err=''): self.returncode, self.stdout, self.stderr = rc, out, err
        calls = []
        def run(cmd, **k):
            calls.append(cmd); return R(0, json.dumps([dict(Path='S1/S1.all.bam', Size=5), dict(Path='S1/S1.all.bam.bai', Size=1), dict(Path='S2/S2.all.bam', Size=7), dict(Path='stray.txt', Size=2)]))
        self.assertEqual(A.cloud_listing('remote:/Bamfiles/G1/', run=run), {'S1': {'S1.all.bam': 5, 'S1.all.bam.bai': 1}, 'S2': {'S2.all.bam': 7}})
        self.assertEqual(calls, [['rclone', 'lsjson', '-R', '--files-only', 'remote:/Bamfiles/G1/']])
        self.assertEqual(A.cloud_listing('x', run=lambda *a, **k: R(3)), {})
        with self.assertRaises(RuntimeError):
            A.cloud_listing('x', run=lambda *a, **k: R(1, '', 'network unreachable'))

    def test_the_local_listing(self):
        root = tempfile.mkdtemp() + '/'; os.makedirs(root + 'S1/sub'); open(root + 'S1/S1.all.bam', 'w').write('abc'); open(root + 'S1/sub/x.txt', 'w').write('z'); open(root + 'stray.txt', 'w').write('q')
        self.assertEqual(A.local_listing(root), {'S1': {'S1.all.bam': 3, os.path.join('sub', 'x.txt'): 1}}); self.assertEqual(A.local_listing('/nonexistent/'), {})


if __name__ == '__main__':
    unittest.main()
