import os, random, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from validate_against_old import fingerprint, TYPES


class TestFingerprint(unittest.TestCase):
    def lines(self, n=500, seed=1):
        rng = random.Random(seed)
        return [f'r{i}\t{rng.choice([99, 147, 1024 + 99])}\tNC_1\t{rng.randint(1, 10**6)}\t60\t150M\t=\t{rng.randint(1, 10**6)}\t300\tACGT\tIIII\tNM:i:0' for i in range(n)]

    def test_order_does_not_matter_content_does(self):
        a = self.lines()
        b = list(reversed(a))
        self.assertEqual(fingerprint(a), fingerprint(b))
        c = list(a); c[7] = c[7].replace('\t60\t', '\t59\t')
        self.assertNotEqual(fingerprint(a), fingerprint(c), 'a changed MAPQ must be detected')
        d = list(a); d[7] = d[7].replace('\t99\t', '\t1123\t', 1)
        self.assertNotEqual(fingerprint(a), fingerprint(d), 'a changed FLAG (e.g. a duplicate mark) must be detected')
        self.assertEqual(fingerprint(a)[0], len(a))

    def test_sequence_qualities_and_tags_are_ignored(self):
        a = self.lines(50)
        b = [l.rsplit('\t', 3)[0] + '\tTTTT\tJJJJ\tNM:i:9' for l in a]
        self.assertEqual(fingerprint(a), fingerprint(b))


class TestTypes(unittest.TestCase):
    def test_only_the_two_bams_the_new_pipeline_writes_are_compared_by_default(self):
        self.assertEqual(TYPES, ('all', 'chimeric'))


if __name__ == '__main__':
    unittest.main()
