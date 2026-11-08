"""LG5 marker discovery on synthetic alignments and pileups: cs parsing, indel left-alignment, sharing between two assemblies, coverage/flank/divergence filters, pileup parsing, selection. No real data, no network."""
import os, random, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from sd_genotypes.lg5inv import markers as K

OFF = 19_400_000


def make_ref(n=6000, seed=3):
    r = random.Random(seed)
    s = [r.choice('ACGT') for _ in range(n)]
    s[3000:3009] = list('GATTTTTCA')                                                  # a homopolymer to test left-alignment
    for i in range(1, n):                                                              # no accidental homopolymers >= 3 outside it
        if not 2995 <= i <= 3012 and s[i] == s[i - 1] == s[i - 2 if i > 1 else 0]:
            s[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[s[i]]
    return ''.join(s)


def make_cs(ref, tstart, tend, snv=None, ins=None, dele=None):
    snv, ins, dele = snv or {}, ins or {}, dele or {}
    cs, run, t = [], 0, tstart
    flush = lambda: (cs.append(f':{run}') if run else None)
    while t < tend:
        if t in dele:
            flush(); run = 0; L = dele[t]; cs.append('-' + ref[t:t + L].lower()); t += L; continue
        if t in snv:
            flush(); run = 0; cs.append('*' + ref[t].lower() + snv[t].lower())
        else:
            run += 1
        t += 1
        if (t - 1) in ins:
            flush(); run = 0; cs.append('+' + ins[t - 1].lower())
    flush()
    return ''.join(cs)


def paf(ref, name='q', tstart=0, tend=None, de=0.01, **ev):
    tend = len(ref) if tend is None else tend
    return '\t'.join(map(str, [name, 7000, 0, tend - tstart, '-', 'LG5_off19400000', len(ref), tstart, tend, tend - tstart, tend - tstart, 60, f'de:f:{de}', 'cs:Z:' + make_cs(ref, tstart, tend, **ev)]))


def alt_of(b):
    return {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[b]


class TestCs(unittest.TestCase):
    def test_short_and_long_form_and_garbage(self):
        self.assertEqual(K.parse_cs(':10*ag:5+acg:3-tt:2'), [('=', 10), ('X', 'A', 'G'), ('=', 5), ('I', 'ACG'), ('=', 3), ('D', 'TT'), ('=', 2)])
        self.assertEqual(K.parse_cs('=ACGT*ag'), [('=', 4), ('X', 'A', 'G')])
        with self.assertRaises(ValueError):
            K.parse_cs(':10?:5')

    def test_indels_are_left_aligned_to_the_same_site_wherever_the_aligner_put_them(self):
        ref = 'GATTTTTCA'
        for anchor in (2, 3, 4, 5, 6):                                                 # insertion of a T anywhere in the run
            self.assertEqual(K.normalise(anchor, 'T', ref), (1, 'T'))
        self.assertEqual(K.normalise(4, 'T', ref), (1, 'T'))                          # deletion of one T of the run
        self.assertEqual(K.normalise(7, 'G', ref), (7, 'G'))                          # not in a repeat: unchanged
        self.assertEqual(K.normalise(0, 'A', 'AAAC')[0], -1)                          # runs off the start: cannot be anchored


class TestSharing(unittest.TestCase):
    def setUp(self):
        self.ref = make_ref()

    def lines(self, *recs):
        return list(recs)

    def test_only_variants_seen_identically_in_both_assemblies_are_candidates(self):
        r = self.ref; a1, b2, ao, bo = r[1000], r[2000], r[500], r[4500]
        A = [paf(r, 'A', snv={1000: alt_of(a1), 2000: alt_of(b2), 500: alt_of(ao)}, ins={3002: 'T'}, dele={4000: 2})]
        B = [paf(r, 'B', snv={1000: alt_of(a1), 2000: alt_of(alt_of(b2)), 4500: alt_of(bo)}, ins={3006: 'T'}, dele={4000: 2})]
        df, s = K.shared_candidates({'A': A, 'B': B}, r, OFF)
        got = {(int(p), ref_, alt) for p, ref_, alt in zip(df.pos, df.ref, df.alt)}
        want = {(OFF + 1001, a1, alt_of(a1)),                                          # the shared SNV
                (OFF + 3001 + 1, 'A', 'AT')}                                           # the same insertion placed at different anchors in the homopolymer: normalised to one site
        self.assertTrue(want <= got, got)
        self.assertEqual(sorted(df.type), sorted(['SNV', 'INS', 'DEL']))
        self.assertNotIn(OFF + 2001, df.pos.tolist())                                  # same site, different alt allele: not shared
        self.assertNotIn(OFF + 501, df.pos.tolist()); self.assertNotIn(OFF + 4501, df.pos.tolist())   # private to one assembly
        self.assertEqual(s['candidates'], 3)

    def test_a_site_covered_twice_in_one_assembly_is_not_a_candidate(self):
        r = self.ref; b = r[2500]
        A = [paf(r, 'A', snv={2500: alt_of(b)})]
        B = [paf(r, 'B1', snv={2500: alt_of(b)}), paf(r, 'B2', tstart=2300, tend=2800)]        # a second alignment over the same stretch (a duplicated copy)
        df, s = K.shared_candidates({'A': A, 'B': B}, r, OFF)
        self.assertEqual(len(df), 0); self.assertEqual(s['shared_before_coverage_filter'], 1)

    def test_variants_without_enough_matching_flank_are_dropped(self):
        r = self.ref
        ev = {1000: alt_of(r[1000]), 1005: alt_of(r[1005]), 2000: alt_of(r[2000]), 5: alt_of(r[5])}   # 1000/1005 are a cluster; 5 is 5 bp from the alignment start
        A = [paf(r, 'A', snv=ev)]; B = [paf(r, 'B', snv=ev)]
        df, _ = K.shared_candidates({'A': A, 'B': B}, r, OFF)
        self.assertEqual(df.pos.tolist(), [OFF + 2001])

    def test_a_divergent_alignment_gives_no_variants_but_still_counts_as_coverage(self):
        r = self.ref; b = r[1000]
        v, cov, s = K.assembly_variants([paf(r, 'A', de=0.2, snv={1000: alt_of(b)})], r, OFF)
        self.assertEqual((len(v), int(cov[1000]), s['too_divergent']), (0, 1, 1))

    def test_a_cs_that_disagrees_with_the_target_is_an_error_not_a_silent_wrong_answer(self):
        r = self.ref; wrong = next(c for c in 'ACGT' if c != r[100]).lower()
        with self.assertRaises(ValueError):
            K.alignment_variants('\t'.join(['q', '1', '0', '200', '+', 'T', '6000', '0', '200', '200', '200', '60', f'cs:Z::100*{wrong}c:99']), r, OFF)
        with self.assertRaises(ValueError):                                            # cs shorter than the PAF coordinates
            K.alignment_variants('\t'.join(['q', '1', '0', '200', '+', 'T', '6000', '0', '200', '200', '200', '60', 'cs:Z::150']), r, OFF)
        with self.assertRaises(ValueError):
            K.alignment_variants('\t'.join(['q', '1', '0', '200', '+', 'T', '6000', '0', '200', '200', '200', '60']), r, OFF)

    def test_end_to_end_cli_writes_the_candidates(self):
        d = tempfile.mkdtemp(); r = self.ref; b = r[1500]
        open(d + '/t.fa', 'w').write('>LG5_off19400000\n' + '\n'.join(r[i:i + 80] for i in range(0, len(r), 80)) + '\n')
        for n in 'AB':
            open(f'{d}/{n}.paf', 'w').write(paf(r, n, snv={1500: alt_of(b)}) + '\n')
        K.main(['candidates', '--target-fa', d + '/t.fa', '--paf', f'A={d}/A.paf', '--paf', f'B={d}/B.paf', '--out', d + '/c.tsv'])
        c = pd.read_csv(d + '/c.tsv', sep='\t')
        self.assertEqual((c.chrom.tolist(), c.pos.tolist(), c.ref.tolist(), c.alt.tolist()), (['NC_135171.1'], [OFF + 1501], [b], [alt_of(b)]))


class TestCounting(unittest.TestCase):
    def test_pileup_string_parsing(self):
        c = K.pileup_counts('.,.+2AT,-1ga^Fc$*T', 'G')
        self.assertEqual((c['G'], c['A'], c['C'], c['T'], c['*']), (4, 1, 1, 1, 1))

    def test_counts_per_site_without_a_base_quality_filter_and_a_wrong_reference_is_refused(self):
        snv = pd.DataFrame(dict(pos=[100, 200, 300], ref=['A', 'C', 'G'], alt=['T', 'G', 'A']))
        seen = {}

        def runner(ref, bam, bed, mapq):
            seen['mapq'] = mapq; seen['bed'] = open(bed).read()
            return ['\t'.join(['NC_135171.1', '100', 'A', '10', '....,,,TTt', 'IIIIIIIIII']), '\t'.join(['NC_135171.1', '200', 'C', '4', '.,.,', 'IIII'])]   # 300 has no coverage
        d, a, rc = K.count_sites(snv, 'x.bam', 'ref.fa', 20, runner)
        self.assertEqual((d.tolist(), a.tolist(), rc.tolist()), ([10, 4, 0], [3, 0, 0], [7, 4, 0]))
        self.assertEqual(seen['mapq'], 20); self.assertEqual(seen['bed'].split('\n')[0], 'NC_135171.1\t99\t100')
        with self.assertRaises(ValueError):
            K.count_sites(snv.iloc[:1], 'x.bam', 'ref.fa', 20, lambda *a: ['\t'.join(['NC_135171.1', '100', 'C', '10', '..........', 'I' * 10])])

    def test_selection_keeps_sites_clean_in_negatives_and_near_half_in_both_positives(self):
        t = pd.DataFrame(dict(pos=[1, 2, 3, 4, 5],
                              BB_depth=[30, 30, 30, 30, 10], BB_alt=[0, 1, 5, 0, 0],                      # 3 has alt in the BB pool, 5 is too shallow
                              OB1_depth=[40, 40, 40, 40, 40], OB1_alt=[20, 18, 20, 3, 20],               # 4 is near 0 in the OB pool
                              OB2_depth=[40, 40, 40, 40, 40], OB2_alt=[15, 20, 20, 3, 20]))
        m = K.select(t, ['BB'], ['OB1', 'OB2'])
        self.assertEqual(m.pos.tolist(), [1, 2])
        self.assertAlmostEqual(m.mean_positive_fraction[0], (0.5 + 0.375) / 2)


if __name__ == '__main__':
    unittest.main()
