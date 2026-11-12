"""region_block (module 06): the statistics on data with a known answer, the VCF reader, and the command end to end on a planted XY block with a female run of homozygosity."""
import os, sys, tempfile, unittest, io, json, shutil
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from sd_regions.region_block import analysis as A


def haplotypes(rng, m, k, freq=0.3):
    return (rng.random((m, k)) < freq).astype(np.int8)


class TestWindows(unittest.TestCase):
    def test_accumulator_equals_a_naive_count_whatever_the_chunking(self):
        rng = np.random.default_rng(1)
        pos = np.sort(rng.integers(1, 1_000_000, 3000)); chrom = np.where(np.arange(3000) < 1800, 'c1', 'c2')
        pos = np.concatenate([np.sort(pos[:1800]), np.sort(pos[1800:])])
        dos = rng.integers(-1, 3, (3000, 5)).astype(np.int8)
        acc1, acc2 = A.WindowAccumulator(5, 100_000), A.WindowAccumulator(5, 100_000)
        acc1.add(chrom, pos, dos)
        for lo, hi in ((0, 700), (700, 1900), (1900, 3000)):                                  # the same data in three chunks, one of them straddling the chromosomes
            acc2.add(chrom[lo:hi], pos[lo:hi], dos[lo:hi])
        k1, c1, h1 = acc1.arrays(); k2, c2, h2 = acc2.arrays()
        self.assertEqual(k1, k2); np.testing.assert_array_equal(c1, c2); np.testing.assert_array_equal(h1, h2)
        for i, (ch, w) in enumerate(k1):
            sel = (chrom == ch) & (pos // 100_000 == w)
            np.testing.assert_array_equal(c1[i], (dos[sel] >= 0).sum(0)); np.testing.assert_array_equal(h1[i], (dos[sel] == 1).sum(0))
        self.assertEqual([k[0] for k in k1], sorted([k[0] for k in k1]))                       # chromosome order is the order of appearance

    def test_an_empty_chunk_is_harmless(self):
        acc = A.WindowAccumulator(3, 1000); acc.add(np.array([]), np.array([], dtype=np.int64), np.zeros((0, 3), np.int8))
        self.assertEqual(acc.arrays()[0], [])

    def test_window_het_is_nan_below_the_minimum_number_of_sites(self):
        h = A.window_het(np.array([[100, 10]]), np.array([[24, 5]]), 30)
        self.assertAlmostEqual(h[0, 0], 0.24); self.assertTrue(np.isnan(h[0, 1]))


class TestRoh(unittest.TestCase):
    def setUp(self):
        self.size = 250_000
        self.keys = [('c1', w) for w in range(40)]
        self.called = np.full((40, 3), 100); self.het = np.full((40, 3), 24)
        self.het[5:15, 0] = 0                                                      # sample 0: ten homozygous windows
        self.het[5:8, 1] = 0                                                       # sample 1: three: too short
        self.het[10:20, 2] = 0; self.het[14, 2] = 30                               # sample 2: broken in the middle: 10-13 and 15-19
        self.called[25, 0] = 10                                                    # a window with too few sites in sample 0 (heterozygous anyway)

    def test_runs_have_the_right_coordinates_and_short_ones_are_dropped(self):
        segs, ok = A.roh_segments(self.keys, self.called, self.het, self.size, het_max=0.03, min_sites=30, min_windows=4)
        got = sorted((g['sample'], g['start'], g['end'], g['n_windows']) for g in segs)
        self.assertEqual(got, [(0, 5 * self.size, 15 * self.size, 10), (2, 10 * self.size, 14 * self.size, 4), (2, 15 * self.size, 20 * self.size, 5)])
        self.assertFalse(ok[25, 0]); self.assertTrue(ok[25, 1])

    def test_a_window_with_too_few_sites_cannot_be_homozygous(self):
        self.het[22:30, 0] = 0                                                     # windows 22-29 have no heterozygotes, but window 25 holds only 10 calls
        segs, ok = A.roh_segments(self.keys, self.called, self.het, self.size)
        s0 = sorted((g['start'] // self.size, g['n_windows']) for g in segs if g['sample'] == 0)
        self.assertEqual(s0, [(5, 10), (26, 4)], '22-24 (3 windows) is too short, 25 cannot be judged, 26-29 is a run')

    def test_froh_is_the_share_of_the_callable_genome(self):
        segs, ok = A.roh_segments(self.keys, self.called, self.het, self.size)
        f = A.froh(segs, ok, self.size, 3)
        self.assertAlmostEqual(f[0], 10 / 39); self.assertAlmostEqual(f[1], 0.0); self.assertAlmostEqual(f[2], 9 / 40)

    def test_runs_do_not_cross_chromosomes(self):
        keys = [('c1', w) for w in range(5)] + [('c2', w) for w in range(5)]
        called, het = np.full((10, 1), 100), np.zeros((10, 1), int)
        segs, _ = A.roh_segments(keys, called, het, self.size, min_windows=4)
        self.assertEqual(sorted((g['chrom'], g['n_windows']) for g in segs), [('c1', 5), ('c2', 5)])

    def test_a_missing_window_breaks_a_run(self):
        keys = [('c1', w) for w in (0, 1, 2, 3, 5, 6, 7, 8)]                       # window 4 holds no SNP at all: it is not in the list
        called, het = np.full((8, 1), 100), np.zeros((8, 1), int)
        segs, _ = A.roh_segments(keys, called, het, self.size, min_windows=4)
        self.assertEqual(sorted((g['start'] // self.size, g['n_windows']) for g in segs), [(0, 4), (5, 4)], 'two runs of four, not one of eight')

    def test_region_summary(self):
        segs, ok = A.roh_segments(self.keys, self.called, self.het, self.size)
        s = A.region_summary(self.keys, self.called, self.het, segs, ok, 'c1', 10 * self.size, 20 * self.size, self.size)
        self.assertEqual(s['n_windows'], 10)
        self.assertAlmostEqual(s['roh_share'][0], 0.5); self.assertAlmostEqual(s['hom_share'][0], 0.5)       # sample 0: windows 10-14 homozygous, 15-19 not
        self.assertAlmostEqual(s['roh_share'][2], 0.9)                                                          # sample 2: 10-13 and 15-19 in runs; 14 is the break
        self.assertAlmostEqual(s['het'][1], 0.24)
        self.assertEqual(A.region_summary(self.keys, self.called, self.het, segs, ok, 'zz', 0, 10, self.size)['n_windows'], 0)

    def test_chromosome_sex_het(self):
        keys = [('c1', 0), ('c1', 1), ('c2', 0)]
        called = np.array([[100, 100, 100, 100]] * 3); het = np.array([[50, 50, 10, 10], [50, 50, 10, 10], [10, 10, 10, 10]])
        male = np.array([True, True, False, False])
        r = A.chrom_sex_het(keys, called, het, male, ~male)
        self.assertAlmostEqual(r['c1'][0], 0.5); self.assertAlmostEqual(r['c1'][1], 0.1); self.assertAlmostEqual(r['c1'][2], 0.4); self.assertAlmostEqual(r['c2'][2], 0.0)


class TestRelatedness(unittest.TestCase):
    def build(self, seed=2):
        """24 fish: 8 unrelated, family A (4 full sibs), family B (4 full sibs), 8 more unrelated; and one inbred fish (both haplotypes identical) in place of the last unrelated one."""
        rng = np.random.default_rng(seed); m = 20000
        pool = haplotypes(rng, m, 80, 0.35); nxt = iter(range(80))
        def ind(): return pool[:, next(nxt)] + pool[:, next(nxt)]
        def fam(k):
            pa, pb = [pool[:, next(nxt)] for _ in range(2)], [pool[:, next(nxt)] for _ in range(2)]
            kids = []
            for _ in range(k):
                ha = pa[rng.integers(0, 2)]; hb = pb[rng.integers(0, 2)]
                kids.append(ha + hb)
            return kids
        # families built from CHROMOSOME-wide parental haplotypes: whole-haplotype transmission makes sibs share about half their alleles
        cols = [ind() for _ in range(8)] + fam(4) + fam(4) + [ind() for _ in range(7)]
        h = pool[:, next(nxt)]; cols.append(2 * h)                                           # the inbred fish
        return np.stack(cols, axis=1).astype(np.int8)

    def test_siblings_are_more_related_and_the_inbred_fish_has_F_near_one(self):
        dos = self.build(); acc = A.GrmAccumulator(24, 0.05); acc.add(dos)
        G = acc.grm()
        self.assertGreater(G[8:12, 8:12][np.triu_indices(4, 1)].mean(), 0.3)                 # family A
        self.assertGreater(G[12:16, 12:16][np.triu_indices(4, 1)].mean(), 0.3)               # family B
        self.assertLess(abs(G[:8, :8][np.triu_indices(8, 1)].mean()), 0.08)                  # unrelated
        self.assertGreater(acc.f_het()[23], 0.8); self.assertLess(abs(acc.f_het()[:8].mean()), 0.1)
        self.assertGreater(G[23, 23], 1.5)                                                   # G_ii = 1 + F

    def test_chunks_give_the_same_matrix(self):
        dos = self.build(); a, b = A.GrmAccumulator(24), A.GrmAccumulator(24); a.add(dos)
        for lo in range(0, 20000, 7000): b.add(dos[lo:lo + 7000])
        np.testing.assert_allclose(a.grm(), b.grm(), atol=1e-9); self.assertEqual(a.n_sites, b.n_sites)

    def test_sex_test_detects_a_family_confounded_with_sex(self):
        dos = self.build(); acc = A.GrmAccumulator(24); acc.add(dos); G = acc.grm()
        sex = np.array(['M'] * 12 + ['F'] * 12)                                              # family A (8-11) males; family B (12-15) females: each family is one sex
        r = A.kinship_sex_test(G, sex, n_perm=2000)
        self.assertLess(r['p'], 0.02); self.assertGreater(r['within_m'], r['between'])
        mixed = np.array((['M', 'F'] * 12)); r2 = A.kinship_sex_test(G, mixed, n_perm=2000)  # sexes alternate: every family has both
        self.assertGreater(r2['p'], 0.1)

    def test_clusters_are_connected_components(self):
        G = np.zeros((6, 6)); np.fill_diagonal(G, 1)
        for i, j in ((0, 1), (1, 2), (4, 5)): G[i, j] = G[j, i] = 0.4
        self.assertEqual(A.related_clusters(G, 0.2), [[0, 1, 2], [4, 5], [3]])
        self.assertEqual(A.related_clusters(G, 0.5), [[0], [1], [2], [3], [4], [5]])


def planted_xy(seed=3, m=400, n_m=12, n_f=12, block_frac=1.0):
    """Region SNPs with a Y haplotype: every female XX (one X haplotype, so homozygous), every male X/Y with the Y differing from the X at the planted SNPs. Returns dos (m x n), male mask."""
    rng = np.random.default_rng(seed)
    X = (rng.random(m) < 0.4).astype(np.int8)
    Y = X.copy(); flip = rng.random(m) < block_frac; Y[flip] = 1 - Y[flip]
    cols = [X + Y for _ in range(n_m)] + [2 * X for _ in range(n_f)]
    dos = np.stack(cols, axis=1).astype(np.int8)
    return dos, np.array([True] * n_m + [False] * n_f)


class TestStructure(unittest.TestCase):
    def test_pca_separates_the_sexes_when_there_is_a_sex_haplotype_and_not_without(self):
        dos, male = planted_xy()
        rng = np.random.default_rng(5); noise = rng.integers(0, 3, dos.shape).astype(np.int8)
        sc, var, used = A.pca_scores(np.vstack([dos, noise[:200]]))
        a, p, pc = A.best_pc_auc(sc, male, n_perm=2000)
        self.assertEqual(a, 1.0); self.assertLess(p, 0.001)
        null = rng.integers(0, 3, (600, 24)).astype(np.int8)
        a0, p0, _ = A.best_pc_auc(A.pca_scores(null)[0], male, n_perm=2000)
        self.assertLess(a0, 0.9)

    def test_auc_extremes_and_fold(self):
        m = np.array([True] * 5 + [False] * 5)
        self.assertEqual(A.auc_sex(np.arange(10.0), m, 500)[0], 1.0); self.assertEqual(A.auc_sex(-np.arange(10.0), m, 500)[0], 1.0)   # either direction is perfect
        a, p = A.auc_sex(np.array([1, 2, 1, 2, 1, 2, 1, 2, 1, 2.0]), m, 500); self.assertTrue(0.5 <= a <= 1)
        self.assertTrue(np.isnan(A.auc_sex(np.arange(4.0), np.array([True] * 4))[0]))

    def test_pca_needs_enough_snps(self):
        with self.assertRaises(ValueError):
            A.pca_scores(np.zeros((5, 24), np.int8))

    def test_pca_from_the_grm_matches_pca_from_genotypes(self):
        dos = TestRelatedness().build(); acc = A.GrmAccumulator(24); acc.add(dos)
        s1, _ = A.pca_from_grm(acc.grm()); s2, _, _ = A.pca_scores(dos)
        self.assertGreater(abs(np.corrcoef(s1[:, 0], s2[:, 0])[0, 1]), 0.95)


class TestLd(unittest.TestCase):
    def test_perfectly_linked_sites_give_r2_one_and_independent_ones_near_zero(self):
        rng = np.random.default_rng(1)
        a = rng.integers(0, 3, 24); dos = np.vstack([a, a, 2 - a, rng.integers(0, 3, 24), rng.integers(0, 3, 24)]).astype(np.int8)
        r2 = A.ld_matrix(dos)
        self.assertAlmostEqual(r2[0, 1], 1.0); self.assertAlmostEqual(r2[0, 2], 1.0); self.assertLess(r2[0, 3], 0.5)
        self.assertTrue(np.isnan(A.ld_matrix(np.array([[1] * 24, a], dtype=np.int8))[0, 1]))                   # a site that does not vary

    def test_missing_genotypes_do_not_break_it(self):
        rng = np.random.default_rng(2); a = rng.integers(0, 3, 24); d = np.vstack([a, a]).astype(np.int8); d[0, :3] = -1
        self.assertGreater(A.ld_matrix(d)[0, 1], 0.7)

    def test_thinning_keeps_even_spacing_and_the_frequency_filter(self):
        rng = np.random.default_rng(3); dos = rng.integers(0, 3, (1000, 24)).astype(np.int8); dos[::2] = 0; pos = np.arange(1000) * 100
        d, p = A.thin_sites(dos, pos, min_maf=0.15, max_sites=100)
        self.assertEqual(len(p), 100); self.assertTrue(np.all(np.diff(p) > 0)); self.assertTrue((p % 200 != 0).all())    # the monomorphic even sites are gone

    def test_a_site_heterozygous_in_every_fish_has_no_ld_to_measure_and_is_not_chosen(self):
        rng = np.random.default_rng(4); dos = rng.integers(0, 3, (60, 24)).astype(np.int8); dos[::3] = 1                      # every third site: all fish heterozygous (frequency 0.5, no variation)
        d, p = A.thin_sites(dos, np.arange(60) * 100, min_maf=0.15, max_sites=100)
        self.assertTrue(all(q % 300 != 0 for q in p)); self.assertGreater(len(p), 20)
        self.assertEqual(len(A.thin_sites(np.ones((50, 24), np.int8), np.arange(50), 0.15, 100)[1]), 0)

    def test_decay_bins(self):
        pos = np.array([0, 100, 5000, 9000]); r2 = np.array([[1, .9, .1, .1], [.9, 1, .2, .1], [.1, .2, 1, .5], [.1, .1, .5, 1]])
        out = A.ld_by_distance(r2, pos, [0, 1000, 6000, 20000])
        self.assertAlmostEqual(out[0], 0.9)                                        # the pair 100 bp apart
        self.assertAlmostEqual(out[1], np.mean([0.1, 0.2, 0.5]))                   # 5000, 4900 and 4000 bp apart
        self.assertAlmostEqual(out[2], 0.1)                                        # 9000 and 8900 bp apart
        self.assertTrue(np.isnan(A.ld_by_distance(r2, pos, [20000, 30000])[0]))


class TestDiagnostic(unittest.TestCase):
    def test_the_rule_depends_on_the_number_of_fish_and_the_strict_rule_does_not_forgive(self):
        # one site: 11 of 12 males heterozygous (the 12th homozygous), 12 of 12 females homozygous reference
        dos = np.array([[1] * 11 + [0] + [0] * 12], dtype=np.int8)
        male = np.array([True] * 12 + [False] * 12)
        self.assertTrue(A.diag_sexloci(dos, male, ~male)[0])                                  # 11/12 = 0.917 >= 0.9: allowed at 12 per sex
        sub_m = np.array([True] * 8 + [False] * 3 + [True] + [False] * 12)                    # 9 males including the exception: 8/9 = 0.889 < 0.9
        sub_f = np.array([False] * 12 + [True] * 7 + [False] * 5)
        self.assertFalse(A.diag_sexloci(dos, sub_m, sub_f)[0])
        self.assertFalse(A.diag_strict(dos, male, ~male)[0])
        # one heterozygous female among 12 is allowed (1/12 = 0.083, maf 0.042), among 7 it is not (1/7 = 0.143)
        d2 = np.array([[1] * 12 + [1] + [0] * 11], dtype=np.int8)
        self.assertTrue(A.diag_sexloci(d2, male, ~male)[0]); self.assertFalse(A.diag_sexloci(d2, male, np.array([False] * 12 + [True] * 7 + [False] * 5))[0])

    def test_strict_rule_and_planted_block(self):
        dos, male = planted_xy(m=300, block_frac=0.5)
        strict = A.diag_strict(dos, male, ~male); loose = A.diag_sexloci(dos, male, ~male)
        self.assertEqual(int(strict.sum()), int(loose.sum())); self.assertGreater(int(strict.sum()), 100)      # a clean Y haplotype: both rules agree
        self.assertLess(int(strict.sum()), 300)

    def test_the_count_is_much_lower_in_the_small_runs_subsets(self):
        rng = np.random.default_rng(7); dos, male = planted_xy(m=500, block_frac=0.6)
        err = rng.random(dos.shape) < 0.012; dos = dos.copy(); dos[err] = rng.integers(0, 3, int(err.sum())).astype(np.int8)        # genotype errors, as at low depth
        full = int(A.diag_sexloci(dos, male, ~male).sum())
        small = A.diag_subsample(dos, male, ~male, 9, 7, n_draws=60, seed=1)
        self.assertEqual(len(small), 60); self.assertLess(np.median(small), full)
        np.testing.assert_array_equal(small, A.diag_subsample(dos, male, ~male, 9, 7, n_draws=60, seed=1))                      # reproducible

    def test_depth_masking(self):
        dos = np.array([[1] * 12 + [0] * 12], dtype=np.int8); dp = np.array([[6.0] * 12 + [9.0] * 12], dtype=np.float32); male = np.array([True] * 12 + [False] * 12)
        r = A.diag_by_depth(dos, dp, male, ~male, thresholds=(3, 5, 8))
        self.assertEqual(r[3], (1, 1)); self.assertEqual(r[5], (1, 1)); self.assertEqual(r[8], (0, 0))                          # males at depth 6 are masked at 8
        dp[:, :12] = 5.0                                                                                                        # a depth EQUAL to the threshold is kept (min_dp 5 keeps depth 5)
        self.assertEqual(A.diag_by_depth(dos, dp, male, ~male, thresholds=(5,))[5], (1, 1)); self.assertEqual(A.diag_by_depth(dos, dp, male, ~male, thresholds=(6,))[6], (0, 0))

    def test_per_sample_het_at_the_diagnostic_sites(self):
        dos, male = planted_xy(m=100, block_frac=1.0); sites = A.diag_sexloci(dos, male, ~male)
        h = A.diag_sample_het(dos, sites); self.assertTrue((h[male] == 1.0).all()); self.assertTrue((h[~male] == 0.0).all())




# ====================================================================== the VCF reader and the command, end to end
from sd_regions.region_block import reader as RD, cli as CLI


def write_planted_vcf(path, seed=11, n_m=12, n_f=12, block=(1_000_000, 3_000_000), c1_len=4_000_000, c2_len=2_000_000, step=1000, low_dp_share=0.02):
    """Two chromosomes, SNPs every `step` bases. c1 block: females XX with ONE X haplotype (so homozygous), males X/Y with the Y differing from the X at every SNP of the block. Elsewhere unrelated fish
    drawn from allele frequencies. The first female (F01) is homozygous over c1 0-1 Mb (a run of homozygosity). Depth is about 8; a few genotypes get depth 1 (masked at min_dp 3). Returns (names, sex)."""
    rng = np.random.default_rng(seed)
    names = [f'M{i + 1:02d}' for i in range(n_m)] + [f'F{i + 1:02d}' for i in range(n_f)]
    sex = np.array(['M'] * n_m + ['F'] * n_f); n = len(names)
    lines = ['##fileformat=VCFv4.2', '##INFO=<ID=MQ,Number=1,Type=Float,Description="mq">', '##INFO=<ID=DP,Number=1,Type=Integer,Description="dp">', '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">',
             '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp">', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(names)]
    for chrom, length in (('c1', c1_len), ('c2', c2_len)):
        for pos in range(step, length, step):
            in_block = chrom == 'c1' and block[0] <= pos < block[1]
            if in_block:
                x = int(rng.random() < 0.5)
                g = np.array([x + (1 - x)] * n_m + [2 * x] * n_f)                       # males X + Y (Y is the other allele): heterozygous; females XX: homozygous
            else:
                q = rng.uniform(0.15, 0.5)
                g = (rng.random((n, 2)) < q).sum(axis=1)
                if chrom == 'c1' and pos < 1_000_000:
                    g[n_m] = 2 * int(rng.random() < q)                                    # F01: homozygous over c1 0-1 Mb
                if g.min() == g.max():
                    continue                                                              # a site the cohort does not vary at is not in the VCF
            dp = rng.poisson(8, n).clip(4, None)
            dp[rng.random(n) < low_dp_share] = 1
            gt = {0: '0/0', 1: '0/1', 2: '1/1'}
            lines.append(f'{chrom}\t{pos}\t.\tA\tC\t60\tPASS\tMQ=60;DP={int(dp.sum())}\tGT:DP\t' + '\t'.join(f'{gt[int(v)]}:{int(d)}' for v, d in zip(g, dp)))
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    return names, sex


class TestRegionAndSamples(unittest.TestCase):
    def test_region_parsing(self):
        self.assertEqual(RD.parse_region('NC_135176.1:10,000,000-30,000,000'), ('NC_135176.1', 10_000_000, 30_000_000))
        for bad in ('c1', 'c1:5', 'c1:9-3', ''):
            with self.assertRaises(ValueError):
                RD.parse_region(bad)

    def test_samples_come_from_the_qc_table_and_other_sexes_are_dropped(self):
        d = tempfile.mkdtemp()
        try:
            p = os.path.join(d, 'sample_qc.tsv'); open(p, 'w').write('sample\tsex\trole\tmean_depth\nA\tM\ta\t5\nB\tF\ta\t6\nC\t\ta\t7\nD\tU\ta\t7\nE\t F \ta\t8\n')
            names, sex, t = CLI.load_samples(qc=p)
            self.assertEqual(names, ['A', 'B', 'E']); self.assertEqual(list(sex), ['M', 'F', 'F'])
            with self.assertRaises(ValueError): CLI.load_samples()
            q = os.path.join(d, 'x.tsv'); open(q, 'w').write('id\tgender\n1\tM\n')
            with self.assertRaises(ValueError): CLI.load_samples(sex_file=q)
        finally:
            shutil.rmtree(d, ignore_errors=True)


class TestStream(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp(); cls.vcf = os.path.join(cls.dir, 'planted.vcf')
        cls.names, cls.sex = write_planted_vcf(cls.vcf)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_the_reader_applies_the_min_dp_mask_and_keeps_the_region_chromosome(self):
        col = RD.stream(self.vcf, self.names, ('c1', 1_000_000, 3_000_000), min_dp=3, backend='builtin', log=lambda *_: None)
        self.assertEqual(col.n_seen, col.n_used if col.n_used == col.n_seen else col.n_seen)
        n_c1 = int(sum(1 for l in open(self.vcf) if l.startswith('c1\t')))
        self.assertLessEqual(len(col.chr_pos), n_c1); self.assertGreater(len(col.chr_pos), n_c1 * 0.9)      # a few sites can fall to the call-rate rule
        d3 = RD.masked(col.chr_raw, col.chr_dp, 3); d5 = RD.masked(col.chr_raw, col.chr_dp, 9)
        self.assertGreater((d5 < 0).sum(), (d3 < 0).sum())                                                 # a stricter depth masks more
        self.assertTrue(np.all(np.nan_to_num(col.chr_dp[d3 >= 0], nan=99) >= 3))
        keys, called, het = col.windows.arrays(); self.assertEqual({k[0] for k in keys}, {'c1', 'c2'})

    def test_sites_with_too_few_genotypes_at_min_dp_are_dropped(self):
        d = tempfile.mkdtemp()
        try:
            names, sex = write_planted_vcf(os.path.join(d, 'a.vcf'), c1_len=100_000, c2_len=50_000)
            lines = open(os.path.join(d, 'a.vcf')).read().rstrip('\n').split('\n')
            head, body = [l for l in lines if l.startswith('#')], [l for l in lines if not l.startswith('#')]
            bad = []
            for pos in (5_001, 5_002, 5_003):                                                         # three sites where 60% of the fish have depth 1: below min_dp 3 they are missing
                f = body[0].split('\t'); f[1] = str(pos)
                for i in range(9, 9 + 15):
                    f[i] = f[i].split(':')[0] + ':1'
                bad.append('\t'.join(f))
            open(os.path.join(d, 'b.vcf'), 'w').write('\n'.join(head + sorted(body + bad, key=lambda l: (l.split('\t')[0], int(l.split('\t')[1])))) + '\n')
            good = RD.stream(os.path.join(d, 'a.vcf'), names, ('c1', 1, 10), backend='builtin', log=lambda *_: None)
            col = RD.stream(os.path.join(d, 'b.vcf'), names, ('c1', 1, 10), backend='builtin', log=lambda *_: None)
            self.assertEqual(col.n_seen, good.n_seen + 3); self.assertEqual(col.n_used, good.n_used)         # seen, but not used
            col0 = RD.stream(os.path.join(d, 'b.vcf'), names, ('c1', 1, 10), min_call_rate=0.3, backend='builtin', log=lambda *_: None)
            self.assertEqual(col0.n_used, good.n_used + 3)                                                    # with a lenient rule they are used
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_the_region_chromosome_stays_out_of_the_relatedness_unless_asked(self):
        n_c2 = int(sum(1 for l in open(self.vcf) if l.startswith('c2\t')))
        col = RD.stream(self.vcf, self.names, ('c1', 1_000_000, 3_000_000), backend='builtin', log=lambda *_: None)
        self.assertLessEqual(col.grm_out.n_sites, n_c2); self.assertGreater(col.grm_out.n_sites, 0.9 * n_c2 * 0.5)
        keep = RD.stream(self.vcf, self.names, ('c1', 1_000_000, 3_000_000), keep_region_chromosome=True, backend='builtin', log=lambda *_: None)
        self.assertGreater(keep.grm_out.n_sites, col.grm_out.n_sites + 1000)                                  # the flanks of c1 are added, the region itself is still out
        ex = RD.stream(self.vcf, self.names, ('c1', 1_000_000, 3_000_000), exclude=[('c2', 0, 1_000_000)], backend='builtin', log=lambda *_: None)
        self.assertLess(ex.grm_out.n_sites, col.grm_out.n_sites)                                              # an excluded stretch is out too

    def test_the_boundaries_of_the_region_are_exact(self):
        base = RD.stream(self.vcf, self.names, ('c1', 1_000_000, 3_000_000), keep_region_chromosome=True, backend='builtin', log=lambda *_: None).grm_out.n_sites
        start_out = RD.stream(self.vcf, self.names, ('c1', 1_000_001, 3_000_000), keep_region_chromosome=True, backend='builtin', log=lambda *_: None).grm_out.n_sites
        end_in = RD.stream(self.vcf, self.names, ('c1', 1_000_000, 3_000_001), keep_region_chromosome=True, backend='builtin', log=lambda *_: None).grm_out.n_sites
        self.assertEqual(start_out, base + 1, 'the SNP at the start coordinate is INSIDE the region; one base later it is outside and joins the relatedness')
        self.assertEqual(end_in, base - 1, 'the SNP at the end coordinate is OUTSIDE the region; one base later it is inside and leaves the relatedness')

    def test_unknown_samples_are_an_error(self):
        with self.assertRaises(ValueError):
            RD.stream(self.vcf, self.names + ['NOPE'], ('c1', 1, 10), backend='builtin', log=lambda *_: None)


class TestCommandEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp(); cls.vcf = os.path.join(cls.dir, 'planted.vcf'); cls.out = os.path.join(cls.dir, 'out')
        cls.names, cls.sex = write_planted_vcf(cls.vcf)
        cls.qc = os.path.join(cls.dir, 'sample_qc.tsv')
        with open(cls.qc, 'w') as fh:
            fh.write('sample\tsex\trole\tmean_depth\n' + ''.join(f'{n}\t{s}\ta\t8\n' for n, s in zip(cls.names, cls.sex)))
        args = CLI.build_parser().parse_args(['--vcf', cls.vcf, '--qc', cls.qc, '--region', 'c1:1000000-3000000', '--core', '1250000-2750000', '--outdir', cls.out, '--backend', 'builtin', '--n-perm', '400',
                                              '--n-draws', '30', '--subregion', '1000000-2000000', '--subregion', '2000000-3000000', '--ld-sites', '120'])
        cls.r = CLI.run(args, log=lambda *_: None)
        cls.st = __import__('pandas').read_csv(os.path.join(cls.out, 'samples.tsv'), sep='\t')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_all_outputs_are_written(self):
        for f in ('summary.md', 'samples.tsv', 'roh_segments.tsv', 'window_het.tsv', 'kinship.tsv', 'kinship_pairs.tsv', 'chromosome_sex_het.tsv', 'pca_region.tsv', 'pca_genome.tsv', 'ld_decay.tsv',
                  'diagnostic.tsv', 'diagnostic_subsample.tsv', 'subregions.tsv', 'run_info.json', 'het_heatmap.png', 'kinship.png', 'pca.png', 'ld.png', 'diagnostic.png'):
            p = os.path.join(self.out, f); self.assertTrue(os.path.isfile(p) and os.path.getsize(p) > 50, f)

    def test_heterozygosity_inside_and_outside_the_block(self):
        st = self.st; m, f = st[st.sex == 'M'], st[st.sex == 'F']
        self.assertGreater(m.region_het.min(), 0.9); self.assertLess(f.region_het.max(), 0.05)             # males X/Y, females XX
        self.assertAlmostEqual(m.genome_het_outside_region.mean(), f.genome_het_outside_region.mean(), delta=0.05)    # the same elsewhere

    def test_the_inbred_female_has_a_run_outside_the_block_and_the_others_do_not(self):
        st = self.st.set_index('sample')
        self.assertGreater(st.loc['F01', 'F_roh_outside_region'], 0.15)                                      # c1 0-1 Mb of the 3.5 Mb outside... about a fifth
        others = st.drop(index='F01'); self.assertLess(others.F_roh_outside_region.max(), 0.02)
        self.assertTrue((st.loc[[n for n in st.index if n.startswith('F')], 'region_in_roh_share'] == 1.0).all())   # every female is homozygous over the block: it is a run in all of them
        self.assertTrue((st.loc[[n for n in st.index if n.startswith('M')], 'region_in_roh_share'] == 0.0).all())
        seg = __import__('pandas').read_csv(os.path.join(self.out, 'roh_segments.tsv'), sep='\t')
        f01 = seg[seg['sample'] == 'F01']; self.assertEqual((int(f01.start.min()), int(f01.end.max())), (0, 3_000_000))     # 0-1 Mb of its own, joined to the block

    def test_relatedness_does_not_line_up_with_sex_in_unrelated_fish(self):
        self.assertGreater(self.r['ksex']['p'], 0.01); self.assertTrue(all(len(c) == 1 for c in self.r['clusters']))

    def test_the_block_separates_the_sexes_and_the_genome_does_not(self):
        auc, p, pc, used = self.r['res']['pca_region']
        self.assertEqual(auc, 1.0); self.assertLess(p, 0.01)
        self.assertLess(self.r['res']['pca_genome'][0], 0.95)

    def test_the_region_is_half_open_start_inclusive_end_exclusive(self):
        pos = [int(l.split('\t')[1]) for l in open(self.vcf) if l.startswith('c1\t')]
        self.assertIn(1_000_000, pos); self.assertIn(3_000_000, pos)                                          # SNPs sit exactly on both boundaries
        self.assertEqual(self.r['res']['n_region_sites'], sum(1 for q in pos if 1_000_000 <= q < 3_000_000))

    def test_diagnostic_counts(self):
        n_block = sum(1 for l in open(self.vcf) if l.startswith('c1\t') and 1_000_000 <= int(l.split('\t')[1]) < 3_000_000)
        self.assertGreater(self.r['diag_loose'], 0.9 * n_block); self.assertLessEqual(self.r['diag_loose'], n_block); self.assertEqual(self.r['diag_strict'], self.r['diag_loose'] if self.r['diag_strict'] == self.r['diag_loose'] else -1)
        self.assertEqual(sorted(self.r['by_depth']), [3, 5, 8])

    def test_ld_is_complete_in_the_block_and_low_elsewhere(self):
        import pandas as pd
        t = pd.read_csv(os.path.join(self.out, 'ld_decay.tsv'), sep='\t')
        reg = t[(t.set == 'region') & (t.group == 'all fish')].iloc[0]; rest = t[(t.set == 'rest of the chromosome') & (t.group == 'all fish')].iloc[0]
        self.assertGreater(reg['r2_1000-2000kb'], 0.9)                                                    # the block: every pair of SNPs is in complete LD (sex alone explains the genotypes)
        self.assertLess(rest['r2_2000-5000kb'], 0.2)                                                       # the flanks are two 1 Mb pieces, so the first pairs far apart are 2-5 Mb apart
        self.assertTrue(np.isnan(rest['r2_1000-2000kb']))                                                  # ... and the 1-2 Mb bin is empty, which the table must show as empty, not as 0
        none = t[(t.set == 'region') & (t.group != 'all fish')]                                            # within the block males are all heterozygous and females all alike: nothing varies
        self.assertEqual(sorted(none.group), ['females', 'males']); self.assertTrue((none.n_sites == 0).all()); self.assertTrue(none[[c for c in t.columns if c.startswith('r2_')]].isna().all().all())

    def test_chromosome_table_and_summary_text(self):
        import pandas as pd
        ch = pd.read_csv(os.path.join(self.out, 'chromosome_sex_het.tsv'), sep='\t').set_index('chrom')
        self.assertGreater(ch.loc['c1', 'male_minus_female'], 0.1); self.assertLess(abs(ch.loc['c2', 'male_minus_female']), 0.05)
        s = self.r['summary']
        for frag in ('# Region block analysis: c1:1,000,000-3,000,000', '12 males and 12 females', 'F_ROH outside the region', 'permutation p', 'not comparable between runs of different size', '## 6. Parts of the region'):
            self.assertIn(frag, s)

    def test_too_few_fish_is_an_error_not_a_crash(self):
        d = tempfile.mkdtemp()
        try:
            q = os.path.join(d, 'q.tsv'); open(q, 'w').write('sample\tsex\nM01\tM\nM02\tM\nF01\tF\n')
            self.assertEqual(CLI.main(['--vcf', self.vcf, '--qc', q, '--region', 'c1:1000000-3000000', '--outdir', os.path.join(d, 'o'), '--backend', 'builtin']), 2)
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_the_command_is_reachable_from_the_module_entry_point(self):
        from sd_regions import cli as top
        d = tempfile.mkdtemp()
        try:
            self.assertEqual(top.main(['region_block', '--vcf', self.vcf, '--qc', self.qc, '--region', 'c1:1000000-3000000', '--outdir', os.path.join(d, 'o'), '--backend', 'builtin', '--n-perm', '50', '--n-draws', '5']), 0)
            self.assertTrue(os.path.isfile(os.path.join(d, 'o', 'summary.md')))
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
