"""Tests of `bc_analysis analyze` (module 07): the functions on exact small examples and on a synthetic cross with planted truths, and the commands end to end."""
import io
import json
import os
import shutil
import sys
import tempfile
import unittest
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from bc_analysis.analyze import cli as AC                                    # noqa: E402
from bc_analysis.analyze import crossovers as XO                             # noqa: E402
from bc_analysis.analyze import data as D                                    # noqa: E402
from bc_analysis.analyze import distortion as TD                             # noqa: E402
from bc_analysis.analyze import recomb as R                                  # noqa: E402
from bc_analysis.analyze import sexlocus as SL                               # noqa: E402

GENOME = 'Mconophoros_GT1_NCBI'
CH = {'NC_140526.1': 300, 'NC_140527.1': 200, 'NC_140528.1': 250, 'NC_140535.1': 370}             # LG1, LG2, LG3, LG10 (bins of 100 kb)
SD_CHROM, SD_BIN = 'NC_140535.1', 275
HOT = (270, 300)                                                                                    # the planted hotspot on LG10, in bins; the sex determiner is inside the 27.5 Mb bin next to it


def setUpModule():
    from bc_analysis.genotyper import cli as G, sites as S
    G.set_lg_names(S.load_config(), GENOME)                                                          # the labels of Mconophoros_GT1_NCBI come from the config


def lg_label(c):
    from bc_analysis.genotyper import cli as G
    return G.lg_label(c)


def tiny_calls(spec, bins=10, chrom='NC_140526.1', bin_bp=100_000):
    """Calls of one chromosome from {fish: [(start state, [bins after which the state changes])]}: exact, for hand-checked results."""
    fish = list(spec)
    X = np.zeros((bins, len(fish)), dtype=np.int8)
    for j, f in enumerate(fish):
        s0, sw = spec[f]
        s = s0
        for b in range(bins):
            X[b, j] = s
            if b in sw:
                s = 1 - s
    meta = pd.DataFrame({'chrom': chrom, 'LG': lg_label(chrom), 'lg': 1, 'pos': np.arange(bins) * bin_bp})
    return D.Calls(meta, X, fish, bin_bp)


def make_cross(n_fish=240, seed=3, p_distorted=0.7, n_unsexed=24):
    """(Calls, fish table, truth). Crossovers: 0.17 cM per 100 kb, 5x that in the hotspot; LG2 is inherited with p_distorted for the YH state; sex = the state at the determiner bin, with two fish whose sex is wrong."""
    rng = np.random.default_rng(seed)
    fish = [f'F{i:03d}' for i in range(n_fish)]
    rows, mats = [], []
    states_by = {}
    for c, nb in CH.items():
        rate = np.full(nb - 1, 0.0017)
        if c == SD_CHROM:
            rate[HOT[0]:HOT[1]] = 0.0085
        X = np.zeros((nb, n_fish), dtype=np.int8)
        for j in range(n_fish):
            s0 = int(rng.random() < (p_distorted if c == 'NC_140527.1' else 0.5))
            sw = np.zeros(nb, dtype=int)
            for b in np.where(rng.random(nb - 1) < rate)[0]:
                sw[b + 1:] ^= 1
            X[:, j] = np.bitwise_xor(s0, sw)
        mats.append(X)
        states_by[c] = X
        rows.append(pd.DataFrame({'chrom': c, 'pos': np.arange(nb) * 100_000}))
    meta = pd.concat(rows, ignore_index=True)
    meta['LG'] = [lg_label(c) for c in meta.chrom]
    meta['lg'] = [int(l[2:]) for l in meta.LG]
    X = np.vstack(mats)
    order = np.lexsort((meta.pos.to_numpy(), meta.lg.to_numpy()))
    calls = D.Calls(meta.iloc[order].reset_index(drop=True), X[order], fish, 100_000)
    st = states_by[SD_CHROM][SD_BIN]
    sex = np.where(st == 1, 'M', 'F').astype(object)
    wrong = [3, 11]                                                                                    # two fish whose recorded sex is wrong
    for w in wrong:
        sex[w] = 'F' if sex[w] == 'M' else 'M'
    unsexed = rng.choice([i for i in range(n_fish) if i not in wrong], n_unsexed, replace=False)
    sex[unsexed] = 'U'
    ft = pd.DataFrame({'SampleID': fish, 'Sex': sex, 'Sire': ['S1' if i < n_fish // 4 else 'S2' for i in range(n_fish)], 'Dam': [f'D{i % 3}' for i in range(n_fish)], 'BroodID': [f'B{i % 5}' for i in range(n_fish)]})
    return calls, D.fish_table(ft), dict(wrong=[fish[w] for w in wrong], states=states_by)


def write_counts(folder, calls, fish_list, chrom, coverage, seed=5, sites_per_bin=24, eps=0.02, q=0.45):
    """Reads at informative sites drawn from the TRUE states of the calls: counts files as `genotype counts` writes them."""
    os.makedirs(folder, exist_ok=True)
    rng = np.random.default_rng(seed)
    ix = calls.idx(chrom)
    X = calls.X[ix]
    for f in fish_list:
        j = calls.fish.index(f)
        rows = []
        for b in range(X.shape[0]):
            pos = np.sort(rng.integers(b * 100_000, (b + 1) * 100_000, sites_per_bin))
            n = rng.poisson(coverage, len(pos))
            y = rng.binomial(n, q if X[b, j] == 1 else eps)
            rows += [(chrom, int(p), int(a - c), int(c), 1) for p, a, c in zip(pos, n, y) if a > 0]
        pd.DataFrame(rows, columns=['chrom', 'pos', 'n_mc', 'n_yh', 'yh_is_ref']).to_csv(os.path.join(folder, f'{f}.counts.tsv.gz'), sep='\t', index=False)


def write_run(root, calls, ft, run='hmm_mask_int2_deep'):
    os.makedirs(f'{root}{run}', exist_ok=True)
    d = pd.DataFrame(calls.X, columns=calls.fish)
    d.insert(0, 'Position', calls.meta.pos.to_numpy())
    d.insert(0, 'Chromosome', calls.meta.chrom.to_numpy())
    d.to_csv(f'{root}{run}/genotypes_calls.csv')
    ft.reset_index().rename(columns={'index': 'SampleID'}).to_csv(f'{root}fish.tsv', sep='\t', index=False)


class Base(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.calls, cls.ft, cls.truth = make_cross()
        cls.sex = D.sex_vector(cls.calls, cls.ft)
        cls.steps = R.Steps(cls.calls)


# ------------------------------------------------------------------ data
class TestData(unittest.TestCase):
    def test_load_orders_by_linkage_group_and_reads_the_labels(self):
        calls, ft, _ = make_cross(n_fish=20, n_unsexed=2)
        d = tempfile.mkdtemp()
        try:
            write_run(d + '/', calls, ft)
            c2 = D.Calls.load(d + '/hmm_mask_int2_deep/genotypes_calls.csv', lg_label)
            self.assertEqual([c2.lg_of(c) for c in c2.chroms], ['LG1', 'LG2', 'LG3', 'LG10'])
            self.assertEqual(c2.chrom_of(10), 'NC_140535.1'); self.assertEqual(c2.chrom_of('LG2'), 'NC_140527.1'); self.assertEqual(c2.length('NC_140535.1'), 37_000_000)
            self.assertTrue((c2.X == calls.X).all())
            with self.assertRaises(ValueError): c2.chrom_of(99)
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_probabilities_are_refused(self):
        d = tempfile.mkdtemp()
        try:
            pd.DataFrame({'Chromosome': ['a', 'a'], 'Position': [0, 100000], 'f1': [0.3, 0.9]}).to_csv(d + '/c.csv')
            with self.assertRaises(ValueError): D.Calls.load(d + '/c.csv')
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_sex_only_m_and_f_count(self):
        ft = D.fish_table(pd.DataFrame({'SampleID': ['a', 'b', 'c', 'd', 'e'], 'Sex': ['M', 'F', 'M?', 'U', '']}))
        calls = tiny_calls({k: (0, []) for k in 'abcde'})
        v = D.sex_vector(calls, ft)
        self.assertEqual(list(v[:2]), [1.0, 0.0]); self.assertTrue(np.isnan(v[2:]).all())


# ------------------------------------------------------------------ marey and rates
class TestRecomb(Base):
    def test_exact_map_of_a_tiny_example(self):
        calls = tiny_calls({'a': (0, [4]), 'b': (1, [2, 6]), 'c': (0, []), 'd': (1, [4])})          # 1 + 2 + 0 + 1 = 4 crossovers in 4 fish
        s = R.Steps(calls)
        mb = R.marey_bins(s)
        self.assertAlmostEqual(mb.cM.iloc[-1], 100.0)                                                 # 4 crossovers / 4 fish = 1 per fish = 100 cM
        lens = R.map_lengths(s, terminal_bp=200_000)
        self.assertEqual((lens.cM.iloc[0], lens.Mb.iloc[0]), (100.0, 1.0))
        self.assertEqual(lens.cM_without_tips.iloc[0], 100.0)                                         # none of the switches is within 2 bins of an end

    def test_a_crossover_belongs_to_the_window_that_contains_its_boundary(self):
        s = R.Steps(tiny_calls({'a': (0, [4])}, bins=10))                                              # the state changes between bins 4 and 5: at 500 kb
        self.assertEqual((s.counts('NC_140526.1', 500_000, 600_000)[0], s.counts('NC_140526.1', 400_000, 500_000)[0]), (1.0, 0.0))

    def test_the_region_ratio_is_exact_and_the_background_excludes_the_region(self):
        calls = tiny_calls({'a': (0, [12, 25]), 'b': (0, [15]), 'c': (0, [30]), 'd': (0, [5, 14])}, bins=40)
        s = R.Steps(calls)
        r = R.region_test(s, 'NC_140526.1', 1_000_000, 2_000_000, terminal_bp=0, n_boot=200)
        self.assertAlmostEqual(r['rate_cM_per_Mb'], 75.0); self.assertAlmostEqual(r['background_cM_per_Mb'], 100 * 0.75 / 2.9, places=6); self.assertAlmostEqual(r['ratio'], 75 / (100 * 0.75 / 2.9), places=6)
        self.assertEqual((r['crossovers'], r['fish_with_a_crossover']), (3, 3))

    def test_the_genome_wide_p_counts_the_observed_window_once(self):
        null = R.window_rates(self.steps, 2_500_000, 500_000, 2_000_000, skip={'LG3', SD_CHROM})
        r = R.region_test(self.steps, SD_CHROM, 27_000_000, 29_500_000, n_boot=100, null_windows=null)
        self.assertAlmostEqual(r['p_genome'], (1 + r['null_as_hot_or_hotter']) / (1 + r['null_windows']))

    def test_tips_and_masked_windows_are_left_out(self):
        calls = tiny_calls({'a': (0, [0]), 'b': (0, [10]), 'c': (0, [4])}, bins=12)                   # crossovers at the first and the last boundary (tips) and one inside
        s = R.Steps(calls)
        self.assertAlmostEqual(R.map_lengths(s, terminal_bp=200_000).cM_without_tips.iloc[0], 100 / 3, places=1)
        w = R.window_rates(s, width_bp=300_000, step_bp=100_000, terminal_bp=200_000)
        self.assertGreater(len(w), 0); self.assertTrue((w.start_Mb >= 0.2 - 1e-9).all()); self.assertTrue((w.start_Mb + 0.3 <= 1.2 - 0.2 + 1e-9).all())
        wm = R.window_rates(s, 300_000, 100_000, 200_000, masked={('NC_140526.1', 4)})
        self.assertLess(len(wm), len(w))

    def test_a_planted_hotspot_is_significant_and_a_control_region_is_not(self):
        steps = self.steps
        null = R.window_rates(steps, 2_500_000, 500_000, 2_000_000, skip={'LG3', SD_CHROM})
        hot = R.region_test(steps, SD_CHROM, 27_000_000, 29_500_000, n_boot=2000, null_windows=null)
        self.assertGreater(hot['ratio'], 2.5); self.assertGreater(hot['ci_lo'], 1.5); self.assertLess(hot['p_bootstrap_ratio_le_1'], 0.01)
        self.assertLessEqual(hot['null_as_hot_or_hotter'], 2); self.assertLess(hot['p_genome'], 0.2)
        cold = R.region_test(steps, SD_CHROM, 8_000_000, 10_500_000, n_boot=2000)
        self.assertLess(cold['ratio'], 1.8); self.assertLess(cold['ci_lo'], 1.0)

    def test_the_hottest_window_of_each_chromosome_is_ranked_and_the_hotspot_wins(self):
        best = R.best_windows(self.steps, 2_500_000, 500_000, 2_000_000)
        self.assertEqual(best.iloc[0].chrom, SD_CHROM); self.assertEqual(list(best['rank']), list(range(1, len(best) + 1)))
        self.assertTrue(26.0 <= best.iloc[0].start_Mb <= 30.0)

    def test_by_group_and_profile(self):
        t = R.by_group(self.steps, SD_CHROM, 27_000_000, 29_500_000, D.groups_of(self.calls, self.ft, 'Sire'))
        self.assertEqual(set(t.group), {'S1', 'S2'}); self.assertEqual(int(t.fish.sum()), 240)
        p = R.rate_profile(self.steps, SD_CHROM, 1_000_000, 300)
        self.assertEqual(len(p), 37); self.assertTrue((p.ci_lo <= p.rate + 1e-9).all() and (p.rate <= p.ci_hi + 1e-9).all())
        self.assertGreater(p[(p.start_Mb >= 27) & (p.start_Mb < 29)].rate.mean(), 3 * p[p.start_Mb < 20].rate.mean())


# ------------------------------------------------------------------ sex locus
class TestSexLocus(Base):
    def test_the_scan_finds_the_planted_determiner(self):
        scan = SL.sexscan(self.calls, self.sex)
        summ = SL.scan_summary(scan)
        self.assertEqual(summ.iloc[0].chrom, SD_CHROM); self.assertLessEqual(abs(summ.iloc[0].best_Mb - SD_BIN / 10), 0.5)
        self.assertGreater(summ.iloc[0].nlogp, 30); self.assertLess(summ.iloc[1].nlogp, 6)

    def test_exact_limits_from_a_hand_made_example(self):
        # the determiner is at bin 5; males (sex 1) carry state 1 at bin 5. Crossovers: m1 after bin 3 (state 0 then 1: SD right of 4), f1 after bin 7 (state 0 then 1: SD left of 8),
        # m2 after bin 4 (1 then 0 at 5? no: state at 5 must be 1) -> use fish with their own bounds
        spec = {'m_far': (1, []), 'f_far': (0, []), 'm_left': (0, [3]), 'f_right': (0, [7]), 'm_left2': (0, [1]), 'f_right2': (0, [8])}
        calls = tiny_calls(spec, bins=12)
        sex = np.array([1, 0, 1, 0, 1, 0], dtype=float)
        r = SL.finemap(calls, sex, 'NC_140526.1', locus_bp=500_000, window_bp=1_000_000)
        lim = r['limits']
        self.assertEqual(int(lim.left_bp.iloc[0]), 400_000); self.assertEqual(int(lim.right_bp.iloc[0]), 800_000)        # the nearest crossovers: 4 | 5 and 7 | 8
        self.assertEqual((int(lim.left_bp.iloc[1]), int(lim.right_bp.iloc[1])), (200_000, 900_000))                            # one nearest fish wrong on each side
        self.assertEqual(r['discordant'], [])

    def test_a_fish_with_the_wrong_sex_is_ignored_and_reported(self):
        r = SL.finemap(self.calls, self.sex, SD_CHROM)
        self.assertTrue(set(self.truth['wrong']) <= set(r['discordant']))
        lim = r['limits'].iloc[0]
        self.assertLessEqual(lim.left_bp, SD_BIN * 100_000); self.assertGreaterEqual(lim.right_bp, (SD_BIN + 1) * 100_000)
        self.assertTrue(set(r['discordant']) <= set(r['excluded'])); self.assertFalse(set(r['recombinants'].fish) & set(r['discordant']))

    def test_a_fish_with_several_crossovers_on_one_side_still_gives_its_nearest_bound(self):
        # the determiner at bin 5 (males carry 1). g crosses over at 7|8 and again at 9|10: its nearest crossover on the right is 7|8, not the later one
        spec = {'g': (0, [7, 9]), 'h': (1, [2, 3]), 'm': (1, []), 'f': (0, [])}
        calls = tiny_calls(spec, bins=12)
        sex = np.array([0, 1, 1, 0], dtype=float)
        r = SL.finemap(calls, sex, 'NC_140526.1', locus_bp=500_000, window_bp=1_000_000)
        rec = r['recombinants'].set_index('fish')
        self.assertEqual(int(rec.loc['g'].boundary_bp), 800_000); self.assertEqual(rec.loc['g'].bound, 'SD left of')
        self.assertEqual(int(rec.loc['h'].boundary_bp), 400_000); self.assertEqual(rec.loc['h'].bound, 'SD right of')                # h: states 1 1 1 0 1 1 ...: its nearest crossover left of the locus is 3|4, the earlier 2|3 does not matter
        self.assertEqual(r['discordant'], [])

    def test_a_crossover_beyond_the_window_is_not_used(self):
        spec = {'m': (0, [0]), 'f': (0, [15]), 'f2': (0, [18]), 'm0': (1, []), 'f0': (0, [])}          # the locus at bin 10 (1.0 Mb): m's crossover is at 0.1 Mb (900 kb away), f's at 1.6 Mb (600 kb), f2's at 1.9 Mb (900 kb)
        calls = tiny_calls(spec, bins=20)
        sex = np.array([1, 0, 0, 1, 0], dtype=float)
        near = SL.finemap(calls, sex, 'NC_140526.1', locus_bp=1_000_000, window_bp=800_000)['recombinants']
        self.assertEqual(list(near.fish), ['f']); self.assertEqual(near.bound.iloc[0], 'SD left of')
        wide = SL.finemap(calls, sex, 'NC_140526.1', locus_bp=1_000_000, window_bp=1_000_000)['recombinants']
        self.assertEqual(sorted(wide.fish), ['f', 'f2', 'm'])

    def test_an_excluded_fish_gives_no_limit_on_either_side(self):
        spec = {'m1': (0, [3]), 'f1': (0, [7]), 'm_far': (1, []), 'f_far': (0, [])}
        calls = tiny_calls(spec, bins=12); sex = np.array([1, 0, 1, 0], dtype=float)
        full = SL.finemap(calls, sex, 'NC_140526.1', locus_bp=500_000, window_bp=1_000_000)
        self.assertEqual(len(full['limits']), 1)
        for gone in ('m1', 'f1'):
            r = SL.finemap(calls, sex, 'NC_140526.1', locus_bp=500_000, window_bp=1_000_000, exclude=[gone])
            self.assertEqual(len(r['limits']), 0, gone); self.assertFalse(r['recombinants'][r['recombinants'].fish == gone].used.any())
            self.assertIn(gone, r['excluded'])

    def test_the_default_locus_is_the_middle_of_the_run_of_fewest_discordant_fish(self):
        spec = {'m_left': (0, [3]), 'f_right': (0, [7]), 'm_far': (1, []), 'f_far': (0, [])}          # no fish disagrees with its sex in bins 4-7
        calls = tiny_calls(spec, bins=12)
        r = SL.finemap(calls, np.array([1, 0, 1, 0], dtype=float), 'NC_140526.1', window_bp=1_000_000)
        self.assertEqual(r['locus_bp'], 600_000); self.assertEqual((int(r['limits'].left_bp.iloc[0]), int(r['limits'].right_bp.iloc[0])), (400_000, 800_000))

    def test_a_fish_that_disagrees_at_the_locus_is_excluded_and_does_not_move_the_limits(self):
        spec = {'m1': (0, [3]), 'f1': (0, [7]), 'x': (1, [8])}                  # x is recorded female but carries the 1 state at the locus
        calls = tiny_calls(spec, bins=12)
        r = SL.finemap(calls, np.array([1, 0, 0], dtype=float), 'NC_140526.1', locus_bp=500_000, window_bp=1_000_000)
        self.assertEqual(r['discordant'], ['x']); self.assertEqual((int(r['limits'].left_bp.iloc[0]), int(r['limits'].right_bp.iloc[0])), (400_000, 800_000))
        r2 = SL.finemap(calls, np.array([1, 0, 0], dtype=float), 'NC_140526.1', locus_bp=500_000, window_bp=1_000_000, exclude=['f1'])
        self.assertEqual(len(r2['limits']), 0)                                  # without f1 there is no crossover on the right


# ------------------------------------------------------------------ crossovers from reads
class TestCrossovers(Base):
    def test_a_perfect_split_is_found_between_the_right_sites(self):
        pos = np.arange(1000, 101000, 1000)
        n = np.full(len(pos), 3)
        y = np.where(pos < 50_000, 0, 1) * 1                                   # all YH reads at one in 3 right of 50 kb
        loc = XO.localize(pos, n, y)
        self.assertEqual(loc['state'], '0->1'); self.assertTrue(49_000 <= loc['best'] <= 51_000); self.assertLessEqual(loc['span1'][1] - loc['span1'][0], 4_000); self.assertFalse(loc['edge'])

    def test_planted_crossovers_are_localized_and_more_reads_narrow_them(self):
        res = SL.finemap(self.calls, self.sex, SD_CHROM)
        roles = XO.roles_from_finemap(res['recombinants'], 2)
        self.assertEqual({'left', 'right'}, set(roles.values()))
        d = tempfile.mkdtemp()
        try:
            write_counts(d + '/after', self.calls, list(roles), SD_CHROM, coverage=2.0, seed=1); write_counts(d + '/before', self.calls, list(roles), SD_CHROM, coverage=0.3, seed=2)
            lo, hi = res['locus_bp'] - 1_500_000, res['locus_bp'] + 1_500_000
            truth = {}
            ix = self.calls.idx(SD_CHROM); pos = self.calls.pos(SD_CHROM)
            for f in roles:
                s = self.calls.X[ix][:, self.calls.fish.index(f)]
                bs = [int(pos[i + 1]) for i in np.where(s[:-1] != s[1:])[0] if lo <= pos[i + 1] < hi]
                truth[f] = bs
            widths = {}
            for tag in ('after', 'before'):
                out = {}
                for f in roles:
                    p, n, y = XO.read_sites(f'{d}/{tag}', f, SD_CHROM, lo, hi)
                    out[f] = XO.localize(p, n, y)
                    self.assertTrue(any(abs(out[f]['best'] - b) <= 30_000 for b in truth[f]), (tag, f, out[f]['best'], truth[f]))
                widths[tag] = np.median([o['span1'][1] - o['span1'][0] for o in out.values()])
                best, cons = XO.interval(out, roles)
                self.assertLessEqual(best[0], best[1]) if tag == 'after' else None
            self.assertLessEqual(widths['after'], widths['before'])
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_interval_and_sensitivity(self):
        res = {'a': dict(best=100.0, span1=(90, 110)), 'b': dict(best=200.0, span1=(190, 215)), 'c': dict(best=50.0, span1=(40, 60)), 'd': dict(best=300.0, span1=(280, 330))}
        roles = {'a': 'left', 'c': 'left', 'b': 'right', 'd': 'right'}
        best, cons = XO.interval(res, roles)
        self.assertEqual(best, (100.0, 200.0)); self.assertEqual(cons, (90, 215))
        with self.assertRaises(ValueError): XO.interval({'a': res['a']}, {'a': 'left'})
        pos = np.arange(1000, 101000, 1000); n = np.full(len(pos), 4); y = np.where(pos < 40_000, 0, 2)
        sw = XO.sweep(pos, n, y); self.assertLess(max(sw.values()) - min(sw.values()), 2000)

    def test_the_two_lod_span_contains_and_exceeds_the_one_lod_span_on_noisy_reads(self):
        rng = np.random.default_rng(4)
        pos = np.arange(1000, 61000, 1000); n = np.ones(len(pos), dtype=int); y = rng.binomial(1, np.where(pos < 30_000, 0.02, 0.45))
        loc = XO.localize(pos, n, y)
        self.assertTrue(loc['span2'][0] <= loc['span1'][0] and loc['span1'][1] <= loc['span2'][1]); self.assertGreater(loc['span2'][1] - loc['span2'][0], loc['span1'][1] - loc['span1'][0])

    def test_a_split_at_the_edge_of_the_data_is_flagged(self):
        pos = np.arange(1000, 31000, 1000); n = np.full(len(pos), 4); y = np.full(len(pos), 2)          # heterozygous all the way: no crossover in the window
        self.assertTrue(XO.localize(pos, n, y)['edge'])
        y2 = np.where(pos < 15_000, 0, 2); self.assertFalse(XO.localize(pos, n, y2)['edge'])

    def test_too_few_sites_is_an_error(self):
        with self.assertRaises(ValueError): XO.localize(np.array([1]), np.array([2]), np.array([1]))


# ------------------------------------------------------------------ distortion and linkage
class TestDistortion(Base):
    def test_the_planted_distortion_is_the_only_significant_region(self):
        sc = TD.trd_scan(self.calls, n_perm=300, seed=1)
        reg = TD.trd_regions(self.calls, sc, 3.0)
        self.assertGreater(len(reg), 0); top = reg.iloc[0]
        self.assertEqual(top.LG, 'LG2'); self.assertGreater(top.mean_share, 0.6); self.assertLess(top.p_genomewide, 0.01)
        self.assertEqual(set(reg[reg.p_genomewide < 0.05].LG), {'LG2'})

    def test_without_distortion_nothing_is_significant(self):
        calls, _, _ = make_cross(p_distorted=0.5, seed=11)
        sc = TD.trd_scan(calls, n_perm=300, seed=3)
        reg = TD.trd_regions(calls, sc, 3.0)
        self.assertEqual(len(reg[reg.p_genomewide < 0.05]), 0)

    def test_runs_of_bins_are_split_by_the_gap_setting(self):
        calls = tiny_calls({f'f{i}': (0, []) for i in range(4)}, bins=30)
        z = np.zeros(30); z[[3, 4, 5, 9, 10]] = 4.0
        sc = dict(z=z, f=0.5 + z / 40, G=np.array([10.0]))
        self.assertEqual(len(TD.trd_regions(calls, sc, 3.0, gap=3)), 1)                               # 5 -> 9 skips 3 bins: one run
        self.assertEqual(len(TD.trd_regions(calls, sc, 3.0, gap=2)), 2)

    def test_the_sign_flip_null_has_the_right_scale(self):
        sc = TD.trd_scan(self.calls, n_perm=300, seed=2)
        self.assertTrue(2.0 < np.percentile(sc['G'], 50) < 4.2)                              # the median maximum |z| over ~1100 correlated bins
        self.assertAlmostEqual(TD.p_genomewide(sc, 1e9), 1 / 301)

    def test_group_checks(self):
        st = np.array([1] * 30 + [0] * 10 + [1] * 10 + [0] * 30)
        g = pd.Series(['a'] * 40 + ['b'] * 40)
        t = TD.group_checks(st, g, 'Sire')
        self.assertEqual(list(t.group), ['a', 'b']); self.assertAlmostEqual(t.share.iloc[0], 0.75); self.assertLess(t.heterogeneity_p.iloc[0], 1e-3)

    def test_linkscan_finds_the_partner_of_a_segment(self):
        rng = np.random.default_rng(1)
        n = 120
        a = rng.integers(0, 2, n)
        c1 = np.tile(a, (20, 1)); c2 = np.tile(rng.integers(0, 2, n), (20, 1)); c2[5:12] = a                          # chromosome 2 bins 5-11 are inherited with chromosome 1
        meta = pd.DataFrame({'chrom': ['NC_140526.1'] * 20 + ['NC_140527.1'] * 20, 'LG': ['LG1'] * 20 + ['LG2'] * 20, 'lg': [1] * 20 + [2] * 20, 'pos': list(np.arange(20) * 100_000) * 2})
        calls = D.Calls(meta, np.vstack([c1, c2]).astype(np.int8), [f'f{i}' for i in range(n)])
        best, own = TD.linkscan(calls, 'NC_140526.1', 0, 1_000_000, smooth=3)
        top = best[best.chrom != 'NC_140526.1'].iloc[0]
        self.assertEqual(top.LG, 'LG2'); self.assertGreater(top.link, 0.6); self.assertTrue(0.4 <= top.best_Mb <= 1.2)


# ------------------------------------------------------------------ the commands
class TestCommands(Base):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.tmp = tempfile.mkdtemp()
        cls.master, cls.cloud = cls.tmp + '/master/', cls.tmp + '/dropbox/'
        cls.root = f'{cls.master}Outputs/BCGenotypes/{GENOME}/'
        os.makedirs(cls.cloud)
        write_run(cls.root, cls.calls, cls.ft)
        res = SL.finemap(cls.calls, cls.sex, SD_CHROM)
        cls.roles = XO.roles_from_finemap(res['recombinants'], 2)
        write_counts(cls.root + 'counts', cls.calls, list(cls.roles), SD_CHROM, 2.0, seed=1)
        write_counts(cls.root + 'before', cls.calls, list(cls.roles), SD_CHROM, 0.3, seed=2)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def twin(self, path, master=None):
        """The Dropbox twin of a local path: the same replacement of the master folders that the FileManager and cichlid_common.runner.to_cloud make."""
        return path.replace(master or self.master, self.cloud)

    def fake_uploader(self, seen):
        def up(path):
            seen.append(path)
            dst = self.twin(path)
            shutil.copytree(path, dst, dirs_exist_ok=True)
            return dst
        return up

    def fake_downloader(self, seen, master):
        def down(path):
            seen.append(path)
            src = self.twin(path, master)
            if not os.path.exists(src):
                raise FileNotFoundError(src)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            shutil.copytree(src, path) if os.path.isdir(src) else shutil.copy2(src, path)
        return down

    def run_cmd(self, *argv, uploader=None, downloader=None, root=None):
        lines, err, old = [], io.StringIO(), sys.stderr
        sys.stderr = err
        try:
            rc = AC.main(['--genome', GENOME, '--out', root or self.root, '--fish-table', self.root + 'fish.tsv', '--no-plots'] + list(argv), out=lines.append, uploader=uploader, downloader=downloader)
        finally:
            sys.stderr = old
        return rc, '\n'.join(lines), err.getvalue()

    def test_marey(self):
        rc, out, err = self.run_cmd('marey'); self.assertEqual(rc, 0, err)
        d = self.root + 'analysis/marey/'
        lens = pd.read_csv(d + 'map_lengths.tsv', sep='\t')
        self.assertEqual(list(lens.LG), ['LG1', 'LG2', 'LG3', 'LG10']); self.assertTrue((lens.cM > 15).all() and (lens.cM < 110).all())
        self.assertTrue(os.path.exists(d + 'marey_bins.tsv') and os.path.exists(d + 'marey_windows.tsv'))
        prov = json.load(open(d + 'provenance.json')); self.assertEqual(prov['step'], 'marey'); self.assertEqual(prov['run'], 'hmm_mask_int2_deep'); self.assertEqual(len(prov['inputs']['calls']['sha1']), 40)
        self.assertIn('map_lengths.tsv', prov['written']); self.assertIn('mean map', out)

    def test_rate(self):
        rc, out, err = self.run_cmd('rate', '--lg', '10', '--region', '27.0-29.5'); self.assertEqual(rc, 0, err)
        d = self.root + 'analysis/rate/'
        r = pd.read_csv(d + 'rate_region.tsv', sep='\t').iloc[0]
        self.assertGreater(r.ratio, 2.5); self.assertEqual(r.LG, 'LG10'); self.assertEqual(int(r.rank_of_chromosome_best_window), 1); self.assertEqual(int(r.null_as_hot_or_hotter), 0, 'the same-width windows are from OTHER chromosomes')
        self.assertEqual(set(pd.read_csv(d + 'rate_by_group.tsv', sep='\t').by), {'Sire', 'BroodID', 'Dam'})
        self.assertIn('same-width windows', out)

    def test_sexscan_and_finemap(self):
        rc, out, err = self.run_cmd('sexscan'); self.assertEqual(rc, 0, err); self.assertIn('LG10', out)
        rc, out, err = self.run_cmd('finemap'); self.assertEqual(rc, 0, err)
        d = self.root + 'analysis/finemap/'
        lim = pd.read_csv(d + 'finemap_limits.tsv', sep='\t')
        self.assertLessEqual(lim.left_Mb.iloc[0], 27.5); self.assertGreaterEqual(lim.right_Mb.iloc[0], 27.6)
        ex = pd.read_csv(d + 'finemap_excluded.tsv', sep='\t'); self.assertTrue(set(self.truth['wrong']) <= set(ex.fish))
        self.assertIn('ignored fish', out)
        rc, out, err = self.run_cmd('finemap', '--lg', '10', '--exclude-fish', 'F100'); self.assertEqual(rc, 0, err); self.assertIn('F100', out)

    def test_crossovers_from_the_finemap_and_with_explicit_fish(self):
        rc, out, err = self.run_cmd('crossovers', '--from-finemap', '--counts', self.root + 'counts', '--counts-before', self.root + 'before'); self.assertEqual(rc, 0, err)
        d = self.root + 'analysis/crossovers/'
        iv = pd.read_csv(d + 'crossovers_interval.tsv', sep='\t'); self.assertEqual(list(iv.when), ['after', 'before'])
        t = pd.read_csv(d + 'crossovers.tsv', sep='\t'); self.assertEqual(len(t), 2 * len(self.roles)); self.assertTrue({'span1_lo', 'fitted_best_bp', 'sweep_min_bp'} <= set(t.columns))
        a = iv[iv.when == 'after'].iloc[0]; self.assertLessEqual(a.best_left_bp, SD_BIN * 100_000 + 100_000); self.assertGreaterEqual(a.best_right_bp, SD_BIN * 100_000)
        fish = ','.join(f'{f}:{r}' for f, r in self.roles.items())
        rc, out, err = self.run_cmd('crossovers', '--fish', fish, '--lg', '10', '--span', '26-29', '--counts', self.root + 'counts'); self.assertEqual(rc, 0, err)
        rc, out, err = self.run_cmd('crossovers', '--counts', self.root + 'counts'); self.assertEqual(rc, 2); self.assertIn('--from-finemap', err)

    def test_trd_and_linkscan(self):
        rc, out, err = self.run_cmd('trd'); self.assertEqual(rc, 0, err)
        d = self.root + 'analysis/trd/'
        reg = pd.read_csv(d + 'trd_regions.tsv', sep='\t'); self.assertEqual(reg.iloc[0].LG, 'LG2'); self.assertIn('genome-wide 5% threshold', out)
        self.assertTrue(os.path.exists(d + 'trd_groups.tsv')); self.assertIn('in_tip', reg.columns); self.assertTrue(bool(reg.iloc[0].in_tip), 'the planted distortion covers all of LG2, ends included'); self.assertIn('chromosome tip', out)
        rc, out, err = self.run_cmd('linkscan', '--lg', '2', '--region', '1-9'); self.assertEqual(rc, 0, err); self.assertIn('inherited together with', out)

    def test_upload_goes_through_the_verified_sync_to_the_twin_path(self):
        seen = []
        rc, out, err = self.run_cmd('marey', '--local-master', self.master, '--cloud-master', self.cloud, '--upload', uploader=self.fake_uploader(seen)); self.assertEqual(rc, 0, err)
        self.assertEqual(seen, [self.root + 'analysis/marey/'])
        twin = f'{self.cloud}Outputs/BCGenotypes/{GENOME}/analysis/marey/'
        self.assertTrue(os.path.exists(twin + 'map_lengths.tsv') and os.path.exists(twin + 'provenance.json')); self.assertIn('uploaded and verified', out)
        rc, out, err = self.run_cmd('marey', '--local-master', self.master, '--cloud-master', self.cloud, '--upload', uploader=lambda p: (_ for _ in ()).throw(RuntimeError('UPLOAD VERIFICATION FAILED for x')))
        self.assertEqual(rc, 2); self.assertIn('UPLOAD VERIFICATION FAILED', err)

    def test_publish_uploads_the_whole_genome_folder_verified(self):
        seen = []
        rc, out, err = self.run_cmd('publish', '--local-master', self.master, '--cloud-master', self.cloud, uploader=self.fake_uploader(seen)); self.assertEqual(rc, 0, err)
        self.assertEqual(seen, [self.root]); self.assertTrue(os.path.exists(f'{self.cloud}Outputs/BCGenotypes/{GENOME}/hmm_mask_int2_deep/genotypes_calls.csv')); self.assertIn('uploaded and verified', out)
        self.assertFalse(os.path.exists(self.root + 'analysis/publish'), 'publish writes no analysis folder')
        rc, out, err = self.run_cmd('publish', '--local-master', self.master, '--cloud-master', self.cloud, '--out', '/tmp/not_in_the_master/', uploader=self.fake_uploader([])); self.assertEqual(rc, 2); self.assertIn('no Dropbox twin', err)

    def test_a_folder_outside_the_local_master_cannot_be_uploaded(self):
        other = tempfile.mkdtemp() + '/'
        try:
            rc, out, err = self.run_cmd('marey', '--local-master', self.master, '--cloud-master', self.cloud, '--analysis-dir', other, '--upload', uploader=self.fake_uploader([]))
            self.assertEqual(rc, 2); self.assertIn('no Dropbox twin', err)
        finally:
            shutil.rmtree(other, ignore_errors=True)

    def test_missing_inputs_are_fetched_with_the_file_manager_and_never_without_permission(self):
        m2 = self.tmp + '/master2/'
        root2 = f'{m2}Outputs/BCGenotypes/{GENOME}/'
        for sub in ('hmm_mask_int2_deep', 'counts'):                                              # on Dropbox only: the calls and the counts of the bounding fish
            shutil.copytree(self.root + sub, f'{self.cloud}Outputs/BCGenotypes/{GENOME}/{sub}', dirs_exist_ok=True)
        os.makedirs(root2)
        shutil.copy2(self.root + 'fish.tsv', root2 + 'fish.tsv')
        seen = []
        rc, out, err = self.run_cmd('marey', '--local-master', m2, '--cloud-master', self.cloud, root=root2, downloader=self.fake_downloader(seen, m2)); self.assertEqual(rc, 0, err)
        self.assertIn(root2 + 'hmm_mask_int2_deep/genotypes_calls.csv', seen); self.assertTrue(os.path.exists(root2 + 'hmm_mask_int2_deep/genotypes_calls.csv'))
        self.assertIn('note: no mask', out)                                                         # the mask is optional: absent on Dropbox too, so none is applied (said, not silent)
        seen.clear()
        rc, out, err = self.run_cmd('crossovers', '--from-finemap', '--counts', root2 + 'counts', '--local-master', m2, '--cloud-master', self.cloud, root=root2, downloader=self.fake_downloader(seen, m2)); self.assertEqual(rc, 0, err)
        self.assertEqual(sorted(os.path.basename(p) for p in seen), sorted(f'{f}.counts.tsv.gz' for f in self.roles))
        shutil.rmtree(root2 + 'hmm_mask_int2_deep')
        rc, out, err = self.run_cmd('marey', '--local-master', m2, '--cloud-master', self.cloud, '--no-download', root=root2, downloader=self.fake_downloader([], m2)); self.assertEqual(rc, 2); self.assertIn('no calls at', err)
        rc, out, err = self.run_cmd('marey', '--local-master', m2, '--cloud-master', self.cloud, '--run', 'not_on_dropbox', root=root2, downloader=self.fake_downloader([], m2)); self.assertEqual(rc, 2); self.assertIn('no calls at', err); self.assertIn('could not be fetched from Dropbox either', err)
        rc, out, err = self.run_cmd('marey', '--calls', '/nowhere/c.csv', '--local-master', m2, '--cloud-master', self.cloud, root=root2, downloader=self.fake_downloader([], m2)); self.assertEqual(rc, 2); self.assertIn('no calls at', err)

    def test_clear_errors(self):
        rc, out, err = self.run_cmd('rate'); self.assertEqual(rc, 2); self.assertIn('--lg and --region', err)
        rc, out, err = self.run_cmd('rate', '--lg', '10', '--region', '30-20'); self.assertEqual(rc, 2); self.assertIn('START must be below END', err)
        rc, out, err = self.run_cmd('marey', '--run', 'nothing_here'); self.assertEqual(rc, 2); self.assertIn('no calls at', err)
        rc, out, err = self.run_cmd('rate', '--lg', '99', '--region', '1-2'); self.assertEqual(rc, 2); self.assertIn('no chromosome', err)

    def test_the_figures_are_drawn_when_matplotlib_is_there(self):
        try:
            import matplotlib  # noqa
        except ImportError:
            self.skipTest('no matplotlib')
        lines, err = [], io.StringIO()
        rc = AC.main(['marey', '--genome', GENOME, '--out', self.root, '--analysis-dir', self.root + 'figtest/'], out=lines.append)
        self.assertEqual(rc, 0); self.assertTrue(os.path.getsize(self.root + 'figtest/marey_maps.png') > 5000)
        rc = AC.main(['finemap', '--genome', GENOME, '--out', self.root, '--fish-table', self.root + 'fish.tsv', '--analysis-dir', self.root + 'figtest2/'], out=lines.append)
        self.assertEqual(rc, 0); self.assertTrue(os.path.getsize(self.root + 'figtest2/finemap.png') > 5000)

    def test_paths_with_doubled_slashes_reach_the_file_manager_clean(self):
        """A HOME that ends in a slash turns ~/x into //x on mzebra; the FileManager maps local to Dropbox paths by string replacement and cannot match that."""
        self.assertEqual(AC.clean_path('/a//b/./c/../d/'), '/a/b/d/'); self.assertEqual(AC.clean_path('/a//b'), '/a/b'); self.assertEqual(AC.clean_path(None), None); self.assertEqual(AC.clean_path(''), '')
        m2 = self.tmp + '/master3/'
        root2 = f'{m2}Outputs/BCGenotypes/{GENOME}/'
        shutil.copytree(self.root + 'counts', f'{self.cloud}Outputs/BCGenotypes/{GENOME}/counts', dirs_exist_ok=True); shutil.copytree(self.root + 'hmm_mask_int2_deep', f'{self.cloud}Outputs/BCGenotypes/{GENOME}/hmm_mask_int2_deep', dirs_exist_ok=True)
        os.makedirs(root2); shutil.copy2(self.root + 'fish.tsv', root2 + 'fish.tsv')
        seen = []
        slashy = lambda p: p.replace('/Outputs/', '//Outputs/')
        rc, out, err = self.run_cmd('crossovers', '--from-finemap', '--counts', slashy(root2) + 'counts', '--local-master', slashy(m2), '--cloud-master', self.cloud, root=slashy(root2), downloader=self.fake_downloader(seen, m2))
        self.assertEqual(rc, 0, err); self.assertTrue(seen)
        self.assertTrue(all('//' not in p for p in seen), seen)

    def test_a_missing_counts_folder_is_named(self):
        rc, out, err = self.run_cmd('crossovers', '--from-finemap', '--counts', self.root + 'counts', '--counts-before', self.root + 'never_copied', '--no-download'); self.assertEqual(rc, 2)
        self.assertIn('no counts file for', err); self.assertIn('the folder does not exist here', err)

    def test_the_group_is_reachable_from_the_top_level_command(self):
        from bc_analysis import cli as TOP
        rc = TOP.main(['analyze', 'marey', '--genome', GENOME, '--out', self.root, '--analysis-dir', self.root + 'top/', '--no-plots'])
        self.assertEqual(rc, 0); self.assertTrue(os.path.exists(self.root + 'top/map_lengths.tsv'))


if __name__ == '__main__':
    unittest.main()
