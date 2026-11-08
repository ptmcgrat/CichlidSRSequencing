"""B-chromosome panel measurement and call on synthetic BAM output (fake samtools): definition build, measurement, definition key, presence calls and flags, the scale validation, and the panel registry."""
import json, os, sys, tempfile, types, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import yaml
from cichlid_common.common_light import config_path
from sd_genotypes.bchrom import definition as D, depth, measure as M, score as S, validate_scale as V
from sd_genotypes.panel import measure as pm

CFG0 = yaml.safe_load(open(config_path('bchrompanel.yaml')))
LENGTHS = {'chrA': 6_000_000, 'chrB': 8_000_000, 'chrC': 7_000_000, 'scaf1': 100_000}


def make_core():
    rows = []
    for i in range(30):
        rows.append(('chrA', 100_000 + 150_000 * i, 'A' if i < 25 else 'B'))
    for i in range(30):
        rows.append(('chrB', 200_000 + 200_000 * i, 'A'))
    for s in (10_000, 40_000, 70_000):
        rows.append(('scaf1', s, 'A'))
    core = pd.DataFrame(rows, columns=['chrom', 'start', 'tier']); core['end'] = core.start + 1000
    core.insert(0, 'core_id', [f'core{i + 1:04d}' for i in range(len(core))]); core['LG'] = 'LGx'
    return core


def make_signal(core):
    sig = pd.DataFrame({'core_id': core.core_id})
    for n in ('N1', 'N2', 'N3'):
        sig[f'lab|RUN_{n}'] = 1.0
        sig.loc[:9, f'lab|RUN_{n}'] = 2.0                                                # the first 10 blocks are elevated in every B-free sample: not clean
    sig['carrier|RUN_C1'] = 5.0
    summ = pd.DataFrame(dict(run=['RUN_N1', 'RUN_N2', 'RUN_N3', 'RUN_C1'], role=['baseline', 'baseline', 'control', 'carrier']))
    return sig, summ


def make_def(n_bins=400):
    tmp = tempfile.mkdtemp(); core = make_core(); sig, summ = make_signal(core)
    D.build(core, LENGTHS, tmp, sig, summ, n_bins=n_bins)
    cfg = dict(CFG0); cfg['definition_dir'] = tmp
    cfg['thresholds'] = dict(cfg['thresholds'], min_scale_bins=int(0.8 * n_bins), min_blocks_used=20)       # the shipped limits are for the real definition (12,000 bins, 836 blocks)
    return cfg, tmp


def fake_run(cfg, depth_of, lengths=LENGTHS, log=None):
    def run(cmd):
        if log is not None:
            log.append(cmd)
        if cmd[:2] == ['samtools', '--version']:
            return types.SimpleNamespace(stdout='samtools 1.x\n')
        if cmd[:3] == ['samtools', 'view', '-H']:
            return types.SimpleNamespace(stdout='@HD\tVN:1.6\n' + ''.join(f'@SQ\tSN:{c}\tLN:{L}\n' for c, L in lengths.items()))
        if cmd[:2] == ['samtools', 'bedcov']:
            lines = []
            for line in open(cmd[-2]):
                c, s, e = line.split('\t'); s, e = int(s), int(e)
                lines.append(f'{c}\t{s}\t{e}\t{depth_of(c, s, e) * (e - s)}')
            return types.SimpleNamespace(stdout='\n'.join(lines) + '\n')
        raise AssertionError(cmd)
    return run


def depth_fn(blocks, base=10.0, fold=1.0, zero_every=50):
    iv = {c: list(zip(g.start, g.end)) for c, g in blocks.groupby('chrom')}

    def f(c, s, e):
        if any(s < be and e > bs for bs, be in iv.get(c, [])):
            return base * fold
        return 0.0 if (s // 50) % zero_every == 0 else base
    return f


class TestDefinition(unittest.TestCase):
    def test_blocks_scale_bins_and_clean_flags(self):
        cfg, tmp = make_def(); blocks, bins, meta = M.load_definition(cfg)
        self.assertEqual((len(blocks), blocks.tier.value_counts().to_dict()), (63, {'A': 58, 'B': 5}))
        self.assertEqual(blocks.clean.tolist()[:12], [False] * 10 + [True] * 2)           # elevated in every B-free sample = not clean
        self.assertEqual(sorted(bins.chrom.unique()), ['chrA', 'chrB', 'chrC'])           # scaf1 (100 kb) is not a placed chromosome
        self.assertTrue(abs(len(bins) - 400) <= 3 and (bins.end - bins.start == 50).all() and (bins.start % 50 == 0).all() and not bins.duplicated().any())
        for c, L in (('chrA', 6e6), ('chrB', 8e6), ('chrC', 7e6)):
            self.assertAlmostEqual((bins.chrom == c).sum(), 400 * L / 21e6, delta=3)      # proportional to length
        for r in bins.itertuples():
            self.assertFalse(((blocks.chrom == r.chrom) & (blocks.start < r.end) & (blocks.end > r.start)).any())   # no scale bin overlaps a block
        g = np.diff(bins[bins.chrom == 'chrB'].start.to_numpy()); self.assertLess(g.max(), 2 * g.mean())           # evenly spread
        self.assertEqual(meta['lengths'], LENGTHS)                                         # the sequences the blocks and bins use, with their lengths

    def test_bins_that_would_land_on_a_block_move_on(self):
        L, n = 6_000_000, 100; step = (L // 50) / n
        first = [int((i + 0.5) * step) * 50 for i in range(n)]                              # where the sampler puts its bins when nothing is in the way
        blocks = pd.DataFrame(dict(chrom='chrA', start=first, end=[f + 50 for f in first]))
        bins = D.sample_scale_bins({'chrA': L}, blocks, n_bins=n)
        self.assertEqual(len(bins), n)
        self.assertFalse(set(bins.start) & set(first))                                      # every bin that would have overlapped a block was moved
        self.assertTrue((bins.start.diff().dropna() > 0).all())

    def test_lengths_from_fai_and_track_json_and_missing_sequences(self):
        d = tempfile.mkdtemp(); open(d + '/g.fa.fai', 'w').write(''.join(f'{c}\t{L}\t0\t60\t61\n' for c, L in LENGTHS.items())); json.dump(dict(lengths=LENGTHS), open(d + '/t.json', 'w'))
        self.assertEqual(D.read_lengths(d + '/g.fa.fai'), LENGTHS); self.assertEqual(D.read_lengths(d + '/t.json'), LENGTHS)
        with self.assertRaises(ValueError):
            D.build(make_core(), {k: v for k, v in LENGTHS.items() if k != 'scaf1'}, d)


class TestMeasure(unittest.TestCase):
    def setUp(self):
        self.cfg, self.tmp = make_def(); self.blocks, self.bins, _ = M.load_definition(self.cfg)

    def test_depths_of_blocks_and_scale_bins(self):
        log = []
        out = M.measure_bam('x.bam', self.cfg, fake_run(self.cfg, depth_fn(self.blocks, 10, 4.0), log=log))
        self.assertEqual((out['n_blocks'], out['n_scale_bins'], set(out['block_depth']), out['definition']), (63, len(self.bins), {40.0}, M.definition_key(self.cfg)))
        sb = np.array(out['scale_bins']); self.assertTrue(set(np.unique(sb)) <= {0.0, 10.0}); self.assertGreater((sb == 0).sum(), 0)
        bc = [c for c in log if c[:2] == ['samtools', 'bedcov']]
        self.assertEqual(len(bc), 2); self.assertNotIn('-Q', bc[0])                         # no MAPQ filter by default: B reads mismap on purpose

    def test_a_mapq_filter_only_when_asked_for(self):
        self.assertEqual(M.bedcov_command('b.bam', 'r.bed'), ['samtools', 'bedcov', 'r.bed', 'b.bam'])
        self.assertEqual(M.bedcov_command('b.bam', 'r.bed', 20), ['samtools', 'bedcov', '-Q', '20', 'r.bed', 'b.bam'])

    def test_a_bam_from_another_reference_or_a_missing_definition_is_refused(self):
        with self.assertRaises(ValueError):
            M.measure_bam('x.bam', self.cfg, fake_run(self.cfg, depth_fn(self.blocks), lengths=dict(LENGTHS, chrB=7_999_999)))
        with self.assertRaises(ValueError):
            M.measure_bam('x.bam', self.cfg, fake_run(self.cfg, depth_fn(self.blocks), lengths={k: v for k, v in LENGTHS.items() if k != 'chrC'}))
        c = dict(self.cfg); c['definition_dir'] = tempfile.mkdtemp()
        with self.assertRaises(FileNotFoundError):
            M.measure_bam('x.bam', c, fake_run(self.cfg, depth_fn(self.blocks)))

    def test_the_definition_key_follows_what_is_measured(self):
        k = M.definition_key(self.cfg)
        c2 = dict(self.cfg); c2['min_mapq'] = 20; self.assertNotEqual(M.definition_key(c2), k)
        c3 = dict(self.cfg); c3['thresholds'] = dict(c3['thresholds'], present_min=9); self.assertEqual(M.definition_key(c3), k)       # a CALL threshold is not part of the definition
        cfg4, _ = make_def(n_bins=300); self.assertNotEqual(M.definition_key(cfg4), k)

    def test_the_panel_registry_runs_it_and_the_json_round_trips(self):
        self.assertTrue({'gsdf', 'lg5inv', 'bchrom'} <= set(pm.MODULES))
        for m in ('lg5inv', 'bchrom'):
            self.assertIsInstance(pm.module_cfg(m), dict)                                  # the shipped config of each registered module loads
        d = tempfile.mkdtemp(); bam = d + '/s.all.bam'; open(bam, 'wb').write(b'x'); open(bam + '.bai', 'wb').write(b'x')
        res = pm.measure(bam, 'S1', modules=('bchrom',), run=fake_run(self.cfg, depth_fn(self.blocks, 10, 4.0)), cfgs={'bchrom': self.cfg})
        pm.write_json(d + '/S1.panel.json', res); rows, skipped = pm.collect(d, 'bchrom')
        self.assertEqual((len(rows), skipped, rows[0]['SampleID'], rows[0]['definition']), (1, [], 'S1', M.definition_key(self.cfg)))
        with self.assertRaises(ValueError):
            pm.measure(bam, 'S1', modules=('nope',), run=fake_run(self.cfg, depth_fn(self.blocks)))


class TestScore(unittest.TestCase):
    def setUp(self):
        self.cfg, _ = make_def(); self.blocks, self.bins, _ = M.load_definition(self.cfg)

    def raw(self, fold, base=10.0, sample='S'):
        out = M.measure_bam('x.bam', self.cfg, fake_run(self.cfg, depth_fn(self.blocks, base, fold)))
        return dict(SampleID=sample, **out)

    def call(self, raw):
        t = S.call_table([raw], self.cfg).iloc[0]; return t

    def test_present_absent_uncertain(self):
        for fold, want in [(4.0, 'B+'), (1.0, 'B-'), (2.0, 'B?')]:
            t = self.call(self.raw(fold)); self.assertEqual((t.bchrom, t.call_flags, round(t.score, 3)), (want, '', fold))

    def test_the_scale_is_the_median_of_the_non_zero_bins(self):
        raw = self.raw(4.0); m = S.metrics(raw, self.cfg, self.blocks)
        self.assertEqual(m['scale_depth'], 10.0)                                           # the zero bins (2%) are not part of the median
        raw['scale_bins'] = [0.0] * (len(raw['scale_bins']) // 2) + raw['scale_bins'][len(raw['scale_bins']) // 2:]
        self.assertEqual(S.metrics(raw, self.cfg, self.blocks)['scale_depth'], 10.0)       # even with half of them zero

    def test_flags(self):
        t = self.call(self.raw(4.0, base=2.0)); self.assertIn('low_depth', t.call_flags)
        raw = self.raw(1.0); n = len(raw['block_depth']); bd = np.array(raw['block_depth']); bd[: int(0.4 * n)] = 30.0; raw['block_depth'] = bd.tolist()   # 40% of blocks at 3x, the median still ~1
        t = self.call(raw); self.assertEqual(t.bchrom, 'B-'); self.assertIn('fraction_disagrees', t.call_flags)
        raw = self.raw(4.0); raw['scale_bins'] = raw['scale_bins'][:100]; self.assertIn('few_scale_bins', self.call(raw).call_flags)
        raw = self.raw(4.0); raw['block_depth'] = raw['block_depth'][:63]; c = dict(self.cfg); c['blocks'] = dict(tiers=['B'], clean_only=False); self.assertIn('few_blocks', S.call_one(S.metrics(raw, c, self.blocks), c)['flags'])

    def test_block_selection_follows_tier_and_clean(self):
        self.assertEqual(int(S.select_blocks(self.blocks, self.cfg).sum()), 58)
        c = dict(self.cfg); c['blocks'] = dict(tiers=['A'], clean_only=True); self.assertEqual(int(S.select_blocks(self.blocks, c).sum()), 58 - 10)
        c['blocks'] = dict(tiers=['A', 'B'], clean_only=False); self.assertEqual(int(S.select_blocks(self.blocks, c).sum()), 63)

    def test_measurements_from_another_definition_are_refused_and_the_cli_calls_a_folder(self):
        raw = self.raw(4.0); raw['definition'] = 'other'
        with self.assertRaises(ValueError):
            S.call_table([raw], self.cfg)
        d = tempfile.mkdtemp()
        for name, fold in (('A', 1.0), ('B', 4.0)):
            r = self.raw(fold, sample=name); r.pop('SampleID')
            pm.write_json(f'{d}/{name}.panel.json', dict(schema=1, sample=name, bam='/x.bam', bam_size=1, panel_version='t', measured='now', samtools='s', modules={'bchrom': r}))
        cfg_path = d + '/bchrompanel.yaml'; yaml.safe_dump(self.cfg, open(cfg_path, 'w'))
        import cichlid_common.common_light as cl
        orig = cl.config_path; cl.config_path = lambda n: cfg_path if n == 'bchrompanel.yaml' else orig(n)
        try:
            S.main(['--panel', d, '--out', d + '/scores.tsv'])
        finally:
            cl.config_path = orig
        t = pd.read_csv(d + '/scores.tsv', sep='\t').set_index('SampleID'); self.assertEqual(t.bchrom.to_dict(), {'A': 'B-', 'B': 'B+'})


class TestValidateScale(unittest.TestCase):
    def test_the_sampled_scale_reproduces_the_exact_score(self):
        cfg, tmp = make_def(n_bins=3000); blocks, bins, _ = M.load_definition(cfg); td = tempfile.mkdtemp(); rng = np.random.default_rng(3); rep = []
        for run, fold in (('R1', 1.0), ('R2', 5.0)):
            arrays = {c: np.clip(rng.normal(10, 2.0, L // 50), 0.5, None).astype(np.float32) for c, L in LENGTHS.items()}
            for a in arrays.values():
                a[rng.random(len(a)) < 0.02] = 0
            for r in blocks.itertuples():
                arrays[r.chrom][r.start // 50: r.end // 50] *= fold
            stats = depth.depth_stats(arrays, 50, 5_000_000); depth.save_track(f'{td}/{run}', arrays, LENGTHS, stats, 50)
            sel = S.select_blocks(blocks, cfg).to_numpy()
            bd = np.array([arrays[r.chrom][r.start // 50: -(-r.end // 50)].mean() for r in blocks.itertuples()])
            rep.append(dict(run=run, tierA_median_scaled=float(np.median(bd[sel] / stats['median']))))
        t = V.compare(td, cfg, pd.DataFrame(rep))
        self.assertEqual(sorted(t.run), ['R1', 'R2'])
        self.assertLess((t.score_ratio - 1).abs().max(), 0.03); self.assertLess((t.scale_ratio - 1).abs().max(), 0.03)
        self.assertLess((t.exact_vs_reported - 1).abs().max(), 1e-6)                       # the exact route reproduces the discovery's own report to float32 rounding


class TestCompareWithMosdepth(unittest.TestCase):
    def test_a_panel_measurement_against_a_mosdepth_track_of_the_same_reads(self):
        import gzip
        cfg, tmp = make_def(n_bins=3000); blocks, bins, _ = M.load_definition(cfg); rng = np.random.default_rng(4)
        arrays = {c: np.clip(rng.normal(10, 2.0, L // 50), 0.5, None).astype(np.float32) for c, L in LENGTHS.items()}
        for a in arrays.values():
            a[rng.random(len(a)) < 0.02] = 0
        for r in blocks.itertuples():
            arrays[r.chrom][r.start // 50: r.end // 50] *= 4.0
        d = tempfile.mkdtemp(); mo = d + '/S.regions.bed.gz'
        with gzip.open(mo, 'wt') as fh:
            for c, a in arrays.items():
                for i, v in enumerate(a):
                    fh.write(f'{c}\t{i * 50}\t{(i + 1) * 50}\t{v}\n')
        depth_of = lambda c, s, e: float(arrays[c][s // 50: -(-e // 50)].mean())          # the panel route reads the same bases
        out = M.measure_bam('x.bam', cfg, fake_run(cfg, depth_of))
        pm.write_json(d + '/S.panel.json', dict(schema=1, sample='S', bam='/x', bam_size=1, panel_version='t', measured='n', samtools='s', modules={'bchrom': out}))
        r = V.compare_bam(d + '/S.panel.json', mo, cfg)
        self.assertLess(abs(r['score_ratio'] - 1), 0.03); self.assertAlmostEqual(r['block_ratio_median'], 1.0, places=2); self.assertEqual(r['block_ratio_within_5pct'], 1.0)
        out['definition'] = 'other'; pm.write_json(d + '/T.panel.json', dict(schema=1, sample='T', bam='/x', bam_size=1, panel_version='t', measured='n', samtools='s', modules={'bchrom': out}))
        with self.assertRaises(ValueError):
            V.compare_bam(d + '/T.panel.json', mo, cfg)


if __name__ == '__main__':
    unittest.main()
