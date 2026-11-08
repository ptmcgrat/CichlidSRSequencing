"""LG5 panel measurement and call on synthetic BAM output (fake samtools): definition build, measurement, definition key, genotype calls and their flags, the panel JSON round trip."""
import json, os, random, sys, tempfile, types, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import yaml
from cichlid_common.common_light import config_path
from sd_genotypes.lg5inv import call as C, definition as D, markers as K, measure as M
from sd_genotypes.panel import measure as pm

CFG0 = yaml.safe_load(open(config_path('lg5panel.yaml')))
rc = K.rc


def base(i):
    return 'ACGT'[(i * 2654435761 >> 7) % 4]


def fetch(a, b):
    out = []
    for i in range(a, b + 1):
        out.append(D.CLIP_R15[i - 19657763] if 19657763 <= i <= 19657777 else base(i))
    return ''.join(out)


def marker_table(n=300):
    pos = [19_700_000 + 20000 * i for i in range(n)]
    ref = [fetch(p, p) for p in pos]
    return pd.DataFrame(dict(pos=pos, ref=ref, alt=[{'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[r] for r in ref], type='SNV'))


def make_def(tmp, n=300):
    mp, jp, k = D.build(fetch, marker_table(n), tmp)
    cfg = dict(CFG0); cfg['definition_dir'] = tmp
    return cfg


def fake_run(cfg, ref_n, alt_n, sam_lines, contig=('NC_135171.1', 42954235), log=None):
    markers, _ = M.load_definition(cfg)

    def run(cmd):
        if log is not None:
            log.append(cmd)
        if cmd[:2] == ['samtools', '--version']:
            return types.SimpleNamespace(stdout='samtools 1.x\n')
        if cmd[:3] == ['samtools', 'view', '-H']:
            return types.SimpleNamespace(stdout=f'@HD\tVN:1.6\n@SQ\tSN:{contig[0]}\tLN:{contig[1]}\n')
        if cmd[:2] == ['samtools', 'mpileup']:
            lines = []
            for p, r, a, nr, na in zip(markers.pos, markers.ref, markers.alt, ref_n, alt_n):
                bases = (r * (nr // 2) + r.lower() * (nr - nr // 2)) + a * na + ('^F' + 'g' + '$' if False else '')
                lines.append('\t'.join(['NC_135171.1', str(p), 'N', str(nr + na), bases or '*', 'I' * (nr + na)]))
            return types.SimpleNamespace(stdout='\n'.join(lines) + '\n')
        if cmd[:2] == ['samtools', 'view']:
            return types.SimpleNamespace(stdout='\n'.join(sam_lines) + '\n')
        raise AssertionError(cmd)
    return run


def sam(seq):
    return '\t'.join(['r', '0', 'NC_135171.1', '26726600', '60', '100M', '*', '0', '0', seq, 'I' * len(seq)])


class TestDefinition(unittest.TestCase):
    def test_build_writes_markers_and_junction_probes_from_the_reference(self):
        tmp = tempfile.mkdtemp(); mp, jp, n = D.build(fetch, marker_table(), tmp)
        m = pd.read_csv(mp, sep='\t'); pr = json.load(open(jp))['probes']
        self.assertEqual((n, list(m.columns), m.pos.is_monotonic_increasing), (300, ['chrom', 'pos', 'ref', 'alt'], True))
        self.assertTrue(pr['A'].startswith(D.CLIP_L15) and pr['B'].endswith(D.CLIP_R15) and len(pr['A']) == len(pr['B']) == len(pr['C']) == 30)
        self.assertEqual(pr['C'], fetch(26726671, 26726700))                           # the standard sequence is contiguous GT3 across the breakpoint

    def test_a_wrong_reference_or_marker_file_is_refused(self):
        bad = marker_table(); bad.loc[3, 'ref'] = next(c for c in 'ACGT' if c != bad.loc[3, 'ref'])
        with self.assertRaises(ValueError):
            D.build(fetch, bad, tempfile.mkdtemp())
        with self.assertRaises(ValueError):
            D.build(lambda a, b: 'N' * (b - a + 1), marker_table(), tempfile.mkdtemp())


class TestMeasure(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(); self.cfg = make_def(self.tmp); _, self.pr = M.load_definition(self.cfg)

    def test_counts_per_marker_and_junction_reads(self):
        ref_n = [7] * 300; alt_n = [0] * 300; alt_n[5] = 3; ref_n[7] = 0; alt_n[7] = 9
        noise = ''.join(random.Random(1).choice('ACGT') for _ in range(100))
        sam_lines = [sam('T' * 20 + self.pr['A'] + 'G' * 50), sam(rc('T' * 30 + self.pr['B'] + 'G' * 40)), sam('A' * 10 + self.pr['C'] + 'C' * 60), sam('C' * 30 + self.pr['C'] + 'A' * 40), sam('G' * 5 + self.pr['C'] + 'T' * 65), sam(noise)]
        log = []
        out = M.measure_bam('x.bam', self.cfg, fake_run(self.cfg, ref_n, alt_n, sam_lines, log=log))
        self.assertEqual((out['ref'][0], out['alt'][5], out['ref'][7], out['alt'][7], out['n_markers']), (7, 3, 0, 9, 300))
        self.assertEqual(out['junction'], dict(A=1, B=1, C=3, reads=6))
        mp = next(c for c in log if c[1] == 'mpileup')
        self.assertTrue(all(x in mp for x in ('-B', '-q', '20', '-r', 'NC_135171.1:19680001-26730000')) and mp[mp.index('-Q') + 1] == '0')    # no base-quality filter
        jv = next(c for c in log if c[1] == 'view' and c[2] == '-F')
        self.assertEqual((jv[3], jv[-1]), ('3840', 'NC_135171.1:26726301-26727060'))

    def test_the_definition_key_changes_with_what_is_measured_and_not_otherwise(self):
        k = M.definition_key(self.cfg)
        self.assertEqual(M.definition_key(dict(self.cfg)), k)
        c2 = dict(self.cfg); c2['min_mapq'] = 30
        self.assertNotEqual(M.definition_key(c2), k)
        c3 = dict(self.cfg); c3['thresholds'] = dict(c3['thresholds'], min_site_depth=9)     # a CALL threshold: not part of the definition
        self.assertEqual(M.definition_key(c3), k)
        tmp2 = tempfile.mkdtemp(); D.build(fetch, marker_table(299), tmp2)                    # a different marker list
        c4 = dict(self.cfg); c4['definition_dir'] = tmp2
        self.assertNotEqual(M.definition_key(c4), k)

    def test_a_bam_from_another_reference_or_without_the_definition_is_refused(self):
        run = lambda contig: fake_run(self.cfg, [5] * 300, [0] * 300, [], contig=contig)
        with self.assertRaises(ValueError):
            M.measure_bam('x.bam', self.cfg, run(('NC_135171.1', 40_000_000)))
        with self.assertRaises(ValueError):
            M.measure_bam('x.bam', self.cfg, run(('chr5', 42954235)))
        c = dict(self.cfg); c['definition_dir'] = tempfile.mkdtemp()
        with self.assertRaises(FileNotFoundError):
            M.measure_bam('x.bam', c, fake_run(self.cfg, [5] * 300, [0] * 300, []))

    def test_it_runs_as_a_panel_module_and_the_json_round_trips(self):
        d = tempfile.mkdtemp(); bam = d + '/s.all.bam'; open(bam, 'wb').write(b'x'); open(bam + '.bai', 'wb').write(b'x')
        pm.MODULES.setdefault('lg5inv', ('sd_genotypes.lg5inv.measure', 'lg5panel.yaml'))
        res = pm.measure(bam, 'S1', modules=('lg5inv',), run=fake_run(self.cfg, [6] * 300, [0] * 300, []), cfgs={'lg5inv': self.cfg})
        pm.write_json(d + '/S1.panel.json', res)
        rows, skipped = pm.collect(d, 'lg5inv')
        self.assertEqual((len(rows), skipped, rows[0]['SampleID'], rows[0]['n_markers'], rows[0]['definition']), (1, [], 'S1', 300, M.definition_key(self.cfg)))


def raw(rng, n=300, depth=10, p=0.5, standing=0.0, jc=(0, 0, 0)):
    d = rng.poisson(depth, n); d = np.maximum(d, 0)
    pp = np.full(n, p); pp[rng.random(n) < standing] = 0.45                                # sites where the alt allele also segregates on standard chromosomes
    alt = rng.binomial(d, pp)
    return dict(definition='x', ref=(d - alt).tolist(), alt=alt.tolist(), junction=dict(A=jc[0], B=jc[1], C=jc[2], reads=sum(jc)))


class TestCall(unittest.TestCase):
    def one(self, r, cfg=CFG0):
        m = C.metrics(r, cfg); return m, C.call_one(m, cfg)

    def test_the_three_genotypes(self):
        rng = np.random.default_rng(5)
        for p, jc, want in [(0.005, (0, 0, 9), 0), (0.5, (4, 4, 4), 1), (0.985, (9, 8, 0), 2)]:
            m, c = self.one(raw(rng, p=p, jc=jc))
            self.assertEqual((c['n_inverted'], c['flags']), (want, []), (p, m))

    def test_standing_variation_does_not_move_a_standard_fish_because_the_call_uses_the_median(self):
        rng = np.random.default_rng(6)
        m, c = self.one(raw(rng, p=0.004, standing=0.30, jc=(0, 0, 8)))                      # 30% of sites carry the alt at ~0.45: the aggregate fraction rises, the median does not
        self.assertGreater(m['aggregate_fraction'], 0.10); self.assertLess(m['median_site_fraction'], 0.1)
        self.assertEqual(c['n_inverted'], 0)

    def test_flags(self):
        rng = np.random.default_rng(7)
        m, c = self.one(raw(rng, depth=1.5, p=0.5)); self.assertIn('low_depth', c['flags'])
        r = raw(rng, p=0.5, jc=(0, 0, 10)); m, c = self.one(r)                                  # markers say heterozygous, junction says standard
        self.assertEqual((c['n_inverted'], 'junction_disagrees' in c['flags']), (1, True))
        m, c = self.one(raw(rng, p=0.5, jc=(1, 0, 1))); self.assertIn('junction_uninformative', c['flags'])
        r = dict(definition='x', ref=[7] * 300, alt=[3] * 300, junction=dict(A=3, B=3, C=3, reads=9)); m, c = self.one(r)    # every site at 0.3: near the first cut
        self.assertEqual(c['n_inverted'], 1); self.assertIn('near_threshold', c['flags'])
        r = dict(definition='x', ref=[8] * 100, alt=[0] * 100, junction=dict(A=0, B=0, C=9, reads=9)); m, c = self.one(r); self.assertIn('few_markers', c['flags'])

    def test_definitions_are_not_mixed_and_the_cli_calls_a_folder(self):
        rng = np.random.default_rng(8); d = tempfile.mkdtemp()
        samples = [('A', 0.004, (0, 0, 9)), ('B', 0.5, (4, 4, 4)), ('C', 0.985, (9, 8, 0))]
        for name, p, jc in samples:
            res = dict(schema=1, sample=name, bam='/x.bam', bam_size=1, panel_version='t', measured='now', samtools='s', modules={'lg5inv': dict(raw(rng, p=p, jc=jc), contig_length=1, region='r', n_markers=300)})
            pm.write_json(f'{d}/{name}.panel.json', res)
        C.main(['--panel', d, '--out', d + '/calls.tsv'])
        t = pd.read_csv(d + '/calls.tsv', sep='\t').set_index('SampleID')
        self.assertEqual(t.genotype.to_dict(), {'A': 'std/std', 'B': 'std/inv', 'C': 'inv/inv'})
        rows, _ = pm.collect(d, 'lg5inv'); rows[0]['definition'] = 'other'
        with self.assertRaises(ValueError):
            C.call_table(rows, CFG0)


if __name__ == '__main__':
    unittest.main()
