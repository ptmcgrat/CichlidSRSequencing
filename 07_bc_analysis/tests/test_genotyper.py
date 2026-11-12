"""bc_analysis genotype (module 07): the sites from a parents VCF, the allele counts, the HMM checked against enumeration of all paths, the simulation, and the commands end to end."""
import gzip, io, itertools, json, os, shutil, sys, tempfile, unittest
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from bc_analysis.genotyper import sites as S, counts as C, hmm as H, simulate as SIM, cli as GC

CFG = S.load_config()
PARENTS = ['YH_008_m', 'MC-1R6R-f', 'MC-3R6R-f', 'MC-1G2G-f', 'MC-5G11G-f', 'MCYHF1-002', 'MCYHF1-003']
DAMS = ['MC-3R6R-f', 'MC-1G2G-f', 'MC-5G11G-f']


# ====================================================================== a parents VCF with one site of every kind
def call(gt, dp=11, ad=None):
    ad = ad or {'0/0': (dp, 0), '0/1': (dp // 2, dp - dp // 2), '1/1': (0, dp)}[gt]
    return f'{gt}:{ad[0]},{ad[1]}:{dp}'


def write_vcf(path, sites, filler=300):
    """sites: list of (chrom, pos, ref, alt, {sample: call string}); samples missing from the dict get a hom-ref call at depth 10. Filler records of depth 10 set the median depth of every parent."""
    lines = ['##fileformat=VCFv4.2', '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">', '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(PARENTS)]
    for i in range(filler):
        lines.append(f'c1\t{1000 + i}\t.\tA\tC\t50\t.\t.\tGT:AD:DP\t' + '\t'.join(call('0/0', 10) for _ in PARENTS))
    for chrom, pos, ref, alt, g in sites:
        lines.append(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t.\t.\tGT:AD:DP\t' + '\t'.join(g.get(p, call('0/0', 10)) for p in PARENTS))
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'wt') as fh:
        fh.write('\n'.join(lines) + '\n')


def parents_case(yh, f1, dams=None, **over):
    """Genotype strings: YH sire `yh`, F1 dam `f1`, the three dams `dams` (a dict or one genotype for all)."""
    d = dict(zip(DAMS, [dams] * 3)) if isinstance(dams, str) else (dams or {})
    g = {'YH_008_m': call(yh), 'MC-1R6R-f': call(f1)}
    g.update({k: call(v) for k, v in d.items()})
    g.update(over)
    return g


def planted():
    return [
        ('c1', 5000, 'A', 'C', parents_case('1/1', '0/0', '0/0')),                                       # YH carries ALT: informative for all three dams, yh_is_ref False
        ('c1', 6000, 'G', 'T', parents_case('0/0', '1/1', '1/1')),                                       # YH carries REF: yh_is_ref True
        ('c1', 7000, 'A', 'C', parents_case('1/1', '0/0', {'MC-3R6R-f': '0/0', 'MC-1G2G-f': '0/1', 'MC-5G11G-f': '0/0'})),   # 1G2G is heterozygous: no site for ITS brood, still one for the others
        ('c1', 8000, 'A', 'C', parents_case('1/1', '0/1', '0/0')),                                       # the F1 dam is heterozygous: nobody
        ('c1', 9000, 'A', 'C', parents_case('1/1', '1/1', '1/1')),                                       # YH and MC carry the same allele: nobody
        ('c1', 10000, 'A', 'C', parents_case('1/1', '0/0', '0/0', **{'YH_008_m': call('1/1', 3)})),      # the YH sire has 3 reads: below the minimum
        ('c1', 11000, 'A', 'C', parents_case('1/1', '0/0', '0/0', **{'YH_008_m': call('1/1', 100)})),    # 100 reads: a collapsed repeat
        ('c1', 12000, 'A', 'C', parents_case('1/1', '0/0', '0/0', **{'YH_008_m': call('1/1', 16, (4, 12))})),   # 4 reads of the other allele: not homozygous enough
        ('c1', 13000, 'A', 'C,G', parents_case('1/1', '0/0', '0/0')),                                    # multi-allelic: skipped
        ('c1', 14000, 'AT', 'A', parents_case('1/1', '0/0', '0/0')),                                     # an indel: skipped
        ('c1', 15000, 'A', 'AT', parents_case('1/1', '0/0', '0/0')),                                     # an insertion with a one-base REF: skipped
        ('c1', 16000, 'A', 'C', parents_case('1/1', '0/0', '0/0', **{'YH_008_m': call('1/1', 12, (2, 10))})),   # 2 reads of the other allele among 12: allowed by the reads rule, refused by the fraction (0.17 > 0.1)
        ('c1', 17000, 'A', 'C', parents_case('1/1', '0/0', '0/0', **{'YH_008_m': call('1/1', 12, (1, 11))})),   # 1 among 12 (0.083): accepted
        ('c1', 18000, 'A', 'C', parents_case('1/1', '0/0', '0/0', **{'YH_008_m': call('1/1', 5, (0, 5))})),     # 5 reads: above half the median (5) but below the absolute minimum of 6
        ('c2', 2000, 'A', 'C', parents_case('1/1', '0/0', {'MC-3R6R-f': '1/1', 'MC-1G2G-f': '0/0', 'MC-5G11G-f': '0/0'})),   # the 3R6R dam carries the YH allele: no site for her brood
    ]


class TestSites(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()
        cls.vcf = os.path.join(cls.dir, 'bc_family.vcf.gz')
        write_vcf(cls.vcf, planted())
        cfg = json.loads(json.dumps(CFG)); cfg['sites']['depth_sample_sites'] = 100000
        cls.cfg = cfg
        cls.res = S.define_sites(cls.vcf, cfg, log=lambda *_: None)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def positions(self, dam):
        return [(c, p) for c, p, *_ in self.res['sites'][dam]]

    def test_each_brood_gets_the_sites_of_its_own_dam(self):
        self.assertEqual(self.positions('MC-3R6R-f'), [('c1', 5000), ('c1', 6000), ('c1', 7000), ('c1', 17000)])
        self.assertEqual(self.positions('MC-1G2G-f'), [('c1', 5000), ('c1', 6000), ('c1', 17000), ('c2', 2000)])                  # her own heterozygous site 7000 is out; 3R6R's carrying the YH allele at c2:2000 does not matter to her
        self.assertEqual(self.positions('MC-5G11G-f'), [('c1', 5000), ('c1', 6000), ('c1', 7000), ('c1', 17000), ('c2', 2000)])

    def test_the_polarity_is_recorded_by_allele_not_by_reference(self):
        sites = {(c, p): y for c, p, _, _, y in self.res['sites']['MC-3R6R-f']}
        self.assertEqual((sites[('c1', 5000)], sites[('c1', 6000)]), (False, True))                                 # the YH allele is ALT at 5000, REF at 6000

    def test_the_depth_and_purity_rules_and_the_skipped_records(self):
        for dam in DAMS:
            for bad in (8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 18000):
                self.assertNotIn(('c1', bad), self.positions(dam), bad)

    def test_attrition_is_monotone_and_the_last_step_is_the_sites(self):
        att = pd.DataFrame(self.res['attrition'])
        for dam, g in att.groupby('dam'):
            v = g.set_index('step').sites.reindex(S.STEPS).tolist()
            self.assertEqual(v, sorted(v, reverse=True)); self.assertEqual(v[-1], len(self.res['sites'][dam]))
        self.assertEqual(att[(att.dam == 'MC-3R6R-f') & (att.step == 'snp_biallelic')].sites.iloc[0], 300 + 12)      # the filler and the twelve biallelic SNPs of the plant (the multi-allelic site, the deletion and the insertion are not SNPs)

    def test_the_median_depths_come_from_the_vcf_not_from_a_constant(self):
        self.assertEqual({p: round(m) for p, m in self.res['medians'].items() if p in ('YH_008_m', 'MC-1R6R-f')}, {'YH_008_m': 10, 'MC-1R6R-f': 10})

    def test_the_files(self):
        out = os.path.join(self.dir, 'out'); S.write_sites(self.res, out, self.vcf, self.cfg, fai=[('c1', 2_000_000), ('c2', 1_000_000)])
        t = pd.read_csv(os.path.join(out, 'sites_MC-5G11G-f.tsv.gz'), sep='\t'); self.assertEqual(list(t.columns), ['chrom', 'pos', 'ref', 'alt', 'yh_is_ref']); self.assertEqual(len(t), 5)
        tg = pd.read_csv(os.path.join(out, 'targets_MC-5G11G-f.tsv'), sep='\t', header=None); self.assertEqual(tg.shape, (5, 2))
        d = json.load(open(os.path.join(out, 'definition.json'))); self.assertEqual(d['key'], self.res['key']); self.assertEqual(d['sites']['MC-1G2G-f'], 4)
        pm = pd.read_csv(os.path.join(out, 'sites_per_mb.tsv'), sep='\t'); self.assertAlmostEqual(float(pm[(pm.dam == 'MC-5G11G-f') & (pm.chrom == 'c2')].sites_per_mb.iloc[0]), 1.0)
        self.assertTrue(os.path.exists(os.path.join(out, 'attrition.tsv')))

    def test_the_definition_key_follows_the_parameters_and_the_vcf(self):
        k = S.definition_key(self.cfg, self.vcf)
        c2 = json.loads(json.dumps(self.cfg)); c2['sites']['min_dp'] = 7
        self.assertNotEqual(k, S.definition_key(c2, self.vcf)); self.assertEqual(k, S.definition_key(json.loads(json.dumps(self.cfg)), self.vcf))
        other = os.path.join(self.dir, 'other.vcf'); write_vcf(other, planted()[:3]); self.assertNotEqual(k, S.definition_key(self.cfg, other))

    def test_a_missing_parent_and_a_missing_depth_are_named(self):
        p = os.path.join(self.dir, 'x.vcf'); write_vcf(p, planted()[:2])
        text = open(p).read().replace('MC-1G2G-f', 'MC-OTHER-f'); open(p, 'w').write(text)
        with self.assertRaises(ValueError) as cm: S.define_sites(p, self.cfg, log=lambda *_: None)
        self.assertIn('MC-1G2G-f', str(cm.exception))
        q = os.path.join(self.dir, 'y.vcf'); write_vcf(q, [], filler=0)
        with self.assertRaises(ValueError) as cm: S.define_sites(q, self.cfg, log=lambda *_: None)
        self.assertIn('no depth', str(cm.exception))


# ====================================================================== counts
SITES = {('c1', 100): ('A', 'C', False), ('c1', 200): ('A', 'C', True), ('c1', 300): ('G', 'T', False)}


class TestCounts(unittest.TestCase):
    def test_counts_are_read_by_allele_name_whatever_the_order(self):
        text = '\n'.join(['c1\t100\tA\tC,<*>\t3,2,0',                 # YH allele is ALT (C): 2 YH reads, 3 MC reads
                          'c1\t200\tA\tC,<*>\t1,4,0',                # YH allele is REF here: 1 YH read, 4 MC (ALT) reads
                          'c1\t300\tG\tA,T,<*>\t0,5,2,0'])           # the first observed alternative is a third allele; T is second
        rows = {(r[0], r[1]): r[2:] for r in C.parse_pileup(text, SITES)}
        self.assertEqual(rows[('c1', 100)], (3, 2, 0)); self.assertEqual(rows[('c1', 200)], (4, 1, 1)); self.assertEqual(rows[('c1', 300)], (0, 2, 0))      # (n_mc, n_yh, yh_is_ref); site 300: the REF G has 0 reads, T has 2

    def test_no_alternative_reads_means_zero_YH_reads_at_an_alt_polarity_site(self):
        rows = C.parse_pileup('c1\t100\tA\t<*>\t6,0', SITES)
        self.assertEqual(rows, [('c1', 100, 6, 0, 0)])
        self.assertEqual(C.parse_pileup('c1\t100\tA\t.\t6', SITES), [('c1', 100, 6, 0, 0)])

    def test_unwanted_lines_are_dropped(self):
        text = '\n'.join(['c1\t999\tA\tC\t5,5',                      # not an informative site
                          'c1\t100\tG\tC\t5,5',                      # the reference base is not the VCF's
                          'c1\t200\tA\tC,<*>\t0,0,0',                # no reads
                          'c1\t300\tG\tT\tx,y', 'short line'])
        self.assertEqual(C.parse_pileup(text, SITES), [])

    def test_the_pileup_command_has_the_filters_the_original_lacked(self):
        s = C.mpileup_script('/g.fa', '/b.bam', '/t.tsv', CFG['counts'])
        for frag in ('-T /t.tsv', '-q 20', '-Q 20', '-I', '--ff UNMAP,SECONDARY,QCFAIL,DUP', '-a FORMAT/AD', '-f /g.fa'):
            self.assertIn(frag, s)


class Fake:
    """rclone, samtools and the pileup, as the build tests do: the BAM 'downloaded' is a file of the size the ledger gave."""
    def __init__(self, size=1000, pileup='c1\t100\tA\tC,<*>\t3,2,0\nc1\t200\tA\tC,<*>\t1,4,0\n', fail=None):
        self.size, self.pileup, self.fail, self.calls = size, pileup, fail, []

    def run(self, argv, input_text=None):
        self.calls.append(argv[0])
        if argv[0] == 'rclone':
            dest = argv[3]; os.makedirs(dest, exist_ok=True)
            s = [a for a in argv if a.endswith('.all.bam')][0].replace('.all.bam', '')
            open(os.path.join(dest, f'{s}.all.bam'), 'wb').write(b'x' * self.size); open(os.path.join(dest, f'{s}.all.bai'), 'wb').write(b'i')
            return 0, '', ''
        return 0, '', ''                                                    # samtools quickcheck

    def script(self, text):
        self.calls.append('pileup')
        return (1, '', 'boom: no index') if self.fail else (0, self.pileup, '')


class TestCountSample(unittest.TestCase):
    def setUp(self):
        from variant_calling import build as B
        self.dir = tempfile.mkdtemp()
        self.env = B.Env(self.dir + '/local', self.dir + '/cloud', 'G', '/ref.fa')
        self.sites = os.path.join(self.dir, 'sites.tsv.gz')
        pd.DataFrame([('c1', 100, 'A', 'C', 0), ('c1', 200, 'A', 'C', 1)], columns=['chrom', 'pos', 'ref', 'alt', 'yh_is_ref']).to_csv(self.sites, sep='\t', index=False)
        self.staged = self.dir + '/staged/'

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def go(self, fake, key='K1', want=1000, keep=False, sample='S1'):
        return C.count_sample(self.env, CFG, sample, 'MC-3R6R-f', self.sites, '/targets.tsv', key, self.dir + '/out/', self.staged, want, keep, fake.run, fake.script)

    def test_a_fish_is_counted_stored_and_its_bam_deleted(self):
        f = Fake(); m = self.go(f)
        self.assertEqual((m['status'], m['n_sites'], m['n_reads'], m['key'], m['bam_bytes'], m['dam']), ('ok', 2, 10, 'K1', 1000, 'MC-3R6R-f'))
        self.assertAlmostEqual(m['yh_fraction'], 0.3)
        cc = C.read_counts(self.dir + '/out/', 'S1')
        np.testing.assert_array_equal(cc['c1'][1], [3, 4]); np.testing.assert_array_equal(cc['c1'][2], [2, 1]); np.testing.assert_array_equal(cc['c1'][3], [False, True])
        self.assertFalse(os.path.exists(self.staged + 'S1')); self.assertEqual(f.calls.count('pileup'), 1)

    def test_current_follows_the_key_and_the_bam(self):
        self.go(Fake()); out = self.dir + '/out/'
        self.assertTrue(C.current(out, 'S1', 'K1', 1000)); self.assertTrue(C.current(out, 'S1', 'K1'))
        self.assertFalse(C.current(out, 'S1', 'K2', 1000), 'another definition of the sites'); self.assertFalse(C.current(out, 'S1', 'K1', 1001), 'a realigned BAM has another size'); self.assertFalse(C.current(out, 'NOPE', 'K1'))

    def test_a_failure_is_recorded_not_raised_and_is_not_current(self):
        m = self.go(Fake(fail=True))
        self.assertEqual(m['status'], 'failed'); self.assertIn('pileup of S1', m['error']); self.assertFalse(C.current(self.dir + '/out/', 'S1', 'K1'))
        self.assertFalse(os.path.exists(self.staged + 'S1'))

    def test_a_fish_whose_later_run_failed_is_not_current_even_though_its_old_counts_are_still_there(self):
        self.go(Fake()); out = self.dir + '/out/'; self.assertTrue(C.current(out, 'S1', 'K1'))
        m = self.go(Fake(fail=True)); self.assertEqual(m['status'], 'failed')
        self.assertTrue(os.path.exists(C.counts_paths(out, 'S1')[0])); self.assertFalse(C.current(out, 'S1', 'K1'))

    def test_keep_bams_keeps_them(self):
        self.go(Fake(), keep=True); self.assertTrue(os.path.exists(self.staged + 'S1/S1.all.bam'))

    def test_a_download_of_the_wrong_size_fails_that_fish_only(self):
        m = self.go(Fake(size=999), want=1000); self.assertEqual(m['status'], 'failed'); self.assertIn('ledger says', m['error'])


# ====================================================================== the HMM, checked against enumeration of every path
def brute_posterior(ll0, ll1, lam, A, young, Wt, term_factor):
    """P(state 1) of each bin by summing over all 2^n paths of the chain (small n): interior hazards by age, a first segment with its own hazards, a last segment weighted by its length."""
    n = len(ll0); W = int(Wt) if Wt else 0; cap = max(A, W, 1)
    h = np.full(cap + 1, min(lam, 0.5)); h[:A] = min(lam * young, 0.5)
    hF = np.full(W, min(lam * term_factor, 0.5)) if W else np.zeros(0)
    if W:
        hF[W - 1] = min(lam, 0.5)
    fend = np.ones(cap + 1)
    if W:
        fend[:max(W - 1, 0)] = term_factor
    m = np.maximum(ll0, ll1); e = np.stack([np.exp(ll0 - m), np.exp(ll1 - m)], 1)
    num = np.zeros(n); den = 0.0
    for path in itertools.product((0, 1), repeat=n):
        w = 0.5 * e[0][path[0]]
        first, b, age = bool(W), 0, cap
        for i in range(1, n):
            sw = path[i] != path[i - 1]
            hz = hF[b] if first else h[age]
            w *= (hz if sw else 1 - hz) * e[i][path[i]]
            if sw:
                first, age = False, 0
            elif first:
                if b + 1 < W:
                    b += 1
                else:
                    first, age = False, cap
            else:
                age = min(age + 1, cap)
        if not first:
            w *= fend[age]
        den += w
        for i in range(n):
            if path[i] == 1:
                num[i] += w
    return num / den


class TestHmm(unittest.TestCase):
    def test_bin_counts(self):
        pos = np.array([50_000, 60_000, 150_000, 150_001, 950_000, 5_000_000]); n_mc = np.array([1, 5, 0, 1, 2, 3]); n_yh = np.array([0, 0, 4, 1, 0, 9]); ref = np.array([True, True, False, False, True, False])
        c = H.bin_counts(pos, n_mc, n_yh, ref, 1_000_000, 100_000, site_cap=2, bin_cap=20)
        self.assertEqual(c.shape, (4, 10))
        self.assertEqual(list(c[:, 0]), [0, 3, 0, 0])               # bin 0: two ref-polarity sites, MC reads 1 and 5 capped at 2, YH reads 0
        self.assertEqual(list(c[:, 1]), [0, 0, 3, 1])               # bin 1: alt-polarity sites: YH reads 4 -> 2 and 1, MC reads 0 and 1
        self.assertEqual(list(c[:, 9]), [0, 2, 2, 2])               # bin 9: the site at 950 kb (ref polarity, 2 MC reads) and the site at 5 Mb, beyond the chromosome, put in the last bin (alt polarity: YH 9 -> 2, MC 3 -> 2)
        self.assertEqual(H.bin_counts([], [], [], [], 1_000_000, 100_000).sum(), 0)

    def test_the_bin_cap_scales_a_bin_down(self):
        c = H.bin_counts(np.arange(10) * 1000 + 1, np.full(10, 2), np.full(10, 2), np.zeros(10, bool), 100_000, 100_000, site_cap=2, bin_cap=20)
        self.assertAlmostEqual(c.sum(), 20.0); self.assertAlmostEqual(c[2, 0] / c[3, 0], 1.0)

    def test_forward_backward_equals_enumeration_of_all_paths(self):
        rng = np.random.default_rng(3)
        for A, young, Wt, tf in ((3, 1e-3, None, 1.0), (3, 1e-3, 2, 1e-2), (1, 1.0, None, 1.0), (4, 0.05, 3, 0.1), (2, 1.0, 2, 0.2), (1, 1.0, 4, 1e-3), (1, 1.0, 1, 1.0), (5, 1e-2, 7, 0.05)):
            ll0, ll1 = rng.normal(0, 2, 9), rng.normal(0, 2, 9)
            got = H.forward_backward(ll0, ll1, 0.15, A, young, Wt, tf)
            np.testing.assert_allclose(got, brute_posterior(ll0, ll1, 0.15, A, young, Wt, tf), atol=1e-9, err_msg=str((A, young, Wt, tf)))

    def evid(self, n, yh_bins, per_bin=6):
        """Counts that say 'heterozygous' on the bins in yh_bins and 'MC' elsewhere (alt-polarity sites): per bin `per_bin` reads, half YH when heterozygous, none when MC."""
        c = np.zeros((4, n)); c[2, yh_bins] = per_bin / 2; c[3, yh_bins] = per_bin / 2; c[3, [i for i in range(n) if i not in set(yh_bins)]] = per_bin
        return c

    def fb(self, c, **kw):
        return H.forward_backward(*H.bin_loglik(c, 0.01, 0.6), 1e-3, kw.pop('A', 20), **kw)

    def test_two_states_are_called_where_the_reads_say(self):
        p = self.fb(self.evid(60, list(range(0, 25))), A=20)
        self.assertTrue((p[:23] > .9).all()); self.assertTrue((p[28:] < .1).all())

    def test_a_short_interior_segment_is_refused_a_long_one_is_called(self):
        short = self.evid(80, list(range(30, 34)))                    # 400 kb block of strong evidence inside an MC chromosome
        self.assertTrue((self.fb(short, A=20) < .5).all(), 'an interior segment shorter than the minimum is not called')
        self.assertTrue((self.fb(short, A=1, young=1.0) > .5)[30:34].all(), 'without the limit it is called')
        long = self.evid(80, list(range(30, 55)))
        self.assertTrue((self.fb(long, A=20) > .5)[32:53].all(), 'a 2.5 Mb segment is called')

    def test_one_switch_near_an_end_costs_an_ordinary_switch_and_the_terminal_limit_stretches_or_drops_it(self):
        n = 80
        strong_left, strong_right = self.evid(n, [0, 1, 2]), self.evid(n, [77, 78, 79])
        # free: a 300 kb terminal segment with strong reads is called exactly where the reads are
        p = self.fb(strong_left, A=20); self.assertTrue((p[0:3] > .5).all() and p[3] < .5)
        p = self.fb(strong_right, A=20); self.assertTrue((p[77:80] > .5).all() and p[76] < .5)
        # with a terminal minimum of 5 bins (500 kb) the same strong reads give a call STRETCHED to the minimum, at both ends
        p = self.fb(strong_left, A=20, Wt=5, term_factor=1e-30); self.assertTrue((p[0:5] > .5).all() and p[5] < .5)
        p = self.fb(strong_right, A=20, Wt=5, term_factor=1e-30); self.assertTrue((p[75:80] > .5).all() and p[74] < .5)
        # weak reads (one YH and one MC read per bin against six MC reads in the neighbours) are called when free and DROPPED under the minimum
        weak = np.zeros((4, n)); weak[3, 3:] = 6; weak[2, 0:3] = 1; weak[3, 0:3] = 1
        self.assertTrue((self.fb(weak, A=20)[0:3] > .5).all()); self.assertTrue((self.fb(weak, A=20, Wt=5, term_factor=1e-30)[0:3] < .5).all())
        # a terminal segment longer than the minimum is untouched
        self.assertTrue((self.fb(self.evid(n, list(range(0, 8))), A=20, Wt=5, term_factor=1e-30)[:8] > .5).all())

    def test_the_limits_mean_what_the_configuration_says_in_bins(self):
        hp = dict(CFG['hmm']); hp['min_terminal_segment_mb'] = 0.5; hp['min_interior_segment_mb'] = 2.0
        A, Wt = H.segment_limits(hp); self.assertEqual((A, Wt), (19, 5))
        n = 90
        # a terminal segment of exactly Wt = 5 bins is allowed and one of 4 bins is stretched to 5, at both ends
        for yh, lo, hi in (([0, 1, 2, 3, 4], 0, 5), ([n - 5, n - 4, n - 3, n - 2, n - 1], n - 5, n)):
            p = self.fb(self.evid(n, yh), A=A, Wt=Wt, term_factor=1e-30); self.assertTrue((p[lo:hi] > .5).all()); self.assertEqual(int((p > .5).sum()), 5)
        for yh in ([0, 1, 2, 3], [n - 4, n - 3, n - 2, n - 1]):
            p = self.fb(self.evid(n, yh), A=A, Wt=Wt, term_factor=1e-30); self.assertEqual(int((p > .5).sum()), 5)
        # an interior block needs about 70 nats beyond the ordinary switches to be called when it is shorter than the minimum, and none when it is as long as the minimum
        self.assertEqual(int((self.fb(self.evid(n, list(range(40, 60))), A=A) > .5).sum()), 20)        # a block of exactly the minimum (20 bins) is called as it is
        self.assertEqual(int((self.fb(self.evid(n, list(range(40, 44))), A=A) > .5).sum()), 0)

    def test_the_terminal_limit_works_with_the_interior_limit_off_at_both_ends(self):
        """The setting of the first real run. (The age dimension used to be capped at the interior limit, so this combination did nothing.)"""
        n = 80; kw = dict(A=1, young=1.0, Wt=5, term_factor=1e-30)
        for yh, expect in (([0, 1, 2], list(range(0, 5))), ([77, 78, 79], list(range(75, 80)))):
            self.assertEqual(np.where(self.fb(self.evid(n, yh), **kw) > .5)[0].tolist(), expect, 'strong reads: stretched to the minimum')
            self.assertEqual(np.where(self.fb(self.evid(n, yh), A=1, young=1.0) > .5)[0].tolist(), yh, 'no limits: exactly where the reads are')
        for lo, hi in ((0, 3), (n - 3, n)):
            weak = np.zeros((4, n)); weak[3, :] = 6; weak[3, lo:hi] = 1; weak[2, lo:hi] = 1
            self.assertEqual(int((self.fb(weak, A=1, young=1.0) > .5).sum()), 3); self.assertEqual(int((self.fb(weak, **kw) > .5).sum()), 0, 'weak reads: dropped')

    def test_the_two_ends_are_treated_alike_at_a_realistic_switch_rate(self):
        rng = np.random.default_rng(6)
        for kw in (dict(A=1, young=1.0, Wt=4, term_factor=1e-6), dict(A=6, young=1e-3, Wt=4, term_factor=1e-6), dict(A=6, young=1e-3, Wt=None)):
            ll0, ll1 = rng.normal(0, 2, 40), rng.normal(0, 2, 40)
            a = H.forward_backward(ll0, ll1, 1e-3, **kw); b = H.forward_backward(ll0[::-1], ll1[::-1], 1e-3, **kw)[::-1]
            self.assertLess(float(np.abs(a - b).max()), 0.02, kw)

    def test_short_interior_segments(self):
        post = np.array([.1] * 5 + [.9] * 3 + [.1] * 20 + [.9] * 30)
        segs = H.short_interior_segments(post, 100_000, 2.0)
        self.assertEqual([(s['start'], s['end'], s['state']) for s in segs], [(500_000, 800_000, 1)])      # the 20-bin segment is exactly 2 Mb: not short; the first and the last segment touch an end
        self.assertAlmostEqual(segs[0]['mean_posterior'], 0.9)
        self.assertEqual(len(H.short_interior_segments(post, 100_000, 3.0)), 2)
        self.assertEqual(H.short_interior_segments(np.array([.1] * 10 + [.9] * 10), 100_000, 5.0), [])      # a single switch makes no interior segment
        self.assertEqual(H.short_interior_segments(np.zeros(0), 100_000, 2.0), [])

    def test_bins_without_reads_give_probabilities_not_nans(self):
        c = np.zeros((4, 30)); c[3, 0:5] = 6; c[2, 25:30] = 3; c[3, 25:30] = 3
        p = self.fb(c); self.assertTrue(np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()); self.assertEqual(len(H.forward_backward(np.zeros(0), np.zeros(0), 1e-3, 5)), 0)

    def test_segment_limits_from_the_parameters(self):
        hp = dict(CFG['hmm']); self.assertEqual(hp['min_interior_segment_mb'], 0)                           # the first run measures the problem: no interior limit by default
        self.assertEqual(H.segment_limits(hp), (1, 2))
        hp['min_interior_segment_mb'] = 2.0; self.assertEqual(H.segment_limits(hp), (19, 2))                  # 2 Mb = 20 bins = A + 1; 0.2 Mb = 2 bins
        hp['min_terminal_segment_mb'] = 0; self.assertEqual(H.segment_limits(hp), (19, None))
        hp['min_interior_segment_mb'] = 0; self.assertEqual(H.segment_limits(hp), (1, None))

    def test_calls_and_breakpoints(self):
        hp = dict(CFG['hmm']); post = np.concatenate([np.full(3, .02), [.3, .7], np.full(20, .98), [.7, .3], np.full(4, .02)])
        bps = H.breakpoints(post, len(post) * 100_000, hp)
        self.assertEqual([(b['bin_lo'], b['bin_hi'], b['to_state']) for b in bps], [(300_000, 400_000, 1), (2_500_000, 2_600_000, 0)])
        self.assertEqual((bps[0]['interval_lo'], bps[0]['interval_hi']), (200_000, 600_000))        # the last bin confidently MC is bin 2; the first confidently YH is bin 5 (so the interval ends where bin 5 does)
        self.assertEqual((bps[0]['p_left'], bps[0]['p_right']), (0.3, 0.7))
        self.assertEqual((bps[0]['distance_to_end'], bps[1]['distance_to_end']), (400_000, 500_000))
        self.assertTrue(bps[0]['terminal_low_confidence'], 'a switch 0.4 Mb from the start is flagged low-confidence')
        self.assertFalse(bps[1]['terminal_low_confidence'], 'one exactly 0.5 Mb from the end is not (the limit is strict)')
        self.assertEqual(list(H.hard_calls([.2, .5, .51])), [0, 0, 1])

    def test_end_support(self):
        hp = dict(CFG['hmm']); c = np.zeros((4, 30)); c[3, 0:3] = 5; c[3, 28:30] = 1; post = np.concatenate([np.full(10, .02), np.full(10, .5), np.full(10, .98)])
        e = H.end_support(c, post, hp, end_mb=1.0)
        self.assertEqual((e['left_reads'], e['right_reads'], e['left_confident'], e['right_confident'], e['n_bins']), (15.0, 2.0, 1.0, 1.0, 30))

    def test_em_recovers_the_error_and_the_bias_and_the_genotypes(self):
        rng = np.random.default_rng(5); L = 40e6; chroms, truths = {}, {}
        for k in range(6):
            xs = SIM.draw_crossovers(rng, L); s0 = int(rng.integers(0, 2))
            p, m, y, r = SIM.simulate_reads(rng, L, xs, s0, 1.0, 1.0, False, True, eps=0.01, q=0.6, n_repeats=0)
            chroms[k] = (p, m, y, r, L); truths[k] = SIM.truth_bins(L, xs, s0)
        posts, cs, eps, q = H.genotype_fish(chroms, CFG['hmm'])
        self.assertTrue(0.005 < eps < 0.02, eps); self.assertTrue(0.54 < q < 0.68, q)
        acc = np.mean([(H.hard_calls(posts[k]) == truths[k]).mean() for k in range(6)]); self.assertGreater(acc, 0.98)


# ====================================================================== the simulation and the claims of the README
class TestSimulation(unittest.TestCase):
    def test_crossovers_respect_the_spacing_and_the_end_bias_moves_them_to_the_ends(self):
        rng = np.random.default_rng(1); allx = [SIM.draw_crossovers(rng, 40e6, True) for _ in range(1500)]
        self.assertTrue(all(all(b - a >= 8e6 for a, b in zip(x, x[1:])) for x in allx))
        near = lambda xs: np.mean([min(x, 40e6 - x) < 3e6 for g in xs for x in g])
        rng = np.random.default_rng(1); flat = [SIM.draw_crossovers(rng, 40e6, False) for _ in range(1500)]
        self.assertGreater(near(allx), near(flat) + 0.15)
        a = SIM.draw_crossovers(np.random.default_rng(9), 40e6, True); b = SIM.draw_crossovers(np.random.default_rng(9), 40e6, True); self.assertEqual(a, b)

    def test_truth_bins_flip_at_the_crossovers(self):
        t = SIM.truth_bins(1_000_000, [250_000, 600_000], 1); self.assertEqual(list(t), [1, 1, 1, 0, 0, 0, 1, 1, 1, 1])

    def test_scoring(self):
        truth = SIM.truth_bins(40e6, [1.0e6, 20.0e6], 0)
        called = truth.copy(); called[300:304] ^= 1                                      # a flip-flop of 400 kb in the interior
        det, spur, flip, term, acc = SIM.score(called, truth, 40e6)
        self.assertEqual([ok for _, ok in det], [True, True]); self.assertEqual((spur, flip, term), (2, 1, 0)); self.assertAlmostEqual(acc, 1 - 4 / 400)
        self.assertEqual(sorted(round(d / 1e6) for d, _ in det), [1, 20])
        missed = np.zeros_like(truth); det, spur, flip, term, acc = SIM.score(missed, truth, 40e6); self.assertEqual([ok for _, ok in det], [False, False])
        term_artifact = truth.copy(); term_artifact[-3:] ^= 1; self.assertEqual(SIM.score(term_artifact, truth, 40e6)[3], 1)

    def test_the_original_cannot_be_run_unchanged_on_a_reference_where_polarity_is_mixed(self):
        rng = np.random.default_rng(2); acc = {'as written': [], 'named': []}
        for _ in range(10):
            xs = SIM.draw_crossovers(rng, 40e6); s0 = int(rng.integers(0, 2)); t = SIM.truth_bins(40e6, xs, s0)
            f = SIM.simulate_reads(rng, 40e6, xs, s0, 1.0, 1.0, False, True, n_repeats=0)
            for name, kw in (('as written', {}), ('named', dict(oriented=True))):
                c = SIM.original_like(*f, 40e6, **kw); n = min(len(c), len(t)); acc[name].append((c[:n] == t[:n]).mean())
        self.assertLess(np.mean(acc['as written']), 0.8); self.assertGreater(np.mean(acc['named']), 0.95)

    def test_interior_limit_removes_flip_flops_and_keeps_the_ends_where_the_original_loses_them(self):
        rows = SIM.run_comparison(dict(CFG['hmm']), coverages=(0.5,), scenario='S2', n_fish=10, seed=21, which=['original (MC/YH alleles named)', 'bins, no segment limits', 'bins, interior limit'])
        d = pd.DataFrame([dict(m=r['method'], flip=r['flip'], det=r['det']) for r in rows])
        flips = d.groupby('m').flip.sum(); self.assertEqual(flips['bins, interior limit'], 0); self.assertGreater(flips['bins, no segment limits'], 0)
        near = lambda m: np.mean([ok for dets in d[d.m == m].det for dist, ok in dets if dist < 1e6])
        self.assertGreater(near('bins, interior limit'), near('original (MC/YH alleles named)') + 0.3)

    def test_the_named_variants_use_the_settings_in_their_names(self):
        hp = dict(CFG['hmm']); ms = SIM.methods(hp); rng = np.random.default_rng(5)
        for _ in range(120):
            xs = SIM.draw_crossovers(rng, 40e6, True); s0 = int(rng.integers(0, 2))
            fish = [(*SIM.simulate_reads(rng, 40e6, xs, s0, 0.5, 1.0, True, True, end_artifact=True), 40e6)]
            direct = lambda **o: H.hard_calls(H.genotype_fish({0: fish[0]}, dict(hp, **o))[0][0])
            none, first, interior, both = direct(min_interior_segment_mb=0, min_terminal_segment_mb=0), direct(min_interior_segment_mb=0, min_terminal_segment_mb=0.2), \
                direct(min_interior_segment_mb=2.0, min_terminal_segment_mb=0), direct(min_interior_segment_mb=2.0, min_terminal_segment_mb=0.2)
            if not np.array_equal(none, interior) and not np.array_equal(interior, both):          # a fish on which each limit changes the calls on its own
                got = lambda name: ms[name](fish)[0]
                np.testing.assert_array_equal(got('bins, no segment limits'), none); np.testing.assert_array_equal(got('bins, terminal limit only (the first run)'), first)
                np.testing.assert_array_equal(got('bins, interior limit'), interior); np.testing.assert_array_equal(got('bins, interior + terminal limits'), both)
                self.assertFalse(np.array_equal(got('bins, terminal limit only (the first run)'), both)); return
        self.fail('no simulated fish on which each limit changes the calls')

    def test_the_summary_tables(self):
        rows = SIM.run_comparison(dict(CFG['hmm']), coverages=(0.5,), scenario='S1', n_fish=2, seed=4, which=['bins, interior limit'])
        recall, false = SIM.summarize(rows)
        self.assertEqual(list(recall.columns), ['<0.5 Mb', '0.5-1 Mb', '1-2 Mb', '>=2 Mb']); self.assertIn('flip_flops_per_100_chr', false.columns)


# ====================================================================== the commands, end to end
import contextlib


def run_cmd(argv):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = GC.main(argv)
    return rc, buf.getvalue()


class TestLabelsAndFish(unittest.TestCase):
    def test_lg_labels_of_the_new_reference(self):
        self.assertEqual([GC.lg_label(c) for c in ('NC_135167.1', 'NC_135176.1', 'NC_135189.1', 'NC_135190.1', 'chr1', 'NW_0001.1')], ['LG1', 'LG10', 'LG23', 'NC_135190.1', 'chr1', 'NW_0001.1'])

    def test_only_fish_of_the_configured_dams_are_taken(self):
        t = pd.DataFrame({'SampleID': ['A', 'B', 'C', 'D', 'E', 'A'], 'BroodID': ['MCYH-BC1-B1', 'MCYH-BC1-B1', 'MCYH-BC1-B2', 'MCYH-BC1-B3', '', 'MCYH-BC1-B1'],
                          'Parents': ['MCYH-F1-003m,MC-3R6R-f', 'MCYH-F1-003m,MC-3R6R-f', 'MCYH-F1-002m,MC-1G2G-f', 'X,MC-OTHER-f', 'MCYH-F1-003m,MC-3R6R-f', 'MCYH-F1-003m,MC-3R6R-f']})
        f = GC.bc_fish(t, CFG)
        self.assertEqual(f.SampleID.tolist(), ['A', 'B', 'C']); self.assertEqual(f.Dam.tolist(), ['MC-3R6R-f', 'MC-3R6R-f', 'MC-1G2G-f']); self.assertEqual(f.Sire.iloc[2], 'MCYH-F1-002m')
        with self.assertRaises(SystemExit): GC.bc_fish(pd.DataFrame({'SampleID': ['A']}), CFG)

    def test_an_f1_brood_of_the_same_mothers_is_not_a_backcross(self):
        t = pd.DataFrame({'SampleID': ['A', 'B', 'F1a', 'F1b', 'Z'], 'BroodID': ['MCYH-BC1-24-01-10-1', 'MCYH-BC1-24-03-27-1', 'MCYH-F1-23-02-21', 'MCYH-F1-23-02-21', 'OTHER-1'],
                          'Parents': ['MCYHF1-002,MC-3R6R-f', 'MCYHF1-003,MC-5G11G-f', 'YH_008_m,MC-3R6R-f', 'YH_008_m,MC-3R6R-f', 'X,MC-OTHER-f']})
        self.assertEqual(GC.bc_fish(t, CFG).SampleID.tolist(), ['A', 'B'])
        lo = GC.left_out(t, CFG); self.assertEqual(lo.to_dict(), {('MCYH-F1-23-02-21', 'YH_008_m', 'MC-3R6R-f'): 2})        # fish of another dam are not listed: they are not ours at all
        cfg = {'sites': dict(CFG['sites'], brood_prefixes=[])}; self.assertEqual(GC.bc_fish(t, cfg).SampleID.tolist(), ['A', 'B', 'F1a', 'F1b']); self.assertEqual(len(GC.left_out(t, cfg)), 0)


class TestCommands(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp(); cls.root = cls.dir + '/G/'
        cls.vcf = os.path.join(cls.dir, 'bc_family.vcf.gz'); write_vcf(cls.vcf, planted())

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_a_sites_then_counts_plan_then_run_then_current(self):
        rc, text = run_cmd(['sites', '--vcf', self.vcf, '--out', self.root]); self.assertEqual(rc, 0)
        for frag in ('definition key', 'snp_biallelic', 'dam_opposite', 'MC-3R6R-f'): self.assertIn(frag, text)
        for dam in DAMS:
            self.assertTrue(os.path.exists(f'{self.root}sites/sites_{dam}.tsv.gz') and os.path.exists(f'{self.root}sites/targets_{dam}.tsv'))
        st = pd.DataFrame({'SampleID': ['F1', 'F2', 'F3', 'F4', 'F5', 'F6'], 'BroodID': ['MCYH-BC1-B1', 'MCYH-BC1-B1', 'MCYH-BC1-B2', 'MCYH-BC1-B9', 'MCYH-BC1-B1', 'MCYH-F1-X'],
                           'Parents': ['S,MC-3R6R-f', 'S,MC-3R6R-f', 'S,MC-1G2G-f', 'S,MC-OTHER-f', 'S,MC-3R6R-f', 'YH_008_m,MC-3R6R-f']})
        al = pd.DataFrame({'SampleID': ['F1', 'F2', 'F3', 'F4', 'F6'], 'BamSize': [1000, 1000, 1000, 1000, 1000], 'GenomeVersion': 'G'})      # F5 has no alignment row; F6 is an F1 fish
        stf, alf = self.dir + '/st.tsv', self.dir + '/al.tsv'; st.to_csv(stf, sep='\t', index=False); al.to_csv(alf, sep='\t', index=False)
        off = ['--sample-table', stf, '--alignment-table', alf, '--local-master', self.dir + '/local', '--cloud-master', self.dir + '/cloud', '--fasta', '/ref.fa', '--out', self.root, '--genome', 'G']
        rc, text = run_cmd(['counts'] + off); self.assertEqual(rc, 0)
        self.assertIn('3 BC fish with an alignment (1 without one)', text); self.assertIn('3 to count', text); self.assertIn('PLAN ONLY', text); self.assertIn('F5', text)
        self.assertRegex(text, r'MCYH-BC1-B1\s+S\s+MC-3R6R-f\s+2'); self.assertRegex(text, r'MCYH-BC1-B2\s+S\s+MC-1G2G-f\s+1')        # the plan shows brood, sire, dam and the number of fish
        self.assertIn('left out, not backcross broods', text); self.assertIn('MCYH-F1-X (sire YH_008_m, dam MC-3R6R-f): 1', text); self.assertNotIn('F6', text.split('left out')[0])
        self.assertFalse(os.path.exists(self.root + 'counts'))                                    # a plan writes nothing
        fake = Fake()
        saved = dict(C.RUNNERS); C.RUNNERS.update(run=fake.run, script=fake.script)
        try:
            rc, text = run_cmd(['counts', '--run', '--jobs', '2'] + off)
            self.assertEqual(rc, 0); self.assertEqual(text.count(' reads at '), 3); self.assertEqual(fake.calls.count('pileup'), 3)
            rc, text = run_cmd(['counts', '--run'] + off); self.assertIn('0 to count', text); self.assertIn('3 current', text); self.assertEqual(fake.calls.count('pileup'), 3)     # nothing is counted twice
            rc, text = run_cmd(['counts', '--run', '--redo', '--samples', 'F1'] + off); self.assertEqual(fake.calls.count('pileup'), 4)
            fake.fail = True; rc, text = run_cmd(['counts', '--run', '--redo', '--samples', 'F2'] + off); self.assertEqual(rc, 3); self.assertIn('FAILED', text)
        finally:
            C.RUNNERS.update(saved)

    def test_b_counts_without_sites_says_so(self):
        d = tempfile.mkdtemp()
        try:
            with self.assertRaises(SystemExit) as cm: run_cmd(['counts', '--sample-table', 'x', '--alignment-table', 'y', '--local-master', d, '--cloud-master', d, '--fasta', '/r.fa', '--out', d + '/G'])
            self.assertIn('genotype sites', str(cm.exception))
        finally:
            shutil.rmtree(d, ignore_errors=True)


class TestHmmCommand(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp(); cls.root = cls.dir + '/G/'
        cls.L = 20e6; cls.chroms = ['NC_135167.1', 'NC_135168.1', 'NC_135176.1']
        fai = cls.dir + '/ref.fa.fai'
        open(fai, 'w').write(''.join(f'{c}\t{int(cls.L)}\t0\t60\t61\n' for c in cls.chroms) + 'NW_0001.1\t1000000\t0\t60\t61\n'); cls.fai = fai
        rng = np.random.default_rng(8); cls.truth = {}
        os.makedirs(cls.root + 'counts')
        for f in ('FISH1', 'FISH2', 'FISH3'):
            rows = []
            for ch in cls.chroms:
                xs = [rng.uniform(3e6, 17e6)] if f != 'FISH3' else [1.3e6]                       # FISH3 has a crossover 1.3 Mb from the start
                s0 = int(rng.integers(0, 2)); p, m, y, r = SIM.simulate_reads(rng, cls.L, xs, s0, 1.0, 1.0, False, True, n_repeats=0)
                cls.truth[(f, ch)] = (xs, SIM.truth_bins(cls.L, xs, s0))
                rows += [(ch, int(a), int(b), int(c), int(d)) for a, b, c, d in zip(p, m, y, r)]
            pd.DataFrame(rows, columns=['chrom', 'pos', 'n_mc', 'n_yh', 'yh_is_ref']).to_csv(f'{cls.root}counts/{f}.counts.tsv.gz', sep='\t', index=False, compression='gzip')
            json.dump(dict(sample=f, status='ok', key='K'), open(f'{cls.root}counts/{f}.counts.json', 'w'))
        json.dump(dict(sample='BAD', status='failed', key='K'), open(cls.root + 'counts/BAD.counts.json', 'w'))
        cls.rc, cls.text = run_cmd(['hmm', '--out', cls.root, '--fai', fai, '--jobs', '1'])
        cls.o = cls.root + 'hmm/'

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_outputs_and_the_grid_cover_the_whole_chromosomes(self):
        self.assertEqual(self.rc, 0); self.assertIn('3 fish, 3 chromosomes', self.text)                      # the failed fish and the unplaced contig are not used
        calls = pd.read_csv(self.o + 'genotypes_calls.csv', index_col=0)
        self.assertEqual(list(calls.columns), ['Chromosome', 'Position', 'FISH1', 'FISH2', 'FISH3']); self.assertEqual(len(calls), 3 * 200)
        self.assertEqual(calls.groupby('Chromosome').Position.max().tolist(), [19_900_000] * 3)               # nothing is cut at the shortest fish: the grid runs to the end of every chromosome
        post = pd.read_csv(self.o + 'genotypes_posterior.csv.gz', index_col=0); self.assertEqual(post.shape, calls.shape); self.assertTrue(((post.iloc[:, 2:] >= 0) & (post.iloc[:, 2:] <= 1)).all().all())

    def test_the_calls_agree_with_the_truth(self):
        calls = pd.read_csv(self.o + 'genotypes_calls.csv', index_col=0)
        for f in ('FISH1', 'FISH2', 'FISH3'):
            for ch in self.chroms:
                c = calls[calls.Chromosome == ch][f].to_numpy(); t = self.truth[(f, ch)][1]; self.assertGreater((c == t).mean(), 0.97, (f, ch))

    def test_breakpoints_end_support_and_qc(self):
        bp = pd.read_csv(self.o + 'breakpoints.tsv', sep='\t')
        self.assertEqual(set(bp.columns) >= {'sample', 'chrom', 'LG', 'bin_lo', 'interval_lo', 'interval_hi', 'p_left', 'p_right', 'to_state', 'distance_to_end', 'terminal_low_confidence'}, True)
        for f in ('FISH1', 'FISH2', 'FISH3'):
            for ch in self.chroms:
                x = self.truth[(f, ch)][0][0]; b = bp[(bp['sample'] == f) & (bp.chrom == ch)]
                self.assertTrue(((b.bin_hi - x).abs() < 1.0e6).any(), (f, ch, x))                                  # every planted crossover is found within 1 Mb
        self.assertEqual(set(bp.LG), {'LG1', 'LG2', 'LG10'})
        f3 = bp[(bp['sample'] == 'FISH3') & (bp.chrom == 'NC_135167.1')]; self.assertTrue(1.0e6 <= f3.distance_to_end.min() <= 1.6e6, f3.distance_to_end.min())        # the crossover 1.3 Mb from the start
        es = pd.read_csv(self.o + 'end_support.tsv', sep='\t'); self.assertEqual(len(es), 9); self.assertTrue((es.left_reads > 0).all() and (es.right_reads > 0).all())
        qc = pd.read_csv(self.o + 'sample_qc.tsv', sep='\t'); self.assertEqual(len(qc), 3); self.assertFalse(qc.flag.any()); self.assertTrue(((qc.eps > 0.002) & (qc.eps < 0.05)).all())
        self.assertEqual(json.load(open(self.o + 'run_info.json'))['n_fish'], 3)

    def test_a_tag_gives_another_folder_and_a_changed_limit_is_recorded(self):
        rc, text = run_cmd(['hmm', '--out', self.root, '--fai', self.fai, '--jobs', '1', '--tag', '_free', '--interior-min-mb', '0', '--terminal-min-mb', '0', '--rho', '2e-8'])
        self.assertEqual(rc, 0); info = json.load(open(self.root + 'hmm_free/run_info.json'))['parameters']
        self.assertEqual((info['min_interior_segment_mb'], info['min_terminal_segment_mb'], info['rho']), (0, 0, 2e-8)); self.assertEqual(json.load(open(self.o + 'run_info.json'))['parameters']['rho'], 1e-8)

    def test_missing_inputs_are_named(self):
        d = tempfile.mkdtemp()
        try:
            with self.assertRaises(SystemExit) as cm: run_cmd(['hmm', '--out', d + '/G', '--fai', self.fai])
            self.assertIn('no counts', str(cm.exception))
            os.makedirs(d + '/G/counts')
            with self.assertRaises(SystemExit) as cm: run_cmd(['hmm', '--out', d + '/G', '--fai', self.fai])
            self.assertIn('no counts', str(cm.exception))                                         # a folder without a finished fish
        finally:
            shutil.rmtree(d, ignore_errors=True)
        with self.assertRaises(SystemExit) as cm: run_cmd(['hmm', '--out', self.root, '--fai', self.dir + '/nope.fai'])
        self.assertIn('reference index', str(cm.exception))


class TestLikelihoodAndPolarity(unittest.TestCase):
    def test_the_likelihood_of_the_four_read_classes(self):
        """YH reads / MC reads on sites where the YH allele is the reference (a) and on sites where it is the alternative allele (b). In a heterozygote a reference-allele read has probability q."""
        c = np.array([[2.0], [3.0], [5.0], [7.0]])                  # yh_a, mc_a, yh_b, mc_b
        l0, l1 = H.bin_loglik(c, 0.01, 0.6)
        self.assertAlmostEqual(float(l0[0]), 7 * np.log(0.01) + 10 * np.log(0.99))
        self.assertAlmostEqual(float(l1[0]), 2 * np.log(0.6) + 3 * np.log(0.4) + 5 * np.log(0.4) + 7 * np.log(0.6))
        one = lambda k: H.bin_loglik(np.eye(4)[:, [k]], 0.01, 0.6)[1][0]
        self.assertGreater(one(0), one(2)); self.assertGreater(one(3), one(1))            # a YH read is likelier in a heterozygote on an `a` site; an MC read on a `b` site

    def reads(self, rng, ref_polarity, q=0.65, L=40e6, n_chr=5):
        chroms, truth = {}, {}
        for k in range(n_chr):
            xs = SIM.draw_crossovers(rng, L); s0 = int(rng.integers(0, 2)); t = SIM.truth_bins(L, xs, s0)
            pos = np.sort(rng.integers(1000, int(L), 40000)); st = t[np.minimum(pos // 100_000, len(t) - 1)]
            p_het = q if ref_polarity else 1 - q
            n = rng.poisson(1.0, len(pos)); y = rng.binomial(n, np.where(st == 1, p_het, 0.01)); keep = n > 0
            chroms[k] = (pos[keep], (n - y)[keep], y[keep], np.full(int(keep.sum()), ref_polarity), L); truth[k] = t
        return chroms, truth

    def test_the_bias_is_recovered_and_the_calls_are_right_whichever_polarity_all_the_sites_have(self):
        for ref_polarity in (True, False):
            chroms, truth = self.reads(np.random.default_rng(12), ref_polarity)
            posts, _, eps, q = H.genotype_fish(chroms, CFG['hmm'])
            self.assertTrue(0.57 < q < 0.7, (ref_polarity, q)); self.assertGreater(np.mean([(H.hard_calls(posts[k]) == truth[k]).mean() for k in truth]), 0.97, ref_polarity)


def write_counts(root, fish, chroms, L, extra=None):
    os.makedirs(root + 'counts', exist_ok=True)
    for f in fish:
        rows = []
        for ch in chroms:
            rows += (extra or (lambda f, ch: None))(f, ch) or [(ch, int(p), 2, 0, 0) for p in np.arange(1000, L, 1000)]
        pd.DataFrame(rows, columns=['chrom', 'pos', 'n_mc', 'n_yh', 'yh_is_ref']).to_csv(f'{root}counts/{f}.counts.tsv.gz', sep='\t', index=False, compression='gzip')
        json.dump(dict(sample=f, status='ok', key='K'), open(f'{root}counts/{f}.counts.json', 'w'))


class TestShortSegmentsAndCompare(unittest.TestCase):
    """The setting of the first real run (no interior limit, terminal limit 0.2 Mb) lists the short interior segments; the 2 Mb sensitivity run removes them; compare says by how much."""
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp(); cls.root = cls.dir + '/G/'; L = 20_000_000; chroms = ['NC_135167.1', 'NC_135168.1', 'NC_135176.1']
        cls.fai = cls.dir + '/ref.fa.fai'; open(cls.fai, 'w').write(''.join(f'{c}\t{L}\t0\t60\t61\n' for c in chroms))

        def extra(f, ch):
            if f == 'ART' and ch == 'NC_135168.1':                                  # MC reads everywhere except 10.0-10.4 Mb, where four sites carry two YH reads each (an artifact of 400 kb)
                pos = [p for p in range(1000, L, 1000) if not (10_000_000 <= p < 10_400_000)]
                return [(ch, p, 2, 0, 0) for p in pos] + [(ch, p, 0, 2, 0) for p in (10_050_000, 10_150_000, 10_250_000, 10_350_000)]
            return None
        write_counts(cls.root, ['ART', 'CLEAN'], chroms, L, extra)
        cls.rc1, cls.t1 = run_cmd(['hmm', '--out', cls.root, '--fai', cls.fai, '--jobs', '1'])
        cls.rc2, cls.t2 = run_cmd(['hmm', '--out', cls.root, '--fai', cls.fai, '--jobs', '1', '--interior-min-mb', '2', '--tag', '_int2'])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_the_first_run_lists_the_artifact_and_says_how_large_the_problem_is(self):
        self.assertEqual((self.rc1, self.rc2), (0, 0)); self.assertIn('interior minimum 0 Mb, terminal minimum 0.2 Mb', self.t1)
        sh = pd.read_csv(self.root + 'hmm/short_segments.tsv', sep='\t')
        self.assertEqual(len(sh), 1); r = sh.iloc[0]
        self.assertEqual((r['sample'], r.chrom, r.LG, r.start, r.end, r.length, r.state), ('ART', 'NC_135168.1', 'LG2', 10_000_000, 10_400_000, 400_000, 1)); self.assertGreater(r.mean_posterior, 0.5)
        self.assertIn('1 interior segments shorter than 2.0 Mb in 1 of 2 fish', self.t1)
        qc = pd.read_csv(self.root + 'hmm/sample_qc.tsv', sep='\t').set_index('sample'); self.assertEqual((int(qc.loc['ART', 'short_interior_segments']), int(qc.loc['CLEAN', 'short_interior_segments'])), (1, 0))

    def test_the_two_megabase_minimum_removes_it(self):
        sh = pd.read_csv(self.root + 'hmm_int2/short_segments.tsv', sep='\t'); self.assertEqual(len(sh), 0)
        self.assertEqual(json.load(open(self.root + 'hmm_int2/run_info.json'))['parameters']['min_interior_segment_mb'], 2.0)
        self.assertEqual(json.load(open(self.root + 'hmm/run_info.json'))['parameters']['min_interior_segment_mb'], 0)

    def test_compare_reports_the_difference(self):
        rc, text = run_cmd(['compare', '--out', self.root, '--b', 'hmm_int2']); self.assertEqual(rc, 0)
        for frag in ('2 fish', 'bins called differently: 4', 'in 1 fish', 'breakpoints: 2 -> 0', 'short interior segments: 1 -> 0'): self.assertIn(frag, text)
        t = pd.read_csv(self.root + 'compare_hmm_vs_hmm_int2.tsv', sep='\t').set_index('sample'); self.assertEqual((int(t.loc['ART', 'bins_differ']), int(t.loc['CLEAN', 'bins_differ'])), (4, 0))
        self.assertEqual((int(t.loc['ART', 'breakpoints_a']), int(t.loc['ART', 'breakpoints_b'])), (2, 0))

    def test_compare_asks_for_what_it_needs(self):
        with self.assertRaises(SystemExit) as cm: run_cmd(['compare', '--out', self.root])
        self.assertIn('--b', str(cm.exception))
        with self.assertRaises(SystemExit) as cm: run_cmd(['compare', '--out', self.root, '--b', 'hmm_nope'])
        self.assertIn('genotypes_calls.csv is missing', str(cm.exception))


class TestSimulateCommand(unittest.TestCase):
    def test_simulate_prints_the_tables_and_writes_them(self):
        d = tempfile.mkdtemp()
        try:
            rc, text = run_cmd(['simulate', '--fish', '1', '--coverage', '0.5', '--scenario', 'S1', '--out', d]); self.assertEqual(rc, 0)
            for frag in ('recall of the true crossovers', 'flip_flops_per_100_chr', 'old-genome-like'): self.assertIn(frag, text)
            self.assertTrue(os.path.exists(d + '/simulation_S1_recall.tsv') and os.path.exists(d + '/simulation_S1_false_calls.tsv'))
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_the_command_is_reachable_from_bc_analysis(self):
        from bc_analysis import cli as top
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = top.main(['genotype', 'simulate', '--fish', '1', '--coverage', '0.5', '--scenario', 'S1'])
        self.assertEqual(rc, 0); self.assertIn('recall of the true crossovers', buf.getvalue())


if __name__ == '__main__':
    unittest.main()
