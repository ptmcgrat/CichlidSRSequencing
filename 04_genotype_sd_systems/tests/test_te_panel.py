"""Transposon-insertion panel module (te): vendored code pinned to sexloci, definition build, classifier equivalence, measurement on simulated read pairs (fake samtools), calls, flags, panel registry."""
import hashlib, json, math, os, sys, tempfile, types, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import yaml
from cichlid_common.common_light import config_path
from sd_genotypes.te import call as C, definition as D, longindel as L, measure as M, sites as K
from sd_genotypes.panel import measure as pm

CFG0 = yaml.safe_load(open(config_path('tepanel.yaml'))); CFG0['backend'] = 'numpy'
LENGTHS = {'NC_135176.1': 34_949_860, 'NC_135173.1': 69_537_369}
ANCHORS = {('NC_135176.1', 24921177): 'T', ('NC_135176.1', 24921393): 'A', ('NC_135173.1', 50765049): 'G', ('NC_135173.1', 50786786): 'G'}      # 0-based positions of the anchor bases


def base(chrom, i):
    return ANCHORS.get((chrom, i)) or 'ACGT'[int(hashlib.md5(f'{chrom}:{i}'.encode()).hexdigest()[:4], 16) % 4]


def fetch(chrom, a, b):
    return ''.join(base(chrom, i) for i in range(max(a, 0), b))


def make_def():
    tmp = tempfile.mkdtemp(); D.build(fetch, LENGTHS, tmp); cfg = dict(CFG0); cfg['definition_dir'] = tmp
    return cfg, tmp


rc = L.revcomp


def sam(name, flag, seq, chrom='NC_135176.1', pos=24921000):
    return '\t'.join([name, str(flag), chrom if not flag & 4 else '*' if False else chrom, str(pos), '60', '100M', '=', str(pos + 200), '300', seq, 'I' * len(seq)])


def pairs_for(site, allele, bps, n, tag, mate_unmapped=False):
    """n read pairs from an allele: mate 1 crosses a breakpoint (50 bp each side), mate 2 lies in the left flank (the same on both alleles)."""
    seq = allele; out = []
    for k in range(n):
        bp = bps[k % len(bps)]; m1 = seq[bp - 50:bp + 50]; m2 = rc(seq[60 + k:160 + k])
        out.append(sam(f'{tag}{k}', 99 | (4 if mate_unmapped else 0), m1, site['chrom'])); out.append(sam(f'{tag}{k}', 147, m2, site['chrom']))
    return out


def site_reads(site, n_ref, n_alt, unmapped_alt=False):
    left, right, ref, alt = site['left'], site['right'], site['ref'], site['alt']
    ref_allele, alt_allele = left + ref + right, left + alt + right
    lines = pairs_for(site, ref_allele, [len(left) + 1], n_ref, f"{site['name']}r")
    lines += pairs_for(site, alt_allele, [len(left) + 1, len(left) + len(alt)], n_alt, f"{site['name']}a", mate_unmapped=unmapped_alt)
    return lines


def fake_run(cfg, plan, lengths=LENGTHS, log=None, fail_site=None):
    """plan: {site name: (n_ref pairs, n_alt pairs)}. samtools view over a window returns the reads of every site whose anchor is inside it."""
    d = M.load_definition(cfg); by = {s['name']: s for s in d['sites']}

    def run(cmd):
        if log is not None:
            log.append(cmd)
        if cmd[:2] == ['samtools', '--version']:
            return types.SimpleNamespace(stdout='samtools 1.x\n')
        if cmd[:3] == ['samtools', 'view', '-H']:
            return types.SimpleNamespace(stdout='@HD\tVN:1.6\n' + ''.join(f'@SQ\tSN:{c}\tLN:{n}\n' for c, n in lengths.items()))
        if cmd[:2] == ['samtools', 'view']:
            chrom, rng = cmd[-1].split(':'); a, b = map(int, rng.split('-'))
            if fail_site and chrom == by[fail_site]['chrom'] and a <= by[fail_site]['pos'] <= b:
                raise RuntimeError('samtools view failed: truncated file')
            lines = []
            for n, v in plan.items():
                s = by[n]; nr, na = v[0], v[1]
                if s['chrom'] == chrom and a <= s['pos'] <= b:
                    lines += site_reads(s, nr, na, unmapped_alt=len(v) > 2)
            return types.SimpleNamespace(stdout='\n'.join(lines) + '\n')
        raise AssertionError(cmd)
    return run


class TestVendored(unittest.TestCase):
    def body(self, f):
        return open(os.path.join(os.path.dirname(M.__file__), f)).read().split('\n', 2)[2]

    def test_the_vendored_files_are_byte_identical_to_sexloci_0_1_26(self):
        self.assertEqual(hashlib.md5(self.body('longindel.py').encode()).hexdigest(), '0ee7541e51610bb101a811ae7e8ac12d')
        self.assertEqual(hashlib.md5(self.body('sites.py').encode()).hexdigest(), '89ad7175765a70a14b9f9cacb4951146')

    def test_allele_fingerprints_are_those_of_VARIANTS_md_appendix_A(self):
        fp = {n: K.fingerprint(a) for g in K.BUILTIN.values() for c, p, n, r, a, i in g}
        self.assertEqual(fp, D.EXPECTED_FINGERPRINTS)
        pos = {n: (c, p, r) for g in K.BUILTIN.values() for c, p, n, r, a, i in g}
        self.assertEqual(pos, {'TE1': ('NC_135176.1', 24921178, 't'), 'TE2': ('NC_135176.1', 24921394, 'a'), 'chr7ins_c1': ('NC_135173.1', 50765050, 'G'), 'chr7ins_c2': ('NC_135173.1', 50786787, 'G')})
        self.assertEqual(K.SITE_GROUPS, {'chr7ins': ('chr7ins_c1', 'chr7ins_c2')})

    def test_the_call_model_is_the_documented_binomial(self):
        def ref_gt(nr, na):
            if nr + na < 5:
                return './.'
            ll = {g: na * math.log(p) + nr * math.log(1 - p) for g, p in (('0/0', 0.02), ('0/1', 0.5), ('1/1', 0.98))}
            return max(ll, key=ll.get)
        for nr in range(0, 21):
            for na in range(0, 21):
                self.assertEqual(L.call_from_counts(nr, na)[0], ref_gt(nr, na), (nr, na))
        gt, gq, _ = L.call_from_counts(4, 4); self.assertEqual(gt, '0/1'); self.assertGreater(gq, 40)


class TestDefinition(unittest.TestCase):
    def test_build_stores_records_flanks_groups_and_fingerprints(self):
        cfg, tmp = make_def(); d = M.load_definition(cfg)
        self.assertEqual([s['name'] for s in d['sites']], ['TE1', 'TE2', 'chr7ins_c1', 'chr7ins_c2'])
        for s in d['sites']:
            self.assertEqual((len(s['left']), len(s['right'])), (250, 250)); self.assertEqual(s['ref'], base(s['chrom'], s['pos'] - 1).upper())
        self.assertEqual(d['groups'], {'chr7ins': ['chr7ins_c1', 'chr7ins_c2']}); self.assertEqual(d['lengths'], LENGTHS)
        self.assertEqual(d['sites'][1]['fingerprint'], 'ALT 9,869 bp, md5 d2b75a3b30e2')

    def test_a_wrong_reference_or_a_contig_edge_is_refused(self):
        bad = lambda c, a, b: ''.join('C' for _ in range(max(a, 0), b))
        with self.assertRaises(ValueError):
            D.build(bad, LENGTHS, tempfile.mkdtemp())                                       # TE1's anchor base is t, the fake reference has C
        short = lambda c, a, b: fetch(c, a, b) if a > 0 else fetch(c, 0, b)[:10]
        with self.assertRaises(ValueError):
            D.build(lambda c, a, b: fetch(c, a, b)[:100] if b - a == 250 else fetch(c, a, b), LENGTHS, tempfile.mkdtemp())
        with self.assertRaises(ValueError):
            D.build(fetch, {'NC_135176.1': 1}, tempfile.mkdtemp())

    def test_a_changed_allele_fails_the_pinned_fingerprint(self):
        orig = K.BUILTIN['chr7_gsdf_ins']; K.BUILTIN['chr7_gsdf_ins'] = [(c, p, n, r, a[:-1] + ('A' if a[-1] != 'A' else 'C'), i) for c, p, n, r, a, i in orig]
        try:
            with self.assertRaises(ValueError):
                D.build(fetch, LENGTHS, tempfile.mkdtemp())
        finally:
            K.BUILTIN['chr7_gsdf_ins'] = orig

    def test_the_classifier_from_stored_flanks_equals_the_one_from_the_reference(self):
        cfg, _ = make_def(); d = M.load_definition(cfg)

        class Ref:
            def __getitem__(self, chrom):
                return types.SimpleNamespace(__getitem__=None) if False else _W(chrom)

        class _W:
            def __init__(self, chrom): self.c = chrom
            def __getitem__(self, sl): return fetch(self.c, sl.start, sl.stop)
        be = L.make_backend('numpy', 5, 2)
        for s in d['sites']:
            a = M.classifier_from_site(s, 250, cfg['classifier'], be)
            b = L.IndelReadClassifier.from_vcf_record(Ref(), s['chrom'], s['pos'], s['ref'], s['alt'], flanking=250, backend=be, **cfg['classifier'])
            self.assertEqual((a.ref_allele, a.alt_allele, a.ref_breakpoints, a.alt_breakpoints, a.fetch_window(), a.insertion_pos), (b.ref_allele, b.alt_allele, b.ref_breakpoints, b.alt_breakpoints, b.fetch_window(), b.insertion_pos), s['name'])

    def test_regions_are_the_merged_read_windows(self):
        cfg, _ = make_def()
        # window = anchor +- (250 + 1 + 200): TE1 24,920,728-24,921,629, TE2 24,920,943-24,921,845 (216 bp apart: merged); the two chr7ins copies are 21.7 kb apart
        self.assertEqual(M.regions(cfg), ['NC_135173.1:50764600-50765501', 'NC_135173.1:50786337-50787238', 'NC_135176.1:24920728-24921845'])      # sorted by sequence name; TE1 and TE2 merged


class TestMeasureAndCall(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg, cls.tmp = make_def(); cls.d = M.load_definition(cls.cfg)

    def measure(self, plan, **kw):
        return M.measure_bam('x.bam', self.cfg, fake_run(self.cfg, plan, **kw))

    def test_pair_counts_per_site(self):
        out = self.measure({'TE1': (6, 0), 'chr7ins_c1': (3, 3), 'chr7ins_c2': (2, 4)})
        self.assertEqual(out['backend'], 'numpy'); self.assertEqual(out['definition'], M.definition_key(self.cfg))
        s = out['sites']
        self.assertEqual((s['TE1']['n_ref'], s['TE1']['n_alt'], s['TE1']['status']), (6, 0, 'OK'))
        self.assertEqual((s['chr7ins_c1']['n_ref'], s['chr7ins_c1']['n_alt']), (3, 3)); self.assertEqual((s['chr7ins_c2']['n_ref'], s['chr7ins_c2']['n_alt']), (2, 4))
        self.assertEqual((s['TE2']['n_ref'], s['TE2']['n_alt']), (0, 0))                    # no reads of TE2 in the plan; TE1's reads are not TE2 evidence

    def test_the_chr10_insertions_TE1_and_TE2_are_classified_at_their_own_site_only(self):
        """The two LG10 insertions are 216 bp apart: each window contains the other's anchor, so reads of one are presented to the other's classifier too and must not count there."""
        out = self.measure({'TE1': (0, 6), 'TE2': (3, 3)})                                   # TE1: six insertion pairs; TE2 (the 9,869 bp allele): three reference and three insertion pairs
        s = out['sites']
        self.assertEqual((s['TE1']['n_ref'], s['TE1']['n_alt']), (0, 6)); self.assertEqual((s['TE2']['n_ref'], s['TE2']['n_alt']), (3, 3))
        t = C.call_table([dict(SampleID='x', **out)], self.cfg).set_index('SampleID')
        self.assertEqual((t.loc['x', 'TE1_GT'], t.loc['x', 'TE2_GT']), ('1/1', '0/1'))        # TE1 homozygous for the insertion, TE2 heterozygous
        both = self.measure({'TE1': (4, 4), 'TE2': (0, 6)})['sites']                          # and the other way round
        self.assertEqual(((both['TE1']['n_ref'], both['TE1']['n_alt']), (both['TE2']['n_ref'], both['TE2']['n_alt'])), ((4, 4), (0, 6)))

    def test_unmapped_mates_and_duplicates(self):
        s = self.d['sites'][2]; lines = site_reads(s, 0, 4, unmapped_alt=True)
        pairs = L.reads_from_samtools_lines(lines)
        self.assertEqual(len(pairs), 8)                                                     # unmapped mates are kept
        dup = [l.split('\t') for l in lines]; dup[0][1] = str(int(dup[0][1]) | 1024)
        self.assertEqual(len(L.reads_from_samtools_lines(['\t'.join(x) for x in dup])), 7)   # a duplicate-flagged read is dropped

    def test_pairs_whose_informative_mate_is_unmapped_count_only_when_unmapped_reads_are_included(self):
        plan = {'chr7ins_c1': (0, 5, 'u'), 'chr7ins_c2': (0, 0)}                           # 5 insertion pairs whose junction mate is stored unmapped at its partner's position
        self.assertEqual(self.measure(plan)['sites']['chr7ins_c1']['n_alt'], 5)
        c = dict(self.cfg); c['include_unmapped'] = False
        out = M.measure_bam('x.bam', c, fake_run(c, plan))
        self.assertEqual(out['sites']['chr7ins_c1']['n_alt'], 0)

    def test_a_site_that_fails_is_recorded_and_the_others_go_on(self):
        out = self.measure({'chr7ins_c1': (0, 5), 'chr7ins_c2': (0, 5)}, fail_site='TE1')
        self.assertEqual(out['sites']['TE1']['status'], 'error'); self.assertIn('truncated', out['sites']['TE1']['message'])
        self.assertEqual(out['sites']['TE2']['status'], 'error')                            # the TE2 window contains the TE1 anchor too: same request fails
        self.assertEqual((out['sites']['chr7ins_c1']['status'], out['sites']['chr7ins_c1']['n_alt']), ('OK', 5))

    def test_a_bam_of_another_reference_is_refused(self):
        with self.assertRaises(ValueError):
            self.measure({}, lengths=dict(LENGTHS, **{'NC_135176.1': 34_000_000}))
        with self.assertRaises(ValueError):
            self.measure({}, lengths={'NC_135173.1': LENGTHS['NC_135173.1']})

    def test_the_definition_key_follows_what_is_measured(self):
        k = M.definition_key(self.cfg)
        c2 = dict(self.cfg); c2['classifier'] = dict(c2['classifier'], score_margin=30); self.assertNotEqual(M.definition_key(c2), k)
        c3 = dict(self.cfg); c3['call'] = dict(c3['call'], min_gq=99); c3['backend'] = 'parasail'; self.assertEqual(M.definition_key(c3), k)    # call settings and the backend are not part of it
        tmp2 = tempfile.mkdtemp(); d = json.load(open(os.path.join(self.tmp, 'te_definition.json'))); d['sites'][0]['left'] = 'A' + d['sites'][0]['left'][1:]
        json.dump(d, open(os.path.join(tmp2, 'te_definition.json'), 'w')); c4 = dict(self.cfg); c4['definition_dir'] = tmp2
        self.assertNotEqual(M.definition_key(c4), k)                                        # a different flank is a different definition

    def rows(self, plans):
        return [dict(SampleID=sid, **self.measure(p)) for sid, p in plans.items()]

    def test_genotypes_and_the_combined_chr7ins_call(self):
        rows = self.rows({'hom_ref': {'chr7ins_c1': (4, 0), 'chr7ins_c2': (4, 0)}, 'het': {'chr7ins_c1': (3, 3), 'chr7ins_c2': (3, 3)}, 'hom_alt': {'chr7ins_c1': (0, 5), 'chr7ins_c2': (0, 5)},
                          'low': {'chr7ins_c1': (1, 1), 'chr7ins_c2': (1, 1)}, 'split': {'chr7ins_c1': (2, 1), 'chr7ins_c2': (2, 2)}})
        t = C.call_table(rows, self.cfg).set_index('SampleID')
        self.assertEqual([t.loc[s, 'chr7ins_GT'] for s in ('hom_ref', 'het', 'hom_alt', 'low')], ['0/0', '0/1', '1/1', './.'])
        self.assertEqual((t.loc['het', 'chr7ins_ref'], t.loc['het', 'chr7ins_alt']), (6, 6))   # the informative pairs of both copies are summed
        self.assertEqual(t.loc['split', 'chr7ins_GT'], C.call_from_counts(4, 3)[0])            # 4 ref + 3 alt in total
        self.assertEqual(t.loc['low', 'chr7ins_call'], ''); self.assertEqual(t.loc['het', 'chr7ins_call'], '0/1')   # the sheet value: confident calls only
        self.assertEqual(t.loc['hom_ref', 'TE1_GT'], './.')                                    # no TE1 reads at all: no call
        self.assertIn('chr7ins: 4 confident calls (GQ >= 10) of 5 measured: 0/0 1, 0/1 2, 1/1 1', C.summary(t.reset_index(), self.cfg))      # hom_ref, het, split (4 + 3 pairs) and hom_alt; `low` has only 4 informative pairs

    def test_the_dilution_flag_needs_a_gsdf_copy_number_above_two_and_a_0_0_with_insertion_pairs(self):
        rows = self.rows({'a': {'chr7ins_c1': (8, 1), 'chr7ins_c2': (8, 0)}, 'b': {'chr7ins_c1': (8, 1), 'chr7ins_c2': (8, 0)}, 'c': {'chr7ins_c1': (8, 0), 'chr7ins_c2': (8, 0)}, 'd': {'chr7ins_c1': (3, 3), 'chr7ins_c2': (3, 3)}})
        t = C.call_table(rows, self.cfg, copies=pd.Series({'a': 3, 'b': 2, 'c': 4, 'd': 4})).set_index('SampleID')
        self.assertEqual(t.chr7ins_flags.to_dict(), {'a': 'dilution_possible', 'b': '', 'c': '', 'd': ''})   # a: 3 copies, 0/0 with an insertion pair; c has no insertion pair; d is called 0/1
        self.assertEqual(C.call_table(rows, self.cfg).chr7ins_flags.unique().tolist(), [''])    # without copy numbers no flag is made

    def test_definitions_are_not_mixed_and_a_sample_twice_is_refused(self):
        rows = self.rows({'a': {'chr7ins_c1': (4, 0), 'chr7ins_c2': (4, 0)}}); bad = [dict(rows[0], definition='other')]
        with self.assertRaises(ValueError):
            C.call_table(bad, self.cfg)
        with self.assertRaises(ValueError):
            C.call_table(rows + rows, self.cfg)

    def test_the_panel_runs_it_and_the_json_round_trips_and_the_cli_calls_a_folder(self):
        self.assertIn('te', pm.MODULES); self.assertIsInstance(pm.module_cfg('te'), dict)
        d = tempfile.mkdtemp(); bam = d + '/s.all.bam'; open(bam, 'wb').write(b'x'); open(bam + '.bai', 'wb').write(b'x')
        res = pm.measure(bam, 'S1', modules=('te',), run=fake_run(self.cfg, {'chr7ins_c1': (0, 6), 'chr7ins_c2': (0, 6)}), cfgs={'te': self.cfg})
        pm.write_json(d + '/S1.panel.json', res)
        import cichlid_common.common_light as cl
        cfg_path = d + '/tepanel.yaml'; yaml.safe_dump(self.cfg, open(cfg_path, 'w')); orig = cl.config_path; cl.config_path = lambda n: cfg_path if n == 'tepanel.yaml' else orig(n)
        try:
            C.main(['--panel', d, '--out', d + '/te_calls.tsv'])
        finally:
            cl.config_path = orig
        t = pd.read_csv(d + '/te_calls.tsv', sep='\t', keep_default_na=False)
        self.assertEqual((t.SampleID.tolist(), t.chr7ins_GT.tolist(), t.chr7ins_call.tolist()), (['S1'], ['1/1'], ['1/1']))


if __name__ == '__main__':
    unittest.main()
