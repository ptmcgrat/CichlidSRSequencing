"""Tandem (double) inversions: the parts, and which assembly carries which (LG11 = 11a + 11b, LG20 = 20a + 20b). Synthetic alignments built like the real LG11 case (selftest_breakpoints):
the double carriers have two inverted segments with a gap between them; the single-part carrier (Rhamphochromis on LG20, Kumar et al. 2025) has one inverted run whose right end is the internal breakpoint."""
import glob, json, os, sys, tempfile, types, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from selftest_breakpoints import line, M
from cichlid_common import paf
from landmarks.inversion_breakpoints import parts as P, steps

T, LT, QL = 'NC_T20', 40 * M, 40.2 * M
PARTS = ['20a', '20b']


def flip(l):
    f = l.split('\t'); qlen, qs, qe = int(f[1]), int(f[2]), int(f[3]); f[2], f[3] = str(qlen - qe), str(qlen - qs); f[4] = '-' if f[4] == '+' else '+'; return '\t'.join(f)


def double(q='d'):
    """A B' C' D: 20a inverted in place (ref 9.0-16.6 Mb), a 0.8 Mb gap, 20b inverted in place (ref 17.4-21.3 Mb)."""
    return [line(q, QL, 0, 8.9 * M, '+', T, LT, 0, 8.9 * M),
            line(q, QL, 9.0 * M, 12.8 * M, '-', T, LT, 12.8 * M, 16.6 * M), line(q, QL, 12.9 * M, 16.7 * M, '-', T, LT, 9.0 * M, 12.8 * M),
            line(q, QL, 17.4 * M, 19.35 * M, '-', T, LT, 19.35 * M, 21.3 * M), line(q, QL, 19.45 * M, 21.4 * M, '-', T, LT, 17.4 * M, 19.35 * M),
            line(q, QL, 21.5 * M, 40.2 * M, '+', T, LT, 21.4 * M, 40.0 * M)]


def first_only(q='r', end=17.2):
    """20a inverted (ref 9.0 to `end` Mb), 20b and the rest forward: the right breakpoint of the single inversion is the internal breakpoint of the double."""
    h = (end - 9.0) / 2
    return [line(q, QL, 0, 8.9 * M, '+', T, LT, 0, 8.9 * M),
            line(q, QL, 9.0 * M, (9.0 + h) * M, '-', T, LT, (9.0 + h) * M, end * M), line(q, QL, (9.1 + h) * M, (9.1 + 2 * h) * M, '-', T, LT, 9.0 * M, (9.0 + h) * M),
            line(q, QL, (end + 0.1 + 0.0) * M, 40.0 * M, '+', T, LT, (end + 0.1) * M, 40.0 * M)]


def second_only(q='s', start=17.2):
    """20b inverted (ref `start` to 21.3 Mb), everything else forward."""
    h = (21.3 - start) / 2
    return [line(q, QL, 0, (start - 0.1) * M, '+', T, LT, 0, (start - 0.1) * M),
            line(q, QL, start * M, (start + h) * M, '-', T, LT, (start + h) * M, 21.3 * M), line(q, QL, (start + h + 0.1) * M, (21.3 + 0.1) * M, '-', T, LT, start * M, (start + h) * M),
            line(q, QL, 21.5 * M, 40.2 * M, '+', T, LT, 21.4 * M, 40.0 * M)]


def collinear(q='n'):
    return [line(q, 40 * M, 0, 40 * M, '+', T, LT, 0, 40 * M)]


def calls_of(**assemblies):
    rows = []
    for name, ls in assemblies.items():
        r = paf.call_method_AC(paf.parse_paf(ls), T)
        if len(r):
            rows.append(r.assign(assembly=name))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def resolve(calls, expect=None, assemblies=None):
    cons = paf.consensus(calls, 150_000, 100_000) if len(calls) else pd.DataFrame()
    return cons, *P.resolve_parts(calls, cons, PARTS, expect, 150_000, 100_000, assemblies)


class TestGeometry(unittest.TestCase):
    def test_the_synthetic_carriers_look_like_the_real_ones(self):
        d = paf.call_method_AC(paf.parse_paf(double()), T)
        self.assertEqual((len(d), d.method.iloc[0], int(d.n_segments.iloc[0])), (1, 'A', 2)); self.assertEqual(d.segments.iloc[0], '9000000-16600000;17400000-21300000')
        self.assertAlmostEqual(int(d.internal_breakpoints.iloc[0]), 17.0 * M, delta=1000)
        r = paf.call_method_AC(paf.parse_paf(first_only()), T)
        self.assertEqual((len(r), r.method.iloc[0], int(r.n_segments.iloc[0])), (1, 'A', 1)); self.assertAlmostEqual(r.right_lo.iloc[0], 17.2 * M, delta=1000)
        self.assertEqual(len(paf.call_method_AC(paf.parse_paf(collinear()), T)), 0)

    def test_the_consensus_alone_lists_the_single_part_carrier_as_a_different_inversion(self):
        calls = calls_of(fDipLim2=double(), fAulStu2=[flip(x) for x in double()], fRhaChi2=first_only())
        cons = paf.consensus(calls, 150_000, 100_000)
        self.assertEqual(len(cons), 2, 'this is why the parts need resolving: the consensus has no link between the two rows')
        self.assertEqual(sorted(cons.n_calls), [1, 2])


class TestResolver(unittest.TestCase):
    def setUp(self):
        self.calls = calls_of(fDipLim2=double(), fAulStu2=[flip(x) for x in double()], fRhaChi2=first_only(), YH=collinear())
        self.cons, self.parts, self.states, self.notes = resolve(self.calls, assemblies=['YH', 'fAulStu2', 'fDipLim2', 'fRhaChi2'])

    def test_both_parts_are_placed_with_their_carriers(self):
        a, b = self.parts.iloc[0], self.parts.iloc[1]
        self.assertEqual((a.part, b.part), ('20a', '20b'))
        self.assertAlmostEqual(a.left_breakpoint, 8.95 * M, delta=100_000); self.assertAlmostEqual(b.right_breakpoint, 21.35 * M, delta=100_000)
        self.assertEqual(a.right_breakpoint, b.left_breakpoint, 'the two parts share the internal breakpoint')
        self.assertEqual((a.carriers, a.n_carriers), ('fAulStu2;fDipLim2;fRhaChi2', 3)); self.assertEqual((b.carriers, b.n_carriers), ('fAulStu2;fDipLim2', 2))
        self.assertAlmostEqual(b.length, 4.1 * M, delta=400_000); self.assertGreater(a.length, 2 * b.length)

    def test_the_single_part_carrier_sharpens_the_internal_breakpoint(self):
        a = self.parts.iloc[0]
        self.assertGreaterEqual(a.right_min, 17.2 * M - 1000); self.assertLessEqual(a.right_max, 17.3 * M + 1000)
        self.assertAlmostEqual(a.right_breakpoint, 17.25 * M, delta=60_000, msg='the segment gap alone would say 17.0 Mb (its middle), 800 kb wide')
        self.assertIn('fRhaChi2', a.evidence); self.assertIn('inside the segment gap', a.evidence)
        self.assertEqual(self.notes, [])

    def test_without_a_single_part_carrier_the_gap_is_the_estimate_and_the_note_says_so(self):
        calls = calls_of(fDipLim2=double(), fAulStu2=[flip(x) for x in double()])
        cons, parts, states, notes = resolve(calls)
        a = parts.iloc[0]
        self.assertAlmostEqual(a.right_breakpoint, 17.0 * M, delta=1000); self.assertEqual((a.right_min, a.right_max), (16_600_000, 17_400_000))
        self.assertTrue(any('no assembly carries a single part' in n and '800000 bp wide' in n for n in notes), str(notes))

    def test_each_assembly_gets_its_state(self):
        s = dict(zip(self.states.assembly, self.states.state))
        self.assertEqual(s, {'YH': 'none', 'fAulStu2': 'both', 'fDipLim2': 'both', 'fRhaChi2': '20a only'})
        got = dict(zip(self.states.assembly, self.states.parts_carried)); self.assertEqual(got, {'YH': '', 'fAulStu2': '20a+20b', 'fDipLim2': '20a+20b', 'fRhaChi2': '20a'})

    def test_a_contig_that_is_reverse_complemented_gives_the_same_answer(self):
        c2 = calls_of(fDipLim2=[flip(x) for x in double()], fRhaChi2=[flip(x) for x in first_only()])
        _, parts, states, _ = resolve(c2)
        self.assertEqual(dict(zip(states.assembly, states.state)), {'fDipLim2': 'both', 'fRhaChi2': '20a only'}); self.assertEqual(parts.iloc[0].n_carriers, 2)

    def test_the_second_part_alone_is_recognised_too(self):
        calls = calls_of(fDipLim2=double(), fAulStu2=[flip(x) for x in double()], other=second_only())
        _, parts, states, notes = resolve(calls)
        self.assertEqual(dict(zip(states.assembly, states.state))['other'], '20b only')
        self.assertEqual(parts.iloc[1].carriers, 'fAulStu2;fDipLim2;other'); self.assertEqual(parts.iloc[0].carriers, 'fAulStu2;fDipLim2')
        self.assertIn('other', parts.iloc[1].evidence)

    def test_an_inversion_with_other_breakpoints_is_not_forced_into_a_part(self):
        odd = [line('o', QL, 0, 12.0 * M, '+', T, LT, 0, 12.0 * M), line('o', QL, 12.1 * M, 16.0 * M, '-', T, LT, 14.1 * M, 18.0 * M), line('o', QL, 16.1 * M, 18.0 * M, '-', T, LT, 12.0 * M, 14.1 * M),
               line('o', QL, 18.1 * M, 40.2 * M, '+', T, LT, 18.1 * M, 40.0 * M)]
        calls = calls_of(fDipLim2=double(), fAulStu2=[flip(x) for x in double()], odd=odd)
        _, parts, states, _ = resolve(calls)
        r = states[states.assembly == 'odd'].iloc[0]
        self.assertEqual((r.state, r.parts_carried), ('another inversion', '?')); self.assertIn('not the breakpoints of 20a or 20b', r.detail)

    def test_without_any_double_carrier_the_parts_cannot_be_placed(self):
        calls = calls_of(fRhaChi2=first_only(), YH=collinear())
        cons, parts, states, notes = resolve(calls, assemblies=['YH', 'fRhaChi2'])
        self.assertEqual(len(parts), 0); self.assertTrue(any('cannot be placed' in n for n in notes))
        self.assertEqual(dict(zip(states.assembly, states.state)), {'YH': 'none', 'fRhaChi2': 'another inversion'})

    def test_both_parts_with_a_single_segment_is_reported_as_unresolved_not_as_both(self):
        merged = [line('m', QL, 0, 8.9 * M, '+', T, LT, 0, 8.9 * M), line('m', QL, 9.0 * M, 15.2 * M, '-', T, LT, 15.1 * M, 21.3 * M), line('m', QL, 15.3 * M, 21.4 * M, '-', T, LT, 9.0 * M, 15.1 * M),
                  line('m', QL, 21.5 * M, 40.2 * M, '+', T, LT, 21.4 * M, 40.0 * M)]
        calls = calls_of(fDipLim2=double(), fAulStu2=[flip(x) for x in double()], m=merged)
        _, _, states, _ = resolve(calls)
        self.assertEqual(dict(zip(states.assembly, states.state))['m'], 'both (parts not resolved)')

    def test_three_parts_are_refused(self):
        with self.assertRaises(ValueError):
            P.resolve_parts(self.calls, self.cons, ['a', 'b', 'c'])


class TestExpectations(unittest.TestCase):
    def setUp(self):
        self.calls = calls_of(fDipLim2=double(), fAulStu2=[flip(x) for x in double()], fRhaChi2=first_only(), YH=collinear())

    def test_what_the_literature_says_is_compared_with_what_the_assemblies_show(self):
        _, parts, states, notes = resolve(self.calls, expect={'fDipLim2': ['20a', '20b'], 'fRhaChi2': ['20a']}, assemblies=['YH', 'fAulStu2', 'fDipLim2', 'fRhaChi2'])
        by = states.set_index('assembly')
        self.assertTrue(by.loc['fDipLim2', 'matches_expectation']); self.assertTrue(by.loc['fRhaChi2', 'matches_expectation'])
        self.assertEqual(by.loc['YH', 'matches_expectation'], '', 'no expectation, no verdict: the benthic assemblies are open'); self.assertEqual(by.loc['fRhaChi2', 'expected'], '20a')
        self.assertNotIn('DIFFERS', P.format_parts(parts, states, notes))

    def test_a_difference_is_flagged_loudly(self):
        _, parts, states, notes = resolve(self.calls, expect={'fRhaChi2': ['20a', '20b'], 'YH': ['20a', '20b']}, assemblies=['YH', 'fAulStu2', 'fDipLim2', 'fRhaChi2'])
        by = states.set_index('assembly')
        self.assertFalse(by.loc['fRhaChi2', 'matches_expectation']); self.assertFalse(by.loc['YH', 'matches_expectation'])
        text = P.format_parts(parts, states, notes)
        self.assertIn('DIFFERS FROM THE LITERATURE', text); self.assertIn("YH is 'none'", text); self.assertIn("fRhaChi2 is '20a only'", text)

    def test_the_shipped_configuration_declares_both_double_inversions(self):
        cfg = steps.cfg_of()['inversions']
        self.assertEqual(cfg['LG20']['parts'], ['20a', '20b']); self.assertEqual(cfg['LG11']['parts'], ['11a', '11b']); self.assertIn('fRhaChi2', cfg['LG20']['assemblies'])
        self.assertEqual(cfg['LG20']['expect'], {'fDipLim2': ['20a', '20b'], 'fRhaChi2': ['20a']}); self.assertTrue(cfg['LG20']['reference_is_standard'])
        self.assertNotIn('YH', cfg['LG20']['expect']); self.assertNotIn('fAulStu2', cfg['LG20']['expect'])
        for lg in ('LG2', 'LG9', 'LG10', 'LG13'):
            self.assertNotIn('parts', cfg[lg])

    def test_the_new_assembly_is_registered_with_its_accession(self):
        import yaml
        from cichlid_common.common_light import config_path
        g = yaml.safe_load(open(config_path('genomes.yaml')))['genomes']['fRhaChi2']
        self.assertEqual((g['source'], g['accession'], g['role']), ('ncbi', 'GCA_963969265.1', 'query'))


class TestBedAndText(unittest.TestCase):
    def test_bed_lines_name_every_breakpoint_of_every_part(self):
        calls = calls_of(fDipLim2=double(), fAulStu2=[flip(x) for x in double()], fRhaChi2=first_only())
        _, parts, states, notes = resolve(calls)
        bed = P.parts_bed(parts, T, 'LG20')
        self.assertEqual([b.split('\t')[3] for b in bed], ['LG20_20a_left', 'LG20_20a_right', 'LG20_20b_left', 'LG20_20b_right'])
        self.assertTrue(all(len(b.split('\t')) == 4 and int(b.split('\t')[1]) <= int(b.split('\t')[2]) for b in bed))
        self.assertEqual(bed[1].split('\t')[1:3], bed[2].split('\t')[1:3], 'the internal breakpoint interval is the same for both parts')


class TestThroughTheRealCallStep(unittest.TestCase):
    """The same data through steps.call_step with a fake project, as selftest_breakpoints does: the files appear next to the consensus."""
    def run_step(self, with_parts=True):
        from cichlid_common.runner import Run
        ref = 'Mzebra_GT3_NCBI'
        adir = tempfile.mkdtemp() + '/'
        data = dict(fDipLim2=double(), fAulStu2=[flip(x) for x in double()], fRhaChi2=first_only(), YH=collinear())
        for name, ls in data.items():
            open(f'{adir}{name}_vs_{ref}.targets.paf', 'w').write('\n'.join(ls) + '\n')
        spec = {'assemblies': list(data), 'reference_is_standard': True, 'label': 'LG20'}
        if with_parts:
            spec.update(parts=PARTS, expect={'fDipLim2': ['20a', '20b'], 'fRhaChi2': ['20a', '20b']})
        cfg = steps.cfg_of(); cfg = dict(cfg, inversions={'LG20': spec})
        real_cfg, real_excl = steps.cfg_of, steps.excluded_regions
        steps.cfg_of, steps.excluded_regions = (lambda: cfg), (lambda p_: {})
        loc = tempfile.mkdtemp() + '/'
        proj = types.SimpleNamespace(fm=types.SimpleNamespace(localMasterDir=loc, cloudMasterDir=loc + 'cloud/'), root=loc + 'P/G/', genome=ref, cfg={'runner': {'upload': False}},
                                     chrom_for_lg=lambda lg: T, chrom_lengths=lambda: {T: LT}, lg={T: 'LG20'}, cores=2, run=None)
        proj.run = lambda a, s, params=None, tar=(): Run(proj, a, s, params, tar)
        os.makedirs(proj.root)
        import contextlib, io
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                steps.call_step(proj, align_dir=adir)
        finally:
            steps.cfg_of, steps.excluded_regions = real_cfg, real_excl
        d = sorted(glob.glob(proj.root + 'inversion_breakpoints/call/*/'))[-1]
        return d, out.getvalue()

    def test_the_step_writes_the_parts_the_states_and_the_bed_and_prints_them(self):
        d, text = self.run_step()
        for f in ('LG20_parts.tsv', 'LG20_assembly_states.tsv', 'LG20_parts_breakpoint_regions.bed', 'LG20_consensus_breakpoints.tsv', 'LG20_strand_runs.tsv', 'LG20_polarity.json', 'summary.json'):
            self.assertTrue(os.path.exists(d + f), f)
        parts = pd.read_csv(d + 'LG20_parts.tsv', sep='\t'); states = pd.read_csv(d + 'LG20_assembly_states.tsv', sep='\t', keep_default_na=False)
        self.assertEqual(list(parts.part), ['20a', '20b']); self.assertEqual(dict(zip(states.assembly, states.state))['fRhaChi2'], '20a only')
        self.assertEqual(len(open(d + 'LG20_parts_breakpoint_regions.bed').read().split()), 16, '4 BED lines of 4 columns')
        info = json.load(open(d + 'LG20_polarity.json'))
        self.assertEqual(info['assembly_states']['fDipLim2'], 'both'); self.assertEqual(info['differs_from_literature'], ['fRhaChi2']); self.assertEqual([p['part'] for p in info['parts']], ['20a', '20b'])
        self.assertIn('-- the parts of LG20 --', text); self.assertIn('DIFFERS FROM THE LITERATURE', text)
        self.assertEqual(json.load(open(d + 'summary.json'))['LG20']['assembly_states']['YH'], 'none')

    def test_a_chromosome_without_parts_is_unchanged(self):
        d, text = self.run_step(with_parts=False)
        self.assertFalse(os.path.exists(d + 'LG20_parts.tsv')); self.assertTrue(os.path.exists(d + 'LG20_consensus_breakpoints.tsv')); self.assertNotIn('the parts of', text)
        self.assertNotIn('parts', json.load(open(d + 'LG20_polarity.json')))


if __name__ == '__main__':
    unittest.main()
