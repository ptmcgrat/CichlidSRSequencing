"""variant_calling panel against module 04's REAL call layer (sd_genotypes.panel.calls) and panel JSONs made under the REAL definitions (the frozen keys), with module 04's own fake-BAM helpers.

Skipped when module 04's call layer or its test helpers are not in the checkout (the git HEAD of 2026-11-17 did not have `panel/calls.py`: it was only on utaka). The other tests of this module use a stand-in
for that layer; this is the one that says the stand-in describes the real thing."""
import io, json, os, shutil, sys, tempfile, unittest
import numpy as np
import pandas as pd
from cichlid_common.common_light import REPO_ROOT
from cichlid_common.fake_sheet import FakeWorksheet
from unittest import mock
from variant_calling import cohort_cli, panelpass as P, panelpass_cli as C

try:
    from sd_genotypes.panel import calls as REAL, measure as pm
    from sd_genotypes.te import measure as TM
    from sd_genotypes.lg5inv import measure as LM
    from sd_genotypes.bchrom import measure as BM
    sys.path.append(os.path.join(REPO_ROOT, '04_genotype_sd_systems', 'tests'))
    import test_panel_calls as TPC, test_te_panel as TT, test_lg5inv_panel as TL, test_bchrom_panel as TB
    AVAILABLE, WHY = True, ''
except Exception as e:                                         # noqa: any import problem means "not here"
    AVAILABLE, WHY = False, f'module 04 call layer or its test helpers are not importable: {type(e).__name__}: {e}'

# name: (gsdf depth ratio, transposon read pairs (ref, alt) per site, LG5 alt fraction, B fold)
SCEN = {'M1': (0.530, {'chr7ins_c1': (0, 0), 'chr7ins_c2': (0, 0), 'TE1': (20, 0), 'TE2': (17, 0)}, 0.0, 1.0),                        # 2 copies, TEs 0/0, chr7ins without pairs (no call)
        'M2': (0.798, {'chr7ins_c1': (8, 1), 'chr7ins_c2': (8, 0), 'TE1': (4, 4), 'TE2': (6, 0)}, 0.5, 1.0),                           # 3 copies, TE1 0/1, chr7ins 0/0 with the dilution flag, LG5 1
        'M3': (1.066, {'chr7ins_c1': (0, 0), 'chr7ins_c2': (0, 0), 'TE1': (2, 0), 'TE2': (2, 0)}, 0.0, 4.0),                            # 4 copies, every transposon site a no_call, a B carrier
        'M4': (0.664, {'chr7ins_c1': (0, 0), 'chr7ins_c2': (0, 0), 'TE1': (15, 0), 'TE2': (0, 14)}, 1.0, 1.0)}                          # gsdf half way between 2 and 3 (flagged), TE2 1/1, LG5 2


def make_blocks(ratio, plan, lg5_p, fold, seed):
    tcfg = pm.module_cfg('te'); tcfg['backend'] = 'numpy'
    te = TM.measure_bam('x.bam', tcfg, TT.fake_run(tcfg, plan))
    lcfg = pm.module_cfg('lg5inv'); markers, _ = LM.load_definition(lcfg); n = len(markers); rng = np.random.default_rng(seed)
    d = rng.poisson(10, n); alt = rng.binomial(d, lg5_p)
    lg = LM.measure_bam('x.bam', lcfg, TL.fake_run(lcfg, (d - alt).tolist(), alt.tolist(), []))
    bcfg = pm.module_cfg('bchrom'); bl, _, meta = BM.load_definition(bcfg)
    bc = BM.measure_bam('x.bam', bcfg, TB.fake_run(bcfg, TB.depth_fn(bl, 10.0, fold), lengths=meta['lengths']))
    return {'gsdf': TPC.gsdf_block(ratio), 'te': te, 'lg5inv': lg, 'bchrom': bc}


@unittest.skipUnless(AVAILABLE, WHY)
class TestRealCallLayer(unittest.TestCase):
    HEAD = ['SampleID', 'Species', 'gsdf number', 'TE1', 'TE2', 'Ecogroup']

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(); cls.jd = cls.tmp + '/json/'; os.makedirs(cls.jd)
        for i, (s, (ratio, plan, p, fold)) in enumerate(SCEN.items()):
            pm.write_json(f'{cls.jd}{s}.panel.json', dict(schema=2, sample=s, bam=f'{s}.all.bam', bam_size_bytes=1000 + i, panel_version='v35-reorg1', measured='2026-11-17T12:00:00', samtools='fake',
                                                          modules=make_blocks(ratio, plan, p, fold, i + 1)))
        ratio, plan, p, fold = SCEN['M2']                                                  # U1: a non-mbuna fish with M2's measurements (the LG5 and B calls must stay out of its row)
        pm.write_json(f'{cls.jd}U1.panel.json', dict(schema=2, sample='U1', bam='U1.all.bam', bam_size_bytes=2000, panel_version='v35-reorg1', measured='2026-11-17T12:00:00', samtools='fake', modules=make_blocks(ratio, plan, p, fold, 9)))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def sheet(self):
        # M2's TE1 and M3's gsdf already hold a DIFFERENT value, M1's gsdf the same one
        return FakeWorksheet(self.HEAD, [['M1', 'Cm', '2', '', '', 'Mbuna'], ['M2', 'Cm', '', '0/0', '', 'Mbuna'], ['M3', 'Lf', '3', '', '', 'Mbuna'], ['M4', 'Lf', '', '', '', 'Mbuna'], ['U1', 'Pm', '', '', '', 'Utaka'], ['OTHER', 'Lf', '5', '0/0', '', 'Utaka']])

    def write(self, ws, *args):
        out, err = io.StringIO(), io.StringIO()
        rc = C.main(['write', '--json-dir', self.jd] + list(args), out=out, err=err, worksheet_fn=lambda: ws)               # no calls_module: the real one is imported
        return rc, out.getvalue(), err.getvalue()

    def test_the_stand_in_describes_the_real_layer(self):
        r = REAL.call_panel(f'{self.jd}M2.panel.json')
        self.assertEqual([c['column'] for c in r['calls']], ['gsdf number', 'chr7ins', 'TE1', 'TE2', 'Inversion5', 'B score'])
        for c in r['calls']:
            for k in ('column', 'flags_column', 'value', 'flags', 'confident', 'value_to_write', 'call_tag', 'details'): self.assertIn(k, c)
        values, flags = REAL.sheet_row(r); self.assertEqual(values['gsdf number'], 3); self.assertEqual(flags['chr7ins flags'], 'pilot-v1;dilution_possible')

    ADDED = ['gsdf flags', 'chr7ins', 'chr7ins flags', 'TE1 flags', 'TE2 flags', 'Inversion5', 'Inversion5 flags', 'B score', 'B flags']

    def got(self, ws):
        return {str(r[0]): {k: str(v) for k, v in zip(ws.grid[0], r)} for r in ws.grid[1:]}

    def test_the_writer_writes_what_the_handoff_says_and_nothing_else(self):
        ws = self.sheet(); rc, out, err = self.write(ws, '--apply'); self.assertEqual(rc, 0, err + out)
        self.assertEqual(ws.grid[0], self.HEAD + self.ADDED); g = self.got(ws)
        want = {'M1': dict(**{'gsdf number': '2', 'TE1': '0/0', 'TE2': '0/0', 'gsdf flags': 'pilot-v1', 'chr7ins': './.', 'chr7ins flags': 'pilot-v1;no_call', 'TE1 flags': 'pilot-v1', 'TE2 flags': 'pilot-v1', 'Inversion5': '0', 'Inversion5 flags': 'pilot-v1;junction_uninformative', 'B score': '1.0', 'B flags': 'pilot-v1'}),
                'M2': {'gsdf number': '3', 'TE1': '0/0', 'TE2': '0/0', 'gsdf flags': 'pilot-v1', 'chr7ins': '0/0', 'chr7ins flags': 'pilot-v1;dilution_possible', 'TE1 flags': '', 'TE2 flags': 'pilot-v1', 'Inversion5': '1', 'Inversion5 flags': 'pilot-v1;junction_uninformative', 'B score': '1.0', 'B flags': 'pilot-v1'},
                'M3': {'gsdf number': '3', 'TE1': './.', 'TE2': './.', 'gsdf flags': '', 'chr7ins': './.', 'chr7ins flags': 'pilot-v1;no_call', 'TE1 flags': 'pilot-v1;no_call', 'TE2 flags': 'pilot-v1;no_call', 'Inversion5': '0', 'Inversion5 flags': 'pilot-v1;junction_uninformative', 'B score': '4.0', 'B flags': 'pilot-v1'}}
        for s, row in want.items():
            for col, v in row.items(): self.assertEqual(g[s][col], v, f'{s} {col}')
        m4 = g['M4']; self.assertEqual((m4['gsdf number'], m4['TE1'], m4['TE2'], m4['chr7ins'], m4['gsdf flags'], m4['B score']), ('', '0/0', '1/1', './.', 'pilot-v1;non_integer', '1.0'))
        if 'unvalidated_value' in m4['Inversion5 flags']: self.assertEqual(m4['Inversion5'], '', 'module 04\'s guard (decision 79): an Inversion5 of 2 is never written')          # the call layer of the checkout has the guard
        else: self.assertEqual(m4['Inversion5'], '2', 'the call layer of this checkout predates decision 79')
        self.assertEqual({k: g['OTHER'][k] for k in ('gsdf number', 'TE1')}, {'gsdf number': '5', 'TE1': '0/0'}); self.assertEqual({g['OTHER'][c] for c in self.ADDED}, {''}, 'a sample with no JSON is untouched')
        u1 = g['U1']; self.assertEqual((u1['gsdf number'], u1['gsdf flags'], u1['chr7ins'], u1['TE1']), ('3', 'pilot-v1', '0/0', '0/1')); self.assertEqual([u1[c] for c in ('Inversion5', 'Inversion5 flags', 'B score', 'B flags')], ['', '', '', ''], 'LG5 and B only for the mbuna')
        self.assertIn('Inversion5: left out for 1 samples whose Ecogroup is not Mbuna', out); self.assertIn('B score: left out for 1 samples', out)
        self.assertIn('conflict 2', out); self.assertIn("sheet '0/0', measured 0/1", out); self.assertIn("sheet '3', measured 4", out); self.assertIn('2 flags cells are left alone', out)
        self.assertIn('every written cell was read back and is in the right row', out)

    def test_module_04s_guard_is_respected_by_the_writer(self):
        orig = REAL.call_panel
        def guarded(path, *a, **k):                                   # decision 79: an Inversion5 of 2 keeps its value but is not confident, value_to_write None, flag unvalidated_value
            r = orig(path, *a, **k)
            for c in r['calls']:
                if c['column'] == 'Inversion5' and c['value'] == 2: c['value_to_write'] = None; c['confident'] = False; c['flags'] = c['flags'] + ['unvalidated_value']
            return r
        ws = self.sheet()
        with mock.patch.object(REAL, 'call_panel', guarded): rc, out, err = self.write(ws, '--apply')
        self.assertEqual(rc, 0, err + out); m4 = self.got(ws)['M4']
        self.assertEqual((m4['Inversion5'], m4['Inversion5 flags']), ('', 'pilot-v1;junction_uninformative;unvalidated_value'), 'no value, no ./., the reason in the flags cell')

    def test_a_tab_without_an_ecogroup_column_is_refused_with_the_real_layer_too(self):
        ws = FakeWorksheet(['SampleID', 'Species', 'gsdf number', 'TE1', 'TE2'], [['M1', 'Cm', '', '', '']]); rc, out, err = self.write(ws, '--apply')
        self.assertEqual(rc, 3); self.assertIn('cannot give it', err); self.assertEqual(list(ws.calls), [])

    def run_readback(self, ws):
        import csv, runpy, contextlib
        sheet = self.tmp + '/sheet_export.tsv'
        with open(sheet, 'w', newline='') as fh: csv.writer(fh, delimiter='\t').writerows(ws.grid)
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts', 'panel_readback.py'); cfg = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config', 'panelpass.yaml')
        buf = io.StringIO()
        with mock.patch.dict(os.environ, {'PANEL': self.jd, 'SHEET': sheet, 'CONFIG': cfg}), contextlib.redirect_stdout(buf): runpy.run_path(script)
        return buf.getvalue()

    def test_the_readback_script_agrees_with_the_writer_and_sees_a_difference(self):
        ws = FakeWorksheet(self.HEAD, [[s, 'x', '', '', '', 'Mbuna' if s.startswith('M') else 'Utaka'] for s in ('M1', 'M2', 'M3', 'M4', 'U1')]); self.write(ws, '--apply'); out = self.run_readback(ws)          # a clean sheet: no cell holds another value
        self.assertIn('0 differences', out); self.assertTrue(out.rstrip().endswith('PASS'), out)
        h = ws.grid[0]; row = [i for i, r in enumerate(ws.grid) if r[0] == 'U1'][0]; ws.grid[row][h.index('B score')] = '4.0'                 # a value where the scope says there must be none
        out = self.run_readback(ws); self.assertIn("('U1', 'B score', '4.0'", out); self.assertTrue(out.rstrip().endswith('FAIL'), out)
        row = [i for i, r in enumerate(ws.grid) if r[0] == 'M3'][0]; ws.grid[row][h.index('TE1')] = ''                                         # a ./. that went missing
        out = self.run_readback(ws); self.assertIn("('M3', 'TE1', ''", out)

    def test_the_readback_reports_the_conflicts_the_writer_left_alone(self):
        ws = self.sheet(); self.write(ws, '--apply'); out = self.run_readback(ws)                # the sheet has M2's TE1 and M3's gsdf with a DIFFERENT value than the panel's
        self.assertIn("('M2', 'TE1', '0/0'", out); self.assertIn("('M3', 'gsdf number', '3'", out); self.assertTrue(out.rstrip().endswith('FAIL'), out)

    def cfg_file(self, **sheet):
        import yaml
        base = P.load_config(); base['sheet'].update(sheet)
        p = tempfile.mkdtemp() + '/cfg.yaml'; yaml.safe_dump({'sheet': base['sheet']}, open(p, 'w')); return p

    def test_the_shipped_config_gives_dot_slash_dot_only_to_genotype_cells_that_could_not_be_called(self):
        ws = self.sheet(); self.write(ws, '--apply')
        got = {r[0]: dict(zip(ws.grid[0], r)) for r in ws.grid[1:]}
        self.assertEqual([got['M1']['chr7ins'], got['M3']['chr7ins'], got['M3']['TE1'], got['M3']['TE2'], got['M4']['chr7ins']], ['./.'] * 5, 'the real no_call cells')
        self.assertEqual(got['M2']['chr7ins'], '0/0', 'dilution_possible is information only: the value is written'); self.assertEqual(got['M4']['gsdf number'], '', 'a flagged number stays empty, never ./.')

    def test_without_the_nocall_block_those_cells_stay_empty(self):
        ws = self.sheet(); rc, out, err = self.write(ws, '--apply', '--config', self.cfg_file(nocall={})); self.assertEqual(rc, 0, err + out)
        got = {r[0]: dict(zip(ws.grid[0], r)) for r in ws.grid[1:]}
        self.assertEqual([got['M1']['chr7ins'], got['M3']['chr7ins'], got['M3']['TE1'], got['M4']['chr7ins']], [''] * 4); self.assertNotIn('./.', [v for r in ws.grid[1:] for v in r])

    def test_flags_in_sheet_false_with_the_real_layer_writes_only_values(self):
        ws = self.sheet(); rc, out, err = self.write(ws, '--apply', '--config', self.cfg_file(flags_in_sheet=False)); self.assertEqual(rc, 0, err + out)
        self.assertEqual(ws.grid[0], self.HEAD + ['chr7ins', 'Inversion5', 'B score']); got = {r[0]: dict(zip(ws.grid[0], r)) for r in ws.grid[1:]}
        self.assertEqual((got['M1']['chr7ins'], got['M2']['chr7ins'], got['M3']['TE1']), ('./.', '0/0', './.'))
        jd = self.jd + 'sheet_writes/'; self.assertTrue(any(f.endswith('.pilot-v1.tsv') for f in os.listdir(jd)))

    def test_a_second_apply_writes_nothing(self):
        ws = self.sheet(); self.write(ws, '--apply'); n = len(ws.calls); rc, out, _ = self.write(ws, '--apply'); self.assertEqual(rc, 0); self.assertIn('nothing to write', out); self.assertEqual(len(ws.calls), n)

    def test_a_measurement_of_another_definition_is_refused_by_the_real_layer(self):
        d = tempfile.mkdtemp()
        try:
            for f in os.listdir(self.jd):
                if f.endswith('.panel.json'): shutil.copy(self.jd + f, d)
            j = json.load(open(d + '/M1.panel.json')); j['modules']['gsdf']['definition'] = 'aaaaaaaaaaaa'; json.dump(j, open(d + '/M1.panel.json', 'w'))
            out, err = io.StringIO(), io.StringIO(); rc = C.main(['write', '--json-dir', d], out=out, err=err, worksheet_fn=lambda: self.sheet())
            self.assertEqual(rc, 0, err.getvalue()); self.assertIn('REFUSED M1: DefinitionMismatch', out.getvalue()); self.assertIn('5 samples with a JSON', out.getvalue()); self.assertIn('4 have something to write', out.getvalue())
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_the_table_has_all_six_columns_and_the_numbers_behind_them(self):
        out, err = io.StringIO(), io.StringIO(); path = self.tmp + '/calls.tsv'
        sample_dt = pd.DataFrame({'SampleID': list(SCEN), 'Species': ['Cm', 'Cm', 'Lf', 'Lf'], 'Ecogroup': 'Mbuna'})
        with mock.patch.object(cohort_cli, 'load_tables', return_value=(sample_dt, pd.DataFrame())):                      # never the real sheets in a test
            rc = C.main(['table', '--json-dir', self.jd, '--out', path], out=out, err=err)
        self.assertEqual(rc, 0, err.getvalue())
        t = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False).set_index('SampleID')
        for col in ('gsdf number', 'chr7ins', 'TE1', 'TE2', 'Inversion5', 'B score', 'gsdf flags', 'B flags', 'Inversion5 flags'): self.assertIn(col, t.columns)
        self.assertEqual((t.loc['M1', 'Species'], t.loc['M3', 'Species']), ('Cm', 'Lf'))
        self.assertEqual((t.loc['M1', 'gsdf number'], t.loc['M2', 'gsdf number'], t.loc['M3', 'gsdf number'], t.loc['M4', 'gsdf number']), ('2', '3', '4', '2'))
        self.assertEqual((t.loc['M4', 'gsdf number confident'], t.loc['M4', 'gsdf flags']), ('False', 'pilot-v1;non_integer')); self.assertEqual((t.loc['M1', 'Inversion5'], t.loc['M2', 'Inversion5'], t.loc['M4', 'Inversion5']), ('0', '1', '2'))
        self.assertEqual((t.loc['M1', 'B score'], t.loc['M3', 'B score']), ('1.0', '4.0')); self.assertEqual((t.loc['M2', 'TE1'], t.loc['M4', 'TE2'], t.loc['M3', 'TE1']), ('0/1', '1/1', ''))
        self.assertTrue(any(c.startswith('Inversion5.') for c in t.columns)); self.assertTrue(any(c.startswith('B score.') for c in t.columns)); self.assertIn('B score: 5 confident', out.getvalue())


if __name__ == '__main__':
    unittest.main()
