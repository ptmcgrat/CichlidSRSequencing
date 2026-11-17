"""The pilot call layer (sd_genotypes.panel.calls): uniform calls from STORED measurements for gsdf, the three transposon sites, the LG5 inversion and the B chromosome; flags and confidence; the call tag;
thresholds that change calls without any BAM; refusal of measurements of another definition; the writer helper; and that the new config keys never touch a definition key."""
import copy, json, os, subprocess, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import yaml
from cichlid_common.common_light import module_dirs, REPO_ROOT
from sd_genotypes.gsdf import measure as GM, steps as GS
from sd_genotypes.lg5inv import measure as LM
from sd_genotypes.bchrom import measure as BM
from sd_genotypes.te import call as TC, measure as TM
from sd_genotypes.panel import calls as PC, measure as pm
import test_te_panel as TT, test_lg5inv_panel as TL, test_bchrom_panel as TB

GCFG = GS.cfg_of()


def gsdf_block(ratio, jr=0.0, flank=15.0):
    n = (GCFG['region']['end'] - GCFG['region']['start']) // GCFG['bin_bp']; b0 = GCFG['region']['start']; bp = GCFG['bin_bp']
    bins = [flank] * n
    for i in range((GCFG['unit_window'][0] - b0) // bp, (GCFG['unit_window'][1] - b0) // bp + 1):         # the bins that START inside the unit window
        bins[i] = flank * ratio
    return dict(definition=GM.definition_key(GCFG), bins=bins, bins_q=list(bins), junction_depth=flank * ratio, junction_depth_q=flank * jr)


def te_block(plan, cfg=None):
    cfg = cfg or TT.make_def()[0]
    return cfg, TM.measure_bam('x.bam', cfg, TT.fake_run(cfg, plan))


def lg5_block(p, depth=10, rng_seed=1):
    cfg = dict(TL.CFG0)
    cfg['definition_dir'] = tempfile.mkdtemp(); cfg = TL.make_def(cfg['definition_dir'])
    b = dict(TL.raw(np.random.default_rng(rng_seed), depth=depth, p=p, jc=(4, 4, 4))); b['definition'] = LM.definition_key(cfg)
    return cfg, b


def bchrom_block(fold, base=10.0):
    cfg, _ = TB.make_def(); blocks, _, _ = BM.load_definition(cfg)
    return cfg, BM.measure_bam('x.bam', cfg, TB.fake_run(cfg, TB.depth_fn(blocks, base, fold)))


by_col = lambda calls: {c['column']: c for c in calls}


class TestGsdf(unittest.TestCase):
    def test_copies_flags_and_confidence(self):
        for ratio, want in [(0.530, 2), (0.798, 3), (1.066, 4), (1.33, 5)]:
            c = PC.call_module('gsdf', gsdf_block(ratio), GCFG)[0]
            self.assertEqual((c['column'], c['flags_column'], c['value'], c['confident'], c['value_to_write'], c['call_tag']), ('gsdf number', 'gsdf flags', want, True, want, 'pilot-v1'), ratio)
        c = PC.call_module('gsdf', gsdf_block(0.664), GCFG)[0]                                        # half way between 2 and 3: reported, flagged, not confident
        self.assertIn('non_integer', c['flags']); self.assertFalse(c['confident']); self.assertIsNone(c['value_to_write']); self.assertIn(c['value'], (2, 3))
        c = PC.call_module('gsdf', gsdf_block(0.53, flank=2.0), GCFG)[0]; self.assertIn('low_flank_depth', c['flags']); self.assertFalse(c['confident'])
        c = PC.call_module('gsdf', gsdf_block(1.60), GCFG)[0]; self.assertEqual((c['value'], 'outside_calibrated_range' in c['flags'], c['confident']), (6, True, False))

    def test_the_junction_check_is_information_only(self):
        c = PC.call_module('gsdf', gsdf_block(0.53, jr=0.5), GCFG)[0]
        self.assertEqual((c['value'], c['flags'], c['confident']), (2, ['junction_in_two_copy'], True))

    def test_the_value_is_the_one_of_the_gsdf_call_step(self):
        from sd_genotypes.gsdf import call as GC
        for ratio in (0.50, 0.55, 0.78, 0.80, 1.0, 1.07):
            b = gsdf_block(ratio); row = GC.call_table([dict(SampleID='x', **b)], GCFG).iloc[0]
            self.assertEqual(PC.call_module('gsdf', b, GCFG)[0]['value'], int(row.gsdf_number))


class TestTransposons(unittest.TestCase):
    def test_three_sites_with_their_columns(self):
        cfg, blk = te_block({'chr7ins_c1': (4, 0), 'chr7ins_c2': (4, 0), 'TE1': (0, 6)})
        c = by_col(PC.call_module('te', blk, cfg))
        self.assertEqual(sorted(c), ['TE1', 'TE2', 'chr7ins'])
        self.assertEqual((c['chr7ins']['value'], c['chr7ins']['confident'], c['chr7ins']['flags_column']), ('0/0', True, 'chr7ins flags'))
        self.assertEqual((c['TE1']['value'], c['TE1']['confident']), ('1/1', True))
        self.assertEqual((c['TE2']['value'], c['TE2']['flags'], c['TE2']['confident'], c['TE2']['value_to_write']), (None, ['no_call'], False, None))   # no reads at TE2: no call, no value

    def test_a_genotype_with_a_low_gq_is_reported_but_not_confident(self):
        cfg = TT.make_def()[0]; low = next((r, a) for r in range(0, 9) for a in range(0, 9) if r + a >= 5 and TC.call_from_counts(r, a)[0] in ('0/0', '0/1', '1/1') and TC.call_from_counts(r, a)[1] < cfg['call']['min_gq'])
        _, blk = te_block({'chr7ins_c1': (low[0], low[1]), 'chr7ins_c2': (0, 0)}, cfg)
        c = by_col(PC.call_module('te', blk, cfg))['chr7ins']
        self.assertEqual((c['flags'], c['confident'], c['value_to_write']), (['low_gq'], False, None)); self.assertIn(c['value'], ('0/0', '0/1', '1/1'))

    def test_dilution_needs_the_gsdf_copy_number(self):
        cfg, blk = te_block({'chr7ins_c1': (8, 1), 'chr7ins_c2': (8, 0)})
        self.assertEqual(by_col(PC.call_module('te', blk, cfg))['chr7ins']['flags'], [])
        c = by_col(PC.call_module('te', blk, cfg, context={'gsdf_copies': 3}))['chr7ins']
        self.assertEqual((c['value'], c['flags'], c['confident']), ('0/0', ['dilution_possible'], True))      # information only: the call stays confident
        self.assertEqual(by_col(PC.call_module('te', blk, cfg, context={'gsdf_copies': 2}))['chr7ins']['flags'], [])

    def test_a_failed_site_is_a_site_error(self):
        cfg = TT.make_def()[0]
        blk = TM.measure_bam('x.bam', cfg, TT.fake_run(cfg, {'chr7ins_c1': (0, 5), 'chr7ins_c2': (0, 5)}, fail_site='TE1'))
        c = by_col(PC.call_module('te', blk, cfg))
        self.assertEqual((c['TE1']['flags'], c['TE1']['confident']), (['site_error'], False)); self.assertEqual((c['chr7ins']['value'], c['chr7ins']['confident']), ('1/1', True))


class TestInversionAndB(unittest.TestCase):
    def test_lg5_zero_one_two(self):
        for p, want in [(0.004, 0), (0.5, 1), (0.985, 2)]:
            cfg, b = lg5_block(p); c = PC.call_module('lg5inv', b, cfg)[0]
            self.assertEqual((c['column'], c['flags_column'], c['value'], c['confident'], type(c['value'])), ('Inversion5', 'Inversion5 flags', want, True, int), p)
        cfg, b = lg5_block(0.5, depth=1.5); c = PC.call_module('lg5inv', b, cfg)[0]
        self.assertIn('low_depth', c['flags']); self.assertFalse(c['confident']); self.assertIsNone(c['value_to_write'])

    def test_lg5_flags_that_are_information_only(self):
        cfg, b = lg5_block(0.5); b = dict(b, junction=dict(A=0, B=0, C=10, reads=10))        # the markers say heterozygous, the junction says standard
        c = PC.call_module('lg5inv', b, cfg)[0]
        self.assertEqual((c['value'], 'junction_disagrees' in c['flags'], c['confident']), (1, True, True))

    def test_the_b_chromosome_is_a_score_and_nothing_else(self):
        cfg, blk = bchrom_block(4.0); c = PC.call_module('bchrom', blk, cfg)[0]
        self.assertEqual((c['column'], c['flags_column'], c['value'], c['confident'], type(c['value'])), ('B score', 'B flags', 4.0, True, float))
        self.assertNotIn(c['value'], ('B+', 'B-', 'B?')); self.assertFalse(any(f in ('B+', 'B-', 'B?') for f in c['flags'])); self.assertNotIn('class', ' '.join(c['details']))
        cfg, blk = bchrom_block(1.0, base=2.0); c = PC.call_module('bchrom', blk, cfg)[0]
        self.assertIn('low_depth', c['flags']); self.assertFalse(c['confident']); self.assertIsNone(c['value_to_write']); self.assertAlmostEqual(c['value'], 1.0, places=2)


class TestPanelAndWriterHelper(unittest.TestCase):
    def full(self):
        tcfg, tblk = te_block({'chr7ins_c1': (8, 1), 'chr7ins_c2': (8, 0), 'TE1': (4, 4), 'TE2': (6, 0)})
        lcfg, lblk = lg5_block(0.5); bcfg, bblk = bchrom_block(4.0)
        panel = dict(schema=2, sample='S1', bam='S1.all.bam', bam_size_bytes=999, panel_version='v', measured='t', samtools='s',
                     modules={'gsdf': gsdf_block(0.798), 'te': tblk, 'lg5inv': lblk, 'bchrom': bblk})
        return panel, {'gsdf': GCFG, 'te': tcfg, 'lg5inv': lcfg, 'bchrom': bcfg}

    def test_every_module_in_one_call_with_the_identity_and_the_gsdf_copies_feeding_the_dilution_flag(self):
        panel, cfgs = self.full(); r = PC.call_panel(panel, cfgs=cfgs)
        self.assertEqual((r['sample'], r['bam'], r['bam_size_bytes']), ('S1', 'S1.all.bam', 999))
        c = by_col(r['calls']); self.assertEqual(list(c), ['gsdf number', 'chr7ins', 'TE1', 'TE2', 'Inversion5', 'B score'])
        self.assertEqual((c['gsdf number']['value'], c['chr7ins']['value'], c['TE1']['value'], c['TE2']['value'], c['Inversion5']['value']), (3, '0/0', '0/1', '0/0', 1))
        self.assertEqual(c['chr7ins']['flags'], ['dilution_possible'])                         # 3 gsdf copies, a 0/0 with an insertion pair
        self.assertEqual(json.loads(json.dumps(r))['calls'][0]['value'], 3)                    # plain Python types: it serialises
        self.assertEqual({k: v['call_tag'] for k, v in c.items()}, {k: 'pilot-v1' for k in c})

    def test_the_writer_gets_values_only_for_confident_calls_and_a_tagged_flags_cell_always(self):
        panel, cfgs = self.full(); panel['modules']['gsdf'] = gsdf_block(0.664)                  # a flagged gsdf call
        r = PC.call_panel(panel, cfgs=cfgs); values, flags = PC.sheet_row(r)
        self.assertNotIn('gsdf number', values); self.assertEqual(values['TE1'], '0/1'); self.assertEqual(values['B score'], 4.0)
        self.assertTrue(flags['gsdf flags'].startswith('pilot-v1;') and 'non_integer' in flags['gsdf flags']); self.assertEqual(flags['TE2 flags'], 'pilot-v1')
        values, _ = PC.sheet_row(r, write_unconfident=True); self.assertIn('gsdf number', values)
        self.assertNotIn('chr7ins', PC.call_panel(panel, modules=['gsdf'], cfgs=cfgs)['calls'][0]['column'])

    def test_it_runs_on_a_json_file_and_from_the_command_line_without_the_bam(self):
        d = tempfile.mkdtemp(); panel = dict(schema=2, sample='S9', bam='S9.all.bam', bam_size_bytes=5, panel_version='v', measured='t', samtools='s', modules={'gsdf': gsdf_block(0.53)})
        pm.write_json(d + '/S9.panel.json', panel)                                             # no BAM anywhere
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(module_dirs()))
        r = subprocess.run([sys.executable, '-m', 'sd_genotypes.panel.calls', '--json', d + '/S9.panel.json'], capture_output=True, text=True, cwd=REPO_ROOT, env=env)
        self.assertEqual(r.returncode, 0, r.stderr[-300:]); self.assertIn('gsdf number', r.stdout); self.assertIn('pilot-v1', r.stdout); self.assertIn('S9.all.bam', r.stdout)


class TestSettings(unittest.TestCase):
    def test_a_threshold_changes_the_call_from_the_same_stored_measurement(self):
        cfg, b = lg5_block(0.5); self.assertEqual(PC.call_module('lg5inv', b, cfg)[0]['value'], 1)
        c2 = copy.deepcopy(cfg); c2['thresholds']['class_cuts'] = [0.6, 0.9]
        self.assertEqual(PC.call_module('lg5inv', b, c2)[0]['value'], 0)                       # the same measurement, other cuts, no BAM
        g = gsdf_block(0.53); g2 = copy.deepcopy(GCFG); g2['calibration'] = dict(g2['calibration'], ratio_at_two_copies=0.80)
        self.assertNotEqual(PC.call_module('gsdf', g, GCFG)[0]['value'], PC.call_module('gsdf', g, g2)[0]['value'])

    def test_the_call_tag_follows_status_and_version_and_a_config_without_them_is_refused(self):
        c2 = copy.deepcopy(GCFG); c2['version'] = 7; c2['status'] = 'provisional'
        self.assertEqual(PC.call_module('gsdf', gsdf_block(0.53), c2)[0]['call_tag'], 'provisional-v7')
        for k in ('status', 'version', 'sheet'):
            bad = {x: v for x, v in GCFG.items() if x != k}
            with self.assertRaises(PC.ConfigError):
                PC.call_module('gsdf', gsdf_block(0.53), bad)

    def test_a_measurement_of_another_definition_is_refused(self):
        b = gsdf_block(0.53); b['definition'] = 'other'
        with self.assertRaises(PC.DefinitionMismatch):
            PC.call_module('gsdf', b, GCFG)
        self.assertEqual(PC.call_module('gsdf', b, GCFG, check_definition=False)[0]['value'], 2)
        cfg, tb = te_block({'chr7ins_c1': (4, 0), 'chr7ins_c2': (4, 0)}); tb['definition'] = 'other'
        with self.assertRaises(PC.DefinitionMismatch):
            PC.call_module('te', tb, cfg)

    def test_the_new_config_keys_are_not_part_of_any_definition_key(self):
        for mod, cfg in (('gsdf', GCFG), ('te', TT.make_def()[0]), ('lg5inv', lg5_block(0.5)[0]), ('bchrom', bchrom_block(1.0)[0])):
            m = __import__(pm.MODULES[mod][0], fromlist=['x']); k = m.definition_key(cfg)
            bare = {x: v for x, v in cfg.items() if x not in ('status', 'version', 'sheet')}
            self.assertEqual(m.definition_key(bare), k, mod)                                   # removing status, version and sheet changes no key
            c2 = copy.deepcopy(cfg); c2['version'] = 99; c2['status'] = 'x'; self.assertEqual(m.definition_key(c2), k, mod)

    def test_the_real_configs_have_status_pilot_version_and_sheet_columns(self):
        for mod in pm.MODULES:
            c = pm.module_cfg(mod); self.assertEqual((c['status'], c['version']), ('pilot', 1), mod)
        self.assertEqual(pm.module_cfg('gsdf')['sheet']['column'], 'gsdf number'); self.assertEqual(sorted(pm.module_cfg('te')['sheet']['sites']), ['TE1', 'TE2', 'chr7ins'])
        self.assertEqual((pm.module_cfg('lg5inv')['sheet']['column'], pm.module_cfg('bchrom')['sheet']['column']), ('Inversion5', 'B score'))


if __name__ == '__main__':
    unittest.main()
