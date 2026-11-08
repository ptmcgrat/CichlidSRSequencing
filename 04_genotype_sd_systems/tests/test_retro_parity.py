"""Retrospective (streamed) measurement and the parity check with sexloci, on fake samtools: region merging, disk-only modules refused, fetch failures classified, a resumable run whose JSONs feed the
usual call steps, and the comparison tables."""
import os, sys, tempfile, threading, types, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import yaml
from cichlid_common.common_light import config_path
from sd_genotypes import parity as P, retro as R
from sd_genotypes.gsdf import steps as GS
from sd_genotypes.panel import measure as pm
import test_te_panel as T

GCFG = GS.cfg_of(); LG = GCFG['region']['end'] - GCFG['region']['start']
RATIO = {'S1': 0.530, 'S2': 0.798, 'S3': 0.530, 'S4': 0.798}                               # planted depth ratios (2 and 3 copies)
FLAKY = {'n': 0}


class FakeServer:
    def __enter__(self):
        return 'http://fake'

    def __exit__(self, *a):
        return False


def fake_run(log=None, lock=threading.Lock()):
    def run(cmd):
        if log is not None:
            with lock:
                log.append(list(cmd))
        if cmd[:2] == ['samtools', '--version']:
            return types.SimpleNamespace(stdout='samtools 1.x\n')
        if cmd[:3] == ['samtools', 'view', '-b']:
            url = cmd[cmd.index('-o') + 2]
            if '/S3/' in url:
                raise RuntimeError('FAILED: samtools view\n[E::hts_open_format] 404 Not Found')
            if '/S4/' in url and FLAKY['n'] < 2:
                FLAKY['n'] += 1; raise RuntimeError('FAILED: samtools view\nconnection reset by peer')
            open(cmd[cmd.index('-o') + 1], 'wb').write(b'raw'); return types.SimpleNamespace(stdout='')
        if cmd[:2] == ['samtools', 'sort']:
            open(cmd[2:][cmd[2:].index('-o') + 1], 'wb').write(b'bam'); return types.SimpleNamespace(stdout='')
        if cmd[:2] == ['samtools', 'index']:
            open(cmd[2] + '.bai', 'wb').write(b'i'); return types.SimpleNamespace(stdout='')
        if cmd[:3] == ['samtools', 'view', '-H']:
            return types.SimpleNamespace(stdout='@HD\tVN:1.6\n@SQ\tSN:NC_135173.1\tLN:68000000\n')
        if cmd[0:2] == ['samtools', 'depth']:
            sample = os.path.basename(cmd[-1]).split('.')[0]; ratio = RATIO[sample]
            r = GCFG['region']; u0, u1 = GCFG['unit_window']; j0, j1 = GCFG['junction_window']; q = '-Q' in cmd; out = []
            for i in range(LG):
                pos = r['start'] + i; d = 15.0 * ratio if u0 <= pos < u1 + 5000 else 15.0
                if q and j0 <= pos < j1:
                    d = 15.0 * (ratio - 0.530) / 0.268 / 2
                out.append(f"{r['chrom']}\t{pos + 1}\t{d:g}\n")
            return types.SimpleNamespace(stdout=''.join(out))
        raise AssertionError(cmd)
    return run


class TestRegions(unittest.TestCase):
    def test_overlapping_and_adjacent_regions_are_merged(self):
        self.assertEqual(R.merge_region_strings(['c:1-10', 'c:5-20', 'd:1-3', 'c:21-30', 'c:40-50']), ['c:1-30', 'c:40-50', 'd:1-3'])
        with self.assertRaises(ValueError):
            R.merge_region_strings(['nonsense'])

    def test_modules_with_regions_and_disk_only_modules(self):
        self.assertEqual(R.module_regions(['gsdf']), ['NC_135173.1:50700001-50830000'])
        cfg, _ = T.make_def()
        regs = R.module_regions(['gsdf', 'te'], cfgs={'te': cfg})
        self.assertEqual(regs, ['NC_135173.1:50700001-50830000', 'NC_135176.1:24920728-24921845'])      # both chr7ins windows lie INSIDE the 130 kb gsdf region and merge into it; TE1 and TE2 merge into one window
        for m in ('bchrom', 'lg5inv'):
            with self.assertRaises(R.DiskOnly):
                R.module_regions(['gsdf', m])
        with self.assertRaises(ValueError):
            R.module_regions(['nope'])


class TestFetch(unittest.TestCase):
    def test_commands_and_cleanup(self):
        log = []; d = tempfile.mkdtemp()
        res = R.measure_sample('S1', 'http://fake/S1/S1.all.bam', ['gsdf'], d, run=fake_run(log))
        v = [c for c in log if c[:3] == ['samtools', 'view', '-b']][0]
        self.assertEqual((v[v.index('-o') + 2], v[v.index('-o') + 3:]), ('http://fake/S1/S1.all.bam', ['NC_135173.1:50700001-50830000']))
        self.assertEqual([c[1] for c in log if c[1] in ('sort', 'index')], ['sort', 'index'])
        self.assertEqual((res['source'], res['retro_regions'], sorted(res['modules'])), ('http://fake/S1/S1.all.bam', ['NC_135173.1:50700001-50830000'], ['gsdf']))
        self.assertEqual(os.listdir(d), [])                                                  # the small BAM, its index and the raw file are removed
        R.measure_sample('S1', 'http://fake/S1/S1.all.bam', ['gsdf'], d, run=fake_run(), keep_bam=True); self.assertIn('S1.regions.bam', os.listdir(d))

    def test_failures_are_classified_and_only_transient_ones_are_retried(self):
        log = []; slept = []
        with self.assertRaises(R.FetchError) as e:
            R.measure_sample('S3', 'http://fake/S3/S3.all.bam', ['gsdf'], tempfile.mkdtemp(), run=fake_run(log), sleep=slept.append)
        self.assertEqual((e.exception.kind, sum(c[:3] == ['samtools', 'view', '-b'] for c in log), slept), ('no_bam', 1, []))
        FLAKY['n'] = 0
        res = R.measure_sample('S4', 'http://fake/S4/S4.all.bam', ['gsdf'], tempfile.mkdtemp(), run=fake_run(), sleep=slept.append)      # two transient failures, then it works
        self.assertEqual((slept, 'gsdf' in res['modules']), ([2, 4], True))
        FLAKY['n'] = -10
        with self.assertRaises(R.FetchError) as e:
            R.measure_sample('S4', 'http://fake/S4/S4.all.bam', ['gsdf'], tempfile.mkdtemp(), run=fake_run(), sleep=lambda s: None)
        self.assertEqual(e.exception.kind, 'error'); FLAKY['n'] = 0

    def test_a_disk_only_module_is_refused_before_any_network_step(self):
        log = []
        with self.assertRaises(R.DiskOnly):
            R.run_samples(['S1'], 'remote:', ['bchrom'], tempfile.mkdtemp(), run=fake_run(log), server_factory=lambda: FakeServer())
        self.assertEqual(log, [])


class TestRunSamples(unittest.TestCase):
    def test_a_resumable_run_whose_json_feeds_the_gsdf_call(self):
        out = tempfile.mkdtemp(); FLAKY['n'] = 0; msgs = []
        r = R.run_samples(['S1', 'S2', 'S3'], 'remote:', ['gsdf'], out, jobs=2, run=fake_run(), server_factory=lambda: FakeServer(), log=msgs.append, sleep=lambda s: None)
        self.assertEqual((sorted(r['done']), sorted(r['failed']), r['skipped']), (['S1', 'S2'], ['S3'], []))
        self.assertTrue(os.path.exists(out + '/S1.panel.json') and not os.path.exists(out + '/S3.panel.json'))
        self.assertIn('no_bam', open(out + '/retro_failures.tsv').read())
        r2 = R.run_samples(['S1', 'S2', 'S3'], 'remote:', ['gsdf'], out, run=fake_run(), server_factory=lambda: FakeServer(), log=msgs.append, sleep=lambda s: None)
        self.assertEqual((r2['skipped'], sorted(r2['failed'])), (['S1', 'S2'], ['S3']))        # finished samples are not fetched again
        r3 = R.run_samples(['S1'], 'remote:', ['gsdf'], out, run=fake_run(), server_factory=lambda: FakeServer(), refresh=True, log=msgs.append)
        self.assertEqual(r3['done'], ['S1'])
        calls_dir = tempfile.mkdtemp(); GS.call_step(out, out_dir=calls_dir)                  # the retro JSONs are ordinary panel JSONs
        c = pd.read_csv(calls_dir + '/gsdf_calls.tsv', sep='\t', keep_default_na=False, na_values=['']).set_index('SampleID')
        self.assertEqual((int(c.loc['S1', 'gsdf_number']), int(c.loc['S2', 'gsdf_number'])), (2, 3))


class TestParity(unittest.TestCase):
    def tables(self):
        ours = pd.DataFrame(dict(SampleID=['a', 'b', 'c', 'd'], TE1_GT=['0/0', '0/1', '1/1', './.'], TE1_GQ=[30, 40, 25, 0], TE1_status=['OK'] * 4, TE1_ref=[9, 5, 0, 2], TE1_alt=[0, 5, 8, 1]))
        theirs = ours.copy(); return ours, theirs

    def test_te_parity_counts_and_genotypes(self):
        ours, theirs = self.tables()
        s, mm = P.compare_te(ours, theirs)
        self.assertEqual((int(s.n_both_measured[0]), int(s.n_counts_identical[0]), int(s.n_both_confident[0]), int(s.n_genotype_agree[0]), len(mm)), (4, 4, 3, 3, 0))
        theirs.loc[1, 'TE1_GT'] = '1/1'; theirs.loc[2, 'TE1_alt'] = 6
        s, mm = P.compare_te(ours, theirs)
        self.assertEqual(sorted(zip(mm.kind, mm.SampleID)), [('counts', 'c'), ('genotype', 'b')]); self.assertEqual(int(s.n_genotype_agree[0]), 2); self.assertEqual(float(s.max_count_difference[0]), 2.0)

    def test_gsdf_parity(self):
        ours = pd.DataFrame(dict(SampleID=['a', 'b', 'c'], ratio=[0.53, 0.80, 1.05], flank_depth=[15, 14, 13], junction_ratio=[0, 0.5, 1.0], gsdf_number=[2, 3, 4], call_flags=['', '', 'low_flank_depth']))
        theirs = pd.DataFrame(dict(SampleID=['a', 'b', 'c'], DepthRatio=[0.531, 0.80, 1.05], FlankDepth=[15, 14, 13], JunctionJ=[0, 0.5, 1.0], CopyNumber=[2, 4, 4], Flag=['', '', 'low depth']))
        s, mm = P.compare_gsdf(ours, theirs)
        self.assertEqual((s['n_compared'], s['n_unflagged_in_both'], s['n_copy_number_agree'], s['n_flagged_in_one_only']), (3, 2, 1, 0))
        self.assertAlmostEqual(s['max_ratio_difference'], 0.001, places=6); self.assertEqual(mm.SampleID.tolist(), ['b'])

    def test_cli_writes_the_summary_and_difference_tables(self):
        d = tempfile.mkdtemp(); ours, theirs = self.tables(); theirs.loc[1, 'TE1_GT'] = '1/1'
        ours.to_csv(d + '/ours.tsv', sep='\t', index=False); theirs.to_csv(d + '/sx.tsv', sep='\t', index=False)
        P.main(['--sexloci', d + '/sx.tsv', '--te-calls', d + '/ours.tsv', '--out-dir', d + '/o'])
        self.assertTrue(os.path.exists(d + '/o/parity_summary.txt')); self.assertIn('genotype', open(d + '/o/parity_te_differences.tsv').read())


if __name__ == '__main__':
    unittest.main()
