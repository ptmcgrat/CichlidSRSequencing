"""lg5inv with fake tools (the fakes of the bchrom self-test): the pools' slices and unmapped reads are kept, every finished pool is recorded as its own run folder, one failed pool does not
stop the others, a re-run resumes, `clean` never deletes the only copy. Nothing here touches the network or a real tool."""
import json, os, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import selftest_bchrom as M
from sd_genotypes.lg5inv import steps as S, workflow as LW

ALIGN = 'align (SamToFastq | minimap2 | MergeBamAlignment)'


class Ex(M.FakeExecutor):
    def run(self, cmd, label, cwd=None):
        super().run(cmd, label, cwd)
        if label == 'unmapped slice':
            self.unmapped_cmd = cmd
            open(cmd[cmd.index('-o') + 1], 'wb').write(b'unm')


def make_cfg(**slices):
    cfg = S.cfg_of(); cfg['download'].update(source='ena', backoff_s=0)
    cfg.update(bin_size=M.BIN, scale_min_chrom_bp=500_000); cfg['resources']['min_free_gb'] = 0
    cfg['slices'].update(slices)
    for r in cfg['runs']:
        r.pop('bases_gb', None)                                                     # the configured sizes are the real ones
    return cfg


class Env:
    def __init__(self, cfg=None, **ex_kw):
        self.cfg = cfg or make_cfg()
        self.tmp = tempfile.mkdtemp() + '/'
        self.proj = M.make_proj(self.tmp)
        glen = sum(M.LENS.values()) * M.BIN
        self.ids = [r['id'] for r in self.cfg['runs']]
        self.reports = {i: M.filereport(i, nbytes=1000, bases=int(12 * glen), alias=i) for i in self.ids}
        self.depth = lambda rid: {c: np.full(n, 12.0, dtype=np.float32) for c, n in M.LENS.items()}
        self.ex = Ex(self.depth, M.BIN, **ex_kw)
        self.downloads = []

    def downloader(self, rec, d):
        os.makedirs(d, exist_ok=True); self.downloads.append(rec['name'])
        open(d + rec['name'], 'wb').write(b'fq'); return d + rec['name']

    def prepare(self, only=None, **kw):
        return S.prepare_step(self.proj, only=only, executor=self.ex, cmds=M.FakeCmds, fetch=lambda i: self.reports[i], downloader=self.downloader, cfg=self.cfg, **kw)

    @property
    def work(self):
        return S.work_dir(self.proj, self.cfg)

    def run_dir(self, rid):
        base = f'{self.proj.root}lg5inv/pool_{rid}/'
        return base + open(base + 'LATEST').read().strip() + '/'


class TestConfig(unittest.TestCase):
    def test_the_five_pools_and_the_slice_definition(self):
        cfg = S.cfg_of(); by = S.validate_config(cfg)
        self.assertEqual(sorted(by), sorted(['SRR18083717', 'SRR18083718', 'SRR18083720', 'SRR18083724', 'SRR18083730']))
        self.assertEqual({r: by[r]['expected_inverted_fraction'] for r in by}, {'SRR18083717': 0.5, 'SRR18083718': 0.0, 'SRR18083720': 0.0, 'SRR18083724': 0.5, 'SRR18083730': 0.0})
        self.assertEqual(cfg['slices']['regions'], [dict(name='LG5_inversion_with_flanks', lg='LG5', start=19400000, end=26980000)])
        self.assertTrue(cfg['slices']['unmapped'])
        self.assertEqual(cfg['resources']['runs_in_flight'], 1)
        self.assertTrue(all(by[r]['sex'] == 'F' for r in by))

    def test_the_inversion_region_contains_both_breakpoints_and_the_flanks(self):
        r = S.cfg_of()['slices']['regions'][0]
        for bp in (19_657_805, 26_726_680):
            self.assertTrue(r['start'] + 200_000 < bp < r['end'] - 200_000)

    def test_bad_configs_are_refused(self):
        cfg = S.cfg_of(); cfg['runs'].append(dict(cfg['runs'][0]))
        with self.assertRaises(ValueError):
            S.validate_config(cfg)
        cfg = S.cfg_of(); cfg['runs'][0]['role'] = 'carrier'
        with self.assertRaises(ValueError):
            S.validate_config(cfg)
        cfg = S.cfg_of(); cfg['slices']['enabled'] = False
        with self.assertRaises(ValueError):
            S.validate_config(cfg)


class TestPrepare(unittest.TestCase):
    def test_a_pool_keeps_the_region_slice_and_the_unmapped_reads_and_is_recorded_at_once(self):
        e = Env(); e.prepare(only=['SRR18083720'])
        sl = e.work + 'slices/'
        self.assertTrue(all(os.path.exists(sl + f'SRR18083720{x}') for x in ('.slice.bam', '.slice.bam.bai', '.unmapped.bam')))
        self.assertFalse(os.path.exists(e.work + 'SRR18083720'))                         # reads, BAMs and scratch are gone
        self.assertIn('chrB:19400001-26980000', e.ex.slice_cmd)                          # 1-based, inclusive, LG5 resolved to the sequence name
        c = e.ex.unmapped_cmd
        self.assertEqual((c[0], c[1], c[c.index('-f') + 1], c[-1].endswith('SRR18083720.all.bam')), ('samtools', 'view', '4', True))
        labels = [l for l, r in e.ex.calls]
        self.assertLess(labels.index('index slice'), labels.index('unmapped slice'))
        self.assertLess(labels.index('unmapped slice'), labels.index('mosdepth'))        # both slices are cut before the BAM's depth is taken and everything is deleted
        d = e.run_dir('SRR18083720')
        self.assertTrue(all(os.path.exists(d + f) for f in ('run.tsv', 'summary.json', 'depth/SRR18083720.npz', 'depth/SRR18083720.json', 'slices/SRR18083720.slice.bam', 'slices/SRR18083720.unmapped.bam')))
        self.assertEqual(json.load(open(d + 'slices/regions.json')), {'LG5_inversion_with_flanks': 'chrB:19400001-26980000'})
        row = open(d + 'run.tsv').read().splitlines()[1].split('\t')
        self.assertIn('0.0', row)                                                         # expected_inverted_fraction of the hard negative

    def test_one_failed_pool_does_not_stop_the_other_and_a_rerun_resumes_where_it_failed(self):
        e = Env(fail_once=(ALIGN, 'SRR18083717'))
        with self.assertRaises(SystemExit):
            e.prepare(only=['SRR18083717', 'SRR18083720'])
        self.assertTrue(os.path.exists(f'{e.proj.root}lg5inv/pool_SRR18083720/LATEST'))
        self.assertFalse(os.path.exists(f'{e.proj.root}lg5inv/pool_SRR18083717/LATEST'))
        n = len(e.downloads)
        e.prepare(only=['SRR18083717', 'SRR18083720'])
        self.assertEqual(len(e.downloads), n)                                             # nothing downloaded again: the fetch of the failed pool was marked, the other is complete
        self.assertEqual([r for l, r in e.ex.calls if l == 'FastqToSam' and r == 'SRR18083717'], ['SRR18083717'])
        self.assertTrue(os.path.exists(f'{e.proj.root}lg5inv/pool_SRR18083717/LATEST'))

    def test_a_finished_pool_is_not_processed_again(self):
        e = Env(); e.prepare(only=['SRR18083718'])
        n = len(e.ex.calls); e.prepare(only=['SRR18083718'])
        self.assertEqual(len(e.ex.calls), n)

    def test_unmapped_reads_can_be_switched_off(self):
        e = Env(cfg=make_cfg(unmapped=False)); e.prepare(only=['SRR18083730'])
        self.assertNotIn('unmapped slice', [l for l, r in e.ex.calls])
        self.assertFalse(os.path.exists(e.work + 'slices/SRR18083730.unmapped.bam'))

    def test_a_pool_whose_unmapped_reads_are_missing_is_not_finished(self):
        e = Env(cfg=make_cfg(unmapped=False)); e.prepare(only=['SRR18083730'])
        e.cfg['slices']['unmapped'] = True                                                # asked for after the run
        rp = LW.Lg5RunProcessor(e.cfg['runs'][2], {'files': []}, e.work, e.cfg, dict(reference='/r', index='/i', slice_regions=['x']))
        self.assertEqual((rp.id, rp.finished()), ('SRR18083730', False))

    def test_unknown_run_ids_are_refused(self):
        e = Env()
        with self.assertRaises(SystemExit):
            e.prepare(only=['SRR0'])

    def test_a_foreign_work_folder_is_refused(self):
        e = Env(); w = e.tmp + 'P/G/tmp/lg5inv/'; os.makedirs(w); open(w + 'precious.txt', 'w').write('x')
        with self.assertRaises(SystemExit):
            S.work_dir(e.proj, e.cfg)


class TestStatusAndClean(unittest.TestCase):
    def test_status_and_a_clean_that_never_deletes_the_only_copy(self):
        e = Env()
        orig = S._record_run
        S._record_run = lambda *a, **k: (_ for _ in ()).throw(OSError('Dropbox is down'))
        try:
            with self.assertRaises(SystemExit):
                e.prepare(only=['SRR18083720'])                                           # finished, but could not be recorded: exit code 1, slices stay
        finally:
            S._record_run = orig
        t = S.status_step(e.proj, e.cfg)
        row = t[t.run == 'SRR18083720'].iloc[0]
        self.assertEqual((row.state, row.on_dropbox), ('complete', '-'))
        with self.assertRaises(SystemExit):
            S.clean_step(e.proj, cfg=e.cfg)                                               # the slices are the only copy
        e.prepare(only=['SRR18083720'])                                                   # the re-run records without redoing the work
        n = len([1 for l, r in e.ex.calls if l == 'FastqToSam'])
        self.assertEqual(n, 1)
        t = S.status_step(e.proj, e.cfg)
        self.assertEqual(t[t.run == 'SRR18083720'].iloc[0].on_dropbox, 'yes')
        S.clean_step(e.proj, cfg=e.cfg)
        self.assertFalse(os.path.exists(e.tmp + 'P/G/tmp/lg5inv/'))


if __name__ == '__main__':
    unittest.main()
