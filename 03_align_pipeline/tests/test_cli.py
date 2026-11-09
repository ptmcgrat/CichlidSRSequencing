"""The command line: what a run is made of (the ledger's per-sample facts, the slices, the panel, the commands) and the gates that stop a run BEFORE any work starts. A fake file manager stands in for the
lab's; no tool is run."""
import io, os, sys, tempfile, unittest
from unittest import mock
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from align_pipeline import cli, panel as pn, seq_pipeline as sp
from align_pipeline.alignment_db import LEDGER_COLUMNS
from test_alignment_db import FakeSheet, OLD

NAN = float('nan')
FAI = {'NC_135168.1': 45_000_000, 'NC_135176.1': 34_949_860, 'NC_135179.1': 30_000_000, 'NC_135171.1': 40_000_000}


SHIPPED = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'slice_regions', 'Mzebra_GT3_NCBI.tsv')


class FakeFM:
    def __init__(self, tmp, header=None, genome='G'):
        self.tmp, self.genome = tmp, genome
        self.localMasterDir, self.localErrorsDir, self.localGenomeDir = tmp + '/master/', tmp + '/errors/', tmp + '/genomes/'
        os.makedirs(self.localGenomeDir)
        self.localGenomeFile, self.localMinimapGenomeFile = self.localGenomeDir + 'G.fna', self.localGenomeDir + 'G.mmi'
        open(self.localGenomeFile, 'w').write('>x\n'); open(self.localMinimapGenomeFile, 'w').write('i')
        open(self.localGenomeFile + '.fai', 'w').write(''.join(f'{k}\t{v}\t0\t60\t61\n' for k, v in FAI.items()))
        open(self.localGenomeDir + 'G.dict', 'w').write('@HD')
        self.sample_dt = pd.DataFrame(dict(SampleID=['S1', 'S2', 'S3'], Species=['Maylandia zebra', 'Aulonocara minutus', NAN], Ecogroup=['Mbuna', 'Deep_Benthic', 'Exclude']))
        self.reads_dt = pd.DataFrame(dict(SampleID=['S1', 'S2', 'S2', 'S3'], RunID=['R1', 'R2b', 'R2a', 'R3'], ReadLength=[250, 250, 250, 150], Instrument=['HiSeq 2500', 'HiSeq 2500', 'HiSeq 2500', NAN],
                                          FileSize=[1e9, 3e9, 2e9, 4e9], ProjectID=['P'] * 4))
        self.merged_dt = self.sample_dt.merge(self.reads_dt, on='SampleID')
        self.genome_version, self.a_ID, self.downloads, self.sheet = None, 'A', [], FakeSheet(header or (OLD + [c for c in LEDGER_COLUMNS if c not in OLD]))
        outer = self
        self.gc = type('GC', (), {'open_by_key': lambda s, k: type('SS', (), {'worksheet': lambda s2, n: outer.sheet})()})()

    def returnOptions(self, kind):
        return {'Genomes': [self.genome], 'Samples': ['S1', 'S2', 'S3'], 'Species': ['Maylandia zebra'], 'ProjectIDs': ['P'], 'Ecogroups': ['Mbuna', 'Deep_Benthic', 'Exclude'], 'Subgroups': ['x']}[kind]

    def setGenome(self, g):
        self.genome_version = g

    def setSamples(self, projectIDs=None, sampleIDs=None, species=None, ecogroups=None, subgroups=None, rerun=False):
        self.samples = list(sampleIDs or ['S1', 'S2', 'S3'])

    def readSampleDatabase(self): pass
    def readAlignmentDatabase(self): pass

    def createSampleFiles(self, s):
        t = self.tmp
        self.localSampleBamDir, self.localSampleTempDir = f'{t}/bam/{s}/', f'{t}/temp/{s}/'
        self.localRawBamFiles = [f'{t}/reads/{r}.bam' for r in self.reads_dt[self.reads_dt.SampleID == s].RunID]
        self.localBamFile, self.localTempSortedBamFile, self.localChimericBamFile = self.localSampleBamDir + s + '.all.bam', self.localSampleTempDir + s + '.sorted.bam', self.localSampleBamDir + s + '.chimeric.bam'

    def uploadData(self, path, parallel=False):
        return ['rclone', 'copy', path, 'cloud:' + path]

    def downloadData(self, path, parallel=False):
        self.downloads.append(path)


class Proc:
    def __init__(self, out, rc=0):
        self.stdout, self.stderr, self.returncode = out, '', rc


def fake_run(cmd, **kw):
    if cmd[0] == 'git':
        return Proc('abc1234\n')
    return Proc({'minimap2': '2.28-r1209\n', 'samtools': 'samtools 1.21\nUsing htslib 1.21\n', 'gatk': 'The Genome Analysis Toolkit (GATK) v4.6.2.0\nHTSJDK Version: 4.1\n'}[cmd[0]])


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.fm = FakeFM(self.tmp)
        self.lines = []

    def main(self, argv, **kw):
        kw.setdefault('which', lambda t: '/usr/bin/' + t); kw.setdefault('run', fake_run); kw.setdefault('pysam_version', '0.24.1'); kw.setdefault('preflight', lambda cfg: (True, 'gsdf 5d71ee8d6c2b\nte 2af781e90119'))
        return cli.main(argv, fm=kw.pop('fm', self.fm), out=self.lines.append, **kw)

    @property
    def text(self):
        return '\n'.join(self.lines)


class TestHelpers(unittest.TestCase):
    def test_the_fingerprint_ignores_order_and_changes_with_the_runs_and_their_files(self):
        a = cli.input_fingerprint([('R1', 1e9), ('R2', 2e9)])
        self.assertEqual(a, cli.input_fingerprint([('R2', 2e9), ('R1', 1e9)]), 'order does not matter')
        self.assertEqual(len(a), 12)
        self.assertNotEqual(a, cli.input_fingerprint([('R1', 1e9)]), 'a run removed'); self.assertNotEqual(a, cli.input_fingerprint([('R1', 1e9), ('R2', 2e9), ('R3', 1)]), 'a run added')
        self.assertNotEqual(a, cli.input_fingerprint([('R1', 1e9), ('R2', 2e9 + 1)]), 'a uBAM re-made under the same RunID has another size')
        self.assertEqual(cli.input_fingerprint([('R1', NAN)]), cli.input_fingerprint([('R1', None)]), 'an unknown size is one value, whatever it is written as')
        self.assertEqual(cli.input_fingerprint([('R1', 1e9)]), cli.input_fingerprint([('R1', '1000000000')]), 'a size read as text or as a float is the same size')

    def test_sample_meta_collects_what_the_ledger_row_needs(self):
        fm = FakeFM(tempfile.mkdtemp())
        m = cli.sample_meta(fm.sample_dt, fm.reads_dt, ['S1', 'S2', 'S3', 'NOPE'])
        self.assertEqual((m['S1']['Species'], m['S1']['Instrument'], m['S1']['ReadLength'], m['S1']['RunIDs'], m['S1']['NRuns']), ('Maylandia zebra', 'HiSeq 2500', 250, 'R1', 1))
        self.assertEqual((m['S2']['RunIDs'], m['S2']['NRuns']), ('R2b,,R2a', 2), 'the existing rows join RunIDs with a double comma')
        self.assertEqual((m['S3']['Species'], m['S3']['Instrument'], m['S3']['ReadLength']), ('', '', 150), 'unknown values are empty, not NaN')
        self.assertEqual((m['NOPE']['RunIDs'], m['NOPE']['NRuns'], m['NOPE']['ReadLength']), ('', 0, ''))
        self.assertEqual(m['S2']['InputFingerprint'], cli.input_fingerprint([('R2a', 2e9), ('R2b', 3e9)]))

    def test_versions_are_short_like_conda_lists_them(self):
        self.assertEqual([cli.short_version(x) for x in ('2.28-r1209', 'samtools 1.21', 'The Genome Analysis Toolkit (GATK) v4.6.2.0', '', 'unknown')], ['2.28', '1.21', '4.6.2.0', '', 'unknown'])
        self.assertEqual(cli.tool_versions(run=fake_run, pysam_version='0.24.1'), {'minimap2': '2.28', 'samtools': '1.21', 'gatk': '4.6.2.0', 'pysam': '0.24.1'})

    def test_the_pipeline_version_carries_the_commit_when_there_is_one(self):
        nogit = tempfile.mkdtemp()                                                       # a folder that is in no checkout
        self.assertEqual(cli.pipeline_version(run=fake_run, path=nogit), 'align_pipeline 3.0 abc1234')
        self.assertEqual(cli.pipeline_version(run=lambda *a, **k: Proc('', 128), path=nogit), 'align_pipeline 3.0')
        def no_git(*a, **k):
            raise FileNotFoundError('git')
        self.assertEqual(cli.pipeline_version(run=no_git, path=nogit), 'align_pipeline 3.0')


class TestReadingTheCommitWithoutGit(unittest.TestCase):
    """An environment may have no git program (yellowhead's `align` environment has none): the commit is then read from the .git folder."""
    SHA = '4d55ab4f0123456789abcdef0123456789abcdef'

    def checkout(self, head, refs=None, packed=None):
        root = tempfile.mkdtemp(); os.makedirs(root + '/.git/refs/heads'); os.makedirs(root + '/03_align_pipeline/align_pipeline')
        open(root + '/.git/HEAD', 'w').write(head)
        for name, sha in (refs or {}).items():
            open(root + '/.git/' + name, 'w').write(sha + '\n')
        if packed:
            open(root + '/.git/packed-refs', 'w').write('# pack-refs with: peeled fully-peeled sorted\n' + packed)
        return root

    def test_a_branch_a_detached_head_and_packed_refs(self):
        root = self.checkout('ref: refs/heads/reorg_spec\n', refs={'refs/heads/reorg_spec': self.SHA})
        self.assertEqual(cli.read_git_head(root + '/03_align_pipeline/align_pipeline'), '4d55ab4', 'found from a folder deep inside the checkout')
        self.assertEqual(cli.read_git_head(self.checkout(self.SHA + '\n')), '4d55ab4', 'a detached HEAD holds the commit itself')
        packed = self.checkout('ref: refs/heads/reorg_spec\n', packed=f'{self.SHA} refs/heads/reorg_spec\n{"1" * 40} refs/heads/other\n')
        self.assertEqual(cli.read_git_head(packed), '4d55ab4', 'a branch that exists only in packed-refs')

    def test_a_linked_worktree_keeps_its_refs_in_the_main_git_folder(self):
        main = self.checkout('ref: refs/heads/main\n', refs={'refs/heads/reorg_spec': self.SHA})
        wt = tempfile.mkdtemp(); gitdir = main + '/.git/worktrees/wt'; os.makedirs(gitdir)
        open(gitdir + '/HEAD', 'w').write('ref: refs/heads/reorg_spec\n'); open(gitdir + '/commondir', 'w').write('../..\n')
        open(wt + '/.git', 'w').write('gitdir: ' + gitdir + '\n')
        self.assertEqual(cli.read_git_head(wt), '4d55ab4')

    def test_no_checkout_a_broken_checkout_or_a_missing_ref_give_an_empty_string(self):
        self.assertEqual(cli.read_git_head(tempfile.mkdtemp()), '')
        self.assertEqual(cli.read_git_head(self.checkout('ref: refs/heads/nowhere\n')), '')
        self.assertEqual(cli.read_git_head(self.checkout('not a head\n')), '')
        self.assertEqual(cli.read_git_head(self.checkout('ref: refs/heads/x\n', refs={'refs/heads/x': 'not a sha'})), '')

    def test_git_is_asked_first_and_the_folder_is_read_when_git_is_missing_or_refuses(self):
        root = self.checkout('ref: refs/heads/reorg_spec\n', refs={'refs/heads/reorg_spec': self.SHA})
        self.assertEqual(cli.git_commit(root, run=lambda *a, **k: Proc('abc1234\n')), 'abc1234', 'git, when there is one, is the authority')
        def no_git(*a, **k):
            raise FileNotFoundError('git')
        self.assertEqual(cli.git_commit(root, run=no_git), '4d55ab4')
        self.assertEqual(cli.git_commit(root, run=lambda *a, **k: Proc('', 128)), '4d55ab4', 'for example "dubious ownership": a checkout made by another user')
        self.assertEqual(cli.pipeline_version(run=no_git, path=root), 'align_pipeline 3.0 4d55ab4')


class TestDefaults(Base):
    def test_the_disk_reservation_default_is_above_the_measured_peak(self):
        args = cli.build_parser(self.fm).parse_args(['run', 'G'])
        self.assertEqual(args.disk_reserve_factor, 5.5); self.assertGreater(args.disk_reserve_factor, 4.6)


class TestSetupDatabase(Base):
    def test_it_adds_the_ledger_columns_once(self):
        self.fm = FakeFM(self.tmp + '/a', header=OLD) if False else FakeFM(tempfile.mkdtemp(), header=OLD)
        self.assertEqual(self.main(['setup-database', 'G']), 0)
        self.assertIn('Species', self.text); self.assertIn('DuplicationRate', self.text); self.assertIn('PanelStatus', self.text)
        self.assertEqual(self.fm.sheet.row_values(1)[:16], OLD)
        self.lines.clear(); self.assertEqual(self.main(['setup-database', 'G']), 0); self.assertIn('Nothing to do', self.text)


class TestListAndPlan(Base):
    def test_list_skips_flagged_samples_and_reports_the_data(self):
        self.assertEqual(self.main(['run', 'G', '--list']), 0)
        self.assertIn('Skipping 1 sample(s) with Ecogroup Error/Exclude', self.text); self.assertIn('S3', self.text.split('\n')[0])
        self.assertIn('2 samples, 3 uBAM files, 6 GB of uBAMs', self.text); self.assertIn('S1,S2', self.text)
        self.lines.clear(); self.main(['run', 'G', '--list', '--include-excluded']); self.assertIn('3 samples, 4 uBAM files', self.text)

    def test_the_commands_of_a_sample_are_chimeric_only_with_the_three_slices_and_the_panel(self):
        self.assertEqual(self.main(['run', 'G', '--print-commands', 'S1', '--slices', SHIPPED]), 0)
        t = self.text
        self.assertEqual([l[4:] for l in t.split('\n') if l.startswith('--- ')][:2], ['align', 'markdup'])
        self.assertIn('bam_split.py', t); self.assertNotIn('discordant', t); self.assertNotIn('unmapped.bam', t)
        for name, chrom, end in (('LG2', 'NC_135168.1', 45_000_000), ('LG10', 'NC_135176.1', 34_949_860), ('LG13', 'NC_135179.1', 30_000_000)):
            self.assertIn(f'--- slice {name}', t); self.assertIn(f'{chrom}:1-{end}', t)
        self.assertIn('-m sd_genotypes panel measure --bam', t); self.assertIn('--sample S1', t)
        self.assertEqual(self.fm.downloads, [], 'printing the commands downloads nothing')

    def test_the_slice_table_of_the_genome_is_found_by_default(self):
        self.fm = FakeFM(tempfile.mkdtemp(), genome='Mzebra_GT3_NCBI')
        self.assertEqual(self.main(['run', 'Mzebra_GT3_NCBI', '--print-commands', 'S1']), 0)
        for name in ('LG2', 'LG10', 'LG13'):
            self.assertIn(f'--- slice {name}', self.text)
        for name in ('LG9', 'LG11', 'LG20'):
            self.assertNotIn(f'--- slice {name}', self.text)                       # decision 42

    def test_no_panel_and_no_slices_remove_those_commands(self):
        self.main(['run', 'G', '--print-commands', 'S1', '--no-panel', '--no-slices'])
        self.assertNotIn('sd_genotypes', self.text); self.assertNotIn('--- slice', self.text)

    def test_the_panel_python_and_checkout_reach_the_command(self):
        self.main(['run', 'G', '--print-commands', 'S1', '--panel-python', '/env/sdpanel/bin/python', '--panel-pythonpath', '/repo/04_genotype_sd_systems'])
        self.assertIn('env PYTHONPATH=/repo/04_genotype_sd_systems', self.text); self.assertIn('/env/sdpanel/bin/python -m sd_genotypes panel measure', self.text)

    def test_an_unknown_sample_is_reported(self):
        self.assertEqual(self.main(['run', 'G', '--print-commands', 'S2', '-s', 'S1']), 1); self.assertIn('S2 is not among the selected samples', self.text)


class TestGates(Base):
    """Everything that can be known before the first alignment is checked before it, so no hours of work are lost to it."""
    def run_main(self, extra=(), **kw):
        with mock.patch.object(sp.SeqPipeline, 'run_batch', autospec=True, return_value=('1 of 2 samples complete', ['S1', 'S2'], {})) as rb:
            rc = self.main(['run', 'G'] + list(extra), **kw)
        return rc, rb

    def test_a_panel_that_is_not_ready_stops_the_run_with_the_reason_and_the_way_out(self):
        rc, rb = self.run_main(preflight=lambda cfg: (False, 'te: definition file missing; build it with `te.definition`'))
        self.assertEqual(rc, 1); rb.assert_not_called()
        self.assertIn('The panel is not ready, so nothing was started', self.text); self.assertIn('te: definition file missing', self.text); self.assertIn('--no-panel', self.text)

    def test_expected_panel_keys_reach_the_preflight_and_a_different_key_stops_the_run(self):
        real = lambda cfg: pn.preflight(cfg, run=lambda cmd, **kw: Proc('gsdf 5d71ee8d6c2b\nte ffffffffffff\n'), timeout=1)
        rc, rb = self.run_main(['--panel-expect', '5d71ee8d6c2b', '2af781e90119'], preflight=real)
        self.assertEqual(rc, 1); rb.assert_not_called(); self.assertIn('not the expected ones', self.text); self.assertIn('2af781e90119', self.text)
        self.lines.clear()
        rc, rb = self.run_main(['--panel-expect', '5d71ee8d6c2b'], preflight=real); self.assertEqual(rc, 0); rb.assert_called_once()

    def test_no_panel_does_not_ask_the_panel_anything(self):
        def must_not_be_called(cfg):
            raise AssertionError('preflight called')
        rc, rb = self.run_main(['--no-panel'], preflight=must_not_be_called)
        self.assertEqual(rc, 0); rb.assert_called_once()

    def test_a_missing_tool_stops_the_run(self):
        rc, rb = self.run_main(which=lambda t: None if t == 'rclone' else '/usr/bin/' + t)
        self.assertEqual(rc, 1); rb.assert_not_called(); self.assertIn('rclone not found on PATH', self.text)

    def test_a_worksheet_without_the_ledger_columns_stops_the_run_with_the_setup_command(self):
        self.fm = FakeFM(tempfile.mkdtemp(), header=OLD)
        rc, rb = self.run_main()
        self.assertEqual(rc, 1); rb.assert_not_called(); self.assertIn('python -m align_pipeline setup-database G', self.text); self.assertIn('Species', self.text)

    def test_a_slice_table_that_does_not_fit_the_reference_stops_the_run(self):
        bad = self.tmp + '/bad.tsv'; open(bad, 'w').write('LG99\tNC_999.1\t1\t100\n')
        rc, rb = self.run_main(['--slices', bad])
        self.assertEqual(rc, 1); rb.assert_not_called(); self.assertIn('does not fit the reference', self.text); self.assertIn('NC_999.1', self.text)
        self.lines.clear(); bad2 = self.tmp + '/bad2.tsv'; open(bad2, 'w').write('LG99\tNC_999.1\t1\tEND\n')
        rc, rb = self.run_main(['--slices', bad2]); self.assertEqual(rc, 1); self.assertIn('not in the reference index', self.text)

    def test_the_happy_path_builds_the_pipeline_with_everything_the_row_needs(self):
        rc, rb = self.run_main(['--no-upload', '--slices', SHIPPED])
        self.assertEqual(rc, 0)
        pipe = rb.call_args[0][0]
        self.assertEqual([s.name for s in pipe.slices], ['LG2', 'LG10', 'LG13']); self.assertIsInstance(pipe.cfg.panel, pn.PanelConfig)
        self.assertEqual(pipe.versions, {'minimap2': '2.28', 'samtools': '1.21', 'gatk': '4.6.2.0', 'pysam': '0.24.1'}); self.assertEqual(pipe.pipeline_version, 'align_pipeline 3.0 abc1234')
        self.assertEqual(pipe.sample_meta['S2']['NRuns'], 2); self.assertEqual(pipe.run_ids, {'S1': 'R1', 'S2': 'R2b,,R2a'}); self.assertEqual(pipe.read_lengths, {'S1': 250, 'S2': 250, 'S3': 150})
        self.assertEqual(pipe.ref_size, sum(FAI.values())); self.assertEqual(pipe.genome_version, 'G')
        self.assertEqual(self.fm.downloads, [self.fm.localGenomeDir], 'the genome data is fetched as before')
        self.assertIn('panel ready; definition keys', self.text); self.assertIn('gsdf 5d71ee8d6c2b', self.text)

    def test_a_failed_sample_makes_the_exit_code_one(self):
        with mock.patch.object(sp.SeqPipeline, 'run_batch', autospec=True, return_value=('summary', ['S1'], {'S2': 'align: exit code 1'})):
            self.assertEqual(self.main(['run', 'G', '--no-upload']), 1)

    def test_no_samples_is_reported(self):
        self.fm.merged_dt = self.fm.merged_dt[self.fm.merged_dt.SampleID == 'S3']
        self.assertEqual(self.main(['run', 'G', '-s', 'S3']), 1); self.assertIn('No samples to run', self.text)


if __name__ == '__main__':
    unittest.main()
