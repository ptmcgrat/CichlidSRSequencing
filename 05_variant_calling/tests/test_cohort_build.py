"""`cohort build` (module 05): command construction, the plan and preflight, and the whole build, run with FAKE bcftools / samtools / rclone (fake_tools.py).
The tests establish that the builder issues the intended command lines and handles staging, resume, failure, renaming, the manifest and its safety stops. They do not establish how the
real tools behave: that is what the first real run on yellowhead, with a small trial, is for."""
import json, os, shutil, stat, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from cohort_fixtures import *
from variant_calling import build as B, cohort_cli

HERE = os.path.dirname(os.path.abspath(__file__))
CFG = B.load_build_config()


class BuildBase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        bindir = os.path.join(self.tmp, 'bin'); os.makedirs(bindir)
        for tool in ('bcftools', 'samtools', 'rclone'):
            p = os.path.join(bindir, tool)
            open(p, 'w').write(f"#!{sys.executable}\nimport sys\nsys.path.insert(0, {HERE!r})\nimport fake_tools\nsys.exit(fake_tools.main({tool!r}, sys.argv[1:]))\n")
            os.chmod(p, os.stat(p).st_mode | stat.S_IEXEC)
        self.saved = {k: os.environ.get(k) for k in ('PATH', 'FAKE_LOG', 'FAKE_NO_GROUPS', 'FAKE_NO_IGNORE_RG', 'FAKE_FAIL_CONTIG', 'FAKE_FAIL_FLAG', 'FAKE_EXTRA_COLUMNS', 'FAKE_RECORDS')}
        os.environ['PATH'] = bindir + os.pathsep + os.environ['PATH']
        self.logfile = os.path.join(self.tmp, 'calls.log'); os.environ['FAKE_LOG'] = self.logfile
        for k in ('FAKE_NO_GROUPS', 'FAKE_NO_IGNORE_RG', 'FAKE_FAIL_CONTIG', 'FAKE_FAIL_FLAG', 'FAKE_EXTRA_COLUMNS', 'FAKE_RECORDS'):
            os.environ.pop(k, None)
        self.local, self.cloud = os.path.join(self.tmp, 'local') + '/', os.path.join(self.tmp, 'cloud') + '/'
        os.makedirs(self.local)
        self.fa = os.path.join(self.tmp, 'ref.fa'); open(self.fa, 'w').write('>NC_1\nACGT\n')
        open(self.fa + '.fai', 'w').write('NC_1\t20000000\t0\t60\t61\nNC_2\t15000000\t0\t60\t61\nNC_3\t8000000\t0\t60\t61\nNC_MT\t16000\t0\t60\t61\nNW_1\t900000\t0\t60\t61\n')
        self.env = B.Env(self.local, self.cloud, G, self.fa)
        rows = brood()
        self.dt = sheet(rows)
        self.ids = [r['SampleID'] for r in rows]
        led = []
        for i, s in enumerate(self.ids):
            d = f'{self.cloud}Bamfiles/{G}/{s}/'; os.makedirs(d)
            open(d + f'{s}.all.bam', 'w').write('BAM' * (100 + 10 * i)); open(d + f'{s}.all.bai', 'w').write('BAI')
            led.append(dict(SampleID=s, GenomeVersion=G, Coverage=20.0, BamSize=os.path.getsize(d + f'{s}.all.bam')))
        self.ledger = pd.DataFrame(led)
        self.spec = family_spec()
        self.msgs, self.uploads = [], []

    def tearDown(self):
        for k, v in self.saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        shutil.rmtree(self.tmp, ignore_errors=True)

    def build(self, **kw):
        kw.setdefault('run_it', True)
        return B.build_cohort('fam', kw.pop('spec', self.spec), kw.pop('dt', self.dt), kw.pop('ledger', self.ledger), self.env, cfg=kw.pop('cfg', CFG), log=self.msgs.append, free_bytes=kw.pop('free_bytes', 1e12),
                              uploader=self.uploads.append, cloud_exists=kw.pop('cloud_exists', lambda p: False), **kw)

    def calls(self, tool, sub=None):
        out = []
        for l in open(self.logfile) if os.path.exists(self.logfile) else []:
            e = json.loads(l)
            if e['tool'] == tool and (sub is None or (e['argv'] and e['argv'][0] == sub)):
                out.append(e['argv'])
        return out


class TestCommands(unittest.TestCase):
    def test_the_pileup_follows_blumers_parameters(self):
        self.assertEqual(B.mpileup_cmd('ref.fa', 'bams.list', 'NC_1', CFG), ['bcftools', 'mpileup', '-f', 'ref.fa', '-b', 'bams.list', '-r', 'NC_1', '--ignore-RG', '-q', '13', '-Q', '20',
                                                                            '-p', '-m', '2', '-F', '0.05', '-a', 'FORMAT/AD,FORMAT/DP', '-Ou'])

    def test_the_caller_has_no_shared_hardy_weinberg_prior(self):
        self.assertEqual(B.call_cmd('out.vcf.gz', CFG), ['bcftools', 'call', '-m', '-v', '-P', '0.001', '-G', '-', '-f', 'GQ', '-Oz', '-o', 'out.vcf.gz'])

    def test_other_group_modes_are_refused_for_now(self):
        with self.assertRaises(B.BuildError) as cm:
            B.call_cmd('o.vcf.gz', CFG, groups='species')
        self.assertIn('per_sample', str(cm.exception))

    def test_a_contig_script_pipes_checks_both_tools_and_indexes(self):
        s = B.contig_script('ref.fa', 'bams.list', 'NC_1:1-100', 'o.vcf.gz', CFG)
        self.assertTrue(s.startswith('set -o pipefail; bcftools mpileup')); self.assertIn(' | bcftools call ', s); self.assertTrue(s.endswith('&& bcftools index -t o.vcf.gz'))

    def test_the_parameters_come_from_the_config_and_can_be_changed(self):
        cfg = B.load_build_config(); cfg['bcftools'].update(prior=0.0001, indel_min_reads=0, mpileup_annotate=[], call_format_fields=[])
        self.assertNotIn('-p', B.mpileup_cmd('r', 'b', 'c', cfg)); self.assertNotIn('-a', B.mpileup_cmd('r', 'b', 'c', cfg))
        self.assertEqual(B.call_cmd('o', cfg), ['bcftools', 'call', '-m', '-v', '-P', '0.0001', '-G', '-', '-Oz', '-o', 'o'])

    def test_regions_and_contigs(self):
        self.assertEqual(B.parse_region('NC_1'), ('NC_1', None, None)); self.assertEqual(B.parse_region('NC_1:1,000-2,000'), ('NC_1', 1000, 2000))
        with self.assertRaises(B.BuildError):
            B.parse_region('NC_1:abc')
        fai = [('NC_1', 20_000_000), ('NC_2', 6_000_000), ('NC_MT', 16_000), ('NW_1', 900_000), ('NC_3', 1_000_000)]
        self.assertEqual([c for c, _, _ in B.choose_contigs(fai, 'chromosomes', 5_000_000)], ['NC_1', 'NC_2'])
        self.assertEqual(len(B.choose_contigs(fai, 'all', 0)), 5)
        self.assertEqual(B.choose_contigs(fai, 'chromosomes', 5_000_000, region='NC_2:10-20'), [('NC_2', 10, 20)])
        self.assertEqual([c for c, _, _ in B.choose_contigs(fai, 'NC_3,NW_1', 0)], ['NC_3', 'NW_1'])
        for bad in (dict(mode='NC_9'), dict(mode='chromosomes', region='NC_9')):
            with self.assertRaises(B.BuildError):
                B.choose_contigs(fai, bad['mode'], 0, region=bad.get('region'))
        with self.assertRaises(B.BuildError):
            B.choose_contigs([('NW_1', 100)], 'chromosomes', 5_000_000)


class TestPlan(BuildBase):
    def test_the_plan_lists_the_bams_contigs_commands_and_outputs_and_writes_nothing(self):
        r = self.build(run_it=False)
        self.assertEqual(r['status'], 'plan'); self.assertEqual(r['problems'], [])
        p = r['plan']; self.assertEqual(len(p.members), 26); self.assertEqual([c[0] for c in p.contigs], ['NC_1', 'NC_2', 'NC_3'], 'the mitochondrion and the unplaced piece are left out')
        self.assertEqual(p.total_bytes, sum(os.path.getsize(f'{self.cloud}Bamfiles/{G}/{s}/{s}.all.bam') for s in self.ids))
        t = r['text']
        for want in ('26 samples (family)', 'caller bcftools, groups per_sample', 'bcftools mpileup', '-G -', 'fam.vcf.gz', 'fam.snps.vcf.gz', 'fam.manifest.json', 'uploaded and verified', 'PLAN ONLY'):
            self.assertIn(want, t)
        self.assertFalse(os.path.exists(p.work), 'a plan stages and writes nothing')

    def test_a_blocked_cohort_is_not_planned(self):
        led = self.ledger[self.ledger.SampleID != 'B1_f1']
        r = self.build(run_it=False, ledger=led)
        self.assertEqual(r['status'], 'blocked'); self.assertIn('no_alignment', r['report'].codes('BLOCK'))

    def test_a_trial_takes_the_smallest_bams_stays_local_and_says_so(self):
        r = self.build(run_it=False, max_samples=4)
        p = r['plan']; self.assertTrue(p.trial); self.assertEqual(len(p.members), 4)
        sizes = sorted(p.sizes.values()); self.assertEqual(sizes, sorted(self.ledger.BamSize)[:4])
        self.assertIn('TRIAL', r['text']); self.assertIn('trial_out', p.out_dir); self.assertTrue(p.files['vcf'].endswith('fam.trial.vcf.gz')); self.assertNotIn('Outputs/PedigreeVCFs', p.out_dir)

    def test_a_region_limits_the_plan_to_that_region(self):
        r = self.build(run_it=False, region='NC_2:1-3000000')
        self.assertEqual(r['plan'].contigs, [('NC_2', 1, 3000000)]); self.assertIn('NC_2:1-3000000', r['text'])

    def test_preflight_problems(self):
        r = self.build(run_it=False, which=lambda t: None); self.assertTrue(any('not on the PATH' in p for p in r['problems']))
        os.environ['FAKE_NO_GROUPS'] = '1'; r = self.build(run_it=False); self.assertTrue(any('group-samples' in p for p in r['problems'])); del os.environ['FAKE_NO_GROUPS']
        os.environ['FAKE_NO_IGNORE_RG'] = '1'; r = self.build(run_it=False); self.assertTrue(any('ignore-RG' in p for p in r['problems'])); del os.environ['FAKE_NO_IGNORE_RG']
        r = self.build(run_it=False, free_bytes=1000); self.assertTrue(any('GB free' in p for p in r['problems']))
        os.remove(self.fa + '.fai'); r = self.build(run_it=False); self.assertTrue(any('.fai is missing' in p for p in r['problems']))
        os.remove(self.fa); r = self.build(run_it=False); self.assertTrue(any('FASTA is not at' in p for p in r['problems']))

    def test_members_without_a_bamsize_are_a_preflight_problem(self):
        led = self.ledger.copy(); led['BamSize'] = led['BamSize'].astype(object); led.loc[led.SampleID == 'B1_m0', 'BamSize'] = ''
        r = self.build(run_it=False, ledger=led); self.assertTrue(any('no BamSize' in p for p in r['problems']))

    def test_a_run_with_preflight_problems_starts_nothing(self):
        with self.assertRaises(B.BuildError) as cm:
            self.build(which=lambda t: None)
        self.assertIn('preflight failed, nothing was started', str(cm.exception)); self.assertFalse(os.path.exists(f'{self.local}Cohorts/'))


class TestRun(BuildBase):
    def test_a_build_stages_calls_merges_renames_writes_the_manifest_and_uploads(self):
        r = self.build()
        self.assertEqual(r['status'], 'built')
        p = r['plan']
        head = open(p.files['vcf']).read().splitlines()
        cols = [l for l in head if l.startswith('#CHROM')][0].split('\t')[9:]
        self.assertEqual(cols, p.members, 'the columns are the SampleIDs, in the order of the cohort')
        self.assertEqual(cols[:2], ['B1_m0', 'B1_m1']); self.assertEqual(cols[-2:], ['S_m', 'D_f'], 'children first, then the parents'); self.assertEqual(sorted(cols), sorted(self.ids))
        self.assertEqual(r['counts'], {'vcf': 9, 'snps': 9}, '3 variants on each of 3 contigs')
        self.assertTrue(os.path.exists(p.files['vcf'] + '.tbi') and os.path.exists(p.files['snps'] + '.tbi'))
        man = json.load(open(p.files['manifest']))
        self.assertEqual((man['cohort'], man['kind'], man['genome'], man['caller'], man['groups'], man['trial']), ('fam', 'family', G, 'bcftools', 'per_sample', False))
        self.assertEqual([m['sample'] for m in man['members']], p.members)
        want = dict(zip(self.ledger.SampleID, self.ledger.BamSize)); self.assertTrue(all(m['bam_bytes'] == want[m['sample']] for m in man['members'])); self.assertEqual(man['members'][0]['coverage'], 20.0)
        self.assertEqual(man['parameters']['bcftools']['prior'], 0.001); self.assertEqual([c['contig'] for c in man['contigs']], ['NC_1', 'NC_2', 'NC_3'])
        self.assertEqual(man['outputs']['vcf']['records'], 9); self.assertFalse(man['check']['blocked']); self.assertIn('bcftools 1.99', man['tools']['bcftools'])
        self.assertEqual([os.path.basename(u) for u in self.uploads], ['fam.vcf.gz', 'fam.vcf.gz.tbi', 'fam.snps.vcf.gz', 'fam.snps.vcf.gz.tbi', 'fam.manifest.json'])
        self.assertIn('Outputs/PedigreeVCFs/', p.out_dir)
        self.assertFalse(os.path.exists(p.staged), 'the staged BAMs are deleted'); self.assertFalse(os.path.exists(p.per_contig))

    def test_the_commands_that_reached_the_tools_are_the_intended_ones(self):
        self.build()
        mp = self.calls('bcftools', 'mpileup'); mp = [a for a in mp if len(a) > 1]
        self.assertEqual(len(mp), 3)
        for a in mp:
            for want in ('--ignore-RG', '-q', '13', '-Q', '20', '-p', '-m', '2', '-F', '0.05', '-a', 'FORMAT/AD,FORMAT/DP', '-f', self.fa):
                self.assertIn(want, a)
        cl = [a for a in self.calls('bcftools', 'call') if len(a) > 1]
        self.assertEqual(len(cl), 3)
        for a in cl:
            self.assertEqual(a[a.index('-G') + 1], '-'); self.assertEqual(a[a.index('-P') + 1], '0.001'); self.assertIn('-m', a); self.assertIn('-v', a)
        self.assertEqual(len(self.calls('rclone', 'copy')), 26)

    def test_a_trial_is_local_only_and_never_reaches_the_analysed_folder(self):
        r = self.build(max_samples=4, region='NC_2:1-3000000')
        self.assertEqual(r['status'], 'built'); self.assertEqual(self.uploads, [])
        p = r['plan']; self.assertTrue(os.path.exists(p.files['vcf'])); self.assertIn('trial_out', p.files['vcf'])
        self.assertFalse(os.path.exists(f'{self.local}Outputs/PedigreeVCFs/'), 'nothing was written where `sexloci run-folder` looks')
        man = json.load(open(p.files['manifest'])); self.assertTrue(man['trial']); self.assertEqual(len(man['members']), 4); self.assertEqual(man['n_members_in_cohort'], 26); self.assertEqual(man['region'], 'NC_2:1-3000000')
        self.assertEqual([c['contig'] for c in man['contigs']], ['NC_2_1-3000000'])

    def test_no_upload_writes_locally_only(self):
        r = self.build(upload=False); self.assertEqual(self.uploads, []); self.assertTrue(os.path.exists(r['plan'].files['vcf']))

    def test_keep_staged_and_keep_work_keep_them(self):
        r = self.build(keep_staged=True, keep_work=True)
        self.assertTrue(os.path.isdir(r['plan'].staged)); self.assertTrue(os.path.isdir(r['plan'].per_contig))

    def test_staged_bams_get_both_index_names(self):
        r = self.build(keep_staged=True)
        d = f"{r['plan'].staged}B1_m0/"
        self.assertTrue(os.path.exists(d + 'B1_m0.all.bai')); self.assertTrue(os.path.lexists(d + 'B1_m0.all.bam.bai'))


class TestResumeAndFailure(BuildBase):
    def test_a_failed_contig_keeps_the_others_and_a_rerun_does_only_the_rest(self):
        flag = os.path.join(self.tmp, 'fail.flag'); open(flag, 'w').write('x')
        os.environ['FAKE_FAIL_CONTIG'], os.environ['FAKE_FAIL_FLAG'] = 'NC_2', flag
        with self.assertRaises(B.BuildError) as cm:
            self.build(keep_staged=True)
        msg = str(cm.exception); self.assertIn('1 of 3 contigs failed', msg); self.assertIn('NC_2', msg); self.assertIn('simulated failure', msg); self.assertIn('rerun does only the rest', msg)
        pc = f'{self.local}Cohorts/fam/per_contig/'
        self.assertTrue(os.path.exists(pc + 'NC_1.vcf.gz.done') and os.path.exists(pc + 'NC_3.vcf.gz.done')); self.assertFalse(os.path.exists(pc + 'NC_2.vcf.gz.done'))
        self.assertFalse(os.path.exists(f'{self.local}Outputs/PedigreeVCFs/fam.vcf.gz'), 'no output from a failed build'); self.assertEqual(self.uploads, [])
        os.remove(flag)
        n_before = len([a for a in self.calls('bcftools', 'mpileup') if len(a) > 1]); copies_before = len(self.calls('rclone', 'copy'))
        r = self.build(); self.assertEqual(r['status'], 'built')
        new = [a for a in self.calls('bcftools', 'mpileup') if len(a) > 1][n_before:]
        self.assertEqual([a[a.index('-r') + 1] for a in new], ['NC_2'], 'only the failed contig was called again')
        self.assertEqual(len(self.calls('rclone', 'copy')), copies_before, 'the staged BAMs were reused')
        self.assertEqual(r['counts']['vcf'], 9); self.assertTrue(any('already done' in m for m in self.msgs))
        self.assertEqual([c['resumed'] for c in json.load(open(r['plan'].files['manifest']))['contigs']], [True, False, True])

    def test_a_wrong_number_of_sample_columns_stops_before_anything_is_written(self):
        os.environ['FAKE_EXTRA_COLUMNS'] = '2'
        with self.assertRaises(B.BuildError) as cm:
            self.build()
        self.assertIn('28 sample columns for 26 BAMs', str(cm.exception)); self.assertFalse(os.path.exists(f'{self.local}Outputs/PedigreeVCFs/fam.vcf.gz')); self.assertEqual(self.uploads, [])

    def test_a_damaged_download_is_named(self):
        open(f'{self.cloud}Bamfiles/{G}/B1_f3/B1_f3.all.bam', 'a').write('CORRUPT')
        led = self.ledger.copy(); led.loc[led.SampleID == 'B1_f3', 'BamSize'] = os.path.getsize(f'{self.cloud}Bamfiles/{G}/B1_f3/B1_f3.all.bam')
        with self.assertRaises(B.BuildError) as cm:
            self.build(ledger=led)
        self.assertIn('B1_f3', str(cm.exception)); self.assertIn('damaged', str(cm.exception))

    def test_a_size_that_differs_from_the_ledger_is_an_error(self):
        led = self.ledger.copy(); led.loc[led.SampleID == 'B1_m2', 'BamSize'] += 5
        with self.assertRaises(B.BuildError) as cm:
            self.build(ledger=led)
        self.assertIn('B1_m2', str(cm.exception)); self.assertIn('ledger says', str(cm.exception))

    def test_a_bam_missing_from_dropbox_is_an_error(self):
        os.remove(f'{self.cloud}Bamfiles/{G}/B1_m4/B1_m4.all.bam')
        with self.assertRaises(B.BuildError) as cm:
            self.build()
        self.assertIn('B1_m4', str(cm.exception))


class TestOverwriteGuards(BuildBase):
    def test_existing_local_outputs_are_not_replaced_without_force(self):
        self.build()
        with self.assertRaises(B.BuildError) as cm:
            self.build()
        self.assertIn('already exists', str(cm.exception)); self.assertIn('--force', str(cm.exception))
        self.assertEqual(self.build(force=True)['status'], 'built')

    def test_an_output_already_on_dropbox_is_not_replaced_without_force(self):
        with self.assertRaises(B.BuildError) as cm:
            self.build(cloud_exists=lambda p: True)
        self.assertIn('already on Dropbox', str(cm.exception)); self.assertFalse(os.path.exists(f'{self.local}Cohorts/fam/staged/')), 'nothing was staged before the refusal'
        self.assertEqual(self.build(cloud_exists=lambda p: True, force=True)['status'], 'built')

    def test_a_blocked_cohort_does_no_work(self):
        led = self.ledger[self.ledger.SampleID != 'B1_f1']
        r = self.build(ledger=led); self.assertEqual(r['status'], 'blocked'); self.assertFalse(os.path.exists(f'{self.local}Cohorts/'))


class TestCli(BuildBase):
    def setUp(self):
        super().setUp()
        self.cohorts = os.path.join(self.tmp, 'c.yaml'); open(self.cohorts, 'w').write("callsets:\n  fam:\n    kind: family\n    broods: [B1]\n")
        self.sample_tsv, self.ledger_tsv = os.path.join(self.tmp, 's.tsv'), os.path.join(self.tmp, 'l.tsv')
        self.dt.to_csv(self.sample_tsv, sep='\t', index=False); self.ledger.to_csv(self.ledger_tsv, sep='\t', index=False)

    def call(self, *extra):
        import contextlib, io
        out, err = io.StringIO(), io.StringIO()
        args = ['build', 'fam', '--cohorts', self.cohorts, '--sample-table', self.sample_tsv, '--alignment-table', self.ledger_tsv, '--local-master', self.local, '--cloud-master', self.cloud, '--fasta', self.fa] + list(extra)
        with contextlib.redirect_stderr(err):
            rc = cohort_cli.main(args, out=out)
        return rc, out.getvalue(), err.getvalue()

    def test_the_default_is_a_plan_and_exits_0(self):
        rc, out, _ = self.call(); self.assertEqual(rc, 0); self.assertIn('PLAN ONLY', out); self.assertFalse(os.path.exists(f'{self.local}Cohorts/'))

    def test_a_plan_with_a_preflight_problem_exits_3(self):
        os.environ['FAKE_NO_GROUPS'] = '1'
        rc, out, _ = self.call(); self.assertEqual(rc, 3); self.assertIn('[PREFLIGHT]', out)

    def test_a_blocked_cohort_exits_1_and_shows_the_report(self):
        pd.read_csv(self.ledger_tsv, sep='\t').iloc[1:].to_csv(self.ledger_tsv, sep='\t', index=False)
        rc, out, _ = self.call(); self.assertEqual(rc, 1); self.assertIn('no_alignment', out)

    def test_run_builds_and_exits_0_and_a_build_error_exits_4(self):
        rc, out, _ = self.call('--run', '--no-upload'); self.assertEqual(rc, 0, out); self.assertIn('[done] fam', out); self.assertTrue(os.path.exists(f'{self.local}Outputs/PedigreeVCFs/fam.vcf.gz'))
        rc, out, err = self.call('--run', '--no-upload'); self.assertEqual(rc, 4); self.assertIn('BUILD ERROR', err); self.assertIn('--force', err)

    def test_the_offline_paths_go_together(self):
        with self.assertRaises(SystemExit):
            cohort_cli.main(['build', 'fam', '--cohorts', self.cohorts, '--sample-table', self.sample_tsv, '--alignment-table', self.ledger_tsv, '--local-master', self.local], out=None)


if __name__ == '__main__':
    unittest.main()
