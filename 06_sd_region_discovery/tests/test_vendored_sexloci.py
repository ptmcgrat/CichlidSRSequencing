"""sexloci 0.1.26, vendored UNCHANGED into this module (SPEC decision 58).
  * the manifest: every file is the one that was copied in, and nothing was added;
  * the package imports and has the version of the manifest;
  * a smoke test with sexloci's own simulator: a pedigree of 3 families with an XY locus planted on chr2 (1.25 to 1.75 Mb) is recovered at tier A, and a wrong genome version is refused;
  * its own pytest suite (105 tests), when pytest is installed (it is not in `envs/analysis.yml`: use `envs/sex_mapping.yml`).
Not exercised here, as in sexloci's own README: the scikit-allel backend, the live rclone and the live Google sheets."""
import contextlib, hashlib, io, os, subprocess, sys, tempfile, unittest

os.environ.setdefault('MPLBACKEND', 'Agg')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def manifest():
    rows = {}
    for line in open(os.path.join(HERE, 'SEXLOCI_MANIFEST.tsv')):
        if line.startswith('#') or line.startswith('path\t') or not line.strip():
            continue
        path, size, md5 = line.rstrip('\n').split('\t')
        rows[path] = (int(size), md5)
    return rows


class TestManifest(unittest.TestCase):
    def test_every_vendored_file_is_the_one_that_was_copied(self):
        for path, (size, md5) in manifest().items():
            p = os.path.join(HERE, path)
            self.assertTrue(os.path.exists(p), f'{path} is missing')
            data = open(p, 'rb').read()
            self.assertEqual((len(data), hashlib.md5(data).hexdigest()), (size, md5), f'{path} was changed: change sexloci in its own project and copy it in again')

    def test_nothing_was_added_to_the_vendored_folders(self):
        listed = set(manifest())
        for top in ('sexloci', 'sexloci_dist'):
            for root, dirs, files in os.walk(os.path.join(HERE, top)):
                dirs[:] = [d for d in dirs if d != '__pycache__']
                for f in files:
                    if f.endswith('.pyc'):
                        continue
                    rel = os.path.relpath(os.path.join(root, f), HERE)
                    self.assertIn(rel, listed, f'{rel} is not in SEXLOCI_MANIFEST.tsv')

    def test_the_manifest_covers_the_package_and_its_documents(self):
        paths = set(manifest())
        for must in ('sexloci/cli.py', 'sexloci/pipeline.py', 'sexloci/vcf.py', 'sexloci_dist/docs/METHODS.md', 'sexloci_dist/README.md', 'sexloci_dist/tests/test_sexloci.py'):
            self.assertIn(must, paths)


class TestPackage(unittest.TestCase):
    def test_it_imports_with_the_manifest_version(self):
        import sexloci
        self.assertEqual(sexloci.__version__, '0.1.26')
        self.assertIn('0.1.26', open(os.path.join(HERE, 'SEXLOCI_MANIFEST.tsv')).readline())

    def test_the_module_command_line_delegates_to_sexloci(self):
        from sd_regions import cli
        out = io.StringIO()
        with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as cm:
            cli.main(['sexloci', '--help'])
        self.assertEqual(cm.exception.code, 0); self.assertIn('run-folder', out.getvalue())


class TestSimulatedLocus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import sexloci.cli as scli
        cls.scli = scli
        cls.dir = tempfile.mkdtemp()
        with contextlib.redirect_stdout(io.StringIO()):
            scli.main(['simulate', '--outdir', cls.dir, '--system', 'XY', '--mode', 'pedigree', '--seed', '3'])

    def args(self, out, genome):
        return ['run', '--vcf', f'{self.dir}/sim.vcf', '--sample-table', f'{self.dir}/sample_table.tsv', '--alignment-table', f'{self.dir}/alignment_table.tsv',
                '--reads-table', f'{self.dir}/reads_table.tsv', '--outdir', out, '--n-perm', '100', '--genome-version', genome]

    def test_a_planted_xy_locus_is_recovered_at_tier_a(self):
        import json, pandas as pd
        truth = json.load(open(f'{self.dir}/truth.json'))
        out = tempfile.mkdtemp()
        with contextlib.redirect_stdout(io.StringIO()):
            rc = self.scli.main(self.args(out, 'Mconophoros_GT1'))
        self.assertEqual(rc, 0)
        reg = pd.read_csv(os.path.join(out, 'sex_linked_regions.tsv'), sep='\t')
        self.assertEqual(len(reg), 1, reg.to_string())
        r = reg.iloc[0]
        self.assertEqual((r.chrom, r.tier, r.model), (truth['chrom'], 'A', 'XY'))
        self.assertTrue(r.start <= truth['sex_pos'] <= r.end, 'the planted position lies inside the called region')
        self.assertEqual((r.units_supporting, r.units_total), (3, 3), 'all three independent families support it')

    def test_a_wrong_genome_version_is_refused_by_the_preflight(self):
        out = tempfile.mkdtemp()
        with contextlib.redirect_stdout(io.StringIO()):
            rc = self.scli.main(self.args(out, 'Mzebra_GT3_NCBI'))
        self.assertEqual(rc, 2)
        self.assertTrue(os.path.exists(os.path.join(out, 'preflight_report.md')))
        self.assertIn('genome_version_mismatch', open(os.path.join(out, 'preflight_report.md')).read())
        self.assertFalse(os.path.exists(os.path.join(out, 'sex_linked_regions.tsv')))


try:
    import pytest as _pytest
except ImportError:
    _pytest = None


@unittest.skipIf(_pytest is None, "pytest is not installed: run `python -m pytest 06_sd_region_discovery/sexloci_dist/tests` in an environment from envs/sex_mapping.yml")
class TestItsOwnSuite(unittest.TestCase):
    def test_sexlocis_own_105_tests_pass(self):
        env = dict(os.environ, PYTHONPATH=os.pathsep.join([HERE] + ([os.environ['PYTHONPATH']] if os.environ.get('PYTHONPATH') else [])), MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1')
        r = subprocess.run([sys.executable, '-m', 'pytest', '-q', '-x', '-p', 'no:cacheprovider', os.path.join(HERE, 'sexloci_dist', 'tests')], capture_output=True, text=True, env=env, cwd=HERE, timeout=1800)
        self.assertEqual(r.returncode, 0, r.stdout[-2500:] + r.stderr[-1500:])


if __name__ == '__main__':
    unittest.main()
