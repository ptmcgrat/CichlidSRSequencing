#!/usr/bin/env python3
"""Replace the vendored sexloci by a new release, unchanged, and rewrite SEXLOCI_MANIFEST.tsv.

    python 06_sd_region_discovery/vendor_sexloci.py ~/Downloads/sexloci-0_1_27.zip [--date 2026-12-01]

The zip holds `sexloci/` with `src/sexloci/` (the package) and the rest of the distribution. `src/sexloci/*.py` becomes `06_sd_region_discovery/sexloci/`; README.md,
CHANGELOG.md, LICENSE, pyproject.toml, environment.yml, config/, docs/, examples/, scripts/ and tests/ become `06_sd_region_discovery/sexloci_dist/`. `.github/` and `.gitignore`
are left out. Afterwards run `python tools/run_tests.py 06` and read the new CHANGELOG before committing. Change sexloci in its own project, never here."""
import datetime, hashlib, os, shutil, sys, tempfile, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOP_FILES = ['README.md', 'CHANGELOG.md', 'LICENSE', 'pyproject.toml', 'environment.yml']
TOP_DIRS = ['config', 'docs', 'examples', 'scripts', 'tests']


def vendor(zip_path, date=None, here=HERE):
    date = date or datetime.date.today().isoformat()
    tmp = tempfile.mkdtemp()
    with zipfile.ZipFile(zip_path) as z:
        z.extractall(tmp)
    roots = [d for d in os.listdir(tmp) if os.path.isdir(os.path.join(tmp, d, 'src', 'sexloci'))]
    if len(roots) != 1:
        raise SystemExit('the zip must hold one folder with src/sexloci/ in it')
    src = os.path.join(tmp, roots[0])
    for d in ('sexloci', 'sexloci_dist'):
        shutil.rmtree(os.path.join(here, d), ignore_errors=True)
    shutil.copytree(os.path.join(src, 'src', 'sexloci'), os.path.join(here, 'sexloci'), ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    os.makedirs(os.path.join(here, 'sexloci_dist'))
    for f in TOP_FILES:
        shutil.copy(os.path.join(src, f), os.path.join(here, 'sexloci_dist', f))
    for d in TOP_DIRS:
        shutil.copytree(os.path.join(src, d), os.path.join(here, 'sexloci_dist', d), ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    version = next(l for l in open(os.path.join(here, 'sexloci', '__init__.py')) if '__version__' in l).split('=')[1].strip().strip('"\'')
    rows = []
    for top in ('sexloci', 'sexloci_dist'):
        for root, _, files in os.walk(os.path.join(here, top)):
            for f in files:
                p = os.path.join(root, f)
                rows.append((os.path.relpath(p, here), os.path.getsize(p), hashlib.md5(open(p, 'rb').read()).hexdigest()))
    rows.sort()
    zmd5 = hashlib.md5(open(zip_path, 'rb').read()).hexdigest()
    with open(os.path.join(here, 'SEXLOCI_MANIFEST.tsv'), 'w') as o:
        o.write(f'# sexloci {version}, copied UNCHANGED on {date} from {os.path.basename(zip_path)} ({os.path.getsize(zip_path)} bytes, md5 {zmd5}).\n')
        o.write('# sexloci/ is src/sexloci/ of the zip; sexloci_dist/ is the rest of it (README, CHANGELOG, LICENSE, docs, config, examples, scripts, tests, pyproject.toml, environment.yml).\n')
        o.write('# Left out: .github/ (its CI workflow) and .gitignore. A change to any listed file is a failure of tests/test_vendored_sexloci.py: change sexloci in its own project and copy it in again.\n')
        o.write('path\tbytes\tmd5\n')
        for r in rows:
            o.write('%s\t%d\t%s\n' % r)
    shutil.rmtree(tmp, ignore_errors=True)
    return version, len(rows)


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    date = sys.argv[sys.argv.index('--date') + 1] if '--date' in sys.argv else None
    if len(args) < 1 or not os.path.exists(args[0]):
        raise SystemExit(__doc__)
    v, n = vendor(args[0], date)
    print(f'vendored sexloci {v}: {n} files; now run `python tools/run_tests.py 06` and read sexloci_dist/CHANGELOG.md')
