#!/usr/bin/env python3
"""Run the tests of every module (each in its own process, all modules on PYTHONPATH).

    python tools/run_tests.py                 all modules
    python tools/run_tests.py 02 04           only these modules (by number)
    python tools/run_tests.py -v              show the output of every module
"""
import glob, os, subprocess, sys, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main(argv):
    verbose = '-v' in argv
    only = [a for a in argv if a != '-v']
    dirs = [os.path.join(REPO, 'common')] + sorted(d for d in glob.glob(os.path.join(REPO, '[0-9][0-9]_*')) if os.path.isdir(d))
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(dirs + ([os.environ['PYTHONPATH']] if os.environ.get('PYTHONPATH') else [])), PYTHONDONTWRITEBYTECODE='1')
    rows, bad = [], 0
    for d in dirs:
        name = os.path.basename(d)
        t = os.path.join(d, 'tests')
        if not os.path.isdir(t) or (only and name[:2] not in only and name not in only):
            continue
        t0 = time.time()
        r = subprocess.run([sys.executable, '-W', 'ignore', '-m', 'unittest', 'discover', '-s', t, '-t', t], capture_output=True, text=True, cwd=d, env=env)
        ran = next((l for l in (r.stderr.strip().splitlines() or ['?']) if l.startswith('Ran ')), '?')
        bad += r.returncode != 0
        rows.append((name, ran, 'OK' if r.returncode == 0 else 'FAILED', round(time.time() - t0, 1)))
        if verbose or r.returncode != 0:
            print(f'---- {name}\n' + r.stderr[-3500:] + (r.stdout[-1500:] if r.returncode else ''))
    print(f'{"module":32s} {"tests":18s} {"result":8s} seconds')
    for n, ran, st, s in rows:
        print(f'{n:32s} {ran:18s} {st:8s} {s}')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
