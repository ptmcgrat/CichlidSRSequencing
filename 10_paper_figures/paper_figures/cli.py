"""python -m paper_figures list | build <figure> [script options]      each figure is a script in paper_figures/figures/ that reads versioned run folders"""
import glob, os, subprocess, sys

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures')


def figures():
    return sorted(os.path.basename(f)[:-3] for f in glob.glob(os.path.join(FIG, 'fig_*.py')))


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__); return 0
    if argv[0] == 'list':
        print('\n'.join(figures())); return 0
    if argv[0] == 'build' and len(argv) >= 2 and argv[1] in figures():
        return subprocess.run([sys.executable, os.path.join(FIG, argv[1] + '.py')] + argv[2:]).returncode
    raise SystemExit(f'usage: python -m paper_figures list | build <figure> [options]; figures: {figures()}')
