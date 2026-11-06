"""python -m download_reads <command> ...

  run        download runs into uBAMs (the former downloadShortReadData_v2.py): --manifest TAB ..., Run Info files, --ena, --samples-file, --add-pending, --source ena|sra|auto
  manifest   the DownloadDatabase sheet (the former downloadManifest.py): import | status | reconcile | release | tabs
  bench      measure how many parallel downloads help on this network (the former ena_bench.py)
  patch-file-manager   add the download-database methods to the lab's file_manager.py (diff preview first)
  selftest   run this module's tests"""
import os, sys


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__); return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == 'run':
        from . import run; return run.main(rest)
    if cmd == 'manifest':
        from . import manifest_tool; return manifest_tool.main(rest)
    if cmd == 'bench':
        from . import bench; return bench.main(rest)
    if cmd == 'patch-file-manager':
        from .tools import patch_file_manager as p; return p.main(rest)
    if cmd == 'selftest':
        from cichlid_common.cli_base import run_module_tests
        return run_module_tests(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tests'))
    raise SystemExit(f'unknown command {cmd!r}\n' + __doc__)
