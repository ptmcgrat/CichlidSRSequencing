"""python -m variant_calling cohort list | members <name> | check <name>   (see cohort_cli)
   python -m variant_calling effcov measure ... | report ...                (effective coverage of BAMs: see effcov_cli)
   python -m variant_calling callset build <name> [--dry-run] | selftest"""
import os
import sys
from cichlid_common import cli_base


def _callset(proj, a, cfg):
    from variant_calling import callset
    if a.step != 'build' or not a.name:
        raise SystemExit('usage: python -m variant_calling callset build <name> [--dry-run]')
    callset.build(proj, a.name, dry_run=a.dry_run)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] == 'cohort':                       # its own small parser: the shared one is a single global parser
        from variant_calling import cohort_cli
        return cohort_cli.main(argv[1:])
    if argv and argv[0] == 'effcov':
        from variant_calling import effcov_cli
        return effcov_cli.main(argv[1:])
    return cli_base.main('python -m variant_calling', ['callset'], {'callset': _callset}, tests_dir=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tests'), argv=argv)
