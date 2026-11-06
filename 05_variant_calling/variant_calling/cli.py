"""python -m variant_calling callset build <name> [--dry-run] | selftest"""
import os
from cichlid_common import cli_base


def _callset(proj, a, cfg):
    from variant_calling import callset
    if a.step != 'build' or not a.name:
        raise SystemExit('usage: python -m variant_calling callset build <name> [--dry-run]')
    callset.build(proj, a.name, dry_run=a.dry_run)


def main(argv=None):
    return cli_base.main('python -m variant_calling', ['callset'], {'callset': _callset}, tests_dir=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tests'), argv=argv)
