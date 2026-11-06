"""python -m landmarks inversion_breakpoints align|call|figure|all | centromeres scan|arrays|lift --table <xlsx> | selftest"""
import os
from cichlid_common import cli_base


def _ib(proj, a, cfg):
    from landmarks.inversion_breakpoints import steps as ib
    if a.step not in ('align', 'call', 'figure', 'all'):
        raise SystemExit('usage: python -m landmarks inversion_breakpoints align|call|figure|all')
    if a.step in ('align', 'all'): ib.align_step(proj)
    if a.step in ('call', 'all'): ib.call_step(proj)
    if a.step in ('figure', 'all'): ib.figure_step(proj)


def _cen(proj, a, cfg):
    from landmarks.centromeres import steps as cs_
    if a.step == 'scan':
        cs_.scan_step(proj)
    elif a.step == 'lift':
        cs_.lift_step(proj, a.table)
    elif a.step == 'arrays':
        cs_.arrays_step(proj)
    else:
        raise SystemExit('usage: python -m landmarks centromeres scan | centromeres arrays | centromeres lift --table <xlsx>')


def main(argv=None):
    return cli_base.main('python -m landmarks', ['inversion_breakpoints', 'centromeres'], {'inversion_breakpoints': _ib, 'centromeres': _cen}, tests_dir=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tests'), argv=argv)
