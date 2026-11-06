import os, sys
from cichlid_common import cli_base


def _lg5inv(proj, a, cfg):
    from sd_genotypes.lg5inv import steps as S
    if a.step == 'prepare':
        S.prepare_step(proj, only=a.only.split(',') if a.only else None, keep_bams=a.keep_bams, runs_in_flight=a.runs_in_flight, threads=a.threads)
    elif a.step == 'status':
        S.status_step(proj)
    elif a.step == 'clean':
        S.clean_step(proj, force=a.force)
    else:
        raise SystemExit('usage: python -m sd_genotypes.lg5inv prepare [--only RUN,...] [--runs-in-flight N] [--threads N] | status | clean [--force]')


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] != 'selftest':
        argv.insert(0, 'lg5inv')                                                      # `python -m sd_genotypes.lg5inv prepare` = analysis lg5inv, step prepare
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    return cli_base.main('python -m sd_genotypes.lg5inv', ['lg5inv'], {'lg5inv': _lg5inv}, tests_dir=os.path.join(here, 'tests'), argv=argv)


if __name__ == '__main__':
    sys.exit(main())
