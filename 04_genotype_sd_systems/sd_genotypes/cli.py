"""python -m sd_genotypes panel measure --bam B --sample S --out-json J [--modules all] | panel check | gsdf call|validate | bchrom prepare|call|all|status|clean|compare | sample_table inversions | selftest"""
import os
from cichlid_common import cli_base


def _panel(a):                                                  # needs only samtools: runs on any machine, in the lab pipeline, without a Project
    from sd_genotypes.panel import measure as pm
    if a.step == 'check':                                       # preflight: definition files and keys of every module; exit code 1 if anything is missing
        rows = pm.check(a.modules)
        for r in rows:
            print(f"[panel] {r['module']:8s} {r['status']:8s} {r['key']:12s} {r['message']}")
        if any(r['status'] != 'OK' for r in rows):
            raise SystemExit(1)
        return True
    if a.step != 'measure' or not (a.bam and a.sample and a.out_json):
        raise SystemExit('usage: python -m sd_genotypes panel measure --bam <sample.all.bam> --sample <ID> --out-json <ID>.panel.json [--modules all | gsdf,te,lg5inv,bchrom]  (default: gsdf)\n'
                         '       python -m sd_genotypes panel check [--modules all]   (preflight: every module has its definition files)')
    res = pm.measure(a.bam, a.sample, pm.resolve_modules(a.modules))
    pm.write_json(a.out_json, res)
    print(f'[panel] {a.sample}: {", ".join(res["modules"])} -> {a.out_json}')
    return True


def _gsdf_offline(a):
    if not a.out_dir:
        return False
    from sd_genotypes.gsdf import steps as gs
    if a.step == 'call' and a.measurements:
        gs.call_step(a.measurements, out_dir=a.out_dir)
    elif a.step == 'validate' and a.table:
        gs.validate_step(a.calls, a.table, out_dir=a.out_dir)
    else:
        raise SystemExit('usage: gsdf call --measurements <dir> --out-dir <dir>   |   gsdf validate --table <tsv> [--calls <tsv>] --out-dir <dir>')
    return True


def _gsdf(proj, a, cfg):
    from sd_genotypes.gsdf import steps as gs
    if a.step == 'call' and a.measurements:
        gs.call_step(a.measurements, proj=proj)
    elif a.step == 'validate':
        gs.validate_step(a.calls, a.table or (_ for _ in ()).throw(SystemExit('gsdf validate needs --table <sample sheet tsv>')), proj=proj)
    else:
        raise SystemExit('usage: gsdf call --measurements <dir> | gsdf validate --table <tsv> [--calls <tsv>]  (add --out-dir for plain files, no Dropbox)')


def _bchrom(proj, a, cfg):
    from sd_genotypes.bchrom import steps as bc
    only = a.only.split(',') if a.only else None
    if a.step in ('prepare', 'all'):
        bc.prepare_step(proj, only=only, keep_bams=a.keep_bams, runs_in_flight=a.runs_in_flight, threads=a.threads)
    if a.step in ('call', 'all'):
        bc.call_step(proj, clean_depth=a.clean_depth)
    elif a.step == 'status':
        bc.status_step(proj)
    elif a.step == 'clean':
        bc.clean_step(proj, force=a.force)
    elif a.step == 'compare':
        bc.compare_step(proj, a.clark_bed, a.paf, a.name_map)
    elif a.step not in ('prepare', 'all'):
        raise SystemExit('usage: python -m sd_genotypes bchrom prepare | call | all | status | clean | compare')


def _sample_table(proj, a, cfg):
    from sd_genotypes.sample_table import steps as st
    if a.step != 'inversions' or not a.blumer_dir:
        raise SystemExit('usage: python -m sd_genotypes sample_table inversions --blumer-dir <dir> [--write-sheet]')
    st.inversions_step(proj, a.blumer_dir, a.write_sheet, a.bc_genotypes, a.bc_legacy_mc)


def main(argv=None):
    return cli_base.main('python -m sd_genotypes', ['panel', 'gsdf', 'bchrom', 'sample_table'], {'gsdf': _gsdf, 'bchrom': _bchrom, 'sample_table': _sample_table},
                         offline={'panel': _panel, 'gsdf': _gsdf_offline}, tests_dir=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tests'), argv=argv)
