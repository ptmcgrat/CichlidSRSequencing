"""python -m sd_regions region_block --vcf V --qc sample_qc.tsv --region CHROM:START-END [--core S-E]   (is a block of sex-linked signal real? ROH, relatedness, PCA, LD, diagnostic SNPs)
   python -m sd_regions sexloci <command> ...   (sexloci 0.1.26, vendored unchanged: list-vcfs | run-folder | run | check | describe | filter-audit | filter-compare | simulate)
   python -m sd_regions parental_sex_region callset|depth|scan|summarize [--species YH,MC] | sex_association fisher [--table T --out D] | selftest"""
import os
import sys
from cichlid_common import cli_base

PSR_STEPS = ['callset', 'depth', 'scan', 'summarize']


def _sa_offline(a):                                             # offline mode: no Google/Dropbox access needed
    if not a.table:
        return False
    from sd_regions.sex_association import steps as sa
    sa.fisher_offline(a.table, a.out, group_by=a.group_by)
    return True


def _sa(proj, a, cfg):
    from sd_regions.sex_association import steps as sa
    if a.step != 'fisher':
        raise SystemExit('usage: python -m sd_regions sex_association fisher [--table <tsv> --out <dir>] [--group-by clade|species]')
    sa.fisher_step(proj, a.group_by)


def _psr(proj, a, cfg):
    from sd_regions.parental_sex_region import steps
    fn = {'callset': steps.callset_step, 'depth': steps.depth_step, 'scan': steps.scan_step, 'summarize': steps.summarize_step}
    if a.step not in fn:
        raise SystemExit(f'steps: {PSR_STEPS}')
    fn[a.step](proj, a.species.split(',') if a.species else None)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] == 'region_block':                  # its own parser as well: python -m sd_regions region_block --vcf ... --region CHROM:START-END
        from sd_regions.region_block.cli import main as rb_main
        return rb_main(argv[1:])
    if argv and argv[0] == 'sexloci':                      # the vendored package has its own parser: the shared one is a single global parser
        from sexloci.cli import main as sexloci_main
        return sexloci_main(argv[1:])
    return cli_base.main('python -m sd_regions', ['parental_sex_region', 'sex_association'], {'parental_sex_region': _psr, 'sex_association': _sa}, offline={'sex_association': _sa_offline}, tests_dir=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tests'), argv=argv)
