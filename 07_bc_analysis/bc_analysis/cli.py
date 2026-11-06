"""python -m bc_analysis bc_genotyping provenance|flagstat|original-vcfs|original|v2-sites|v2-counts|v2-hmm|validate|compare|stage-legacy | selftest"""
import os
from cichlid_common import cli_base

BC_STEPS = ['provenance', 'flagstat', 'original-vcfs', 'original', 'v2-sites', 'v2-counts', 'v2-hmm', 'validate', 'compare', 'stage-legacy']

# what each BC step publishes to Dropbox (relative to the genome project folder); `tar` = directories of many small files
BC_PUBLISH = {
    'provenance': (['qc/provenance.txt', 'qc/chromosomes.tsv'], ()),
    'flagstat': (['qc/flagstat_and_mapq.tsv', 'qc/flagstat_summary.tsv', 'qc/mapq20_fraction_by_chrom.tsv'], ()),
    'original-vcfs': (['original/AllBroodParents_Unfiltered.vcf.gz*', 'original/AllBroodParents_Filtered.vcf.gz*', 'original/BC_Pileup.vcf.gz*',
                       'original/BC_Variants.vcf.gz*', 'original/BC_Pileup_Samples.txt', 'original/filter_expression.txt'], ()),
    'original': (['original/OutputGenotypes_original.csv', 'original/markers_per_brood_chrom.tsv', 'original/per_fish_chrom_reads_switches.tsv'], ()),
    'v2-sites': (['v2/sites', 'v2/bc_family.xcheck.vcf.gz*', 'qc/v2_site_attrition.tsv', 'qc/v2_sites_per_mb.tsv'], ()),
    'v2-counts': (['v2/counts'], ('counts',)),
    'v2-hmm': (['v2/OutputPosteriors{tag}.csv', 'v2/breakpoints{tag}.tsv', 'v2/sample_qc{tag}.tsv'], ()),
    'validate': (['qc/validation_downsample.tsv'], ()),
    'compare': (['comparison'], ()),
}


def publish_bc(proj, step, tag=''):
    from cichlid_common.runner import publish
    pats, tar = BC_PUBLISH.get(step, ([], ()))
    paths = []
    for p in pats:
        paths += glob.glob(proj.root + p.format(tag=tag))
    if paths:
        publish(proj, 'bc_genotyping', step + tag, paths, tar=tar)




def _bc(proj, a, cfg):
    from bc_analysis.bc_genotyping import stage, provenance, flagstat, original, v2, validate, compare
    fn = {'provenance': lambda: provenance.main(proj), 'flagstat': lambda: flagstat.main(proj),
          'original-vcfs': lambda: original.make_vcfs(proj), 'original': lambda: original.main(proj),
          'v2-sites': lambda: v2.define_sites(proj), 'v2-counts': lambda: v2.count_reads(proj),
          'v2-hmm': lambda: v2.run_hmm(proj, a.rho, a.tag), 'validate': lambda: validate.main(proj),
          'compare': lambda: compare.main(proj), 'stage-legacy': lambda: stage.main(proj, a.execute)}
    if a.step not in fn:
        raise SystemExit(f'steps: {BC_STEPS}')
    fn[a.step]()
    if cfg['runner']['upload']:
        publish_bc(proj, a.step, a.tag if a.step == 'v2-hmm' else '')


def main(argv=None):
    return cli_base.main('python -m bc_analysis', ['bc_genotyping'], {'bc_genotyping': _bc}, tests_dir=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tests'), argv=argv)
