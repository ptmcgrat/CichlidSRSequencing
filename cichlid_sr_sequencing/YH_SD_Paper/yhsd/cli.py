"""python -m yhsd.cli <analysis> <step> [options]      run from cichlid_sr_sequencing/YH_SD_Paper/   (yhsd.cli needs no double-underscore files)

  python -m yhsd.cli selftest
  python -m yhsd.cli chroms
  python -m yhsd.cli callset build <name> [--dry-run]            # bc_family | yh_broods | mc_broods | ... (config/callsets.yaml)
  python -m yhsd.cli bc_genotyping <step>                        # provenance flagstat original-vcfs original v2-sites v2-counts v2-hmm validate compare
  python -m yhsd.cli parental_sex_region <step> [--species YH,MC]  # callset depth scan summarize
  python -m yhsd.cli sample_table inversions --blumer-dir <dir> [--bc-genotypes legacy|<csv> [--bc-legacy-mc]] [--write-sheet]   # add Blumer 2025 inversion genotypes to the SampleDatabase
Every step records a manifest and uploads a verified run folder to Dropbox (disable with --no-upload; --test never uploads)."""
import argparse, glob, os

BC_STEPS = ['provenance', 'flagstat', 'original-vcfs', 'original', 'v2-sites', 'v2-counts', 'v2-hmm', 'validate', 'compare', 'stage-legacy']
PSR_STEPS = ['callset', 'depth', 'scan', 'summarize']

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
    from .runner import publish
    pats, tar = BC_PUBLISH.get(step, ([], ()))
    paths = []
    for p in pats:
        paths += glob.glob(proj.root + p.format(tag=tag))
    if paths:
        publish(proj, 'bc_genotyping', step + tag, paths, tar=tar)


def main():
    ap = argparse.ArgumentParser(prog='python -m yhsd.cli')
    ap.add_argument('analysis', choices=['selftest', 'chroms', 'callset', 'bc_genotyping', 'parental_sex_region', 'sample_table'])
    ap.add_argument('step', nargs='?')
    ap.add_argument('name', nargs='?', help='callset build <name>')
    ap.add_argument('--config', default=None, help='default: config/paper.yaml next to this package')
    ap.add_argument('--execute', action='store_true', help='stage-legacy: actually copy (default dry run)')
    ap.add_argument('--dry-run', action='store_true', help='callset build: list the samples and stop')
    ap.add_argument('--rho', type=float, default=None, help='v2-hmm: crossover rate per bp (sensitivity runs)')
    ap.add_argument('--tag', default='', help='v2-hmm: output suffix, e.g. _rho2e-8')
    ap.add_argument('--species', default=None, help='parental_sex_region: comma-separated, e.g. YH,MC')
    ap.add_argument('--test', action='store_true', help="write to YH_SD_Paper_test/ (never uploaded) so a trial cannot contaminate the real run")
    ap.add_argument('--no-upload', action='store_true', help='skip the Dropbox upload and verification')
    ap.add_argument('--blumer-dir', default=None, help='sample_table inversions: folder with Table_S1.tsv, Table_S7.tsv, Data_S1.nwk, Data_S2.nwk')
    ap.add_argument('--write-sheet', action='store_true', help='sample_table inversions: append the new columns to the Google SampleDatabase (default: dry run)')
    ap.add_argument('--bc-genotypes', default=None, help="sample_table inversions: BC haplotype-call grid to infer chr2/chr13 of BC fish: a path (downloaded via the FileManager if missing), or 'legacy' = 2026_Nikesh/OutputGenotypes.csv")
    ap.add_argument('--bc-legacy-mc', action='store_true', help='the BC grid uses the old Mconophoros_GT1 chromosome names (NC_036781.1 = LG2, NC_036792.1 = LG13)')
    ap.add_argument('--chroms', default=None, help='comma-separated LG names, e.g. LG10,LG3 (use with --test)')
    ap.add_argument('--samples', default=None, help='comma-separated BC sample IDs (use with --test)')
    a = ap.parse_args()
    if a.analysis == 'selftest':
        from . import selftest
        return selftest.main()
    from .common import load_config, Project
    cfg = load_config(a.config)
    if a.test:
        cfg['project_dir_name'] += '_test'
    if a.chroms:
        cfg['chrom_subset'] = a.chroms.split(',')
    if a.samples:
        cfg['sample_subset'] = a.samples.split(',')
    if (a.chroms or a.samples) and not a.test:
        raise SystemExit('--chroms/--samples produce partial outputs; add --test so they go to a separate directory')
    cfg.setdefault('runner', {})['upload'] = not (a.test or a.no_upload)
    proj = Project(cfg)
    species = a.species.split(',') if a.species else None

    if a.analysis == 'chroms':
        print(proj.chrom_table.to_string(index=False))
    elif a.analysis == 'callset':
        from . import callset
        if a.step != 'build' or not a.name:
            raise SystemExit('usage: python -m yhsd.cli callset build <name>')
        callset.build(proj, a.name, dry_run=a.dry_run)
    elif a.analysis == 'bc_genotyping':
        from .analyses.bc_genotyping import stage, provenance, flagstat, original, v2, validate, compare
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
    elif a.analysis == 'sample_table':
        from .analyses.sample_table import steps as st
        if a.step != 'inversions' or not a.blumer_dir:
            raise SystemExit('usage: python -m yhsd.cli sample_table inversions --blumer-dir <dir> [--write-sheet]')
        st.inversions_step(proj, a.blumer_dir, a.write_sheet, a.bc_genotypes, a.bc_legacy_mc)
    elif a.analysis == 'parental_sex_region':
        from .analyses.parental_sex_region import steps
        fn = {'callset': steps.callset_step, 'depth': steps.depth_step, 'scan': steps.scan_step, 'summarize': steps.summarize_step}
        if a.step not in fn:
            raise SystemExit(f'steps: {PSR_STEPS}')
        fn[a.step](proj, species)


if __name__ == '__main__':
    main()
