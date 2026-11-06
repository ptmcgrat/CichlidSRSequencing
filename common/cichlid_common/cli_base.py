"""Shared command-line machinery of the module CLIs (python -m landmarks | sd_genotypes | variant_calling | sd_regions | bc_analysis | paper_figures | cichlid_common).

Every module keeps the options of the original single CLI, so a command line that worked before only changes its first words:
    python -m yhsd.cli centromeres arrays   ->   python -m landmarks centromeres arrays
A module lists the analyses it owns, the ones that need no Project (offline), and a handler per analysis."""
import argparse, os, sys


def build_parser(prog, analyses):
    ap = argparse.ArgumentParser(prog=prog)
    ap.add_argument('analysis', choices=analyses)
    ap.add_argument('step', nargs='?')
    ap.add_argument('name', nargs='?', help='callset build <name>')
    ap.add_argument('--config', default=None, help='default: config/paper.yaml next to this package')
    ap.add_argument('--execute', action='store_true', help='stage-legacy: actually copy (default dry run)')
    ap.add_argument('--dry-run', action='store_true', help='callset build: list the samples and stop')
    ap.add_argument('--rho', type=float, default=None, help='v2-hmm: crossover rate per bp (sensitivity runs)')
    ap.add_argument('--tag', default='', help='v2-hmm: output suffix, e.g. _rho2e-8')
    ap.add_argument('--species', default=None, help='parental_sex_region: comma-separated, e.g. YH,MC')
    ap.add_argument('--test', action='store_true', help="write to the _test project folder/ (never uploaded) so a trial cannot contaminate the real run")
    ap.add_argument('--no-upload', action='store_true', help='skip the Dropbox upload and verification')
    ap.add_argument('--blumer-dir', default=None, help='sample_table inversions: folder with Table_S1.tsv, Table_S7.tsv, Data_S1.nwk, Data_S2.nwk')
    ap.add_argument('--write-sheet', action='store_true', help='sample_table inversions: append the new columns to the Google SampleDatabase (default: dry run)')
    ap.add_argument('--bc-genotypes', default=None, help="sample_table inversions: BC haplotype-call grid to infer chr2/chr13 of BC fish: a path (downloaded via the FileManager if missing), or 'legacy' = 2026_Nikesh/OutputGenotypes.csv")
    ap.add_argument('--bc-legacy-mc', action='store_true', help='the BC grid uses the old Mconophoros_GT1 chromosome names (NC_036781.1 = LG2, NC_036792.1 = LG13)')
    ap.add_argument('--table', default=None, help='sex_association fisher: run OFFLINE on this tab-separated sample table (no FileManager, no upload); centromeres lift: the published centromere workbook (xlsx)')
    ap.add_argument('--out', default='sex_association_out', help='sex_association fisher --table: output folder')
    ap.add_argument('--group-by', default=None, choices=['clade', 'species'], help='sex_association fisher: group by clade (default) or species')
    ap.add_argument('--only', default=None, help='genomes fetch: comma-separated registry keys, e.g. LatrZW,MezeOBm')
    ap.add_argument('--chroms', default=None, help='comma-separated LG names, e.g. LG10,LG3 (use with --test)')
    ap.add_argument('--samples', default=None, help='comma-separated BC sample IDs (use with --test)')
    ap.add_argument('--keep-bams', action='store_true', help='bchrom prepare: keep the duplicate-marked BAM of every run (default: delete everything large once the depth track exists)')
    ap.add_argument('--runs-in-flight', type=int, default=None, help='bchrom prepare: runs processed at the same time (config: resources.runs_in_flight)')
    ap.add_argument('--threads', type=int, default=None, help='bchrom prepare: minimap2 threads per run (config: resources.align_threads)')
    ap.add_argument('--clean-depth', action='store_true', help='bchrom call: delete the depth tracks and the work folder once the call run is saved')
    ap.add_argument('--force', action='store_true', help='bchrom clean: delete even if no `bchrom call` run exists')
    ap.add_argument('--clark-bed', default=None, help='bchrom compare: one of the Core_blocks/*.bed files of Clark et al. 2018 (Appendix A)')
    ap.add_argument('--paf', default=None, help='bchrom compare: PAF of the Clark assembly aligned to the reference (default: align the registered Mzebra_UMD2a)')
    ap.add_argument('--name-map', default=None, help='bchrom compare: tab-separated Clark sequence name -> name in the aligned FASTA')
    ap.add_argument('--bam', default=None, help='panel measure: a coordinate-sorted, indexed BAM aligned to the project reference')
    ap.add_argument('--sample', default=None, help='panel measure: sample id stored in the JSON')
    ap.add_argument('--out-json', default=None, help='panel measure: output file (<ID>.panel.json)')
    ap.add_argument('--modules', default='gsdf', help='panel measure: comma-separated modules (default gsdf)')
    ap.add_argument('--measurements', default=None, help='gsdf call: a folder of *.panel.json (or one file)')
    ap.add_argument('--calls', default=None, help='gsdf validate: the gsdf_calls.tsv of `gsdf call` (omit to look at the existing `gsdf number` values only)')
    ap.add_argument('--out-dir', default=None, help='gsdf call / validate: write plain files here instead of an immutable run folder on Dropbox (also needs no Project)')
    return ap


def open_project(a):
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
    return Project(cfg), cfg


def run_module_tests(tests_dir):
    """Run a module's tests (unittest discovery in <module>/tests) with every module on the path."""
    import unittest
    from .common_light import module_dirs
    for d in module_dirs() + [tests_dir]:
        if d not in sys.path:
            sys.path.insert(0, d)
    suite = unittest.defaultTestLoader.discover(tests_dir, top_level_dir=tests_dir)
    res = unittest.TextTestRunner(verbosity=1).run(suite)
    return 0 if res.wasSuccessful() else 1


def main(prog, analyses, handlers, offline=None, tests_dir=None, argv=None):
    """handlers[analysis](proj, a, cfg); offline[analysis](a) -> True if it handled the command without a Project."""
    a = build_parser(prog, analyses + ['selftest']).parse_args(argv)
    if a.analysis == 'selftest':
        return run_module_tests(tests_dir)
    off = (offline or {}).get(a.analysis)
    if off and off(a):
        return 0
    proj, cfg = open_project(a)
    handlers[a.analysis](proj, a, cfg)
    return 0
