"""python -m variant_calling effcov measure (--cohort NAME | --samples-file F | -s ID [ID ...] | --all [--exclude-ecogroup E]) [--run] [options]
   python -m variant_calling effcov report [--tsv F] [selection]

    effcov measure ...   the plan (default): which samples, how much to download, the definition and any problem; with --run it downloads each BAM, measures the effective coverage (see
                         variant_calling/effcov.py), deletes the BAM and appends one row per sample to the output file (default <local master>/Outputs/EffectiveCoverage/effective_coverage.<genome>.tsv).
                         A rerun skips samples already measured with the same definition (--redo measures them again); a failed sample is recorded and does not stop the others.
    effcov report ...    the table of what is in the file: effective against table coverage, duplicates, MAPQ and base survival, the samples below the minimum depth of the sex analysis.

Exit codes: 0 ok, 1 some samples failed, 2 usage, 3 the plan has problems. Data: the lab sheets through the shared file manager (read-only), or offline from TSV exports (--sample-table, --alignment-table).
This command does not write to any Google sheet."""
import argparse
import os
import sys

DEFAULT_GENOME = 'Mzebra_GT3_NCBI'


def build_parser():
    ap = argparse.ArgumentParser(prog='python -m variant_calling effcov', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('command', choices=['measure', 'report'])
    s = ap.add_argument_group('which samples (exactly one)')
    s.add_argument('--cohort', default=None, help='the members of a cohort of config/callsets.yaml (what `cohort members` prints)')
    s.add_argument('--samples-file', default=None, help='a text file of SampleIDs, one per line (# starts a comment)')
    s.add_argument('-s', '--samples', nargs='+', default=None, help='SampleIDs')
    s.add_argument('--all', action='store_true', help='every sample of the genome\'s AlignmentDatabase')
    s.add_argument('--exclude-ecogroup', action='append', default=[], help='with --all (or any selection): leave out samples of this SampleDatabase Ecogroup (repeatable, e.g. Mbuna)')
    ap.add_argument('--genome', default=None, help=f"default: the cohort's genome, else {DEFAULT_GENOME}")
    ap.add_argument('--run', action='store_true', help='measure: really download and measure (default: print the plan only)')
    ap.add_argument('--redo', action='store_true', help='measure again samples already measured with this definition')
    ap.add_argument('--jobs', type=int, default=None, help='samples measured at the same time (default: config/effcov.yaml)')
    ap.add_argument('--threads', type=int, default=None, help='samtools threads per sample (default: config/effcov.yaml)')
    ap.add_argument('--keep-bams', action='store_true', help='keep the downloaded BAMs (a rerun reuses them)')
    ap.add_argument('--out', default=None, help='the output file (default: <local master>/Outputs/EffectiveCoverage/effective_coverage.<genome>.tsv)')
    ap.add_argument('--tsv', default=None, help='report: the file to report (default: as --out)')
    ap.add_argument('--min-depth', type=float, default=None, help='report: the depth below which a sample is flagged (default: config/effcov.yaml, sexloci\'s min_dp)')
    ap.add_argument('--cohorts', default=None, help='a cohort file (default: config/callsets.yaml)')
    ap.add_argument('--build-config', default=None, help='the caller\'s parameters (default: config/cohort_build.yaml)')
    ap.add_argument('--config', default=None, help='this command\'s parameters (default: config/effcov.yaml)')
    o = ap.add_argument_group('offline / tests')
    o.add_argument('--sample-table', default=None, help='a tab-separated export of the SampleDatabase')
    o.add_argument('--alignment-table', default=None, help='a tab-separated export of the genome\'s AlignmentDatabase')
    o.add_argument('--local-master', default=None)
    o.add_argument('--cloud-master', default=None)
    o.add_argument('--fasta', default=None, help='the reference FASTA (its .fai is read)')
    return ap


def read_ids(path):
    out = []
    for line in open(path):
        line = line.split('#')[0].strip()
        if line:
            out.append(line)
    return out


def select_members(a, sample_dt, ledger, spec=None, notes=None):
    """The SampleIDs to measure, in the order given, without repeats."""
    from variant_calling import cohort
    notes = notes if notes is not None else []
    if a.cohort:
        members, more = cohort.resolve_members(sample_dt, spec)
        notes += more
    elif a.samples_file:
        members = read_ids(a.samples_file)
    elif a.samples:
        members = list(a.samples)
    else:
        members = [s for s in ledger['SampleID'].tolist() if str(s).strip()] if ledger is not None and len(ledger) else []
    if a.exclude_ecogroup and sample_dt is not None and 'Ecogroup' in sample_dt.columns:
        gone = set(sample_dt.loc[sample_dt['Ecogroup'].isin(a.exclude_ecogroup), 'SampleID'])
        left_out = [s for s in members if s in gone]
        if left_out:
            notes.append(f"{len(left_out)} samples of Ecogroup {', '.join(a.exclude_ecogroup)} are left out")
        members = [s for s in members if s not in gone]
    known = set(sample_dt['SampleID']) if sample_dt is not None and len(sample_dt) else set()
    unknown = [s for s in members if known and s not in known]
    if unknown:
        notes.append(f'{len(unknown)} samples are not in the SampleDatabase: ' + ', '.join(unknown[:5]))
    return list(dict.fromkeys(members)), notes


def main(argv=None, out=None, err=None, run=None, script=None, which=None, free_bytes=None):
    out, err = out or sys.stdout, err or sys.stderr
    a = build_parser().parse_args(argv)
    from variant_calling import build as B, cohort, cohort_cli, effcov as E
    chosen = [bool(a.cohort), bool(a.samples_file), bool(a.samples), bool(a.all)]
    if a.command == 'measure' and sum(chosen) != 1:
        print('choose exactly one of --cohort, --samples-file, -s and --all', file=err)
        return 2
    spec, genome = None, a.genome
    if a.cohort:
        cohorts = cohort.load_cohorts(a.cohorts)
        if a.cohort not in cohorts:
            print(f"unknown cohort {a.cohort!r}; the cohorts are: {', '.join(cohorts)}", file=err)
            return 2
        spec = dict(cohorts[a.cohort])
        problems = cohort.validate(spec)
        if problems:
            print('the definition is not usable: ' + '; '.join(problems), file=err)
            return 2
        genome = genome or spec.get('genome')
    genome = genome or DEFAULT_GENOME
    cfg, ecfg = B.load_build_config(a.build_config), E.load_config(a.config)
    min_depth = a.min_depth if a.min_depth is not None else ecfg['min_depth']
    if a.local_master or a.cloud_master or a.fasta:
        if not (a.local_master and a.cloud_master and a.fasta):
            print('--local-master, --cloud-master and --fasta go together', file=err)
            return 2
        env = B.Env(a.local_master, a.cloud_master, genome, a.fasta)
    else:
        env = B.Env.from_file_manager(genome)
    path = a.out or a.tsv or f"{env.local_master}{ecfg['output_subdir']}effective_coverage.{genome}.tsv"
    try:
        sample_dt, ledger = cohort_cli.load_tables(a, genome)
    except SystemExit as e:
        print(str(e), file=err)
        return 2
    sex = dict(zip(sample_dt['SampleID'], sample_dt['Sex'])) if sample_dt is not None and 'Sex' in getattr(sample_dt, 'columns', []) else {}
    if a.command == 'report':
        members = None
        if any(chosen):
            members = set(select_members(a, sample_dt, ledger, spec)[0])
        rows = E.read_rows(path)
        if not rows:
            print(f'{path} has no rows', file=err)
            return 2
        key = None
        if env.genome_fa and os.path.exists(env.genome_fa + '.fai'):
            fai = B.read_fai(env.genome_fa + '.fai')
            key = E.definition_key(cfg, ecfg, E.windows_from_fai(fai, ecfg['window_bp'], ecfg['position_fraction'], cfg['min_contig_bp']))
            if not any(r.get('DefinitionKey') == key for r in rows):
                key = None                                   # nothing under the current definition: show what there is
        print(E.report_text(rows, key=key, sex=sex, min_depth=min_depth, members=members), file=out)
        return 0
    members, notes = select_members(a, sample_dt, ledger, spec)
    for n in notes:
        print('note: ' + n, file=out)
    plan = E.make_plan(env, cfg, ecfg, members, ledger, path, redo=a.redo, jobs=a.jobs, threads=a.threads)
    kw = {k: v for k, v in dict(which=which, run=run).items() if v is not None}
    problems = E.preflight(env, plan, free_bytes=free_bytes, **kw)
    print(E.format_plan(plan, problems, cfg, running=a.run), file=out)
    if not a.run:
        return 3 if problems else 0
    if problems:
        return 3
    kw = {k: v for k, v in dict(run=run, script=script).items() if v is not None}
    res = E.run_effcov(env, cfg, plan, keep_bams=a.keep_bams, log=lambda m: print(m, file=out, flush=True), **kw)
    print('', file=out)
    print(E.report_text(E.read_rows(path), key=plan.key, sex=sex, min_depth=min_depth, members=set(plan.members)), file=out)
    if res['failed']:
        print(f"\n{len(res['failed'])} samples failed (see the Error column of {path}); run the same command again to retry them", file=out)
    return 1 if res['failed'] else 0
