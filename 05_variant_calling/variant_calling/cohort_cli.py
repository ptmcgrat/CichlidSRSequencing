"""python -m variant_calling cohort list | members <name> | check <name>

    cohort list                       the cohorts of config/callsets.yaml: kind, genome, caller, groups
    cohort members <name>             the SampleIDs the definition selects, one per line (notes go to stderr)
    cohort check <name> [--json F]    the rules of its kind against the SampleDatabase and the AlignmentDatabase; exit 0 ok, 1 BLOCKED, 2 usage
    cohort build <name> [--run ...]   the plan (default) or the build of the VCF with bcftools: see variant_calling/build.py; exit 0, 1 BLOCKED, 3 preflight problems, 4 build error

Data: the lab sheets through the shared file manager, or offline from TSV exports (--sample-table, --alignment-table)."""
import argparse
import json
import sys

import pandas as pd


def build_parser():
    ap = argparse.ArgumentParser(prog='python -m variant_calling cohort', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('command', choices=['list', 'members', 'check', 'build'])
    ap.add_argument('name', nargs='?', help='the cohort (an entry of callsets.yaml)')
    ap.add_argument('--cohorts', default=None, help='a cohort file (default: config/callsets.yaml)')
    ap.add_argument('--rules', default=None, help='a rules file (default: config/cohort_rules.yaml)')
    ap.add_argument('--genome', default=None, help="override the cohort's genome (default: its own, else Mzebra_GT3_NCBI)")
    ap.add_argument('--sample-table', default=None, help='OFFLINE: a tab-separated export of the SampleDatabase')
    ap.add_argument('--alignment-table', default=None, help='OFFLINE: a tab-separated export of the genome\'s AlignmentDatabase')
    ap.add_argument('--json', default=None, help='check: also write the report as JSON to this file')
    b = ap.add_argument_group('build')
    b.add_argument('--run', action='store_true', help='build: really stage the BAMs, call and write (default: print the plan only)')
    b.add_argument('--region', default=None, help='build: call only CONTIG or CONTIG:START-END (a trial of the plumbing)')
    b.add_argument('--contigs', default=None, help='build: chromosomes | all | a comma-separated list of contigs (default: config/cohort_build.yaml)')
    b.add_argument('--max-samples', type=int, default=None, help='build: a TRIAL with the N smallest BAMs; written locally, never uploaded')
    b.add_argument('--workers', type=int, default=None, help='build: contigs called at the same time')
    b.add_argument('--keep-staged', action='store_true', help='build: keep the staged BAMs afterwards (a later build of the same cohort reuses them)')
    b.add_argument('--keep-work', action='store_true', help='build: keep the per-contig calls afterwards')
    b.add_argument('--no-upload', action='store_true', help='build: write locally only')
    b.add_argument('--force', action='store_true', help='build: replace outputs that already exist locally or on Dropbox')
    b.add_argument('--build-config', default=None, help='build: a parameter file (default: config/cohort_build.yaml)')
    b.add_argument('--local-master', default=None, help='build, offline/tests: the local master directory')
    b.add_argument('--cloud-master', default=None, help='build, offline/tests: the Dropbox master directory')
    b.add_argument('--fasta', default=None, help='build, offline/tests: the reference FASTA (its .fai is read)')
    return ap


def _read_tsv(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def load_tables(a, genome):
    """(SampleDatabase, AlignmentDatabase) as DataFrames: the TSV exports when given, else the lab sheets (read-only)."""
    if a.sample_table or a.alignment_table:
        if not (a.sample_table and a.alignment_table):
            raise SystemExit('offline mode needs both --sample-table and --alignment-table')
        return _read_tsv(a.sample_table), _read_tsv(a.alignment_table)
    from cichlid_common.file_manager import FileManager
    fm = FileManager()
    fm.readSampleDatabase()
    fm.setGenome(genome)
    fm.readAlignmentDatabase()
    return fm.sample_dt, fm.alignment_dt


def main(argv=None, out=None):
    out = out or sys.stdout
    a = build_parser().parse_args(argv)
    from variant_calling import cohort, rules
    cohorts = cohort.load_cohorts(a.cohorts)
    if a.command == 'list':
        for n, s in cohorts.items():
            print(f"{n:28s} {str(s.get('kind')):11s} {s.get('genome'):18s} caller={s.get('caller'):9s} groups={s.get('groups'):11s} {s.get('description', '')}", file=out)
        return 0
    if not a.name or a.name not in cohorts:
        print(f"unknown cohort {a.name!r}; the cohorts are: {', '.join(cohorts)}", file=sys.stderr)
        return 2
    spec = dict(cohorts[a.name])
    if a.genome:
        spec['genome'] = a.genome
    sample_dt, ledger = load_tables(a, spec['genome'])
    if a.command == 'members':
        problems = cohort.validate(spec)
        if problems:
            print('the definition is not usable: ' + '; '.join(problems), file=sys.stderr)
            return 2
        members, notes = cohort.resolve_members(sample_dt, spec)
        for n in notes:
            print('note: ' + n, file=sys.stderr)
        print('\n'.join(members), file=out)
        return 0
    if a.command == 'build':
        return _build(a, spec, sample_dt, ledger, out)
    report = rules.check_cohort(a.name, spec, sample_dt, ledger, rules.load_rules(a.rules))
    print(report.text(), file=out)
    if a.json:
        with open(a.json, 'w') as f:
            json.dump(report.to_dict(), f, indent=1)
    return 1 if report.blocked else 0


def _build(a, spec, sample_dt, ledger, out):
    from variant_calling import build as B
    cfg = B.load_build_config(a.build_config)
    if a.local_master or a.cloud_master or a.fasta:
        if not (a.local_master and a.cloud_master and a.fasta):
            raise SystemExit('--local-master, --cloud-master and --fasta go together')
        env = B.Env(a.local_master, a.cloud_master, spec['genome'], a.fasta)
    else:
        env = B.Env.from_file_manager(spec['genome'])
    try:
        res = B.build_cohort(a.name, spec, sample_dt, ledger, env, cfg=cfg, region=a.region, contigs=a.contigs, max_samples=a.max_samples, run_it=a.run, workers=a.workers,
                             keep_staged=a.keep_staged, keep_work=a.keep_work, upload=not a.no_upload, force=a.force, log=lambda m: print(m, file=out, flush=True))
    except B.BuildError as e:
        print('BUILD ERROR: ' + str(e), file=sys.stderr)
        return 4
    if res['status'] == 'blocked':
        print(res['report'].text(), file=out)
        return 1
    if res['status'] == 'plan':
        print(res['text'], file=out)
        return 3 if res['problems'] else 0
    return 0
