"""python -m variant_calling cohort list | members <name> | check <name>

    cohort list                       the cohorts of config/callsets.yaml: kind, genome, caller, groups
    cohort members <name>             the SampleIDs the definition selects, one per line (notes go to stderr)
    cohort check <name> [--json F]    the rules of its kind against the SampleDatabase and the AlignmentDatabase; exit 0 ok, 1 BLOCKED, 2 usage

Data: the lab sheets through the shared file manager, or offline from TSV exports (--sample-table, --alignment-table)."""
import argparse
import json
import sys

import pandas as pd


def build_parser():
    ap = argparse.ArgumentParser(prog='python -m variant_calling cohort', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('command', choices=['list', 'members', 'check'])
    ap.add_argument('name', nargs='?', help='the cohort (an entry of callsets.yaml)')
    ap.add_argument('--cohorts', default=None, help='a cohort file (default: config/callsets.yaml)')
    ap.add_argument('--rules', default=None, help='a rules file (default: config/cohort_rules.yaml)')
    ap.add_argument('--genome', default=None, help="override the cohort's genome (default: its own, else Mzebra_GT3_NCBI)")
    ap.add_argument('--sample-table', default=None, help='OFFLINE: a tab-separated export of the SampleDatabase')
    ap.add_argument('--alignment-table', default=None, help='OFFLINE: a tab-separated export of the genome\'s AlignmentDatabase')
    ap.add_argument('--json', default=None, help='check: also write the report as JSON to this file')
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
    report = rules.check_cohort(a.name, spec, sample_dt, ledger, rules.load_rules(a.rules))
    print(report.text(), file=out)
    if a.json:
        with open(a.json, 'w') as f:
            json.dump(report.to_dict(), f, indent=1)
    return 1 if report.blocked else 0
