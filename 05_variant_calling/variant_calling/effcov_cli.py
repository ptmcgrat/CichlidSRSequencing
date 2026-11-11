"""python -m variant_calling effcov measure (--cohort NAME | --samples-file F | -s ID [ID ...] | --all [--exclude-ecogroup E]) [--run] [options]
   python -m variant_calling effcov report [--tsv F] [--by COLUMN [--summary-only]] [selection]
   python -m variant_calling effcov write [--apply] [--overwrite] [--limit N] [--tsv F] [selection]

    effcov measure ...   the plan (default): which samples, how much to download, the definition and any problem; with --run it downloads each BAM, measures the effective coverage (see
                         variant_calling/effcov.py), deletes the BAM and appends one row per sample to the output file (default <local master>/Outputs/EffectiveCoverage/effective_coverage.<genome>.tsv).
                         A rerun skips samples already measured with the same definition (--redo measures them again); a failed sample is recorded and does not stop the others.
    effcov write ...     put the measured EffectiveCoverage, MapqPassFraction and BaseSurvivalFraction into the AlignmentDatabase (see cichlid_common/sheet_writer.py): the default is a PLAN that writes
                         nothing; --apply adds the three columns at the END of the header if they are missing, fills EMPTY cells of the named samples only, re-reads the sheet and checks every cell, and keeps a
                         journal. A cell that already holds a different value is left alone unless --overwrite. --limit N writes only the first N samples (a first trial on the live sheet).
    effcov report ...    the table of what is in the file; --by COLUMN (a column of the SampleDatabase or the AlignmentDatabase: Species, RunIDs, BroodID, Ecogroup ...) adds one line per group with the medians of the
                         quality measures (duplicate rate and its range, MAPQ pass, bases kept, effective/table), --summary-only prints only those lines. Read-only.
                         The table of what is in the file: effective against table coverage, duplicates, MAPQ and base survival, the samples below the minimum depth of the sex analysis.

Exit codes: 0 ok, 1 some samples failed, 2 usage, 3 the plan has problems. Data: the lab sheets through the shared file manager (read-only), or offline from TSV exports (--sample-table, --alignment-table).
This command does not write to any Google sheet."""
import argparse
import os
import sys

DEFAULT_GENOME = 'Mzebra_GT3_NCBI'
TRUE_WORDS = ('1', '1.0', 'true', 'yes')          # as the cohort rules read LabReared


def build_parser():
    ap = argparse.ArgumentParser(prog='python -m variant_calling effcov', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('command', choices=['measure', 'report', 'write'])
    s = ap.add_argument_group('which samples (exactly one)')
    s.add_argument('--cohort', default=None, help='the members of a cohort of config/callsets.yaml (what `cohort members` prints)')
    s.add_argument('--samples-file', default=None, help='a text file of SampleIDs, one per line (# starts a comment)')
    s.add_argument('-s', '--samples', nargs='+', default=None, help='SampleIDs')
    s.add_argument('--all', action='store_true', help='every sample of the genome\'s AlignmentDatabase')
    s.add_argument('--lab-reared', action='store_true', help='every sample of the AlignmentDatabase that the SampleDatabase marks LabReared')
    s.add_argument('--exclude-ecogroup', action='append', default=[], help='with --all (or any selection): leave out samples of this SampleDatabase Ecogroup (repeatable, e.g. Mbuna)')
    ap.add_argument('--genome', default=None, help=f"default: the cohort's genome, else {DEFAULT_GENOME}")
    ap.add_argument('--run', action='store_true', help='measure: really download and measure (default: print the plan only)')
    ap.add_argument('--apply', action='store_true', help='write: really write to the sheet (default: print the plan only)')
    ap.add_argument('--overwrite', action='store_true', help='write: replace a different value already in a cell (default: leave it and say so)')
    ap.add_argument('--limit', type=int, default=None, help='write: only the first N samples, alphabetically')
    ap.add_argument('--by', default=None, help='report: group the samples by this column of the SampleDatabase or the AlignmentDatabase (Species, RunIDs, BroodID, Ecogroup, ...)')
    ap.add_argument('--summary-only', action='store_true', help='report with --by: only the group lines')
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


def _text(v):
    """A cell as text: the lab sheets give NaN (a float) for a blank cell, which is not an empty string and is truthy."""
    return v.strip() if isinstance(v, str) else ''


def sex_map(sample_dt):
    """{SampleID: Sex} for the samples that have both."""
    if sample_dt is None or 'Sex' not in getattr(sample_dt, 'columns', []):
        return {}
    return {_text(s): _text(x) for s, x in zip(sample_dt['SampleID'], sample_dt['Sex']) if _text(s) and _text(x)}


def group_map(column, sample_dt, ledger):
    """{SampleID: label} from `column` of the SampleDatabase, else of the AlignmentDatabase; None when neither has it."""
    for table in (sample_dt, ledger):
        if table is not None and column in getattr(table, 'columns', []) and 'SampleID' in table.columns:
            return {_text(s): _text(v) for s, v in zip(table['SampleID'], table[column]) if _text(s)}
    return None


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
    elif a.lab_reared:
        in_ledger = [_text(s) for s in ledger['SampleID'].tolist() if _text(s)] if ledger is not None and len(ledger) else []
        if sample_dt is None or 'LabReared' not in getattr(sample_dt, 'columns', []):
            notes.append('the SampleDatabase has no LabReared column')
            members = []
        else:
            flagged = {_text(s) for s, v in zip(sample_dt['SampleID'], sample_dt['LabReared']) if _text(s) and _text(v).lower() in TRUE_WORDS}
            members = [s for s in in_ledger if s in flagged]
            if len(flagged) > len(set(members)):
                notes.append(f'{len(flagged - set(members))} samples marked LabReared have no row in the AlignmentDatabase and are skipped')
    else:
        members = [_text(s) for s in ledger['SampleID'].tolist() if _text(s)] if ledger is not None and len(ledger) else []
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


def write_to_sheet(a, env, cfg, ecfg, rows, members, path, genome, out, err, worksheet_fn):
    """`effcov write`: the plan, and with --apply the guarded write into the AlignmentDatabase."""
    from variant_calling import build as B, effcov as E
    from cichlid_common import sheet_writer as SW
    if not (env.genome_fa and os.path.exists(env.genome_fa + '.fai')):
        print('the reference index is needed to tell which rows were measured with the current definition', file=err)
        return 2
    key = E.definition_key(cfg, ecfg, E.windows_from_fai(B.read_fai(env.genome_fa + '.fai'), ecfg['window_bp'], ecfg['position_fraction'], cfg['min_contig_bp']))
    ok = E.latest_ok(rows, key)
    if not ok:
        print(f'{path} has no successful measurement under the current definition ({key}): measure first (older rows are not written)', file=err)
        return 2
    updates = E.sheet_updates(ok, members)
    if not updates:
        print('no measured sample is selected', file=err)
        return 2
    if worksheet_fn is None:
        if a.local_master or a.cloud_master or a.fasta:
            print('writing needs the live sheet: do not combine it with --local-master, --cloud-master or --fasta', file=err)
            return 2
        from cichlid_common.file_manager import FileManager
        fm = FileManager()
        fm.setGenome(genome)
        worksheet_fn = lambda: fm.gc.open_by_key(fm.a_ID).worksheet(fm.genome_version)         # exactly how the alignment pipeline opens it
    cu = SW.ColumnUpdater(worksheet_fn, key='SampleID', log=lambda m: print(m, file=out))
    try:
        plan = cu.plan(updates, overwrite=a.overwrite)
    except SW.SheetWriteError as e:
        print(f'the AlignmentDatabase cannot be written: {e}', file=err)
        return 3
    print(f'AlignmentDatabase, worksheet {genome}: {len(updates)} samples measured with definition {key}', file=out)
    print(plan.format(), file=out)
    if a.limit is not None:
        first = sorted({c['sample'] for c in plan.writable})[:a.limit]
        print(f"--limit {a.limit}: this run writes only {len(first)} samples ({', '.join(first[:6])}{', ...' if len(first) > 6 else ''}): {sum(1 for c in plan.writable if c['sample'] in first)} cells", file=out)
    if not a.apply:
        print('\nnothing was written: add --apply to write the cells above', file=out)
        return 0
    if not plan.writable:
        print('\nnothing to write', file=out)
        return 0
    journal = os.path.join(os.path.dirname(path), 'sheet_writes', f"alignment.{genome}.{__import__('datetime').datetime.now().strftime('%Y%m%d-%H%M%S')}.tsv")
    try:
        rep = cu.apply(updates, overwrite=a.overwrite, limit=a.limit, journal=journal)
    except SW.SheetWriteError as e:
        print(f'the write stopped: {e}', file=err)
        return 3
    print(f"\nwrote {rep.written} cells" + (f" (columns added: {', '.join(rep.added_columns)})" if rep.added_columns else '') + f"; journal {journal}", file=out)
    if not rep.ok:
        print(f'\nVERIFICATION FAILED for {len(rep.failures)} cells: the sheet changed while writing, or a write did not take. Nothing was undone automatically; see the journal (old values) and:', file=out)
        for f in rep.failures[:10]:
            print(f"  {SW.a1(f['row'], f['col'])} {f['sample']} {f['column']}: {f['problem']}", file=out)
        return 1
    print('every written cell was read back and is in the right row', file=out)
    return 0


def main(argv=None, out=None, err=None, run=None, script=None, which=None, free_bytes=None, worksheet_fn=None):
    out, err = out or sys.stdout, err or sys.stderr
    a = build_parser().parse_args(argv)
    from variant_calling import build as B, cohort, cohort_cli, effcov as E
    chosen = [bool(a.cohort), bool(a.samples_file), bool(a.samples), bool(a.all), bool(a.lab_reared)]
    if (a.command == 'measure' and sum(chosen) != 1) or sum(chosen) > 1:
        print('choose exactly one of --cohort, --samples-file, -s, --all and --lab-reared', file=err)
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
    sex = sex_map(sample_dt)
    if a.command in ('report', 'write'):
        members = None
        if any(chosen):
            members = set(select_members(a, sample_dt, ledger, spec)[0])
        rows = E.read_rows(path)
        if not rows:
            print(f'{path} has no rows', file=err)
            return 2
        key = None
        if a.command == 'write':
            return write_to_sheet(a, env, cfg, ecfg, rows, members, path, genome, out, err, worksheet_fn)
        if env.genome_fa and os.path.exists(env.genome_fa + '.fai'):
            fai = B.read_fai(env.genome_fa + '.fai')
            key = E.definition_key(cfg, ecfg, E.windows_from_fai(fai, ecfg['window_bp'], ecfg['position_fraction'], cfg['min_contig_bp']))
            if not any(r.get('DefinitionKey') == key for r in rows):
                key = None                                   # nothing under the current definition: show what there is
        groups = None
        if a.by:
            groups = group_map(a.by, sample_dt, ledger)
            if groups is None:
                print(f"no column {a.by!r} in the SampleDatabase or the AlignmentDatabase", file=err)
                return 2
        print(E.report_text(rows, key=key, sex=sex, min_depth=min_depth, members=members, groups=groups, summary_only=a.summary_only, group_label=a.by or 'group'), file=out)
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
