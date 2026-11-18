"""python -m variant_calling panel measure (--samples-file F | -s S [S ...] | --ecogroup E) [--run] [--limit N] [--jobs N] [--redo] [--no-upload] [--keep-bams] [--panel-python P] [--panel-pythonpath D]
   python -m variant_calling panel write   [selection] [--tab T] [--json-dir D] [--apply] [--limit N]
   python -m variant_calling panel table   [selection] [--json-dir D] [--out F]

  measure   module 04's panel over existing BAMs: download, `panel measure`, JSON kept (locally and on Dropbox next to the BAM), BAM deleted; a plan unless --run. A stopped run resumes.
  write     the calls of the JSONs into the SampleDatabase tab (default SampleDatabase; --tab names a copy for trials) through the guarded writer: empty cells only, conflicts reported and never
            overwritten, a plan unless --apply, --limit N writes the first N samples only. Writes gsdf number, chr7ins, TE1, TE2 and, unless sheet.flags_in_sheet is false, their flags columns; sheet.nocall puts './.' in a genotype cell whose call has too little data (config/panelpass.yaml).
  table     every call (all six columns) and the numbers behind them in one TSV, with a summary: for reading the distributions.
Without a selection, `write` and `table` take every JSON of the folder. See panelpass.py."""
import argparse
import os
import sys

DEFAULT_GENOME = 'Mzebra_GT3_NCBI'


def build_parser():
    ap = argparse.ArgumentParser(prog='python -m variant_calling panel', description='module 04 panel over existing BAMs and its calls into the SampleDatabase')
    ap.add_argument('command', choices=['measure', 'write', 'table'])
    s = ap.add_argument_group('which samples (measure: exactly one; write and table: none = every JSON of the folder)')
    s.add_argument('--samples-file', default=None, help='a text file of SampleIDs, one per line (# starts a comment)')
    s.add_argument('-s', '--samples', nargs='+', default=None, help='SampleIDs')
    s.add_argument('--ecogroup', default=None, help='every sample of the SampleDatabase with this Ecogroup (e.g. Mbuna) that has a BAM in the AlignmentDatabase of the genome')
    ap.add_argument('--genome', default=DEFAULT_GENOME, help=f'the genome the BAMs are aligned to (default {DEFAULT_GENOME}: the panel is defined on it)')
    ap.add_argument('--run', action='store_true', help='measure: really download and measure (default: print the plan only)')
    ap.add_argument('--apply', action='store_true', help='write: really write to the sheet (default: print the plan only)')
    ap.add_argument('--limit', type=int, default=None, help='measure: only the first N samples still to do; write: only the first N samples')
    ap.add_argument('--jobs', type=int, default=None, help='measure: BAMs on disk at a time (default: config/panelpass.yaml)')
    ap.add_argument('--redo', action='store_true', help='measure: measure again samples already measured')
    ap.add_argument('--no-upload', action='store_true', help='measure: do not copy the JSONs to Dropbox')
    ap.add_argument('--keep-bams', action='store_true', help='measure: keep the downloaded BAMs (a rerun reuses them)')
    ap.add_argument('--panel-python', default=None, help='measure: the interpreter that can import sd_genotypes (default: the one running this command)')
    ap.add_argument('--panel-pythonpath', default=None, help='measure: where sd_genotypes is checked out, only if `source tools/env.sh` was not run')
    ap.add_argument('--json-dir', default=None, help='the folder of the <sample>.panel.json files (default: <local master>/Outputs/PanelJSON/<genome>/)')
    ap.add_argument('--tab', default=None, help='write: the worksheet of the SampleDatabase spreadsheet (default: config sheet.tab)')
    ap.add_argument('--out', default=None, help='table: the output file (default: <json dir>/panel_calls.tsv)')
    ap.add_argument('--config', default=None, help='default: config/panelpass.yaml')
    o = ap.add_argument_group('offline / tests')
    o.add_argument('--sample-table', default=None, help='a tab-separated export of the SampleDatabase')
    o.add_argument('--alignment-table', default=None, help='a tab-separated export of the genome\'s AlignmentDatabase')
    o.add_argument('--local-master', default=None)
    o.add_argument('--cloud-master', default=None)
    return ap


def select(a, sample_dt, ledger, genome, notes):
    """(the SampleIDs, or None when no selection was made); a selected sample without a BAM in the ledger is left out and said so."""
    from variant_calling import effcov_cli
    if a.samples_file:
        ids = effcov_cli.read_ids(a.samples_file)
    elif a.samples:
        ids = list(a.samples)
    elif a.ecogroup:
        if sample_dt is None or 'Ecogroup' not in sample_dt.columns:
            raise SystemExit('--ecogroup needs the SampleDatabase with an Ecogroup column')
        ids = sorted(str(x) for x in sample_dt.loc[sample_dt['Ecogroup'].astype(str) == a.ecogroup, 'SampleID'])
        notes.append(f'{len(ids)} samples of Ecogroup {a.ecogroup} in the SampleDatabase')
    else:
        return None
    ids = list(dict.fromkeys(ids))
    if ledger is not None and len(ledger) and a.command == 'measure':
        have = {str(x) for x in ledger['SampleID']}
        missing = [s for s in ids if s not in have]
        if missing:
            notes.append(f"{len(missing)} selected samples have no row in the AlignmentDatabase of {genome} (no BAM there): left out ({', '.join(missing[:6])}{', ...' if len(missing) > 6 else ''})")
        ids = [s for s in ids if s in have]
    return ids


def main(argv=None, out=None, err=None, run=None, runner=None, which=None, free_bytes=None, worksheet_fn=None, calls_module=None, fm_factory=None):
    out, err = out or sys.stdout, err or sys.stderr
    a = build_parser().parse_args(argv)
    from variant_calling import build as B, cohort_cli, panelpass as P
    chosen = [bool(a.samples_file), bool(a.samples), bool(a.ecogroup)]
    if sum(chosen) > 1 or (a.command == 'measure' and sum(chosen) != 1):
        print('choose exactly one of --samples-file, -s and --ecogroup' + ('' if a.command == 'measure' else ' (or none: every JSON of the folder)'), file=err)
        return 2
    genome = a.genome
    try:
        pcfg = P.load_config(a.config)
    except ValueError as e:
        print(str(e), file=err)
        return 2
    env = None
    if a.local_master or a.cloud_master:
        if not (a.local_master and a.cloud_master):
            print('--local-master and --cloud-master go together', file=err)
            return 2
        env = B.Env(a.local_master, a.cloud_master, genome, '')
    elif a.command == 'measure' or not a.json_dir:
        env = B.Env.from_file_manager(genome)
    json_dir = (a.json_dir.rstrip('/') + '/') if a.json_dir else f"{env.local_master}{pcfg['json_subdir']}{genome}/"
    sample_dt = ledger = None
    notes = []
    if a.command == 'measure' or a.ecogroup or a.command == 'table':
        try:
            sample_dt, ledger = cohort_cli.load_tables(a, genome)
        except (SystemExit, Exception) as e:                          # noqa: the sheets cannot be read (offline, no credentials, no gspread ...)
            if a.command == 'table' and not a.ecogroup:
                sample_dt = ledger = None                              # the table is useful without the species columns
                print(f'note: the SampleDatabase could not be read, so the table has no Species or Ecogroup column ({type(e).__name__}: {str(e)[:120]})', file=out)
            else:
                print(str(e), file=err)
                return 2
    try:
        members = select(a, sample_dt, ledger, genome, notes)
    except SystemExit as e:
        print(str(e), file=err)
        return 2
    for n in notes:
        print('note: ' + n, file=out)
    kw_run = {k: v for k, v in dict(runner=runner).items() if v is not None}

    if a.command == 'measure':
        try:
            panel = P.make_panel_config(pcfg, a.panel_python, a.panel_pythonpath)
            ok, text, keys = P.check_panel(panel, **kw_run)
        except P.PanelPassError as e:
            print(str(e), file=err)
            return 2
        plan = P.make_plan(env, pcfg, members, ledger, keys, redo=a.redo, jobs=a.jobs, limit=a.limit, json_dir=json_dir)
        kw = {k: v for k, v in dict(which=which, free_bytes=free_bytes).items() if v is not None}
        problems = P.preflight(env, pcfg, plan, ok, text, upload=not a.no_upload, **kw)
        print(P.format_plan(plan, problems, pcfg, running=a.run, upload=not a.no_upload), file=out)
        if not a.run:
            return 3 if problems else 0
        if problems:
            return 3
        kw = {k: v for k, v in dict(run=run, runner=runner, free_bytes=free_bytes).items() if v is not None}
        res = P.run_pass(env, pcfg, plan, panel, upload=not a.no_upload, keep_bams=a.keep_bams, log=lambda m: print(m, file=out, flush=True), **kw)
        if res['failed']:
            print(f"\n{len(res['failed'])} samples failed (see the Error column of {plan.log_path}); run the same command again to retry them", file=out)
        return 1 if res['failed'] else 0

    on_disk = P.list_json_samples(json_dir)
    samples = [s for s in (members if members is not None else on_disk)]
    if not on_disk:
        print(f'{json_dir} has no *.panel.json files', file=err)
        return 2
    try:
        calls_mod = calls_module or P.load_calls_module()
    except P.PanelPassError as e:
        print(str(e), file=err)
        return 2

    if a.command == 'table':
        info = {}
        if sample_dt is not None:
            keep = [c for c in ('Species', 'Ecogroup', 'Sex') if c in sample_dt.columns]
            info = {str(r['SampleID']): {c: r[c] for c in keep} for _, r in sample_dt.iterrows()} if keep else {}
        rows, cols = P.table_rows(json_dir, samples, calls_mod, info)
        if not rows:
            print('no JSON of the selection is in the folder', file=err)
            return 2
        path = a.out or f'{json_dir}panel_calls.tsv'
        import csv
        with open(path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore')
            w.writeheader()
            w.writerows(rows)
        print(P.summary_text(rows, cols), file=out)
        print(f'\nwrote {path}', file=out)
        return 0

    from cichlid_common import sheet_writer as SW
    tab = a.tab or pcfg['sheet']['tab']
    if worksheet_fn is None:
        if a.sample_table or a.local_master or a.cloud_master:
            print('writing needs the live sheet: do not combine it with offline options', file=err)
            return 2
        if fm_factory is None:
            from cichlid_common.file_manager import FileManager as fm_factory
        fm = fm_factory()
        worksheet_fn = lambda: fm.gc.open_by_key(fm.s_ID).worksheet(tab)
    eco_of = None
    if pcfg['sheet'].get('scope'):
        try:
            eco_of = P.ecogroups_of(worksheet_fn(), pcfg['sheet']['key'])
        except (ValueError, SW.SheetWriteError) as e:
            print(f'sheet.scope needs the Ecogroup of every sample, and the {tab} tab cannot give it: {e}', file=err)
            return 3
    updates, flag_of, info = P.collect_updates(json_dir, samples, calls_mod, pcfg['sheet']['write_columns'], pcfg['sheet'], eco_of)
    print(f"{len(samples) - len(info['no_json'])} samples with a JSON in {json_dir}; {len(updates)} have something to write"
          + (f"; {len(info['no_json'])} selected samples have no JSON ({', '.join(info['no_json'][:6])}{', ...' if len(info['no_json']) > 6 else ''})" if info['no_json'] else ''), file=out)
    for s, why in info['refused'][:10]:
        print(f'  REFUSED {s}: {why}', file=out)
    if len(info['refused']) > 10:
        print(f"  ... and {len(info['refused']) - 10} more refused", file=out)
    if not updates:
        print('nothing to write', file=err)
        return 2
    cu = SW.ColumnUpdater(worksheet_fn, key=pcfg['sheet']['key'], log=lambda m: print(m, file=out), column_order=P.column_order(pcfg['sheet']['write_columns'], flag_of))
    try:
        plan, dropped = P.drop_flags_of_conflicts(cu, updates, flag_of)
    except SW.SheetWriteError as e:
        print(f'the {tab} tab cannot be written: {e}', file=err)
        return 3
    print(f'SampleDatabase, tab {tab!r}: calls of {len(updates)} samples', file=out)
    print(plan.format(), file=out)
    for col, n in sorted(info['scoped_out'].items()):
        print(f"  {col}: left out for {n} samples whose Ecogroup is not {' or '.join(pcfg['sheet']['scope'][col])} (sheet.scope: module 04 validated it only there)", file=out)
    nc = pcfg['sheet'].get('nocall') or {}
    if nc:
        n_nc = sum(1 for u in updates.values() for v in u.values() if v == nc['value'])
        if n_nc:
            print(f"  {n_nc} genotype cells of the calls get {nc['value']!r} (no call or low GQ: too little data; a site_error stays empty)", file=out)
    if not pcfg['sheet'].get('flags_in_sheet', True):
        print('  the flags are not written to the sheet (sheet.flags_in_sheet: false): they stay in the JSONs and in `panel table`', file=out)
    if dropped:
        print(f'  {dropped} flags cells are left alone because the value next to them conflicts with what the sheet holds', file=out)
    if a.limit is not None:
        first = sorted({c['sample'] for c in plan.writable})[:a.limit]
        print(f"--limit {a.limit}: this run writes only {len(first)} samples ({', '.join(first[:6])}{', ...' if len(first) > 6 else ''}): {sum(1 for c in plan.writable if c['sample'] in first)} cells", file=out)
    if not a.apply:
        print('\nnothing was written: add --apply to write the cells above', file=out)
        return 0
    if not plan.writable:
        print('\nnothing to write', file=out)
        return 0
    import datetime
    tags = '+'.join(info.get('tags') or [])                       # the call settings that made the cells: with flags_in_sheet false the journal is where they are recorded
    journal = os.path.join(json_dir, 'sheet_writes', f"{tab.replace(' ', '_')}.{datetime.datetime.now().strftime('%Y%m%d-%H%M%S')}{'.' + tags if tags else ''}.tsv")
    try:
        rep = cu.apply(updates, overwrite=False, limit=a.limit, journal=journal)
    except SW.SheetWriteError as e:
        print(f'the write stopped: {e}', file=err)
        return 3
    print(f"\nwrote {rep.written} cells" + (f" (columns added: {', '.join(rep.added_columns)})" if rep.added_columns else '') + f'; journal {journal}', file=out)
    if not rep.ok:
        print(f'\nVERIFICATION FAILED for {len(rep.failures)} cells: the sheet changed while writing, or a write did not take. Nothing was undone automatically; see the journal (old values) and:', file=out)
        for f in rep.failures[:10]:
            print(f"  {SW.a1(f['row'], f['col'])} {f['sample']} {f['column']}: {f['problem']}", file=out)
        return 1
    print('every written cell was read back and is in the right row', file=out)
    return 0
