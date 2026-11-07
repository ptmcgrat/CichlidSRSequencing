"""Download short-read data from ENA into the McGrath sequencing folder as unmapped BAMs.  Replaces downloadShortReadData.py (which stays untouched; its --Local mode is not in v2).

Where the run list comes from (any combination):
  Run_Info_File           an SRA Run Selector CSV on Dropbox, as before
  --ena ACC [ACC ...]     ENA accessions (project PRJEB..., sample SAMEA..., run ERR...): the run list is taken from ENA, no Run Selector file needed
  --samples-file FILE     TSV/CSV with a Biosample (or SampleID) column [and optionally a Bioproject column], e.g. blumer_missing_from_sheet.tsv

Every run is classified BEFORE anything is downloaded (ok / already in the database / no FASTQ on ENA / single end / missing mate / not on ENA) and written to a plan file, so a run can no longer
vanish silently; --dry-run stops there. At the end a status table records what happened to every run and why. A run is added to the sample database only if its uBAM is verified on the cloud.
  python3 downloadShortReadData_v2.py --samples-file blumer_missing_from_sheet.tsv --dry-run
  python3 downloadShortReadData_v2.py --samples-file blumer_missing_from_sheet.tsv --workers 4
  python3 downloadShortReadData_v2.py --samples-file blumer_missing_from_sheet.tsv --no-db        # while another job may overwrite the sample database: rows are kept in pending_db_rows.jsonl
  python3 downloadShortReadData_v2.py --add-pending ~/Temp/.../pending_db_rows.jsonl             # later, from ONE machine, when nothing else writes the database
  python3 downloadShortReadData_v2.py --retry download_status_20261106-0900.tsv          # only the runs that failed in a way a retry may fix

Speed: downloads and conversions run in SEPARATE pools. --dl-workers (default 3) is the most downloads at the same time: one Aspera transfer is built to fill the link, so more of them
compete and get throttled; it shrinks by itself when the server pushes back and grows again when it does not. --workers (default 4) is the number of conversions (FastqToSam, adapter marking, upload),
which start as soon as a run's FASTQs are on disk. --ascp-rate is the target rate PER download (aggregate = rate x downloads). See the README for how to find the best values with the benchmark."""
import argparse, json, os, platform, shutil, subprocess, sys, time
import pandas as pd
from helper_modules import ena_fetch as E, ena_scheduler as SCH, download_manifest as MF

RUNINFO_COLUMNS = ['Run', 'AvgSpotLen', 'Bases', 'BioProject', 'BioSample', 'Experiment', 'Instrument', 'Library Name', 'LibraryLayout', 'LibrarySelection', 'LibrarySource', 'Organism', 'Platform', 'SRA Study']


from helper_modules.ena_fetch import find_first, ASCP_CANDIDATES, KEY_CANDIDATES      # shared with unit_scripts/ena_bench.py


def validate_args(args):
    if not (args.Run_Info_File or args.ena or args.samples_file or args.add_pending or args.manifest):
        raise SystemExit('give a Run_Info_File, --manifest TAB ..., --ena ACC ..., --samples-file FILE or --add-pending JSONL')


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('Run_Info_File', nargs='*', help='one or more Run Info CSVs (in the read-download folder on Dropbox); put them BEFORE any option that takes several values (--ena)')
    p.add_argument('--ena', nargs='+', default=[], metavar='ACC', help='ENA project / sample / run accessions to list runs from')
    p.add_argument('--samples-file', help='TSV/CSV with a Biosample or SampleID column (and optionally Bioproject)')
    p.add_argument('--manifest', nargs='+', metavar='TAB', help='take the runs from these tabs of the DownloadDatabase spreadsheet (see downloadManifest.py) and record status, source, time and failures there')
    p.add_argument('--retry-failed', choices=['no', 'retryable', 'all'], default='retryable', help='manifest mode: which failed runs are tried again (default: those that failed in a retryable way)')
    p.add_argument('--claim-hours', type=float, default=6.0, help='manifest mode: a claim older than this is considered abandoned')
    p.add_argument('--adopt-cloud', action='store_true', help='for runs whose uBAM is already on the cloud but which have no database row: add the row (no download). Without this they are only reported')
    p.add_argument('--no-db', action='store_true', help='do NOT write the sample database: download, upload and verify the uBAMs, and append each verified run to pending_db_rows.jsonl in the temp folder (add them later, from ONE machine, with --add-pending)')
    p.add_argument('--add-pending', metavar='JSONL', help='add the rows of a pending_db_rows.jsonl to the sample database (only for uBAMs that are on the cloud and runs not yet in the database), then exit')
    p.add_argument('--retry', metavar='STATUS_TSV', help='only the runs that failed transiently in a previous status table')
    p.add_argument('--dry-run', action='store_true', help='classify the runs and write the plan; download nothing')
    p.add_argument('--dl-workers', type=int, default=3, help='most downloads at the same time (default 3; adapts downward when the server pushes back)')
    p.add_argument('--workers', type=int, default=4, help='conversions at the same time (FastqToSam, adapter marking, upload; default 4)')
    p.add_argument('--max-ready', type=int, default=None, help='most runs whose FASTQs wait on disk for a free conversion slot (default: --workers); bounds the disk use')
    p.add_argument('--ascp-rate', default='300m', help='Aspera target rate per download (ENA suggests 300m); aggregate = this x the number of downloads')
    p.add_argument('--source', choices=['ena', 'sra', 'auto'], default='auto', help="where the reads come from: ena (Aspera / https), sra (NCBI's copy on AWS, converted with fasterq-dump) or auto: the Aspera probe decides, "
                   "ENA first if it moves data, otherwise SRA first; the other one is the fallback for a run that fails")
    p.add_argument('--sra-connections', type=int, default=8, help='parallel connections of one SRA download (aria2c)')
    p.add_argument('--fq-threads', type=int, default=8, help='fasterq-dump threads per conversion')
    p.add_argument('--sra-require-check', action='store_true', help='SRA route: refuse a run whose first reads cannot be compared with ENA (ENA unreachable) instead of recording "skipped"')
    p.add_argument('--max-runs', type=int, default=None, help='process only the N SMALLEST runs that need downloading (a cheap first real test)')
    p.add_argument('--bulk', action='store_true', help='with --add-pending: add all rows with one read and one write of the two sheets (backed up and verified); only when nothing else writes the database. --dry-run shows what would be added')
    p.add_argument('--no-probe', action='store_true', help='skip the 20 s Aspera probe before the downloads start (the probe decides between Aspera and https)')
    p.add_argument('--probe-seconds', type=float, default=20.0, help='length of the Aspera probe')
    p.add_argument('--stall-seconds', type=float, default=90.0, help='an Aspera transfer that receives no data for this long is stopped and the run falls back to https')
    p.add_argument('--stagger', type=float, default=8.0, help='seconds between starts of downloads (avoids a burst of simultaneous connections)')
    p.add_argument('--max-attempts', type=int, default=3, help='attempts per run before it is recorded as failed')
    p.add_argument('--backoff', type=float, default=60.0, help='seconds before a retry, times the attempt number')
    p.add_argument('--no-one-session', action='store_true', help='download the two mates in separate Aspera sessions (default: one session for both)')
    p.add_argument('--min-free-gb', type=float, default=100.0, help='keep this much free space on the temp disk')
    p.add_argument('--disk-factor', type=float, default=4.0, help='peak temp disk of a run = this x its FASTQ size')
    p.add_argument('--ascp', default=None, help='Aspera client; default: $ASCP, `which ascp`, ~/bin/ascp, ~/.aspera/sdk/ascp, ... (https is the fallback if none is found)')
    p.add_argument('--key', default=None, help='Aspera key; default: ~/.aspera/sdk/aspera_bypass_rsa.pem, ~/.aspera/ascli/sdk/..., ~/.aspera/connect/etc/asperaweb_id_dsa.openssh')
    return p.parse_args(argv)


def read_table(path):
    return pd.read_csv(path, sep=None, engine='python', dtype=str, keep_default_na=False)


def sample_ids_of(df):
    for c in ('Biosample', 'BioSample', 'SampleID', 'sample_accession'):
        if c in df.columns:
            return df[c].astype(str).str.strip(), c
    return df.iloc[:, 0].astype(str).str.strip(), df.columns[0]


# ---------------------------------------------------------------------------------------------------------------- run table
def build_records(args, fm, fetch=E.fetch_runs, log=print):
    """ENA records (one per run) for everything that was asked for. Queries ENA by project where possible (one call for hundreds of runs), by sample otherwise; every accession once."""
    cache = {}
    def cfetch(acc):
        if acc not in cache:
            cache[acc] = fetch(acc)
        return cache[acc]
    frames, wanted_samples = [], None
    args.run_tab, args.run_attempts, args.run_spotlen = {}, {}, {}
    sources = []
    for name in args.Run_Info_File:
        fm.downloadData(fm.localReadDownloadDir + name)
        sources.append((name, pd.read_csv(fm.localReadDownloadDir + name, dtype=str, keep_default_na=False)))
    for tab, mf in (getattr(args, 'manifests', None) or {}).items():
        todo_rows = mf.runs_to_download(retry_failed={'no': False, 'retryable': True, 'all': 'all'}[args.retry_failed], ttl_s=args.claim_hours * 3600)
        sources.append((f'tab {tab}', todo_rows[MF.RUNINFO_COLUMNS].reset_index(drop=True)))
        args.run_tab.update({r: tab for r in todo_rows['Run']})
        args.run_attempts.update({r: int(a) if str(a).isdigit() else 0 for r, a in zip(todo_rows['Run'], todo_rows['Attempts'])})
    for name, ri in sources:
        if 'AvgSpotLen' in ri.columns:
            args.run_spotlen.update(dict(zip(ri['Run'], ri['AvgSpotLen'])))
        bad = [c for c in ('Run', 'BioSample') if c not in ri.columns]
        if bad:
            raise SystemExit(f'{name}: columns {bad} are required')
        proj = ri['BioProject'] if 'BioProject' in ri.columns else ri.get('SRA Study', pd.Series(dtype=str))
        wanted_runs = set(ri['Run'])
        part = [cfetch(acc) for acc in (sorted(set(proj[proj != ''])) or sorted(wanted_runs))]
        df = pd.concat(part, ignore_index=True) if part else pd.DataFrame()
        df = df[df['run_accession'].isin(wanted_runs)] if len(df) else df
        missing = sorted(wanted_runs - set(df['run_accession'])) if len(df) else sorted(wanted_runs)
        frames.append(df)
        log(f'  {name}: {len(wanted_runs)} runs, {len(df)} found on ENA' + (f', {len(missing)} not returned by the project query (asked for one by one)' if missing else ''))
        for r in missing:                                               # runs the project query did not return: ask for them one by one
            try:
                frames.append(cfetch(r))
            except E.EnaError as e:
                log(f'  {r}: not retrievable from ENA ({str(e)[-120:]})')
    for acc in args.ena:
        frames.append(cfetch(acc))
    if args.samples_file:
        sf = read_table(args.samples_file)
        ids, col = sample_ids_of(sf)
        wanted_samples = set(ids[ids != ''])
        projects = sorted(set(sf['Bioproject'][sf['Bioproject'] != ''])) if 'Bioproject' in sf.columns else []
        got = pd.DataFrame()
        for pr in projects:
            part = cfetch(pr)
            frames.append(part[part['sample_accession'].isin(wanted_samples) | part['secondary_sample_accession'].isin(wanted_samples)])
        have = set(pd.concat(frames, ignore_index=True)['sample_accession']) if frames else set()
        for s in sorted(wanted_samples - have):
            try:
                frames.append(cfetch(s))
            except E.EnaError as e:
                log(f'  {s}: not retrievable from ENA ({str(e)[-120:]})')
    frames = [f for f in frames if len(f)]
    out = pd.concat(frames, ignore_index=True).drop_duplicates('run_accession').reset_index(drop=True) if frames else pd.DataFrame(columns=E.RUN_FIELDS)
    if len(out):
        out['avg_spot_len'] = out['run_accession'].map(args.run_spotlen)
    return out, wanted_samples


def classify(records, reads_in_db, wanted_samples=None, retry_runs=None, on_cloud=()):
    """The plan: one row per run with a status. Nothing is dropped silently."""
    rows = []
    for _, rec in records.iterrows():
        sel = E.select_fastqs(rec)
        status = sel['status']
        if rec['run_accession'] in reads_in_db:
            status = 'in_database'
        elif rec['run_accession'] in on_cloud:
            status = 'on_cloud'                                       # uBAM already on Dropbox, no database row: not downloaded again
        if retry_runs is not None and rec['run_accession'] not in retry_runs:
            status = 'not_in_retry_list'
        rows.append(dict(run=rec['run_accession'], sample=rec['sample_accession'], project=rec['study_accession'], organism=rec['scientific_name'], platform=rec['instrument_platform'],
                         layout=rec['library_layout'], status=status, note=sel['note'], fastq_gb=(sel['r1']['bytes'] + sel['r2']['bytes']) / 1e9 if sel['status'] == 'ok' else 0.0, sel=sel, rec=rec))
    plan = pd.DataFrame(rows)
    if wanted_samples:                                                  # samples that were asked for but have no run on ENA at all
        seen = set(records['sample_accession']) | set(records['secondary_sample_accession']) if len(records) else set()
        for s in sorted(wanted_samples - seen):
            plan = pd.concat([plan, pd.DataFrame([dict(run='', sample=s, project='', organism='', platform='', layout='', status='not_on_ena', note='no run found on ENA for this sample', fastq_gb=0.0, sel=None, rec=None)])], ignore_index=True)
    return plan


# ---------------------------------------------------------------------------------------------------------------- running
def database_rows(fm, t):
    rec = t['rec']
    file_locations = f"{t['project']}/{t['run']}.unmapped_marked_adapters.bam"
    read_data = {'SampleID': t['sample'], 'ProjectID': t['project'], 'RunID': t['run'], 'ReadLength': _spot_len(rec), 'TotalBases': rec.get('base_count'),
                 'Instrument': rec.get('instrument_model'), 'LibraryID': t['library'], 'LibraryLayout': rec.get('library_layout'), 'LibrarySource': rec.get('library_source'),
                 'Platform': t['platform'], 'FileLocations': file_locations}
    try:
        read_data['FileSize'] = fm.returnFileSize(fm.localReadsDir + file_locations)          # rclone size on the cloud copy; raises while Dropbox has not listed the file yet
    except Exception:
        read_data['FileSize'] = ''                                                           # filled in by --add-pending
    sample_data = {'SampleID': t['sample'], 'Sex': '', 'Species': t['organism'], 'DoB': '', 'BroodID': '', 'Parents': '', 'Ecogroup': '', 'LabReared': '', 'Inversion10': ''}
    return read_data, sample_data


def add_to_database(fm, t, st, no_db=False, log=print):
    """The row of a run whose uBAM is verified on the cloud: into the sample database, or (no_db, or if the write fails) appended to pending_db_rows.jsonl for --add-pending.
    FileManager.addDNAReadRow re-reads both sheets and rewrites them completely, and asserts that the run is not yet in the database: any exception of that is caught here, so a run
    whose uBAM is on the cloud never loses its row."""
    read_data, sample_data = database_rows(fm, t)
    st['ubam_bytes'], st['ubam_path'], st['db_written'] = read_data['FileSize'], read_data['FileLocations'], False
    if not no_db and read_data['FileSize'] != '':
        try:
            fm.addDNAReadRow(read_data, sample_data)
            st['db_written'] = True
            return
        except Exception as e:
            st['db_pending'] = f'{type(e).__name__}: {str(e)[-200:]}'
            log(f"  {t['run']}: database write failed ({st['db_pending']}); the row is kept in pending_db_rows.jsonl")
    elif not no_db:
        st['db_pending'] = 'file size not available yet'
    if st.get('db_pending'):
        st['stage'] = 'done (row pending)'
    with open(fm.localTempDir + 'pending_db_rows.jsonl', 'a') as f:
        f.write(json.dumps(dict(read_data=read_data, sample_data=sample_data), default=str) + '\n')


READ_KEYS = {'SampleID', 'ProjectID', 'RunID', 'ReadLength', 'TotalBases', 'Instrument', 'LibraryID', 'LibraryLayout', 'LibrarySource', 'Platform', 'FileLocations', 'FileSize'}
SAMPLE_KEYS = {'SampleID', 'Sex', 'Species', 'DoB', 'BroodID', 'Parents', 'Ecogroup', 'LabReared', 'Inversion10'}


def _ensure_rows(fm, tab, need, log):
    """set_with_dataframe does not enlarge the grid of a worksheet: make sure it has `need` rows (header included)."""
    try:
        ws = fm.gc.open_by_key(fm.s_ID).worksheet(tab)
        if ws.row_count < need:
            old = ws.row_count; ws.resize(rows=need + 50); log(f'  {tab}: grid enlarged from {old} to {need + 50} rows')
    except Exception as e:
        log(f'  {tab}: could not check the grid size ({type(e).__name__}: {str(e)[:100]}); continuing')


def _cloud_listing(fm, folder, run=None):
    """{file name: size in bytes} of a folder on Dropbox with ONE `rclone lsjson`. (FileManager.returnCloudFiles is an uncached `rclone lsf` of the whole folder and returnFileSize an `rclone size`: one call of each per run,
    which is many minutes for a few hundred runs.)"""
    cloud = folder.replace(fm.localMasterDir, fm.cloudMasterDir, 1)
    r = (run or subprocess.run)(['rclone', 'lsjson', '--files-only', cloud], capture_output=True, text=True)
    if r.returncode == 3:                                                                     # rclone: directory not found
        return {}
    if r.returncode != 0:
        raise RuntimeError(f'rclone lsjson {cloud} failed (exit {r.returncode}): {(r.stderr or "")[-200:]}')
    return {x['Name']: int(x['Size']) for x in json.loads(r.stdout or '[]')}


def _default_listing(fm):
    if hasattr(fm, 'cloudMasterDir') and hasattr(fm, 'localMasterDir'):
        return lambda folder: _cloud_listing(fm, folder)
    return lambda folder: {n: None for n in fm.returnCloudFiles(folder)}                    # a file manager without cloud paths (the tests' fakes): sizes are asked per run


def add_pending_bulk(fm, path, log=print, dry_run=False, listing=None):
    """Add the rows of a pending file with ONE read and ONE write of the two sheets (add_pending adds run by run: two sheet reads and two whole-sheet writes per run, which hits Google's 60 reads per minute).
    The rules are those of FileManager.addDNAReadRow: a run is added only if it is not in DNAReadsDatabase and its uBAM is on the cloud; a sample is added to SampleDatabase only if it is not there yet.
    Safe only when nothing else writes the two sheets meanwhile. Both sheets are saved to CSV first, nothing is written if any line is malformed, and the result is read back and checked."""
    listing = listing or _default_listing(fm)
    fm.readSampleDatabase()
    old_reads, old_samples = fm.reads_dt.copy(), fm.sample_dt.copy()
    in_db, have_samples = set(old_reads['RunID']), set(old_samples['SampleID'])
    new_reads, new_samples, skipped, missing, warned, bad, listed = [], [], 0, [], [], [], {}
    for n, line in enumerate(open(path), 1):
        if not line.strip():
            continue
        try:
            row = json.loads(line); rd, sd = dict(row['read_data']), dict(row['sample_data'])
        except (ValueError, KeyError, TypeError):
            bad.append(f'line {n}: not a pending row'); continue
        if set(rd) != READ_KEYS or set(sd) != SAMPLE_KEYS:
            bad.append(f"line {n} ({rd.get('RunID', '?')}): unexpected fields {sorted(set(rd) ^ READ_KEYS)} {sorted(set(sd) ^ SAMPLE_KEYS)}"); continue
        if rd['RunID'] in in_db:
            skipped += 1; continue
        folder = fm.localReadsDir + os.path.dirname(rd['FileLocations']) + '/'
        if folder not in listed:
            listed[folder] = listing(folder); log(f'  Reads/{os.path.dirname(rd["FileLocations"])}/: {len(listed[folder])} files on Dropbox')
        names, base = listed[folder], os.path.basename(rd['FileLocations'])
        if base not in names:
            missing.append(rd['RunID']); continue
        if rd.get('FileSize') in ('', None) and not dry_run:                                   # a dry run needs no sizes
            rd['FileSize'] = names[base] if names[base] is not None else fm.returnFileSize(fm.localReadsDir + rd['FileLocations'])
        new_reads.append(rd); in_db.add(rd['RunID'])
        if sd['SampleID'] in have_samples:
            if sd['SampleID'] not in {x['SampleID'] for x in new_samples}:
                warned.append(sd['SampleID'])                                                  # already in the sheet before this run: addDNAReadRow only warns
        else:
            new_samples.append(sd); have_samples.add(sd['SampleID'])
    if bad:
        raise SystemExit('nothing was written; the pending file has malformed lines:\n  ' + '\n  '.join(bad[:10]))
    log(f'pending rows (bulk): {len(new_reads)} to add, {skipped} already in the database, {len(missing)} not added because the uBAM is not on the cloud {missing[:10]}; '
        f'{len(new_samples)} new samples ({len(warned)} runs belong to samples that are already in SampleDatabase)')
    if dry_run or not new_reads:
        log('dry run: nothing written' if dry_run else 'nothing to add'); return len(new_reads), skipped, missing
    stamp = time.strftime('%Y%m%d-%H%M%S'); bdir = f'{fm.localTempDir}database_backup_{stamp}/'; os.makedirs(bdir, exist_ok=True)
    old_samples.to_csv(bdir + 'SampleDatabase.csv', index=False); old_reads.to_csv(bdir + 'DNAReadsDatabase.csv', index=False)
    log(f'  both sheets saved to {bdir} (Google Sheets also keeps its version history)')
    fm.reads_dt = pd.concat([old_reads, pd.DataFrame(new_reads)], ignore_index=True)
    fm.sample_dt = pd.concat([old_samples, pd.DataFrame(new_samples)], ignore_index=True) if new_samples else old_samples
    assert fm.reads_dt['RunID'].is_unique and list(fm.reads_dt['RunID'][:len(old_reads)]) == list(old_reads['RunID'])
    assert list(fm.sample_dt['SampleID'][:len(old_samples)]) == list(old_samples['SampleID'])
    _ensure_rows(fm, 'SampleDatabase', len(fm.sample_dt) + 1, log); _ensure_rows(fm, 'DNAReadsDatabase', len(fm.reads_dt) + 1, log)
    if not fm.setSampleDatabase():
        raise SystemExit(f'the write to the sheets failed after retries; they may be partly written: restore from {bdir} (or the sheet\'s version history), then run again')
    want_r, want_s = list(fm.reads_dt['RunID']), list(fm.sample_dt['SampleID'])
    fm.readSampleDatabase()                                                                    # read back and check
    if list(fm.reads_dt['RunID']) != want_r or list(fm.sample_dt['SampleID']) != want_s:
        raise SystemExit(f'the sheets do not hold what was written (expected {len(want_r)} runs and {len(want_s)} samples, found {len(fm.reads_dt)} and {len(fm.sample_dt)}): see {bdir}')
    log(f'  verified: DNAReadsDatabase {len(old_reads)} -> {len(fm.reads_dt)} rows, SampleDatabase {len(old_samples)} -> {len(fm.sample_dt)} rows')
    return len(new_reads), skipped, missing


def add_pending(fm, path, log=print):
    """Add the rows of a pending file to the database. A row is added only if its uBAM is listed on the cloud and its run is not yet in the database; the rest is reported."""
    in_db = set(fm.reads_dt['RunID']); added, skipped, missing = 0, 0, []
    for line in open(path):
        if not line.strip():
            continue
        row = json.loads(line); rd, sd = row['read_data'], row['sample_data']
        if rd['RunID'] in in_db:
            skipped += 1; continue
        folder = fm.localReadsDir + os.path.dirname(rd['FileLocations']) + '/'
        if os.path.basename(rd['FileLocations']) not in set(fm.returnCloudFiles(folder)):
            missing.append(rd['RunID']); continue
        if rd.get('FileSize') in ('', None):
            rd['FileSize'] = fm.returnFileSize(fm.localReadsDir + rd['FileLocations'])
        fm.addDNAReadRow(rd, sd); in_db.add(rd['RunID']); added += 1
    log(f'pending rows: {added} added, {skipped} already in the database, {len(missing)} not added because the uBAM is not on the cloud {missing[:10]}')
    return added, skipped, missing


def _spot_len(rec):
    """ReadLength of the database row: the Run Info's AvgSpotLen when there is one (it is what the old script used). ENA's read_count means reads for some runs and pairs for others, so no length is guessed from it."""
    v = rec.get('avg_spot_len') if hasattr(rec, 'get') else None
    try:
        return int(float(v)) if v not in (None, '') and v == v else ''
    except (TypeError, ValueError):
        return ''


def report(plan, results, path, log=print):
    rows = []
    for _, p in plan.iterrows():
        r = results.get(p['run'])
        if r:
            rows.append(dict(run=p['run'], sample=p['sample'], project=p['project'], plan=p['status'], state=r['state'], stage=r['stage'], kind=r.get('kind', ''), attempts=r.get('attempts', ''), error=(r.get('error') or r.get('db_pending') or '')[:300].replace('\n', ' ')))
        else:
            rows.append(dict(run=p['run'], sample=p['sample'], project=p['project'], plan=p['status'], state='skipped' if p['status'] in ('in_database', 'not_in_retry_list', 'on_cloud') else ('planned' if p['status'] == 'ok' else 'not_downloadable'),
                             stage='', kind='', attempts='', error=p['note']))
    t = pd.DataFrame(rows)
    t.to_csv(path, sep='\t', index=False)
    log('\nsummary (state / plan status / stage):')
    log(t.groupby(['state', 'plan', 'stage']).size().to_string())
    return t


def manifest_events(manifests, now=time.time):
    """Scheduler events -> manifest updates (queued, written in batches). `Source` is the site the reads came from: ENA-aspera / ENA-https (SRA-aws when that source exists)."""
    host = platform.node().split('.')[0]
    def ev(name, t, st):
        mf = manifests.get(t.get('tab'))
        if mf is None:
            return
        run = t['run']; A = lambda n: t.get('base_attempts', 0) + (n or 0)
        if name == 'fetch_start':
            mf.update(run, DownloadStatus='downloading', Attempts=A(t.get('attempts', 0) + 1), Machine=host)
        elif name == 'fetched':
            d = st.get('download') or {}
            mf.update(run, Source=d.get('source') or 'ENA-' + str(d.get('method', '?')), DownloadSeconds=d.get('seconds'), MBPerSec=d.get('mb_per_s'), Verified=True,
                      QualityCheck=d.get('quality_check', 'md5 of the ENA files' if str(d.get('source', 'ENA')).startswith('ENA') else ''), ReadCountCheck=d.get('count_check', ''),
                      ReadsInUBAM=d.get('reads_fastq', ''))
        elif name == 'convert_start':
            mf.update(run, DownloadStatus='converting')
        elif name == 'done':
            mf.update(run, DownloadStatus='in_database' if st.get('db_written') else 'uploaded', InDatabase=bool(st.get('db_written')), DownloadedOn=MF.utc_iso(now()), UBamBytes=st.get('ubam_bytes'),
                      UBamPath=st.get('ubam_path'), ConvertSeconds=round(sum(v for k, v in (st.get('timings') or {}).items() if k in ('ubam_s', 'adapters_s')), 1) if st.get('timings') else '',
                      UploadSeconds=(st.get('timings') or {}).get('upload_s', ''), Attempts=A(st.get('attempts')), FailureKind='', LastError=str(st.get('db_pending') or ''), ClaimedBy='', ClaimedAt='', Machine=host)
        elif name == 'retry':
            mf.update(run, Attempts=A(st.get('attempts')), FailureKind=st.get('kind'), LastError=str(st.get('error', ''))[:250])
        elif name == 'failed':
            mf.update(run, DownloadStatus='failed', FailureKind=st.get('kind'), LastError=f"{st.get('stage', '')}: {str(st.get('error', ''))[:250]}", Attempts=A(st.get('attempts')), ClaimedBy='', ClaimedAt='')
        mf.flush()
    return ev


def sample_summary(plan):
    """Per BioSample: how many of its runs are in the database / on the cloud / to download / not downloadable, and whether the sample is complete."""
    p = plan[plan.run != '']
    g = p.groupby('sample').agg(project=('project', 'first'), organism=('organism', 'first'), runs=('run', 'size'), in_database=('status', lambda x: int((x == 'in_database').sum())),
                                on_cloud=('status', lambda x: int((x == 'on_cloud').sum())), to_download=('status', lambda x: int((x == 'ok').sum())), fastq_gb=('fastq_gb', 'sum'))
    g['not_downloadable'] = g.runs - g.in_database - g.on_cloud - g.to_download
    g['state'] = ['complete' if r.in_database == r.runs else ('none in database' if r.in_database == 0 else 'partial') for r in g.itertuples()]
    return g.reset_index()


def cloud_runs(fm, projects, log=print):
    """Run ids whose uBAM is listed on the cloud (one rclone lsf per project folder)."""
    out = set()
    for pr in sorted(projects):
        for n in fm.returnCloudFiles(fm.localReadsDir + pr + '/'):
            if n.endswith('.unmapped_marked_adapters.bam'):
                out.add(n.split('.')[0])
    return out


def run(args, fm, fetch=E.fetch_runs, launch=None, sleep=time.sleep, free_fn=None, log=print, now=time.time, probe=E.probe_aspera, which=shutil.which):
    fm.readSampleDatabase()
    args.manifests = {}
    if args.manifest:
        who = f'{platform.node().split(".")[0]}-{os.getpid()}'
        for tab in args.manifest:
            args.manifests[tab] = fm.downloadManifest(tab, who)
            args.manifests[tab].setup()
        log(f'manifest tabs {args.manifest} as {who}')
    if args.add_pending:
        if args.bulk or args.dry_run:
            add_pending_bulk(fm, args.add_pending, log, dry_run=args.dry_run)           # --dry-run never writes, with or without --bulk
        else:
            add_pending(fm, args.add_pending, log)
        if args.manifests and not args.dry_run:                                                  # runs added to the database now: their rows in the sheet become in_database
            fm.readSampleDatabase(); in_db_now = set(fm.reads_dt['RunID'])
            for tab, mf in args.manifests.items():
                n = mf.reconcile(in_db_now, set()); mf.flush(force=True); log(f'  manifest {tab}: {n} rows are now in_database')
        return 0
    if args.no_db:
        log('--no-db: the sample database will NOT be written; verified runs are appended to ' + fm.localTempDir + 'pending_db_rows.jsonl')
    records, wanted = build_records(args, fm, fetch, log)
    if not len(records) and not wanted:
        log('nothing to download: no pending runs' + (' in the manifest tabs' if args.manifests else '')); return 0
    in_db = set(fm.reads_dt['RunID'])
    retry_runs = None
    if args.retry:
        prev = read_table(args.retry)
        retry_runs = set(prev.loc[(prev['state'] == 'failed') & (prev['kind'] != 'permanent'), 'run'])
        log(f'retrying {len(retry_runs)} runs that failed in a retryable way (not permanent) in {args.retry}')
    on_cloud = cloud_runs(fm, set(records['study_accession']) if len(records) else set(), log) - in_db
    plan = classify(records, in_db, wanted, retry_runs, on_cloud)
    stamp = time.strftime('%Y%m%d-%H%M%S')
    os.makedirs(fm.localTempDir, exist_ok=True)
    plan.drop(columns=['sel', 'rec']).to_csv(f'{fm.localTempDir}ena_plan_{stamp}.tsv', sep='\t', index=False)
    log('plan: ' + ', '.join(f'{k}: {v}' for k, v in plan.status.value_counts().items()))
    for st in ('no_fastq', 'single_end', 'missing_mate', 'bad_listing', 'not_on_ena'):
        sub = plan[plan.status == st]
        if len(sub):
            log(f'  {st}: {len(sub)} (e.g. {sub.iloc[0]["run"] or sub.iloc[0]["sample"]}: {sub.iloc[0]["note"][:150]})')
    log(f'plan written to {fm.localTempDir}ena_plan_{stamp}.tsv')
    todo = plan[plan.status == 'ok']
    if args.max_runs:
        todo = todo.sort_values('fastq_gb').head(args.max_runs)
        log(f'--max-runs {args.max_runs}: only the {len(todo)} smallest runs are processed: {list(todo.run)}')
    if args.manifests and not args.dry_run:
        for tab, mf in args.manifests.items():
            mine = plan[plan.run.map(lambda r: args.run_tab.get(r) == tab)]
            for _, r in mine.iterrows():
                if r.status in ('no_fastq', 'single_end', 'missing_mate', 'bad_listing', 'not_on_ena'):
                    mf.update(r.run, DownloadStatus='skipped', FailureKind='permanent', LastError=f'{r.status}: {r.note}'[:250])
                if r.rec is not None and r.status == 'ok':
                    mf.update(r.run, ReadsSource=r.rec.get('read_count', ''))
            n = mf.reconcile(in_db, on_cloud); mf.flush(force=True)
            log(f'  manifest {tab}: {len(mine)} runs planned, {n} rows reconciled with the database / cloud')
    ss = sample_summary(plan); ss.to_csv(f'{fm.localTempDir}ena_plan_samples_{stamp}.tsv', sep='\t', index=False)
    log('samples: ' + ', '.join(f'{k}: {v}' for k, v in ss.state.value_counts().items()) + f' (of {len(ss)}); per-sample table {fm.localTempDir}ena_plan_samples_{stamp}.tsv')
    log(f'to download: {len(todo)} runs, {todo.fastq_gb.sum():.0f} GB of FASTQ' + (f'; {int((plan.status == "on_cloud").sum())} runs have a uBAM on the cloud but no database row ' +
        ('(adding their rows now: --adopt-cloud)' if args.adopt_cloud else '(use --adopt-cloud to add their rows)') if (plan.status == 'on_cloud').any() else ''))
    ascp = find_first(args.ascp, ASCP_CANDIDATES, which='ascp', env='ASCP'); key = find_first(args.key, KEY_CANDIDATES)
    if ascp and key:
        log(f'Aspera: {ascp} (key {key}), {args.ascp_rate} per transfer')
    else:
        log(f"Aspera NOT used ({'ascp not found' if not ascp else 'ascp found at ' + ascp + ' but no key'}): downloads go over https (md5-verified, slower). Install ascp / pass --ascp and --key.")
        ascp = None
    if ascp and key and len(todo) and not args.dry_run and not args.no_probe:
        ok, mb, rate = probe(todo.iloc[0].sel['r1'], fm.localTempDir + 'aspera_probe/', ascp, key, args.ascp_rate, seconds=args.probe_seconds)
        if ok:
            log(f'Aspera probe: {mb:.0f} MB in {args.probe_seconds:.0f} s ({rate:.1f} MB/s): Aspera works on this network')
        else:
            log(f'Aspera probe: only {mb:.1f} MB in {args.probe_seconds:.0f} s: Aspera moves no data on this network (TCP 33001 connects, the UDP data channel is probably blocked): using https. '
                f'(--no-probe to force Aspera; ask IT to allow outbound UDP 33001 to fasp.sra.ebi.ac.uk)')
            ascp = None
    sources = ['ena'] if args.source == 'ena' else (['sra'] if args.source == 'sra' else (['ena', 'sra'] if ascp else ['sra', 'ena']))
    if len(todo) and not args.dry_run:
        log(f"sources, in the order they are tried per run: {sources}" + (' (Aspera did not move data in the probe, so SRA comes first)' if args.source == 'auto' and sources[0] == 'sra' else ''))
        if 'sra' in sources and launch is None:
            missing = [t for t in ('aria2c', 'vdb-validate', 'fasterq-dump') if not which(t)]
            if missing:
                raise SystemExit(f"the SRA source needs {missing} on PATH: conda install -c conda-forge -c bioconda sra-tools aria2 (or --source ena)")
        if 'sra' in sources:
            args.disk_factor = max(args.disk_factor, 8.0)                       # .sra + uncompressed FASTQ (~4x the gz) + fasterq-dump scratch + the uBAM
    tasks = []
    for _, r in todo.iterrows():
        out = f"{fm.localReadsDir}{r.project}/{r.run}.unmapped_marked_adapters.bam"
        spec = dict(run=r.run, sample=r['sample'], library=r.rec['library_name'] or r.run, platform=r.platform, sel=r.sel, output_bam=out, temp_dir=fm.localTempDir,
                    ascp=ascp, key=key, ascp_rate=args.ascp_rate, one_session=not args.no_one_session, stall_s=args.stall_seconds,
                    source=sources[0], ena_read_count=int(r.rec['read_count']) if str(r.rec.get('read_count', '')).isdigit() else None, paired=str(r.rec.get('library_layout', '')).upper() == 'PAIRED', ena_url_r1=E.full_url(r.sel['r1']['url']),
                    instrument=r.rec.get('instrument_model', ''), fq_threads=args.fq_threads, sra_connections=args.sra_connections, require_quality_check=args.sra_require_check)
        tasks.append(dict(run=r.run, sample=r['sample'], project=r.project, organism=r.organism, platform=r.platform, library=spec['library'], sel=r.sel, rec=r.rec, fastq_gb=r.fastq_gb,
                          output_bam=out, spec=spec))
    results = {}
    if args.adopt_cloud and not args.dry_run:
        for _, r in plan[plan.status == 'on_cloud'].iterrows():
            t = dict(run=r.run, sample=r['sample'], project=r.project, organism=r.organism, platform=r.platform, library=r.rec['library_name'] or r.run, rec=r.rec)
            st = dict(state='done', stage='adopted from cloud', kind='')
            add_to_database(fm, t, st, args.no_db, log); results[r.run] = st
            mf = args.manifests.get(args.run_tab.get(r.run))
            if mf:
                mf.update(r.run, DownloadStatus='in_database' if st.get('db_written') else 'uploaded', InDatabase=bool(st.get('db_written')), Source='legacy', Notes='adopted from cloud')
    if not args.dry_run and tasks:
        here = os.path.dirname(os.path.abspath(__file__))
        launch = launch or (lambda kind, spec, logf: subprocess.Popen([sys.executable, os.path.join(here, 'helper_modules', 'grab_ena_v2.py'), spec], stdout=open(logf, 'w'), stderr=subprocess.STDOUT))
        for t in tasks:
            os.makedirs(fm.localReadsDir + t['project'], exist_ok=True)
        cfg = dict(temp_dir=fm.localTempDir, dl_max=args.dl_workers, workers=args.workers, max_ready=args.max_ready or args.workers, stagger_s=args.stagger, max_attempts=args.max_attempts,
                   backoff_s=args.backoff, min_free_gb=args.min_free_gb, disk_factor=args.disk_factor, status_every_s=60, poll_s=5)
        for t in tasks:
            t['sources'] = sources
            t['tab'] = args.run_tab.get(t['run']); t['base_attempts'] = args.run_attempts.get(t['run'], 0)           # Attempts in the sheet add up over invocations
        results.update(SCH.run_tasks(tasks, cfg, launch, lambda t, st: add_to_database(fm, t, st, args.no_db, log), sleep=sleep, now=now, free_fn=free_fn, log=log,
                                     claim=(lambda t: args.manifests[t['tab']].claim(t['run'], args.claim_hours * 3600) if t.get('tab') else True) if args.manifests else None,
                                     on_event=manifest_events(args.manifests, now) if args.manifests else None))
    for mf in args.manifests.values():
        mf.flush(force=True)
    table = report(plan, results, f'{fm.localTempDir}download_status_{stamp}.tsv', log)
    log(f'status table: {fm.localTempDir}download_status_{stamp}.tsv')
    pend = [k for k, v in results.items() if v.get('db_pending')]
    if pend:
        log(f'{len(pend)} runs are on the cloud but their database rows are pending: python3 downloadShortReadData_v2.py --add-pending {fm.localTempDir}pending_db_rows.jsonl')
    bad = table[(table.state == 'failed')]
    return 1 if len(bad) else 0


def main():
    args = parse_args()
    validate_args(args)
    from helper_modules.file_manager import FileManager as FM
    fm = FM()
    sys.exit(run(args, fm))


if __name__ == '__main__':
    main()
