"""Manage the DownloadDatabase spreadsheet (one tab per Run Info file). The downloader (downloadShortReadData_v2.py --manifest TAB ...) reads the runs from it and writes status back.

  python3 downloadManifest.py import Blumer BlumerPaperData.csv        # create / extend the tab 'Blumer' from a Run Info file (a file name in .../Reads/ReadDownloadFiles/ on Dropbox, or a local path)
  python3 downloadManifest.py import Malinsky MalinskyPaperData.csv
  python3 downloadManifest.py status                                    # per tab: how many runs are pending / downloading / uploaded / in_database / failed / skipped, stale claims, sources
  python3 downloadManifest.py reconcile Blumer Malinsky                 # after --add-pending (or any change to the database): runs now in DNAReadsDatabase become in_database
  python3 downloadManifest.py release Blumer --hours 6                  # give runs whose claim is older than 6 h back (after a crash)
  python3 downloadManifest.py tabs

import never touches a row that is already in the tab (so a person's notes, Skip flags and status edits survive a re-import): it only appends new runs. A run that is already in
DNAReadsDatabase starts as in_database (Source legacy), one whose uBAM is on Dropbox without a database row as uploaded, every other as pending."""
import argparse, os, sys, time
import pandas as pd
from download_reads import manifest as MF


def read_runinfo(fm, src):
    path = src if os.path.exists(src) else fm.localReadDownloadDir + src
    if not os.path.exists(path):
        fm.downloadData(fm.localReadDownloadDir + src)                       # from Dropbox, as the old script did
    ri = pd.read_csv(path, dtype=str, keep_default_na=False)
    bad = [c for c in ('Run', 'BioSample') if c not in ri.columns]
    if bad:
        raise SystemExit(f'{src}: columns {bad} are required (found {list(ri.columns)})')
    return ri


def cloud_runs(fm, projects):
    out = set()
    for pr in sorted(projects):
        out |= {n.split('.')[0] for n in fm.returnCloudFiles(fm.localReadsDir + pr + '/') if n.endswith('.unmapped_marked_adapters.bam')}
    return out


def cmd_import(fm, tab, src, log=print):
    ri = read_runinfo(fm, src)
    fm.readSampleDatabase()
    in_db = set(fm.reads_dt['RunID'])
    on_cloud = cloud_runs(fm, set(ri['BioProject'])) - in_db if 'BioProject' in ri.columns else set()
    mf = fm.downloadManifest(tab)
    res = mf.import_runinfo(ri, in_db, on_cloud)
    mf.flush(force=True)
    log(f"{tab}: {res['added']} runs added, {res['existing']} were already in the tab; of the {len(ri)} runs {len(set(ri['Run']) & in_db)} are in the database and {len(set(ri['Run']) & on_cloud)} are on the cloud without a row")
    return res


def cmd_status(fm, tabs=None, hours=6.0, log=print, now=time.time):
    tabs = tabs or [t for t in fm.downloadManifest('_').b.tabs()]
    out = {}
    for tab in tabs:
        mf = fm.downloadManifest(tab)
        if 'Run' not in mf.b.header(tab):                                   # e.g. the empty default 'Sheet1' of a new spreadsheet
            log(f'{tab}: not a manifest tab (no Run column), skipped'); continue
        df = mf.load()
        c = df['DownloadStatus'].replace('', 'pending').value_counts().to_dict()
        stale = [r.Run for r in df.itertuples() if r.DownloadStatus in MF.ACTIVE and (MF.parse_iso(r.ClaimedAt) is None or now() - MF.parse_iso(r.ClaimedAt) > hours * 3600)]
        log(f'{tab}: {len(df)} runs | ' + ', '.join(f'{k} {c.get(k, 0)}' for k in MF.STATUSES) + (f' | {len(stale)} STALE claims (older than {hours:.0f} h; `release` gives them back)' if stale else ''))
        src = df.loc[df['Source'] != '', 'Source'].value_counts().to_dict()
        if src:
            log('   sources: ' + ', '.join(f'{k} {v}' for k, v in src.items()))
        fk = df.loc[df['DownloadStatus'] == 'failed', 'FailureKind'].value_counts().to_dict()
        if fk:
            log('   failed by kind: ' + ', '.join(f'{k or "?"} {v}' for k, v in fk.items()))
        out[tab] = dict(counts=c, stale=stale)
    return out


def cmd_reconcile(fm, tabs, log=print):
    """Bring the sheet up to date with the database and Dropbox: runs in DNAReadsDatabase become in_database (also `uploaded` ones after an --add-pending); uBAMs on Dropbox become uploaded."""
    fm.readSampleDatabase(); in_db = set(fm.reads_dt['RunID']); out = {}
    for tab in tabs:
        mf = fm.downloadManifest(tab); df = mf.load()
        on_cloud = cloud_runs(fm, set(df['BioProject']) - {''}) - in_db if len(df) else set()
        n = mf.reconcile(in_db, on_cloud); mf.flush(force=True); out[tab] = n
        log(f'{tab}: {n} rows updated from the database / Dropbox')
    return out


def cmd_release(fm, tab, hours=6.0, log=print, now=time.time):
    mf = fm.downloadManifest(tab); df = mf.load(); n = 0
    for r in df.itertuples():
        ca = MF.parse_iso(r.ClaimedAt)
        if r.DownloadStatus in MF.ACTIVE and (ca is None or now() - ca > hours * 3600):
            mf.update(r.Run, DownloadStatus='pending', ClaimedBy='', ClaimedAt=''); n += 1
    mf.flush(force=True)
    log(f'{tab}: {n} stale claims released')
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m download_reads manifest', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('import'); p.add_argument('tab'); p.add_argument('source')
    p = sub.add_parser('status'); p.add_argument('tabs', nargs='*'); p.add_argument('--hours', type=float, default=6.0)
    p = sub.add_parser('reconcile'); p.add_argument('tabs', nargs='+')
    p = sub.add_parser('release'); p.add_argument('tab'); p.add_argument('--hours', type=float, default=6.0)
    sub.add_parser('tabs')
    a = ap.parse_args(argv)
    from cichlid_common.file_manager import FileManager as FM
    fm = FM()
    if not hasattr(fm, 'downloadManifest'):
        raise SystemExit('file_manager.py has no downloadManifest: run `python -m download_reads patch-file-manager` first')
    if a.cmd == 'import':
        cmd_import(fm, a.tab, a.source)
    elif a.cmd == 'status':
        cmd_status(fm, a.tabs, a.hours)
    elif a.cmd == 'reconcile':
        cmd_reconcile(fm, a.tabs)
    elif a.cmd == 'release':
        cmd_release(fm, a.tab, a.hours)
    else:
        print('\n'.join(fm.downloadManifest('_').b.tabs()))


if __name__ == '__main__':
    main()
