"""Compare the first reads of a FINISHED uBAM (on Dropbox) with the first reads of ENA's FASTQ, for runs whose SRA quality check was skipped because ENA did not answer at the time.

    python3 unit_scripts/verify_ubam.py --skipped [--tabs Blumer Malinsky] [--write]      every SRA-aws run of the sheet whose QualityCheck starts with "skipped"
    python3 unit_scripts/verify_ubam.py ERR14734399 ERR... --tab Blumer [--write]

Only the first ~1 MB of the uBAM is read (rclone cat | samtools view -f 0x40, stopped early) and the first 3 MB of ENA's FASTQ_1. Read 1 of the uBAM is compared with FASTQ_1 read by read (sequence and quality). A sheet is written
only with --write, and only the QualityCheck cell: `identical-to-ENA (rechecked)`, or `MISMATCH (rechecked)` with a note. A mismatch never changes the status of the run: you decide what to do with it."""
import argparse, os, subprocess, sys
try:
    from helper_modules import ena_fetch as E, sra_fetch as SF
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from helper_modules import ena_fetch as E, sra_fetch as SF


def sam_lines(remote, popen=subprocess.Popen):
    """The read-1 records of a uBAM on Dropbox as SAM lines, streamed; the caller stops reading when it has enough."""
    p1 = popen(['rclone', 'cat', remote], stdout=subprocess.PIPE)
    p2 = popen(['samtools', 'view', '-f', '0x40', '-'], stdin=p1.stdout, stdout=subprocess.PIPE, universal_newlines=True)
    try:
        for line in p2.stdout:
            yield line
    finally:
        for p in (p2, p1):
            try:
                p.kill()
            except Exception:
                pass


def ubam_reads(remote, n=1000, lines=sam_lines):
    out = []
    for line in lines(remote):
        f = line.rstrip('\n').split('\t')
        if len(f) > 10:
            out.append((f[9], f[10]))
        if len(out) >= n:
            break
    return out


def ena_reads(run, n=1000, fetch=E.fetch_runs, head=SF.ena_head):
    df = fetch(run)
    if df is None or not len(df):
        return []
    sel = E.select_fastqs(df.iloc[0])
    if not sel or not sel.get('r1'):
        return []
    return SF.records(head(E.full_url(sel['r1']['url'])), n)


def compare(mine, theirs, min_reads=100):
    k = min(len(mine), len(theirs))
    if k < min_reads:
        return 'skipped (ENA unreachable)', f'only {len(theirs)} ENA reads and {len(mine)} uBAM reads could be compared'
    bad_seq = sum(1 for a, b in zip(mine[:k], theirs[:k]) if a[0] != b[0]); bad_q = sum(1 for a, b in zip(mine[:k], theirs[:k]) if a[1] != b[1])
    if bad_seq == 0 and bad_q == 0:
        return 'identical-to-ENA', f'{k} reads compared'
    return 'MISMATCH', f'{bad_seq} sequences and {bad_q} quality strings differ in the first {k} reads'


def main(argv=None, fm=None, ubam_reader=ubam_reads, ena_reader=ena_reads, out=print):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('runs', nargs='*'); ap.add_argument('--skipped', action='store_true'); ap.add_argument('--tabs', nargs='+', default=['Blumer', 'Malinsky']); ap.add_argument('--tab', default=None)
    ap.add_argument('--write', action='store_true'); ap.add_argument('--n', type=int, default=1000)
    a = ap.parse_args(argv)
    if not (a.runs or a.skipped):
        raise SystemExit('give run accessions (with --tab) or --skipped')
    if fm is None:
        try:
            from helper_modules.file_manager import FileManager
        except ImportError:
            from cichlid_common.file_manager import FileManager
        fm = FileManager()
    todo = []
    for tab in ([a.tab] if a.tab else a.tabs):
        d = fm.readDownloadTab(tab)
        if not len(d):
            continue
        sel = d[d.Run.isin(a.runs)] if a.runs else d[(d.Source == 'SRA-aws') & d.QualityCheck.str.startswith('skipped')]
        todo += [(tab, r) for r in sel.itertuples()]
    out(f'{len(todo)} run(s) to check')
    results = []
    for tab, r in todo:
        remote = fm.cloudMasterDir + 'Reads/' + r.UBamPath
        verdict, detail = compare(ubam_reader(remote, a.n), ena_reader(r.Run, a.n))
        out(f'{tab} {r.Run}: {verdict} ({detail})' + ('   <-- ACTION: this run does not match ENA; do not use it until it is resolved' if verdict == 'MISMATCH' else ''))
        results.append((tab, r.Run, verdict, detail))
        if a.write and not verdict.startswith('skipped'):
            mf = fm.downloadManifest(tab)
            mf.update(r.Run, QualityCheck=verdict + ' (rechecked)', **({'Notes': f'rechecked from the uBAM: {detail}'} if verdict == 'MISMATCH' else {}))
            mf.flush(force=True)
    return results


if __name__ == '__main__':
    main()
