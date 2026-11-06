"""ENA helpers for downloadShortReadData_v2.py / grab_ena_v2.py. Every external call (curl, ascp) goes through an injectable `run`, so the logic is unit tested without a network.

What the old scripts did that lost samples, and what is done here instead:
  * FASTQ files were taken by POSITION from fastq_ftp (fq1 = first, fq2 = second). ENA lists an unpaired file (<run>.fastq.gz) FIRST for runs that have orphan reads, so such runs got the orphans
    and mate 1 and then looked for <run>_2.fastq.gz, which was never downloaded. Files are now chosen by NAME (_1 / _2).
  * ENA answers were fetched with three bare retries without a pause; ENA rate-limits. Retries now back off, and a failure is a recorded status, not a crash of the whole batch.
  * Downloads were checked only by the exit code of ascp. ENA's md5 and size are now verified; a truncated mate file can no longer be silently zipped against a complete one.
  * Run lists came only from an SRA Run Selector file; ENA-only projects (PRJEB...) can be listed directly from ENA."""
import glob, gzip, hashlib, io, os, re, shutil, subprocess, time
import pandas as pd

RUN_FIELDS = ['run_accession', 'sample_accession', 'secondary_sample_accession', 'experiment_accession', 'study_accession', 'library_name', 'library_layout', 'library_source',
              'library_selection', 'instrument_platform', 'instrument_model', 'scientific_name', 'read_count', 'base_count', 'fastq_ftp', 'fastq_aspera', 'fastq_md5', 'fastq_bytes', 'submitted_ftp', 'submitted_format']


class EnaError(RuntimeError):
    """kind: 'permanent' (a retry will not help), 'throttle' (the server or the network is saying slow down: lower the download concurrency) or 'transient' (try again)."""
    def __init__(self, msg, kind='transient'):
        super().__init__(msg)
        self.kind = kind


def filereport_url(accession, fields=RUN_FIELDS):
    return f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields={','.join(fields)}&format=tsv&limit=0"


def http_get(url, run=None, attempts=5, backoff=15.0, sleep=time.sleep):
    """GET through curl with a User-Agent, a time limit and a growing pause between attempts (ENA answers 429 / 5xx under load)."""
    run = run or subprocess.run
    last = ''
    for k in range(1, attempts + 1):
        r = run(['curl', '-fsSL', '-A', 'Mozilla/5.0', '--max-time', '180', url], capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout
        last = f'curl exit {r.returncode}: {(r.stderr or "")[-200:].strip()}' if r.returncode else 'empty answer'
        if k < attempts:
            sleep(backoff * k)
    raise EnaError(f'{url}: {last} after {attempts} attempts')


def fetch_runs(accession, **kw):
    """All runs of a project / sample / run accession as a DataFrame of strings (one row per run)."""
    text = http_get(filereport_url(accession), **kw)
    return pd.read_csv(io.StringIO(text), sep='\t', dtype=str, keep_default_na=False)


# ----------------------------------------------------------------------------------------------------------------- which files are the mates
def select_fastqs(row):
    """Classify one ENA run record. Returns dict(status, r1, r2, orphans, note). status is one of
       ok | no_fastq (ENA has no FASTQ: note says what it has instead) | single_end | missing_mate (only one of _1/_2) | bad_listing (md5 / size lists do not line up)."""
    layout = str(row.get('library_layout', '')).upper()
    urls = [u for u in str(row.get('fastq_ftp', '')).split(';') if u]
    md5s = [m for m in str(row.get('fastq_md5', '')).split(';') if m]
    sizes = [s for s in str(row.get('fastq_bytes', '')).split(';') if s]
    if not urls:
        sub = str(row.get('submitted_ftp', '')); fmt = str(row.get('submitted_format', ''))
        return dict(status='no_fastq', r1=None, r2=None, orphans=[], note=(f'ENA has the submitted file(s) only ({fmt}): {sub[:200]}' if sub else 'ENA has no files for this run'))
    if not (len(urls) == len(md5s) == len(sizes)):
        return dict(status='bad_listing', r1=None, r2=None, orphans=[], note=f'{len(urls)} urls / {len(md5s)} md5 / {len(sizes)} sizes')
    asp = [a for a in str(row.get('fastq_aspera', '')).split(';') if a]
    files = [dict(name=os.path.basename(u), url=u, md5=m, bytes=int(b), aspera=(asp[i] if len(asp) == len(urls) else u.replace('ftp.sra.ebi.ac.uk/', 'fasp.sra.ebi.ac.uk:/', 1)))
             for i, (u, m, b) in enumerate(zip(urls, md5s, sizes))]
    r1 = [f for f in files if re.search(r'_1\.f(ast)?q\.gz$', f['name'])]
    r2 = [f for f in files if re.search(r'_2\.f(ast)?q\.gz$', f['name'])]
    orphans = [f for f in files if f not in r1 and f not in r2]
    if layout == 'SINGLE':
        return dict(status='single_end', r1=None, r2=None, orphans=files, note='library layout is SINGLE')
    if len(r1) == 1 and len(r2) == 1:
        return dict(status='ok', r1=r1[0], r2=r2[0], orphans=orphans, note=(f'{len(orphans)} unpaired file(s) ignored' if orphans else ''))
    return dict(status='missing_mate', r1=None, r2=None, orphans=files, note=f'files: {[f["name"] for f in files]}')


def ena_spots(row):
    """Read PAIRS (spots) of a run. ENA's `read_count` counts READS: both mates of a pair (measured: ERR13389718 has read_count 109,512,218 and fasterq-dump finds 54,756,109 spots). None if unknown."""
    try:
        n = int(row['read_count'])
    except (KeyError, TypeError, ValueError):
        return None
    return n // 2 if str(row.get('library_layout', '')).upper() == 'PAIRED' else n


def runinfo_from_ena(df):
    """ENA run records -> the columns of an SRA Run Selector file that downloadShortReadData expects. AvgSpotLen = bases per spot (both mates), as in Run Selector."""
    spots = pd.to_numeric(df['read_count'], errors='coerce') / df['library_layout'].str.upper().map(lambda x: 2 if x == 'PAIRED' else 1); bases = pd.to_numeric(df['base_count'], errors='coerce')
    return pd.DataFrame({'Run': df['run_accession'], 'AvgSpotLen': (bases / spots).round().astype('Int64'), 'Bases': bases.astype('Int64'), 'BioProject': df['study_accession'],
                         'BioSample': df['sample_accession'], 'Experiment': df['experiment_accession'], 'Instrument': df['instrument_model'], 'Library Name': df['library_name'],
                         'LibraryLayout': df['library_layout'], 'LibrarySelection': df['library_selection'], 'LibrarySource': df['library_source'], 'Organism': df['scientific_name'],
                         'Platform': df['instrument_platform'], 'SRA Study': df['study_accession']})


# ----------------------------------------------------------------------------------------------------------------- downloads
def full_url(u, scheme='https'):
    """ENA lists paths such as ftp.sra.ebi.ac.uk/vol1/fastq/...; they are served over https as well as ftp."""
    return u if '://' in u else f'{scheme}://{u}'


def md5_file(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def _complete(path, rec, md5=md5_file):
    return os.path.exists(path) and os.path.getsize(path) == rec['bytes'] and md5(path) == rec['md5']


def https_download(rec, dest_dir, attempts=4, backoff=30.0, run=None, sleep=time.sleep, md5=md5_file, log=print):
    """Resumable curl download verified by size and md5; a complete-but-wrong file is deleted and fetched again."""
    run = run or subprocess.run
    os.makedirs(dest_dir, exist_ok=True)
    final = os.path.join(dest_dir, rec['name']); part = final + '.partial'
    if _complete(final, rec, md5):
        return final
    url = rec['url'] if '://' in rec['url'] else 'https://' + rec['url']
    last = ''
    for k in range(1, attempts + 1):
        r = run(['curl', '-fL', '--retry', '3', '--retry-delay', '10', '-C', '-', '-o', part, url], capture_output=True, text=True)
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have == rec['bytes'] and md5(part) == rec['md5']:
            os.replace(part, final); return final
        if have >= rec['bytes']:
            last = f'size/md5 wrong ({have} bytes)'; os.remove(part)
        else:
            last = f'curl exit {r.returncode}: {(r.stderr or "")[-200:].strip()} ({have}/{rec["bytes"]} bytes)'
        log(f'    {rec["name"]}: https attempt {k}/{attempts} failed: {last}')
        if k < attempts:
            sleep(backoff * k)
    raise EnaError(f'{rec["name"]}: {last}')


ASCP_PERMANENT = re.compile(r'No such file|not found|Permission denied|Access denied|Invalid argument|does not exist', re.I)
ASCP_THROTTLE = re.compile(r'too many|limit|refused|busy|server aborted|timed? ?out|timeout|disconnect|connection (reset|refused)|failed to authenticate|session stop|unable to (connect|establish)', re.I)


def classify_ascp_error(stderr):
    """Best-effort reading of ascp's stderr (the patterns are my reading of typical ascp messages, NOT verified against ENA: adjust them when you see the real texts in the logs).
    permanent = the file is not there / not allowed; throttle = the session was refused, aborted or timed out (what too many parallel sessions looks like); otherwise transient."""
    t = str(stderr or '')
    if ASCP_PERMANENT.search(t):
        return 'permanent'
    if ASCP_THROTTLE.search(t):
        return 'throttle'
    return 'transient'


def ascp_command(recs, dest_dir, ascp, key, rate='300m', port=33001, resume=1):
    """ascp for one or several files of the same folder (several sources = ONE session: one SSH/UDP handshake for both mates).
    -Q fair policy, -T no encryption (the data are public), -l target rate PER TRANSFER, -k 1 = resume from the size of a partial file (the old call restarted from zero after every failure)."""
    return [ascp, '-QT', '-l', str(rate), '-P', str(port), '-k', str(resume), '-i', key] + ['era-fasp@' + r['aspera'] for r in recs] + [dest_dir]


def run_monitored(cmd, watch, stall_s=90.0, timeout=None, popen=subprocess.Popen, sleep=time.sleep, now=time.time, poll=5.0):
    """Run `cmd` and watch its progress. watch() returns the number of bytes present so far; if that does not grow for `stall_s` seconds the process is terminated and `stalled` is True.
    A silent ascp is what a firewall that lets the TCP control connection through but drops the UDP data channel looks like: it neither fails nor transfers. `timeout` is a hard limit.
    Returns (returncode, stderr text, stalled)."""
    p = popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    t0 = last_change = now(); last = watch(); stalled = False
    while p.poll() is None:
        sleep(poll)
        cur = watch()
        if cur != last:
            last, last_change = cur, now()
        if now() - last_change > stall_s or (timeout and now() - t0 > timeout):
            stalled = now() - last_change > stall_s
            p.terminate()
            try:
                p.wait(timeout=10)
            except Exception:
                p.kill()
            break
    err = ''
    try:
        err = (p.stderr.read() or b'').decode('utf-8', 'replace')
    except Exception:
        pass
    return p.returncode, err, stalled


def _present(dest_dir, recs):
    tot = 0
    for r in recs:
        for f in glob.glob(os.path.join(dest_dir, r['name']) + '*'):
            if not f.endswith(('.aspx', '.aspera-ckpt')):
                try:
                    tot += os.path.getsize(f)
                except OSError:
                    pass
    return tot


def aspera_download(recs, dest_dir, ascp, key, rate='300m', port=33001, run=None, md5=md5_file, timeout=None, stall_s=90.0, monitor=run_monitored):
    """One ascp session for the files in `recs`. Returns (ok, kind, message): ok only if every file is complete and its md5 matches. A file that is complete but wrong is deleted.
    kind 'stall' = no data for `stall_s` seconds (the session is killed). If `run` is given it is used instead of the monitored process (tests)."""
    if not (ascp and os.path.exists(ascp)):
        return False, 'transient', 'ascp not available'
    os.makedirs(dest_dir, exist_ok=True)
    cmd = ascp_command(recs, dest_dir, ascp, key, rate, port)
    if run is not None:
        try:
            r = run(cmd, capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return False, 'throttle', f'ascp timed out after {timeout} s'
        rc, err, stalled = r.returncode, (r.stderr.decode('utf-8', 'replace') if isinstance(r.stderr, (bytes, bytearray)) else str(r.stderr or '')), False
    else:
        rc, err, stalled = monitor(cmd, lambda: _present(dest_dir, recs), stall_s=stall_s, timeout=timeout)
    if stalled:
        return False, 'stall', f'no data received for {stall_s:.0f} s (TCP 33001 connects; the UDP data channel is probably blocked by a firewall)'
    bad = []
    for rec in recs:
        path = os.path.join(dest_dir, rec['name'])
        if not _complete(path, rec, md5):
            bad.append(rec['name'])
            if os.path.exists(path) and os.path.getsize(path) >= rec['bytes']:
                os.remove(path)                                            # complete-looking but wrong: do not let -k 1 'resume' it
    if not bad and rc == 0:
        return True, '', ''
    kind = classify_ascp_error(err) if rc else 'transient'
    return False, kind, f'ascp exit {rc}, incomplete/wrong: {bad}: {err[-300:].strip()}'


def probe_aspera(rec, dest_dir, ascp, key, rate='300m', seconds=20.0, min_mb=10.0, monitor=run_monitored, now=time.time):
    """Does Aspera actually move data on THIS network? Starts one transfer of `rec`, stops it after `seconds`, and reports (ok, MB received, MB/s). Not ok = nothing (or almost nothing)
    arrived: use https. The partial file is deleted."""
    os.makedirs(dest_dir, exist_ok=True)
    t0 = now()
    monitor(ascp_command([rec], dest_dir, ascp, key, rate), lambda: _present(dest_dir, [rec]), stall_s=seconds, timeout=seconds)
    mb = _present(dest_dir, [rec]) / 1e6
    for f in glob.glob(os.path.join(dest_dir, rec['name']) + '*'):
        try:
            os.remove(f)
        except OSError:
            pass
    return mb >= min_mb, round(mb, 1), round(mb / max(now() - t0, 1e-6), 1)


def download_run(sel, dest_dir, ascp=None, key=None, rate='300m', both_in_one_session=True, log=print, timeout=None, run=None, stall_s=90.0, ascp_monitor=None, **kw):
    """Both mate files. Aspera first (one session for both mates), https as the fallback (UDP/33001 is blocked on some networks, and a refused session should not end the run).
    Returns (path1, path2, stats). Raises EnaError(kind=...) when both routes fail; the kind is the Aspera one if it was throttling, so that the scheduler can slow down."""
    recs = [sel['r1'], sel['r2']]
    t0 = time.time(); method = 'https'; akind = ''; stalled = False
    if ascp:
        groups = [recs] if both_in_one_session else [[recs[0]], [recs[1]]]
        ok_all, msgs = True, []
        for g in groups:
            ok, kind, msg = aspera_download(g, dest_dir, ascp, key, rate, run=(run if ascp_monitor is None else None), timeout=timeout, stall_s=stall_s, monitor=ascp_monitor or run_monitored)
            if not ok:
                ok_all = False; akind = kind; stalled = stalled or kind == 'stall'; msgs.append(msg); log(f'    aspera failed ({kind}): {msg[:200]}')
        if ok_all:
            method = 'aspera'
        elif akind == 'permanent':
            raise EnaError('; '.join(msgs), 'permanent')
    if method != 'aspera':
        try:
            for rec in recs:
                https_download(rec, dest_dir, log=log, run=run, **kw)
        except EnaError as e:
            raise EnaError(str(e), 'throttle' if akind == 'throttle' else 'transient')
    secs = max(time.time() - t0, 1e-6); nbytes = sum(r['bytes'] for r in recs)
    return os.path.join(dest_dir, recs[0]['name']), os.path.join(dest_dir, recs[1]['name']), dict(method=method, seconds=round(secs, 1), bytes=nbytes, mb_per_s=round(nbytes / 1e6 / secs, 1), aspera_stalled=stalled)


# ----------------------------------------------------------------------------------------------------------------- FASTQ clean-up (only when FastqToSam rejects the files)
def sanitize_pair(fq1, fq2, out1, out2):
    """Stream both gzipped FASTQs and keep a pair only if BOTH mates have a name, a sequence and a quality of equal, non-zero length (what FastqToSam chokes on). Names get the same ' 1:N:0:2' /
    ' 2:N:0:2' suffix as in the old script. Unlike zip(), raises if one mate file has more reads than the other: a truncated download must not pass as a shorter, valid sample.
    Returns (pairs read, pairs kept)."""
    n = kept = 0
    with gzip.open(fq1, 'rt') as a, gzip.open(fq2, 'rt') as b, open(out1, 'w') as oa, open(out2, 'w') as ob:
        while True:
            ra = [a.readline() for _ in range(4)]; rb = [b.readline() for _ in range(4)]
            if not ra[0] and not rb[0]:
                break
            if not ra[0] or not rb[0]:
                raise EnaError(f'{os.path.basename(fq1)} and {os.path.basename(fq2)} have different numbers of reads (read {n} pairs) although both passed their md5: inconsistent files on ENA', 'permanent')
            n += 1
            na, nb = ra[0][1:].split()[0] if len(ra[0]) > 1 else '', rb[0][1:].split()[0] if len(rb[0]) > 1 else ''
            sa, sb, qa, qb = ra[1].strip(), rb[1].strip(), ra[3].strip(), rb[3].strip()
            if not (na and nb and sa and sb and qa and qb) or len(sa) != len(qa) or len(sb) != len(qb):
                continue
            base = re.sub(r'/[12]$', '', na)
            oa.write(f'@{base} 1:N:0:2\n{sa}\n+\n{qa}\n'); ob.write(f'@{base} 2:N:0:2\n{sb}\n+\n{qb}\n'); kept += 1
    return n, kept


EMPTY_READ_ERRORS = re.compile(r'empty|zero|length|differ|quality|Fastq record|Sequence', re.I)


ASCP_CANDIDATES = ['~/bin/ascp', '~/.aspera/sdk/ascp', '~/.aspera/ascli/sdk/ascp', '~/.aspera/connect/bin/ascp', '/home/mcgrath-lab/bin/ascp']
KEY_CANDIDATES = ['~/.aspera/sdk/aspera_bypass_rsa.pem', '~/.aspera/ascli/sdk/aspera_bypass_rsa.pem', '~/.aspera/connect/etc/asperaweb_id_dsa.openssh']


def find_first(explicit, candidates, which=None, exists=os.path.exists, env=None):
    """The Aspera client / key: the explicit path if it exists, else $ASCP (client only), `which ascp`, then the usual places. None if nothing is found (the caller then uses https and says so)."""
    paths = [explicit] if explicit else []
    if env:
        paths.append(os.environ.get(env))
    if which:
        paths.append(shutil.which(which))
    paths += candidates
    for p in paths:
        if p and exists(os.path.expanduser(p)):
            return os.path.expanduser(p)
    return None
