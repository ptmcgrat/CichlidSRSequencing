"""The SRA source: reads from NCBI's copy of the data on AWS (public bucket sra-pub-run-odp, plain https, many parallel connections with aria2c) converted to FASTQ with fasterq-dump.
Used when ENA is slow or unreachable. The run list and the metadata still come from ENA's portal API; only the bytes of the reads come from here.

  aria2c (8 connections, resume) -> vdb-validate -> fasterq-dump --split-3 -> checks -> <run>_1.fastq / <run>_2.fastq (uncompressed) for the same conversion stage as the ENA route

Every run is CHECKED, never trusted:
  * integrity: vdb-validate on the .sra file;
  * completeness: fasterq-dump's `spots read` must equal ENA's read_count (read pairs);
  * quality scores: NCBI also distributes "SRA Lite" files whose base qualities are simplified (not usable for variant calling). The first reads of the converted FASTQ are therefore COMPARED with the first
    reads of ENA's FASTQ (sequence and quality). Identical: recorded `identical-to-ENA`. Different: the run is rejected as `source_unavailable` (another source is tried). ENA not answering (it is slow right
    now): `skipped (ENA unreachable)` is recorded in the manifest, and a run with only ONE or TWO distinct quality symbols (what SRA Lite looks like) is rejected. Binned qualities are NOT a sign of Lite: the
    original HiSeq 2500 data of ERR13389718 have 6 distinct symbols at ENA and at NCBI alike, which is why no instrument-based rule is used.
Failure kinds: transient (try again), throttle, permanent, and `source_unavailable` = this source cannot deliver this run correctly (object missing, counts or qualities differ): the next source is tried."""
import os, re, subprocess, time, zlib
try:
    from helper_modules.ena_fetch import EnaError, run_monitored
except ImportError:
    from ena_fetch import EnaError, run_monitored

S3 = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}'


def sra_url(run):
    return S3.format(run=run)


def aria2_command(url, dest_dir, name, connections=8):
    """Resumable multi-connection download. The log file (warnings only) is where the reason of a failure is read from."""
    return ['aria2c', '-x', str(connections), '-s', str(connections), '-k', '32M', '-c', '--file-allocation=none', '--max-tries=3', '--retry-wait=10', '--console-log-level=error',
            '--summary-interval=0', '--download-result=hide', '--log=' + os.path.join(dest_dir, name + '.aria2.log'), '--log-level=warn', '-d', dest_dir, '-o', name, url]


def validate_command(sra_path):
    return ['vdb-validate', sra_path]


def fasterq_command(sra_path, out_dir, tmp_dir, threads=8):
    return ['fasterq-dump', '--split-3', '-e', str(threads), '-t', tmp_dir, '-O', out_dir, '--force', sra_path]


def parse_counts(text):
    """fasterq-dump's summary: spots read / reads read / reads written (with thousands separators)."""
    out = {}
    for key, label in (('spots', 'spots read'), ('reads_read', 'reads read'), ('reads_written', 'reads written')):
        m = re.search(label + r'\s*:\s*([\d,]+)', text)
        if m:
            out[key] = int(m.group(1).replace(',', ''))
    return out


def gunzip_partial(data):
    """Whatever can be decompressed of a TRUNCATED gzip stream (a range request of the first bytes of an ENA FASTQ)."""
    try:
        return zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(data)
    except zlib.error:
        return b''


def records(text, n):
    """The first n complete FASTQ records of a text as (sequence, quality) pairs."""
    lines = text.split('\n'); out = []
    for i in range(0, len(lines) - 3, 4):
        if len(out) >= n:
            break
        out.append((lines[i + 1], lines[i + 3]))
    return out


def quality_alphabet(path, n_reads=20000):
    seen = set()
    with open(path, errors='replace') as f:
        for i, line in enumerate(f):
            if i % 4 == 3:
                seen.update(line.rstrip('\n'))
            if i >= n_reads * 4:
                break
    return seen


def ena_head(url, nbytes=3_000_000, timeout=60, run=subprocess.run):
    """The first reads of ENA's FASTQ through a range request; [] if ENA does not deliver in time."""
    try:
        r = run(['curl', '-sL', '-r', f'0-{nbytes}', '-m', str(timeout), url], capture_output=True)
        return gunzip_partial(r.stdout or b'').decode('utf-8', 'replace')
    except Exception:
        return ''


def quality_check(fq1_path, ena_url, instrument='', n=1000, fetch=ena_head):
    """Returns (verdict, detail). verdict: 'identical-to-ENA' | 'MISMATCH' | 'skipped (ENA unreachable)' | 'simplified-qualities-suspected'. (instrument is accepted for compatibility and not used.)"""
    theirs = records(fetch(ena_url), n)
    mine = records(open(fq1_path, errors='replace').read(n * 4 * 400 + 4000), n)          # enough text for n reads even at 400 characters per line
    k = min(len(theirs), len(mine))
    if k >= 100:
        bad_seq = sum(1 for a, b in zip(theirs[:k], mine[:k]) if a[0] != b[0]); bad_q = sum(1 for a, b in zip(theirs[:k], mine[:k]) if a[1] != b[1])
        if bad_seq == 0 and bad_q == 0:
            return 'identical-to-ENA', f'{k} reads compared'
        return 'MISMATCH', f'{bad_seq} sequences and {bad_q} quality strings differ in the first {k} reads'
    alpha = quality_alphabet(fq1_path)
    if len(alpha) <= 2:
        return 'simplified-qualities-suspected', f'only {len(alpha)} distinct quality symbol(s): SRA Lite files have a constant or two-level quality'
    return 'skipped (ENA unreachable)', f'{len(alpha)} distinct quality symbols'


def download_sra(run_id, dest_dir, connections=8, stall_s=120.0, timeout=None, monitor=run_monitored):
    path = os.path.join(dest_dir, run_id + '.sra'); logp = path + '.aria2.log'
    os.makedirs(dest_dir, exist_ok=True)
    t0 = time.time()
    watch = lambda: os.path.getsize(path) if os.path.exists(path) else 0
    rc, err, stalled = monitor(aria2_command(sra_url(run_id), dest_dir, run_id + '.sra', connections), watch, stall_s=stall_s, timeout=timeout)
    log = ''
    try:
        log = open(logp, errors='replace').read()[-600:]
    except OSError:
        pass
    if stalled:
        raise EnaError(f'SRA download stalled for {stall_s:.0f} s', 'throttle')
    if rc == 3 or re.search(r'\b404\b|not found|NoSuchKey', err + log, re.I):
        raise EnaError(f'{sra_url(run_id)} does not exist (aria2c exit {rc})', 'source_unavailable')
    if rc != 0 or os.path.exists(path + '.aria2') or not os.path.exists(path):
        raise EnaError(f'aria2c exit {rc}: {(err + log)[-300:].strip()}', 'transient')
    secs = max(time.time() - t0, 1e-6); size = os.path.getsize(path)
    return path, dict(source='SRA-aws', method='aria2', seconds=round(secs, 1), bytes=size, mb_per_s=round(size / 1e6 / secs, 1))


def fetch_sra(spec, dest_dir, run=subprocess.run, monitor=run_monitored, qcheck=quality_check, log=print):
    """spec: run, reads_expected (ENA read_count), ena_url_r1, instrument, fq_threads, sra_connections, stall_s. Returns (paths, stats). Raises EnaError with a kind."""
    rid = spec['run']
    sra_path, stats = download_sra(rid, dest_dir, spec.get('sra_connections', 8), spec.get('stall_s', 120.0), spec.get('ascp_timeout'), monitor)
    tmpdir = os.path.join(dest_dir, rid + '.fqtmp'); os.makedirs(tmpdir, exist_ok=True)
    created = [sra_path, sra_path + '.aria2.log', os.path.join(dest_dir, f'{rid}_1.fastq'), os.path.join(dest_dir, f'{rid}_2.fastq'), os.path.join(dest_dir, f'{rid}.fastq')]

    def cleanup(keep_fastq=False):
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)
        for f in created[:2] + ([] if keep_fastq else created[2:]):
            try:
                os.remove(f)
            except OSError:
                pass
    try:
        v = run(validate_command(sra_path), capture_output=True, text=True)
        if v.returncode != 0:
            raise EnaError(f'vdb-validate failed: {(v.stdout + v.stderr)[-300:].strip()}', 'transient')
        t0 = time.time()
        r = run(fasterq_command(sra_path, dest_dir, tmpdir, spec.get('fq_threads', 8)), capture_output=True, text=True)
        if r.returncode != 0:
            raise EnaError(f'fasterq-dump exit {r.returncode}: {(r.stderr or "")[-300:].strip()}', 'transient')
        counts = parse_counts((r.stderr or '') + (r.stdout or ''))
        f1, f2 = os.path.join(dest_dir, f'{rid}_1.fastq'), os.path.join(dest_dir, f'{rid}_2.fastq')
        if not (os.path.exists(f1) and os.path.exists(f2)):
            raise EnaError(f'fasterq-dump wrote no _1/_2 pair for {rid} (single-end or technical reads only)', 'permanent')
        exp = spec.get('reads_expected')
        if exp and counts.get('spots') is not None and counts['spots'] != int(exp):
            raise EnaError(f"read count differs: {counts['spots']} spots from SRA, {exp} reads on ENA", 'source_unavailable')
        verdict, detail = qcheck(f1, spec.get('ena_url_r1', ''), spec.get('instrument', ''))
        if verdict in ('MISMATCH', 'simplified-qualities-suspected'):
            raise EnaError(f'{verdict}: {detail} (SRA Lite?)', 'source_unavailable')
        if verdict.startswith('skipped') and spec.get('require_quality_check'):
            raise EnaError(f'quality check against ENA not possible ({detail}) and --sra-require-check is set', 'transient')
        stats.update(quality_check=verdict, quality_detail=detail, reads_fastq=counts.get('spots', ''), convert_seconds=round(time.time() - t0, 1))
        cleanup(keep_fastq=True)
        return [f1, f2], stats
    except BaseException:
        cleanup()
        raise
