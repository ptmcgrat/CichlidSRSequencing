"""ENA (the EBI mirror of NCBI SRA) metadata and resumable, md5-verified downloads. Pure helpers: every external call goes through an injectable function, so the logic is unit tested
without a network."""
import hashlib, os, subprocess, time

FIELDS = ['run_accession', 'fastq_ftp', 'fastq_md5', 'fastq_bytes', 'base_count', 'read_count', 'library_layout', 'sample_alias']


class DownloadError(RuntimeError):
    pass


def filereport_url(accession, fields=FIELDS):
    return f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields={','.join(fields)}&format=tsv"


def parse_filereport(text):
    """ENA filereport TSV -> {run_accession: {run, alias, layout, base_count, read_count, files: [{name, url, md5, bytes}]}}.
    fastq_ftp / fastq_md5 / fastq_bytes are ';'-separated lists that must line up."""
    lines = [l for l in text.splitlines() if l.strip()]
    if not lines:
        return {}
    header = lines[0].split('\t')
    if 'run_accession' not in header:
        raise ValueError('not an ENA filereport (no run_accession column)')
    out = {}
    for l in lines[1:]:
        cells = l.split('\t')
        row = dict(zip(header, cells + [''] * (len(header) - len(cells))))
        urls = [u for u in row.get('fastq_ftp', '').split(';') if u]
        md5s = [m for m in row.get('fastq_md5', '').split(';') if m]
        sizes = [s for s in row.get('fastq_bytes', '').split(';') if s]
        if not (len(urls) == len(md5s) == len(sizes)):
            raise ValueError(f"{row['run_accession']}: fastq_ftp/fastq_md5/fastq_bytes lists differ in length ({len(urls)}/{len(md5s)}/{len(sizes)})")
        files = [dict(name=os.path.basename(u), url=u, md5=m, bytes=int(s)) for u, m, s in zip(urls, md5s, sizes)]
        out[row['run_accession']] = dict(run=row['run_accession'], alias=row.get('sample_alias', ''), layout=row.get('library_layout', ''),
                                         base_count=int(row['base_count']) if row.get('base_count') else None,
                                         read_count=int(row['read_count']) if row.get('read_count') else None, files=files)
    return out


def paired_files(rec):
    """(read 1, read 2) file records of a paired run. A third file without _1/_2 (reads whose mate was lost) is ignored and reported by the caller via `unpaired_files`."""
    r1 = [f for f in rec['files'] if f['name'].endswith('_1.fastq.gz')]
    r2 = [f for f in rec['files'] if f['name'].endswith('_2.fastq.gz')]
    if len(r1) != 1 or len(r2) != 1:
        raise ValueError(f"{rec['run']}: expected one _1 and one _2 FASTQ, found {[f['name'] for f in rec['files']]}")
    return r1[0], r2[0]


def unpaired_files(rec):
    return [f for f in rec['files'] if not (f['name'].endswith('_1.fastq.gz') or f['name'].endswith('_2.fastq.gz'))]


def full_url(u, scheme='https'):
    """ENA lists paths such as ftp.sra.ebi.ac.uk/vol1/fastq/...; they are served over https as well as ftp."""
    return u if '://' in u else f'{scheme}://{u}'


def md5_file(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def curl_command(url, dest):
    """Resumable download (-C -) into `dest`; redirects followed; transient errors retried by curl itself, and again by download_file."""
    return ['curl', '-fL', '--retry', '3', '--retry-delay', '10', '-C', '-', '-o', dest, url]


def _complete(path, rec, md5=md5_file):
    return os.path.exists(path) and os.path.getsize(path) == rec['bytes'] and md5(path) == rec['md5']


def download_file(rec, dest_dir, scheme='https', attempts=4, backoff=30.0, run=subprocess.run, sleep=time.sleep, md5=md5_file, log=print):
    """Download one FASTQ to dest_dir/<name>, verifying size and md5. Resumes a partial file; a file that is complete but wrong is deleted and fetched again.
    Returns the final path. Raises DownloadError after `attempts` failures."""
    os.makedirs(dest_dir, exist_ok=True)
    final = os.path.join(dest_dir, rec['name']); part = final + '.partial'
    if _complete(final, rec, md5):
        return final
    last = ''
    for k in range(1, attempts + 1):
        r = run(curl_command(full_url(rec['url'], scheme), part), capture_output=True, text=True)
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have == rec['bytes']:                                  # also covers "curl complained that there was nothing left to resume"
            if md5(part) == rec['md5']:
                os.replace(part, final)
                return final
            last = f"md5 mismatch after {have} bytes"
            os.remove(part)
        elif have > rec['bytes']:
            last = f'file larger than expected ({have} > {rec["bytes"]})'
            os.remove(part)
        else:
            last = f'curl exit {r.returncode}: {(r.stderr or "")[-200:].strip()} ({have}/{rec["bytes"]} bytes)'
        log(f"    download {rec['name']}: attempt {k}/{attempts} failed: {last}")
        if k < attempts:
            sleep(backoff * k)
    raise DownloadError(f"{rec['name']}: {last}")


def fetch_filereport(accession, run=subprocess.run):
    r = run(['curl', '-fsSL', '--retry', '3', filereport_url(accession)], capture_output=True, text=True)
    if r.returncode:
        raise DownloadError(f'ENA filereport for {accession} failed: curl exit {r.returncode} {(r.stderr or "")[-200:]}')
    return r.stdout
