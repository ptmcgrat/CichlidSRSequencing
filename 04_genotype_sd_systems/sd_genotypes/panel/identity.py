"""Identity of the BAM a measurement was made on: its file name and its size in bytes, stored in every <sample>.panel.json (schema 2) as `bam` and `bam_size_bytes`.
A realigned sample usually keeps its file name (<sample>.all.bam) but not its size, so a pipeline tells a measurement of the old BAM from one of the new BAM by comparing both.
For a sample measured through `retro` the identity is that of the SOURCE BAM on Dropbox (the regional BAM that is measured is temporary), read with an HTTP HEAD request.

  python -m sd_genotypes.panel.identity --json S.panel.json --bam S.all.bam
      exit 0: the measurement was made on this BAM (same name and size); 1: it was made on another BAM (realigned? re-measure); 2: the JSON does not record a size (made before the identity was stored)."""
import argparse, json, os, sys, urllib.parse, urllib.request


def bam_identity(path):
    """{'bam': file name, 'bam_size_bytes': size} of a local BAM."""
    return dict(bam=os.path.basename(path), bam_size_bytes=os.path.getsize(path))


def http_size(url, timeout=60):
    """Size in bytes of the file at a URL: HTTP HEAD (Content-Length), or a one-byte range request (Content-Range) if the server does not answer HEAD with a length."""
    with urllib.request.urlopen(urllib.request.Request(url, method='HEAD'), timeout=timeout) as r:
        n = r.headers.get('Content-Length')
    if n is not None:
        return int(n)
    req = urllib.request.Request(url, headers={'Range': 'bytes=0-0'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        cr = r.headers.get('Content-Range', '')
    if '/' in cr and cr.rsplit('/', 1)[1].isdigit():
        return int(cr.rsplit('/', 1)[1])
    raise ValueError(f'{url}: the server reports no size (no Content-Length, no Content-Range)')


def source_identity(url, head=http_size):
    """Identity of the file a URL (or a local path) points to. `head(url)` returns the size of a URL (replaced in tests)."""
    if os.path.exists(url):
        return bam_identity(url)
    return dict(bam=urllib.parse.unquote(urllib.parse.urlsplit(url).path.rsplit('/', 1)[-1]), bam_size_bytes=int(head(url)))


def normalize(d):
    """(bam file name, size in bytes or None) of a parsed panel JSON, schema 2 or the older schema 1. In schema 1 `bam` was an absolute path and the size was `bam_size`; for a streamed sample
    (it has `source`) they were those of a temporary regional BAM, so the name is taken from the source URL and the size is unknown."""
    if 'bam_size_bytes' in d:
        return d['bam'], d['bam_size_bytes']
    if 'source' in d:
        return os.path.basename(urllib.parse.urlsplit(d['source']).path), None
    return os.path.basename(d.get('bam', '')), d.get('bam_size')


def same_bam(panel, bam):
    """('same' | 'different' | 'unknown', reason). `panel`: a parsed panel JSON or its path; `bam`: a local path or an identity dict."""
    d = json.load(open(panel)) if isinstance(panel, str) else panel
    name, size = normalize(d)
    cur = bam if isinstance(bam, dict) else bam_identity(bam)
    if size is None:
        return 'unknown', f'the measurement does not record the size of {name} (made before the BAM identity was stored)'
    if name != cur['bam']:
        return 'different', f"measured on {name}, the BAM is {cur['bam']}"
    if size != cur['bam_size_bytes']:
        return 'different', f"{name} was {size} bytes when it was measured and is {cur['bam_size_bytes']} now (realigned?)"
    return 'same', ''


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.panel.identity', description=__doc__.split('\n')[0])
    ap.add_argument('--json', required=True); ap.add_argument('--bam', required=True)
    a = ap.parse_args(argv)
    status, why = same_bam(a.json, a.bam)
    print({'same': 'CURRENT: measured on this BAM', 'different': 'STALE: ' + why, 'unknown': 'UNKNOWN: ' + why}[status])
    raise SystemExit({'same': 0, 'different': 1, 'unknown': 2}[status])


if __name__ == '__main__':
    main()
