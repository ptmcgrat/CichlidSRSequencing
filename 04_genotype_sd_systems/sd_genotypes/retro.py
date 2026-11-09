"""Retrospective measurement: samples whose BAM is not on local disk are measured by streaming ONLY the regions the modules need from the BAM on Dropbox (a read-only rclone web server on the Bamfiles folder and one
`samtools view` per sample, as sexloci does), then running the panel on that small BAM. Writes <sample>.panel.json exactly as `panel measure` does, so the same call steps read them. Modules that read the whole
BAM (bchrom) or megabases (lg5inv) have no regions(cfg) and are refused: they are measured on disk.

  python -m sd_genotypes.retro --remote 'ptm_dropbox:.../Bamfiles/Mzebra_GT3_NCBI/' --samples S1,S2 (or --samples-file F) --modules gsdf,te --out-dir DIR [--jobs 4] [--keep-bams] [--refresh]

As a library: measure_sample(sample, url, modules, workdir) for one sample (url may also be a local path), run_samples(...) for a list (resumable: existing JSONs are kept unless refresh)."""
import argparse, importlib, os, re, socket, subprocess, threading, time, urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from sd_genotypes.panel import identity as ID, measure as pm


class DiskOnly(Exception):
    pass


class FetchError(Exception):
    def __init__(self, kind, message):
        super().__init__(f'{kind}: {message}')
        self.kind, self.message = kind, message


def classify_failure(err):
    e = err.lower()
    if '404' in e or 'not found' in e or 'no such file' in e:
        return 'no_bam'
    if 'index' in e or '.bai' in e or '.csi' in e:
        return 'no_index'
    return 'error'


def merge_region_strings(regions):
    """['c:1-10', 'c:5-20', 'd:1-3'] -> ['c:1-20', 'd:1-3'] (1-based, inclusive). Overlapping regions are merged because samtools would return their reads twice."""
    by = {}
    for r in regions:
        m = re.fullmatch(r'(.+):(\d+)-(\d+)', r)
        if not m:
            raise ValueError(f'not a region: {r!r}')
        by.setdefault(m.group(1), []).append((int(m.group(2)), int(m.group(3))))
    out = []
    for c, iv in by.items():
        iv.sort(); cur = list(iv[0])
        for a, b in iv[1:]:
            if a <= cur[1] + 1:
                cur[1] = max(cur[1], b)
            else:
                out.append(f'{c}:{cur[0]}-{cur[1]}'); cur = [a, b]
        out.append(f'{c}:{cur[0]}-{cur[1]}')
    return out


def module_regions(modules, cfgs=None):
    regs = []
    for m in modules:
        if m not in pm.MODULES:
            raise ValueError(f'unknown panel module {m!r}; available: {sorted(pm.MODULES)}')
        mod = importlib.import_module(pm.MODULES[m][0])
        if not hasattr(mod, 'regions'):
            raise DiskOnly(f'module {m!r} reads the whole BAM (or megabases of it) and cannot be measured retrospectively: measure it while the BAM is on disk')
        regs += mod.regions((cfgs or {}).get(m) or pm.module_cfg(m))
    return merge_region_strings(regs)


class HttpServer:
    """A read-only rclone web server on the cloud Bamfiles folder (the only network step)."""

    def __init__(self, remote_path, port=None):
        self.remote_path, self.port, self.proc = remote_path, port, None

    def __enter__(self):
        if self.port is None:
            with socket.socket() as s:
                s.bind(('127.0.0.1', 0))
                self.port = s.getsockname()[1]
        try:
            self.proc = subprocess.Popen(['rclone', 'serve', 'http', '--read-only', '--addr', f'127.0.0.1:{self.port}', self.remote_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except FileNotFoundError:
            raise RuntimeError('rclone is not installed or not on PATH') from None
        for _ in range(100):
            if self.proc.poll() is not None:
                raise RuntimeError(f"rclone serve exited immediately; check that {self.remote_path} exists (rclone lsf '{self.remote_path}')")
            with socket.socket() as s:
                if s.connect_ex(('127.0.0.1', self.port)) == 0:
                    return f'http://127.0.0.1:{self.port}'
            time.sleep(0.3)
        raise RuntimeError('rclone serve did not start within 30 s')

    def __exit__(self, *exc):
        if self.proc is not None:
            self.proc.terminate()
            try:
                self.proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.proc.kill()


def fetch_region_bam(url, regions, workdir, sample, run, retries=2, sleep=time.sleep):
    """samtools view -b of the regions from a URL or path into a sorted, indexed BAM in workdir. Returns (bam, [intermediate files])."""
    raw, bam = os.path.join(workdir, f'{sample}.raw.bam'), os.path.join(workdir, f'{sample}.regions.bam')
    err = ''
    for attempt in range(retries + 1):
        try:
            run(['samtools', 'view', '-b', '-o', raw, url] + regions)
            break
        except Exception as e:
            err = str(e); kind = classify_failure(err)
            if kind != 'error' or attempt == retries:
                raise FetchError(kind, (err.strip().splitlines() or [''])[-1][:160]) from None
            sleep(2 * (attempt + 1))
    for cmd in (['samtools', 'sort', '-o', bam, raw], ['samtools', 'index', bam]):
        try:
            run(cmd)
        except Exception as e:
            raise FetchError('error', f'{cmd[1]} failed: {str(e).strip()[-120:]}') from None
    return bam, [raw, bam + '.bai']


def measure_sample(sample, url, modules, workdir, run=None, cfgs=None, keep_bam=False, retries=2, sleep=time.sleep, head=None):
    """The panel dict of one sample measured from the regions streamed from `url`. Its `bam` and `bam_size_bytes` are those of the SOURCE BAM at `url` (read first, with an HTTP HEAD request; `head(url)`
    returns a size and is replaced in tests), not of the temporary regional BAM that is measured. Raises FetchError (no_bam, no_index, error) or DiskOnly."""
    from cichlid_common.common_light import sh
    run = run or sh
    regions = module_regions(modules, cfgs)
    try:
        ident = ID.source_identity(url, head or ID.http_size)
    except urllib.error.HTTPError as e:
        raise FetchError('no_bam' if e.code == 404 else 'error', f'cannot read the size of {url}: HTTP {e.code}') from None
    except Exception as e:
        raise FetchError('error', f'cannot read the size of {url}: {str(e)[:120]}') from None
    os.makedirs(workdir, exist_ok=True)
    bam, tmp = fetch_region_bam(url, regions, workdir, sample, run, retries, sleep)
    try:
        res = pm.measure(bam, sample, modules, run=run, cfgs=cfgs, identity=ident)
        res['source'] = url; res['retro_regions'] = regions
        return res
    finally:
        for p in tmp + ([] if keep_bam else [bam]):
            if os.path.exists(p):
                os.remove(p)


def run_samples(samples, remote, modules, out_dir, jobs=4, workdir=None, run=None, cfgs=None, keep_bams=False, refresh=False, server_factory=None, log=print, sleep=time.sleep, head=None):
    """Measure every sample not yet in out_dir (or all with refresh) and write <sample>.panel.json; failures go to retro_failures.tsv and never stop the others. Returns {'done', 'skipped', 'failed'}."""
    bad = [x for x in samples if not x or re.search(r'\s', x)]
    if bad:
        raise ValueError(f'sample names must not contain whitespace (tabs turned into spaces when pasting?): {bad[:3]}')
    os.makedirs(out_dir, exist_ok=True); workdir = workdir or os.path.join(out_dir, 'region_bams')
    module_regions(modules, cfgs)                                                         # refuse disk-only modules before any network step
    todo = [s for s in samples if refresh or not os.path.exists(os.path.join(out_dir, f'{s}.panel.json'))]
    skipped = [s for s in samples if s not in todo]; done, failed = [], {}; lock = threading.Lock()
    log(f'[retro] {len(samples)} samples: {len(skipped)} already measured, {len(todo)} to fetch')
    if todo:
        factory = server_factory or (lambda: HttpServer(remote))
        with factory() as base:
            def job(s):
                return measure_sample(s, f'{base}/{s}/{s}.all.bam', modules, workdir, run, cfgs, keep_bams, sleep=sleep, head=head)
            with ThreadPoolExecutor(max_workers=max(1, jobs)) as ex:
                futs = {ex.submit(job, s): s for s in todo}
                for i, f in enumerate(as_completed(futs), 1):
                    s = futs[f]
                    try:
                        pm.write_json(os.path.join(out_dir, f'{s}.panel.json'), f.result()); done.append(s); log(f'[retro] [{i}/{len(todo)}] {s}: measured')
                    except Exception as e:
                        with lock:
                            failed[s] = str(e)[:300]
                        log(f'[retro] [{i}/{len(todo)}] {s}: FAILED {str(e)[:160]}')
    if failed:
        with open(os.path.join(out_dir, 'retro_failures.tsv'), 'w') as fh:
            fh.write('sample\tmessage\n' + ''.join(f'{k}\t{v.replace(chr(9), " ").replace(chr(10), " ")}\n' for k, v in failed.items()))
    return dict(done=done, skipped=skipped, failed=failed)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.retro', description=__doc__.split('\n')[0])
    ap.add_argument('--remote', required=True, help='rclone path of the Bamfiles folder, e.g. ptm_dropbox:.../Bamfiles/Mzebra_GT3_NCBI/ (holding <sample>/<sample>.all.bam)')
    g = ap.add_mutually_exclusive_group(required=True); g.add_argument('--samples'); g.add_argument('--samples-file')
    ap.add_argument('--modules', default='gsdf,te'); ap.add_argument('--out-dir', required=True); ap.add_argument('--jobs', type=int, default=4); ap.add_argument('--keep-bams', action='store_true'); ap.add_argument('--refresh', action='store_true')
    a = ap.parse_args(argv)
    samples = a.samples.split(',') if a.samples else [x.strip() for x in open(a.samples_file) if x.strip() and not x.startswith('#')]
    r = run_samples(samples, a.remote, a.modules.split(','), a.out_dir, a.jobs, keep_bams=a.keep_bams, refresh=a.refresh)
    print(f"[retro] measured {len(r['done'])}, skipped {len(r['skipped'])}, failed {len(r['failed'])}")
    raise SystemExit(1 if r['failed'] else 0)


if __name__ == '__main__':
    main()
