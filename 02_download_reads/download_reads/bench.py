"""Find how many parallel Aspera downloads from ENA actually help, on YOUR network, before a big download campaign.

  python3 unit_scripts/ena_bench.py ERR715541 ERR13389714 ERR13389716 ERR13389717 --concurrency 1 2 3 4 --seconds 60 --rate 300m

For each concurrency c it starts c ascp transfers (one file each, different runs; give at least max(c) runs, preferably large ones so that no file finishes during the test), lets them run for
--seconds, stops them, and reports the aggregate MB/s and how many sessions died early (refused / aborted: a sign of throttling). The right --dl-workers is the smallest c after which the aggregate
stops growing (or the sessions start dying); the right --ascp-rate is about (link capacity / that c). The files are deleted afterwards. Needs the same ascp and key as the downloader."""
import argparse, os, shutil, subprocess, sys, tempfile, time
from download_reads.sources import ena as E


def dir_bytes(d):
    return sum(os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(d) for f in fs if os.path.exists(os.path.join(r, f)))


def bench(files, c, seconds, rate, ascp, key, popen=subprocess.Popen, sleep=time.sleep, now=time.time):
    d = tempfile.mkdtemp(prefix='ena_bench_')
    procs = [popen(E.ascp_command([f], d, ascp, key, rate), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE) for f in files[:c]]
    t0 = now(); sleep(seconds)
    died = [p for p in procs if p.poll() is not None and p.returncode != 0]
    errs = [(p.stderr.read() or b'').decode('utf-8', 'replace')[-120:] for p in died]
    for p in procs:
        if p.poll() is None:
            p.terminate()
    for p in procs:
        try:
            p.wait(timeout=10)
        except Exception:
            p.kill()
    mb = dir_bytes(d) / 1e6; secs = max(now() - t0, 1e-6)
    shutil.rmtree(d, ignore_errors=True)
    return dict(concurrency=c, aggregate_mb_s=round(mb / secs, 1), per_transfer_mb_s=round(mb / secs / c, 1), died_early=len(died), errors=errs)


def bench_https(files, c, seconds, popen=subprocess.Popen, sleep=time.sleep, now=time.time):
    """The same measurement over https: c curl downloads of c different files, each stopped by curl's own time limit."""
    d = tempfile.mkdtemp(prefix='ena_bench_')
    procs = [popen(['curl', '-sL', '-m', str(int(seconds)), '-o', os.path.join(d, f['name']), f['url'] if '://' in f['url'] else 'https://' + f['url']], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE) for f in files[:c]]
    t0 = now(); sleep(seconds + 2)
    for p in procs:
        if p.poll() is None:
            p.terminate()
    died = [p for p in procs if p.poll() not in (None, 0, 28)]                    # 28 = curl's own time limit, expected
    mb = dir_bytes(d) / 1e6; secs = max(now() - t0, 1e-6)
    shutil.rmtree(d, ignore_errors=True)
    return dict(concurrency=c, aggregate_mb_s=round(mb / secs, 1), per_transfer_mb_s=round(mb / secs / c, 1), died_early=len(died), errors=[])


def resolve(ascp=None, key=None, exists=os.path.exists):
    a = E.find_first(ascp, E.ASCP_CANDIDATES, which='ascp', env='ASCP', exists=exists); k = E.find_first(key, E.KEY_CANDIDATES, exists=exists)
    if not (a and k):
        raise SystemExit(f"ascp {'found at ' + a if a else 'NOT found'}, key {'found at ' + k if k else 'NOT found'}: give --ascp and --key (searched {E.ASCP_CANDIDATES} and {E.KEY_CANDIDATES})")
    return a, k


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m download_reads bench', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('runs', nargs='+'); ap.add_argument('--concurrency', nargs='+', type=int, default=[1, 2, 3, 4]); ap.add_argument('--seconds', type=float, default=60)
    ap.add_argument('--rate', default='300m'); ap.add_argument('--https', action='store_true', help='measure https instead of Aspera'); ap.add_argument('--ascp', default=None); ap.add_argument('--key', default=None)
    a = ap.parse_args(argv)
    a.ascp, a.key = resolve(a.ascp, a.key)
    print(f'ascp {a.ascp}, key {a.key}')
    files = []
    for r in a.runs:
        s = E.select_fastqs(E.fetch_runs(r).iloc[0])
        if s['status'] == 'ok':
            files.append(s['r1'])
    if len(files) < max(a.concurrency):
        raise SystemExit(f'need at least {max(a.concurrency)} downloadable runs, got {len(files)}')
    print(('https' if a.https else f'Aspera, rate per transfer {a.rate}') + f', {a.seconds:.0f} s per setting')
    rows = []
    for c in a.concurrency:
        r = bench_https(files, c, a.seconds) if a.https else bench(files, c, a.seconds, a.rate, a.ascp, a.key); rows.append(r)
        print(f"  {c} parallel: {r['aggregate_mb_s']:7.1f} MB/s aggregate, {r['per_transfer_mb_s']:6.1f} per transfer, {r['died_early']} sessions died early {r['errors'][:1]}")
        time.sleep(10)                                                    # let ENA release the sessions
    best = max(rows, key=lambda r: r['aggregate_mb_s'])
    if best['aggregate_mb_s'] < 1.0 and not a.https:
        print('-> Aspera moved almost no data at any setting: the TCP control connection works but the UDP data channel is probably blocked. Run again with --https, and ask IT to allow outbound UDP 33001 to fasp.sra.ebi.ac.uk.')
        return rows
    knee = next((r for r in rows if r['aggregate_mb_s'] >= 0.9 * best['aggregate_mb_s'] and r['died_early'] == 0), best)
    print(f"-> best aggregate at {best['concurrency']} parallel; the smallest setting within 10 % of it without dying sessions is {knee['concurrency']}: use --dl-workers {knee['concurrency']}")
    return rows


if __name__ == '__main__':
    main()
