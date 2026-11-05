"""Two-pool scheduler for the ENA downloads: a SMALL download pool that the server and the network tolerate, and a larger conversion pool for the CPU-heavy stages (FastqToSam, adapter marking, upload).

Why: one Aspera/FASP transfer is built to fill the link by itself, so more parallel downloads do not add speed, they compete and get throttled. The old script ran download + conversion in one worker,
so keeping the CPUs busy meant more simultaneous downloads. Here the download pool is capped (and shrinks by itself when the server pushes back: AIMD, like TCP), while converting runs as soon as
its FASTQs are on disk, on as many workers as the machine has. Backpressure (max_ready) and a disk ledger keep the FASTQs from piling up. Everything time-related is injectable for tests."""
import json, os, shutil, threading, time
from . import ena_fetch as E


class Governor:
    """AIMD limit on concurrent downloads: a throttle signal halves it (not below `lo`) and starts a cool-down; `up_after` clean downloads add one (not above `hi`)."""
    def __init__(self, lo=1, hi=3, start=None, up_after=3, cooldown=60.0, now=time.time):
        self.lo, self.hi, self.up_after, self.cooldown, self.now = lo, hi, up_after, cooldown, now
        self.allowed = min(hi, max(lo, start if start is not None else hi)); self.ok = 0; self.quiet_until = 0.0; self.lock = threading.Lock(); self.events = []

    def on_success(self):
        with self.lock:
            self.ok += 1
            if self.ok >= self.up_after and self.allowed < self.hi and self.now() >= self.quiet_until:
                self.allowed += 1; self.ok = 0; self.events.append(('up', self.allowed))

    def on_throttle(self):
        with self.lock:
            self.allowed = max(self.lo, self.allowed // 2); self.ok = 0; self.quiet_until = self.now() + self.cooldown; self.events.append(('down', self.allowed))

    def may_start(self, n_running):
        return n_running < self.allowed and self.now() >= self.quiet_until


def usage(folder, run):
    tot = 0
    try:
        for f in os.listdir(folder):
            if run in f:
                try:
                    tot += os.path.getsize(os.path.join(folder, f))
                except OSError:
                    pass
    except OSError:
        pass
    return tot


def read_status(path, log_path, rc):
    try:
        st = json.load(open(path))
    except (OSError, ValueError):
        tail = ''
        try:
            tail = open(log_path, errors='replace').read()[-600:]
        except OSError:
            pass
        return dict(state='failed', stage='worker crashed', kind='transient', error=f'no status file (exit {rc}): {tail}')
    if rc != 0 and st.get('state') in ('done', 'fetched'):
        st.update(state='failed', stage='exit status', kind='transient', error=f'worker reported {st.get("state")} but exited {rc}')
    return st


def run_tasks(tasks, cfg, launch, on_done, sleep=time.sleep, now=time.time, free_fn=None, log=print):
    """tasks: dicts with run, fastq_gb, and whatever the workers need. cfg: temp_dir, dl_max, workers, max_ready, stagger_s, max_attempts, backoff_s, min_free_gb, disk_factor, status_every_s.
    launch(kind, spec_path, log_path) -> process with .poll(); kind is 'fetch' or 'convert'. on_done(task, status) is called for a finished conversion that reported done (database row);
    it may raise, which turns the run into a failure. Returns {run: final status}."""
    free_fn = free_fn or (lambda p: shutil.disk_usage(p).free)
    tmp = cfg['temp_dir']; gov = Governor(lo=1, hi=cfg['dl_max'], now=now)
    queue = [dict(t, attempts=0, not_before=0.0) for t in tasks]
    fetching, ready, converting, results = [], [], [], {}
    last_start, last_status, bytes_done, t_begin = -1e9, now(), 0.0, now()

    def spec_for(t, kind):
        spec = dict(t['spec'], stage=kind, status_json=f"{tmp}{t['run']}.{kind}.status.json")
        path = f"{tmp}{t['run']}.{kind}.spec.json"
        json.dump(spec, open(path, 'w'), indent=1)
        if os.path.exists(spec['status_json']):
            os.remove(spec['status_json'])
        return path, spec['status_json']

    def retry_or_fail(t, st):
        t['attempts'] += 1
        if st.get('kind') != 'permanent' and t['attempts'] < cfg['max_attempts']:
            t['not_before'] = now() + cfg['backoff_s'] * t['attempts']; queue.append(t)
            log(f"  {t['run']}: {st['stage']} failed ({st.get('kind')}), retry {t['attempts']}/{cfg['max_attempts'] - 1} in {cfg['backoff_s'] * t['attempts']:.0f} s")
        else:
            results[t['run']] = dict(st, state='failed', attempts=t['attempts'])
            log(f"  {t['run']}: FAILED ({st['stage']}, {st.get('kind')}): {str(st.get('error', ''))[:150].replace(chr(10), ' ')}")

    while queue or fetching or ready or converting:
        for r in list(fetching):
            rc = r['proc'].poll()
            if rc is None:
                continue
            fetching.remove(r); st = read_status(r['status_json'], r['log'], rc)
            if st['state'] == 'fetched':
                gov.on_success(); bytes_done += r['task']['fastq_gb'] * 1e9; r['task']['dl'] = st.get('download'); ready.append(r['task'])
                d = st.get('download') or {}
                log(f"  {r['task']['run']}: downloaded {r['task']['fastq_gb']:.1f} GB via {d.get('method', '?')} at {d.get('mb_per_s', '?')} MB/s")
            else:
                if st.get('kind') == 'throttle':
                    gov.on_throttle(); log(f'  server pushing back: downloads limited to {gov.allowed}, cool-down {gov.cooldown:.0f} s')
                retry_or_fail(r['task'], st)
        for r in list(converting):
            rc = r['proc'].poll()
            if rc is None:
                continue
            converting.remove(r); st = read_status(r['status_json'], r['log'], rc)
            if st['state'] == 'done':
                try:
                    on_done(r['task'], st); results[r['task']['run']] = dict(st, attempts=r['task']['attempts'] + 1)
                    log(f"  {r['task']['run']}: done")
                except Exception as e:
                    results[r['task']['run']] = dict(st, state='failed', stage='database', kind='permanent', error=f'uBAM uploaded but the database row failed: {e}')
                    log(f"  {r['task']['run']}: uBAM uploaded but the database row failed: {e}")
            else:
                retry_or_fail(dict(r['task']), st) if st.get('kind') != 'permanent' else results.__setitem__(r['task']['run'], dict(st, state='failed'))
                if st.get('kind') == 'permanent':
                    log(f"  {r['task']['run']}: FAILED ({st['stage']}): {str(st.get('error', ''))[:150].replace(chr(10), ' ')}")
        while ready and len(converting) < cfg['workers']:                       # drain first: a run whose FASTQs are on disk converts before a new download starts
            t = ready.pop(0); path, sj = spec_for(t, 'convert')
            converting.append(dict(task=t, status_json=sj, log=f"{tmp}{t['run']}.convert.log", proc=launch('convert', path, f"{tmp}{t['run']}.convert.log")))
        while queue and gov.may_start(len(fetching)) and len(ready) < cfg['max_ready'] and now() - last_start >= cfg['stagger_s']:
            cand = next((t for t in queue if t['not_before'] <= now()), None)
            if cand is None:
                break
            need = cfg['disk_factor'] * cand['fastq_gb'] * 1e9
            inflight = [x['task'] for x in fetching + converting] + ready
            reserved = sum(max(0.0, cfg['disk_factor'] * t['fastq_gb'] * 1e9 - usage(tmp, t['run'])) for t in inflight)
            if free_fn(tmp) - reserved < need + cfg['min_free_gb'] * 1e9:
                if not inflight:
                    queue.remove(cand)
                    results[cand['run']] = dict(state='failed', stage='disk', kind='transient', error=f"not enough free disk for {need / 1e9:.0f} GB (free {free_fn(tmp) / 1e9:.0f} GB, margin {cfg['min_free_gb']:.0f} GB)", attempts=cand['attempts'])
                    continue
                break
            queue.remove(cand); path, sj = spec_for(cand, 'fetch'); last_start = now()
            log(f"Starting download of {cand['run']} ({cand['sample']}, {cand['fastq_gb']:.1f} GB); downloads {len(fetching) + 1}/{gov.allowed}")
            fetching.append(dict(task=cand, status_json=sj, log=f"{tmp}{cand['run']}.fetch.log", proc=launch('fetch', path, f"{tmp}{cand['run']}.fetch.log")))
        if now() - last_status >= cfg['status_every_s']:
            last_status = now(); el = max(now() - t_begin, 1e-6)
            log(f"[{time.strftime('%H:%M:%S')}] downloading {len(fetching)}/{gov.allowed}, ready {len(ready)}, converting {len(converting)}/{cfg['workers']}, waiting {len(queue)}, done {sum(1 for v in results.values() if v['state'] == 'done')}, "
                f"failed {sum(1 for v in results.values() if v['state'] == 'failed')}, average download {bytes_done / 1e6 / el:.0f} MB/s, free {free_fn(tmp) / 1e9:.0f} GB")
        if fetching or converting or queue:
            sleep(cfg.get('poll_s', 5))
    return results
