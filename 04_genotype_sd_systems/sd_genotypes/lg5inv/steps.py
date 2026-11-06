"""lg5inv steps: prepare (per pool: download -> align -> keep the LG5 slices), status, clean. Reuses bchrom's tested pieces; own config (config/lg5inv.yaml) and work folder."""
import json, os, shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import yaml
from cichlid_common.common_light import config_path
from sd_genotypes.bchrom import ena, steps as BS, workflow as BW
from . import workflow as LW

MARKER = '.lg5inv_work'


def cfg_of(path=None):
    return yaml.safe_load(open(path or config_path('lg5inv.yaml')))


def validate_config(cfg):
    if cfg['download'].get('source', 'ena') not in ('ena', 'sra', 'auto'):
        raise SystemExit(f"download.source must be ena, sra or auto, not {cfg['download'].get('source')!r}")
    ids = [r['id'] for r in cfg['runs']]
    if len(set(ids)) != len(ids):
        raise ValueError('duplicate run ids in config/lg5inv.yaml')
    for r in cfg['runs']:
        if r.get('role') != 'pool':
            raise ValueError(f"{r['id']}: role must be pool")
    sl = cfg.get('slices') or {}
    if not (sl.get('enabled') and sl.get('regions')):
        raise ValueError('lg5inv needs slices.enabled and at least one slices.regions entry: the slices are the point of this analysis')
    return {r['id']: r for r in cfg['runs']}


def work_dir(proj, cfg):
    """A new or empty folder is stamped with a marker; a non-empty folder without it is not ours and is refused (clean deletes it)."""
    w = cfg.get('work_dir') or (proj.dirs['tmp'] + 'lg5inv/')
    w = w if w.endswith('/') else w + '/'
    if os.path.isdir(w) and os.listdir(w) and not os.path.exists(w + MARKER):
        raise SystemExit(f'work folder {w} exists, is not empty and was not created by lg5inv (no {MARKER} file); choose another work_dir in config/lg5inv.yaml')
    os.makedirs(w + 'depth/', exist_ok=True)
    open(w + MARKER, 'a').close()
    return w


def recorded(proj, rid):
    return os.path.exists(f'{proj.root}lg5inv/pool_{rid}/LATEST')


def _record_run(proj, cfg, work, rec, spec, summ):
    """One immutable run folder per pool, written as soon as the pool is finished: summary, depth track (a few tens of MB), the slice BAMs."""
    rid = spec['id']
    row = dict(run=rid, label=spec.get('label'), population=spec.get('population'), species=spec.get('species'), sex=spec.get('sex'), morph=spec.get('morph'), n_fish=spec.get('n_fish'),
               expected_inverted_fraction=spec.get('expected_inverted_fraction'), bases=summ.get('bases'), observed_mean_depth=summ.get('observed_mean_depth'),
               expected_mean_depth=summ.get('expected_mean_depth'), depth_ratio=summ.get('depth_ratio'), median_depth=summ.get('median_depth'), zero_fraction=summ.get('zero_fraction'),
               percent_duplication=(summ.get('dup_metrics') or {}).get('PERCENT_DUPLICATION'), flags=' | '.join(summ.get('flags') or []))
    with proj.run('lg5inv', f'pool_{rid}', params=dict(cfg=cfg, run=rid)) as run:
        pd.DataFrame([row]).to_csv(run.dir + 'run.tsv', sep='\t', index=False)
        shutil.copy2(f'{work}depth/{rid}.summary.json', run.dir + 'summary.json')
        if os.path.exists(f'{work}ena/{rid}.tsv'):
            shutil.copy2(f'{work}ena/{rid}.tsv', run.dir + 'ena.tsv')
        os.makedirs(run.dir + 'depth/', exist_ok=True)
        for ext in ('.npz', '.json'):
            shutil.copy2(f'{work}depth/{rid}{ext}', run.dir + 'depth/' + rid + ext)
        if (cfg.get('slices') or {}).get('upload') and os.path.isdir(work + 'slices/'):
            os.makedirs(run.dir + 'slices/', exist_ok=True)
            for f in sorted(os.listdir(work + 'slices/')):
                if f.startswith(rid + '.'):
                    shutil.copy2(work + 'slices/' + f, run.dir + 'slices/' + f)
            json.dump(dict(zip([r['name'] for r in cfg['slices']['regions']], BS.slice_regions(proj, cfg))), open(run.dir + 'slices/regions.json', 'w'), indent=1)


def prepare_step(proj, only=None, keep_bams=False, runs_in_flight=None, threads=None, executor=None, cmds=None, fetch=None, downloader=None, cfg=None, sra_fetcher=None):
    cfg = cfg or cfg_of(); validate_config(cfg)
    if threads:
        cfg['resources']['align_threads'] = int(threads)
    work = work_dir(proj, cfg)
    specs = cfg['runs']
    if only:
        unknown = set(only) - {s['id'] for s in specs}
        if unknown:
            raise SystemExit(f'unknown run ids {sorted(unknown)}; see config/lg5inv.yaml')
        specs = [s for s in specs if s['id'] in only]
    recs = BS.get_records(work, specs, fetch or ena.fetch_filereport)
    g = BS.reference_info(proj, cfg)
    if executor is None:
        problems = BS.preflight(g, cfg=cfg)
        if problems:
            raise SystemExit('lg5inv prepare cannot start:\n  ' + '\n  '.join(problems))
    ledger = BW.SpaceLedger(work, cfg['resources']['min_free_gb'])
    cmds = cmds or BW.LabCommands()
    n = int(runs_in_flight or cfg['resources']['runs_in_flight'])
    total_gb = sum(sum(f['bytes'] for f in recs[s['id']]['files']) for s in specs) / 1e9
    print(f'[lg5inv] {len(specs)} runs, {total_gb:.0f} GB of FASTQ, {n} in flight, work folder {work}; slices: {BS.slice_regions(proj, cfg)} + unmapped reads: {bool(cfg["slices"].get("unmapped"))}')
    results, failures = {}, {}

    def job(spec):
        return LW.Lg5RunProcessor(spec, recs[spec['id']], work, cfg, g, executor=executor, cmds=cmds, ledger=ledger, downloader=downloader, keep_bams=keep_bams, sra_fetcher=sra_fetcher).process()

    with ThreadPoolExecutor(max_workers=n) as ex:
        futs = {ex.submit(job, s): s for s in specs}
        for f in as_completed(futs):
            s = futs[f]
            try:
                results[s['id']] = f.result()
                print(f"[lg5inv] {s['id']} {s.get('species', '')} {s.get('morph', '')}: done")
            except Exception as e:                                                      # one failed run must not stop the others
                failures[s['id']] = str(e)[-600:]
                print(f"[lg5inv] {s['id']} FAILED: {str(e)[-600:]}")
                continue
            try:
                _record_run(proj, cfg, work, recs[s['id']], s, results[s['id']])
                print(f"[lg5inv] {s['id']}: recorded (lg5inv/pool_{s['id']})")
            except Exception as e:
                failures[s['id']] = f'finished, but the run folder was not written: {str(e)[-500:]}; the slices are still in {work}slices/ (re-run prepare to record it)'
                print(f"[lg5inv] {s['id']} RECORD FAILED: {str(e)[-500:]}")
            if results[s['id']].get('flags'):
                print(f"[lg5inv] {s['id']} warnings: {' | '.join(results[s['id']]['flags'])}")
    if failures:
        print('[lg5inv] failed runs (re-run `prepare`; finished runs and finished stages are skipped):', sorted(failures))
        raise SystemExit(1)


def status_step(proj, cfg=None):
    cfg = cfg or cfg_of(); validate_config(cfg)
    work = work_dir(proj, cfg)
    print(f'work folder {work}   free {shutil.disk_usage(work).free / 1e9:.0f} GB   in use {BW.dir_size(work) / 1e9:.1f} GB')
    rows = []
    for s in cfg['runs']:
        rp = BW.RunPaths(work, s['id'])
        lp = LW.Lg5RunProcessor.__new__(LW.Lg5RunProcessor)                           # only for finished(); no processing object is needed
        lp.p, lp.g, lp.cfg, lp.id = rp, dict(slice_regions=['x']), cfg, s['id']
        lp.unmapped_bam = os.path.join(os.path.dirname(rp.slice_bam), s['id'] + '.unmapped.bam')
        done = lp.finished()
        stages = [st for st in BW.STAGES if os.path.exists(rp.stages + st + '.done')]
        rows.append(dict(run=s['id'], sample=f"{s.get('species', '')} {s.get('sex', '')} {s.get('morph', '')}"[:34], state='complete' if done else (stages[-1] + ' done' if stages else ('started' if os.path.isdir(rp.dir) else 'not started')),
                         on_dropbox='yes' if recorded(proj, s['id']) else '-', scratch_gb=round(BW.dir_size(rp.dir) / 1e9, 1) if os.path.isdir(rp.dir) else 0))
    t = pd.DataFrame(rows)
    print(t.to_string(index=False)); print(f"\n{(t.state == 'complete').sum()} of {len(t)} runs complete, {(t.on_dropbox == 'yes').sum()} on Dropbox")
    return t


def clean_step(proj, force=False, cfg=None):
    """Delete the work folder. Refuses while a finished run has no run folder on Dropbox (the slices in the work folder are then the only copy), unless --force."""
    cfg = cfg or cfg_of(); validate_config(cfg); work = work_dir(proj, cfg)
    if not force:
        missing = [s['id'] for s in cfg['runs'] if os.path.exists(BW.RunPaths(work, s['id']).slice_bam) and not recorded(proj, s['id'])]
        if missing:
            raise SystemExit(f'runs with slices but no lg5inv/pool_<run> run folder: {missing}; run `prepare` again to record them (or use --force)')
    if not os.path.exists(work + MARKER):
        raise SystemExit(f'{work} was not created by lg5inv (no {MARKER} file): refusing to delete it')
    n = BW.dir_size(work)
    shutil.rmtree(work, ignore_errors=True)
    print(f'[lg5inv] deleted {work} ({n / 1e9:.1f} GB)')
