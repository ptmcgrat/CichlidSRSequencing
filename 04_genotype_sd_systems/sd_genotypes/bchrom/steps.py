"""bchrom steps: prepare (download -> align -> depth track), call (blocks, core, evaluation), status, clean, compare."""
import json, os, shutil, time
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd
import yaml
from cichlid_common.common_light import config_path
from . import ena, depth as D, calling as C, compare as K, workflow as W


def cfg_of(path=None):
    return yaml.safe_load(open(path or config_path('bchrom.yaml')))


def validate_config(cfg):
    if cfg['download'].get('source', 'ena') not in ('ena', 'sra', 'auto'):
        raise SystemExit(f"download.source must be ena, sra or auto, not {cfg['download'].get('source')!r}")
    runs = cfg['runs']; ids = [r['id'] for r in runs]
    if len(set(ids)) != len(ids):
        raise ValueError('duplicate run ids in config/bchrom.yaml')
    by = {r['id']: r for r in runs}
    for r in runs:
        if r['role'] not in ('carrier', 'baseline', 'control'):
            raise ValueError(f"{r['id']}: role must be carrier, baseline or control")
        if r['role'] in ('carrier', 'control'):
            b = by.get(r.get('baseline'))
            if b is None or b['role'] != 'baseline':
                raise ValueError(f"{r['id']}: baseline {r.get('baseline')} is not a run with role baseline")
    n = sum(r['role'] == 'carrier' for r in runs)
    if cfg['core']['min_individuals'] > n:
        raise ValueError(f"core.min_individuals {cfg['core']['min_individuals']} exceeds the {n} carriers")
    return by


REQUIRED_TOOLS = ['curl', 'gatk', 'minimap2', 'samtools', 'mosdepth', 'bash']
SRA_TOOLS = ['aria2c', 'vdb-validate', 'fasterq-dump']                       # needed when download.source is sra or auto


def preflight(g, which=shutil.which, exists=os.path.exists, cfg=None):
    """Everything the real run needs, checked BEFORE the first download (a missing tool would otherwise show up after 30 GB of data). Returns a list of problems."""
    tools = REQUIRED_TOOLS + (SRA_TOOLS if cfg and 'sra' in W.sources_of(cfg) else [])
    bad = [f'`{t}` is not on PATH' for t in tools if not which(t)]
    for k in ('reference', 'index'):
        if not (g.get(k) and exists(g[k])):
            bad.append(f'{k} file not found: {g.get(k)}')
    return bad


MARKER = '.bchrom_work'


def work_dir(proj, cfg):
    """The work folder. A new or empty folder is stamped with a marker file; a non-empty folder without the marker is not ours and is refused (clean deletes this folder)."""
    w = cfg.get('work_dir') or (proj.dirs['tmp'] + 'bchrom/')
    w = w if w.endswith('/') else w + '/'
    if os.path.isdir(w) and os.listdir(w) and not os.path.exists(w + MARKER):
        raise SystemExit(f'work folder {w} exists, is not empty and was not created by bchrom (no {MARKER} file); choose another bchrom.work_dir in config/bchrom.yaml')
    os.makedirs(w + 'depth/', exist_ok=True)
    open(w + MARKER, 'a').close()
    return w


def slice_regions(proj, cfg):
    """samtools region strings (1-based, inclusive) of config slices.regions; LG labels are resolved to the reference's sequence names."""
    sl = cfg.get('slices') or {}
    if not sl.get('enabled'):
        return []
    return [f"{proj.chrom_for_lg(r['lg'])}:{r['start'] + 1}-{r['end']}" for r in sl.get('regions', [])]


def reference_info(proj, cfg=None):
    idx = getattr(proj.fm, 'localMinimapGenomeFile', None)
    return dict(reference=proj.genome_fa, index=idx if idx and os.path.exists(idx) else proj.genome_fa, slice_regions=slice_regions(proj, cfg) if cfg else [])


def get_records(work, specs, fetch=ena.fetch_filereport):
    """ENA file lists (names, sizes, md5, base counts) of the runs, cached in <work>/ena/ so that an interrupted run does not ask again."""
    os.makedirs(work + 'ena/', exist_ok=True); recs = {}
    for s in specs:
        f = f"{work}ena/{s['id']}.tsv"
        if not os.path.exists(f):
            text = fetch(s['id'])
            open(f + '.tmp', 'w').write(text); os.replace(f + '.tmp', f)
        parsed = ena.parse_filereport(open(f).read())
        if s['id'] not in parsed:
            os.remove(f)
            raise ValueError(f"ENA returned no record for {s['id']}")
        recs[s['id']] = parsed[s['id']]
    return recs


# =================================================================================================== prepare
def prepare_step(proj, only=None, keep_bams=False, runs_in_flight=None, threads=None, executor=None, cmds=None, fetch=None, downloader=None, cfg=None, sra_fetcher=None):
    cfg = cfg or cfg_of(); validate_config(cfg)
    if threads:
        cfg['resources']['align_threads'] = int(threads)
    work = work_dir(proj, cfg)
    specs = cfg['runs']
    if only:
        unknown = set(only) - {s['id'] for s in specs}
        if unknown:
            raise SystemExit(f'unknown run ids {sorted(unknown)}; see config/bchrom.yaml')
        specs = [s for s in specs if s['id'] in only]
    recs = get_records(work, specs, fetch or ena.fetch_filereport)
    g = reference_info(proj, cfg)
    if executor is None:                                                                 # real run (tests pass a fake executor)
        problems = preflight(g, cfg=cfg)
        if problems:
            raise SystemExit('bchrom prepare cannot start:\n  ' + '\n  '.join(problems))
    ledger = W.SpaceLedger(work, cfg['resources']['min_free_gb'])
    cmds = cmds or W.LabCommands()
    n = int(runs_in_flight or cfg['resources']['runs_in_flight'])
    total_gb = sum(sum(f['bytes'] for f in recs[s['id']]['files']) for s in specs) / 1e9
    print(f'[bchrom] {len(specs)} runs, {total_gb:.0f} GB of FASTQ, {n} in flight, work folder {work}')
    order = sorted(specs, key=lambda s: -(s.get('bases_gb') or 0))
    results, failures = {}, {}

    def job(spec):
        return W.RunProcessor(spec, recs[spec['id']], work, cfg, g, executor=executor, cmds=cmds, ledger=ledger, downloader=downloader, keep_bams=keep_bams, sra_fetcher=sra_fetcher).process()

    with ThreadPoolExecutor(max_workers=n) as ex:
        futs = {ex.submit(job, s): s for s in order}
        for f in as_completed(futs):
            s = futs[f]
            try:
                results[s['id']] = f.result()
                print(f"[bchrom] {s['id']} {s.get('label', '')}: done")
            except Exception as e:                                                       # one failed run must not stop the others
                failures[s['id']] = str(e)[-600:]
                print(f"[bchrom] {s['id']} FAILED: {str(e)[-600:]}")
    if results:
        _record_prepare(proj, cfg, work, recs, results, failures)
    if failures:
        print('[bchrom] failed runs (re-run `bchrom prepare`; finished runs and finished stages are skipped):', sorted(failures))
        raise SystemExit(1)


def g_regions(cfg, proj):
    return dict(zip([r['name'] for r in cfg['slices']['regions']], slice_regions(proj, cfg)))


def _record_prepare(proj, cfg, work, recs, results, failures):
    """Provenance of the data that was used (the reads themselves are deleted): accessions, file md5s, base counts, depth checks, flags."""
    rows = []
    for rid, s in results.items():
        dm = s.get('dup_metrics') or {}
        rows.append(dict(run=rid, role=s['role'], label=s.get('label'), alias=s.get('alias'), baseline=s.get('baseline'), bases=s.get('bases'), observed_mean_depth=s.get('observed_mean_depth'),
                         expected_mean_depth=s.get('expected_mean_depth'), depth_ratio=s.get('depth_ratio'), median_depth=s.get('median_depth'), zero_fraction=s.get('zero_fraction'),
                         percent_duplication=dm.get('PERCENT_DUPLICATION'), flags=' | '.join(s.get('flags') or [])))
    with proj.run('bchrom', 'prepare', params=dict(cfg=cfg, failed=sorted(failures))) as run:
        pd.DataFrame(rows).to_csv(run.dir + 'runs.tsv', sep='\t', index=False)
        os.makedirs(run.dir + 'summaries/', exist_ok=True); os.makedirs(run.dir + 'ena/', exist_ok=True)
        for rid in results:
            shutil.copy2(f'{work}depth/{rid}.summary.json', run.dir + f'summaries/{rid}.json')
            if os.path.exists(f'{work}ena/{rid}.tsv'):
                shutil.copy2(f'{work}ena/{rid}.tsv', run.dir + f'ena/{rid}.tsv')
        if (cfg.get('slices') or {}).get('enabled') and cfg['slices'].get('upload') and os.path.isdir(work + 'slices/'):
            os.makedirs(run.dir + 'slices/', exist_ok=True)
            for f in sorted(os.listdir(work + 'slices/')):
                if any(f.startswith(rid + '.') for rid in results):
                    shutil.copy2(work + 'slices/' + f, run.dir + 'slices/' + f)
            json.dump(g_regions(cfg, proj), open(run.dir + 'slices/regions.json', 'w'), indent=1)
        warn = [r for r in rows if r['flags']]
        if warn:
            print('[bchrom] runs with warnings:'); [print('  ', r['run'], r['flags']) for r in warn]


# =================================================================================================== call
def _excluded_mask(proj, arrays, bin_size):
    try:
        from cichlid_common.genomes_registry import excluded_regions
        ex = excluded_regions(proj)
    except Exception:
        ex = {}
    mask = {}
    for c, regs in ex.items():
        if c in arrays:
            m = np.zeros(len(arrays[c]), dtype=bool)
            for s, e in regs:
                m[s // bin_size: -(-e // bin_size)] = True
            mask[c] = m
    return mask, ex


def load_tracks(work, cfg):
    tracks = {}
    for s in cfg['runs']:
        p = f"{work}depth/{s['id']}"
        if os.path.exists(p + '.npz') and os.path.exists(p + '.json'):
            a, meta = D.load_track(p)
            if meta['bin_size'] != cfg['bin_size']:
                raise SystemExit(f"{s['id']}: track has bin_size {meta['bin_size']} but config says {cfg['bin_size']}; run `bchrom prepare` again (or restore the setting)")
            tracks[s['id']] = (a, meta)
    return tracks


def call_step(proj, clean_depth=False, cfg=None):
    cfg = cfg or cfg_of(); by = validate_config(cfg)
    work = work_dir(proj, cfg); b = cfg['bin_size']; thr = cfg['thresholds']
    tr = load_tracks(work, cfg)
    carriers = [i for i, s in by.items() if s['role'] == 'carrier' and i in tr and by[i]['baseline'] in tr]
    controls = [i for i, s in by.items() if s['role'] == 'control' and i in tr and by[i]['baseline'] in tr]
    need = cfg['core']['min_individuals']
    if len(carriers) < need:
        raise SystemExit(f'only {len(carriers)} carriers have a depth track and a baseline (need {need}); run `bchrom prepare` / `bchrom status`')
    arrays = {i: t[0] for i, t in tr.items()}
    sf = {i: D.scale_factor(t[1]['stats'], cfg['scale_by']) for i, t in tr.items()}
    mean = {i: t[1]['stats']['mean'] for i, t in tr.items()}
    lengths = {i: t[1]['lengths'] for i, t in tr.items()}
    mask, excl = _excluded_mask(proj, arrays[carriers[0]], b)
    print(f'[bchrom] {len(carriers)} carriers, {len(controls)} controls, {len(tr)} tracks; scaling by the {cfg["scale_by"]}; excluded regions: {excl or "none"}')
    blocks = {}
    for i in carriers + controls:
        bl = by[i]['baseline']; p_exp = mean[i] / (mean[i] + mean[bl])
        blocks[i] = C.call_blocks(arrays[i], arrays[bl], sf[i], sf[bl], p_exp, b, thr, mask=mask, lengths=lengths[i])
        print(f"  {by[i]['label']}: {len(blocks[i])} blocks, {blocks[i].length.sum() / 1e6:.2f} Mb")
    core = C.core_regions({i: blocks[i] for i in carriers}, b, need, cfg['core']['min_len_bp'])
    dropped = 0
    if cfg['core'].get('exclude_blocks_in_controls') and len(core):
        allc = pd.concat([blocks[i] for i in controls]) if controls else pd.DataFrame(columns=core.columns)
        keep = [not ((allc.chrom == r.chrom) & (allc.start < r.end) & (allc.end > r.start)).any() for r in core.itertuples()]
        dropped = int(len(core) - sum(keep)); core = core[keep].reset_index(drop=True)
    core['length'] = core.end - core.start
    core.insert(0, 'core_id', [f'core{k + 1:04d}' for k in range(len(core))])
    core.insert(1, 'LG', core.chrom.map(getattr(proj, 'lg', {}) or {}))
    print(f'[bchrom] core set: {len(core)} blocks, {core.length.sum() / 1e6:.2f} Mb (Clark: 619 blocks in the files, 1.37 Mb){"; " + str(dropped) + " dropped as present in a control" if dropped else ""}')

    sig = C.block_signal(core, arrays, sf, b)
    base_ids = [i for i in by if by[i]['role'] == 'baseline' and i in tr]
    for col, ids in (('median_carriers', carriers), ('median_controls', controls), ('median_baselines', base_ids)):
        core[col] = sig[ids].median(axis=1).values if ids else np.nan
    core['tier'] = np.where(core.median_baselines >= cfg['core']['tier_a_min_baseline'], 'A', 'B')      # B: the male pools have little coverage there: a high ratio may be a low denominator
    summ = C.signal_summary(sig).rename(columns={'sample': 'run', 'n_blocks': 'n_core_blocks'})      # n_blocks (below) = blocks called in that sample
    sA = C.signal_summary(sig[(core.tier == 'A').values]).rename(columns={'sample': 'run', 'n_blocks': 'tierA_n_blocks', 'median_scaled': 'tierA_median_scaled', 'frac_ge2': 'tierA_frac_ge2', 'frac_ge3': 'tierA_frac_ge3'})
    summ = summ.merge(sA[['run', 'tierA_n_blocks', 'tierA_median_scaled', 'tierA_frac_ge2', 'tierA_frac_ge3']], on='run')
    summ.insert(1, 'role', summ.run.map(lambda i: by[i]['role'])); summ.insert(2, 'label', summ.run.map(lambda i: by[i].get('label')))
    # in-sample blocks per carrier / control; leave-one-out evaluation of the carriers (the core set was chosen because the carriers are high there: that inflates their own signal)
    extra = []
    loo_min = max(1, int(round(need * (len(carriers) - 1) / len(carriers))))
    for i in summ.run:
        r = dict(run=i)
        if i in blocks:
            a, bsp = C.b_space_mb(blocks[i]); r.update(n_blocks=len(blocks[i]), a_space_mb=a, b_space_mb=bsp, published_a_space_mb=by[i].get('published_A_mb'), published_b_space_mb=by[i].get('published_B_mb'))
        if i in carriers:
            core_i = C.core_regions({j: blocks[j] for j in carriers if j != i}, b, loo_min, cfg['core']['min_len_bp'])
            s_i = C.signal_summary(C.block_signal(core_i, {i: arrays[i]}, {i: sf[i]}, b)).iloc[0] if len(core_i) else None
            r.update(loo_core_blocks=len(core_i), loo_core_mb=float((core_i.end - core_i.start).sum() / 1e6) if len(core_i) else 0.0)
            if s_i is not None:
                r.update(loo_median_scaled=s_i.median_scaled, loo_mean_scaled=s_i.mean_scaled, loo_frac_ge2=s_i.frac_ge2, loo_frac_ge3=s_i.frac_ge3)
        extra.append(r)
    summ = summ.merge(pd.DataFrame(extra), on='run', how='left')
    sig_out = sig.copy(); sig_out.insert(0, 'core_id', core.core_id.values)
    text = _summary_text(cfg, by, core, summ, carriers, controls, dropped)
    print(text)
    with proj.run('bchrom', 'call', params=dict(cfg=cfg, carriers=carriers, controls=controls, excluded=excl, loo_min=loo_min)) as run:
        core[['chrom', 'start', 'end', 'core_id', 'n_carriers']].to_csv(run.dir + f'core_blocks_{proj.genome}.bed', sep='\t', header=False, index=False)
        core.to_csv(run.dir + 'core_blocks.tsv', sep='\t', index=False)
        sig_out.rename(columns={i: f"{by[i]['label']}|{i}" for i in sig.columns}).to_csv(run.dir + 'block_signal.tsv', sep='\t', index=False)
        summ.to_csv(run.dir + 'sample_summary.tsv', sep='\t', index=False)
        for i in carriers + controls:
            blocks[i].to_csv(run.dir + f"blocks_{by[i]['label']}_{i}.tsv", sep='\t', index=False)
        open(run.dir + 'summary.txt', 'w').write(text)
    if clean_depth:
        clean_step(proj, force=True, cfg=cfg)
    return core, summ


def _fmt_range(s):
    s = s.dropna()
    return f'{s.min():.2f}-{s.max():.2f}' if len(s) else 'n/a'


def _summary_text(cfg, by, core, summ, carriers, controls, dropped):
    ca, co = summ[summ.run.isin(carriers)], summ[summ.run.isin(controls)]
    ba = summ[summ.role == 'baseline']
    L = [f'B chromosome core set on the project reference: {len(core)} blocks, {core.length.sum() / 1e6:.2f} Mb (A space; Clark 2018: 1.37 Mb, 619 blocks in the published files)',
         f'  blocks per carrier: {int(ca.n_blocks.min())}-{int(ca.n_blocks.max())}; A space {_fmt_range(ca.a_space_mb)} Mb (Clark: 3.51-14.06); B space {_fmt_range(ca.b_space_mb)} Mb (Clark: 23.19-99.69)',
         f'  controls: {", ".join(str(int(x)) for x in co.n_blocks)} blocks, A space {_fmt_range(co.a_space_mb)} Mb (Clark 0.28-1.52): the false-positive load of a single B-free sample',
         '', 'signal over the core blocks (scaled coverage = depth / genome-wide depth; a B-free sample is ~1, a carrier is well above 2):',
         f'  carriers, in sample        median {_fmt_range(ca.median_scaled)}   fraction of blocks >= 2x: {_fmt_range(ca.frac_ge2)}   >= 3x: {_fmt_range(ca.frac_ge3)}',
         f'  carriers, leave-one-out    median {_fmt_range(ca.loo_median_scaled)}   fraction of blocks >= 2x: {_fmt_range(ca.loo_frac_ge2)}   >= 3x: {_fmt_range(ca.loo_frac_ge3)}   '
         '(core defined WITHOUT that carrier: what a new carrier would show)',
         f'  B-free controls            median {_fmt_range(co.median_scaled)}   fraction >= 2x: {_fmt_range(co.frac_ge2)}   >= 3x: {_fmt_range(co.frac_ge3)}',
         f'  baselines (male pools)     median {_fmt_range(ba.median_scaled)}   fraction >= 2x: {_fmt_range(ba.frac_ge2)}', '']
    nA, nB = int((core.tier == 'A').sum()), int((core.tier == 'B').sum()); mbA, mbB = core[core.tier == 'A'].length.sum() / 1e6, core[core.tier == 'B'].length.sum() / 1e6
    L += [f'tiers of the core set: A (male-pool coverage >= {cfg["core"]["tier_a_min_baseline"]}, the ratio is a high numerator): {nA} blocks, {mbA:.2f} Mb; B (poorly covered in the male pools, the ratio may be a low denominator): {nB} blocks, {mbB:.2f} Mb',
          f'  tier A only: carriers median {_fmt_range(ca.tierA_median_scaled)} (fraction >= 2x {_fmt_range(ca.tierA_frac_ge2)}); controls median {_fmt_range(co.tierA_median_scaled)} (fraction >= 2x {_fmt_range(co.tierA_frac_ge2)})', '']
    checks = []
    if len(ca) and len(co):
        gap = ca.loo_median_scaled.min() - co.median_scaled.max()
        checks.append(('every carrier (leave-one-out) has a median scaled coverage of at least 2', bool(ca.loo_median_scaled.min() >= 2)))
        checks.append(('every control has a median scaled coverage below 1.5', bool(co.median_scaled.max() < 1.5)))
        checks.append((f'carriers and controls are separated by the median (gap {gap:.2f}) and by the fraction of blocks >= 2x', bool(gap > 1 and ca.loo_frac_ge2.min() > co.frac_ge2.max() + 0.3)))
    for t, ok in checks:
        L.append(('OK   ' if ok else 'CHECK') + ' ' + t)
    return '\n'.join(L) + '\n'


# =================================================================================================== status / clean
def status_step(proj, cfg=None):
    cfg = cfg or cfg_of(); validate_config(cfg)
    work = work_dir(proj, cfg)
    free = shutil.disk_usage(work).free / 1e9
    print(f'work folder {work}   free {free:.0f} GB   in use {W.dir_size(work) / 1e9:.1f} GB')
    rows = []
    for s in cfg['runs']:
        rp = W.RunPaths(work, s['id'])
        done = os.path.exists(rp.summary) and os.path.exists(rp.track + '.npz')
        stages = [st for st in W.STAGES if os.path.exists(rp.stages + st + '.done')]
        rows.append(dict(run=s['id'], role=s['role'], label=s.get('label', '')[:34], state='complete' if done else (stages[-1] + ' done' if stages else ('started' if os.path.isdir(rp.dir) else 'not started')),
                         scratch_gb=round(W.dir_size(rp.dir) / 1e9, 1) if os.path.isdir(rp.dir) else 0))
    t = pd.DataFrame(rows)
    print(t.to_string(index=False)); print(f"\n{(t.state == 'complete').sum()} of {len(t)} runs complete")
    return t


def clean_step(proj, force=False, cfg=None):
    """Delete the work folder (depth tracks, scratch, cached ENA listings). Refuses unless `bchrom call` has produced a run (its outputs are on Dropbox), or --force."""
    cfg = cfg or cfg_of(); work = work_dir(proj, cfg)
    if not force and not os.path.exists(f'{proj.root}bchrom/call/LATEST'):
        raise SystemExit('no `bchrom call` run found: the depth tracks are the only copy of the data; run `bchrom call` first (or use --force)')
    if not os.path.exists(work + MARKER):
        raise SystemExit(f'{work} was not created by bchrom (no {MARKER} file): refusing to delete it')
    n = W.dir_size(work)
    shutil.rmtree(work, ignore_errors=True)
    print(f'[bchrom] deleted {work} ({n / 1e9:.1f} GB)')


# =================================================================================================== compare
def compare_step(proj, clark_bed, paf=None, name_map=None, slop=0, cfg=None):
    from cichlid_common.runner import latest_dir
    cfg = cfg or cfg_of()
    if not clark_bed:
        raise SystemExit('usage: python -m yhsd.cli bchrom compare --clark-bed <one of Appendix A Core_blocks/*.bed> [--paf <Clark assembly vs reference>] [--name-map <tsv>]')
    cdir = latest_dir(proj, 'bchrom', 'call')
    ours = pd.read_csv(cdir + 'core_blocks.tsv', sep='\t')
    if paf is None:
        from cichlid_common.genome_align import ensure_paf
        mm = yaml.safe_load(open(config_path('inversion_breakpoints.yaml')))['minimap2']
        paf = ensure_paf(proj, cfg['clark']['assembly_key'], proj.genome, mm)
    from cichlid_common import paf as P
    df = P.parse_paf(paf)
    blocks = K.read_clark_bed(clark_bed)
    nm = K.read_name_map(name_map) if name_map else None
    if not blocks.chrom.isin(set(df.qname.unique()) | set(nm or {})).any():
        print('none of Clark\'s sequence names is in the alignment. Names in his file:', blocks.chrom.unique()[:5].tolist(), '| names in the alignment:', df.qname.unique()[:5].tolist(),
              '\nGive --name-map (two tab-separated columns: his name, the name in the aligned FASTA).')
        raise SystemExit(1)
    lifted, summ = K.compare(K.lift_blocks(blocks, df, nm), ours, slop)
    with proj.run('bchrom', 'compare', params=dict(clark_bed=os.path.basename(clark_bed), paf=os.path.basename(paf), slop=slop)) as run:
        lifted.to_csv(run.dir + 'clark_blocks_lifted.tsv', sep='\t', index=False)
        json.dump(summ, open(run.dir + 'comparison.json', 'w'), indent=1)
        txt = []
        for k in ('all', 'anchored', 'unplaced'):
            s = summ[k]
            txt.append(f"{k:9s} {s['n_clark_blocks']:4d} Clark blocks ({s['bp'] / 1e6:.2f} Mb): {s['n_lifted']} lift ({s['frac_lifted']:.0%}); of those {s['n_recovered']} are inside our core set ({s['frac_recovered_of_lifted']:.0%})")
        open(run.dir + 'comparison.txt', 'w').write('\n'.join(txt) + '\n'); print('\n'.join(txt))
