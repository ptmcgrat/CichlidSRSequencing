"""inversion_breakpoints analysis. Steps (python -m yhsd.cli inversion_breakpoints <step>):
  align    minimap2 of every inversion-carrying assembly to the reference; keeps alignments on the target chromosomes and records which reference chromosome each contig belongs to
  call     per chromosome: breakpoints per assembly (methods A, B, C), consensus, strand polarity, dotplots, BED
  figure   overview table + draft diagram of the inversions (and centromeres if `centromeres scan` has been run), all in reference coordinates
  all      align, call, figure"""
import json, os
import numpy as np
import pandas as pd
import yaml
from cichlid_common.common_light import sh, config_path
from cichlid_common.genomes_registry import fasta_path, entry_dir, excluded_regions
from cichlid_common.runner import latest_dir
from cichlid_common import paf as P
from cichlid_common.genome_align import ensure_paf
from ..centromeres import cen as cenmod
from . import parts as parts_mod


def cfg_of():
    return yaml.safe_load(open(config_path('inversion_breakpoints.yaml')))


def assemblies_of(cfg):
    return list(dict.fromkeys(a for spec in cfg['inversions'].values() for a in spec['assemblies']))


def targets(proj, cfg):
    return {lg: proj.chrom_for_lg(lg) for lg in cfg['inversions']}


def _genome_md5(proj, key):
    return json.load(open(entry_dir(proj, key) + 'manifest.json')).get('md5_fasta')


def align_step(proj):
    cfg = cfg_of(); ref = cfg['reference']
    if proj.genome != ref:
        raise SystemExit(f"config genome_version is {proj.genome} but inversion_breakpoints.reference is {ref}")
    tch = targets(proj, cfg); keep = set(tch.values())
    ref_fa = fasta_path(proj, ref)
    queries = assemblies_of(cfg)
    with proj.run('inversion_breakpoints', 'align', params=dict(cfg=cfg, targets=tch, md5={k: _genome_md5(proj, k) for k in queries + [ref] if os.path.exists(entry_dir(proj, k) + 'manifest.json')})) as run:
        for q in queries:
            try:
                qfa = fasta_path(proj, q)
            except FileNotFoundError as e:
                print(f'[align] {q}: skipped ({e})'); continue
            full = ensure_paf(proj, q, ref, cfg['minimap2'])
            n = 0; share = {}
            with open(full) as fi, open(f'{run.dir}{q}_vs_{ref}.targets.paf', 'w') as fo:
                for line in fi:
                    f = line.split('\t', 11)
                    share[(f[0], f[5])] = share.get((f[0], f[5]), 0) + int(f[10])          # aligned bases per (contig, reference chromosome)
                    if f[5] in keep:
                        fo.write(line); n += 1
            cm = pd.DataFrame([(c, t, b) for (c, t), b in share.items()], columns=['contig', 'reference_chrom', 'aligned_bp'])
            cm['contig_total'] = cm.groupby('contig').aligned_bp.transform('sum'); cm['share'] = cm.aligned_bp / cm.contig_total
            cm = cm.sort_values(['contig', 'aligned_bp'], ascending=[True, False]).groupby('contig').head(1)
            cm['lg'] = cm.reference_chrom.map(proj.lg)
            cm.to_csv(f'{run.dir}{q}_contig_map.tsv', sep='\t', index=False)
            print(f'[align] {q}: kept {n} blocks on {sorted(keep)}; {len(cm)} contigs mapped')


def dotplot(df, tname, contigs, consensus_df, title, path_base, excluded=()):
    """One panel per contig (shared reference axis): contigs no longer overlap on a common y-axis. Dashed lines = consensus breakpoints with high/medium confidence."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.fonttype': 'none', 'font.size': 6.5})
    contigs = list(contigs)[:4]
    n = max(1, len(contigs))
    fig, axs = plt.subplots(n, 1, figsize=(4.4, 1.9 * n + 0.7), sharex=True, squeeze=False)
    xs = []
    for q in contigs:
        g = df[(df.qname == q) & (df.tname == tname)]
        xs += [g.tstart.min(), g.tend.max()]
    bps = []
    for c in (consensus_df.itertuples() if consensus_df is not None and len(consensus_df) else []):
        if getattr(c, 'confidence', 'high') in ('high', 'medium'):
            bps += [c.left_breakpoint, c.right_breakpoint]
    for ax, q in zip(axs.flat, contigs):
        g = df[(df.qname == q) & (df.tname == tname)]
        for r in g.itertuples():
            if r.strand == '+': ax.plot([r.tstart / 1e6, r.tend / 1e6], [r.qstart / 1e6, r.qend / 1e6], color='#2a6fb0', lw=1.3)
            else: ax.plot([r.tstart / 1e6, r.tend / 1e6], [r.qend / 1e6, r.qstart / 1e6], color='#d9531e', lw=1.3)
        for x in bps:
            ax.axvline(x / 1e6, color='#555555', lw=0.6, ls='--')
        for s, e in excluded:                                              # known misassembly of the reference: shown, but ignored by the caller
            ax.axvspan(s / 1e6, e / 1e6, color='#999999', alpha=0.25, lw=0)
        ax.set_ylabel(f'{q}\nposition (Mb)', fontsize=5.5)
        ax.tick_params(labelsize=5.5)
    axs.flat[-1].set_xlabel(f'reference {tname} (Mb)')
    axs.flat[0].set_title(title + '   (blue = same strand, orange = reverse)', fontsize=6.5)
    for ext in ('svg', 'png'):
        fig.savefig(f'{path_base}.{ext}', dpi=250, bbox_inches='tight')
    plt.close(fig)


def haplotig_table(df, tname):
    rows = []
    for q, g in df[df.tname == tname].groupby('qname'):
        plus = int(g[g.strand == '+'].alen.sum()); minus = int(g[g.strand == '-'].alen.sum())
        rows.append(dict(contig=q, qlen=int(g.qlen.iloc[0]), aligned_plus=plus, aligned_minus=minus, ref_start=int(g.tstart.min()), ref_end=int(g.tend.max()),
                         dominant_strand='+' if plus >= minus else '-', main_target=P.main_target(df, q) == tname))
    return pd.DataFrame(rows)


def call_step(proj, align_dir=None):
    cfg = cfg_of(); ref = cfg['reference']; cc = cfg['calling']
    adir = align_dir or latest_dir(proj, 'inversion_breakpoints', 'align')
    L = proj.chrom_lengths()
    excl = excluded_regions(proj)
    if excl: print('[call] ignoring alignments in known misassembled regions of the reference:', {k: v for k, v in excl.items()})
    with proj.run('inversion_breakpoints', 'call', params=dict(cfg=cfg, align_run=adir, excluded_regions=excl)) as run:
        summary = {}
        for lg, spec in cfg['inversions'].items():
            tname = proj.chrom_for_lg(lg)
            allcalls, hap_rows, pafs = [], [], {}
            for q in spec['assemblies']:
                f = f'{adir}{q}_vs_{ref}.targets.paf'
                if not os.path.exists(f):
                    print(f'[call] {lg}: no alignment file for {q}; run `align` first'); continue
                run.add_input(f)
                df_all = P.parse_paf(f, {tname}); pafs[q] = df_all                          # the dotplot shows everything
                df = P.drop_blocks_in_regions(df_all, excl.get(tname, []))                  # the caller does not see the misassembled region
                if not len(df):
                    print(f'[call] {lg}: no alignments of {q} on {tname}'); continue
                ac = P.call_method_AC(df, tname, cc['min_block'], cc['min_inversion'], None, cc.get('despeckle_run', 300000), cc.get('max_breakpoint_interval', 1000000)).assign(assembly=q)
                b = P.discordance_regions(df, tname, L[tname], cc['window'], cc['min_window_bp'], cc['min_block'], cc['min_inversion']).assign(assembly=q)
                allcalls += [x for x in (ac, b) if len(x)]
                h = haplotig_table(df, tname); h['assembly'] = q; hap_rows.append(h)
            calls = pd.concat(allcalls, ignore_index=True) if allcalls else pd.DataFrame()
            calls.to_csv(f'{run.dir}{lg}_breakpoints_per_assembly.tsv', sep='\t', index=False)
            if hap_rows: pd.concat(hap_rows).to_csv(f'{run.dir}{lg}_haplotigs.tsv', sep='\t', index=False)
            run_rows = []
            for q, df in pafs.items():
                if len(df):
                    df = P.drop_blocks_in_regions(df, excl.get(tname, []))
                    top = df.groupby('qname').alen.sum().sort_values(ascending=False).index[:6].tolist()
                    rt = P.runs_table(df, tname, top, cc['min_block'], cc.get('despeckle_run', 300000)); rt.insert(0, 'assembly', q); run_rows.append(rt)
            if run_rows: pd.concat(run_rows).to_csv(f'{run.dir}{lg}_strand_runs.tsv', sep='\t', index=False)      # what the caller sees: strand runs per contig after despeckling
            cons = P.consensus(calls, cc['consensus_tolerance'], cc.get('consensus_gap', 100_000)) if len(calls) else pd.DataFrame()
            cons.to_csv(f'{run.dir}{lg}_consensus_breakpoints.tsv', sep='\t', index=False)
            parts_info, parts_text = {}, ''
            if spec.get('parts') and len(calls) and len(cons):                       # a tandem (double) inversion: where each part lies and who carries which (parts.py)
                pt, st, pnotes = parts_mod.resolve_parts(calls, cons, spec['parts'], spec.get('expect'), cc['consensus_tolerance'], cc.get('consensus_gap', 100_000), spec['assemblies'])
                pt.to_csv(f'{run.dir}{lg}_parts.tsv', sep='\t', index=False); st.to_csv(f'{run.dir}{lg}_assembly_states.tsv', sep='\t', index=False)
                with open(f'{run.dir}{lg}_parts_breakpoint_regions.bed', 'w') as fb:
                    fb.write('\n'.join(parts_mod.parts_bed(pt, tname, lg)) + ('\n' if len(pt) else ''))
                parts_info = dict(parts=pt.to_dict('records'), assembly_states=dict(zip(st.assembly, st.state)), parts_notes=pnotes,
                                  differs_from_literature=[r.assembly for r in st.itertuples() if r.matches_expectation is False])
                parts_text = parts_mod.format_parts(pt, st, pnotes)
            info = dict(lg=lg, reference_chrom=tname, assemblies=spec['assemblies'], n_consensus=int(len(cons)),
                        reference_is_standard=spec.get('reference_is_standard'), polarity=(P.strand_polarity(spec.get('reference_is_standard')) if len(cons) else 'no inversion called: polarity not applicable'))
            info.update(parts_info)
            json.dump(info, open(f'{run.dir}{lg}_polarity.json', 'w'), indent=1, default=str)
            with open(f'{run.dir}{lg}_breakpoint_regions.bed', 'w') as b1, open(f'{run.dir}{lg}_breakpoint_flanks_50kb.bed', 'w') as b2:
                for i, c in enumerate(cons.itertuples() if len(cons) else []):
                    for side, lo, hi in (('left', c.left_min, c.left_max), ('right', c.right_min, c.right_max)):
                        b1.write(f'{tname}\t{lo}\t{hi}\t{lg}_inv{i + 1}_{side}\n'); b2.write(f'{tname}\t{max(0, lo - 50000)}\t{hi + 50000}\t{lg}_inv{i + 1}_{side}_flank\n')
            for q, df in pafs.items():
                if len(df):
                    top = df.groupby('qname').alen.sum().sort_values(ascending=False).index[:6].tolist()
                    dotplot(df, tname, top, cons, f'{lg}: {q} vs {ref}', f'{run.dir}{lg}_{q}_dotplot', excl.get(tname, []))
            summary[lg] = info
            print(f'\n== {lg} ({tname}) ==')
            print(cons.to_string(index=False) if len(cons) else 'no inversion found on this chromosome')
            print('polarity:', info['polarity'])
            if parts_text:
                print(f'-- the parts of {lg} --\n' + parts_text)
        json.dump(summary, open(run.dir + 'summary.json', 'w'), indent=1, default=str)


# ------------------------------------------------------------------ overview for the figure
def overview_table(proj_lengths, lg_to_chrom, cons_by_lg, centromeres, order, excluded=None):
    """One row per chromosome: reference length, the top consensus inversion and the centromere (if known). centromeres: DataFrame(lg, ref_start, ref_end, confidence)."""
    rows = []
    for lg in order:
        chrom = lg_to_chrom[lg]; r = dict(lg=lg, chrom=chrom, length=int(proj_lengths[chrom]))
        c = cons_by_lg.get(lg)
        if c is not None and len(c):
            top = c.iloc[0]; r.update(inv_n_segments=int(top.get('n_segments', 1)), inv_segments=top.get('segments', '') if isinstance(top.get('segments', ''), str) else '', inv_start=int(top.left_breakpoint), inv_end=int(top.right_breakpoint), inv_left_min=int(top.left_min), inv_left_max=int(top.left_max),
                                      inv_right_min=int(top.right_min), inv_right_max=int(top.right_max), inv_n_calls=int(top.n_calls), inv_sources=top.sources)
        if centromeres is not None and len(centromeres):
            m = centromeres[centromeres.lg == lg]
            if len(m):
                m = m.iloc[0]; r.update(cen_start=int(m.ref_start), cen_end=int(m.ref_end), cen_confidence=m.confidence, cen_size=m.get('size', np.nan))
                if 'inv_start' in r:
                    r['cen_relation'] = cenmod.relation_to_inversion(r['cen_start'], r['cen_end'], r['inv_start'], r['inv_end'])
        for s, e in (excluded or {}).get(chrom, [])[:1]:
            r.update(excl_start=int(s), excl_end=int(e))
        rows.append(r)
    return pd.DataFrame(rows)


def draw_overview(ov, path_base, label_of=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Polygon
    plt.rcParams.update({'svg.fonttype': 'none', 'font.size': 7})
    fig, ax = plt.subplots(figsize=(6.5, 0.55 * len(ov) + 0.6))
    ax.set_xlim(-6, ov.length.max() / 1e6 + 2); ax.set_ylim(len(ov), -0.8); ax.axis('off')
    for i, r in enumerate(ov.itertuples()):
        ax.add_patch(Rectangle((0, i - 0.12), r.length / 1e6, 0.24, fc='#e6e6e6', ec='#888888', lw=0.5))
        ax.text(-0.6, i, (label_of or {}).get(r.lg, r.lg), ha='right', va='center', fontsize=7, weight='bold')
        if hasattr(r, 'excl_start') and not pd.isna(r.excl_start):
            ax.add_patch(Rectangle((r.excl_start / 1e6, i - 0.12), (r.excl_end - r.excl_start) / 1e6, 0.24, fc='none', ec='#888888', hatch='////', lw=0.4))
        if hasattr(r, 'inv_start') and not pd.isna(r.inv_start):
            ax.add_patch(Rectangle((r.inv_start / 1e6, i - 0.2), (r.inv_end - r.inv_start) / 1e6, 0.4, fc='#e76f51', ec='none', alpha=0.9))
            ax.text((r.inv_start + r.inv_end) / 2e6, i - 0.28, f'{(r.inv_end - r.inv_start) / 1e6:.1f} Mb', ha='center', va='bottom', fontsize=5.5)
        if hasattr(r, 'cen_start') and not pd.isna(r.cen_start):
            x = (r.cen_start + r.cen_end) / 2e6
            if (r.cen_end - r.cen_start) > 150_000:                              # the interval itself (the triangle marks its middle)
                ax.add_patch(Rectangle((r.cen_start / 1e6, i + 0.12), (r.cen_end - r.cen_start) / 1e6, 0.07, fc='#333333', ec='none', alpha=0.55))
            ax.add_patch(Polygon([[x - 0.6, i + 0.38], [x + 0.6, i + 0.38], [x, i + 0.14]], closed=True, fc='#333333', ec='none'))
    ax.text(0, len(ov) - 0.1, 'orange: inversion (consensus breakpoints, reference coordinates)    black triangle: centromere    hatched: reference misassembly (ignored)', fontsize=5.5, color='#555555', va='top')
    for ext in ('svg', 'png'):
        fig.savefig(f'{path_base}.{ext}', dpi=250, bbox_inches='tight')
    plt.close(fig)


def figure_step(proj, call_dir=None):
    cfg = cfg_of()
    cdir = call_dir or latest_dir(proj, 'inversion_breakpoints', 'call')
    L = proj.chrom_lengths(); lg2c = {lg: proj.chrom_for_lg(lg) for lg in cfg['inversions']}
    cons = {}
    for lg in cfg['inversions']:
        f = f'{cdir}{lg}_consensus_breakpoints.tsv'
        if os.path.exists(f) and os.path.getsize(f) > 5:
            c = pd.read_csv(f, sep='\t')
            if 'confidence' in c.columns:
                c = c[c.confidence.isin(['high', 'medium'])].reset_index(drop=True)          # calls that only assume the reference orientation are never drawn
            cons[lg] = c
    cen = None
    try:
        cen_dir = latest_dir(proj, 'centromeres', 'scan'); c = pd.read_csv(cen_dir + 'centromeres.tsv', sep='\t')
        c = c[c.confidence.isin(['high', 'medium']) & c.ref_start.notna() & c.lg.notna()].copy()
        c['is_ref'] = (c.genome == cfg['reference']).astype(int)
        cen = c.sort_values('is_ref', ascending=False).drop_duplicates('lg')[['lg', 'ref_start', 'ref_end', 'confidence', 'genome']]    # prefer the reference's own call
    except FileNotFoundError:
        print('[figure] no centromere scan found')
    try:                                                                                  # published centromeres lifted into the reference win over the satellite scan
        pub = pd.read_csv(latest_dir(proj, 'centromeres', 'lift') + 'published_centromeres_lifted.tsv', sep='\t')
        for c_ in ('bounded_start', 'bounded_end'):
            if c_ not in pub: pub[c_] = np.nan                                                  # lift tables from before v23 have no bounded columns
        pub['ref_start'] = pub.gt3_start.fillna(pub.bounded_start); pub['ref_end'] = pub.gt3_end.fillna(pub.bounded_end)           # clean lift, else the region between the flank landings
        pub = pub[pub.ref_start.notna()][['lg', 'ref_start', 'ref_end', 'lift_status']].assign(confidence=lambda d: np.where(d.lift_status == 'ok', 'published', 'published-bounded'), genome='published').drop(columns='lift_status')
        cen = pd.concat([pub, cen[~cen.lg.isin(pub.lg)]], ignore_index=True) if cen is not None else pub
        print(f'[figure] {len(pub)} centromeres from the published assembly (lifted), the rest from the satellite scan')
    except (FileNotFoundError, KeyError):
        pass
    try:                                                                                  # the array-level call (largest family array > 500 kb from both ends) is the primary definition
        arr = pd.read_csv(latest_dir(proj, 'centromeres', 'arrays') + 'centromeres_arrays.tsv', sep='\t')
        arr = arr[(arr.genome == cfg['reference']) & arr.start.notna() & arr.lg.notna()].copy()
        arr['ref_start'] = arr.start; arr['ref_end'] = arr.end; arr['genome'] = 'arrays'
        arr['confidence'] = np.where(arr.ambiguous == True, 'array-ambiguous', 'array')
        arr = arr[['lg', 'ref_start', 'ref_end', 'confidence', 'genome', 'size']]
        cen = pd.concat([arr, cen[~cen.lg.isin(arr.lg)]], ignore_index=True) if cen is not None else arr
        print(f'[figure] {len(arr)} centromeres from the satellite-array call; the rest from the published lift / scan')
    except (FileNotFoundError, KeyError):
        pass
    if cen is None:
        print('[figure] no centromere information; the overview shows inversions only')
    ov = overview_table(L, lg2c, cons, cen, list(cfg['inversions']), excluded_regions(proj))
    labels = {lg: spec.get('label', lg) for lg, spec in cfg['inversions'].items()}
    with proj.run('inversion_breakpoints', 'figure', params=dict(call_run=cdir)) as run:
        ov.to_csv(run.dir + 'inversion_overview.tsv', sep='\t', index=False)
        draw_overview(ov, run.dir + 'inversion_overview', labels)
        print(ov.to_string(index=False))
