"""centromeres analysis: python -m yhsd.cli centromeres scan
For every genome in config/centromeres.yaml that has a repeat annotation: satellite density per window, centromere call per chromosome (confidence-graded), assembly gaps,
profile plots, and (for non-reference genomes) the called interval lifted to the reference through the assembly alignments of inversion_breakpoints (bounded by the nearest alignable flanks)."""
import os
import numpy as np
import pandas as pd
import yaml
from ...common_light import config_path
from ...genomes_registry import fasta_path, annotation_for, excluded_regions
from ...runner import latest_dir
from ..inversion_breakpoints import paf as P
from . import cen


def cfg_of():
    return yaml.safe_load(open(config_path('centromeres.yaml')))


def scan_genome(name, fa, annotation, cfg, name_map=None, exclude=None):
    """Returns (calls DataFrame, dict of per-sequence profiles)."""
    seqs = cen.fasta_lengths_and_gaps(fa)
    ann = cen.parse_annotation(annotation, cfg['classes'])
    if name_map:
        ann['chrom'] = ann.chrom.map(name_map); ann = ann.dropna(subset=['chrom'])        # annotation names -> genome sequence names
    rows, prof = [], {}
    for seq, (length, gaps) in seqs.items():
        if length < cfg['min_chrom_length']:
            continue
        sub = ann[ann.chrom == seq]
        dens = cen.window_coverage(sub.start.values, sub.end.values, length, cfg['window'])
        gcov = cen.window_coverage([a for a, b in gaps], [b for a, b in gaps], length, cfg['window'])
        for s, e in (exclude or {}).get(seq, []):                              # known misassembly of the reference: no centromere evidence from there
            dens[int(s // cfg['window']):int(e // cfg['window']) + 1] = 0.0; gcov[int(s // cfg['window']):int(e // cfg['window']) + 1] = 0.0
        r = cen.call_centromere(dens, cfg['window'], cfg['smooth'], cfg['confidence'], gcov)
        r.update(genome=name, sequence=seq, length=length, n_satellite_features=int(len(sub)), masked_regions=';'.join(f'{s}-{e}' for s, e in (exclude or {}).get(seq, [])),
                 top_families=cen.top_families(sub, r['start'], r['end']) if r['confidence'] != 'none' else '')
        rows.append(r); prof[seq] = (dens, gcov, r)
    return pd.DataFrame(rows), prof


def plot_profiles(prof, window, smooth_bp, title, path_base):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.fonttype': 'none', 'font.size': 6})
    n = len(prof); ncol = 4; nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(7.2, 1.05 * nrow + 0.4), squeeze=False)
    colors = {'high': '#d9531e', 'medium': '#e9a23b', 'low': '#bbbbbb', 'none': '#dddddd'}
    for ax, (seq, (dens, g, r)) in zip(axs.flat, prof.items()):
        x = np.arange(len(dens)) * window / 1e6
        ax.plot(x, cen.smooth(dens, smooth_bp // window), color='#333333', lw=0.7)
        if r['confidence'] != 'none':
            ax.axvspan(r['start'] / 1e6, r['end'] / 1e6, color=colors[r['confidence']], alpha=0.5, lw=0)
        ax.set_title(f"{seq} ({r['confidence']})", fontsize=5.5, pad=1); ax.tick_params(labelsize=5); ax.set_ylim(0, max(0.05, dens.max()))
    for ax in list(axs.flat)[n:]:
        ax.axis('off')
    fig.suptitle(title + '  (line: smoothed satellite fraction)', fontsize=6.5)
    fig.tight_layout()
    for ext in ('svg', 'png'):
        fig.savefig(f'{path_base}.{ext}', dpi=220, bbox_inches='tight')
    plt.close(fig)


def scan_step(proj, align_dir=None):
    cfg = cfg_of(); ref = proj.genome
    try:
        adir = align_dir or latest_dir(proj, 'inversion_breakpoints', 'align')
    except FileNotFoundError:
        adir = None
    with proj.run('centromeres', 'scan', params=dict(cfg=cfg)) as run:
        all_calls = []
        for key, g in cfg['genomes'].items():
            spec = g.get('annotation')
            if not spec:
                print(f'[centromeres] {key}: no annotation configured; skipped'); continue
            try:
                fa = fasta_path(proj, key)
                if spec == 'registered':
                    apath, nmap = annotation_for(proj, key)
                else:
                    apath, nmap = os.path.expanduser(spec), None
            except FileNotFoundError as e:
                print(f'[centromeres] {key}: skipped ({e})'); continue
            calls, prof = scan_genome(key, fa, apath, cfg, nmap, excluded_regions(proj) if key == ref else None)
            if not len(calls):
                print(f'[centromeres] {key}: no chromosome-scale sequences with annotation'); continue
            # LG labels
            if key == ref:
                calls['lg'] = calls.sequence.map(proj.lg)
            elif adir and os.path.exists(f'{adir}{key}_contig_map.tsv'):
                cm = pd.read_csv(f'{adir}{key}_contig_map.tsv', sep='\t').set_index('contig')
                calls['lg'] = calls.sequence.map(cm.lg); calls['lg_share'] = calls.sequence.map(cm.share)
            else:
                calls['lg'] = np.nan
            # reference coordinates
            calls['ref_start'] = np.nan; calls['ref_end'] = np.nan; calls['lift_status'] = ''
            if key == ref:
                calls['ref_start'], calls['ref_end'], calls['lift_status'] = calls.start, calls.end, 'reference'
            elif adir and os.path.exists(f"{adir}{key}_vs_{ref}.targets.paf"):
                df = P.parse_paf(f"{adir}{key}_vs_{ref}.targets.paf")
                for i, r in calls.iterrows():
                    if r.confidence in ('high', 'medium') and not pd.isna(r.start):
                        lf = P.lift_interval(df, r.sequence, int(r.start), int(r.end))
                        if lf is None: calls.loc[i, 'lift_status'] = 'no alignable flanks (chromosome not in the alignment targets?)'
                        else:
                            calls.loc[i, 'lift_status'] = lf['status']
                            if lf['status'] == 'ok': calls.loc[i, ['ref_start', 'ref_end']] = [lf['tstart'], lf['tend']]
            all_calls.append(calls)
            plot_profiles(prof, cfg['window'], cfg['smooth'], key, f'{run.dir}{key}_centromere_profiles')
            pd.DataFrame({f'{s}': pd.Series(v[0]) for s, v in prof.items()}).to_csv(f'{run.dir}{key}_satellite_density_{cfg["window"]}bp.tsv.gz', sep='\t', index=False)
            print(f'\n== {key} ==\n' + calls[['sequence', 'lg', 'start', 'end', 'confidence', 'ratio', 'dominance', 'top_families', 'ref_start', 'ref_end']].round(2).to_string(index=False))
        if all_calls:
            pd.concat(all_calls, ignore_index=True)[['genome', 'sequence', 'lg', 'length', 'start', 'end', 'mid', 'confidence', 'peak', 'ratio', 'dominance', 'gap_fraction',
                                                    'n_satellite_features', 'masked_regions', 'top_families', 'ref_start', 'ref_end', 'lift_status']].to_csv(run.dir + 'centromeres.tsv', sep='\t', index=False)


def lift_step(proj, table):
    """centromeres lift --table <xlsx>: align the published near-complete assembly to the reference, lift every published centromere interval and centroid core into the
    reference's coordinates (bounded by the nearest alignable flanks), and compare with the satellite scan."""
    from . import published as pb
    from ..inversion_breakpoints.steps import ensure_paf
    if not table:
        raise SystemExit('usage: python -m yhsd.cli centromeres lift --table <10577_2026_9819_MOESM2_ESM.xlsx>')
    cfg = cfg_of(); pc = cfg['published']; ref = proj.genome
    mm = yaml.safe_load(open(config_path('inversion_breakpoints.yaml')))['minimap2']
    lg_to_chrom = {lg: c for c, lg in proj.lg.items()}
    with proj.run('centromeres', 'lift', params=dict(published=pc, table=os.path.basename(table))) as run:
        run.add_input(table)                                                          # the paper's table is archived in the run folder
        tbl = pb.read_published(table, pc['sheet'], pc['genome'])
        full = ensure_paf(proj, pc['assembly'], ref, mm)
        share = {}
        with open(full) as f:                                                         # pass 1: which sequence of the published assembly is which reference chromosome
            for line in f:
                fl = line.split('\t', 11)
                share[(fl[0], fl[5])] = share.get((fl[0], fl[5]), 0) + int(fl[10])
        cm = pd.DataFrame([(c, tg, b) for (c, tg), b in share.items()], columns=['contig', 'reference_chrom', 'aligned_bp'])
        cm['contig_total'] = cm.groupby('contig').aligned_bp.transform('sum'); cm['share'] = cm.aligned_bp / cm.contig_total
        cm = cm.sort_values(['contig', 'aligned_bp'], ascending=[True, False]).groupby('contig').head(1)
        mapping = pb.map_published_chr(cm, lg_to_chrom, list(tbl.lg.unique()))
        win = pc['window']; wins = {}
        for r in tbl.itertuples():
            m = mapping.get(r.lg)
            if m:
                wins.setdefault(m['seq'], []).append((r.start - win, r.end + win))
        keep = []
        with open(full) as f:                                                         # pass 2: only the alignments near the published centromeres
            for line in f:
                fl = line.split('\t', 11)
                if fl[0] in wins:
                    qs, qe = int(fl[2]), int(fl[3])
                    if any(qs < e and qe > s for s, e in wins[fl[0]]):
                        keep.append(line)
        df = P.parse_paf(keep)
        lifted = pb.lift_all(df, tbl, mapping, lg_to_chrom, win)
        scan = None
        try:
            s = pd.read_csv(latest_dir(proj, 'centromeres', 'scan') + 'centromeres.tsv', sep='\t')
            scan = s[s.genome == ref][['lg', 'start', 'end', 'confidence']]
        except FileNotFoundError:
            print('[lift] no satellite scan found: run `centromeres scan` to compare')
        out = pb.compare_with_scan(lifted, scan, pc['agree_distance'])
        out.to_csv(run.dir + 'published_centromeres_lifted.tsv', sep='\t', index=False)
        cm.to_csv(run.dir + 'published_assembly_contig_map.tsv', sep='\t', index=False)
        show = out[['lg', 'published_start', 'published_end', 'gt3_start', 'gt3_end', 'gt3_width', 'strand', 'lift_status', 'scan_start', 'scan_end', 'scan_confidence', 'scan_vs_published']]
        print(show.round(0).to_string(index=False))
        bad = out[out.lift_status != 'ok']
        print(f'\n{len(out) - len(bad)} of {len(out)} published centromeres lifted cleanly; check the others: {sorted(bad.lg)}')
        if scan is not None:
            print('scan agrees for', int((out.scan_vs_published == 'agrees').sum()), 'chromosomes; disagrees for', sorted(out[out.scan_vs_published.str.startswith('disagrees', na=False)].lg))


# ------------------------------------------------------------------ array-level centromeres (family-specific)
def arrays_for_genome(name, fa, annotation, cfg, name_map=None, exclude=None):
    """Arrays of the centromeric satellite family per chromosome-scale sequence. Returns (calls DataFrame, all-arrays DataFrame)."""
    ac = cfg['arrays']
    seqs = cen.fasta_lengths_and_gaps(fa)
    ann = cen.parse_annotation(annotation, cfg['classes'])
    if name_map:
        ann['chrom'] = ann.chrom.map(name_map); ann = ann.dropna(subset=['chrom'])
    ann = ann[ann.name == ac['family']]
    calls, allrows = [], []
    for seq, (length, _) in seqs.items():
        if length < cfg['min_chrom_length']:
            continue
        sub = ann[ann.chrom == seq]
        for s, e in (exclude or {}).get(seq, []):
            sub = sub[~((sub.start >= s) & (sub.end <= e))]                      # features inside a known misassembly are ignored
        res = cen.call_by_arrays(cen.merge_arrays(sub.start.values, sub.end.values, ac['max_gap']), length, ac['min_end_distance'], ac['min_array'])
        b, o = res['best'], res['best_overall']
        row = dict(genome=name, sequence=seq, length=length, family=ac['family'], n_family_features=int(len(sub)), n_arrays=res['n_arrays'], n_internal_arrays=res['n_internal'],
                   start=b['start'] if b else np.nan, end=b['end'] if b else np.nan, size=b['size'] if b else np.nan, span=b['span'] if b else np.nan, n_features=b['n_features'] if b else np.nan,
                   end_distance=b['end_distance'] if b else np.nan, second_internal_size=res['second_internal_size'], ambiguous=res['ambiguous'],
                   largest_array_is_terminal=res['overall_is_terminal'], largest_array_start=o['start'] if o else np.nan, largest_array_end=o['end'] if o else np.nan,
                   largest_array_size=o['size'] if o else np.nan)
        calls.append(row)
        for a in res['arrays']:
            allrows.append(dict(genome=name, sequence=seq, **a))
    return pd.DataFrame(calls), pd.DataFrame(allrows)


def compare_calls(calls, pub, scan, inv):
    """Add, to the array calls of the REFERENCE genome, the published (lifted) interval, the satellite-scan interval and the relation to the inversion of the same chromosome.
    Sequences without an LG label (unplaced scaffolds, lg = NaN) are left alone: NaN must never be matched against the tables (a NaN row in the scan table once crashed the step)."""
    calls = calls.copy()
    for col in ('published_start', 'published_end', 'scan_start', 'scan_end', 'inversion_start', 'inversion_end'):
        calls[col] = np.nan
    calls['array_vs_published'] = ''
    calls['relation_to_inversion'] = ''
    pub_by = pub[pub.lg.notna()].drop_duplicates('lg').set_index('lg') if pub is not None and len(pub) else None
    scan_by = scan[scan.lg.notna()].drop_duplicates('lg').set_index('lg') if scan is not None and len(scan) else None
    for i, r in calls.iterrows():
        if pd.isna(r.lg) or pd.isna(r.start):
            continue
        if pub_by is not None and r.lg in pub_by.index:
            p_ = pub_by.loc[r.lg]
            ps, pe = (p_.gt3_start, p_.gt3_end) if pd.notna(p_.gt3_start) else (p_.get('bounded_start', np.nan), p_.get('bounded_end', np.nan))
            if pd.notna(ps):
                gap = max(0.0, max(r.start, ps) - min(r.end, pe))
                calls.loc[i, ['published_start', 'published_end']] = [ps, pe]
                calls.loc[i, 'array_vs_published'] = 'overlaps' if gap == 0 else f'{gap / 1e6:.2f} Mb apart'
        if scan_by is not None and r.lg in scan_by.index:
            calls.loc[i, ['scan_start', 'scan_end']] = [scan_by.loc[r.lg, 'start'], scan_by.loc[r.lg, 'end']]
        if r.lg in inv:
            calls.loc[i, ['inversion_start', 'inversion_end']] = [inv[r.lg][0], inv[r.lg][1]]
            calls.loc[i, 'relation_to_inversion'] = cen.relation_to_inversion(r.start, r.end, inv[r.lg][0], inv[r.lg][1])
    return calls


def arrays_step(proj, align_dir=None):
    """centromeres arrays: the centromere is the largest array of the satellite family (config arrays.family) more than 500 kb from either chromosome end, at feature resolution.
    Compared with the published lifted intervals and the satellite scan, and classified against each inversion (pericentric / paracentric)."""
    cfg = cfg_of(); ref = proj.genome
    with proj.run('centromeres', 'arrays', params=dict(cfg=cfg, excluded_regions=excluded_regions(proj))) as run:
        pub = scan = None
        try:
            pub = pd.read_csv(latest_dir(proj, 'centromeres', 'lift') + 'published_centromeres_lifted.tsv', sep='\t')
        except FileNotFoundError:
            pass
        try:
            s_ = pd.read_csv(latest_dir(proj, 'centromeres', 'scan') + 'centromeres.tsv', sep='\t'); scan = s_[s_.genome == ref]
        except FileNotFoundError:
            pass
        inv = {}
        try:
            cdir = latest_dir(proj, 'inversion_breakpoints', 'call')
            for lg in proj.lg.values():
                f = f'{cdir}{lg}_consensus_breakpoints.tsv'
                if os.path.exists(f) and os.path.getsize(f) > 5:
                    c = pd.read_csv(f, sep='\t')
                    if 'confidence' in c:
                        c = c[c.confidence.isin(['high', 'medium'])]
                    if len(c): inv[lg] = (c.iloc[0].left_breakpoint, c.iloc[0].right_breakpoint)
        except FileNotFoundError:
            pass
        out = []
        for key, g in cfg['genomes'].items():
            spec = g.get('annotation')
            if not spec:
                continue
            try:
                fa = fasta_path(proj, key)
                apath, nmap = annotation_for(proj, key) if spec == 'registered' else (os.path.expanduser(spec), None)
            except FileNotFoundError as e:
                print(f'[arrays] {key}: skipped ({e})'); continue
            calls, allarr = arrays_for_genome(key, fa, apath, cfg, nmap, excluded_regions(proj) if key == ref else None)
            if not len(calls):
                continue
            calls['lg'] = calls.sequence.map(proj.lg) if key == ref else np.nan
            if key != ref:
                print(f'[arrays] {key}: arrays are reported in the genome\'s own coordinates (lift to the reference is done by `centromeres scan` for annotated non-reference genomes)')
            if key == ref:
                calls = compare_calls(calls, pub, scan, inv)                                # published lift, satellite scan and inversion
            out.append(calls)
            allarr.to_csv(f'{run.dir}{key}_satellite_arrays.tsv', sep='\t', index=False)
            show = calls[['sequence', 'lg', 'start', 'end', 'size', 'end_distance', 'n_internal_arrays', 'ambiguous', 'largest_array_is_terminal']].copy()
            for c in ('start', 'end', 'size', 'end_distance'): show[c] = (show[c] / 1e6).round(3)
            print(f'\n== {key}: largest {cfg["arrays"]["family"]} array > {cfg["arrays"]["min_end_distance"] / 1e3:.0f} kb from both ends (Mb) ==')
            print(show.to_string(index=False))
            odd = calls[calls.largest_array_is_terminal == True]
            if len(odd): print('chromosomes whose largest array is TERMINAL and therefore not chosen:', sorted(odd.lg.dropna()))
        if out:
            pd.concat(out, ignore_index=True).to_csv(run.dir + 'centromeres_arrays.tsv', sep='\t', index=False)
