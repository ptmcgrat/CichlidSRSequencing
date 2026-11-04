"""centromeres analysis: python -m yhsd.cli centromeres scan
For every genome in config/centromeres.yaml that has a repeat annotation: satellite density per window, centromere call per chromosome (confidence-graded), assembly gaps,
profile plots, and (for non-reference genomes) the called interval lifted to the reference through the assembly alignments of inversion_breakpoints (bounded by the nearest alignable flanks)."""
import os
import numpy as np
import pandas as pd
import yaml
from ...common_light import config_path
from ...genomes_registry import fasta_path, annotation_for
from ...runner import latest_dir
from ..inversion_breakpoints import paf as P
from . import cen


def cfg_of():
    return yaml.safe_load(open(config_path('centromeres.yaml')))


def scan_genome(name, fa, annotation, cfg, name_map=None):
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
        r = cen.call_centromere(dens, cfg['window'], cfg['smooth'], cfg['confidence'], gcov)
        r.update(genome=name, sequence=seq, length=length, n_satellite_features=int(len(sub)),
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
            calls, prof = scan_genome(key, fa, apath, cfg, nmap)
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
                                                    'n_satellite_features', 'top_families', 'ref_start', 'ref_end', 'lift_status']].to_csv(run.dir + 'centromeres.tsv', sep='\t', index=False)
