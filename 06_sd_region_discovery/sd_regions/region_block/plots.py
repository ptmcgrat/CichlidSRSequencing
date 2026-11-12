"""Figures of the region block analysis (matplotlib, Agg). Each function writes one PNG and returns its path."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

MALE, FEMALE = '#1f77b4', '#d6604d'


def _sex_colors(sex):
    return [MALE if s == 'M' else FEMALE for s in sex]


def het_heatmap(path, keys, het_win, names, sex, chrom, region, size, core=None, roh_segments=None, title=''):
    """Samples (rows, males above females) x windows of one chromosome (columns), coloured by heterozygosity, and the mean by sex above. A homozygous stretch is a dark run in a row: a run of
    homozygosity. The region is marked by vertical lines."""
    cols = [i for i, k in enumerate(keys) if k[0] == chrom]
    x = np.array([keys[i][1] for i in cols]) * size / 1e6
    H = het_win[cols].T                                                  # samples x windows
    sex = np.asarray(sex)
    reg = [i for i, w in enumerate(x) if region[1] / 1e6 <= w + size / 2e6 < region[2] / 1e6]
    mean_in_region = np.nanmean(H[:, reg], axis=1) if reg else np.zeros(len(names))
    order = []
    for s in ('M', 'F'):
        idx = [i for i in range(len(names)) if sex[i] == s]
        order += sorted(idx, key=lambda i: (np.nan_to_num(mean_in_region[i], nan=9), names[i]))
    fig, (a0, a1) = plt.subplots(2, 1, figsize=(13, 0.28 * len(names) + 3.4), sharex=True, gridspec_kw=dict(height_ratios=[1, max(3, 0.22 * len(names))]))
    for s, col in (('M', MALE), ('F', FEMALE)):
        sel = sex == s
        if sel.any():
            m = np.nanmean(H[sel], axis=0)
            a0.plot(x, m, color=col, lw=1.6, label=f'{"males" if s == "M" else "females"} (n={int(sel.sum())})')
    a0.set_ylabel('mean het'); a0.legend(frameon=False, ncol=2, loc='lower right', bbox_to_anchor=(1.0, 1.0), fontsize=9); a0.set_title(title or f'Heterozygosity along {chrom}', loc='left', fontweight='bold')
    im = a1.imshow(H[order], aspect='auto', cmap='viridis', vmin=0, vmax=0.45, extent=[x.min(), x.max() + size / 1e6, len(order), 0], interpolation='nearest')
    a1.set_yticks(np.arange(len(order)) + 0.5); a1.set_yticklabels([names[i] for i in order], fontsize=7)
    for t, i in zip(a1.get_yticklabels(), order):
        t.set_color(MALE if sex[i] == 'M' else FEMALE)
    a1.set_xlabel(f'{chrom} (Mb)')
    for a in (a0, a1):
        for v in (region[1] / 1e6, region[2] / 1e6):
            a.axvline(v, color='white' if a is a1 else '#444', lw=1.2, ls='--')
        if core:
            for v in core:
                a.axvline(v / 1e6, color='#ffcc00', lw=1.0, ls=':')
    if roh_segments:
        pos_in_order = {s: r for r, s in enumerate(order)}
        for g in roh_segments:
            if g['chrom'] == chrom and g['sample'] in pos_in_order:
                r = pos_in_order[g['sample']]
                a1.plot([g['start'] / 1e6, g['end'] / 1e6], [r + 0.92, r + 0.92], color='#ff4d4d', lw=1.6)
    cb = fig.colorbar(im, ax=[a0, a1], pad=0.01, fraction=0.025); cb.set_label('heterozygosity per window')
    fig.savefig(path, dpi=140, bbox_inches='tight'); plt.close(fig)
    return path


def kinship_plot(path, G, names, sex, stats=None):
    order = sorted(range(len(names)), key=lambda i: (sex[i], names[i]))
    M = G[np.ix_(order, order)].copy()
    np.fill_diagonal(M, np.nan)                                          # the diagonal (1 + F) would swamp the colour scale; F is in samples.tsv
    fig, ax = plt.subplots(figsize=(8.2, 7.2))
    cmap = plt.get_cmap('RdBu_r').copy(); cmap.set_bad('#f2f2f2')
    im = ax.imshow(M, cmap=cmap, vmin=-0.5, vmax=0.5)
    ax.set_xticks(range(len(order))); ax.set_xticklabels([names[i] for i in order], rotation=90, fontsize=7)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([names[i] for i in order], fontsize=7)
    for tl in (ax.get_xticklabels(), ax.get_yticklabels()):
        for t, i in zip(tl, order):
            t.set_color(MALE if sex[i] == 'M' else FEMALE)
    ax.set_title('Relatedness between fish (genome outside the region)\n' + (stats or ''), loc='left', fontsize=10, fontweight='bold')
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label='relationship coefficient (0.5 = full sibs, 0.25 = half sibs)')
    fig.savefig(path, dpi=140, bbox_inches='tight'); plt.close(fig)
    return path


def pca_plot(path, region_scores, region_var, genome_scores, genome_var, names, sex, region_note='', genome_note=''):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
    for ax, sc, var, note, ttl in ((axes[0], region_scores, region_var, region_note, 'PCA of the SNPs in the region'), (axes[1], genome_scores, genome_var, genome_note, 'PCA of the rest of the genome')):
        ax.scatter(sc[:, 0], sc[:, 1], c=_sex_colors(sex), s=48, edgecolors='k', linewidths=.4)
        for i, nm in enumerate(names):
            ax.annotate(nm.replace('ASG_', '').replace('_1', ''), (sc[i, 0], sc[i, 1]), fontsize=6.5, xytext=(3, 3), textcoords='offset points')
        ax.set_xlabel(f'PC1 ({100 * var[0]:.0f}%)'); ax.set_ylabel(f'PC2 ({100 * var[1]:.0f}%)' if len(var) > 1 else 'PC2')
        ax.set_title(f'{ttl}\n{note}', loc='left', fontsize=10, fontweight='bold'); ax.grid(color='#eee')
    axes[0].scatter([], [], c=MALE, label='male'); axes[0].scatter([], [], c=FEMALE, label='female'); axes[0].legend(frameon=False)
    fig.savefig(path, dpi=140, bbox_inches='tight'); plt.close(fig)
    return path


def ld_plot(path, mats, title='', ncols=3):
    """mats = [(label, r2 matrix, positions)]: one panel each (SNPs in position order), ncols per row."""
    nrows = int(np.ceil(len(mats) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 5.0 * nrows), squeeze=False)
    im = None
    for ax, (lab, r2, pos) in zip(axes.ravel(), mats):
        ax.set_title(lab, fontsize=9.5, fontweight='bold')
        if r2 is None:
            ax.text(0.5, 0.5, 'no SNP varies among these fish\n(no linkage disequilibrium to measure)', ha='center', va='center', fontsize=10, color='#555', transform=ax.transAxes)
            ax.set_xticks([]); ax.set_yticks([])
            continue
        im = ax.imshow(r2, vmin=0, vmax=1, cmap='magma_r', interpolation='nearest')
        ticks = np.linspace(0, len(pos) - 1, 5).astype(int)
        ax.set_xticks(ticks); ax.set_xticklabels([f'{pos[t] / 1e6:.1f}' for t in ticks], fontsize=8); ax.set_yticks(ticks); ax.set_yticklabels([f'{pos[t] / 1e6:.1f}' for t in ticks], fontsize=8)
        ax.set_xlabel('position (Mb)')
    for ax in axes.ravel()[len(mats):]:
        ax.axis('off')
    if im is not None:
        fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02, label='r$^2$')
    fig.suptitle(title, x=0.01, ha='left', fontsize=11, fontweight='bold')
    fig.savefig(path, dpi=130, bbox_inches='tight'); plt.close(fig)
    return path


def diag_plot(path, by_depth, subsamples, full_count, n_full, names=None, sex=None, loo=None):
    """Left: the count of sex-diagnostic SNPs under each rule and minimum depth. Middle: the count in random subsets of fewer fish against all of them. Right (when given): the strict count with each fish left
    out in turn: the fish that contradict the block stand out."""
    three = loo is not None
    fig, axes = plt.subplots(1, 3 if three else 2, figsize=(18.5 if three else 12.5, 4.8))
    a, b = axes[0], axes[1]
    ts = sorted(by_depth)
    w = 0.38
    a.bar(np.arange(len(ts)) - w / 2, [by_depth[t][0] for t in ts], w, label="sexloci's rule (up to 10% exceptions)", color='#6a3d9a')
    a.bar(np.arange(len(ts)) + w / 2, [by_depth[t][1] for t in ts], w, label='strict (no exceptions)', color='#999')
    a.set_xticks(range(len(ts))); a.set_xticklabels([f'depth >= {t}' for t in ts]); a.set_ylabel('male-diagnostic SNPs in the region'); a.legend(frameon=False, fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)
    a.set_title('Does the count depend on the depth cut-off?', loc='left', fontweight='bold', fontsize=10)
    labels, data = [], []
    for (nm, nf), arr in sorted(subsamples.items()):
        labels.append(f'{nm}M + {nf}F'); data.append(arr)
    labels.append(f'all ({n_full})'); data.append(np.array([full_count]))
    rng = np.random.default_rng(0)
    for i, arr in enumerate(data[:-1], start=1):                           # every random draw as a point, and the median
        b.scatter(i + rng.uniform(-0.15, 0.15, len(arr)), arr, s=9, color='#999', alpha=.7, zorder=2)
        b.plot([i - 0.25, i + 0.25], [np.median(arr)] * 2, color='k', lw=2, zorder=3)
    b.scatter([len(data)], [full_count], color='#6a3d9a', zorder=3, s=70); b.set_xticks(range(1, len(data) + 1)); b.set_xticklabels(labels); b.set_xlim(0.4, len(data) + 0.6)
    b.set_ylim(bottom=0); b.set_ylabel('male-diagnostic SNPs (sexloci rule)'); b.set_title('... or on how many fish are analysed?', loc='left', fontweight='bold', fontsize=10)
    if three:
        c = axes[2]
        order = np.argsort(np.array(sex) != 'M')                           # males first, then females
        order = list(sorted(range(len(names)), key=lambda i: (sex[i], names[i])))
        c.bar(range(len(order)), [loo[i, 1] for i in order], color=[MALE if sex[i] == 'M' else FEMALE for i in order])
        c.set_xticks(range(len(order))); c.set_xticklabels([names[i] for i in order], rotation=90, fontsize=7)
        for t, i in zip(c.get_xticklabels(), order):
            t.set_color(MALE if sex[i] == 'M' else FEMALE)
        c.set_ylabel('strict diagnostic SNPs without this fish'); c.set_title('Which fish decide it? (one left out at a time)', loc='left', fontweight='bold', fontsize=10)
    fig.savefig(path, dpi=140, bbox_inches='tight'); plt.close(fig)
    return path


def null_plot(path, rows, region_auc, rs, re_, chrom):
    """The sex AUC of the region PCA on every window of the genome: the region (red line) against the distribution of the others, and the AUC along the genome."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(14, 4.6), gridspec_kw=dict(width_ratios=[1, 2]))
    other = [r for r in rows if not r['overlaps_region']]
    aucs = np.array([r['auc'] for r in other])
    a.hist(aucs, bins=np.linspace(0.5, 1.0, 26), color='#999', edgecolor='white')
    a.axvline(region_auc, color='#d6604d', lw=2.5); a.text(region_auc, a.get_ylim()[1] * 0.92, ' region', color='#d6604d', ha='right' if region_auc > 0.9 else 'left', fontweight='bold')
    a.set_xlabel('AUC for sex of the best of three PCs'); a.set_ylabel('windows'); a.set_title(f'{len(aucs)} windows elsewhere in the genome', loc='left', fontweight='bold', fontsize=10)
    chroms = list(dict.fromkeys(r['chrom'] for r in rows)); offs, tot = {}, 0
    for c in chroms:
        offs[c] = tot; tot += max(r['end'] for r in rows if r['chrom'] == c) + 2_000_000
    for i, c in enumerate(chroms):
        rr = [r for r in rows if r['chrom'] == c]
        b.scatter([offs[c] + (r['start'] + r['end']) / 2 for r in rr], [r['auc'] for r in rr], s=14, color='#d6604d' if c == chrom else ('#555' if i % 2 == 0 else '#999'), zorder=3)
    b.axhline(region_auc, color='#d6604d', lw=1, ls='--')
    b.set_xticks([offs[c] + max(r['end'] for r in rows if r['chrom'] == c) / 2 for c in chroms]); b.set_xticklabels([c.replace('NC_1351', '').replace('.1', '') for c in chroms], fontsize=7, rotation=90)
    b.set_ylabel('AUC for sex'); b.set_title('along the genome (windows of the same length as the region; red = the region\'s chromosome)', loc='left', fontweight='bold', fontsize=10); b.set_ylim(0.45, 1.03)
    fig.savefig(path, dpi=140, bbox_inches='tight'); plt.close(fig)
    return path
