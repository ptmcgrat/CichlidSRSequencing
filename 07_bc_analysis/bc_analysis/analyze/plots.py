"""The figures of the analyses. matplotlib is imported when a figure is drawn: without it the commands still write every table and say so."""
import numpy as np
import pandas as pd


def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def plot_marey(mb, lens, masked, terminal_bp, path, bin_bp=100_000):
    plt = _plt()
    chroms = list(dict.fromkeys(mb.chrom))
    ncol = 6
    nrow = int(np.ceil(len(chroms) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 2.6 * nrow), squeeze=False)
    for ax, c in zip(axes.ravel(), chroms):
        d = mb[mb.chrom == c]
        ax.plot(d.Mb, d.cM, color='#d62728', lw=1.4)
        L = d.Mb.max()
        ax.axvspan(0, terminal_bp / 1e6, color='#e8e8e8', alpha=0.8)
        ax.axvspan(L - terminal_bp / 1e6, L, color='#e8e8e8', alpha=0.8)
        for (cc, b) in masked:
            if cc == c:
                ax.axvline(b * bin_bp / 1e6 + bin_bp / 2e6, color='#1f77b4', lw=0.4, alpha=0.5)
        row = lens[lens.chrom == c].iloc[0]
        ax.set_title(f'LG{row.LG} {row.cM:.0f} cM, {row.Mb:.0f} Mb' if str(row.LG).isdigit() else f'{row.LG} {row.cM:.0f} cM, {row.Mb:.0f} Mb', fontsize=8, loc='left')
        ax.tick_params(labelsize=6)
    for ax in axes.ravel()[len(chroms):]:
        ax.axis('off')
    fig.supxlabel('Mb (grey: the first and last 2 Mb, low confidence; blue ticks: masked bins)')
    fig.supylabel('cumulative cM')
    plt.tight_layout()
    plt.savefig(path, dpi=110)
    plt.close()


def plot_rate(profile, region, null_rates, path):
    plt = _plt()
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.2), gridspec_kw={'width_ratios': [3, 1.2]})
    ax[0].bar(profile.start_Mb + 0.5, profile.rate, width=0.9, color=['#d62728' if region['start_Mb'] <= m < region['end_Mb'] else '#777777' for m in profile.start_Mb + 0.5])
    ax[0].errorbar(profile.start_Mb + 0.5, profile.rate, yerr=[profile.rate - profile.ci_lo, profile.ci_hi - profile.rate], fmt='none', ecolor='k', lw=0.6, capsize=1)
    ax[0].axhline(region['background_cM_per_Mb'], color='#2ca02c', ls=':', label=f"background of the chromosome ({region['background_cM_per_Mb']:.1f})")
    ax[0].set_xlabel('Mb'); ax[0].set_ylabel('cM per Mb (95% CI over fish)'); ax[0].legend(fontsize=8)
    ax[0].set_title(f"{region['chrom']} {region['start_Mb']:.2f}-{region['end_Mb']:.2f} Mb: {region['rate_cM_per_Mb']:.1f} cM/Mb, {region['ratio']:.1f}x the rest (CI {region['ci_lo']:.1f}-{region['ci_hi']:.1f})", loc='left', fontsize=9)
    if null_rates is not None and len(null_rates):
        ax[1].hist(null_rates, bins=30, color='#bbbbbb'); ax[1].axvline(region['rate_cM_per_Mb'], color='#d62728'); ax[1].set_xlabel('cM/Mb of same-width windows elsewhere'); ax[1].set_title('the genome', fontsize=9, loc='left')
    plt.tight_layout(); plt.savefig(path, dpi=125); plt.close()


def plot_sexscan(scan, path):
    plt = _plt()
    fig, ax = plt.subplots(figsize=(14, 3.8))
    off, ticks = 0.0, []
    for i, (c, g) in enumerate(scan.groupby('chrom', sort=False)):
        ax.scatter(off + g.pos / 1e6, g.nlogp, s=3, color=['#444444', '#8899aa'][i % 2])
        ticks.append((off + g.pos.max() / 2e6, str(g.LG.iloc[0]).replace('LG', '')))
        off += g.pos.max() / 1e6 + 3
    ax.set_xticks([t for t, _ in ticks]); ax.set_xticklabels([l for _, l in ticks], fontsize=7); ax.set_ylabel('-log10 p (sex vs state)'); ax.set_xlabel('linkage group')
    plt.tight_layout(); plt.savefig(path, dpi=125); plt.close()


def plot_finemap(res, calls, chrom, sex, path, flank_bp=1_500_000):
    plt = _plt()
    rec, lim, prof = res['recombinants'], res['limits'], res['profile']
    lo, hi = res['locus_bp'] - flank_bp, res['locus_bp'] + flank_bp
    fish = list(dict.fromkeys(rec[(rec.boundary_bp >= lo) & (rec.boundary_bp <= hi)].fish)) + [f for f in res['discordant'] if f not in set(rec.fish)]
    fig, ax = plt.subplots(2, 1, figsize=(11, 3 + 0.3 * len(fish)), gridspec_kw={'height_ratios': [1, max(2, 0.28 * len(fish))]}, sharex=True)
    m = (prof.pos >= lo) & (prof.pos <= hi)
    ax[0].step(prof.pos[m] / 1e6, prof.discordant[m], where='post', color='k'); ax[0].set_ylabel('fish discordant\nwith their sex')
    pos = calls.pos(chrom); ix = calls.idx(chrom)
    names = {f: i for i, f in enumerate(calls.fish)}
    for y, f in enumerate(fish):
        s = calls.X[ix][:, names[f]]
        for k in range(len(pos)):
            if lo <= pos[k] < hi:
                ax[1].barh(y, calls.bin_bp / 1e6, left=pos[k] / 1e6, color='#d62728' if s[k] == 1 else '#9aa0a6', height=0.8, edgecolor='none')
        sx = 'M' if sex[names[f]] == 1 else 'F'
        ax[1].text(lo / 1e6 - 0.01, y, f"{f} {sx}" + ('  (excluded)' if f in res['excluded'] else ''), ha='right', va='center', fontsize=6, color='#b22222' if f in res['excluded'] else 'k')
    if len(lim):
        for a_ in ax:
            a_.axvspan(lim.left_bp.iloc[0] / 1e6, lim.right_bp.iloc[0] / 1e6, color='#7fbf7f', alpha=0.35, zorder=0)
    ax[1].set_ylim(len(fish), -1); ax[1].set_yticks([]); ax[1].set_xlim(lo / 1e6, hi / 1e6); ax[1].set_xlabel('Mb (red = the YH-carrying state, grey = MC; green = the limits)')
    plt.tight_layout(); plt.subplots_adjust(left=0.16); plt.savefig(path, dpi=125); plt.close()


def plot_crossovers(rows, lod_profiles, bounds, path):
    """rows: list of dicts (fish, role, best, best_before or None); lod_profiles: {fish: (x_after, lod_after, x_before or None, lod_before or None)}; bounds: (left_bp, right_bp) of the interval."""
    plt = _plt()
    n = len(rows)
    fig, ax = plt.subplots(2, n, figsize=(3.6 * n, 6.4), squeeze=False)
    for j, r in enumerate(rows):
        xa, la, xb, lb = lod_profiles[r['fish']]
        z = 60_000
        a = ax[0, j]
        if xb is not None:
            m = (xb > r['best'] - z) & (xb < r['best'] + z); a.plot(xb[m] / 1e6, lb[m], color='#9a9a9a', label='before')
        m = (xa > r['best'] - z) & (xa < r['best'] + z); a.plot(xa[m] / 1e6, la[m], color='#d62728', label='after')
        a.axhline(-1, color='k', ls=':', lw=0.8); a.set_ylim(-8, 0.4); a.set_title(f"{r['fish']} ({r['role']}) {r['best'] / 1e6:.3f} Mb", fontsize=9, loc='left'); a.legend(fontsize=7)
        ax[1, j].axis('off')
    a = fig.add_subplot(2, 1, 2)
    for i, r in enumerate(rows):
        a.plot([r['best'] / 1e6], [i], 'o', color='#1f77b4' if r['role'] == 'left' else '#d62728')
        if r.get('best_before') is not None:
            a.plot([r['best_before'] / 1e6], [i - 0.2], 'o', color='#aaaaaa', ms=4)
        a.text(r['best'] / 1e6, i + 0.15, r['fish'], fontsize=7)
    a.axvspan(bounds[0] / 1e6, bounds[1] / 1e6, color='#7fbf7f', alpha=0.4); a.set_yticks([]); a.set_xlabel('Mb: the crossovers and the interval they give (grey = before)')
    plt.tight_layout(); plt.savefig(path, dpi=125); plt.close()


def plot_trd(calls, sc, band, path):
    plt = _plt()
    fig, ax = plt.subplots(figsize=(14, 3.8))
    off, ticks = 0.0, []
    meta = calls.meta
    for i, c in enumerate(calls.chroms):
        ix = calls.idx(c); pos = calls.pos(c) / 1e6
        ax.scatter(off + pos, sc['f'][ix], s=2, color=['#444444', '#8899aa'][i % 2])
        ticks.append((off + pos.max() / 2, str(calls.lg_of(c)).replace('LG', ''))); off += pos.max() + 3
    ax.axhline(0.5, color='k', lw=0.6); ax.axhspan(0.5 - band, 0.5 + band, color='#ffe9a8', alpha=0.6, zorder=0)
    ax.set_xticks([t for t, _ in ticks]); ax.set_xticklabels([l for _, l in ticks], fontsize=7); ax.set_ylabel('share with the YH state'); ax.set_ylim(0.2, 0.85)
    plt.tight_layout(); plt.savefig(path, dpi=125); plt.close()
