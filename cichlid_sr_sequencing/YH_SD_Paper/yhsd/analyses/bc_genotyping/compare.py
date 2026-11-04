"""Compare ORIGINAL (rerun) vs V2 calls on the same fish and grid; write tables and one summary figure."""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def _switches(M, chrom):
    out = {}
    for ch in pd.unique(chrom):
        idx = np.where(chrom == ch)[0]
        out[ch] = (M[idx] != np.roll(M[idx], 1, axis=0))[1:]    # (n-1, fish)
    return out


def summarize(calls, chrom, pos, label, step):
    rows = []
    for ch in pd.unique(chrom):
        idx = np.where(chrom == ch)[0]; M = calls[idx]; p = pos[idx]
        sw = (np.diff(M, axis=0) != 0)
        n_sw = sw.sum(0)
        term, inter = [], []
        L = p.max() + step
        for j in range(M.shape[1]):
            w = np.where(sw[:, j])[0]
            if len(w):
                term += [(p[w[0]] + step) / 1e6, (L - p[w[-1]] - step) / 1e6]
                inter += list(np.diff(p[w]) / 1e6)
        rows.append(dict(method=label, chrom=ch, state1_freq=M.mean(), cM=100 * n_sw.mean(),
                         frac_ge2_switch=(n_sw >= 2).mean(), frac_terminal_lt2Mb=np.mean(np.array(term) < 2) if term else np.nan,
                         median_interswitch_Mb=np.median(inter) if inter else np.nan))
    return pd.DataFrame(rows)


def sex_scan(calls, chrom, pos, samples, pheno, ch, exclude=()):
    ph = pheno.set_index('samplename').sex_100
    keep = [i for i, s in enumerate(samples) if ph.get(s, '-') in ('M', 'F') and s not in exclude]
    y = np.array([ph[samples[i]] == 'M' for i in keep])
    idx = np.where(chrom == ch)[0]
    res = []
    for k in idx:
        x = calls[k, keep] == 1
        t = [[(x & y).sum(), (x & ~y).sum()], [(~x & y).sum(), (~x & ~y).sum()]]
        res.append((pos[k], -np.log10(stats.fisher_exact(t)[1]), int((x != y).sum())))
    return pd.DataFrame(res, columns=['pos', 'nlp', 'discordant']), len(keep)


def main(proj):
    cfg = proj.cfg
    o = pd.read_csv(proj.dirs['original'] + 'OutputGenotypes_original.csv', index_col=0)
    v = pd.read_csv(proj.dirs['v2'] + 'OutputPosteriors.csv')
    fish = [c for c in o.columns[2:] if c in v.columns]
    m = o[['Chromosome', 'Position'] + fish].merge(v[['Chromosome', 'Position'] + fish], on=['Chromosome', 'Position'], suffixes=('_o', '_v'))
    chrom, pos = m.Chromosome.values, m.Position.values
    O = (m[[f + '_o' for f in fish]].values >= 0.5).astype(int)
    V = m[[f + '_v' for f in fish]].values
    Vh = (V >= 0.5).astype(int)
    step = cfg['v2']['grid_step']
    lg = proj.lg
    conc = pd.DataFrame([dict(chrom=ch, LG=lg[ch], agreement=float((O[chrom == ch] == Vh[chrom == ch]).mean()),
                              v2_uncertain_frac=float(((V[chrom == ch] > 0.1) & (V[chrom == ch] < 0.9)).mean()))
                         for ch in pd.unique(chrom)])
    s_o = summarize(O, chrom, pos, 'original', step); s_v = summarize(Vh, chrom, pos, 'v2', step)
    tab = pd.concat([s_o, s_v]); tab['LG'] = tab.chrom.map(lg)
    out = proj.dirs['comparison']
    conc.to_csv(out + 'agreement_by_chrom.tsv', sep='\t', index=False)
    tab.to_csv(out + 'recombination_and_distortion_by_chrom.tsv', sep='\t', index=False)
    # chr10 sex scan, both methods, with and without the two post-mortem phenotyped fish
    ch10 = proj.chrom_for_lg(cfg['sex_lg'])
    pheno = proj.load_pheno()
    rows, scans = [], {}
    for label, calls in (('original', O), ('v2', Vh)):
        for ex_label, ex in (('all', ()), ('excl_566_760', ('MCYHBC1-566-1', 'MCYHBC1-760-1'))):
            sc, n = sex_scan(calls, chrom, pos, fish, pheno, ch10, ex)
            scans[(label, ex_label)] = sc
            b = sc.loc[sc.nlp.idxmax()]
            rows.append(dict(method=label, fish=ex_label, n=n, peak_pos=int(b.pos), max_nlp=float(b.nlp), discordant_at_peak=int(b.discordant)))
    pd.DataFrame(rows).to_csv(out + 'chr10_sex_scan_summary.tsv', sep='\t', index=False)
    print(pd.DataFrame(rows).to_string(index=False)); print(conc.agreement.describe())
    # figure
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    a, b = s_o.set_index('chrom'), s_v.set_index('chrom')
    ax[0].scatter(a.state1_freq, b.state1_freq.reindex(a.index)); ax[0].plot([.4, .8], [.4, .8], 'k--', lw=.7)
    for ch in a.index:
        if abs(a.state1_freq[ch] - .5) > .08 or abs(b.state1_freq[ch] - .5) > .08:
            ax[0].annotate(lg[ch], (a.state1_freq[ch], b.state1_freq[ch]))
    ax[0].set_xlabel('YH-state frequency, original'); ax[0].set_ylabel('v2'); ax[0].set_title('Transmission ratio by chromosome')
    ax[1].scatter(a.cM, b.cM.reindex(a.index)); ax[1].plot([20, 80], [20, 80], 'k--', lw=.7)
    ax[1].set_xlabel('map length cM, original'); ax[1].set_ylabel('v2'); ax[1].set_title('Chromosome map length')
    for (lab, ex), sc in scans.items():
        if ex == 'all':
            ax[2].plot(sc.pos / 1e6, sc.nlp, label=lab)
    ax[2].set_xlabel(f'{lg[ch10]} position (Mb)'); ax[2].set_ylabel('-log10 p (sex)'); ax[2].legend(frameon=False)
    plt.tight_layout(); plt.savefig(out + 'comparison_summary.png', dpi=200)
