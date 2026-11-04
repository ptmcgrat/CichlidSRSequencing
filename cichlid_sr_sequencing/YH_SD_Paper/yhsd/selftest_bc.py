"""Self-test with synthetic data; needs only numpy/pandas/scipy/matplotlib. `python -m yhsd selftest`"""
import os, sys, tempfile, types
import numpy as np
import pandas as pd


def _ok(name, cond, extra=''):
    print(('PASS ' if cond else 'FAIL ') + name, extra)
    if not cond:
        raise SystemExit(1)


def main():
    from .bc_ancestry_hmm import forward_backward, fit_sample, with_grid, call_breakpoints
    from .analyses.bc_genotyping.original import stretch, viterbi_numpy, align_rows
    from .analyses.bc_genotyping import compare
    from .analyses.bc_genotyping.v2 import hmm_sample
    rng = np.random.default_rng(11)

    # 1. allele-aware HMM recovers a breakpoint and the reference-bias parameter with mixed polarity
    n = 8000; pos = np.sort(rng.choice(40_000_000, n, replace=False)); st = (pos > 15_000_000).astype(int)
    flag = rng.random(n) < .5; q = 0.56
    u = rng.poisson(.5, n); yh = rng.binomial(u, np.where(st == 1, np.where(flag, q, 1 - q), 0.005)); obs = u > 0
    ap, ar, aa, g, fl = with_grid(pos[obs], (u - yh)[obs], yh[obs], 40_000_000, flag=flag[obs])
    eps, qh = fit_sample([ap], [ar], [aa], flags=[fl]); post = forward_backward(ap, ar, aa, eps, qh, yh_is_ref=fl)
    bp = call_breakpoints(ap, post)[1]
    _ok('hmm breakpoint within 100 kb', len(bp) == 1 and abs(bp[0][0] - 15_000_000) < 100_000, str(bp))
    _ok('hmm recovers q', abs(qh - q) < 0.03, f'q_hat={qh:.3f}')

    # 2. original-approach helpers
    o, c, p = stretch(np.array([[3, 0], [2, 1], [0, 4]]), np.array(['x'] * 3), np.array([1, 2, 3]), rng)
    _ok('stretch read count', len(o) == 10 and (o == 1).sum() == 5)
    s0 = (np.arange(6000) >= 3000).astype(int); ob = (rng.random(6000) < np.where(s0 == 1, .4, .03)).astype(np.uint8)
    path = viterbi_numpy(ob, [[.97, .03], [.6, .4]], 1e-22)
    _ok('numpy viterbi', (path == s0).mean() > 0.99)
    a = dict(chrom=np.array(['a', 'a', 'b']), pos=np.array([1, 2, 3])); b = dict(chrom=np.array(['a', 'b']), pos=np.array([2, 3]))
    ia, ib = align_rows(a, b); _ok('align_rows', list(ia) == [1, 2] and list(ib) == [0, 1])

    # 3. v2 per-sample wrapper + compare on a toy 2-chromosome BC panel
    cfg = dict(v2=dict(rho=1e-8, cap=2, grid_step=100000), sex_lg='LG10')
    lengths = {'chrA': 20_000_000, 'chrB': 30_000_000}
    fish = [f'F{i}' for i in range(40)]
    tmp = tempfile.mkdtemp() + '/'
    for d in ('original', 'v2', 'comparison'):
        os.makedirs(tmp + d)
    ocol, vcol, sex = {}, {}, {}
    for f in fish:
        male = rng.random() < .5; sex[f] = 'M' if male else 'F'
        rows = []
        for ch, L in lengths.items():
            m = int(L / 2000); ps = np.sort(rng.choice(L, m, replace=False)); fgl = rng.random(m) < .5
            s = np.full(m, 1 if (ch == 'chrB' and male) else 0)
            uu = rng.poisson(.55, m); y = rng.binomial(uu, np.where(s == 1, np.where(fgl, .56, .44), .005))
            rows.append(pd.DataFrame(dict(chrom=ch, pos=ps, mc=uu - y, yh=y, yh_is_ref=fgl)))
        cnt = pd.concat(rows); cnt = cnt[(cnt.mc + cnt.yh) > 0]
        res, _ = hmm_sample(cnt, lengths, cfg['v2'])
        vcol[f] = np.concatenate([res[ch]['post'] for ch in lengths]); ocol[f] = np.round(vcol[f])
    grid = pd.concat([pd.DataFrame({'Chromosome': ch, 'Position': np.arange(0, L, 100000)}) for ch, L in lengths.items()], ignore_index=True)
    pd.concat([grid, pd.DataFrame(ocol)], axis=1).to_csv(tmp + 'original/OutputGenotypes_original.csv')
    pd.concat([grid, pd.DataFrame(vcol)], axis=1).to_csv(tmp + 'v2/OutputPosteriors.csv', index=False)
    pheno = pd.DataFrame(dict(samplename=fish, sex_100=[sex[f] for f in fish]))
    proj = types.SimpleNamespace(cfg=cfg, lg={'chrA': 'LG1', 'chrB': 'LG10'}, chrom_for_lg=lambda x: 'chrB',
                                 dirs=dict(original=tmp + 'original/', v2=tmp + 'v2/', comparison=tmp + 'comparison/'),
                                 load_pheno=lambda: pheno)
    compare.main(proj)
    sc = pd.read_csv(tmp + 'comparison/chr10_sex_scan_summary.tsv', sep='\t')
    _ok('compare finds the planted sex locus', (sc.discordant_at_peak == 0).all() and (sc.max_nlp > 8).all())

