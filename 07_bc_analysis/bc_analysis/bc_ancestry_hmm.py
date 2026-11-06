"""
bc_ancestry_hmm.py
Site-level ancestry HMM for low-coverage backcross (BC) offspring of an F1 sire.

State 0 = offspring inherited the MC haplotype from the F1 sire (homozygous MC allele at informative sites)
State 1 = offspring inherited the YH-derived haplotype (het: MC allele + YH allele at informative sites)

ALLELE POLARITY: counts are (n_mc, n_yh) = reads for the MC-fixed allele and the YH-specific allele, NOT ref/alt.
With a non-MC reference (e.g. Mzebra) the YH allele is the REFERENCE allele at ~half the sites, and mapping bias then
raises the YH fraction in hets above 0.5 instead of lowering it. Pass `yh_is_ref` (bool per site) so the het emission is
P(YH read | het) = q if yh_is_ref else 1-q, with q = P(reference-allele read | het) fitted per sample.
If yh_is_ref is None the old behaviour is used (constant p1 = P(YH read | het)).

Design choices (and why):
  * Observations are per-SITE read counts (ref, alt), not one observation per read.
    PCR duplicates / mismapped read clusters therefore cannot be counted as independent evidence;
    counts are additionally capped at `cap` per allele.
  * Transition probabilities depend on the physical distance between observed sites
    (rate `rho` crossovers per bp), so the crossover penalty is a few nats, not ~50 nats.
  * Emission error (eps) and het alt-fraction (p1) are re-estimated per sample by EM
    because error and reference bias vary with coverage and library quality.
  * Output is the posterior P(state = 1) at every site AND at a regular grid of zero-count
    "pseudo-sites", so the grid values are true marginal posteriors rather than interpolations.
Only numpy is required.
"""
import numpy as np
from math import exp, log


def _prep(pos, ref, alt, cap):
    ref = np.minimum(ref, cap)
    alt = np.minimum(alt, cap)
    return ref.astype(float), alt.astype(float)


def forward_backward(pos, ref, alt, eps=0.005, p1=0.45, rho=1e-8, cap=2, yh_is_ref=None):
    """Return P(state=1) at each entry of `pos` (sorted, int bp). ref=n_mc, alt=n_yh. Zero-count entries allowed."""
    n = len(pos)
    ref, alt = _prep(pos, ref, alt, cap)
    p1v = p1 if yh_is_ref is None else np.where(yh_is_ref, p1, 1.0 - p1)
    l0 = ref * log(1 - eps) + alt * log(eps)
    l1 = ref * np.log(1 - p1v) + alt * np.log(p1v)
    m = np.maximum(l0, l1)
    e0 = np.exp(l0 - m).tolist()
    e1 = np.exp(l1 - m).tolist()
    d = np.diff(pos).astype(float)
    sw = (0.5 * (1.0 - np.exp(-2.0 * rho * d))).tolist()  # symmetric 2-state switch prob

    a0 = [0.0] * n
    a1 = [0.0] * n
    x0, x1 = 0.5 * e0[0], 0.5 * e1[0]
    s = x0 + x1
    a0[0], a1[0] = x0 / s, x1 / s
    for i in range(1, n):
        t = sw[i - 1]
        y0 = (a0[i - 1] * (1 - t) + a1[i - 1] * t) * e0[i]
        y1 = (a1[i - 1] * (1 - t) + a0[i - 1] * t) * e1[i]
        s = y0 + y1
        a0[i], a1[i] = y0 / s, y1 / s
    b0 = [1.0] * n
    b1 = [1.0] * n
    post = np.empty(n)
    post[n - 1] = a1[n - 1]
    for i in range(n - 2, -1, -1):
        t = sw[i]
        z0 = (1 - t) * e0[i + 1] * b0[i + 1] + t * e1[i + 1] * b1[i + 1]
        z1 = (1 - t) * e1[i + 1] * b1[i + 1] + t * e0[i + 1] * b0[i + 1]
        s = z0 + z1
        b0[i], b1[i] = z0 / s, z1 / s
        q0, q1 = a0[i] * b0[i], a1[i] * b1[i]
        post[i] = q1 / (q0 + q1)
    return post


def fit_sample(pos, ref, alt, rho=1e-8, cap=2, n_iter=3, eps0=0.005, p10=0.45, flags=None):
    """EM for eps and p1 (per sample; pass lists with one entry per chromosome).
    Without flags p1 = P(YH read | het). With flags (yh_is_ref per site) the fitted parameter is
    q = P(reference-allele read | het) and the returned p1 is that q (start value 0.55)."""
    eps = eps0
    p1 = p10 if flags is None else 0.55
    for _ in range(n_iter):
        num_e = den_e = num_p = den_p = 0.0
        for k, (P, R, A) in enumerate(zip(pos, ref, alt)):
            fl = None if flags is None else flags[k]
            post = forward_backward(P, R, A, eps, p1, rho, cap, fl)
            Rc, Ac = np.minimum(R, cap), np.minimum(A, cap)
            n = Rc + Ac
            num_e += ((1 - post) * Ac).sum(); den_e += ((1 - post) * n).sum()
            if fl is None:
                num_p += (post * Ac).sum()
            else:
                num_p += (post * np.where(fl, Ac, Rc)).sum()     # reference-allele reads among hets
            den_p += (post * n).sum()
        eps = float(np.clip(num_e / max(den_e, 1), 0.002, 0.05))
        p1 = float(np.clip(num_p / max(den_p, 1), 0.30, 0.50 if flags is None else 0.70))
    return eps, p1


def with_grid(pos, ref, alt, length, step=100_000, offset=None, flag=None):
    """Merge observed sites with zero-count grid points; return merged arrays + mask of grid entries.
    offset=None -> step//2 (simulation default); offset=0 reproduces the original 0,100000,... grid."""
    grid = np.arange(step // 2 if offset is None else offset, length, step)
    allpos = np.concatenate([pos, grid])
    allref = np.concatenate([ref, np.zeros(len(grid), dtype=ref.dtype)])
    allalt = np.concatenate([alt, np.zeros(len(grid), dtype=alt.dtype)])
    isgrid = np.concatenate([np.zeros(len(pos), bool), np.ones(len(grid), bool)])
    o = np.argsort(allpos, kind='stable')
    if flag is not None:
        allflag = np.concatenate([np.asarray(flag, bool), np.zeros(len(grid), bool)])
        return allpos[o], allref[o], allalt[o], isgrid[o], allflag[o]
    return allpos[o], allref[o], allalt[o], isgrid[o]


def call_breakpoints(pos, post, hi=0.95):
    """Hard calls from posteriors plus a breakpoint interval for each switch.
    Interval = (last position confidently in the old state, first position confidently in the new state)."""
    state = (post > 0.5).astype(int)
    out = []
    sw = np.where(np.diff(state) != 0)[0]
    for k in sw:
        old = state[k]
        left = pos[k]
        for j in range(k, -1, -1):
            if (old == 1 and post[j] >= hi) or (old == 0 and post[j] <= 1 - hi):
                left = pos[j]; break
        right = pos[k + 1]
        for j in range(k + 1, len(pos)):
            if (old == 0 and post[j] >= hi) or (old == 1 and post[j] <= 1 - hi):
                right = pos[j]; break
        out.append((int(pos[k]), int(pos[k + 1]), int(left), int(right)))
    return state, out
