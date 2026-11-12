"""The ancestry HMM of the BC genotyper. Pure numpy; every function is tested on data with a known answer (tests/test_genotyper.py) and the design on simulated crosses (genotype simulate).

Data. One BC fish inherited the MC allele from its MC dam and either the MC or the YH-derived allele from its F1 sire at every informative site, so a site has two classes of reads, MC-allele and YH-allele, and the
fish is in state 0 (homozygous MC: YH-looking reads are errors) or state 1 (heterozygous: about half of the reads are YH, shifted by the reference bias). The Markov chain runs over BINS of `bin_bp` (100 kb) because
the coverage is low and the resolution of the mapping is 100 kb; reads are summed per bin, after a cap per site and per bin so that a cluster of mismapped reads cannot count as many independent observations.

Switches. A crossover is a switch of state. The chain has states (genotype, age), the age being the number of bins since the last switch, up to the `interior` limit A. A switch has the ordinary hazard `rho * bin_bp` when the segment is
old and a hazard smaller by `young` while it is shorter than A: a short INTERIOR segment needs two switches, so it is forbidden unless the reads outweigh a very large penalty; this is the belief that two crossovers a
few hundred kb apart do not happen, made explicit, and it is what removes the flip-flops that a plain distance-based model makes. A TERMINAL segment has only one switch, so it costs an ordinary switch: the start of the chromosome
is given the age A - Wt + 1 (the first segment may end after Wt bins) and the last segment is down-weighted only when it is shorter than Wt bins. The effect of this limit is to STRETCH or DROP: a terminal
segment shorter than Wt bins whose reads are strong is called with the minimum length Wt (the extra bins cost less than the penalty), and one whose reads are weak is not called at all. A penalty of 1e-22 per switch, as in the original script, treats the two cases alike
and misses crossovers within about 1 Mb of a chromosome end when the markers are sparse there (see the simulation study in the README).

Reference bias. Each site has a polarity: the YH allele is the reference allele (a) or the alternative allele (b). In a heterozygote a reference-allele read has probability q > 0.5 (mapping bias towards the reference), so the
YH-read probability is q at `a` sites and 1 - q at `b` sites. With the old MC reference every site was `b`; with Mzebra_GT3 about half are `a`, so counts are kept per polarity.
"""
import numpy as np


def bin_counts(pos, n_mc, n_yh, yh_is_ref, length, bin_bp, site_cap=2, bin_cap=20):
    """Counts per bin as a (4, n_bins) array: YH reads at sites where the YH allele is the reference, MC reads there, YH reads where it is the alternative allele, MC reads there. Each site contributes at most
    `site_cap` reads of each allele and each bin at most `bin_cap` reads in all (the counts of a bin are scaled down to it)."""
    nb = int(np.ceil(length / bin_bp))
    out = np.zeros((4, nb))
    if len(pos) == 0:
        return out
    pos, yh_is_ref = np.asarray(pos), np.asarray(yh_is_ref, bool)
    b = np.minimum(pos // bin_bp, nb - 1).astype(int)
    mc, yh = np.minimum(np.asarray(n_mc), site_cap), np.minimum(np.asarray(n_yh), site_cap)
    for k, (arr, m) in enumerate(((yh, yh_is_ref), (mc, yh_is_ref), (yh, ~yh_is_ref), (mc, ~yh_is_ref))):
        np.add.at(out[k], b[m], arr[m])
    tot = out.sum(0)
    return out * np.minimum(1.0, bin_cap / np.maximum(tot, 1e-12))


def bin_loglik(c, eps, q):
    """Log-likelihood of the reads of each bin in state 0 (homozygous MC) and state 1 (heterozygous)."""
    yh_a, mc_a, yh_b, mc_b = c
    l0 = (yh_a + yh_b) * np.log(eps) + (mc_a + mc_b) * np.log(1 - eps)
    l1 = yh_a * np.log(q) + mc_a * np.log(1 - q) + yh_b * np.log(1 - q) + mc_b * np.log(q)
    return l0, l1


def segment_limits(hp, bin_bp=None):
    """(A, Wt) from the parameters. The chain protects a segment while its age is below A, so the shortest unprotected interior segment is A + 1 bins: A is therefore the minimum interior segment in bins
    minus one. Wt is the minimum terminal segment in bins (None = no limit)."""
    bin_bp = bin_bp or hp['bin_bp']
    A = int(round(hp['min_interior_segment_mb'] * 1e6 / bin_bp)) - 1
    Wt = int(round(hp['min_terminal_segment_mb'] * 1e6 / bin_bp))
    return max(A, 1), (Wt if Wt > 0 else None)


def forward_backward(ll0, ll1, lam, A, young=1e-30, Wt=None, term_factor=1e-30):
    """P(state = 1) of each bin, from a chain whose states are (genotype, bins since the last switch) plus, for the first segment of the chromosome, (genotype, bins since the start).

    Interior segments: the hazard of a switch is `lam` per bin when the segment is old and `lam * young` while its age is below A, so a segment shorter than A + 1 bins pays `young` once (at the switch that ends it).
    First segment (only when Wt is given): its own hazard `lam * term_factor` while it is shorter than Wt bins and `lam` from Wt bins on, independent of `young`; the interior limit can be off (young = 1) and the
    terminal limit on. Last segment: down-weighted by `term_factor` when it is shorter than Wt bins (the age dimension is at least Wt long, so this also works when A is 1). Wt = None: both ends are free.
    A = 1 with young = 1 and Wt = None is the plain two-state chain with switch probability lam."""
    n = len(ll0)
    if n == 0:
        return np.zeros(0)
    W = int(Wt) if Wt else 0
    cap = max(A, W, 1)                                                      # ages 0..cap, cap absorbing
    h = np.full(cap + 1, min(lam, 0.5))
    h[:A] = min(lam * young, 0.5)
    hF = np.full(W, min(lam * term_factor, 0.5)) if W else np.zeros(0)
    if W:
        hF[W - 1] = min(lam, 0.5)
    fend = np.ones(cap + 1)
    if W:
        fend[:max(W - 1, 0)] = term_factor                                  # last segment of 1 .. Wt - 1 bins
    m = np.maximum(ll0, ll1)
    e = np.stack([np.exp(ll0 - m), np.exp(ll1 - m)], 1)
    am, af = np.zeros((n, 2, cap + 1)), np.zeros((n, 2, W))
    a_m, a_f = np.zeros((2, cap + 1)), np.zeros((2, W))
    if W:
        a_f[:, 0] = 0.5
    else:
        a_m[:, cap] = 0.5
    a_m, a_f = a_m * e[0][:, None], a_f * e[0][:, None]
    s0 = a_m.sum() + a_f.sum()
    a_m, a_f = a_m / s0, a_f / s0
    am[0], af[0] = a_m, a_f
    for i in range(1, n):
        stay_m = a_m * (1 - h)[None, :]
        sw = (a_m * h[None, :]).sum(1)
        new_m = np.zeros_like(a_m)
        new_m[:, 1:cap] = stay_m[:, :cap - 1]
        new_m[:, cap] = stay_m[:, cap - 1] + stay_m[:, cap]
        new_f = np.zeros_like(a_f)
        if W:
            stay_f = a_f * (1 - hF)[None, :]
            sw = sw + (a_f * hF[None, :]).sum(1)
            new_f[:, 1:W] = stay_f[:, :W - 1]
            new_m[:, cap] += stay_f[:, W - 1]                              # a first segment of Wt bins without a switch becomes an ordinary old segment
        new_m[0, 0] = sw[1]
        new_m[1, 0] = sw[0]
        new_m *= e[i][:, None]
        new_f *= e[i][:, None]
        s0 = new_m.sum() + new_f.sum()
        a_m, a_f = new_m / s0, new_f / s0
        am[i], af[i] = a_m, a_f
    b_m, b_f = np.ones((2, cap + 1)) * fend[None, :], np.ones((2, W))
    post = np.empty(n)
    qm, qf = am[n - 1] * b_m, af[n - 1] * b_f
    post[n - 1] = (qm[1].sum() + qf[1].sum()) / (qm.sum() + qf.sum())
    for i in range(n - 2, -1, -1):
        eb_m, eb_f = e[i + 1][:, None] * b_m, e[i + 1][:, None] * b_f
        nxt = np.empty_like(b_m)
        nxt[:, :cap] = eb_m[:, 1:cap + 1]
        nxt[:, cap] = eb_m[:, cap]
        swv = np.stack([eb_m[1, 0], eb_m[0, 0]])
        nb_m = nxt * (1 - h)[None, :] + swv[:, None] * h[None, :]
        if W:
            nxt_f = np.empty_like(b_f)
            nxt_f[:, :W - 1] = eb_f[:, 1:W]
            nxt_f[:, W - 1] = eb_m[:, cap]
            nb_f = nxt_f * (1 - hF)[None, :] + swv[:, None] * hF[None, :]
        else:
            nb_f = b_f
        s0 = nb_m.sum() + nb_f.sum()
        b_m, b_f = nb_m / s0, nb_f / s0
        gm, gf = am[i] * b_m, af[i] * b_f
        post[i] = (gm[1].sum() + gf[1].sum()) / (gm.sum() + gf.sum())
    return post


def short_interior_segments(post, bin_bp, short_mb, threshold=0.5):
    """Interior segments (touching neither end) of the hard calls that are shorter than `short_mb`: dicts with the start and end (bp), the length, the state and the mean posterior of the state they are in.
    These are the segments a minimum interior segment would remove: counting them in a run WITHOUT the minimum measures how large the problem is."""
    st = (np.asarray(post) > threshold).astype(int)
    starts = np.concatenate([[0], np.where(np.diff(st) != 0)[0] + 1, [len(st)]])
    out = []
    for k in range(1, len(starts) - 2):
        a, b = int(starts[k]), int(starts[k + 1])
        if (b - a) * bin_bp < short_mb * 1e6:
            seg = np.asarray(post)[a:b]
            out.append(dict(start=a * bin_bp, end=b * bin_bp, length=(b - a) * bin_bp, state=int(st[a]), mean_posterior=float(seg.mean() if st[a] == 1 else 1 - seg.mean())))
    return out


def genotype_fish(chroms, hp):
    """chroms = {chrom: (pos, n_mc, n_yh, yh_is_ref, length)} of one fish. eps and q are fitted by EM over all its chromosomes together, then each chromosome gets its posteriors.
    Returns ({chrom: posterior per bin}, {chrom: counts (4, n_bins)}, eps, q)."""
    bp = hp['bin_bp']
    cs = {ch: bin_counts(p, m, y, r, L, bp, hp['site_cap'], hp['bin_cap']) for ch, (p, m, y, r, L) in chroms.items()}
    lam = hp['rho'] * bp
    A, Wt = segment_limits(hp)
    kw = dict(young=hp.get('interior_hazard_factor', 1e-30) if hp['min_interior_segment_mb'] > 0 else 1.0, Wt=Wt, term_factor=hp.get('terminal_factor', 1e-30))
    eps, q = 0.01, 0.55
    elo, ehi = hp['eps_bounds']
    qlo, qhi = hp['q_bounds']
    for _ in range(hp['em_iterations']):
        ne = de = nq = dq = 0.0
        for c in cs.values():
            if c.sum() == 0:
                continue
            post = forward_backward(*bin_loglik(c, eps, q), lam, A, **kw)
            tot = c.sum(0)
            ne += ((1 - post) * (c[0] + c[2])).sum()
            de += ((1 - post) * tot).sum()
            nq += (post * (c[0] + c[3])).sum()                      # reference-allele reads among the heterozygous: YH reads on `a` sites, MC reads on `b` sites
            dq += (post * tot).sum()
        eps = float(np.clip(ne / max(de, 1e-9), elo, ehi))
        q = float(np.clip(nq / max(dq, 1e-9), qlo, qhi))
    posts = {ch: forward_backward(*bin_loglik(c, eps, q), lam, A, **kw) for ch, c in cs.items()}
    return posts, cs, eps, q


def hard_calls(post, threshold=0.5):
    return (np.asarray(post) > threshold).astype(int)


def breakpoints(post, length, hp):
    """One row per switch of the hard calls: bin positions (bp, the START of the bin on each side), the interval between the last bin confidently in the old state and the first confidently in the new one,
    the posteriors on both sides, the distance to the nearest chromosome end and a flag for terminal calls shorter than `terminal_low_confidence_mb`."""
    bp = hp['bin_bp']
    lo_t, hi_t = hp['uncertain_band']
    state = hard_calls(post, hp['call_threshold'])
    rows = []
    for k in np.where(np.diff(state) != 0)[0]:
        old = state[k]
        left = k
        for j in range(k, -1, -1):
            if (old == 1 and post[j] >= hi_t) or (old == 0 and post[j] <= lo_t):
                left = j
                break
        right = k + 1
        for j in range(k + 1, len(post)):
            if (old == 0 and post[j] >= hi_t) or (old == 1 and post[j] <= lo_t):
                right = j
                break
        x = (k + 1) * bp
        d = min(x, length - x)
        rows.append(dict(bin_lo=int(k * bp), bin_hi=int((k + 1) * bp), interval_lo=int(left * bp), interval_hi=int((right + 1) * bp), p_left=float(post[k]), p_right=float(post[k + 1]),
                         to_state=int(state[k + 1]), distance_to_end=int(d), terminal_low_confidence=bool(d < hp['terminal_low_confidence_mb'] * 1e6)))
    return rows


def end_support(counts, post, hp, end_mb=1.0):
    """What the ends of a chromosome rest on: reads in the first and last `end_mb`, and the share of the first and last bins (that long) with a confident posterior."""
    nb = counts.shape[1]
    k = max(1, int(round(end_mb * 1e6 / hp['bin_bp'])))
    lo, hi = hp['uncertain_band']
    tot = counts.sum(0)
    conf = (post <= lo) | (post >= hi)
    return dict(left_reads=float(tot[:k].sum()), right_reads=float(tot[-k:].sum()), left_confident=float(conf[:k].mean()), right_confident=float(conf[-k:].mean()), n_bins=nb)
