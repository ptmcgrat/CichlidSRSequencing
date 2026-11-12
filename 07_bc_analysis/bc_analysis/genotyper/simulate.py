"""Simulated backcrosses with known crossovers, to test the genotyper and to compare it with the original script. `python -m bc_analysis genotype simulate` reproduces the table of the README.

What is simulated. A chromosome of one BC fish: 0 to 3 crossovers in the F1 sire's meiosis (a minimum spacing stands for interference; optionally 40% of them near a chromosome end), informative sites at a given density,
reads per site Poisson(coverage), an error rate for YH-looking reads in the MC state, a reference bias in the heterozygous state, a polarity per site (the YH allele is the reference allele or not), markers depleted near the ends,
and ARTIFACTS: repeat regions of 50 to 400 kb in which MC-state sites show extra YH-looking reads (mismapping), also at the very end (subtelomeric repeats). The artifacts are what a high switch penalty was invented against.
"""
import numpy as np
import pandas as pd

from . import hmm as H

BIN = 100_000
SCENARIOS = {
    'S1': dict(label='old-genome-like: YH allele always ALT, even marker density, no end artifacts', mixed=False, end_depleted=False, end_artifact=False, end_bias=True),
    'S2': dict(label='new-genome-like: mixed polarity, markers sparse near the ends, subtelomeric artifacts', mixed=True, end_depleted=True, end_artifact=True, end_bias=True),
}


# ------------------------------------------------------------------ truth
def draw_crossovers(rng, length, end_bias=False, dmin=8e6):
    n = rng.choice([0, 1, 2, 3], p=[0.50, 0.40, 0.09, 0.01])
    xs = []
    for _ in range(n):
        for _try in range(50):
            if end_bias and rng.random() < 0.4:
                d = min(rng.exponential(1.5e6), length / 2)
                x = d if rng.random() < 0.5 else length - d
            else:
                x = rng.uniform(0, length)
            if all(abs(x - y) >= dmin for y in xs):
                xs.append(x)
                break
    return sorted(xs)


def truth_bins(length, xs, s0, bin_bp=BIN):
    nb = int(np.ceil(length / bin_bp))
    centers = (np.arange(nb) + 0.5) * bin_bp
    st = np.full(nb, s0, int)
    for x in xs:
        st[centers > x] ^= 1
    return st


# ------------------------------------------------------------------ reads
def simulate_reads(rng, length, xs, s0, cov, density=1.0, end_depleted=False, mixed=False, eps=0.01, q=0.6, end_artifact=False, n_repeats=2):
    """(pos, n_mc, n_yh, yh_is_ref) of the sites with at least one read."""
    n_sites = int(length / 1000 * density)
    pos = np.sort(rng.integers(1000, int(length), n_sites))
    if end_depleted:
        keep = np.ones(len(pos), bool)
        near = (pos < 1.5e6) | (pos > length - 1.5e6)
        keep[near] = rng.random(near.sum()) < 0.25
        pos = pos[keep]
    crossed = (pos[:, None] > np.array(xs)[None, :]).sum(1) % 2 if len(xs) else np.zeros(len(pos), int)
    state = np.where(crossed == 0, s0, 1 - s0)
    yh_is_ref = (rng.random(len(pos)) < 0.5) if mixed else np.zeros(len(pos), bool)
    p_yh = np.where(state == 0, eps, np.where(yh_is_ref, q, 1 - q))
    for k in range(n_repeats + (1 if end_artifact else 0)):
        if rng.random() < 0.7 or (end_artifact and k == n_repeats):
            ln = rng.uniform(5e4, 4e5)
            a = rng.uniform(0, length - ln) if k < n_repeats else (0 if rng.random() < 0.5 else length - ln)
            inside = (pos >= a) & (pos < a + ln) & (state == 0) & (rng.random(len(pos)) < 0.5)
            p_yh = np.where(inside, 0.4, p_yh)
    n = rng.poisson(cov, len(pos))
    keep = n > 0
    n_yh = rng.binomial(n, p_yh)
    return pos[keep], (n - n_yh)[keep], n_yh[keep], yh_is_ref[keep]


# ------------------------------------------------------------------ the original script's model, as written
def original_like(pos, n_mc, n_yh, yh_is_ref, length, trans=1e-22, oriented=False, bin_bp=BIN):
    """GenotypeBCs._hmmVector: Viterbi, two states, switch probability 1e-22 (a penalty of about 50 nats), emissions [.97 .03] and [.6 .4] for (REF read, ALT read), then linear interpolation to a 100 kb grid
    from 0 to the last read. Reads of one site are exchangeable, so the per-read chain equals this per-site one. oriented=False reads REF as 0 and ALT as 1 whatever the polarity (the script on a new reference);
    oriented=True names the alleles (MC = 0, YH = 1), which is what the script did on the MC reference."""
    if oriented:
        n_ref, n_alt = np.asarray(n_mc), np.asarray(n_yh)
    else:
        n_ref = np.where(yh_is_ref, n_yh, n_mc)
        n_alt = np.where(yh_is_ref, n_mc, n_yh)
    l0 = n_ref * np.log(.97) + n_alt * np.log(.03)
    l1 = n_ref * np.log(.6) + n_alt * np.log(.4)
    pen = -np.log(trans)
    v0, v1 = l0[0], l1[0]
    b0 = np.zeros(len(pos), bool)
    b1 = np.zeros(len(pos), bool)
    for i in range(1, len(pos)):
        s0 = v0 >= v1 - pen
        s1 = v1 >= v0 - pen
        b0[i], b1[i] = not s0, not s1
        v0, v1 = (v0 if s0 else v1 - pen) + l0[i], (v1 if s1 else v0 - pen) + l1[i]
    st = 0 if v0 >= v1 else 1
    path = np.empty(len(pos), np.uint8)
    for i in range(len(pos) - 1, -1, -1):
        path[i] = st
        if i and (b0[i] if st == 0 else b1[i]):
            st = 1 - st
    grid = np.arange(0, pos.max(), bin_bp)
    return (np.interp(grid, pos, path) >= 0.5).astype(int)


def v2_like(fish, rho=1e-8, bin_bp=BIN):
    """The repo's earlier v2: site-level forward-backward with switch probabilities that depend on the distance between sites."""
    from ..bc_ancestry_hmm import forward_backward, fit_sample, with_grid
    per = []
    for (pos, n_mc, n_yh, yh_is_ref, length) in fish:
        per.append(with_grid(pos, n_mc, n_yh, length, bin_bp, offset=0, flag=yh_is_ref))
    eps, q = fit_sample([p[0] for p in per], [p[1] for p in per], [p[2] for p in per], rho=rho, cap=2, flags=[p[4] for p in per])
    return [(forward_backward(ap, ar, aa, eps, q, rho, 2, fl)[isg] > .5).astype(int) for ap, ar, aa, isg, fl in per]


# ------------------------------------------------------------------ scoring
def score(called, truth_state, length, tol=1.0e6, bin_bp=BIN):
    """The truth is the set of switches of the true bin states (what a 100 kb grid can represent). Returns (detections [(distance of the true switch to the nearest end, found within tol)], spurious switches,
    flip-flops (consecutive called switches less than 2 Mb apart), spurious switches within 2 Mb of an end, bin accuracy)."""
    nb = len(truth_state)
    called = np.asarray(called)[:nb]
    if len(called) < nb:
        called = np.concatenate([called, np.full(nb - len(called), called[-1])])
    tsw = (np.where(np.diff(truth_state) != 0)[0] + 1) * bin_bp
    sp = (np.where(np.diff(called) != 0)[0] + 1) * bin_bp
    matched, det = set(), []
    for x in tsw:
        d = np.abs(sp - x) if len(sp) else np.array([np.inf])
        j = int(np.argmin(d)) if len(sp) else -1
        ok = len(sp) > 0 and d[j] <= tol and j not in matched
        if ok:
            matched.add(j)
        det.append((float(min(x, length - x)), bool(ok)))
    spurious = [i for i in range(len(sp)) if i not in matched]
    flip = sum(1 for i in range(len(sp) - 1) if sp[i + 1] - sp[i] < 2e6)
    term = sum(1 for i in spurious if min(sp[i], length - sp[i]) < 2e6)
    return det, len(spurious), flip, term, float((called == truth_state).mean())


def methods(hp):
    """name -> function(fish list) -> list of called-state arrays. The HMM variants differ only in the two segment limits."""
    def hmm_variant(**over):
        h = dict(hp, **over)

        def run(fish):
            out = []
            posts, _, _, _ = H.genotype_fish({k: f for k, f in enumerate(fish)}, h)
            for k in range(len(fish)):
                out.append(H.hard_calls(posts[k], h['call_threshold']))
            return out
        return run
    return {
        'original (as written)': lambda fish: [original_like(*f) for f in fish],
        'original (MC/YH alleles named)': lambda fish: [original_like(*f, oriented=True) for f in fish],
        'v2': v2_like,
        'bins, no segment limits': hmm_variant(min_interior_segment_mb=0, min_terminal_segment_mb=0),
        'bins, terminal limit only (the first run)': hmm_variant(min_interior_segment_mb=0, min_terminal_segment_mb=0.2),
        'bins, interior limit': hmm_variant(min_interior_segment_mb=2.0, min_terminal_segment_mb=0),
        'bins, interior + terminal limits': hmm_variant(min_interior_segment_mb=2.0, min_terminal_segment_mb=0.2),
    }


def run_comparison(hp, coverages=(0.3, 0.5, 1.0), scenario='S2', n_fish=20, n_chr=6, length=40e6, seed=1, which=None):
    """One row per (coverage, method, chromosome): detections, spurious switches, flip-flops, terminal spurious switches, accuracy."""
    sc = SCENARIOS[scenario]
    rows = []
    for cov in coverages:
        rng = np.random.default_rng(seed + int(cov * 100))
        meths = {k: v for k, v in methods(hp).items() if which is None or k in which}
        for _ in range(n_fish):
            fish, truths = [], []
            for _c in range(n_chr):
                xs = draw_crossovers(rng, length, sc['end_bias'])
                s0 = int(rng.integers(0, 2))
                p, m, y, r = simulate_reads(rng, length, xs, s0, cov, 1.0, sc['end_depleted'], sc['mixed'], end_artifact=sc['end_artifact'] and rng.random() < 0.4)
                fish.append((p, m, y, r, length))
                truths.append(truth_bins(length, xs, s0))
            for name, fn in meths.items():
                calls = fn(fish)
                for k in range(n_chr):
                    det, nsp, flip, term, acc = score(calls[k], truths[k], length)
                    rows.append(dict(scenario=scenario, cov=cov, method=name, det=det, spurious=nsp, flip=flip, term_spurious=term, acc=acc))
    return rows


def summarize(rows):
    """(recall table by distance to the nearest end, false-call table)."""
    recs = []
    for r in rows:
        for d, ok in r['det']:
            b = '<0.5 Mb' if d < 5e5 else '0.5-1 Mb' if d < 1e6 else '1-2 Mb' if d < 2e6 else '>=2 Mb'
            recs.append(dict(cov=r['cov'], method=r['method'], bucket=b, ok=ok))
    order = ['<0.5 Mb', '0.5-1 Mb', '1-2 Mb', '>=2 Mb']
    x = pd.DataFrame(recs)
    rec = x.groupby(['cov', 'method', 'bucket']).ok.agg(['mean', 'size']).reset_index()
    t = rec.pivot_table(index=['cov', 'method'], columns='bucket', values='mean').reindex(columns=order).round(2)
    n = rec.pivot_table(index=['cov', 'method'], columns='bucket', values='size').reindex(columns=order).fillna(0).astype(int)
    recall = (t.astype(str) + ' (' + n.astype(str) + ')')
    d = pd.DataFrame([dict(cov=r['cov'], method=r['method'], spurious=r['spurious'], flip=r['flip'], term=r['term_spurious'], acc=r['acc']) for r in rows])
    false = d.groupby(['cov', 'method']).agg(spurious_switches_per_100_chr=('spurious', lambda s: 100 * s.mean()), flip_flops_per_100_chr=('flip', lambda s: 100 * s.mean()),
                                              terminal_spurious_per_100_chr=('term', lambda s: 100 * s.mean()), bin_accuracy=('acc', 'mean')).round(3)
    return recall, false
