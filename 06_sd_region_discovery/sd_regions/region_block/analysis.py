"""Statistics of the region block analysis. Pure numpy/scipy: genotype matrices in, numbers out, so that every function can be tested on data with a known answer.

Conventions. A genotype matrix `dos` is (sites x samples) int8, the alternative-allele dosage 0/1/2, -1 = missing (the masking by depth has already been done). `male` and `female` are boolean vectors
over the samples. A "window" is a stretch of `size` bases of one chromosome; its heterozygosity in a sample is the heterozygous calls divided by the calls in that window, over the SNPs the cohort is variable at.

What each part is for (the questions it answers):
  * WindowAccumulator + roh_segments + region_summary      Are the females (or the males) HOMOZYGOUS over the block, which individuals, and how long are their runs of homozygosity?
  * GrmAccumulator + kinship_sex_test + related_clusters   How related are the fish, how inbred, and does relatedness line up with sex (a family confounded with sex mimics a sex locus)?
  * pca_scores + auc_sex                                   Do the SNPs of the block split the fish into groups that follow sex, and do those groups differ from the genome-wide structure?
  * ld_matrix + ld_by_distance                             Is there a block of linkage disequilibrium (an inversion or a non-recombining haplotype), in which sex?
  * diag_* + diag_subsample                                How many sex-diagnostic SNPs are there under stricter rules and other depth cut-offs, and how does the count depend on the number of fish?
"""
import numpy as np
from scipy import stats as sps


# ====================================================================== windows
class WindowAccumulator:
    """Per-sample counts of called and heterozygous genotypes per (chromosome, window), filled chunk by chunk from a sorted VCF."""

    def __init__(self, n_samples, size):
        self.n, self.size = n_samples, int(size)
        self.called, self.het, self.order = {}, {}, []

    def add(self, chrom, pos, dos):
        if len(pos) == 0:
            return
        wid = np.asarray(pos) // self.size
        chrom = np.asarray(chrom)
        change = np.ones(len(pos), bool)
        change[1:] = (chrom[1:] != chrom[:-1]) | (wid[1:] != wid[:-1])
        starts = np.flatnonzero(change)
        C = np.add.reduceat((dos >= 0).astype(np.int32), starts, axis=0)
        H = np.add.reduceat((dos == 1).astype(np.int32), starts, axis=0)
        for s, c, h in zip(starts, C, H):
            k = (str(chrom[s]), int(wid[s]))
            if k in self.called:
                self.called[k] += c
                self.het[k] += h
            else:
                self.called[k], self.het[k] = c.copy(), h.copy()
                if k[0] not in self.order:
                    self.order.append(k[0])

    def arrays(self):
        """(keys, called, het): keys = [(chrom, window index)] in chromosome order then position; called and het are (windows x samples)."""
        keys = sorted(self.called, key=lambda k: (self.order.index(k[0]), k[1]))
        if not keys:
            return [], np.zeros((0, self.n), int), np.zeros((0, self.n), int)
        return keys, np.vstack([self.called[k] for k in keys]), np.vstack([self.het[k] for k in keys])


def window_het(called, het, min_sites):
    """Heterozygosity per window and sample, NaN where fewer than min_sites genotypes were called."""
    with np.errstate(invalid='ignore', divide='ignore'):
        h = het / called
    h[called < min_sites] = np.nan
    return h


def roh_segments(keys, called, het, size, het_max=0.03, min_sites=30, min_windows=4):
    """Runs of homozygosity from window heterozygosity: a window is homozygous in a sample when at least `min_sites` genotypes were called and at most `het_max` of them are heterozygous; a run is
    min_windows or more ADJACENT homozygous windows of one chromosome. Returns (segments, callable), segments = [dict(sample, chrom, start, end, n_windows)], callable = (windows x samples) bool.
    This is deliberately simple and transparent (a threshold on window heterozygosity, no model): the resolution is one window, and a window that holds few SNPs cannot be judged."""
    callable_ = called >= min_sites
    with np.errstate(invalid='ignore', divide='ignore'):
        hom = callable_ & (het <= het_max * np.maximum(called, 1))
    chroms = [k[0] for k in keys]
    widx = np.array([k[1] for k in keys]) if keys else np.zeros(0, int)
    segs = []
    for s in range(called.shape[1]):
        i = 0
        while i < len(keys):
            if hom[i, s]:
                j = i
                while j + 1 < len(keys) and hom[j + 1, s] and chroms[j + 1] == chroms[i] and widx[j + 1] == widx[j] + 1:
                    j += 1
                if j - i + 1 >= min_windows:
                    segs.append(dict(sample=s, chrom=chroms[i], start=int(widx[i] * size), end=int((widx[j] + 1) * size), n_windows=j - i + 1))
                i = j + 1
            else:
                i += 1
    return segs, callable_


def froh(segs, callable_, size, n_samples):
    """Fraction of the callable genome in a run of homozygosity, per sample."""
    total = callable_.sum(axis=0) * size
    inroh = np.zeros(n_samples)
    for g in segs:
        inroh[g['sample']] += g['end'] - g['start']
    with np.errstate(invalid='ignore', divide='ignore'):
        return inroh / total


def region_windows(keys, chrom, start, end, size):
    """Indices of the windows that overlap chrom:start-end."""
    return np.array([i for i, (c, w) in enumerate(keys) if c == chrom and w * size < end and (w + 1) * size > start], dtype=int)


def region_summary(keys, called, het, segs, callable_, chrom, start, end, size, het_max=0.03):
    """Per sample, over the windows of the region: heterozygosity, the share of callable windows that are homozygous, and the share inside a run of homozygosity (>= min_windows adjacent)."""
    idx = region_windows(keys, chrom, start, end, size)
    n = called.shape[1]
    if len(idx) == 0:
        nan = np.full(n, np.nan)
        return dict(het=nan, hom_share=nan, roh_share=nan, n_windows=0)
    C, H = called[idx].sum(axis=0), het[idx].sum(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        reg_het = H / C
    ok = callable_[idx]
    with np.errstate(invalid='ignore', divide='ignore'):
        hom = ok & (het[idx] <= het_max * np.maximum(called[idx], 1))
        hom_share = hom.sum(axis=0) / ok.sum(axis=0)
    inroh = np.zeros((len(idx), n), bool)
    wmap = {keys[i]: r for r, i in enumerate(idx)}
    for g in segs:
        if g['chrom'] != chrom:
            continue
        for w in range(g['start'] // size, g['end'] // size):
            r = wmap.get((chrom, w))
            if r is not None:
                inroh[r, g['sample']] = True
    with np.errstate(invalid='ignore', divide='ignore'):
        roh_share = (inroh & ok).sum(axis=0) / ok.sum(axis=0)
    return dict(het=reg_het, hom_share=hom_share, roh_share=roh_share, n_windows=len(idx))


def chrom_sex_het(keys, called, het, male, female):
    """{chromosome: (male het, female het, difference)} from the window counts."""
    out = {}
    for c in dict.fromkeys(k[0] for k in keys):
        idx = [i for i, k in enumerate(keys) if k[0] == c]
        m = het[idx][:, male].sum() / max(1, called[idx][:, male].sum())
        f = het[idx][:, female].sum() / max(1, called[idx][:, female].sum())
        out[c] = (float(m), float(f), float(m - f))
    return out


# ====================================================================== relatedness and inbreeding
class GrmAccumulator:
    """Genomic relationship matrix (VanRaden/Yang) and per-sample inbreeding from heterozygosity, accumulated chunk by chunk.

    A site counts when it is called in at least 2 samples and its minor allele frequency (in the cohort) is at least min_maf. G_ij = sum_s z_is z_js / (sites called in both), z = (x - 2p) / sqrt(2p(1-p)).
    G_ii is about 1 + F; G_ij about the relationship coefficient (1 identical or duplicate, 0.5 full sibs or parent-offspring, 0.25 half sibs, 0 unrelated). Allele frequencies come from the same
    few fish, so the values are biased by about -1/(n-1) off the diagonal: use them to compare pairs and groups, not as absolute pedigree coefficients."""

    def __init__(self, n_samples, min_maf=0.05):
        self.n, self.min_maf = n_samples, min_maf
        self.ZZ, self.W = np.zeros((n_samples, n_samples)), np.zeros((n_samples, n_samples))
        self.het_obs, self.het_exp, self.n_sites = np.zeros(n_samples), np.zeros(n_samples), 0

    def add(self, dos):
        called = dos >= 0
        nc = called.sum(axis=1)
        p = np.where(called, dos, 0).sum(axis=1) / (2.0 * np.maximum(nc, 1))
        maf = np.minimum(p, 1 - p)
        keep = (nc >= 2) & (maf >= self.min_maf)
        if not keep.any():
            return
        x, c, p, nck = dos[keep].astype(float), called[keep], p[keep], nc[keep]
        Z = np.where(c, (x - 2 * p[:, None]) / np.sqrt(2 * p * (1 - p))[:, None], 0.0)
        self.ZZ += Z.T @ Z
        Cf = c.astype(float)
        self.W += Cf.T @ Cf
        he = 2 * p * (1 - p) * (2.0 * nck) / np.maximum(2.0 * nck - 1, 1)        # unbiased expected heterozygosity of the site
        self.het_obs += ((x == 1) & c).sum(axis=0)
        self.het_exp += (he[:, None] * Cf).sum(axis=0)
        self.n_sites += int(keep.sum())

    def grm(self):
        with np.errstate(invalid='ignore', divide='ignore'):
            return self.ZZ / self.W

    def f_het(self):
        """F = 1 - observed / expected heterozygosity, per sample."""
        with np.errstate(invalid='ignore', divide='ignore'):
            return 1 - self.het_obs / self.het_exp


def kinship_sex_test(G, sex, n_perm=10000, seed=1):
    """Does relatedness line up with sex? Statistic = mean G over pairs of the SAME sex minus the mean over pairs of DIFFERENT sex; the p-value permutes the sex labels (one-sided: same-sex pairs
    more related). `sex` is an array of 'M'/'F'. Returns dict(within_f, within_m, between, stat, p)."""
    sex = np.asarray(sex)
    iu = np.triu_indices(len(sex), 1)
    g = G[iu]
    ok = np.isfinite(g)
    si, sj = sex[iu[0]], sex[iu[1]]

    def stat(a, b):
        same = (a == b) & ok
        diff = (a != b) & ok
        return g[same].mean() - g[diff].mean()

    obs = stat(si, sj)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(sex)
        count += stat(perm[iu[0]], perm[iu[1]]) >= obs
    wf = g[(si == 'F') & (sj == 'F') & ok]
    wm = g[(si == 'M') & (sj == 'M') & ok]
    bt = g[(si != sj) & ok]
    return dict(within_f=float(wf.mean()) if len(wf) else np.nan, within_m=float(wm.mean()) if len(wm) else np.nan, between=float(bt.mean()), stat=float(obs), p=float((count + 1) / (n_perm + 1)))


def related_clusters(G, threshold=0.2):
    """Connected components of the graph 'G_ij >= threshold' (union-find): [[sample indices], ...], singletons included, largest first."""
    n = G.shape[0]
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for i in range(n):
        for j in range(i + 1, n):
            if np.isfinite(G[i, j]) and G[i, j] >= threshold:
                parent[find(i)] = find(j)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return sorted(groups.values(), key=lambda g: (-len(g), g[0]))


# ====================================================================== structure: PCA
def pca_scores(dos, min_maf=0.05, n_comp=3):
    """Principal components of the samples from the SNPs of `dos` (sites x samples). Each site is standardised by its allele frequency, missing genotypes are set to the mean (0).
    Returns (scores (n x k), variance explained (k,), sites used)."""
    called = dos >= 0
    nc = called.sum(axis=1)
    p = np.where(called, dos, 0).sum(axis=1) / (2.0 * np.maximum(nc, 1))
    keep = (nc >= 2) & (np.minimum(p, 1 - p) >= min_maf)
    if keep.sum() < 10:
        raise ValueError(f'only {int(keep.sum())} usable SNPs for the PCA (need at least 10)')
    x, c, pk = dos[keep].astype(float), called[keep], p[keep]
    Z = np.where(c, (x - 2 * pk[:, None]) / np.sqrt(2 * pk * (1 - pk))[:, None], 0.0)
    U, S, _ = np.linalg.svd(Z.T, full_matrices=False)
    k = min(n_comp, len(S))
    return U[:, :k] * S[:k], (S ** 2 / (S ** 2).sum())[:k], int(keep.sum())


def pca_from_grm(G, n_comp=3):
    """The same kind of scores from a relationship matrix (eigenvectors scaled by the root of their eigenvalue); missing pairs are set to 0."""
    A = np.nan_to_num(G, nan=0.0)
    w, v = np.linalg.eigh((A + A.T) / 2)
    order = np.argsort(w)[::-1][:n_comp]
    return v[:, order] * np.sqrt(np.maximum(w[order], 0)), w[order] / np.maximum(w[w > 0].sum(), 1e-12)


def auc_sex(score, is_male, n_perm=5000, seed=1):
    """How well does one score separate the sexes? AUC (Mann-Whitney) folded to >= 0.5, with a permutation p-value of the folded AUC. 1.0 = perfect separation, 0.5 = none."""
    score, is_male = np.asarray(score, float), np.asarray(is_male, bool)
    n1, n0 = int(is_male.sum()), int((~is_male).sum())
    if n1 == 0 or n0 == 0:
        return np.nan, np.nan
    ranks = sps.rankdata(score)

    def folded(lbl):
        u = ranks[lbl].sum() - n1 * (n1 + 1) / 2
        a = u / (n1 * n0)
        return max(a, 1 - a)
    obs = folded(is_male)
    rng = np.random.default_rng(seed)
    cnt = sum(folded(rng.permutation(is_male)) >= obs - 1e-12 for _ in range(n_perm))
    return float(obs), float((cnt + 1) / (n_perm + 1))


def best_pc_auc(scores, is_male, n_perm=5000, seed=1):
    """(AUC, p, pc number from 1) of the best of the given components at separating the sexes (the p-value is NOT corrected for choosing among components)."""
    best = (0.0, 1.0, 1)
    for k in range(scores.shape[1]):
        a, p = auc_sex(scores[:, k], is_male, n_perm, seed)
        if a > best[0]:
            best = (a, p, k + 1)
    return best


# ====================================================================== linkage disequilibrium
def ld_matrix(dos):
    """r^2 between all pairs of the sites of `dos` (sites x samples). Missing genotypes are replaced by the site's mean; a site that does not vary gives NaN."""
    x = dos.astype(float)
    x[dos < 0] = np.nan
    mu = np.nanmean(x, axis=1)
    x = np.where(np.isnan(x), mu[:, None], x)
    x -= mu[:, None]
    sd = np.sqrt((x ** 2).sum(axis=1))
    with np.errstate(invalid='ignore', divide='ignore'):
        r = (x @ x.T) / np.outer(sd, sd)
    r[(sd == 0), :] = np.nan
    r[:, (sd == 0)] = np.nan
    return r ** 2


def thin_sites(dos, pos, min_maf=0.15, max_sites=300):
    """Up to max_sites sites evenly spaced in the list, with minor allele frequency >= min_maf that VARY among the fish given (a site at which every fish is heterozygous has frequency 0.5 and no variation,
    and gives no linkage disequilibrium to measure): (dos subset, pos subset)."""
    called = dos >= 0
    p = np.where(called, dos, 0).sum(axis=1) / (2.0 * np.maximum(called.sum(axis=1), 1))
    varies = np.where(called, dos, -1).max(axis=1) > np.where(called, dos, 9).min(axis=1)
    ok = np.flatnonzero((np.minimum(p, 1 - p) >= min_maf) & varies)
    if len(ok) > max_sites:
        ok = ok[np.linspace(0, len(ok) - 1, max_sites).round().astype(int)]
    return dos[ok], pos[ok]


def ld_by_distance(r2, pos, edges):
    """Mean r^2 of the pairs whose distance falls in each bin [edges[i], edges[i+1]) -> array (NaN for an empty bin)."""
    d = np.abs(pos[:, None] - pos[None, :])
    iu = np.triu_indices(len(pos), 1)
    dd, rr = d[iu], r2[iu]
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (dd >= lo) & (dd < hi) & np.isfinite(rr)
        out.append(rr[sel].mean() if sel.any() else np.nan)
    return np.array(out)


# ====================================================================== sex-diagnostic SNPs
def _wc():
    from sexloci.popgen import wc_fst_components            # the vendored definition: the same function sexloci uses
    return wc_fst_components


def diag_sexloci(dos, male, female, tol=0.1, min_called=3):
    """Boolean per site: sexloci's rule for a male-diagnostic ('XY') SNP: at least (1 - tol) of the called males heterozygous, at most tol of the called females heterozygous, and the females' minor
    allele frequency at most tol / 2. With 12 fish of a sex this allows one exception per sex, with 9 or 7 none: the count is NOT comparable between runs of different size."""
    _, _, _, n1, n2, p1, p2, h1, h2 = _wc()(dos[:, male], dos[:, female])
    ok = (n1 >= min_called) & (n2 >= min_called)
    maf2 = np.minimum(p2, 1 - p2)
    return ok & (h1 >= 1 - tol) & (h2 <= tol) & (maf2 <= tol / 2)


def diag_strict(dos, male, female, min_called=8):
    """Boolean per site: EVERY called male heterozygous and EVERY called female homozygous for the same allele (no exceptions), with at least min_called of each sex called."""
    dm, df = dos[:, male], dos[:, female]
    nm, nf = (dm >= 0).sum(1), (df >= 0).sum(1)
    all_het = ((dm == 1) | (dm < 0)).all(axis=1)
    f0 = ((df == 0) | (df < 0)).all(axis=1)
    f2 = ((df == 2) | (df < 0)).all(axis=1)
    return (nm >= min_called) & (nf >= min_called) & all_het & (f0 | f2)


def diag_by_depth(dos_raw, dp, male, female, thresholds=(3, 5, 8), tol=0.1):
    """{min depth: (sexloci-rule count, strict count)}: genotypes with depth below the threshold are masked before the rule is applied."""
    out = {}
    for t in thresholds:
        d = dos_raw.copy()
        d[np.nan_to_num(dp, nan=np.inf) < t] = -1
        out[t] = (int(diag_sexloci(d, male, female, tol).sum()), int(diag_strict(d, male, female).sum()))
    return out


def diag_subsample(dos, male, female, n_m, n_f, n_draws=200, tol=0.1, seed=1):
    """The sexloci-rule count in random subsets of n_m males and n_f females: array of n_draws counts. Shows how much the count depends on the number of fish and on which fish are included."""
    rng = np.random.default_rng(seed)
    mi, fi = np.flatnonzero(male), np.flatnonzero(female)
    out = []
    for _ in range(n_draws):
        m = np.zeros(len(male), bool)
        f = np.zeros(len(male), bool)
        m[rng.choice(mi, n_m, replace=False)] = True
        f[rng.choice(fi, n_f, replace=False)] = True
        out.append(int(diag_sexloci(dos, m, f, tol).sum()))
    return np.array(out)


def diag_sample_het(dos, sites_mask):
    """Per sample: the share of its called genotypes at the diagnostic sites that are heterozygous (males should be near 1 and females near 0 under a real XY haplotype)."""
    d = dos[sites_mask]
    called = (d >= 0).sum(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        return (d == 1).sum(axis=0) / called
