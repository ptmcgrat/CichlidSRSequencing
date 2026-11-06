"""Parent-offspring (family-based) sex-linkage tests.
XY: the FATHER is heterozygous and the MOTHER homozygous -> the paternal allele an offspring inherited is known; at a site linked
to the sex-determining locus sons carry one paternal allele and daughters the other. ZW: roles of the parents swap.
Works on hard genotype calls with depth/GQ filters (offspring are 3-13x, so single sites are noisy; windows aggregate sites)."""
from functools import lru_cache
import numpy as np
import pandas as pd
from scipy import stats


@lru_cache(maxsize=None)
def _fisher(a, b, c, d):
    return stats.fisher_exact([[a, b], [c, d]])[1]


def alt_dosage(gt, dp, gq, min_dp, min_gq):
    """gt: (sites, samples, 2) allele array (-1 missing). Returns alt allele count 0/1/2, NaN if missing or low quality."""
    miss = (gt < 0).any(-1)
    dose = np.where(miss, np.nan, gt.sum(-1).astype(float))
    bad = miss | (dp < min_dp) | (gq < min_gq)
    return np.where(bad, np.nan, dose)


def transmitted_alt(parent_dose, other_dose, off_dose):
    """Alt allele transmitted by `parent` (het) to each offspring, given the other parent is hom (dose 0 or 2).
    Returns (sites, offspring) array of 0/1/NaN; NaN where uninformative or inconsistent (error / de novo)."""
    informative = (parent_dose == 1) & ((other_dose == 0) | (other_dose == 2))                 # (sites,)
    contrib_other = (other_dose / 2.0)[:, None]
    t = off_dose - contrib_other
    t = np.where((t == 0) | (t == 1), t, np.nan)
    return np.where(informative[:, None], t, np.nan)


def site_scan(trans, is_male, min_n=12):
    """Per-site association of transmitted allele with offspring sex. Returns DataFrame (one row per informative site)."""
    male = np.asarray(is_male, bool)
    ok = ~np.isnan(trans)
    n = ok.sum(1)
    t1 = np.nan_to_num(trans, nan=0) == 1
    a = (t1 & ok & male).sum(1); b = (t1 & ok & ~male).sum(1)
    c = (~t1 & ok & male).sum(1); d = (~t1 & ok & ~male).sum(1)
    keep = n >= min_n
    idx = np.where(keep)[0]
    p = np.array([_fisher(int(a[i]), int(b[i]), int(c[i]), int(d[i])) for i in idx]) if len(idx) else np.array([])
    conc = np.maximum(a + d, b + c)[idx] / n[idx] if len(idx) else np.array([])
    return pd.DataFrame(dict(site=idx, n=n[idx], concordance=conc, p=p, nlp=-np.log10(np.maximum(p, 1e-300))))


def window_summary(pos, scan, size, perfect=0.95):
    """Aggregate per-site results into windows: informative sites, median concordance, fraction near-perfect, best -log10 p."""
    if len(scan) == 0:
        return pd.DataFrame(columns=['start', 'n_inf', 'med_conc', 'frac_perfect', 'max_nlp'])
    scan = scan.copy(); scan['start'] = (pos[scan.site.values] // size) * size
    g = scan.groupby('start')
    return pd.DataFrame(dict(n_inf=g.size(), med_conc=g.concordance.median(), frac_perfect=g.concordance.apply(lambda x: (x >= perfect).mean()),
                             max_nlp=g.nlp.max())).reset_index()


def family_scan(pos, dose, parents, sexes, size, min_n=12):
    """dose: (sites, samples) alt dosage; parents: dict(father=col, mother=col, offspring=[cols]); sexes: array over offspring (True=male).
    Runs the XY test (father het) and the ZW test (mother het). Returns {'XY': (site_df, window_df), 'ZW': (...)}."""
    f, m, o = parents['father'], parents['mother'], parents['offspring']
    res = {}
    for mode, (het, hom) in {'XY': (f, m), 'ZW': (m, f)}.items():
        tr = transmitted_alt(dose[:, het], dose[:, hom], dose[:, o])
        sc = site_scan(tr, sexes, min_n)
        res[mode] = (sc, window_summary(pos, sc, size))
    return res
