"""Population-genetic summaries by sex: Hudson Fst, heterozygosity, SNP density, and normalised depth ratios in windows."""
import numpy as np
import pandas as pd
from scipy import stats


def group_counts(dose, cols):
    """Alt allele count and number of sampled chromosomes among `cols` (hard calls; NaN = missing)."""
    d = dose[:, cols]
    called = ~np.isnan(d)
    return np.nansum(d, 1), 2 * called.sum(1)


def hudson_terms(alt1, n1, alt2, n2, min_n=4):
    """Per-site numerator/denominator of Hudson's Fst estimator (Bhatia et al. 2013). NaN where a group has < min_n chromosomes."""
    ok = (n1 >= min_n) & (n2 >= min_n)
    p1 = np.where(ok, alt1 / np.maximum(n1, 1), np.nan); p2 = np.where(ok, alt2 / np.maximum(n2, 1), np.nan)
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / np.maximum(n1 - 1, 1) - p2 * (1 - p2) / np.maximum(n2 - 1, 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return num, den


def window_index(pos, size):
    return (pos // size) * size


def windowed_fst(pos, num, den, size):
    ok = ~np.isnan(num) & ~np.isnan(den)
    df = pd.DataFrame(dict(w=window_index(pos[ok], size), num=num[ok], den=den[ok]))
    g = df.groupby('w').agg(num=('num', 'sum'), den=('den', 'sum'), n_sites=('num', 'size'))
    g['fst'] = g.num / g.den.replace(0, np.nan)
    return g[['fst', 'n_sites']].reset_index().rename(columns={'w': 'start'})


def het_by_group(dose, cols, min_called=4):
    """Observed heterozygosity per site among called individuals; NaN if fewer than min_called called."""
    d = dose[:, cols]; called = ~np.isnan(d)
    n = called.sum(1)
    h = np.where(n >= min_called, (d == 1).sum(1) / np.maximum(n, 1), np.nan)
    return h


def snp_density(dose, cols):
    """Per-site indicator that the group is polymorphic (has both alleles among called genotypes)."""
    alt, n = group_counts(dose, cols)
    return (alt > 0) & (alt < n)


def windowed_mean(pos, values, size, name):
    ok = ~np.isnan(values)
    df = pd.DataFrame(dict(start=window_index(pos[ok], size), v=values[ok]))
    return df.groupby('start').v.mean().rename(name).reset_index()


def sex_summary(pos, dose, is_male, size, min_called=4):
    """Window table: male/female heterozygosity, delta, SNP density, and male-female Fst. dose: (sites, samples)."""
    m = np.where(is_male)[0]; f = np.where(~np.asarray(is_male))[0]
    am, nm = group_counts(dose, m); af, nf = group_counts(dose, f)
    num, den = hudson_terms(am, nm, af, nf, min_called * 2)
    out = windowed_fst(pos, num, den, size)
    for tag, cols in (('M', m), ('F', f)):
        out = out.merge(windowed_mean(pos, het_by_group(dose, cols, min_called), size, f'het_{tag}'), on='start', how='left')
        poly = snp_density(dose, cols).astype(float)
        out = out.merge(windowed_mean(pos, poly, size, f'poly_{tag}'), on='start', how='left')
    out['het_diff_M_minus_F'] = out.het_M - out.het_F
    return out


def depth_ratio(depth, is_male):
    """depth: (windows, samples) mean depth. Normalise each sample by its own median, then log2 of the male/female mean ratio
    per window plus a Mann-Whitney p-value. Detects Y-specific sequence (ratio>0), X-only sequence (<0) and male-specific duplications."""
    norm = depth / np.nanmedian(depth, axis=0, keepdims=True)
    male = np.asarray(is_male, bool)
    rm, rf = np.nanmean(norm[:, male], 1), np.nanmean(norm[:, ~male], 1)
    lr = np.log2(np.maximum(rm, 1e-3) / np.maximum(rf, 1e-3))
    p = np.array([stats.mannwhitneyu(norm[i, male], norm[i, ~male]).pvalue if np.isfinite(norm[i]).all() else np.nan for i in range(len(norm))])
    return lr, p
