"""python -m sd_regions region_block --vcf V --qc <sexloci output>/sample_qc.tsv --region CHROM:START-END [--core START-END] [--outdir D] [options]

Is a block of sex-linked signal a real sex haplotype (an inversion, a Y region), or could it be made by low coverage, by inbreeding or by family structure that happens to line up with sex?
One pass over the VCF through sexloci's own reader (the same filters, depth masking and call-rate rule), then, in `--outdir`:

  summary.md            the numbers in words, with how to read them
  samples.tsv           per sample: sex, heterozygosity inside and outside the region, runs of homozygosity (F_ROH, and F_ROH_outside_region), inbreeding (F), the region's share in a run, the diagnostic-SNP heterozygous share
  roh_segments.tsv      the runs of homozygosity (>= 4 adjacent 250 kb windows with almost no heterozygotes) of every sample
  het_heatmap.png       samples x windows of the chromosome, coloured by heterozygosity; red bars = runs of homozygosity; the mean by sex above
  kinship.png, kinship.tsv, kinship_pairs.tsv     relatedness outside the region and whether it lines up with sex
  pca.png, pca_region.tsv, pca_genome.tsv         the structure of the fish from the SNPs in the region and from the rest of the genome, and how well each separates the sexes
  ld.png, ld_decay.tsv  linkage disequilibrium in the region and in the rest of the chromosome, for all fish, the males and the females
  diagnostic.png, diagnostic.tsv, diagnostic_subsample.tsv, diagnostic_leave_one_out.tsv    male-diagnostic SNPs under sexloci's rule and a strict one, at several depth cut-offs, in random subsets of the fish, with each fish left out in turn and without the fish that contradict the block
  subregions.tsv        the whole region, the --core and every --subregion: PCA AUC, diagnostic counts, heterozygosity by sex, females and males with >= 20% of the part inside a run of homozygosity
  null_scan.png, null_scan.tsv    the region PCA repeated on every window of the genome of the same length: how unusual is the region's AUC for sex in THESE fish (skip with --no-null-scan)
  chromosome_sex_het.tsv    male and female heterozygosity of every chromosome (is the region the only one where the sexes differ?)
  window_het.tsv, run_info.json

It does not decide whether there is an inversion. It tests the three alternative explanations that would make a block of male-diagnostic SNPs misleading, and shows what the block looks like."""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

from . import analysis as A
from . import plots as P
from . import reader as R

DEFAULT_SUBSAMPLES = ['9,7', '8,8']


def parse_range(text):
    a, b = text.replace(',', '').split('-')
    return int(a), int(b)


def load_samples(qc=None, sex_file=None):
    """(names, sex array of 'M'/'F', qc DataFrame or None). A sample with another sex value is left out."""
    if not qc and not sex_file:
        raise ValueError('give --qc (sexloci\'s sample_qc.tsv) or --sex-file (a table with the columns sample and sex)')
    path = qc or sex_file
    t = pd.read_csv(path, sep='\t', dtype=str)
    t.columns = [c.strip().lower() for c in t.columns]
    if 'sample' not in t.columns or 'sex' not in t.columns:
        raise ValueError(f'{path} needs the columns sample and sex; it has {list(t.columns)}')
    t = t[t.sex.str.strip().str.upper().isin(['M', 'F'])].copy()
    t['sex'] = t.sex.str.strip().str.upper()
    return t['sample'].tolist(), t['sex'].to_numpy(), t


def _fmt(x, d=3):
    return 'n/a' if x is None or (isinstance(x, float) and not np.isfinite(x)) else f'{x:.{d}f}'


def run(a, log=print):
    t0 = time.time()
    names, sex, qc = load_samples(a.qc, a.sex_file)
    region = R.parse_region(a.region)
    rc, rs, re_ = region
    core = parse_range(a.core) if a.core else None
    male, female = sex == 'M', sex == 'F'
    if male.sum() < 3 or female.sum() < 3:
        raise ValueError(f'need at least 3 males and 3 females, have {int(male.sum())} and {int(female.sum())}')
    os.makedirs(a.outdir, exist_ok=True)
    out = lambda f: os.path.join(a.outdir, f)
    exclude = [R.parse_region(x) for x in (a.kinship_exclude or [])]
    log(f'reading {a.vcf}: {len(names)} samples ({int(male.sum())} M, {int(female.sum())} F), region {a.region}')
    col = R.stream(a.vcf, names, region, a.config, a.min_dp, a.min_call_rate, a.window, a.min_maf, exclude, a.keep_region_chromosome, a.backend, log, keep_genome=not a.no_null_scan)
    n = len(names)
    log(f'{col.n_used:,} of {col.n_seen:,} sites passed the filters and the call-rate rule ({col.grm_out.n_sites:,} sites for relatedness)')
    size = a.window

    # ------------------------------------------------------------------ windows, heterozygosity, runs of homozygosity
    keys, called, het = col.windows.arrays()
    het_win = A.window_het(called, het, a.roh_min_sites)
    segs, callable_ = A.roh_segments(keys, called, het, size, a.roh_het, a.roh_min_sites, a.roh_min_windows)
    froh = A.froh(segs, callable_, size, n)
    out_idx = [i for i in range(len(keys)) if not (keys[i][0] == rc and keys[i][1] * size < re_ and (keys[i][1] + 1) * size > rs)]
    segs_out, callable_out = A.roh_segments([keys[i] for i in out_idx], called[out_idx], het[out_idx], size, a.roh_het, a.roh_min_sites, a.roh_min_windows)
    froh_out = A.froh(segs_out, callable_out, size, n)                      # inbreeding outside the region: not contaminated by the block itself
    reg = A.region_summary(keys, called, het, segs, callable_, rc, rs, re_, size, a.roh_het)
    inside = set(A.region_windows(keys, rc, rs, re_, size).tolist())
    outside = [i for i in range(len(keys)) if i not in inside]
    genome_het = het[outside].sum(axis=0) / np.maximum(called[outside].sum(axis=0), 1)
    chr_het = A.chrom_sex_het(keys, called, het, male, female)

    # ------------------------------------------------------------------ relatedness and inbreeding (outside the region)
    G = col.grm_out.grm()
    f_grm = np.diag(G) - 1
    f_het = col.grm_out.f_het()
    ksex = A.kinship_sex_test(G, sex, a.n_perm, a.seed)
    clusters = A.related_clusters(G, a.cluster_threshold)
    cl_of = {i: k for k, g in enumerate(clusters) for i in g}

    # ------------------------------------------------------------------ the sites of the region
    in_reg = (col.chr_pos >= rs) & (col.chr_pos < re_)
    pos_reg = col.chr_pos[in_reg]
    dos_reg = R.masked(col.chr_raw[in_reg], col.chr_dp[in_reg], a.min_dp)
    res = dict(n_region_sites=int(in_reg.sum()))
    pca = {}
    if in_reg.sum() >= 10:
        try:
            sc_r, var_r, used_r = A.pca_scores(dos_reg, a.min_maf)
            res['pca_region'] = A.best_pc_auc(sc_r, male, a.n_perm, a.seed) + (used_r,)
            pca['region'] = (sc_r, var_r)
        except ValueError as e:
            log(f'region PCA skipped: {e}')
    sc_g, var_g = A.pca_from_grm(G)
    res['pca_genome'] = A.best_pc_auc(sc_g, male, a.n_perm, a.seed)
    pca['genome'] = (sc_g, var_g)
    res['pc_corr'] = None
    if 'region' in pca:
        from scipy.stats import spearmanr
        res['pc_corr'] = max(abs(spearmanr(pca['region'][0][:, 0], sc_g[:, k])[0]) for k in range(sc_g.shape[1]))
    parts = [('whole region', rs, re_)] + ([('core', core[0], core[1])] if core else []) + [(None, *parse_range(x)) for x in (a.subregion or [])]
    subregion_rows = []
    for name, s0, s1 in parts:
        sel = (pos_reg >= s0) & (pos_reg < s1)
        row = dict(part=name or f'{rc}:{s0}-{s1}', start=s0, end=s1, sites=int(sel.sum()))
        if sel.sum() >= 10:
            try:
                sc, var, used = A.pca_scores(dos_reg[sel], a.min_maf)
                auc, p, pc = A.best_pc_auc(sc, male, a.n_perm, a.seed)
                row.update(pca_auc=auc, pca_p=p, pca_pc=pc)
            except ValueError:
                pass
            row.update(diag_sexloci=int(A.diag_sexloci(dos_reg[sel], male, female, a.diag_tol).sum()), diag_strict=int(A.diag_strict(dos_reg[sel], male, female).sum()))
        rs_ = A.region_summary(keys, called, het, segs, callable_, rc, s0, s1, size, a.roh_het)
        if rs_['n_windows']:
            row.update(het_male=float(np.nanmean(rs_['het'][male])), het_female=float(np.nanmean(rs_['het'][female])),
                       roh_share_male=float(np.nanmean(rs_['roh_share'][male])), roh_share_female=float(np.nanmean(rs_['roh_share'][female])),
                       females_roh20=int((rs_['roh_share'][female] >= 0.2).sum()), males_roh20=int((rs_['roh_share'][male] >= 0.2).sum()))
        subregion_rows.append(row)

    # ------------------------------------------------------------------ linkage disequilibrium: the region against the rest of the chromosome
    edges = [0, 100_000, 250_000, 500_000, 1_000_000, 2_000_000, 5_000_000, 10_000_000, 1e12]
    dos_all = R.masked(col.chr_raw, col.chr_dp, a.min_dp)
    ctrl_sel = (col.chr_pos < rs) | (col.chr_pos >= re_)
    ld_rows, ld_mats = [], []
    for label, d, p in (('region', dos_reg, pos_reg), ('rest of the chromosome (flanks joined)', dos_all[ctrl_sel], col.chr_pos[ctrl_sel])):
        if len(p) < 20:
            continue
        for grp, mask in (('all fish', np.ones(n, bool)), ('males', male), ('females', female)):
            ds, ps = A.thin_sites(d[:, mask], p, a.ld_min_maf, a.ld_sites)
            if len(ps) < 10:                                                  # nothing varies among these fish (e.g. males all heterozygous, females all alike): say so
                ld_rows.append(dict(set=label, group=grp, n_sites=len(ps)))
                ld_mats.append((f'{label}: {grp}', None, np.zeros(0)))
                continue
            r2 = A.ld_matrix(ds)
            dec = A.ld_by_distance(r2, ps, edges)
            ld_rows.append(dict(set=label, group=grp, n_sites=len(ps), **{f'r2_{int(lo / 1000)}-{("inf" if hi > 1e11 else str(int(hi / 1000)))}kb': v for lo, hi, v in zip(edges[:-1], edges[1:], dec)}))
            ld_mats.append((f'{label}: {grp} ({len(ps)} SNPs)', r2, ps))

    # ------------------------------------------------------------------ the diagnostic SNPs
    diag_loose = A.diag_sexloci(dos_reg, male, female, a.diag_tol)
    diag_str = A.diag_strict(dos_reg, male, female)
    thresholds = [t for t in (3, 5, 8) if t >= a.min_dp] or [a.min_dp]
    by_depth = A.diag_by_depth(col.chr_raw[in_reg], col.chr_dp[in_reg], male, female, thresholds, a.diag_tol)
    subs = {}
    for text in (a.subsample or DEFAULT_SUBSAMPLES):
        nm, nf = (int(x) for x in text.split(','))
        if nm <= male.sum() and nf <= female.sum() and (nm, nf) != (int(male.sum()), int(female.sum())):
            subs[(nm, nf)] = A.diag_subsample(dos_reg, male, female, nm, nf, a.n_draws, a.diag_tol, a.seed)
    diag_het = A.diag_sample_het(dos_reg, diag_loose) if diag_loose.any() else np.full(n, np.nan)
    loo = A.diag_leave_one_out(dos_reg, male, female, a.diag_tol)
    exc = A.exception_fish(diag_het, male, female)
    without = A.diag_without(dos_reg, male, female, exc, a.diag_tol) if len(exc) else None
    res['exceptions'] = [names[i] for i in exc]
    res['without'] = without

    # ------------------------------------------------------------------ the null: the same region PCA on every window of the genome
    null_rows = []
    if not a.no_null_scan and col.genome:
        length = re_ - rs
        null_rows = A.null_scan(col.genome, male, length, a.scan_step or max(length // 4, size), a.min_maf, a.scan_min_sites)
        for r in null_rows:
            r['dhet'] = A.scan_dhet(keys, called, het, male, female, r['chrom'], r['start'], r['end'], size)
            r['overlaps_region'] = bool(r['chrom'] == rc and r['start'] < re_ and r['end'] > rs)
        log(f'null scan: {len(null_rows)} windows of {length / 1e6:.0f} Mb')
    res['null'] = null_rows

    # ------------------------------------------------------------------ tables
    longest = np.zeros(n)
    nroh = np.zeros(n, int)
    for g in segs:
        longest[g['sample']] = max(longest[g['sample']], (g['end'] - g['start']) / 1e6)
        nroh[g['sample']] += 1
    st = pd.DataFrame(dict(sample=names, sex=sex, genome_het_outside_region=genome_het, region_het=reg['het'], region_homozygous_share=reg['hom_share'], region_in_roh_share=reg['roh_share'],
                           F_het=f_het, F_grm=f_grm, F_roh=froh, F_roh_outside_region=froh_out, n_roh=nroh, longest_roh_mb=longest, diag_het_share=diag_het, relatedness_cluster=[cl_of[i] for i in range(n)]))
    if qc is not None and 'mean_depth' in qc.columns:
        st['mean_depth'] = pd.to_numeric(qc.set_index('sample').loc[names, 'mean_depth'], errors='coerce').to_numpy()
    st.round(5).to_csv(out('samples.tsv'), sep='\t', index=False)
    pd.DataFrame([dict(sample=names[g['sample']], sex=sex[g['sample']], chrom=g['chrom'], start=g['start'], end=g['end'], length_mb=(g['end'] - g['start']) / 1e6, n_windows=g['n_windows']) for g in segs],
                 columns=['sample', 'sex', 'chrom', 'start', 'end', 'length_mb', 'n_windows']).to_csv(out('roh_segments.tsv'), sep='\t', index=False)
    wh = pd.DataFrame(het_win, columns=names)
    wh.insert(0, 'start', [k[1] * size for k in keys]); wh.insert(0, 'chrom', [k[0] for k in keys])
    wh.round(4).to_csv(out('window_het.tsv'), sep='\t', index=False)
    pd.DataFrame(G, index=names, columns=names).round(4).to_csv(out('kinship.tsv'), sep='\t')
    iu = np.triu_indices(n, 1)
    pairs = pd.DataFrame(dict(a=[names[i] for i in iu[0]], b=[names[j] for j in iu[1]], sex_a=sex[iu[0]], sex_b=sex[iu[1]], relatedness=G[iu])).sort_values('relatedness', ascending=False)
    pairs.round(4).to_csv(out('kinship_pairs.tsv'), sep='\t', index=False)
    pd.DataFrame([dict(chrom=c, male_het=v[0], female_het=v[1], male_minus_female=v[2]) for c, v in chr_het.items()]).round(4).to_csv(out('chromosome_sex_het.tsv'), sep='\t', index=False)
    for nm, key in (('pca_region.tsv', 'region'), ('pca_genome.tsv', 'genome')):
        if key in pca:
            sc = pca[key][0]
            pd.DataFrame({'sample': names, 'sex': sex, **{f'PC{k + 1}': sc[:, k] for k in range(sc.shape[1])}}).round(5).to_csv(out(nm), sep='\t', index=False)
    pd.DataFrame(ld_rows).round(4).to_csv(out('ld_decay.tsv'), sep='\t', index=False)
    if null_rows:
        pd.DataFrame(null_rows).round(4).to_csv(out('null_scan.tsv'), sep='\t', index=False)
    drows = [dict(rule=rule, min_depth=t, count=v[i]) for t, v in by_depth.items() for i, rule in enumerate(('sexloci', 'strict'))]
    drows.append(dict(rule='sexloci', min_depth=f'{a.min_dp} (as run), all fish', count=int(diag_loose.sum())))
    if without is not None:
        drows += [dict(rule='sexloci', min_depth=f'{a.min_dp}, without {", ".join(names[i] for i in exc)}', count=without[0]), dict(rule='strict', min_depth=f'{a.min_dp}, without {", ".join(names[i] for i in exc)}', count=without[1])]
    pd.DataFrame(drows).to_csv(out('diagnostic.tsv'), sep='\t', index=False)
    pd.DataFrame([dict(males=k[0], females=k[1], draws=len(v), median=float(np.median(v)), p05=float(np.quantile(v, .05)), p95=float(np.quantile(v, .95)), all_fish=int(diag_loose.sum())) for k, v in subs.items()]
                 ).round(1).to_csv(out('diagnostic_subsample.tsv'), sep='\t', index=False)
    pd.DataFrame(dict(sample=names, sex=sex, diag_het_share=diag_het, sexloci_count_without_this_fish=loo[:, 0], strict_count_without_this_fish=loo[:, 1])).round(4).to_csv(out('diagnostic_leave_one_out.tsv'), sep='\t', index=False)
    pd.DataFrame(subregion_rows).round(4).to_csv(out('subregions.tsv'), sep='\t', index=False)

    # ------------------------------------------------------------------ figures
    P.het_heatmap(out('het_heatmap.png'), keys, het_win, names, sex, rc, region, size, core, segs, f'Heterozygosity along {rc} ({a.region.split(":")[1]} marked)')
    P.kinship_plot(out('kinship.png'), G, names, sex, f'within-female {_fmt(ksex["within_f"])}, within-male {_fmt(ksex["within_m"])}, between sexes {_fmt(ksex["between"])}; permutation p = {_fmt(ksex["p"], 4)}')
    if 'region' in pca:
        ra = res.get('pca_region')
        P.pca_plot(out('pca.png'), pca['region'][0], pca['region'][1], sc_g, var_g, names, sex, f'best PC{ra[2]} separates the sexes: AUC {_fmt(ra[0], 2)} (p = {_fmt(ra[1], 4)})',
                   f'best PC{res["pca_genome"][2]}: AUC {_fmt(res["pca_genome"][0], 2)} (p = {_fmt(res["pca_genome"][1], 4)})')
    if ld_mats:
        P.ld_plot(out('ld.png'), ld_mats, 'Linkage disequilibrium (r$^2$): the region, then the rest of the chromosome; all fish, males, females')
    P.diag_plot(out('diagnostic.png'), by_depth, subs, int(diag_loose.sum()), n, names, sex, loo)
    if null_rows and res.get('pca_region'):
        P.null_plot(out('null_scan.png'), null_rows, res['pca_region'][0], rs, re_, rc)

    # ------------------------------------------------------------------ summary
    summary = write_summary(a, names, sex, st, segs, ksex, clusters, res, pairs, chr_het, by_depth, subs, diag_loose, diag_str, ld_rows, subregion_rows, col, rc, rs, re_, core, loo)
    open(out('summary.md'), 'w').write(summary)
    json.dump(dict(args={k: v for k, v in vars(a).items() if k != 'func'}, n_samples=n, sites_seen=col.n_seen, sites_used=col.n_used, seconds=round(time.time() - t0, 1)), open(out('run_info.json'), 'w'), indent=1, default=str)
    log(f'done in {time.time() - t0:.0f} s: {a.outdir}')
    return dict(loo=loo, parts=subregion_rows, st=st, ksex=ksex, res=res, by_depth=by_depth, subs=subs, diag_loose=int(diag_loose.sum()), diag_strict=int(diag_str.sum()), segs=segs, summary=summary, clusters=clusters, ld_rows=ld_rows)


def write_summary(a, names, sex, st, segs, ksex, clusters, res, pairs, chr_het, by_depth, subs, diag_loose, diag_str, ld_rows, subregion_rows, col, rc, rs, re_, core, loo):
    male, female = sex == 'M', sex == 'F'
    L = []
    w = L.append
    nm, nf = int(male.sum()), int(female.sum())
    w(f'# Region block analysis: {rc}:{rs:,}-{re_:,}\n')
    w(f'{nm} males and {nf} females from `{os.path.basename(a.vcf)}`; genotypes at depth >= {a.min_dp}, sites called in >= {a.min_call_rate:.0%} of the fish ({col.n_used:,} of {col.n_seen:,} sites). '
      f'This tests three explanations that would make a block of male-diagnostic SNPs misleading (low coverage, inbreeding, family structure aligned with sex) and shows what the block looks like. It does not by itself show an inversion.\n')
    # --- homozygosity
    fm, ff = st[st.sex == 'M'], st[st.sex == 'F']
    w('## 1. Heterozygosity and runs of homozygosity (inbreeding)\n')
    w(f'| | males | females |\n|---|---|---|\n| heterozygosity outside the region (mean) | {_fmt(fm.genome_het_outside_region.mean())} | {_fmt(ff.genome_het_outside_region.mean())} |\n'
      f'| heterozygosity inside the region | {_fmt(fm.region_het.mean())} | {_fmt(ff.region_het.mean())} |\n'
      f'| F_ROH outside the region (share of the genome in a run of homozygosity), mean (range) | {_fmt(fm.F_roh_outside_region.mean())} ({_fmt(fm.F_roh_outside_region.min())}-{_fmt(fm.F_roh_outside_region.max())}) | {_fmt(ff.F_roh_outside_region.mean())} ({_fmt(ff.F_roh_outside_region.min())}-{_fmt(ff.F_roh_outside_region.max())}) |\n'
      f'| inbreeding F from heterozygosity, mean | {_fmt(fm.F_het.mean())} | {_fmt(ff.F_het.mean())} |\n'
      f'| share of the region inside a run of homozygosity, mean (times the same fish\' F_ROH outside it) | {_fmt(fm.region_in_roh_share.mean())} (x{_fmt(fm.region_in_roh_share.mean() / fm.F_roh_outside_region.mean(), 1)}) | {_fmt(ff.region_in_roh_share.mean())} (x{_fmt(ff.region_in_roh_share.mean() / ff.F_roh_outside_region.mean(), 1)}) |\n'
      f'| fish with >= 20% of the region inside a run of homozygosity | {int((fm.region_in_roh_share >= 0.2).sum())} of {nm} | {int((ff.region_in_roh_share >= 0.2).sum())} of {nf} |\n')
    w(f'A run of homozygosity here is {a.roh_min_windows} or more adjacent {a.window // 1000} kb windows in which at most {a.roh_het:.0%} of the calls are heterozygous (a threshold, not a model; resolution one window). '
      f'There are {len(segs)} runs in all; the longest is {st.longest_roh_mb.max():.1f} Mb ({st.loc[st.longest_roh_mb.idxmax(), "sample"]}).\n')
    w('How to read it: compare the share of the region inside a run with the F_ROH the same fish have outside it (the multiplier in the table): a multiplier well above 1 in one sex only, with F_ROH outside the region alike in both sexes, means the homozygosity belongs to the region and not to inbreeding in general. A block where the females (or the males) are homozygous in **many** fish while F_ROH is low elsewhere points at the block itself (a haplotype of low diversity, an inversion '
      'homozygote, a sex chromosome). High F_ROH and F in both sexes everywhere points at inbreeding or a small founder group, and then homozygosity in a block alone is weak evidence. '
      'The picture is `het_heatmap.png`: dark runs in a row are runs of homozygosity, red bars mark the called runs.\n')
    # --- relatedness
    w('## 2. Relatedness and family structure\n')
    w(f'Relatedness is estimated from the genome outside {rc if not a.keep_region_chromosome else "the region"} ({col.grm_out.n_sites:,} SNPs).\n')
    w(f'| mean relationship | value |\n|---|---|\n| between two females | {_fmt(ksex["within_f"])} |\n| between two males | {_fmt(ksex["within_m"])} |\n| between a male and a female | {_fmt(ksex["between"])} |\n')
    w(f'Same-sex pairs minus mixed-sex pairs = {_fmt(ksex["stat"], 4)}; permutation p = {_fmt(ksex["p"], 4)} (small = relatedness lines up with sex, which would confound any sex contrast).\n')
    cl_big = [c for c in clusters if len(c) > 1]
    if cl_big:
        w(f'Groups of related fish (relationship >= {a.cluster_threshold}): ' + '; '.join(f'{len(c)} fish ({sum(sex[i] == "M" for i in c)} M, {sum(sex[i] == "F" for i in c)} F)' for c in cl_big) + '. ')
        w('A group made of one sex only is the warning sign.\n')
    else:
        w(f'No two fish reach a relationship of {a.cluster_threshold}.\n')
    top = pairs.head(5)
    w('Most related pairs: ' + '; '.join(f'{r.a}-{r.b} ({r.sex_a}{r.sex_b}) {r.relatedness:.2f}' for r in top.itertuples()) + '. (The estimate is biased by about -1/(n-1) off the diagonal because the allele frequencies come from the same fish: compare pairs, '
      'do not read pedigree coefficients.)\n')
    # --- pca
    w('## 3. Structure: does the block split the fish in a way the rest of the genome does not?\n')
    if res.get('pca_region'):
        ra, rg = res['pca_region'], res['pca_genome']
        w(f'| | best PC | AUC for sex | permutation p |\n|---|---|---|---|\n| SNPs in the region ({ra[3]:,} SNPs) | PC{ra[2]} | {_fmt(ra[0], 2)} | {_fmt(ra[1], 4)} |\n'
          f'| rest of the genome | PC{rg[2]} | {_fmt(rg[0], 2)} | {_fmt(rg[1], 4)} |\n')
        w(f'The region\'s PC1 correlates with the genome-wide PCs at |rho| = {_fmt(res["pc_corr"], 2)}. The p-values choose the best of three components, so they are optimistic.\n')
        w('How to read it: sexes separate in the region but not in the genome-wide PCA = a locus-specific difference, as a sex haplotype would give. Sexes separate in both = population structure aligned with sex. '
          'Groups in the region that follow neither sex nor the genome-wide structure would be an inversion (or other haplotype) genotype that is not sex-linked.\n')
    else:
        w('Too few SNPs in the region for a PCA.\n')
    # --- ld
    w('## 4. Linkage disequilibrium\n')
    if ld_rows:
        t = pd.DataFrame(ld_rows)
        cols = [c for c in t.columns if c.startswith('r2_')]
        pick = [c for c in cols if c in ('r2_1000-2000kb', 'r2_2000-5000kb', 'r2_5000-10000kb')]
        w('Mean r^2 between SNPs far apart (`ld.png` shows the matrices; `ld_decay.tsv` all distances):\n')
        w('| set | group | SNPs | ' + ' | '.join(c.replace('r2_', '') for c in pick) + ' |\n|---|---|---|' + '---|' * len(pick) + '\n' + '\n'.join(f'| {r["set"]} | {r["group"]} | {int(r["n_sites"])} | ' + ' | '.join(_fmt(float(r[c]), 3) for c in pick) + ' |' for _, r in t.iterrows()) + '\n')
        nn = len(names)
        w(f'Baseline: SNPs with no linkage at all give r^2 of about 1/n = {1 / nn:.2f} for all {nn} fish, {1 / nm:.2f} for the {nm} males and {1 / nf:.2f} for the {nf} females. If the rest of the chromosome is far above that too, the fish share long haplotypes everywhere (a small or structured population), and then the region\'s LD cannot be told from the background.\n')
        w('How to read it: an inversion or a non-recombining haplotype shows as a block of high r^2 across Mb, in the fish that carry both arrangements. A block in the males but not the females fits a Y haplotype; in neither, '
          'LD that is no higher than elsewhere on the chromosome argues against a long haplotype.\n')
    # --- diagnostic
    w('## 5. The male-diagnostic SNPs\n')
    w(f'In the region: {int(diag_loose.sum()):,} under sexloci\'s rule (at least {1 - a.diag_tol:.0%} of the males heterozygous, at most {a.diag_tol:.0%} of the females, females fixed for one allele) and {int(diag_str.sum()):,} under a strict rule (every called male heterozygous, every called female homozygous).\n')
    w('| minimum depth | sexloci rule | strict |\n|---|---|---|\n' + '\n'.join(f'| {t} | {v[0]:,} | {v[1]:,} |' for t, v in by_depth.items()) + '\n')
    if subs:
        w('Random subsets of the fish (`diagnostic_subsample.tsv`): ' + '; '.join(f'{k[0]} males + {k[1]} females: median {int(np.median(v))} (5-95%: {int(np.quantile(v, .05))}-{int(np.quantile(v, .95))})' for k, v in sorted(subs.items())) + f'; all {nm}+{nf}: {int(diag_loose.sum()):,}. '
          'The rule allows one exception per sex with 12 fish and none with 9 or 7, so the count is **not comparable between runs of different size**; use the heterozygosity and Fst profiles for that.\n')
    dh = st.diag_het_share
    if np.isfinite(dh).any():
        exc = st[st['sample'].isin(res.get('exceptions', []))]
        w('Heterozygous share at the diagnostic sites: males ' + _fmt(dh[male].mean(), 2) + ' (range ' + _fmt(dh[male].min(), 2) + '-' + _fmt(dh[male].max(), 2) + '), females ' + _fmt(dh[female].mean(), 2) + ' (range ' + _fmt(dh[female].min(), 2) + '-' + _fmt(dh[female].max(), 2) + '). '
          + ('Exceptions: ' + ', '.join(f'{r.sample} ({r.sex}, {r.diag_het_share:.2f})' for r in exc.itertuples()) + '.\n' if len(exc) else 'No fish is an exception.\n'))
    if len(res.get('exceptions', [])):
        wo = res['without']
        w(f'**The fish that contradict the block decide the count.** Without {", ".join(res["exceptions"])}: {wo[0]:,} SNPs under sexloci\'s rule and {wo[1]:,} under the strict rule (with all fish: {int(diag_loose.sum()):,} and {int(diag_str.sum()):,}). '
          'Leaving one fish out at a time (`diagnostic_leave_one_out.tsv`), the strict count is highest without ' + ', '.join(f'{names[i]} ({int(loo[i, 1]):,})' for i in np.argsort(-loo[:, 1])[:3]) + '. '
          'A fish that is homozygous at SNPs where the other fish of its sex are heterozygous (or the reverse) is either a recombinant, a different haplotype class, or a sex that was recorded wrongly: check its other markers (TE1, TE2) and its record.\n')
    w('Low coverage can only hide heterozygotes, so it works against a male-diagnostic count; it could create one only if the females\' heterozygotes were all missed at the same sites, which needs very low female depth there.\n')
    if subregion_rows:
        w('## 6. Parts of the region\n')
        w('| part | SNPs | best PC AUC (p) | diagnostic: sexloci / strict | male het | female het | females / males with >= 20% in a run |\n|---|---|---|---|---|---|---|\n' + '\n'.join(
            f'| {r["part"]} ({r["start"] / 1e6:.1f}-{r["end"] / 1e6:.1f} Mb) | {r["sites"]:,} | {_fmt(r.get("pca_auc"), 2)} ({_fmt(r.get("pca_p"), 4)}) | {r.get("diag_sexloci", "n/a")} / {r.get("diag_strict", "n/a")} | '
            f'{_fmt(r.get("het_male"), 2)} | {_fmt(r.get("het_female"), 2)} | {r.get("females_roh20", "n/a")} / {r.get("males_roh20", "n/a")} |' for r in subregion_rows) + '\n')
    # --- chromosomes
    w('## 7. Where else do the sexes differ in heterozygosity?\n')
    top = sorted(chr_het.items(), key=lambda kv: -abs(kv[1][2]))[:6]
    w('| chromosome | male het | female het | male - female |\n|---|---|---|---|\n' + '\n'.join(f'| {c} | {v[0]:.3f} | {v[1]:.3f} | {v[2]:+.3f} |' for c, v in top) + '\n')
    w('A real sex locus is expected on one chromosome; several chromosomes with large differences of both signs point at structure or at chance blocks in a small sample.\n')
    nulls = [r for r in res.get('null', []) if not r['overlaps_region']]
    if nulls and res.get('pca_region'):
        auc_r = res['pca_region'][0]
        aucs = np.array([r['auc'] for r in nulls]); dh = np.array([r['dhet'] for r in nulls if np.isfinite(r['dhet'])])
        reg_dh = next((r.get('het_male', np.nan) - r.get('het_female', np.nan) for r in subregion_rows if r['part'] == 'whole region'), np.nan)
        w('## 8. Is this unusual for these fish? (the same analysis on every other stretch of the genome)\n')
        w(f'The region PCA was repeated on {len(nulls)} windows of the same length ({(re_ - rs) / 1e6:.0f} Mb, every {(a.scan_step or max((re_ - rs) // 4, a.window)) / 1e6:.1f} Mb) that do not overlap the region (`null_scan.tsv`, `null_scan.png`). '
          f'The best of three PCs separates the sexes with an AUC of {np.median(aucs):.2f} at the median and {np.quantile(aucs, .95):.2f} at the 95th percentile of these windows; the region reaches {auc_r:.2f}, and {int((aucs >= auc_r - 1e-9).sum())} of the {len(aucs)} windows reach it or more '
          f'({(aucs >= auc_r - 1e-9).mean():.1%}). ' + (f'Male minus female heterozygosity: median {np.median(dh):+.3f}, 95th percentile {np.quantile(dh, .95):+.3f}, the region {reg_dh:+.3f} ({int((dh >= reg_dh).sum())} of {len(dh)} windows at or above it). ' if len(dh) and np.isfinite(reg_dh) else '')
          + 'The windows overlap, so read this as an approximate rank, not as a p-value; and this sample has long shared haplotypes everywhere, which is why a chance association with sex is possible at all.\n')
    w('## What this does not show\n')
    w('**Every p-value above is descriptive.** The region was chosen because it stood out in these same fish, so a test on the same fish is optimistic. Only an independent sample (another pool, or a pedigree with both sexes in each family) can test it.\n')
    w('That the block is an **inversion**: that needs breakpoints. If the block survives the checks above, the next steps are (1) male-only discordant read pairs and split reads at its edges in the chimeric BAMs, '
      '(2) the same region in an independent pedigree with both sexes in every family (for the YH lineage the chr10 Y inversion is a known system: the same coordinates would mean the same rearrangement, different ones a nested or new one), '
      'and (3) a male assembly. Nor does it replace replication: with one pool of unrelated fish the sexloci tiers cannot go above C.\n')
    return '\n'.join(L)


def build_parser():
    ap = argparse.ArgumentParser(prog='python -m sd_regions region_block', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--qc', help="sexloci's sample_qc.tsv (columns sample and sex); the samples and sexes analysed")
    ap.add_argument('--sex-file', help='alternatively a table with the columns sample and sex')
    ap.add_argument('--region', required=True, help='CHROM:START-END, the block to examine (and to leave out of the genome-wide relatedness)')
    ap.add_argument('--core', help='START-END: the called block, drawn on the figures')
    ap.add_argument('--subregion', action='append', help='START-END of a part of the region, analysed separately (repeatable)')
    ap.add_argument('--outdir', help='default: ~/Temp/SexDeterminationAnalysis/region_block_<chrom>_<start>-<end>')
    ap.add_argument('--config', help="a sexloci YAML (the same file as for the run)")
    ap.add_argument('--min-dp', type=int, default=3, help='genotypes below this depth are missing (the value of the sexloci run)')
    ap.add_argument('--min-call-rate', type=float, default=0.8)
    ap.add_argument('--window', type=int, default=250_000)
    ap.add_argument('--min-maf', type=float, default=0.05, help='minor allele frequency for relatedness and PCA')
    ap.add_argument('--roh-het', type=float, default=0.03, help='a window is homozygous when at most this share of its calls are heterozygous')
    ap.add_argument('--roh-min-sites', type=int, default=30)
    ap.add_argument('--roh-min-windows', type=int, default=4, help='adjacent homozygous windows that make a run')
    ap.add_argument('--cluster-threshold', type=float, default=0.2, help='relationship that joins two fish in a group')
    ap.add_argument('--kinship-exclude', action='append', help='CHROM:START-END left out of the genome-wide relatedness too (repeatable); the region\'s chromosome is always left out')
    ap.add_argument('--keep-region-chromosome', action='store_true', help='use the rest of the region\'s chromosome for relatedness')
    ap.add_argument('--ld-sites', type=int, default=300)
    ap.add_argument('--ld-min-maf', type=float, default=0.15)
    ap.add_argument('--diag-tol', type=float, default=0.1, help="sexloci's diag_tol")
    ap.add_argument('--subsample', action='append', help='M,F: count the diagnostic SNPs in random subsets of this many males and females (repeatable; default 9,7 and 8,8)')
    ap.add_argument('--n-draws', type=int, default=200)
    ap.add_argument('--no-null-scan', action='store_true', help='skip the scan of every window of the genome (it keeps every genotype of the genome in memory, one byte each: 3 million SNPs x 24 fish is about 80 MB)')
    ap.add_argument('--scan-step', type=int, default=None, help='distance between the windows of the null scan (default a quarter of the region)')
    ap.add_argument('--scan-min-sites', type=int, default=2000, help='windows with fewer SNPs are skipped')
    ap.add_argument('--n-perm', type=int, default=5000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--backend', choices=['auto', 'builtin', 'allel'], default=None, help="sexloci's VCF backend (default: allel when installed)")
    return ap


def main(argv=None):
    a = build_parser().parse_args(argv)
    if not a.outdir:
        c, s, e = R.parse_region(a.region)
        a.outdir = os.path.expanduser(f'~/Temp/SexDeterminationAnalysis/region_block_{c}_{s}-{e}')
    try:
        run(a)
    except ValueError as e:
        print(f'error: {e}', file=sys.stderr)
        return 2
    return 0
