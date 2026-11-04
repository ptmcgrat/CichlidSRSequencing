"""Rewritten pipeline (v2): deep-sample joint genotyping -> informative sites -> BC allele counts -> site-level HMM.
ALLELE-AWARE: the reference is Mzebra, so MC/YH alleles are polarized per site (MC-fixed allele vs YH-specific allele).
PEDIGREE: F1 002 & 003 are full sibs (YH_008_m x MC-1R6R-f); F1 010-017 are a different family (YH_005_m x MC-3R6R-f)."""
import os
import numpy as np
import pandas as pd
from ...common import sh, parallel
from ...bc_ancestry_hmm import forward_backward, fit_sample, with_grid, call_breakpoints
from ... import callset as cs


def deep_samples(cfg):
    s = [cfg['yh_sire'], cfg['f1_dam']] + cfg['bc_dams'] + cfg['mc_panel'] + list(cfg['f1_sire_names'].values())
    if cfg.get('use_extra_family_in_sites'):
        s += [cfg['extra_f1_family']['yh_sire']] + cfg['extra_f1_family']['f1']
    return list(dict.fromkeys(s))


# ------------------------------------------------------------------ 0/1. GVCFs and joint genotyping
def _gvcf_one(args):
    proj, s = args
    return proj.ensure_gvcf(s)


def joint_call(proj, samples, outdir, name):
    """Per-chromosome CombineGVCFs + GenotypeGVCFs, then concat -> <outdir>/<name>.vcf.gz (+tbi)."""
    out = f'{outdir}{name}.vcf.gz'
    if os.path.exists(out + '.tbi'):
        print('exists:', out); return out
    os.makedirs(outdir + 'per_chrom/', exist_ok=True)
    res, errs = parallel(_gvcf_one, [(proj, s) for s in samples], proj.workers('gvcf'), 'gvcf')
    assert not errs, errs
    gvcfs = [res[i] for i in range(len(samples))]
    items = [(proj, c, gvcfs, outdir, name) for c in proj.chroms]
    r, e = parallel(_joint_chrom, items, proj.workers('joint'), 'joint-genotype')
    assert not e, e
    files = [f'{outdir}per_chrom/{name}.{c}.vcf.gz' for c in proj.chroms]
    sh(['bcftools', 'concat', '-Oz', '-o', out] + files)
    sh(['tabix', '-f', '-p', 'vcf', out])
    return out


def _joint_chrom(args):
    proj, chrom, gvcfs, outdir, name = args
    comb = f"{proj.dirs['tmp']}{name}.comb.{chrom}.g.vcf.gz"
    out = f'{outdir}per_chrom/{name}.{chrom}.vcf.gz'
    if os.path.exists(out + '.tbi'):
        return out
    mem = proj.cfg.get('gatk_mem_gb', 8)
    sh(['gatk', '--java-options', f'-Xmx{mem}g', 'CombineGVCFs', '-R', proj.genome_fa, '-L', chrom, '-O', comb] +
       [x for g in gvcfs for x in ('-V', g)])
    sh(['gatk', '--java-options', f'-Xmx{mem}g', 'GenotypeGVCFs', '-R', proj.genome_fa, '-V', comb, '-L', chrom, '-O', out])
    for f in (comb, comb + '.tbi'):
        if os.path.exists(f):
            os.remove(f)
    return out


def joint_genotype(proj):
    joint_call(proj, deep_samples(proj.cfg), proj.dirs['v2'], 'deep_joint')


# ------------------------------------------------------------------ 2. informative sites (allele-polarized)
def define_sites(proj):
    """Informative site (per BC sire): all deep MC samples carry (almost) only ONE allele (the MC allele); the F1 sire is a
    balanced het carrying the OTHER allele (the YH allele, necessarily YH-derived because the F1 dam is MC); YH_008 carries that
    YH allele; depth in a normal window; if YH_008 is homozygous for the YH allele the other full-sib sire must be het too;
    optional repeat-mask BED. Works for either polarity (MC allele can be REF or ALT in Mzebra)."""
    import allel
    cfg, v = proj.cfg, proj.cfg['v2']
    sires = list(cfg['f1_sire_names'].values())
    os.makedirs(proj.dirs['v2'] + 'sites', exist_ok=True)
    fam_samples = cs.load_callsets()['callsets']['bc_family']['samples']
    fam_vcf = cs.vcf_path(proj, 'bc_family')                       # 7-sample GATK callset (see config/callsets.yaml)
    joint = f"{proj.dirs['v2']}bc_family.xcheck.vcf.gz"            # GATK sites that bcftools recounts agree with
    if not os.path.exists(joint + '.tbi'):
        print('bcftools cross-check dropped', cs.crosscheck(proj, fam_vcf, fam_samples, joint), 'sites')
    panel = list(dict.fromkeys(cfg['bc_dams'] + [cfg['f1_dam']]))   # the MC females inside the family VCF
    keep = {s: [] for s in sires}; attr = []
    for chrom in proj.chroms:
        d = allel.read_vcf(joint, region=chrom, fields=['variants/POS', 'variants/REF', 'variants/ALT', 'calldata/GT', 'calldata/AD', 'samples'])
        if d is None:
            continue
        ix = {s: i for i, s in enumerate(d['samples'])}
        GT = allel.GenotypeArray(d['calldata/GT'])
        AD = np.clip(d['calldata/AD'].astype(float), 0, None)
        ref, alt = AD[:, :, 0], AD[:, :, 1]; DP = ref + alt
        A = d['variants/ALT']; R = d['variants/REF']
        biallelic = (A[:, 1] == '') if A.shape[1] > 1 else np.ones(len(A), bool)
        snp = (np.char.str_len(R.astype(str)) == 1) & (np.char.str_len(A[:, 0].astype(str)) == 1) & biallelic
        meandp = DP[snp].mean(0)
        okdp = lambda s: (DP[:, ix[s]] >= v['min_dp_frac'] * meandp[ix[s]]) & (DP[:, ix[s]] <= v['max_dp_frac'] * meandp[ix[s]])
        ab = lambda s: alt[:, ix[s]] / np.maximum(DP[:, ix[s]], 1)
        het = lambda s: GT[:, ix[s]].is_het() & (ab(s) >= v['ab_lo']) & (ab(s) <= v['ab_hi']) & okdp(s)
        p_ref = sum(ref[:, ix[s]] for s in panel); p_alt = sum(alt[:, ix[s]] for s in panel); p_dp = p_ref + p_alt
        p_n = sum((DP[:, ix[s]] >= v['panel_min_dp']).astype(int) for s in panel)
        fixed_ref = (p_alt / np.maximum(p_dp, 1) <= v['panel_max_alt_frac'])        # MC carries REF only
        fixed_alt = (p_ref / np.maximum(p_dp, 1) <= v['panel_max_alt_frac'])        # MC carries ALT only
        panel_ok = (p_n >= v['panel_min_samples']) & (p_dp > 0) & (fixed_ref ^ fixed_alt)
        mc_is_ref = fixed_ref
        yh = cfg['yh_sire']
        yh_g = GT[:, ix[yh]]
        yh_carries_yh_allele = np.where(mc_is_ref, yh_g.is_het() | yh_g.is_hom_alt(), yh_g.is_het() | yh_g.is_hom_ref()) & okdp(yh)
        yh_hom = np.where(mc_is_ref, yh_g.is_hom_alt(), yh_g.is_hom_ref())
        for sire in sires:
            other = [x for x in sires if x != sire][0]
            steps = [('snp', snp), ('sire_het_balanced', het(sire)), ('mc_panel_fixed', panel_ok),
                     ('yh_carries_yh_allele', yh_carries_yh_allele), ('sib_consistent', (~yh_hom) | het(other))]
            m = np.ones(len(snp), bool); row = dict(chrom=chrom, sire=sire, start=len(m))
            for name, f in steps:
                m &= f; row[name] = int(m.sum())
            attr.append(row)
            keep[sire].append(pd.DataFrame(dict(chrom=chrom, pos=d['variants/POS'][m], ref=R[m], alt=A[m, 0], mc_is_ref=mc_is_ref[m])))
    pd.DataFrame(attr).to_csv(proj.dirs['qc'] + 'v2_site_attrition.tsv', sep='\t', index=False)
    L = proj.chrom_lengths(); summ = []
    extra_females = [s for s in cfg['mc_panel'] if s not in panel]
    for sire in sires:
        df = pd.concat(keep[sire])
        if v.get('repeat_mask_bed'):
            df = _mask(df, v['repeat_mask_bed'], proj.dirs['tmp'])
        _write_sites(proj, sire, df)
        n0 = len(df)
        if extra_females:                                           # extra MC females: recount only, they are not in the VCF
            df = panel_recount(proj, sire, df, extra_females)
            _write_sites(proj, sire, df)
        print(f'{sire}: {n0} sites -> {len(df)} after extra-MC-female recount')
        base = f"{proj.dirs['v2']}sites/sites_{sire}"
        for c, g in df.groupby('chrom'):
            summ.append(dict(sire=sire, chrom=c, LG=proj.lg[c], n_sites=len(g), sites_per_mb=len(g) / (L[c] / 1e6),
                             frac_mc_is_ref=float(g.mc_is_ref.mean())))
    s = pd.DataFrame(summ); s.to_csv(proj.dirs['qc'] + 'v2_sites_per_mb.tsv', sep='\t', index=False)
    print(s.groupby('sire').sites_per_mb.describe())


def _write_sites(proj, sire, df):
    base = f"{proj.dirs['v2']}sites/sites_{sire}"
    df.assign(a=df.ref + ',' + df.alt)[['chrom', 'pos', 'a']].to_csv(base + '.tsv', sep='\t', header=False, index=False)
    sh(['bgzip', '-f', base + '.tsv']); sh(['tabix', '-f', '-s1', '-b2', '-e2', base + '.tsv.gz'])
    df.to_csv(base + '.info.tsv.gz', sep='\t', index=False, compression='gzip')


def _panel_one(args):
    proj, s, sire = args
    out = f"{proj.dirs['tmp']}panel.{sire}.{s}.tsv.gz"
    count_bam(proj, proj.ensure_bam(s), sire, out)
    return pd.read_csv(out, sep='\t')


def panel_recount(proj, sire, df, females):
    """Drop sites where the YH allele appears in the extra deep MC females (they should carry only the MC allele).
    Uses the sites files just written for `sire`; counts come from the same filtered bcftools pileup as the BC fish."""
    v = proj.cfg['v2']
    res, errs = parallel(_panel_one, [(proj, s, sire) for s in females], proj.workers('count'), 'panel-recount')
    assert not errs, errs
    allc = pd.concat(res.values())
    g = allc.groupby(['chrom', 'pos']).agg(mc=('mc', 'sum'), yh=('yh', 'sum')).reset_index()
    g['bad'] = (g.yh / (g.mc + g.yh).clip(lower=1) > v['panel_max_alt_frac']) & (g.mc + g.yh >= 6)
    bad = set(zip(g[g.bad].chrom, g[g.bad].pos))
    return df[[(c, p) not in bad for c, p in zip(df.chrom, df.pos)]]


def _mask(df, bed, tmp):
    b = tmp + 'sites.bed'
    df.assign(s=df.pos - 1)[['chrom', 's', 'pos']].to_csv(b, sep='\t', header=False, index=False)
    out = sh(f'bedtools intersect -v -a {b} -b {bed}').stdout.split('\n')
    kept = {(l.split('\t')[0], int(l.split('\t')[2])) for l in out if l}
    return df[[(c, p) in kept for c, p in zip(df.chrom, df.pos)]]


# ------------------------------------------------------------------ 3. allele counts
def count_bam(proj, bam, sire, out):
    """Counts of the MC allele and the YH allele at every informative site of `sire` (deduplicated, MAPQ/BQ filtered)."""
    v = proj.cfg['v2']
    base = f"{proj.dirs['v2']}sites/sites_{sire}"
    tmp = out + '.raw.gz'
    sh(['bash', '-o', 'pipefail', '-c',
        f"bcftools mpileup -f {proj.genome_fa} -T {base}.tsv.gz -q {v['mapq']} -Q {v['baseq']} -I -d 200 "
        f"--ff UNMAP,SECONDARY,QCFAIL,DUP -a FORMAT/AD -Ou {bam} | "
        f"bcftools query -f '%CHROM\\t%POS\\t%ALT\\t[%AD]\\n' | gzip > {tmp}"])
    info = pd.read_csv(base + '.info.tsv.gz', sep='\t')
    raw = pd.read_csv(tmp, sep='\t', header=None, names=['chrom', 'pos', 'alts', 'ad'], dtype=str)
    raw['pos'] = raw.pos.astype(int)
    raw = raw.merge(info, on=['chrom', 'pos'], how='inner')
    ad = raw.ad.str.split(',', expand=True).apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    alts = raw.alts.str.split(',', expand=True)
    n_alt = np.zeros(len(raw), int)
    for k in range(alts.shape[1]):
        m = (alts[k] == raw.alt).values
        if k + 1 < ad.shape[1]:
            n_alt[m] = ad.iloc[:, k + 1].values[m]
    n_ref = ad.iloc[:, 0].values
    mcr = raw.mc_is_ref.astype(bool).values
    res = pd.DataFrame(dict(chrom=raw.chrom, pos=raw.pos, mc=np.where(mcr, n_ref, n_alt), yh=np.where(mcr, n_alt, n_ref), yh_is_ref=~mcr))
    res = res[(res.mc + res.yh) > 0]
    res.to_csv(out, sep='\t', index=False, compression='gzip')
    os.remove(tmp)
    return len(res)


def _count_one(args):
    proj, s, sire = args
    out = f"{proj.dirs['v2']}counts/{s}.tsv.gz"
    if os.path.exists(out):
        return out
    return count_bam(proj, proj.ensure_bam(s), sire, out)


def count_reads(proj):
    os.makedirs(proj.dirs['v2'] + 'counts', exist_ok=True)
    bc = proj.bc_table()
    assert bc.SireSample.notna().all(), 'unmapped sire names; edit f1_sire_names in config'
    res, errs = parallel(_count_one, [(proj, r.SampleID, r.SireSample) for r in bc.itertuples()], proj.workers('count'), 'count-reads')
    pd.Series(errs).to_csv(proj.dirs['qc'] + 'v2_count_errors.tsv', sep='\t')


# ------------------------------------------------------------------ 4. HMM
def hmm_sample(counts, lengths, vcfg, rho=None, step=None):
    """counts: DataFrame chrom,pos,mc,yh,yh_is_ref. Returns {chrom: grid posteriors, breakpoints}, (eps, q)."""
    rho = rho or vcfg['rho']; step = step or vcfg['grid_step']
    per, P, R, A, F = {}, [], [], [], []
    for ch, L in lengths.items():
        c = counts[counts.chrom == ch].sort_values('pos')
        ap, ar, aa, isg, fl = with_grid(c.pos.values, c.mc.values, c.yh.values, L, step, offset=0, flag=c.yh_is_ref.values)
        per[ch] = (ap, ar, aa, isg, fl); P.append(ap); R.append(ar); A.append(aa); F.append(fl)
    eps, q = fit_sample(P, R, A, rho=rho, cap=vcfg['cap'], flags=F)
    out = {}
    for ch, (ap, ar, aa, isg, fl) in per.items():
        post = forward_backward(ap, ar, aa, eps, q, rho, vcfg['cap'], fl)
        _, bps = call_breakpoints(ap, post)
        out[ch] = dict(pos=ap[isg], post=post[isg], bps=bps, n_obs=int(((ar + aa)[~isg] > 0).sum()))
    return out, (eps, q)


def _hmm_one(args):
    proj, s, rho, lengths = args
    c = pd.read_csv(f"{proj.dirs['v2']}counts/{s}.tsv.gz", sep='\t')
    return hmm_sample(c, lengths, proj.cfg['v2'], rho)


def run_hmm(proj, rho=None, tag=''):
    lengths = proj.chrom_lengths()
    bc = proj.bc_table()
    samples = [s for s in bc.SampleID if os.path.exists(f"{proj.dirs['v2']}counts/{s}.tsv.gz")]
    res, errs = parallel(_hmm_one, [(proj, s, rho, lengths) for s in samples], proj.workers('hmm'), 'v2-hmm')
    done = {samples[i]: r for i, r in res.items()}
    bp_rows, qc_rows = [], []
    for s, (out, (eps, q)) in done.items():
        qc_rows.append(dict(sample=s, eps=eps, q_ref_in_het=q, n_obs_sites=sum(o['n_obs'] for o in out.values()),
                            n_switch=sum(len(o['bps']) for o in out.values())))
        for ch, o in out.items():
            for lo, hi, cl, cr in o['bps']:
                bp_rows.append(dict(sample=s, chrom=ch, bp_lo=lo, bp_hi=hi, conf_left=cl, conf_right=cr))
    first = next(iter(done.values()))[0]
    grid = pd.concat([pd.DataFrame({'Chromosome': ch, 'Position': o['pos']}) for ch, o in first.items()], ignore_index=True)
    wide = {s: np.concatenate([done[s][0][ch]['post'] for ch in first]) for s in done}
    out = pd.concat([grid, pd.DataFrame(wide)], axis=1)
    out.to_csv(f"{proj.dirs['v2']}OutputPosteriors{tag}.csv", index=False)
    pd.DataFrame(bp_rows).to_csv(f"{proj.dirs['v2']}breakpoints{tag}.tsv", sep='\t', index=False)
    qc = pd.DataFrame(qc_rows)
    qc['flag'] = (qc.eps > 0.03) | (qc.n_switch > 12)
    qc.to_csv(f"{proj.dirs['v2']}sample_qc{tag}.tsv", sep='\t', index=False)
    print('v2 HMM done:', out.shape, 'flagged samples:', int(qc.flag.sum()), 'errors:', len(errs))
