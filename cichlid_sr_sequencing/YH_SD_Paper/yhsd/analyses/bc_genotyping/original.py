"""ORIGINAL approach, regenerated from scratch for the chosen genome (clean VCFs), then run.
Algorithm unchanged: markers = YH sire homozygous for one allele & F1 dam and brood dam homozygous for the other (depth > 12);
every read is one observation; 2-state CategoricalHMM, fixed emissions, trans_prob=1e-22, Viterbi; interpolation to a 100 kb grid.
Deliberate changes: (1) the old filter used GT[4] (a BC dam) where the YH sire was meant - fixed by looking samples up by name;
(2) ALLELE POLARIZATION (config original.polarity): with a non-MC reference, MC carries the ALT allele at ~half the sites, so
'both' flips the observation coding there ('alt_only' = the exact old site rule); (3) seeded RNG; (4) parallel; (5) extra bookkeeping."""
import os
from functools import reduce
import numpy as np
import pandas as pd
from ...common import sh, parallel
from .v2 import joint_call

_S = {}   # read-only state inherited by forked workers


# ------------------------------------------------------------------ A. regenerate the VCFs
def parent_samples(cfg):
    return list(dict.fromkeys([cfg['yh_sire'], cfg['f1_dam']] + cfg['bc_dams']))


def make_parents_vcfs(proj):
    cfg, oc = proj.cfg, proj.cfg['original']
    d = proj.dirs['original']
    unf = joint_call(proj, parent_samples(cfg), d, 'AllBroodParents_Unfiltered')
    filt = d + 'AllBroodParents_Filtered.vcf.gz'
    names = sh(['bcftools', 'query', '-l', unf]).stdout.split()
    iy, idm = names.index(cfg['yh_sire']), names.index(cfg['f1_dam'])      # look up by NAME (old code used GT[4])
    expr = f'GT[{iy}]="AA" && GT[{idm}]="RR"'
    if oc['polarity'] == 'both':
        expr = f'(({expr}) || (GT[{iy}]="RR" && GT[{idm}]="AA"))'
    expr += f' && FMT/DP[{iy}]>{oc["min_dp"]} && FMT/DP[{idm}]>{oc["min_dp"]}'
    sh(['bcftools', 'view', '-m2', '-M2', '-v', 'snps', '-i', expr, '-Oz', '-o', filt, unf])
    sh(['tabix', '-f', '-p', 'vcf', filt])
    open(d + 'filter_expression.txt', 'w').write(expr + '\n')
    n = sh(f'bcftools view -H {filt} | wc -l').stdout.strip()
    print('filtered parent sites:', n)


def _bam_one(args):
    proj, s = args
    return proj.ensure_bam(s)


def _pileup_chrom(args):
    proj, chrom, listfile = args
    d = proj.dirs['original']; q = proj.cfg['original']['mpileup_q']
    pile, var = f'{d}per_chrom/BC_Pileup.{chrom}.vcf.gz', f'{d}per_chrom/BC_Variants.{chrom}.vcf.gz'
    if os.path.exists(var + '.tbi'):
        return chrom
    sh(['bash', '-o', 'pipefail', '-c',
        f"ulimit -n 8192 && bcftools mpileup -f {proj.genome_fa} -r {chrom} -R {d}AllBroodParents_Filtered.vcf.gz -b {listfile} "
        f"-a AD -Q {q} -Oz -o {pile} --threads 2 && tabix -f -p vcf {pile} && "
        f"bcftools call -c -Oz -o {var} {pile} && tabix -f -p vcf {var}"])
    return chrom


def make_pileup_vcfs(proj):
    cfg = proj.cfg; d = proj.dirs['original']
    os.makedirs(d + 'per_chrom/', exist_ok=True)
    bc = proj.bc_table().SampleID.tolist()
    samples = list(dict.fromkeys(parent_samples(cfg) + list(cfg['f1_sire_names'].values()) + bc))
    res, errs = parallel(_bam_one, [(proj, s) for s in samples], proj.workers('flagstat'), 'bams')
    assert not errs, 'BAM download/index failed (see above)'
    listfile = d + 'BC_Pileup_Samples.txt'
    open(listfile, 'w').write('\n'.join(proj.bam(s) for s in samples) + '\n')
    r, e = parallel(_pileup_chrom, [(proj, c, listfile) for c in proj.chroms], proj.workers('pileup'), 'pileup')
    assert not e, e
    for kind, out in (('BC_Pileup', 'BC_Pileup.vcf.gz'), ('BC_Variants', 'BC_Variants.vcf.gz')):
        sh(['bcftools', 'concat', '-Oz', '-o', d + out] + [f'{d}per_chrom/{kind}.{c}.vcf.gz' for c in proj.chroms])
        sh(['tabix', '-f', '-p', 'vcf', d + out])
    print('wrote', d + 'BC_Pileup.vcf.gz', 'and BC_Variants.vcf.gz')


def make_vcfs(proj):
    make_parents_vcfs(proj)
    make_pileup_vcfs(proj)


# ------------------------------------------------------------------ B. the HMM
def load_vcf(path):
    import allel
    v = allel.read_vcf(path, fields=['variants/CHROM', 'variants/POS', 'calldata/GT', 'calldata/AD', 'samples'])
    return dict(samples=np.array(v['samples']), chrom=v['variants/CHROM'], pos=v['variants/POS'],
                ad=v['calldata/AD'][:, :, 0:2], gt=allel.GenotypeArray(v['calldata/GT']))


def align_rows(a, b):
    """Return row indices of a and b for the (chrom,pos) keys they share, in a's order."""
    ka = pd.MultiIndex.from_arrays([a['chrom'], a['pos']]); kb = pd.MultiIndex.from_arrays([b['chrom'], b['pos']])
    common = ka.intersection(kb)
    ia = ka.get_indexer(common); ib = kb.get_indexer(common)
    o = np.argsort(ia)
    return ia[o], ib[o]


def brood_masks(var, cfg):
    """dam -> (marker mask, mc_is_alt). Uses genotype calls of the YH sire, F1 dam and the brood dam from BC_Variants."""
    idx = lambda s: int((var['samples'] == s).argmax())
    g = var['gt']
    yh, fd = g[:, idx(cfg['yh_sire'])], g[:, idx(cfg['f1_dam'])]
    out = {}
    for m in cfg['bc_dams']:
        md = g[:, idx(m)]
        mask = yh.is_hom_alt() & fd.is_hom_ref() & md.is_hom_ref()
        if cfg['original']['polarity'] == 'both':
            mask = mask | (yh.is_hom_ref() & fd.is_hom_alt() & md.is_hom_alt())
        out[m] = (mask, np.asarray(fd.is_hom_alt()))
    return out


def stretch(ad, chrom, pos, rng):
    """One observation per read (0=MC allele, 1=YH allele), order shuffled within each site (as in the original)."""
    ad = np.clip(ad, 0, None).astype(np.int64)
    n = len(ad)
    obs = np.repeat(np.tile([0, 1], n), ad.ravel()).astype(np.uint8)
    site = np.repeat(np.repeat(np.arange(n), 2), ad.ravel())
    order = np.lexsort((rng.random(len(obs)), site))
    return obs[order], np.repeat(chrom, ad.sum(1)), np.repeat(pos, ad.sum(1))


def viterbi_numpy(obs, emission, trans_prob):
    """Same model as hmmlearn.CategoricalHMM(2 states) with start=.5; Viterbi path."""
    le = np.log(np.array(emission)); pen = -np.log(trans_prob); stay = np.log1p(-trans_prob)
    n = len(obs)
    v0, v1 = np.log(.5) + le[0][obs[0]], np.log(.5) + le[1][obs[0]]
    b0 = np.zeros(n, bool); b1 = np.zeros(n, bool)
    for i in range(1, n):
        o = obs[i]
        s0 = v0 + stay >= v1 - pen; s1 = v1 + stay >= v0 - pen
        n0 = (v0 + stay if s0 else v1 - pen) + le[0][o]
        n1 = (v1 + stay if s1 else v0 - pen) + le[1][o]
        b0[i], b1[i] = not s0, not s1
        v0, v1 = n0, n1
    st = 0 if v0 >= v1 else 1
    path = np.empty(n, np.uint8)
    for i in range(n - 1, -1, -1):
        path[i] = st
        if i and (b0[i] if st == 0 else b1[i]):
            st = 1 - st
    return path


def hmm_path(obs, oc):
    if oc['engine'] == 'numpy':
        return viterbi_numpy(obs, oc['emission'], oc['trans_prob'])
    from hmmlearn.hmm import CategoricalHMM
    tp = oc['trans_prob']
    m = CategoricalHMM(n_components=2, n_features=2)
    m.startprob_ = np.array([.5, .5]); m.transmat_ = np.array([[1 - tp, tp], [tp, 1 - tp]])
    m.emissionprob_ = np.array(oc['emission'])
    return m.predict(obs.reshape(-1, 1))


def _fish(i):
    s = _S
    brood_mask, mc_is_alt, offspring = s['jobs'][i]
    oc = s['cfg']['original']
    rng = np.random.default_rng(oc['seed'] + i)
    j = int((s['pile']['samples'] == offspring).argmax())
    ad = s['pile']['ad'][brood_mask, j]
    ad = np.where(mc_is_alt[brood_mask][:, None], ad[:, ::-1], ad)        # col0 = MC allele, col1 = YH allele
    chrom = s['pile']['chrom'][brood_mask]; pos = s['pile']['pos'][brood_mask]
    obs, ochrom, opos = stretch(ad, chrom, pos, rng)
    frames, stats = [], []
    for ch in dict.fromkeys(chrom):
        if ch not in s['chroms']:
            continue
        m = ochrom == ch
        if m.sum() < 2:
            continue
        path = hmm_path(obs[m], oc)
        p = opos[m]
        grid = np.arange(0, p.max(), oc['grid_step'])
        frames.append(pd.DataFrame({'Chromosome': ch, 'Position': grid, offspring: np.interp(grid, p, path)}))
        stats.append((offspring, ch, int(m.sum()), int((path[1:] != path[:-1]).sum())))
    return pd.concat(frames), stats


def main(proj):
    cfg = proj.cfg; oc = cfg['original']
    d = proj.dirs['original']
    pile = load_vcf(d + 'BC_Pileup.vcf.gz'); var = load_vcf(d + 'BC_Variants.vcf.gz')
    ia, ib = align_rows(pile, var)
    if len(ia) != len(pile['pos']) or len(ib) != len(var['pos']):
        print(f'NOTE: pileup/variant VCF rows differ ({len(pile["pos"])} vs {len(var["pos"])}); using {len(ia)} shared sites')
    pile = {k: (v[ia] if k != 'samples' else v) for k, v in pile.items()}
    var = {k: (v[ib] if k != 'samples' else v) for k, v in var.items()}
    masks = brood_masks(var, cfg)
    bc = proj.bc_table()
    jobs, rows = [], []
    for brood in cfg['broods']:
        b = bc[bc.BroodID == brood]
        if len(b) == 0:
            continue
        assert b.Dam.nunique() == 1, f'{brood} has several dams: {b.Dam.unique()}'
        mk, mc_is_alt = masks[b.Dam.iloc[0]]
        for ch in dict.fromkeys(pile['chrom'][mk]):
            rows.append((brood, b.Dam.iloc[0], ch, proj.lg.get(ch, ch), int(((pile['chrom'] == ch) & mk).sum()), float(mc_is_alt[mk & (pile['chrom'] == ch)].mean())))
        jobs += [(mk, mc_is_alt, s) for s in b.SampleID if s in pile['samples']]
    pd.DataFrame(rows, columns=['brood', 'dam', 'chrom', 'LG', 'n_markers', 'frac_mc_is_alt']).to_csv(d + 'markers_per_brood_chrom.tsv', sep='\t', index=False)
    _S.update(pile=pile, jobs=jobs, cfg=cfg, chroms=set(proj.chroms))
    res, errs = parallel(_fish, list(range(len(jobs))), proj.workers('hmm'), 'original-hmm')
    frames = [res[i][0] for i in sorted(res)]
    merged = reduce(lambda l, r: pd.merge(l, r, on=['Chromosome', 'Position'], how='inner'), frames)
    merged.to_csv(d + 'OutputGenotypes_original.csv')
    pd.DataFrame([x for i in sorted(res) for x in res[i][1]], columns=['sample', 'chrom', 'n_reads', 'n_switch']).to_csv(d + 'per_fish_chrom_reads_switches.tsv', sep='\t', index=False)
    print('original approach done:', merged.shape, 'failed fish:', len(errs))
