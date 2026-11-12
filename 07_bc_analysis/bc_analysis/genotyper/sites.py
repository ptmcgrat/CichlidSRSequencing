"""Informative sites of the BC cross, defined from the parents' cohort VCF (module 05, cohort `bc_family`, called with bcftools).

The rule (the original script's, written with the alleles named and the depth made relative):
  a site is informative for a BROOD when the YH grandsire is homozygous for one allele, and the F1 sires' MC grandmother (`f1_dam`) and THE BROOD'S OWN MC DAM are both homozygous for the other allele.
  Then every F1 sire of the cross is heterozygous (YH allele from his sire, MC allele from his dam) and the BC fish inherited the MC allele from its dam, so its genotype says which allele it got from the F1 sire.
  A heterozygous site in an MC mother who is not the dam of this brood does not matter, so the sites are defined per dam; the three dams give three site sets.
  Both polarities are kept: the YH allele can be the reference allele (yh_is_ref) or the alternative allele, because the reference (Mzebra) is neither parent.
The depth thresholds are relative to each parent's own median depth: the old filter (more than 12 reads) was written for another alignment; the parents here have 9 to 13x effective depth, and a collapsed repeat
shows as too MUCH depth, so there is an upper bound as well. A homozygous call may carry only a few reads of the other allele.
"""
import gzip
import hashlib
import json
import os
import statistics
from collections import Counter

import pandas as pd
import yaml

CONFIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config', 'genotyper.yaml')


def load_config(path=None):
    return yaml.safe_load(open(path or CONFIG))


# ------------------------------------------------------------------ reading the VCF (no dependency: plain text, gzip or bgzip)
def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_header(path):
    with _open(path) as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                return line.rstrip('\n').split('\t')[9:]
    raise ValueError(f'{path} has no #CHROM line')


def _genotype(field, keys):
    """(gt 0/1/2 = hom ref/het/hom alt, -1 missing; ad_ref, ad_alt, dp) of one sample column."""
    parts = field.split(':')
    d = dict(zip(keys, parts))
    g = d.get('GT', './.').replace('|', '/')
    gt = {'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2}.get(g, -1)
    ad_ref = ad_alt = 0
    if 'AD' in d and d['AD'] not in ('.', ''):
        a = d['AD'].split(',')
        try:
            ad_ref, ad_alt = int(a[0]), int(a[1])
        except (ValueError, IndexError):
            pass
    try:
        dp = int(d['DP']) if d.get('DP', '.') not in ('.', '') else ad_ref + ad_alt
    except ValueError:
        dp = ad_ref + ad_alt
    return gt, ad_ref, ad_alt, dp


def iter_snps(path, samples, limit=None):
    """Biallelic SNP records as (chrom, pos, ref, alt, {sample: (gt, ad_ref, ad_alt, dp)}); only the named samples are parsed."""
    names = read_header(path)
    idx = {s: names.index(s) for s in samples}               # a ValueError names the missing sample
    n = 0
    with _open(path) as fh:
        for line in fh:
            if line[0] == '#':
                continue
            f = line.rstrip('\n').split('\t')
            ref, alt = f[3], f[4]
            if len(ref) != 1 or len(alt) != 1 or alt == '.':
                continue
            keys = f[8].split(':')
            if 'GT' not in keys:
                continue
            yield f[0], int(f[1]), ref, alt, {s: _genotype(f[9 + i], keys) for s, i in idx.items()}
            n += 1
            if limit and n >= limit:
                return


def median_depths(path, samples, n_sites):
    dps = {s: [] for s in samples}
    for _, _, _, _, g in iter_snps(path, samples, limit=n_sites):
        for s, (gt, _, _, dp) in g.items():
            if gt >= 0 and dp > 0:
                dps[s].append(dp)
    return {s: (statistics.median(v) if v else 0.0) for s, v in dps.items()}


# ------------------------------------------------------------------ the rule
def homozygous(call, med, sc):
    """(allele 0/2, ok): a homozygous call with the depth and the purity the parameters ask for."""
    gt, ad_ref, ad_alt, dp = call
    if gt not in (0, 2):
        return None
    lo = max(sc['min_dp'], sc['dp_frac_lo'] * med)
    hi = sc['dp_frac_hi'] * med
    if not (lo <= dp <= hi):
        return None
    minor = ad_alt if gt == 0 else ad_ref
    if minor > sc['hom_max_minor_reads'] or minor > sc['hom_max_minor_frac'] * dp:
        return None
    return gt


STEPS = ['snp_biallelic', 'yh_homozygous', 'f1_dam_opposite', 'dam_opposite']


def define_sites(vcf, cfg, log=print):
    """Returns dict(sites={dam: [(chrom, pos, ref, alt, yh_is_ref)]}, attrition=[rows], medians={parent: depth}, key=str). One pass over the VCF after a pass over its first `depth_sample_sites` SNPs."""
    sc = cfg['sites']
    yh, f1, dams = sc['yh_sire'], sc['f1_dam'], list(sc['dams'])
    parents = list(dict.fromkeys([yh, f1] + dams))
    med = median_depths(vcf, parents, sc['depth_sample_sites'])
    for p in parents:
        if med[p] <= 0:
            raise ValueError(f'no depth could be read for {p} in the first {sc["depth_sample_sites"]} SNPs of {vcf}')
    log(f'median depth of the parents: ' + ', '.join(f'{p} {med[p]:.0f}' for p in parents))
    sites = {d: [] for d in dams}
    counts = {d: Counter() for d in dams}
    per_chrom = {d: Counter() for d in dams}
    for chrom, pos, ref, alt, g in iter_snps(vcf, parents):
        for d in dams:
            counts[d]['snp_biallelic'] += 1
        ay = homozygous(g[yh], med[yh], sc)
        if ay is None:
            continue
        for d in dams:
            counts[d]['yh_homozygous'] += 1
        af = homozygous(g[f1], med[f1], sc)
        if af is None or af == ay:
            continue
        for d in dams:
            counts[d]['f1_dam_opposite'] += 1
            ad = homozygous(g[d], med[d], sc)
            if ad is None or ad != af:
                continue
            counts[d]['dam_opposite'] += 1
            per_chrom[d][chrom] += 1
            sites[d].append((chrom, pos, ref, alt, ay == 0))
    attrition = [dict(dam=d, step=s, sites=counts[d][s]) for d in dams for s in STEPS]
    key = definition_key(cfg, vcf)
    return dict(sites=sites, attrition=attrition, per_chrom=per_chrom, medians=med, key=key)


def definition_key(cfg, vcf):
    blob = json.dumps(dict(sites=cfg['sites'], vcf=os.path.basename(vcf), size=os.path.getsize(vcf)), sort_keys=True)
    return hashlib.sha1(blob.encode()).hexdigest()[:12]


# ------------------------------------------------------------------ outputs
def write_sites(result, out_dir, vcf, cfg, fai=None):
    """sites_<dam>.tsv.gz (chrom, pos, ref, alt, yh_is_ref), targets_<dam>.tsv (chrom, pos: for bcftools mpileup -T), attrition.tsv, sites_per_mb.tsv, definition.json."""
    os.makedirs(out_dir, exist_ok=True)
    for d, rows in result['sites'].items():
        df = pd.DataFrame(rows, columns=['chrom', 'pos', 'ref', 'alt', 'yh_is_ref'])
        df['yh_is_ref'] = df.yh_is_ref.astype(int)
        df.to_csv(os.path.join(out_dir, f'sites_{d}.tsv.gz'), sep='\t', index=False)
        df[['chrom', 'pos']].to_csv(os.path.join(out_dir, f'targets_{d}.tsv'), sep='\t', index=False, header=False)
    pd.DataFrame(result['attrition']).to_csv(os.path.join(out_dir, 'attrition.tsv'), sep='\t', index=False)
    lengths = dict(fai) if fai else {}
    rows = []
    for d, c in result['per_chrom'].items():
        for ch, n in c.items():
            rows.append(dict(dam=d, chrom=ch, sites=n, sites_per_mb=(n / (lengths[ch] / 1e6)) if ch in lengths else None))
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, 'sites_per_mb.tsv'), sep='\t', index=False)
    json.dump(dict(key=result['key'], vcf=os.path.basename(vcf), vcf_bytes=os.path.getsize(vcf), parameters=cfg['sites'], median_depth=result['medians'],
                   sites={d: len(r) for d, r in result['sites'].items()}), open(os.path.join(out_dir, 'definition.json'), 'w'), indent=1)
