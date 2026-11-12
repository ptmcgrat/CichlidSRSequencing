"""Allele counts of one BC fish at the informative sites of ITS dam, from its BAM, stored so that the HMM can be rerun without the BAM.

`bcftools mpileup -T <targets>` over the whole BAM with a mapping-quality and a base-quality filter, duplicates and secondary alignments excluded, indels skipped (-I), AD per allele. The counts are read by ALLELE NAME,
never by position in the list: the alleles of a pileup line are REF followed by the observed alternative alleles (and <*>), in an order that depends on the data. The YH allele of a site is the one the parents'
VCF says (the reference allele when yh_is_ref, else the alternative allele of the site); the MC allele is the other one. Reads of any third allele are not counted.

Output per fish: <out>/counts/<sample>.counts.tsv.gz (chrom, pos, n_mc, n_yh, yh_is_ref; sites with reads only) and <sample>.counts.json (the definition key of the sites, the BAM's name and size, the totals).
A fish is measured again only when the key or the BAM changed (name and size, as the SD panel does)."""
import datetime
import gzip
import json
import os
import shutil
import socket
import time

import pandas as pd

from variant_calling import build as B


RUNNERS = dict(run=B.run_cmd, script=B.run_script)                 # what the command line runs tools with (the tests replace them)


def mpileup_script(fasta, bam, targets, ccfg):
    return (f"bcftools mpileup -f {fasta} -T {targets} -q {ccfg['mapq']} -Q {ccfg['baseq']} -I -d {ccfg['max_depth']} --ff {ccfg['exclude_flags']} -a FORMAT/AD -Ou {bam} | "
            f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AD]\\n'")


def load_sites(path):
    """{(chrom, pos): (ref, alt, yh_is_ref)}"""
    df = pd.read_csv(path, sep='\t', dtype={'chrom': str})
    return {(c, p): (r, a, bool(y)) for c, p, r, a, y in zip(df.chrom, df.pos, df.ref, df.alt, df.yh_is_ref)}


def parse_pileup(text, sites):
    """Rows (chrom, pos, n_mc, n_yh, yh_is_ref) of the lines of a pileup query that are informative sites and have at least one MC or YH read."""
    rows = []
    for line in text.splitlines():
        f = line.split('\t')
        if len(f) < 5:
            continue
        key = (f[0], int(f[1]))
        info = sites.get(key)
        if info is None:
            continue
        ref, alt, yh_is_ref = info
        if f[2] != ref:
            continue                                                  # the reference base differs from the VCF's: a different reference, not our site
        alleles = [f[2]] + ([] if f[3] in ('.', '') else f[3].split(','))
        try:
            ad = [int(x) for x in f[4].split(',')]
        except ValueError:
            continue
        count = {a: (ad[i] if i < len(ad) else 0) for i, a in enumerate(alleles)}
        yh_allele, mc_allele = (ref, alt) if yh_is_ref else (alt, ref)
        n_yh, n_mc = count.get(yh_allele, 0), count.get(mc_allele, 0)
        if n_yh + n_mc > 0:
            rows.append((key[0], key[1], n_mc, n_yh, int(yh_is_ref)))
    return rows


def counts_paths(out_dir, sample):
    d = os.path.join(out_dir, 'counts')
    return os.path.join(d, f'{sample}.counts.tsv.gz'), os.path.join(d, f'{sample}.counts.json')


def current(out_dir, sample, key, bam_bytes=None):
    """True when the counts of this fish exist for this definition key (and for this BAM, when its size is known)."""
    tsv, js = counts_paths(out_dir, sample)
    if not (os.path.exists(tsv) and os.path.exists(js)):
        return False
    try:
        meta = json.load(open(js))
    except ValueError:
        return False
    return meta.get('key') == key and meta.get('status') == 'ok' and (bam_bytes is None or meta.get('bam_bytes') == int(bam_bytes))


def count_sample(env, cfg, sample, dam, sites_path, targets_path, key, out_dir, staged, want_bytes=None, keep_bam=False, run=B.run_cmd, script=B.run_script):
    """Download the BAM, count, delete the BAM. A failure is returned (status 'failed' and the reason), never raised."""
    t0 = time.time()
    tsv, js = counts_paths(out_dir, sample)
    os.makedirs(os.path.dirname(tsv), exist_ok=True)
    meta = dict(sample=sample, dam=dam, key=key, machine=socket.gethostname(), measured=datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), bam=f'{sample}.all.bam')
    d = f'{staged}{sample}/'
    try:
        B.stage_bam(env, staged, sample, want_bytes, run)
        bam = f'{d}{sample}.all.bam'
        meta['bam_bytes'] = os.path.getsize(bam)
        out = B.must(script(mpileup_script(env.genome_fa, bam, targets_path, cfg['counts'])), f'the pileup of {sample} at the sites of {dam}')
        rows = parse_pileup(out, load_sites(sites_path))
        df = pd.DataFrame(rows, columns=['chrom', 'pos', 'n_mc', 'n_yh', 'yh_is_ref'])
        df.to_csv(tsv, sep='\t', index=False, compression='gzip')
        meta.update(status='ok', error='', n_sites=len(df), n_reads=int(df.n_mc.sum() + df.n_yh.sum()), yh_fraction=float(df.n_yh.sum() / max(1, df.n_mc.sum() + df.n_yh.sum())))
    except Exception as e:                                          # noqa: one fish must never stop the others
        meta.update(status='failed', error=f'{type(e).__name__}: {e}'[:400])
    finally:
        if not keep_bam:
            shutil.rmtree(d, ignore_errors=True)
    meta['seconds'] = round(time.time() - t0, 1)
    json.dump(meta, open(js, 'w'), indent=1)
    return meta


def read_counts(out_dir, sample):
    """{chrom: (pos, n_mc, n_yh, yh_is_ref)} of one fish (arrays), from its counts file."""
    import numpy as np
    tsv, _ = counts_paths(out_dir, sample)
    df = pd.read_csv(tsv, sep='\t', dtype={'chrom': str})
    out = {}
    for ch, g in df.groupby('chrom', sort=False):
        g = g.sort_values('pos')
        out[ch] = (g.pos.to_numpy(), g.n_mc.to_numpy(), g.n_yh.to_numpy(), g.yh_is_ref.to_numpy().astype(bool))
    return out
