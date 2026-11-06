"""Generic callset builder. A callset = named sample set + GATK joint calling + documented filters -> one VCF.
Defined in config/callsets.yaml; built with `python -m yhsd callset build <name>` (or from an analysis step).
Calling follows Kumar et al. 2025 (HaplotypeCaller GVCFs -> GenomicsDBImport -> GenotypeGVCFs, --heterozygosity 0.00175) with
hard filters from Kumar / the GBE sex-chromosome paper (QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum), genotype filters
(DP, GQ), biallelic SNPs only and no SNPs within 10 bp of an indel. A bcftools allele-recount cross-check can drop sites where
the two callers disagree (`crosscheck`)."""
import os
import numpy as np
import pandas as pd
import yaml
from cichlid_common.common_light import sh, config_path
from cichlid_common.common import parallel


def load_callsets(path=None):
    return yaml.safe_load(open(path or config_path('callsets.yaml')))


# ------------------------------------------------------------------ sample selection
def resolve_samples(dt, spec):
    """dt: SampleDatabase DataFrame (SampleID, Species, BroodID, Parents, Category). Returns ordered unique sample list."""
    out = list(spec.get('samples', []))
    if spec.get('broods'):
        b = dt[dt.BroodID.isin(spec['broods'])]
        out += b.SampleID.tolist()
        if spec.get('include_parents', True):
            out += sorted({p.strip() for ps in b.Parents.dropna() for p in ps.split(',')})
    if spec.get('categories'):
        m = dt.Category.isin(spec['categories'])
        if spec.get('species'):
            m &= dt.Species == spec['species']
        out += dt[m].SampleID.tolist()
    seen, res = set(), []
    for s in out:
        if s not in seen:
            seen.add(s); res.append(s)
    return res


# ------------------------------------------------------------------ calling
def _gdb_genotype(args):
    proj, chrom, mapfile, outdir, name, g = args
    out = f'{outdir}per_chrom/{name}.{chrom}.raw.vcf.gz'
    if os.path.exists(out + '.tbi'):
        return out
    mem = proj.cfg.get('gatk_mem_gb', 8)
    gdb = f"{proj.dirs['tmp']}gdb.{name}.{chrom}"
    if os.path.exists(gdb):
        sh(['rm', '-rf', gdb])
    sh(['gatk', '--java-options', f'-Xmx{mem}g', 'GenomicsDBImport', '--genomicsdb-workspace-path', gdb, '--intervals', chrom,
        '--sample-name-map', mapfile, '--tmp-dir', proj.dirs['tmp'], '--reader-threads', '4'])
    ann = [x for a in g['annotations'] for x in ('-A', a)] + [x for a in g['annotation_groups'] for x in ('-G', a)]
    sh(['gatk', '--java-options', f'-Xmx{mem}g', 'GenotypeGVCFs', '-R', proj.genome_fa, '-V', f'gendb://{gdb}', '-O', out,
        '--heterozygosity', str(g['heterozygosity'])] + ann)
    sh(['rm', '-rf', gdb])
    return out


def _filter_chrom(args):
    proj, chrom, outdir, name, f = args
    raw = f'{outdir}per_chrom/{name}.{chrom}.raw.vcf.gz'; out = f'{outdir}per_chrom/{name}.{chrom}.vcf.gz'
    if os.path.exists(out + '.tbi'):
        return out
    site = ' || '.join([f"INFO/QD<{f['qd']}", f"INFO/FS>{f['fs']}", f"INFO/SOR>{f['sor']}", f"INFO/MQ<{f['mq']}",
                        f"INFO/MQRankSum<{f['mqranksum']}", f"INFO/ReadPosRankSum<{f['readposranksum']}"])
    gt = f"FMT/DP<{f['min_dp']} | FMT/GQ<{f['min_gq']}"
    cmd = (f"bcftools filter -g {f['snpgap']} -e '{site}' -Ou {raw} | "                      # SnpGap needs indels still present
           f"bcftools view -f PASS,. -m2 -M2 -v snps -Ou | "
           f"bcftools filter -S . -e '{gt}' -Oz -o {out} && tabix -f -p vcf {out}")
    sh(['bash', '-o', 'pipefail', '-c', cmd])
    return out


def build(proj, name, spec=None, dry_run=False):
    cs = load_callsets()
    spec = spec or cs['callsets'][name]
    g = {**cs['defaults']['gatk'], **spec.get('gatk', {})}
    f = {**cs['defaults']['filters'], **spec.get('filters', {})}
    dt = proj.sample_table()
    samples = resolve_samples(dt, spec)
    if dry_run:
        print(name, len(samples), 'samples:', samples); return samples
    with proj.run('callsets', name, params=dict(spec=spec, gatk=g, filters=f, samples=samples), tar=()) as run:
        outdir = run.dir
        os.makedirs(outdir + 'per_chrom/', exist_ok=True)
        gvcfs = [proj.ensure_gvcf(s) for s in samples]                  # cloud GVCFs, or HaplotypeCaller from the BAM
        mapfile = outdir + f'{name}.sample_map.txt'
        open(mapfile, 'w').write(''.join(f'{s}\t{p}\n' for s, p in zip(samples, gvcfs)))
        run.add_input(mapfile)
        r, e = parallel(_gdb_genotype, [(proj, c, mapfile, outdir, name, g) for c in proj.chroms], proj.workers('joint'), 'genotype')
        assert not e, e
        r, e = parallel(_filter_chrom, [(proj, c, outdir, name, f) for c in proj.chroms], proj.workers('joint'), 'filter')
        assert not e, e
        files = [f'{outdir}per_chrom/{name}.{c}.vcf.gz' for c in proj.chroms]
        vcf = f'{outdir}{name}.vcf.gz'
        sh(['bcftools', 'concat', '-Oz', '-o', vcf] + files); sh(['tabix', '-f', '-p', 'vcf', vcf])
        sh(f'bcftools stats -s - {vcf} > {outdir}{name}.stats.txt')
        sh(f'rm -rf {outdir}per_chrom')                                  # per-chromosome intermediates are not kept
        open(outdir + f'{name}.samples.txt', 'w').write('\n'.join(samples) + '\n')
    return run.dir + f'{name}.vcf.gz'


def vcf_path(proj, name):
    """VCF of the most recent run of callset `name` (pulled from Dropbox if not on local disk)."""
    from .runner import latest_dir
    return latest_dir(proj, 'callsets', name) + f'{name}.vcf.gz'


# ------------------------------------------------------------------ bcftools cross-check
def classify_counts(ref, alt, min_dp=8, het_lo=0.25, het_hi=0.75, hom_max=0.1):
    """Genotype class from allele counts: 0 hom-ref, 1 het, 2 hom-alt, -1 ambiguous/low depth."""
    ref, alt = np.asarray(ref, float), np.asarray(alt, float)
    dp = ref + alt; f = alt / np.maximum(dp, 1)
    out = np.full(ref.shape, -1)
    ok = dp >= min_dp
    out[ok & (f <= hom_max)] = 0
    out[ok & (f >= het_lo) & (f <= het_hi)] = 1
    out[ok & (f >= 1 - hom_max)] = 2
    return out


def crosscheck(proj, vcf, check_samples, out_vcf, min_dp=8, max_disagree=0):
    """Recount alleles with bcftools mpileup in `check_samples` (deep samples) at every VCF site and drop sites where the genotype
    class disagrees with GATK in more than `max_disagree` samples. Writes out_vcf (+ list of dropped sites)."""
    v = proj.cfg.get('v2', {})
    sites = out_vcf + '.sites.tsv'
    sh(f"bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' {vcf} > {sites} && bgzip -f {sites} && tabix -f -s1 -b2 -e2 {sites}.gz")
    gatk = {}
    for s in check_samples:
        raw = sh(f"bcftools query -s {s} -f '%CHROM\\t%POS\\t[%GT]\\n' {vcf}").stdout.strip().split('\n')
        gatk[s] = {(l.split('\t')[0], int(l.split('\t')[1])): l.split('\t')[2] for l in raw if l}
    bad = {}
    for s in check_samples:
        bam = proj.ensure_bam(s)
        txt = sh(['bash', '-o', 'pipefail', '-c',
                  f"bcftools mpileup -f {proj.genome_fa} -T {sites}.gz -q {v.get('mapq', 20)} -Q {v.get('baseq', 20)} -I -B -d 400 "
                  f"--ff UNMAP,SECONDARY,QCFAIL,DUP -a FORMAT/AD -Ou {bam} | bcftools query -f '%CHROM\\t%POS\\t%ALT\\t[%AD]\\n'"]).stdout
        for line in txt.strip().split('\n'):
            c, p, alts, ad = line.split('\t'); ad = [int(x) for x in ad.split(',')]
            key = (c, int(p)); g = gatk[s].get(key)
            if g is None:
                continue
            cls = classify_counts([ad[0]], [ad[1] if len(ad) > 1 else 0], min_dp)[0]
            gcls = {'0/0': 0, '0|0': 0, '0/1': 1, '0|1': 1, '1|0': 1, '1/1': 2, '1|1': 2}.get(g, -1)
            if cls != -1 and gcls != -1 and cls != gcls:
                bad[key] = bad.get(key, 0) + 1
    drop = sorted(k for k, n in bad.items() if n > max_disagree)
    pd.DataFrame(drop, columns=['chrom', 'pos']).to_csv(out_vcf + '.dropped.tsv', sep='\t', index=False, header=False)
    sh(f"bgzip -f {out_vcf}.dropped.tsv && tabix -f -s1 -b2 -e2 {out_vcf}.dropped.tsv.gz")
    sh(f"bcftools view -T ^{out_vcf}.dropped.tsv.gz -Oz -o {out_vcf} {vcf} && tabix -f -p vcf {out_vcf}")
    return len(drop)
