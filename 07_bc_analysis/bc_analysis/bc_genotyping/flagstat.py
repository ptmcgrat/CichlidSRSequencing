"""Download BAMs (only <sample>.all.bam) and re-run flagstat + per-chromosome MAPQ/duplicate summaries."""
import json, os
import pandas as pd
from cichlid_common.common import sh, parallel

_P = None


def _one(args):
    proj, s = args
    bam = proj.ensure_bam(s)
    fs = json.loads(sh(['samtools', 'flagstat', '-@', '2', '-O', 'json', bam]).stdout)['QC-passed reads']
    rows = [dict(sample=s, chrom='ALL', total=fs['total'], mapped=fs['mapped'], duplicates=fs['duplicates'],
                 properly_paired=fs['properly paired'], nondup_primary=None, nondup_mapq20=None)]
    for c in proj.chroms:   # 0xD04 = unmapped+secondary+supplementary+duplicate
        n = int(sh(['samtools', 'view', '-c', '-F', '0xD04', bam, c]).stdout)
        q = int(sh(['samtools', 'view', '-c', '-F', '0xD04', '-q', '20', bam, c]).stdout)
        rows.append(dict(sample=s, chrom=c, nondup_primary=n, nondup_mapq20=q))
    return rows


def main(proj, samples=None):
    if samples is None:
        bc = proj.bc_table().SampleID.tolist()
        cfg = proj.cfg
        deep = [cfg['yh_sire'], cfg['f1_dam']] + cfg['bc_dams'] + cfg['mc_panel'] + list(cfg['f1_sire_names'].values())
        samples = list(dict.fromkeys(bc + deep))
    res, errs = parallel(_one, [(proj, s) for s in samples], proj.workers('flagstat'), 'flagstat')
    df = pd.DataFrame([r for rows in res.values() for r in rows])
    df['LG'] = df.chrom.map(proj.lg)
    df.to_csv(proj.dirs['qc'] + 'flagstat_and_mapq.tsv', sep='\t', index=False)
    a = df[df.chrom == 'ALL'].copy()
    a['dup_rate'] = a.duplicates / a.total
    a['mapped_rate'] = a.mapped / a.total
    c = df[df.chrom != 'ALL'].copy()
    c['frac_mapq20'] = c.nondup_mapq20 / c.nondup_primary
    a.to_csv(proj.dirs['qc'] + 'flagstat_summary.tsv', sep='\t', index=False)
    c.groupby('chrom').frac_mapq20.describe().to_csv(proj.dirs['qc'] + 'mapq20_fraction_by_chrom.tsv', sep='\t')
    print(a[['dup_rate', 'mapped_rate']].describe())
    print(c.groupby('chrom').frac_mapq20.median().sort_values().head(5), '\n(lowest median MAPQ>=20 fraction by chromosome)')
    if errs:
        pd.Series(errs).to_csv(proj.dirs['qc'] + 'flagstat_errors.tsv', sep='\t')
