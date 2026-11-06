"""Ground-truth validation with the deep BC fish (PM26, ~7-17x): call at full depth, downsample to target coverages,
call again, and measure disagreement. Gives REAL error rates for this data type (replaces simulated numbers)."""
import os
import numpy as np
import pandas as pd
from cichlid_common.common import sh, parallel
from .v2 import count_bam, hmm_sample

_deep = {}


def _calls(proj, bam, sire, tag, lengths):
    out = f"{proj.dirs['tmp']}val.{tag}.tsv.gz"
    count_bam(proj, bam, sire, out)
    res, _ = hmm_sample(pd.read_csv(out, sep='\t'), lengths, proj.cfg['v2'])
    os.remove(out)
    return res


def _one(args):
    proj, s, sire, cov, t, seed, lengths = args
    bam = proj.ensure_bam(s)
    key = f'{s}.cov{t}.seed{seed}'
    ds = f"{proj.dirs['tmp']}{key}.bam"
    frac = min(t / cov, 0.999)                      # samtools -s INT.FRAC : INT = seed, FRAC = fraction kept
    sh(['samtools', 'view', '-b', '-@', '2', '-s', f'{seed + frac:.5f}', '-o', ds, bam])
    sh(['samtools', 'index', ds])
    full = _calls(proj, bam, sire, s + '.full', lengths)
    low = _calls(proj, ds, sire, key, lengths)
    os.remove(ds); os.remove(ds + '.bai')
    rows = []
    for ch in full:
        a = full[ch]['post'] > 0.5; b = low[ch]['post'] > 0.5
        fb = np.array([x[0] + x[1] for x in full[ch]['bps']]) / 2; lb = np.array([x[0] + x[1] for x in low[ch]['bps']]) / 2
        dist = [np.abs(lb - x).min() if len(lb) else np.nan for x in fb]
        rows.append(dict(sample=s, target_cov=t, seed=seed, chrom=ch, grid_agreement=float((a == b).mean()),
                         n_bp_full=len(fb), n_bp_low=len(lb), n_full_bp_missed=int(sum(np.isnan(dist) or d > 1e6 for d in dist)),
                         med_bp_shift_kb=float(np.nanmedian(dist) / 1e3) if len(fb) and not np.all(np.isnan(dist)) else np.nan))
    return rows


def main(proj):
    v = proj.cfg['v2']['validation']
    proj.fm.readAlignmentDatabase()
    cov = proj.fm.alignment_dt.set_index('SampleID').Coverage
    bc = proj.bc_table().set_index('SampleID')
    deep = [s for s in bc.index if s in cov.index and cov[s] >= v['min_cov']]
    print('deep BC fish for validation:', len(deep), deep)
    lengths = proj.chrom_lengths()
    jobs = [(proj, s, bc.SireSample[s], float(cov[s]), t, seed, lengths) for s in deep for t in v['target_cov'] for seed in v['seeds']]
    res, errs = parallel(_one, jobs, proj.workers('count'), 'validate')
    df = pd.DataFrame([r for rows in res.values() for r in rows])
    df['LG'] = df.chrom.map(proj.lg)
    df.to_csv(proj.dirs['qc'] + 'validation_downsample.tsv', sep='\t', index=False)
    print(df.groupby('target_cov')[['grid_agreement', 'n_full_bp_missed']].agg(['mean', 'sum']))
