"""Parental-species sex-region analysis. Steps (python -m yhsd parental_sex_region <step>):
  callset   build the per-species GATK callsets (yh_broods, mc_broods) -> Dropbox
  depth     window depth per sample (mosdepth) for male/female depth ratios (copy-number: Y-specific sequence, gsdf duplication)
  scan      per-window inheritance (XY and ZW, per family), male-female Fst, heterozygosity, SNP density, depth ratio
  summarize candidate regions (windows with >= 2 independent signals, merged) + figures
All outputs are run folders (manifest + Dropbox upload) written through the runner."""
import os
import numpy as np
import pandas as pd
import yaml
from cichlid_common.common import sh, parallel
from cichlid_common.common_light import config_path
from variant_calling import callset as cs
from cichlid_common.runner import latest_dir
from ..inheritance import alt_dosage, family_scan
from ..popgen import sex_summary, depth_ratio


def cfg_of(proj):
    return yaml.safe_load(open(config_path('parental_sex_region.yaml')))


def sample_meta(proj, scfg):
    """Offspring (sex, brood) and family structure from SampleDatabase. First parent = father, second = mother."""
    dt = proj.sample_table()
    off = dt[dt.BroodID.isin(scfg['broods'])][['SampleID', 'Sex', 'BroodID', 'Parents']].copy()
    off = off[off.Sex.isin(['M', 'F'])]
    fam = {}
    for b, g in off.groupby('BroodID'):
        f, m = [x.strip() for x in g.Parents.iloc[0].split(',')[:2]]
        fam[b] = dict(father=f, mother=m, offspring=g.SampleID.tolist(), male=(g.Sex == 'M').tolist())
    return off, fam


# ------------------------------------------------------------------ steps
def callset_step(proj, species=None):
    cfg = cfg_of(proj)
    for sp, scfg in cfg['species'].items():
        if species and sp not in species:
            continue
        cs.build(proj, scfg['callset'])


def _depth_one(args):
    proj, s, window, outdir = args
    pre = f'{proj.dirs["tmp"]}mosdepth.{s}'
    sh(['mosdepth', '-n', '-t', '2', '--by', str(window), pre, proj.ensure_bam(s)])
    df = pd.read_csv(pre + '.regions.bed.gz', sep='\t', header=None, names=['chrom', 'start', 'end', s])
    df.to_csv(f'{outdir}{s}.depth.tsv.gz', sep='\t', index=False)
    for ext in ('.regions.bed.gz', '.regions.bed.gz.csi', '.mosdepth.global.dist.txt', '.mosdepth.summary.txt'):
        if os.path.exists(pre + ext):
            os.remove(pre + ext)
    return s


def depth_step(proj, species=None):
    cfg = cfg_of(proj)
    for sp, scfg in cfg['species'].items():
        if species and sp not in species:
            continue
        dt = proj.sample_table()
        samples = cs.resolve_samples(dt, cs.load_callsets()['callsets'][scfg['callset']])
        with proj.run('parental_sex_region', f'depth_{sp}', params=dict(window=cfg['depth_window'], samples=samples), tar=('per_sample',)) as run:
            os.makedirs(run.dir + 'per_sample/', exist_ok=True)
            res, errs = parallel(_depth_one, [(proj, s, cfg['depth_window'], run.dir + 'per_sample/') for s in samples], proj.workers('count'), 'depth')
            assert not errs, errs
            mats = [pd.read_csv(f'{run.dir}per_sample/{s}.depth.tsv.gz', sep='\t') for s in samples]
            m = mats[0]
            for x in mats[1:]:
                m = m.merge(x, on=['chrom', 'start', 'end'])
            m.to_csv(run.dir + 'depth_matrix.tsv.gz', sep='\t', index=False)


def scan_arrays(chrom, pos, dose, names, off, fam, cfg):
    """Window table for one chromosome from alt-dosage arrays (sites x samples). Pure numpy/pandas: unit-tested."""
    col = {s: i for i, s in enumerate(names)}
    size = cfg['window']
    ocols = [col[s] for s in off.SampleID if s in col]
    is_male = (off.set_index('SampleID').loc[[s for s in off.SampleID if s in col]].Sex == 'M').values
    out = sex_summary(pos, dose[:, ocols], is_male, size).rename(columns={'fst': 'fst_all'})
    for b, f in fam.items():                                  # per-family inheritance (XY and ZW) and per-family Fst
        oc = [col[s] for s in f['offspring'] if s in col]; mk = np.array([m for s, m in zip(f['offspring'], f['male']) if s in col])
        res = family_scan(pos, dose, dict(father=col[f['father']], mother=col[f['mother']], offspring=oc), mk, size, cfg['min_informative_offspring'])
        for mode, (_, w) in res.items():
            w = w.rename(columns={c: f'{mode}_{b}_{c}' for c in w.columns if c != 'start'})
            out = out.merge(w, on='start', how='outer')
        fs = sex_summary(pos, dose[:, oc], mk, size)[['start', 'fst']].rename(columns={'fst': f'fst_{b}'})
        out = out.merge(fs, on='start', how='left')
    out.insert(0, 'chrom', chrom)
    return out


def _scan_chrom(args):
    vcf, chrom, off, fam, cfg = args
    import allel
    d = allel.read_vcf(vcf, region=chrom, fields=['variants/POS', 'calldata/GT', 'calldata/DP', 'calldata/GQ', 'samples'])
    if d is None:
        return None
    names = list(d['samples'])
    parents = {p for f in fam.values() for p in (f['father'], f['mother'])}
    mind = np.array([cfg['min_dp_parent'] if s in parents else cfg['min_dp_offspring'] for s in names])
    dose = alt_dosage(d['calldata/GT'], d['calldata/DP'], d['calldata/GQ'], mind, cfg['min_gq'])
    return scan_arrays(chrom, d['variants/POS'], dose, names, off, fam, cfg)


def scan_step(proj, species=None):
    cfg = cfg_of(proj)
    for sp, scfg in cfg['species'].items():
        if species and sp not in species:
            continue
        vcf = cs.vcf_path(proj, scfg['callset'])
        off, fam = sample_meta(proj, scfg)                      # sample database is read once, in the parent process
        with proj.run('parental_sex_region', f'scan_{sp}', params=dict(cfg=cfg, species=sp, vcf=vcf)) as run:
            run.add_input(vcf)
            res, errs = parallel(_scan_chrom, [(vcf, c, off, fam, cfg) for c in proj.chroms], proj.workers('count'), f'scan-{sp}')
            assert not errs, errs
            w = pd.concat([r for r in res.values() if r is not None], ignore_index=True)
            w['LG'] = w.chrom.map(proj.lg)
            w['het_absdiff'] = w.het_diff_M_minus_F.abs()
            try:                                                    # depth ratio, if the depth step has been run
                dm = pd.read_csv(latest_dir(proj, 'parental_sex_region', f'depth_{sp}') + 'depth_matrix.tsv.gz', sep='\t')
                sx = off.set_index('SampleID').Sex
                cols = [c for c in dm.columns[3:] if c in sx.index]
                lr, p = depth_ratio(dm[cols].values, (sx[cols] == 'M').values)
                dd = dm[['chrom', 'start']].assign(depth_log2_MF=lr, depth_p=p)
                dd['start'] = (dd.start // cfg['window']) * cfg['window']
                w = w.merge(dd, on=['chrom', 'start'], how='left')
            except FileNotFoundError:
                print('[scan] no depth run found; skipping depth ratio')
            w.to_csv(run.dir + 'windows.tsv.gz', sep='\t', index=False)


# ------------------------------------------------------------------ summarize (pure functions are unit-tested)
def robust_z(x):
    x = np.asarray(x, float); med = np.nanmedian(x); mad = np.nanmedian(np.abs(x - med)) * 1.4826
    return (x - med) / (mad if mad > 0 else np.nanstd(x) or 1.0)


def call_signals(w, rc):
    """Boolean signal columns per window and their count."""
    w = w.copy()
    sig = pd.DataFrame(index=w.index)
    sig['fst'] = robust_z(w.fst_all) > rc['fst_z']
    sig['het'] = robust_z(w.het_absdiff) > rc['het_z']
    inh = [c for c in w.columns if c.endswith('_frac_perfect')]
    sig['inheritance'] = (w[inh].max(axis=1) >= rc['inheritance_frac_perfect']) if inh else False
    sig['depth'] = (w.depth_p < rc['depth_p']) & (w.depth_log2_MF.abs() > 0.3) if 'depth_p' in w else False
    w = pd.concat([w, sig.add_prefix('sig_')], axis=1)
    w['n_signals'] = sig.sum(axis=1)
    return w


def merge_regions(w, size, min_signals, gap):
    """Merge windows with n_signals >= min_signals (per chromosome, gaps <= `gap` windows) into regions."""
    rows = []
    for ch, g in w[w.n_signals >= min_signals].groupby('chrom'):
        starts = np.sort(g.start.values)
        if len(starts) == 0:
            continue
        s0 = prev = starts[0]
        for s in list(starts[1:]) + [None]:
            if s is None or s - prev > gap * size:
                sub = g[(g.start >= s0) & (g.start <= prev)]
                rows.append(dict(chrom=ch, start=int(s0), end=int(prev + size), n_windows=len(sub), max_signals=int(sub.n_signals.max()),
                                 max_fst=float(sub.fst_all.max()), signals=','.join(sorted({c[4:] for c in sub.columns if c.startswith('sig_') and sub[c].any()}))))
                if s is not None:
                    s0 = s
            if s is not None:
                prev = s
    return pd.DataFrame(rows)


def summarize_step(proj, species=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cfg = cfg_of(proj); rc = cfg['region_calling']
    for sp, scfg in cfg['species'].items():
        if species and sp not in species:
            continue
        wdir = latest_dir(proj, 'parental_sex_region', f'scan_{sp}')
        w = pd.read_csv(wdir + 'windows.tsv.gz', sep='\t')
        w = call_signals(w, rc)
        reg = merge_regions(w, cfg['window'], rc['min_signals'], rc['merge_gap_windows'])
        if len(reg):
            reg['LG'] = reg.chrom.map(proj.lg)
        with proj.run('parental_sex_region', f'summary_{sp}', params=dict(region_calling=rc, scan_run=wdir)) as run:
            run.add_input(wdir + 'windows.tsv.gz')
            w.to_csv(run.dir + 'windows_with_signals.tsv.gz', sep='\t', index=False)
            reg.to_csv(run.dir + 'candidate_regions.tsv', sep='\t', index=False)
            panels = [('fst_all', 'Fst (M vs F)'), ('het_diff_M_minus_F', 'heterozygosity M - F')]
            inh = [c for c in w.columns if c.startswith('XY_') and c.endswith('_frac_perfect')]
            panels += [(c, c.replace('_frac_perfect', '')) for c in inh]
            if 'depth_log2_MF' in w:
                panels.append(('depth_log2_MF', 'depth log2(M/F)'))
            lgs = list(dict.fromkeys(w.LG)); off = {}; acc = 0
            chrom_len = w.groupby('LG').start.max()
            for lg in lgs:
                off[lg] = acc; acc += chrom_len[lg] + 2e6
            x = w.start + w.LG.map(off)
            fig, axs = plt.subplots(len(panels), 1, figsize=(11, 1.9 * len(panels)), sharex=True)
            for ax, (c, lab) in zip(np.atleast_1d(axs), panels):
                col = np.where(w.LG.map(lambda z: lgs.index(z) % 2), '#7f8c8d', '#34495e')
                ax.scatter(x / 1e6, w[c], s=3, c=col); ax.set_ylabel(lab, fontsize=7)
                exp = scfg['expected_lg']
                if exp in off:
                    ax.axvspan(off[exp] / 1e6, (off[exp] + chrom_len[exp]) / 1e6, color='#f1c40f', alpha=.25, lw=0)
            np.atleast_1d(axs)[-1].set_xticks([(off[l] + chrom_len[l] / 2) / 1e6 for l in lgs]); np.atleast_1d(axs)[-1].set_xticklabels(lgs, fontsize=6, rotation=90)
            np.atleast_1d(axs)[0].set_title(f'{sp}: window statistics (shaded: expected {scfg["expected_lg"]})', fontsize=9)
            plt.tight_layout(); plt.savefig(run.dir + f'{sp}_genome_scan.png', dpi=200); plt.close()
            print(sp, 'candidate regions:\n', reg.to_string(index=False) if len(reg) else 'none')
