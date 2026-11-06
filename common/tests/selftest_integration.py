"""python -m yhsd selftest : synthetic tests for the whole library (no data, no network, no external tools)."""
import os, subprocess, sys, tempfile, types
import numpy as np
import pandas as pd


def _ok(name, cond, extra=''):
    print(('PASS ' if cond else 'FAIL ') + name, extra)
    if not cond:
        raise SystemExit(1)


def simulate_family(rng, L=30_000_000, locus=12_000_000, depth_off=None, depth_par=15.0, n_sites=6000, het_frac=0.35,
                    n_m=12, n_f=12, rate=1e-8):
    """One XY family: father het at het_frac of sites (X and Y haplotypes differ), mother hom. Sons inherit the Y haplotype at the sex
    locus, daughters the X; crossovers (Poisson, 1 cM/Mb) move the inherited haplotype away from the locus. Genotypes are called
    from simulated reads at the given depths. Returns pos, dose (sites x [father, mother, offspring...]), is_male."""
    pos = np.sort(rng.choice(L, n_sites, replace=False))
    n_off = n_m + n_f
    hx = rng.integers(0, 2, n_sites); hy = np.where(rng.random(n_sites) < het_frac, 1 - hx, hx)         # X / Y alleles of father
    mom = rng.integers(0, 2, n_sites)                                                                    # mother homozygous
    is_male = np.array([True] * n_m + [False] * n_f)
    dose_true = np.zeros((n_sites, 2 + n_off))
    dose_true[:, 0] = hx + hy; dose_true[:, 1] = 2 * mom
    for j in range(n_off):
        start = hy if is_male[j] else hx                                      # haplotype at the locus
        other = hx if is_male[j] else hy
        n_x = rng.poisson(rate * L)
        cps = np.sort(rng.uniform(0, L, n_x))
        flips = np.searchsorted(cps, pos) - np.searchsorted(cps, locus)
        use_start = (flips % 2 == 0)
        pat = np.where(use_start, start, other)
        dose_true[:, 2 + j] = pat + mom
    d_off = depth_off if depth_off is not None else np.full(n_off, 7.0)
    depths = np.concatenate([[depth_par, depth_par], d_off])
    called = np.full_like(dose_true, np.nan)
    for s in range(dose_true.shape[1]):
        dp = rng.poisson(depths[s], n_sites)
        alt = rng.binomial(dp, dose_true[:, s] / 2.0)
        ref = dp - alt
        g = np.where((alt > 0) & (ref > 0), 1.0, np.where(alt == 0, 0.0, 2.0))
        g[dp < 3] = np.nan
        called[:, s] = g
    return pos, called, is_male


def main():
    from sd_regions.inheritance import family_scan
    from sd_regions.popgen import sex_summary, depth_ratio
    from variant_calling.callset import classify_counts, resolve_samples
    rng = np.random.default_rng(21)

    # --- inheritance / popgen power at the real offspring depths (YH-like: median ~7x, 14/48 below 5x)
    yh_depths = np.array([7.6, 4.4, 4.5, 8.2, 7.4, 3.7, 9.5, 7.4, 6.5, 8.7, 12.3, 6.8, 9.6, 4.8, 7.7, 11.6, 5.7, 8.7, 9.3, 3.9, 10.4, 3.5, 10.9, 12.1])
    L, locus, size = 30_000_000, 15_000_000, 100_000
    nears, fars = [], []
    for rep in range(8):                                      # average a few families: one family is noisy
        pos, dose, is_male = simulate_family(rng, L, locus, depth_off=yh_depths)
        res = family_scan(pos, dose, dict(father=0, mother=1, offspring=list(range(2, dose.shape[1]))), is_male, size)
        w = res['XY'][1].copy(); w['dist'] = np.abs(w.start + size / 2 - locus)
        nears.append(w[w.dist <= 1_000_000].med_conc.mean()); fars.append(w[w.dist >= 10_000_000].med_conc.mean())
    # NOTE: with 24 meioses at ~1 cM/Mb the test localises a locus only to roughly +/-5-10 Mb; depth is not the limit.
    _ok('XY inheritance: concordance falls with distance from the locus', np.mean(nears) > np.mean(fars) + 0.06,
        f'median concordance near {np.mean(nears):.3f} vs far (>=10 Mb) {np.mean(fars):.3f}')
    zw = res['ZW'][1]
    zmax = float(np.nan_to_num(zw.max_nlp.max())) if len(zw) else 0.0       # the simulated mother is homozygous everywhere: no ZW-informative sites
    _ok('ZW test finds nothing for an XY system', zmax < w.max_nlp.max() - 1, f'ZW max nlp {zmax:.1f} vs XY {w.max_nlp.max():.1f}')
    tied = w[w.frac_perfect >= w.frac_perfect.max() - 0.02]            # many windows tie: report the centre of the plateau, not one window
    top = tied.start.median()
    _ok('centre of the XY plateau is within 5 Mb of the planted locus', abs(top - locus) <= 5_000_000, f'plateau centre {top:.0f} ({len(tied)} tied windows)')

    # --- male/female Fst among offspring: higher near the locus (sons carry Y, daughters X haplotypes)
    ss = sex_summary(pos, dose[:, 2:], is_male, size)
    ss['dist'] = np.abs(ss.start + size / 2 - locus)
    nearf = ss[ss.dist <= 1_000_000].fst.mean(); farf = ss[ss.dist >= 10_000_000].fst.mean()
    _ok('windowed M-F Fst higher near locus', nearf > farf + 0.05, f'Fst near {nearf:.2f} vs far {farf:.2f}')

    # --- depth ratio detects a male-specific duplication (gsdf-like, 1.5x in males)
    n_w = 300; depth = rng.gamma(40, 0.2, (n_w, 24)) * rng.uniform(0.7, 1.3, 24); depth[150, :12] *= 1.5
    lr, p = depth_ratio(depth, np.array([True] * 12 + [False] * 12))
    _ok('depth ratio finds the male duplication', np.argmax(lr) == 150 and p[150] < 1e-3, f'log2 ratio {lr[150]:.2f}, p={p[150]:.1e}')

    # --- analysis step on arrays: scan_arrays -> call_signals -> merge_regions finds the planted region
    from sd_regions.parental_sex_region.steps import scan_arrays, call_signals, merge_regions, robust_z
    names = ['father', 'mother'] + [f'o{i}' for i in range(24)]
    off = pd.DataFrame(dict(SampleID=names[2:], Sex=['M'] * 12 + ['F'] * 12, BroodID='B1'))
    fam = {'B1': dict(father='father', mother='mother', offspring=names[2:], male=[True] * 12 + [False] * 12)}
    cfgs = dict(window=size, min_informative_offspring=12)
    pos2, dose2, ism2 = simulate_family(rng, L, locus, depth_off=yh_depths)
    wt = scan_arrays('chrA', pos2, dose2, names, off, fam, cfgs)
    _ok('scan_arrays produces inheritance, Fst and heterozygosity columns',
        {'fst_all', 'het_M', 'XY_B1_frac_perfect', 'ZW_B1_n_inf', 'fst_B1'} <= set(wt.columns), f'{wt.shape[1]} columns, {len(wt)} windows')
    wt['het_absdiff'] = wt.het_diff_M_minus_F.abs()
    rc = dict(fst_z=3.0, het_z=3.0, inheritance_frac_perfect=0.8, depth_p=1e-3, min_signals=1, merge_gap_windows=2)
    ws = call_signals(wt, rc)
    reg = merge_regions(ws, size, 1, 2)
    hit = reg[(reg.start <= locus) & (reg.end >= locus)]
    _ok('candidate region covers the planted locus', len(hit) == 1 and (hit.end - hit.start).iloc[0] <= 15_000_000,
        f'region {int(hit.start.iloc[0])}-{int(hit.end.iloc[0])}' if len(hit) else 'none')
    _ok('robust_z flags an outlier', robust_z(np.r_[np.random.default_rng(0).normal(0, 1, 500), 12.0])[-1] > 6)

    # --- sample_table: Blumer join (key, coding, validation of Inversion10 and Sex)
    from sd_genotypes.sample_table.inversions import annotate
    db = pd.DataFrame(dict(SampleID=['SAMEA1', 'SAMEA2', 'LAB1', 'Unnamed_ok'], Sex=['M', 'F', 'M', 'U'], Species=['a', 'b', 'c', 'd'],
                           Inversion10=[1.0, 0.0, 2.0, np.nan]))
    db['Unnamed: 14'] = np.nan
    s1 = pd.DataFrame(dict(sequence_id=['cichlid1', 'cichlid2', 'cichlid3'], specimen_id=['x', 'y', 'z'], Biosample=['SAMEA1', 'SAMEA2', 'SAMEA9'], sex=['M', 'F', 'U']))
    s7 = pd.DataFrame(dict(sample_id=['cichlid1', 'cichlid2', 'cichlid3'], chr2=[0, 2, 0], chr9=[2, 2, 0], chr10=[1, 0, 0], chr11=[0, 1, 0], chr13=[0, 0, 2]))
    tab, rep = annotate(db, s1, s7)
    _ok('annotate: matched, validated, lab sample left blank', rep['n_matched'] == 2 and rep['inv10_discordant'] == 0 and rep['sex_discordant'] == 0
        and tab.loc[tab.SampleID == 'SAMEA2', 'Inversion2'].iloc[0] == 2 and np.isnan(tab.loc[tab.SampleID == 'LAB1', 'Inversion9'].iloc[0])
        and 'Unnamed: 14' not in tab.columns and len(rep['s1_not_in_db']) == 1)
    db.loc[0, 'Inversion10'] = 0.0
    _, rep2 = annotate(db, s1, s7)
    _ok('annotate: flags an Inversion10 discordance', rep2['inv10_discordant'] == 1 and rep2['inv10_discordant_ids'] == ['SAMEA1'])

    # --- BC inversion inference from haplotype calls
    from sd_genotypes.sample_table import bc_inversions as bci
    n_fish, grid_pos = 60, np.arange(0, 30_000_000, 100_000)
    states = rng.integers(0, 2, n_fish)
    R = np.zeros((len(grid_pos), n_fish), int)
    for j in range(n_fish):
        cp = rng.uniform(0, 3e6) if rng.random() < .5 else rng.uniform(26e6, 30e6)       # crossovers only outside the block (6-24 Mb)
        R[:, j] = np.where(grid_pos < cp, 1 - states[j], states[j]) if cp < 3e6 else np.where(grid_pos > cp, 1 - states[j], states[j])
    R[(grid_pos >= 6e6) & (grid_pos <= 24e6), :] = states[None, :]                      # block: one haplotype per fish
    R[(grid_pos >= 10e6) & (grid_pos <= 11e6), 5] = 1 - states[5]                       # fish 5: ambiguous inside the block
    R[(grid_pos >= 10e6) & (grid_pos <= 17e6), 5] = 1 - states[5]
    grid = pd.DataFrame(R, columns=[f'f{i}' for i in range(n_fish)]); grid.insert(0, 'Position', grid_pos); grid.insert(0, 'Chromosome', 'chrX')
    calls, blocks = bci.infer(grid, {'LG2': 'chrX'}, {'LG2': False}, max_rate=0.03, min_len=10e6, min_consistency=0.8)
    b = blocks.iloc[0]
    _ok('BC inversion block detected (contains the true block, excludes the crossover-rich ends)', b.found and b.block_start <= 8e6 and b.block_end >= 22e6 and (b.block_end - b.block_start) < 28e6, f'{b.block_start}-{b.block_end}')
    ok = [i for i in range(n_fish) if i != 5]
    _ok('BC inversion genotype: YH haplotype -> 1, MC haplotype -> 0', all(calls.loc[f'f{i}', 'Inversion2'] == states[i] for i in ok))
    _ok('BC inversion: inconsistent fish left blank', np.isnan(calls.loc['f5', 'Inversion2']) or calls.loc['f5', 'consistency_LG2'] >= 0.8)

    # --- pedigree-expected values and the blank-cell sheet writer
    from sd_genotypes.sample_table.steps import fill_blank_cells
    tb = pd.DataFrame(dict(SampleID=['b1', 'b2', 'w1'], Species=['YHxMC x MC BC', 'YHxMC x MC BC', 'wild'], Inversion9=[np.nan, np.nan, np.nan],
                           Inversion10=[np.nan, np.nan, 1.0], Inversion11=[np.nan, np.nan, np.nan], InversionSource=['inferred_BC_haplotype:chr2,chr13', '', 'Blumer2025_TableS7']))
    out, filled = bci.apply_expected(tb, {'Inversion9': 2, 'Inversion10': 0, 'Inversion11': 2}, {'LG9': True, 'LG10': True, 'LG11': False},
                                     {'Inversion9': 'LG9', 'Inversion10': 'LG10', 'Inversion11': 'LG11'})
    _ok('apply_expected fills BC rows only, respects the data check, labels the source',
        filled == ['Inversion9', 'Inversion10'] and out.loc[0, 'Inversion9'] == 2 and np.isnan(out.loc[0, 'Inversion11']) and np.isnan(out.loc[2, 'Inversion9'])
        and out.loc[2, 'Inversion10'] == 1.0 and out.loc[0, 'InversionSource'] == 'inferred_BC_haplotype:chr2,chr13;pedigree_expected:chr9,chr10'
        and out.loc[1, 'InversionSource'] == 'pedigree_expected:chr9,chr10' and out.loc[2, 'InversionSource'] == 'Blumer2025_TableS7')
    class _WS:
        def __init__(s, v): s.v, s.batches = v, []
        def get_all_values(s): return s.v
        def batch_update(s, data, value_input_option=None): s.batches.append(data)
    ws = _WS([['SampleID', 'Species', 'Inversion10'], ['b1', 'x', ''], ['w1', 'x', '1'], ['b2', 'x', '']])
    pj = types.SimpleNamespace(fm=types.SimpleNamespace(s_ID='x', gc=types.SimpleNamespace(open_by_key=lambda k: types.SimpleNamespace(worksheet=lambda n: ws))))
    tb2 = out.copy(); tb2.loc[tb2.SampleID == 'w1', 'Inversion10'] = 0.0                         # would overwrite the non-empty cell if the writer were careless
    n = fill_blank_cells(pj, tb2, 'Inversion10', tempfile.mkdtemp() + '/', execute=True)
    _ok('fill_blank_cells writes only blank cells (never overwrites)', n == 2 and ws.batches[0] == [{'range': 'C2', 'values': [[0]]}, {'range': 'C4', 'values': [[0]]}], str(ws.batches))

    # --- BC genotype file is fetched through the FileManager when missing
    from sd_genotypes.sample_table.steps import resolve_bc_genotypes
    base = tempfile.mkdtemp() + '/'
    calls_dl = []
    class _FM: localMasterDir = base; localNikeshDir = base + '2026_Nikesh/'
    class _PJ:
        fm = _FM()
        def ensure_local(self, path):
            calls_dl.append(path); os.makedirs(os.path.dirname(path), exist_ok=True); open(path, 'w').write('Chromosome,Position\n')
    got = resolve_bc_genotypes(_PJ(), 'legacy')
    _ok("resolve_bc_genotypes('legacy') downloads via the FileManager when missing", got == base + '2026_Nikesh/OutputGenotypes.csv' and calls_dl == [got])
    resolve_bc_genotypes(_PJ(), got)
    _ok('resolve_bc_genotypes does not re-download an existing file', len(calls_dl) == 1)
    try:
        resolve_bc_genotypes(_PJ(), '/nonexistent/elsewhere.csv'); ok_err = False
    except FileNotFoundError:
        ok_err = True
    _ok('resolve_bc_genotypes gives a clear error outside the master directory', ok_err)

    # --- sex_association: planted XY (gsdf dup in Mbuna), planted ZW (LG5 in Mbuna), fixed marker, underpowered clade, brood and uncertain-sex exclusions
    from sd_regions.sex_association import steps as sa
    rg = np.random.default_rng(3)
    def mk(group, n_m, n_f, car_m, car_f, col, val_car, val_non, extra=None):
        rows = []
        for sx, n, nc in (('M', n_m, car_m), ('F', n_f, car_f)):
            for i in range(n):
                rows.append({'SampleID': f'{group}{sx}{i}', 'Sex': sx, 'Subgroup': group, col: val_car if i < nc else val_non, **(extra or {})})
        return pd.DataFrame(rows)
    d1 = mk('Mbuna', 20, 20, 18, 1, 'gsdf number', 3.0, 2.0)                       # carriers male: should be supported
    d2 = mk('Utaka', 20, 20, 0, 0, 'gsdf number', 3.0, 2.0)                        # no carriers: not testable
    d3 = mk('Eukambuzi', 3, 3, 3, 0, 'gsdf number', 3.0, 2.0)                      # tiny: cannot reach significance
    d4 = mk('SB_Other', 20, 20, 10, 10, 'gsdf number', 3.0, 2.0)                   # carriers unrelated to sex: no evidence
    tab = pd.concat([d1, d2, d3, d4], ignore_index=True)
    tab['Ecogroup'] = ''
    z = mk('Mbuna', 20, 20, 1, 18, 'Inversion5', 1.0, 0.0)                         # W-like carriers female: LG5
    z2 = tab[tab.Subgroup == 'Mbuna'].reset_index(drop=True).copy()
    tab['Inversion5'] = np.nan; tab.loc[tab.Subgroup == 'Mbuna', 'Inversion5'] = z['Inversion5'].values
    extra = mk('Mbuna', 10, 0, 10, 0, 'gsdf number', 3.0, 2.0); extra['BroodID'] = 'B1'; extra['SampleID'] = ['b' + str(i) for i in range(len(extra))]   # lab siblings: excluded
    unc = mk('Mbuna', 0, 5, 0, 5, 'gsdf number', 3.0, 2.0); unc['Sex'] = 'F?'; unc['SampleID'] = ['u' + str(i) for i in range(len(unc))]                 # uncertain sex: excluded
    full = pd.concat([tab, extra, unc], ignore_index=True)
    cfg = sa.load_config(); cfg['markers'] = {k: v for k, v in cfg['markers'].items() if k in ('gsdf_dup', 'lg5')}
    res, flow = sa.run_tests(full, cfg)
    g = lambda m, grp: res[(res.marker == m) & (res.group == grp)].iloc[0]
    _ok('sex_association: planted XY marker is supported (FDR)', g('gsdf_dup', 'Mbuna').evidence == 'supports sex determiner' and g('gsdf_dup', 'Mbuna').n_used == 40, f"q={g('gsdf_dup', 'Mbuna').q:.1e}")
    _ok('sex_association: no carriers -> not testable', g('gsdf_dup', 'Utaka').status == 'not testable: no variation')
    _ok('sex_association: tiny group is flagged, not claimed', g('gsdf_dup', 'Eukambuzi').status.startswith('not testable') or g('gsdf_dup', 'Eukambuzi').status == 'underpowered')
    _ok('sex_association: unrelated marker shows no evidence', g('gsdf_dup', 'Other shallow benthics').evidence.startswith('no evidence'))
    _ok('sex_association: ZW marker (W carriers female) is supported', g('lg5', 'Mbuna').evidence == 'supports sex determiner' and g('lg5', 'Mbuna').carrier_sex_expected == 'F')
    _ok('sex_association: sibling offspring and uncertain sexes are excluded and counted', set(flow[flow.marker == 'gsdf_dup'].dropped_because) >= {'sex not M/F', 'lab-reared offspring (BroodID set)'}
        and int(flow[(flow.marker == 'gsdf_dup') & (flow.dropped_because == 'lab-reared offspring (BroodID set)')].n_samples.iloc[0]) == 10)
    sp = full.copy(); sp['LabReared'] = np.where(sp.SampleID.str.contains('F1|F2|F3|F4|F5'), 1, np.nan)
    cfg2 = dict(cfg, split_lab_wild=True); res2, _ = sa.run_tests(sp, cfg2)
    _ok('sex_association: lab and wild samples are tested separately', {'Mbuna (lab)', 'Mbuna (wild)'} <= set(res2.group), str(sorted(set(res2.group))[:4]))
    gt = sa.to_number(pd.Series(['0/0', '0/1', '1/1', './.', np.nan, '0|1']), 'gt')
    _ok('sex_association: genotype strings -> alt dosage', list(gt.fillna(-1)) == [0, 1, 2, -1, -1, 1])
    _ok('sex_association: BH q-values are monotone and <= 1', (lambda q: np.all(np.diff(q[np.argsort([0.001, 0.01, 0.04, 0.5])]) >= -1e-12) and q.max() <= 1)(sa.bh([0.001, 0.01, 0.04, 0.5])))

    # --- genome registry helpers and the inversion-breakpoint pipeline on synthetic alignments

    # --- callset helpers
    _ok('classify_counts', list(classify_counts([20, 10, 0, 15, 3], [0, 10, 20, 2, 0])) == [0, 1, 2, -1, -1])
    dt = pd.DataFrame(dict(SampleID=['o1', 'o2', 'o3', 'p1', 'p2'], Species=['YH'] * 5, BroodID=['B1', 'B1', 'B2', None, None],
                           Parents=['p1,p2', 'p1,p2', 'p1,p2', None, None], Category=[None] * 5))
    _ok('resolve_samples', resolve_samples(dt, {'broods': ['B1', 'B2']}) == ['o1', 'o2', 'o3', 'p1', 'p2'])

    # --- runner with a fake rclone (local directory as the "Dropbox")
    fake = tempfile.mkdtemp(); cloud = tempfile.mkdtemp() + '/'
    open(fake + '/rclone', 'w').write('''#!/usr/bin/env python3
import sys, os, shutil, filecmp
a=[x for x in sys.argv[1:] if not x.startswith('--')]
cmd,src,dst=a[0],a[1],a[2]
if cmd=='copy':
    if os.path.isdir(src): shutil.copytree(src,dst,dirs_exist_ok=True)
    else: os.makedirs(dst,exist_ok=True); shutil.copy2(src,dst)
elif cmd=='check':
    bad=0
    if os.path.isdir(src):
        for r,_,fs in os.walk(src):
            for f in fs:
                s=os.path.join(r,f); d=os.path.join(dst,os.path.relpath(s,src))
                if not os.path.exists(d) or not filecmp.cmp(s,d,shallow=False): bad+=1
    sys.exit(1 if bad else 0)
''')
    os.chmod(fake + '/rclone', 0o755)
    old = os.environ['PATH']; os.environ['PATH'] = fake + ':' + old
    try:
        from cichlid_common.runner import Run, latest_dir
        loc = tempfile.mkdtemp() + '/'
        proj = types.SimpleNamespace(fm=types.SimpleNamespace(localMasterDir=loc, cloudMasterDir=cloud), root=loc + 'P/G/', genome='G', cfg={})
        os.makedirs(proj.root)
        with Run(proj, 'toy', 'step') as r:
            open(r.dir + 'out.txt', 'w').write('x')
        import shutil; shutil.rmtree(proj.root + 'toy')
        _ok('runner uploads, verifies and restores from the cloud', os.path.exists(latest_dir(proj, 'toy', 'step') + 'out.txt'))
    finally:
        os.environ['PATH'] = old
    print('ALL SELFTESTS PASSED')
