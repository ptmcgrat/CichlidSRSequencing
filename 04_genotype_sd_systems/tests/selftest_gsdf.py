"""Synthetic tests of the gsdf copy-number analysis and the panel (called from selftest). samtools is replaced by a fake executable that prints the depth a sample with a given copy number
would have, so the measurement code, the JSON, the command line, the calling and the validation are exercised end to end."""
import json, os, stat, subprocess, sys, tempfile, types
import numpy as np
import pandas as pd
import yaml


def _ok(name, cond, extra=''):
    print(('PASS ' if cond else 'FAIL ') + name, extra)
    if not cond:
        raise SystemExit(1)


def _raises(fn, exc=Exception):
    try:
        fn()
    except exc:
        return True
    return False


FAKE_SAMTOOLS = r'''#!/usr/bin/env python3
import os, sys
a = sys.argv[1:]
cfg_chrom, start, end, ul, ur, jl, jr_ = "NC_135173.1", 50700000, 50830000, 50755000, 50785000, 50772650, 50773150
if a[0] == "--version":
    print("samtools 1.19 (fake)"); sys.exit(0)
if a[:2] == ["view", "-H"]:
    print("@HD\tVN:1.6\n@SQ\tSN:" + os.environ.get("FAKE_CONTIG", cfg_chrom) + "\tLN:68000000"); sys.exit(0)
if a[0] == "depth":
    q = "-Q" in a
    ratio = float(os.environ["FAKE_RATIO"]); jr = float(os.environ["FAKE_JUNCTION"]); flank = float(os.environ.get("FAKE_FLANK", "15"))
    for pos in range(start + 1, end + 1):
        d = flank * ratio if ul <= pos - 1 < ur + 5000 else flank
        if q and jl <= pos - 1 < jr_:
            d = flank * jr
        print(f"{cfg_chrom}\t{pos}\t{d:g}")
    sys.exit(0)
sys.exit(2)
'''


def run():
    from sd_genotypes.gsdf import measure as M, call as C, validate as V, steps as S
    from sd_genotypes.panel import measure as PM
    cfg = S.cfg_of()
    r = cfg['region']; L = r['end'] - r['start']
    u0, u1 = [x - r['start'] for x in cfg['unit_window']]; j0, j1 = [x - r['start'] for x in cfg['junction_window']]

    def depth_arrays(ratio, jr, flank=15.0):
        d = np.full(L, float(flank)); d[u0:u1 + 5000] = flank * ratio                       # the real block runs to 50,794,777: the bin starting at 50,785,000 is inside it too
        q = d.copy(); q[j0:j1] = flank * jr
        return d, q

    def raw_of(sid, ratio, jr, flank=15.0):
        d, q = depth_arrays(ratio, jr, flank)
        x = M.raw_measurements(d, q, cfg); x.update(SampleID=sid, definition=M.definition_key(cfg), bam_size=1, measured='t', panel_version='t'); return x

    # ======================================================================================== measurement, pure
    _ok('region string is 1-based and inclusive', M.region_string(cfg) == 'NC_135173.1:50700001-50830000')
    _ok('depth command: zeros included, region, optional MAPQ', M.depth_command('x.bam', 'c:1-2') == ['samtools', 'depth', '-a', '-r', 'c:1-2', 'x.bam']
        and M.depth_command('x.bam', 'c:1-2', 30)[-3:-1] == ['-Q', '30'])
    arr = M.parse_depth('c\t11\t5\nc\t12\t7\nc\t99999\t3\n', 10, 20)
    _ok('samtools depth parsed into a 0-based array (positions without a line are 0, out of range ignored)', arr[0] == 5 and arr[1] == 7 and arr[2:].sum() == 0 and len(arr) == 10)
    bins = np.repeat(np.arange(26, dtype=float), 5000)
    bm = M.bin_means(bins, 5000)
    unit = M.window_mean_of_bins(bm, r['start'], 5000, [cfg['unit_window']]); flank = M.window_mean_of_bins(bm, r['start'], 5000, cfg['flank_windows'])
    _ok('a bin belongs to a window when it STARTS inside it (sexloci, the calibration): unit = 7 bins (11-17), flank = 8 + 4 bins (0-7, 20-23)', unit == np.mean(range(11, 18)) and flank == np.mean(list(range(0, 8)) + [20, 21, 22, 23]), f'{unit} {flank}')
    inside = M.mean_of_bins_inside(bm, r['start'], 5000, cfg['junction_flank_window'])
    _ok('the junction ratio divides by the bins ENTIRELY inside the first flank window (bins 0-6)', inside == np.mean(range(0, 7)))
    _ok('the unit window holds none of the flank bins (no overlap)', not set(range(11, 18)) & set(list(range(0, 8)) + [20, 21, 22, 23]))
    k = M.definition_key(cfg)
    c2 = yaml.safe_load(yaml.safe_dump(cfg)); c2['calibration']['ratio_per_copy'] = 0.3; c3 = yaml.safe_load(yaml.safe_dump(cfg)); c3['unit_window'] = [50755000, 50790000]
    _ok('definition key: unchanged by calibration or thresholds, changed by anything that changes what is measured', M.definition_key(c2) == k and M.definition_key(c3) != k)

    # ======================================================================================== calling
    cal = cfg['calibration']
    _ok('copies = 2 + (ratio - 0.530) / 0.268', abs(C.estimate_copies(0.530, cal) - 2) < 1e-9 and abs(C.estimate_copies(0.530 + 0.268, cal) - 3) < 1e-9 and abs(C.estimate_copies(0.530 + 2 * 0.268, cal) - 4) < 1e-9)
    cases = {'two': (0.530, 0.0, 15.0), 'three': (0.798, 0.5, 15.0), 'four': (1.066, 1.0, 15.0), 'five': (1.33, 1.0, 15.0), 'lowdepth': (0.798, 0.5, 2.0), 'half': (0.664, 0.25, 15.0), 'depth3_junction2': (0.798, 0.0, 15.0), 'two_junction': (0.530, 0.5, 15.0), 'two_small_junction': (0.530, 0.10, 15.0), 'six': (1.60, 1.0, 15.0), 'seven': (1.87, 1.0, 15.0)}
    t = C.call_table([raw_of(k, *v) for k, v in cases.items()], cfg).set_index('SampleID')
    _ok('2, 3 and 4 copies are called and named -/-, dup/-, dup/dup, unflagged', [(int(t.loc[k].gsdf_number), t.loc[k].genotype, t.loc[k].call_flags) for k in ('two', 'three', 'four')] == [(2, '-/-', ''), (3, 'dup/-', ''), (4, 'dup/dup', '')])
    _ok('5 copies (they exist in the sheet): reported as 5, named CN5, NOT flagged (sexloci flags only estimates below 0.6 or above 5.4)', t.loc['five'].gsdf_number == 5 and t.loc['five'].genotype == 'CN5' and t.loc['five'].call_flags == '')
    _ok('6 copies is flagged outside the calibrated range; an estimate of ~7 is reported as 6 (clipped) and flagged', t.loc['six'].gsdf_number == 6 and t.loc['six'].call_flags == 'outside_calibrated_range' and t.loc['seven'].gsdf_number == 6 and 'outside_calibrated_range' in t.loc['seven'].call_flags)
    _ok('flank depth below 3x is flagged (the call itself is still made)', t.loc['lowdepth'].call_flags == 'low_flank_depth' and t.loc['lowdepth'].gsdf_number == 3)
    _ok('an estimate halfway between integers is flagged non_integer', 'non_integer' in t.loc['half'].call_flags and t.loc['half'].gsdf_number in (2, 3))
    _ok('the junction check is one-sided: a unique-read junction in a sample called with two copies is reported in check_flags (the call stays unflagged); a small ratio (0.10) is not',
        t.loc['two_junction'].check_flags == 'junction_in_two_copy' and t.loc['two_junction'].call_flags == '' and t.loc['two_small_junction'].check_flags == '')
    _ok('a three-copy sample WITHOUT junction reads is not flagged (true three-copy samples are often at 0 at 5-6x depth)', t.loc['depth3_junction2'].check_flags == '' and t.loc['depth3_junction2'].call_flags == '')
    _ok('no column holds a junction-implied copy number any more (its formula assumed 0.5 and 1.0)', 'junction_implied_copies' not in t.columns)
    _ok('the junction ratio is 0, 0.5, 1 for 2, 3, 4 copies', [round(t.loc[k].junction_ratio, 2) for k in ('two', 'three', 'four')] == [0.0, 0.5, 1.0])
    z = C.call_table([raw_of('zero', 0.5, 0.0, flank=0.0)], cfg).iloc[0]
    _ok('no coverage at all gives no estimate and flags, without crashing', np.isnan(z.gsdf_number) and 'low_flank_depth' in z.call_flags and 'no_estimate' in z.call_flags)
    _ok('the call table has no column called "flags" (it would collide with the pandas attribute)', 'flags' not in t.columns and 'call_flags' in t.columns)

    # ======================================================================================== calibration from anchors
    ratios = [0.50 + 0.25 * (n - 2) + e for n, e in zip([2, 2, 2, 3, 3, 3, 4, 4], [0.01, -0.01, 0.0, 0.01, -0.01, 0.0, 0.0, 0.0])]
    cal2 = C.fit_anchors(ratios, [2, 2, 2, 3, 3, 3, 4, 4])
    _ok('fit_anchors recovers the line through samples of known copy number (here 0.50 and 0.25)', abs(cal2['ratio_at_two_copies'] - 0.5) < 0.02 and abs(cal2['ratio_per_copy'] - 0.25) < 0.02, str(cal2))
    _ok('fit_anchors needs two different copy numbers', _raises(lambda: C.fit_anchors([0.5, 0.51], [2, 2]), ValueError))
    t2 = C.call_table([raw_of('th3', 0.75, 0.5)], cfg, calibration=C.THEORY_CALIBRATION).iloc[0]
    _ok('a calibration passed to call_table replaces the one in the config (theory: ratio 0.75 = 3 copies)', t2.gsdf_number == 3 and t2.call_flags == '')

    # ======================================================================================== measuring a BAM (fake samtools)
    log = []
    def fake_run(ratio, jr, contig='NC_135173.1'):
        def run_(cmd):
            log.append(cmd)
            if cmd[:3] == ['samtools', 'view', '-H']:
                return types.SimpleNamespace(stdout=f'@HD\tVN:1.6\n@SQ\tSN:{contig}\tLN:68000000\n')
            d, q = depth_arrays(ratio, jr)
            arr_ = q if '-Q' in cmd else d
            return types.SimpleNamespace(stdout=''.join(f'{r["chrom"]}\t{r["start"] + 1 + i}\t{arr_[i]:g}\n' for i in range(L)))
        return run_
    m = M.measure_bam('x.bam', cfg, fake_run(0.798, 0.5))
    mm = M.metrics(m, cfg); cc = C.call_one(mm, cfg)
    _ok('measure_bam: two depth passes (all reads, MAPQ >= 30) and a call of 3 copies', sum(c[1] == 'depth' for c in log) == 2 and any('-Q' in c and '30' in c for c in log) and cc['gsdf_number'] == 3 and m['definition'] == k)
    _ok('measure_bam: a BAM without the contig is refused (wrong reference)', _raises(lambda: M.measure_bam('x.bam', cfg, fake_run(0.5, 0, contig='chr7')), ValueError))

    # ======================================================================================== panel + command line, through a real subprocess with a fake samtools
    tmp = tempfile.mkdtemp() + '/'
    os.makedirs(tmp + 'bin'); open(tmp + 'bin/samtools', 'w').write(FAKE_SAMTOOLS); os.chmod(tmp + 'bin/samtools', 0o755)
    from cichlid_common.common_light import module_dirs, REPO_ROOT
    env = dict(os.environ, PATH=tmp + 'bin:' + os.environ['PATH'], PYTHONPATH=os.pathsep.join(module_dirs()))
    pkg = REPO_ROOT
    def cli(*args, **e):
        return subprocess.run([sys.executable, '-m', 'sd_genotypes'] + list(args), capture_output=True, text=True, cwd=pkg, env={**env, **{k: str(v) for k, v in e.items()}})
    open(tmp + 'S1.all.bam', 'wb').write(b'bam')
    r_ = cli('panel', 'measure', '--bam', tmp + 'S1.all.bam', '--sample', 'S1', '--out-json', tmp + 'm/S1.panel.json', FAKE_RATIO=0.53, FAKE_JUNCTION=0)
    _ok('panel measure refuses a BAM without an index', r_.returncode != 0 and 'no index' in (r_.stderr + r_.stdout))
    open(tmp + 'S1.all.bam.bai', 'wb').write(b'i')
    r_ = cli('panel', 'measure', '--bam', tmp + 'S1.all.bam', '--sample', 'S1', '--out-json', tmp + 'm/S1.panel.json', FAKE_RATIO=0.53, FAKE_JUNCTION=0)
    _ok('panel measure (command line, real subprocess) writes the JSON', r_.returncode == 0 and os.path.exists(tmp + 'm/S1.panel.json'), r_.stderr[-300:])
    j = json.load(open(tmp + 'm/S1.panel.json'))
    _ok('JSON: sample, BAM size, panel version, samtools version, and the gsdf module with its definition key', j['sample'] == 'S1' and j['bam_size'] == 3 and j['panel_version'].startswith('v') and 'samtools 1.19' in j['samtools']
        and j['modules']['gsdf']['definition'] == k and len(j['modules']['gsdf']['bins']) == 26)
    r_ = cli('panel', 'measure', '--bam', tmp + 'S1.all.bam', '--sample', 'S1', '--out-json', tmp + 'm/x.json', '--modules', 'lg5', FAKE_RATIO=0.53, FAKE_JUNCTION=0)
    _ok('an unknown module is refused with the list of the available ones', r_.returncode != 0 and 'gsdf' in (r_.stderr + r_.stdout))
    r_ = cli('panel', 'measure', '--bam', tmp + 'S1.all.bam', '--sample', 'S9', '--out-json', tmp + 'm/S9.panel.json', FAKE_RATIO=0.53, FAKE_JUNCTION=0, FAKE_CONTIG='chr7')
    _ok('a BAM of the wrong reference is refused', r_.returncode != 0 and 'NC_135173.1' in (r_.stderr + r_.stdout))
    rows, skipped = PM.collect(tmp + 'm/')
    _ok('collect reads a folder of panel JSONs', len(rows) == 1 and rows[0]['SampleID'] == 'S1' and skipped == [])
    PM.write_json(tmp + 'm/other.panel.json', dict(schema=1, sample='O', bam='x', bam_size=1, panel_version='v', measured='t', samtools='s', modules={}))
    rows, skipped = PM.collect(tmp + 'm/')
    _ok('a panel file without a gsdf measurement is skipped, not an error', len(rows) == 1 and skipped == ['other.panel.json'])
    PM.write_json(tmp + 'm/bad.panel.json', dict(schema=99))
    _ok('a panel file of another schema is refused', _raises(lambda: PM.collect(tmp + 'm/'), ValueError)); os.remove(tmp + 'm/bad.panel.json')

    # ======================================================================================== a synthetic cohort: females 2 copies, males 3 (XY-like species), plus noise cases
    rng = np.random.default_rng(5)
    os.makedirs(tmp + 'cohort')
    sheet = []
    def add(sid, species, sex, ratio, jr, flank=15.0, eco='Shallow_Benthic', sheet_value=None):
        x = raw_of(sid, ratio + rng.normal(0, 0.02), jr, flank + rng.normal(0, 0.5))
        PM.write_json(f'{tmp}cohort/{sid}.panel.json', dict(schema=1, sample=sid, bam='x', bam_size=x.pop('bam_size'), panel_version='v34', measured=x.pop('measured'), samtools='s', modules={'gsdf': {k_: v for k_, v in x.items() if k_ not in ('SampleID', 'panel_version')}}))
        sheet.append(dict(SampleID=sid, Species=species, Sex=sex, Ecogroup=eco, **{'gsdf number': sheet_value}))
    for i in range(12):
        add(f'XYF{i}', 'Xy species', 'F', 0.530, 0.0, sheet_value=2)
        add(f'XYM{i}', 'Xy species', 'M', 0.798, 0.5, sheet_value=3)
    for i in range(6):
        add(f'NDF{i}', 'No dup species', 'F', 0.530, 0.0, sheet_value=2); add(f'NDM{i}', 'No dup species', 'M', 0.530, 0.0, sheet_value=2)
    for i in range(6):                                                                   # duplication present, equally in both sexes
        add(f'DPF{i}', 'Dup both sexes', 'F', 0.798 if i % 2 else 0.530, 0.5 if i % 2 else 0.0, sheet_value=3 if i % 2 else 2)
        add(f'DPM{i}', 'Dup both sexes', 'M', 0.798 if i % 2 else 0.530, 0.5 if i % 2 else 0.0, sheet_value=3 if i % 2 else 2)
    add('XYF_exc', 'Xy species', 'F', 0.798, 0.5, sheet_value=3)                           # a female with the male-specific copy
    add('XYM_four', 'Xy species', 'M', 1.066, 1.0, sheet_value=4)                              # a 4-copy male
    add('LOW1', 'Xy species', 'M', 0.798, 0.5, flank=2.0, sheet_value=3)                   # low depth: flagged
    add('NEWM', 'New species', 'M', 0.798, 0.5, sheet_value=np.nan)                        # not in the sheet yet
    pd.DataFrame(sheet).to_csv(tmp + 'sheet.tsv', sep='\t', index=False)
    r_ = cli('gsdf', 'call', '--measurements', tmp + 'cohort', '--out-dir', tmp + 'calls')
    _ok('gsdf call (command line, offline) runs', r_.returncode == 0 and os.path.exists(tmp + 'calls/gsdf_calls.tsv'), r_.stderr[-400:])
    calls = pd.read_csv(tmp + 'calls/gsdf_calls.tsv', sep='\t', keep_default_na=False, na_values=['']).fillna({'call_flags': ''})
    truth = {row['SampleID']: row['gsdf number'] for row in sheet if row['gsdf number'] == row['gsdf number']}
    good = calls[(calls.call_flags == '') & calls.SampleID.isin(truth)]
    _ok('every unflagged call equals the planted copy number', len(good) > 40 and all(int(r.gsdf_number) == truth[r.SampleID] for r in good.itertuples()), f'{len(good)} unflagged')
    _ok('the low-depth sample is flagged, the 4-copy male is called 4', calls.set_index('SampleID').loc['LOW1'].call_flags == 'low_flank_depth' and calls.set_index('SampleID').loc['XYM_four'].gsdf_number == 4)
    txt = open(tmp + 'calls/summary.txt').read()
    _ok('summary: counts per copy number, flags, and the drift check against the calibration', 'copies (rounded)' in txt and 'low_flank_depth: 1' in txt and 'expected_ratio' in txt, txt.split('\n')[1])
    cl = C.cluster_check(calls, cfg).set_index('copies')
    _ok('drift check: unflagged 2- and 3-copy samples sit at the calibrated ratios (0.530 / 0.798)', abs(cl.loc[2].median_ratio - 0.530) < 0.02 and abs(cl.loc[3].median_ratio - 0.798) < 0.02)
    su = pd.read_csv(tmp + 'calls/sheet_update.tsv', sep='\t')
    _ok('sheet_update.tsv has the columns to paste into the sheet (`gsdf number`, `gsdf flags`)', list(su.columns) == ['SampleID', 'gsdf number', 'gsdf flags'] and len(su) == len(calls))
    PM.write_json(tmp + 'cohort/DUP.panel.json', json.load(open(tmp + 'cohort/XYF0.panel.json')) | dict(sample='XYF1'))
    _ok('the same sample measured twice is refused', _raises(lambda: S.call_step(tmp + 'cohort', out_dir=tmp + 'c2'), SystemExit)); os.remove(tmp + 'cohort/DUP.panel.json')
    old = json.load(open(tmp + 'cohort/XYF0.panel.json')); old['sample'] = 'OLD'; old['modules']['gsdf']['definition'] = 'abc123'
    PM.write_json(tmp + 'cohort/OLD.panel.json', old)
    _ok('measurements made under another definition are refused (not silently mixed)', _raises(lambda: S.call_step(tmp + 'cohort', out_dir=tmp + 'c3'), SystemExit)); os.remove(tmp + 'cohort/OLD.panel.json')

    # ======================================================================================== validation
    sh = V.sheet_frame(pd.DataFrame(dict(SampleID=['a', 'b', 'c'], Sex=['M', 'F?', 'U'], Species=['s'] * 3, Ecogroup=['e'] * 3, **{'gsdf number': ['3', '2', '']})))
    _ok('sheet frame: uncertain and unknown sexes become missing, values numeric', sh.Sex.tolist()[0] == 'M' and sh.Sex.isna().tolist()[1:] == [True, True] and sh.sheet_gsdf.iloc[0] == 3 and np.isnan(sh.sheet_gsdf.iloc[2]))
    sheet_df = pd.read_csv(tmp + 'sheet.tsv', sep='\t')
    summ, mism = V.compare_with_sheet(calls, V.sheet_frame(sheet_df))
    _ok('agreement with the existing column: all unflagged calls agree; the flagged and the new sample are left out', summ['frac_agree'] == 1.0 and summ['n_flagged_among_with_sheet'] == 1 and len(mism) == 0, str(summ))
    bad = calls.copy(); bad.loc[bad.SampleID == 'XYM0', 'gsdf_number'] = 2
    s2, m2 = V.compare_with_sheet(bad, V.sheet_frame(sheet_df))
    _ok('a disagreement is listed', len(m2) == 1 and m2.iloc[0].SampleID == 'XYM0' and s2['n_agree'] == s2['n_compared'] - 1)
    use = calls[calls.call_flags == ''].merge(V.sheet_frame(sheet_df)[['SampleID', 'Sex', 'Species', 'Ecogroup']], on='SampleID').rename(columns={'gsdf_number': 'copies'})
    tab, disc = V.sex_table(use[['SampleID', 'Sex', 'Species', 'Ecogroup', 'copies']], 3, 3)
    pat = dict(zip(tab.Species, tab.pattern))
    _ok('sex x species: male-specific, no duplication, and duplication in both sexes are told apart',
        pat == {'Xy species': 'male-specific (XY-like)', 'No dup species': 'no duplication seen', 'Dup both sexes': 'duplication, not male-specific'}, str(pat))
    _ok('the exceptions of the male-specific species are listed (the female with the male-specific copy; a 4-copy male is not an exception)', disc.SampleID.tolist() == ['XYF_exc'] and disc.why.iloc[0].startswith('female'), str(disc.SampleID.tolist()))
    _ok('a species with fewer than 3 sexed males or females is not tabulated', 'New species' not in pat)
    r_ = cli('gsdf', 'validate', '--table', tmp + 'sheet.tsv', '--calls', tmp + 'calls/gsdf_calls.tsv', '--out-dir', tmp + 'val')
    _ok('gsdf validate (command line, offline) writes the agreement, the species table and the exceptions', r_.returncode == 0 and all(os.path.exists(tmp + 'val/' + f) for f in ('sex_by_species.tsv', 'sex_exceptions.tsv', 'mismatches_vs_sheet.tsv', 'summary.txt')), r_.stderr[-300:])
    r_ = cli('gsdf', 'validate', '--table', tmp + 'sheet.tsv', '--out-dir', tmp + 'val2')
    _ok('validate without calls tabulates the existing values of the sheet', r_.returncode == 0 and 'EXISTING' in r_.stdout and 'male-specific (XY-like)' in r_.stdout)

    # ======================================================================================== run-folder mode (with a Project)
    from selftest_bchrom import make_proj
    proj = make_proj(tempfile.mkdtemp() + '/')
    S.call_step(tmp + 'cohort', proj=proj)
    rid = open(proj.root + 'gsdf/call/LATEST').read().strip()
    _ok('with a Project the calls go to an immutable run folder with a manifest', os.path.exists(f'{proj.root}gsdf/call/{rid}/gsdf_calls.tsv') and os.path.exists(f'{proj.root}gsdf/call/{rid}/manifest.json'))
    print('gsdf / panel selftests passed')
