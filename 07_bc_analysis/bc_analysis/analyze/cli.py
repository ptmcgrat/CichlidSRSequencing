"""python -m bc_analysis analyze marey|rate|sexscan|finemap|crossovers|trd|linkscan  [--genome G] [--run hmm_mask_int2_deep] [--fish-table T.tsv] [--upload] ...

   marey                       genetic against physical position of every chromosome (Marey maps), map lengths, rates in windows
   rate --lg 10 --region 27.5-30.0   is the recombination rate of a region unusual: against its chromosome (bootstrap over fish), against same-width windows elsewhere, the hottest window of every chromosome
   sexscan                     the association of every bin's state with the recorded sex, along the genome
   finemap [--lg 10]           the limits of the sex determiner from the fish with a crossover near it (a table of the recombinants, the limits if 0-3 nearest fish on each side were wrong)
   crossovers [--from-finemap | --fish F:left,F:right ...] [--counts DIR] [--counts-before DIR]   the bounding crossovers localized from the READS, before and after more sequencing, and the interval they give
   trd                         transmission ratio distortion with a sign-flip permutation null (genome-wide p), its regions and the checks by sire, brood, dam and sex
   linkscan --lg 3 --region 8-11   which part of the genome a segment is inherited with (finds a segment that an assembly attached to the wrong chromosome)
   publish                     upload the WHOLE genome folder (sites, counts, HMM runs, mask, analyses) to its Dropbox twin, verified: the record of the BC genotyping must not live only on one machine

Reads <genome folder>/<run>/genotypes_calls.csv (default run: config analysis.run), the mask <genome folder>/mask/mask_bins.tsv, the sample table for sex, sire, dam and brood (the lab sheets, or --fish-table).
Writes <genome folder>/analysis/<step>/ : tables (.tsv), figures (.png, when matplotlib is there) and provenance.json (the command, the inputs with their sha1, versions); --upload copies that folder to Dropbox
(the same path under the cloud master, <cloud master>/Outputs/BCGenotypes/<genome>/analysis/<step>/). Inputs that are not on this machine are fetched from Dropbox with the FileManager (unless --no-download). <genome folder> = <local master>/Outputs/BCGenotypes/<genome>/ unless --out is given."""
import argparse
import datetime
import json
import os
import platform
import socket
import subprocess
import sys

import numpy as np
import pandas as pd

from . import crossovers as X
from . import data as D
from . import distortion as T
from . import plots as P
from . import recomb as R
from . import sexlocus as SL

STEPS = ['marey', 'rate', 'sexscan', 'finemap', 'crossovers', 'trd', 'linkscan', 'publish']
DEFAULTS = dict(run='hmm_mask_int2_deep', terminal_mb=2.0, marey_window_mb=2.5, bootstrap=5000, null_skip_lg=[3], rate_step_mb=0.5, finemap_window_mb=5.0, finemap_flank_mb=1.5, crossover_flank_mb=1.5, crossovers_per_side=2,
                trd_permutations=2000, trd_zmin=3.0, linkscan_smooth=10)


def build_parser():
    ap = argparse.ArgumentParser(prog='python -m bc_analysis analyze', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('step', choices=STEPS)
    ap.add_argument('--genome', default='Mconophoros_GT1_NCBI')
    ap.add_argument('--config', default=None)
    ap.add_argument('--out', default=None, help='the genome folder of the genotype outputs (default <local master>/Outputs/BCGenotypes/<genome>/)')
    ap.add_argument('--run', default=None, help='the HMM run folder to analyse (default: config analysis.run)')
    ap.add_argument('--calls', default=None, help='a genotypes_calls.csv (default <genome folder>/<run>/genotypes_calls.csv)')
    ap.add_argument('--mask', default=None, help="the mask (default <genome folder>/mask/mask_bins.tsv if it exists; 'none' for no mask)")
    ap.add_argument('--fish-table', default=None, help='a tab-separated table with SampleID, Sex (M/F), and Sire, Dam, BroodID (default: from the lab sheets)')
    ap.add_argument('--analysis-dir', default=None, help='where to write (default <genome folder>/analysis/<step>/)')
    ap.add_argument('--upload', action='store_true', help='copy the output folder to its Dropbox twin and verify it (the lab FileManager paths, cichlid_common.runner.sync)')
    ap.add_argument('--no-download', action='store_true', help='do not fetch missing inputs (the calls, the mask, counts files) from Dropbox')
    ap.add_argument('--no-plots', action='store_true')
    ap.add_argument('--lg', default=None, help='a linkage group: 10 or LG10')
    ap.add_argument('--region', default=None, help='START-END in Mb')
    ap.add_argument('--locus-mb', type=float, default=None, help='finemap: where to look (default: the best bin of the sex scan)')
    ap.add_argument('--exclude-fish', default='', help='finemap: fish to ignore, comma-separated')
    ap.add_argument('--counts', default=None, help='crossovers: the counts folder (default <genome folder>/counts/)')
    ap.add_argument('--counts-before', default=None, help='crossovers: a folder with the counts before more sequencing, for the comparison')
    ap.add_argument('--fish', default=None, help='crossovers: SampleID:left or SampleID:right, comma-separated')
    ap.add_argument('--from-finemap', action='store_true', help='crossovers: take the nearest crossovers of the fine-mapping on each side')
    ap.add_argument('--span', default=None, help='crossovers: START-END in Mb of the window with the reads (default: the locus +- config analysis.crossover_flank_mb)')
    ap.add_argument('--group-by', default='Sire,BroodID,Dam', help='the columns of the fish table to summarise by')
    ap.add_argument('--seed', type=int, default=1)
    o = ap.add_argument_group('offline (tests, no lab sheets)')
    o.add_argument('--sample-table', default=None)
    o.add_argument('--alignment-table', default=None)
    o.add_argument('--local-master', default=None)
    o.add_argument('--cloud-master', default=None)
    o.add_argument('--fasta', default=None)
    return ap


def make_env(a):
    """Local and Dropbox master folders: the lab file manager (live), or --local-master and --cloud-master (offline, tests)."""
    from variant_calling import build as B
    if a.local_master or a.cloud_master:
        if not (a.local_master and a.cloud_master):
            raise ValueError('--local-master and --cloud-master go together')
        return B.Env(a.local_master, a.cloud_master, a.genome, a.fasta or '')
    return B.Env.from_file_manager(a.genome)


class Ctx:
    def __init__(self, a, cfg, out, uploader=None, downloader=None):
        from bc_analysis.genotyper import cli as G
        self.a, self.cfg, self.out, self.uploader, self.downloader = a, cfg, out, uploader, downloader
        self.G = G
        self.an = {**DEFAULTS, **(cfg.get('analysis') or {})}
        G.set_lg_names(cfg, a.genome)
        self._env = None
        if a.out:
            self.root = a.out.rstrip('/') + '/'
        elif a.local_master:
            self.root = f"{a.local_master.rstrip('/')}/Outputs/BCGenotypes/{a.genome}/"
        else:
            self.root = G.genome_dir(a)
        self.run = a.run or self.an['run']
        self.dir = (a.analysis_dir.rstrip('/') + '/') if a.analysis_dir else f'{self.root}analysis/{a.step}/'
        if a.step != 'publish':
            os.makedirs(self.dir, exist_ok=True)
        self.inputs = {}
        self.written = []

    def env(self):
        if self._env is None:
            self._env = make_env(self.a)
        return self._env

    def need(self, path):
        """True when `path` is on this machine, fetching it from Dropbox with the FileManager first when it is not (unless --no-download). A path outside the local master has no Dropbox twin.
        Never raises: a failed fetch leaves its reason in self.fetch_error for the message of the caller."""
        self.fetch_error = None
        path = clean_path(path)
        if os.path.exists(path):
            return True
        if self.a.no_download:
            return False
        try:
            env = self.env()
            if not os.path.abspath(path).startswith(os.path.abspath(env.local_master)):
                return False
            self.out(f'fetching {path} from Dropbox')
            dl = self.downloader
            if dl is None:
                from cichlid_common.file_manager import FileManager
                dl = FileManager().downloadData
            dl(path)
        except Exception as e:                                                  # noqa: the file is not in Dropbox either, or the file manager is not there
            self.fetch_error = f'{type(e).__name__}: {str(e)[:120]}'
            return False
        return os.path.exists(path)

    def why(self):
        return f'; it could not be fetched from Dropbox either ({self.fetch_error})' if self.fetch_error else ''

    def calls(self):
        path = self.a.calls or f'{self.root}{self.run}/genotypes_calls.csv'
        if not self.need(path):
            raise ValueError(f'no calls at {path}: run `genotype hmm` first (or give --calls / --run)' + self.why())
        self.inputs['calls'] = path
        return D.Calls.load(path, self.G.lg_label)

    def mask(self):
        p = self.a.mask
        if p == 'none':
            return set()
        p = p or f'{self.root}mask/mask_bins.tsv'
        if self.need(p):
            self.inputs['mask'] = p
        else:
            self.out(f'note: no mask at {p}: none is applied')
        return D.read_mask(p)

    def fish(self, need=True):
        if self.a.fish_table:
            self.inputs['fish_table'] = self.a.fish_table
            return D.fish_table(pd.read_csv(self.a.fish_table, sep='\t', dtype=str, keep_default_na=False))
        try:
            from variant_calling import cohort_cli
            sample_dt, _ = cohort_cli.load_tables(self.a, self.a.genome)
            bc = self.G.bc_fish(sample_dt, self.cfg)
            sex = sample_dt.drop_duplicates('SampleID')[['SampleID', 'Sex']] if 'Sex' in sample_dt.columns else bc[['SampleID']].assign(Sex='')
            return D.fish_table(bc.merge(sex, on='SampleID', how='left'))
        except (SystemExit, Exception) as e:                                   # noqa: the sheets cannot be read
            if need:
                raise ValueError(f'the fish table is needed (sex, sire, dam, brood) and the lab sheets cannot be read ({type(e).__name__}: {str(e)[:100]}): give --fish-table') from e
            return None

    def write(self, name, df):
        path = f'{self.dir}{name}'
        df.to_csv(path, sep='\t', index=False)
        self.written.append(path)
        return path

    def figure(self, name, fn, *args, **kw):
        if self.a.no_plots:
            return None
        path = f'{self.dir}{name}'
        try:
            fn(*args, path, **kw)
        except ImportError:
            self.out('note: matplotlib is not available here, so no figure was drawn (the tables are written)')
            return None
        self.written.append(path)
        return path


def mb_pair(text, what='--region'):
    try:
        a, b = (float(x) for x in str(text).split('-'))
    except (ValueError, TypeError):
        raise ValueError(f'{what} must be START-END in Mb, e.g. 27.5-30.0') from None
    if not a < b:
        raise ValueError(f'{what}: START must be below END')
    return int(round(a * 1e6)), int(round(b * 1e6))


# ------------------------------------------------------------------ the steps
def step_marey(x):
    calls, masked = x.calls(), x.mask()
    steps = R.Steps(calls)
    tip = int(x.an['terminal_mb'] * 1e6)
    mb = R.marey_bins(steps)
    lens = R.map_lengths(steps, tip, masked)
    win = R.windows_table(steps, int(x.an['marey_window_mb'] * 1e6), tip, masked)
    x.write('marey_bins.tsv', mb.round(4)); x.write('map_lengths.tsv', lens); x.write('marey_windows.tsv', win.round(3))
    x.figure('marey_maps.png', P.plot_marey, mb, lens, masked, tip, bin_bp=calls.bin_bp)
    x.out(f'{len(lens)} chromosomes, {len(calls.fish)} fish: mean map {lens.cM.mean():.1f} cM per chromosome ({lens.cM_without_tips.mean():.1f} without the first and last {x.an["terminal_mb"]:g} Mb); '
          f'longest {lens.sort_values("cM").iloc[-1].LG} {lens.cM.max():.0f} cM, shortest {lens.sort_values("cM").iloc[0].LG} {lens.cM.min():.0f} cM')
    return calls


def step_rate(x):
    a = x.a
    if not a.lg or not a.region:
        raise ValueError('rate needs --lg and --region START-END (Mb), e.g. --lg 10 --region 27.5-30.0')
    calls, masked = x.calls(), x.mask()
    steps = R.Steps(calls)
    chrom = calls.chrom_of(a.lg)
    lo, hi = mb_pair(a.region)
    tip = int(x.an['terminal_mb'] * 1e6)
    skip = {f'LG{g}' for g in x.an['null_skip_lg']} | {chrom}
    null = R.window_rates(steps, hi - lo, int(x.an['rate_step_mb'] * 1e6), tip, masked, skip)
    res = R.region_test(steps, chrom, lo, hi, tip, int(x.an['bootstrap']), a.seed, masked, null)
    best = R.best_windows(steps, hi - lo, int(x.an['rate_step_mb'] * 1e6), tip, masked, {f'LG{g}' for g in x.an['null_skip_lg']})
    res['rank_of_chromosome_best_window'] = int(best[best.chrom == chrom]['rank'].iloc[0]) if len(best) and (best.chrom == chrom).any() else None
    res['chromosomes_ranked'] = int(len(best))
    prof = R.rate_profile(steps, chrom, 1_000_000, 2000, a.seed)
    x.write('rate_region.tsv', pd.DataFrame([res]).round(4)); x.write('rate_best_windows.tsv', best.round(3)); x.write('rate_profile.tsv', prof.round(3))
    ft = x.fish(need=False)
    if ft is not None:
        rows = []
        for col in [c for c in a.group_by.split(',') if c in ft.columns]:
            t = R.by_group(steps, chrom, lo, hi, D.groups_of(calls, ft, col)); t.insert(0, 'by', col); rows.append(t)
        if rows:
            x.write('rate_by_group.tsv', pd.concat(rows, ignore_index=True).round(3))
    x.figure('rate.png', P.plot_rate, prof, res, null.rate.to_numpy() if len(null) else None)
    x.out(f"{res['LG']} {a.region} Mb: {res['rate_cM_per_Mb']:.2f} cM/Mb against {res['background_cM_per_Mb']:.2f} for the rest of the chromosome (ratio {res['ratio']:.2f}, 95% CI {res['ci_lo']:.2f}-{res['ci_hi']:.2f}); "
          + (f"{res['null_as_hot_or_hotter']} of {res['null_windows']} same-width windows elsewhere are as hot (p = {res['p_genome']:.4f}); " if 'p_genome' in res else '') + f"rank {res['rank_of_chromosome_best_window']} of {res['chromosomes_ranked']} chromosomes by their hottest window")
    return calls


def _sex(x, calls):
    ft = x.fish()
    return ft, D.sex_vector(calls, ft)


def step_sexscan(x):
    calls = x.calls()
    ft, sex = _sex(x, calls)
    scan = SL.sexscan(calls, sex)
    summ = SL.scan_summary(scan)
    x.write('sexscan_bins.tsv', scan.round(3)); x.write('sexscan_summary.tsv', summ.round(3))
    x.figure('sexscan.png', P.plot_sexscan, scan)
    b = summ.iloc[0]
    x.out(f'{int((~np.isnan(sex)).sum())} sexed fish. The strongest association is {b.LG} at {b.best_Mb:.1f} Mb (-log10 p = {b.nlogp:.1f}; {b.M_with_YH} males with the YH state, {b.F_with_MC} females with the MC state); '
          f'the next chromosome reaches {summ.iloc[1].nlogp:.1f}' if len(summ) > 1 else '')
    return calls


def do_finemap(x, calls, sex):
    a = x.a
    if a.lg:
        chrom = calls.chrom_of(a.lg)
    else:
        scan = SL.sexscan(calls, sex)
        chrom = SL.scan_summary(scan).iloc[0].chrom
    res = SL.finemap(calls, sex, chrom, int(a.locus_mb * 1e6) if a.locus_mb else None, int(x.an['finemap_window_mb'] * 1e6), [s for s in a.exclude_fish.split(',') if s])
    return chrom, res


def step_finemap(x):
    calls = x.calls()
    ft, sex = _sex(x, calls)
    chrom, res = do_finemap(x, calls, sex)
    rec, lim = res['recombinants'], res['limits'].copy()
    lim['left_Mb'], lim['right_Mb'], lim['width_kb'] = lim.left_bp / 1e6, lim.right_bp / 1e6, lim.width_bp / 1e3
    x.write('finemap_recombinants.tsv', rec); x.write('finemap_limits.tsv', lim.round(4)); x.write('finemap_discordance.tsv', res['profile'])
    pd.DataFrame({'fish': res['excluded']}).assign(reason='discordant at the best bin or --exclude-fish').to_csv(f'{x.dir}finemap_excluded.tsv', sep='\t', index=False); x.written.append(f'{x.dir}finemap_excluded.tsv')
    x.figure('finemap.png', P.plot_finemap, res, calls, chrom, sex, flank_bp=int(x.an['finemap_flank_mb'] * 1e6))
    if len(lim):
        l0 = lim.iloc[0]
        x.out(f"{calls.lg_of(chrom)}: the determiner lies between {l0.left_Mb:.3f} and {l0.right_Mb:.3f} Mb ({l0.width_kb:.0f} kb), best bin {res['locus_bp'] / 1e6:.1f} Mb; "
              f"if one nearest fish on each side were wrong: {lim.iloc[1].left_Mb:.3f}-{lim.iloc[1].right_Mb:.3f} Mb" if len(lim) > 1 else f"{calls.lg_of(chrom)}: limits {l0.left_Mb:.3f}-{l0.right_Mb:.3f} Mb")
    else:
        x.out(f'{calls.lg_of(chrom)}: no crossover on both sides of the best bin within {x.an["finemap_window_mb"]:g} Mb')
    x.out(f"ignored fish: {', '.join(res['excluded']) or 'none'}" + (f"; of these {len(res['discordant'])} disagree with their recorded sex at the locus (mis-sexed, or a second factor)" if res['discordant'] else ''))
    return calls


def step_crossovers(x):
    a = x.a
    calls = x.calls()
    counts = (a.counts or f'{x.root}counts/').rstrip('/') + '/'
    roles, chrom, locus = {}, None, None
    if a.fish:
        roles = {f.split(':')[0]: f.split(':')[1] for f in a.fish.split(',')}
        if not a.lg or not a.span:
            raise ValueError('crossovers --fish needs --lg and --span START-END (Mb)')
        chrom = calls.chrom_of(a.lg)
        lo, hi = mb_pair(a.span, '--span')
    elif a.from_finemap:
        ft, sex = _sex(x, calls)
        chrom, res = do_finemap(x, calls, sex)
        roles = X.roles_from_finemap(res['recombinants'], int(x.an['crossovers_per_side']))
        locus = res['locus_bp']
        flank = int(x.an['crossover_flank_mb'] * 1e6)
        lo, hi = (mb_pair(a.span, '--span') if a.span else (max(locus - flank, 0), locus + flank))
    else:
        raise ValueError('crossovers needs --from-finemap or --fish SampleID:left,SampleID:right,...')
    if not any(r == 'left' for r in roles.values()) or not any(r == 'right' for r in roles.values()):
        raise ValueError(f'the crossovers need a left and a right fish; got {roles}')
    x.inputs['counts'] = counts
    rows, tab, profs = [], {}, {}
    sets = [('after', counts)] + ([('before', a.counts_before.rstrip('/') + '/')] if a.counts_before else [])
    for tag, folder in sets:
        res_t = {}
        for f in roles:
            if not x.need(f'{folder}{f}.counts.tsv.gz'):
                raise ValueError(f'no counts file for {f} in {folder}' + (' (the folder does not exist here: was it copied into the genome folder?)' if not os.path.isdir(folder) else '') + x.why())
            pos, n, y = X.read_sites(folder, f, chrom, lo, hi)
            loc = X.localize(pos, n, y)
            e, q = X.fitted_rates(pos, n, y, loc['split_index'], loc['state'])
            fit = X.localize(pos, n, y, e, q)
            sw = X.sweep(pos, n, y)
            tab[(tag, f)] = dict(fish=f, role=roles[f], when=tag, reads=loc['reads'], sites=loc['sites'], best_bp=loc['best'], gap_bp=loc['gap'], span1_lo=loc['span1'][0], span1_hi=loc['span1'][1], span2_lo=loc['span2'][0], span2_hi=loc['span2'][1],
                                 state=loc['state'], edge=loc['edge'], fitted_eps=round(e, 4), fitted_q=round(q, 3), fitted_best_bp=fit['best'], sweep_min_bp=min(sw.values()), sweep_max_bp=max(sw.values()))
            res_t[f] = loc
            profs[(tag, f)] = (pos, n, y, e, q)
        best, cons = X.interval(res_t, roles)
        rows.append(dict(when=tag, best_left_bp=best[0], best_right_bp=best[1], best_width_bp=best[1] - best[0], conservative_left_bp=cons[0], conservative_right_bp=cons[1], conservative_width_bp=cons[1] - cons[0]))
    t = pd.DataFrame(list(tab.values())); iv = pd.DataFrame(rows)
    x.write('crossovers.tsv', t.round(3)); x.write('crossovers_interval.tsv', iv.round(1))
    def lod(tag, f):
        pos, n, y, e, q = profs[(tag, f)]
        l0, l1 = X._ll(n, y, e), X._ll(n, y, q); best = None
        for s, la, lb in ((0, l0, l1), (1, l1, l0)):
            tot = np.concatenate([[0], np.cumsum(la)]) + np.concatenate([np.cumsum(lb[::-1])[::-1], [0]])
            if best is None or tot.max() > best[0].max():
                best = (tot, s)
        tot = best[0]
        xs = np.array([(pos[max(k - 1, 0)] + pos[min(k, len(pos) - 1)]) / 2 for k in range(len(pos) + 1)])
        return xs, (tot - tot.max()) / np.log(10)
    ev = []
    lods = {}
    for f in roles:
        row = tab[('after', f)]
        xa, la = lod('after', f)
        xb, lb = lod('before', f) if ('before', f) in profs else (None, None)
        ev.append(dict(fish=f, role=roles[f], best=row['best_bp'], best_before=tab[('before', f)]['best_bp'] if ('before', f) in tab else None)); lods[f] = (xa, la, xb, lb)
    ev.sort(key=lambda r: r['best'])
    final = iv[iv.when == 'after'].iloc[0]
    x.figure('crossovers.png', P.plot_crossovers, ev, lods, (final.best_left_bp, final.best_right_bp))
    for _, r in iv.iterrows():
        x.out(f"{r['when']:6s}: the determiner lies between {r.best_left_bp / 1e6:.3f} and {r.best_right_bp / 1e6:.3f} Mb ({r.best_width_bp / 1e3:.0f} kb; conservative {r.conservative_left_bp / 1e6:.3f}-{r.conservative_right_bp / 1e6:.3f} Mb)")
    for _, r in t[t.when == 'after'].iterrows():
        x.out(f"  {r.fish} ({r.role}): crossover {r.best_bp / 1e6:.3f} Mb, 1-LOD span {(r.span1_hi - r.span1_lo) / 1e3:.0f} kb, {r.reads} reads at {r.sites} sites" + ('  EDGE: it may lie outside the window' if r.edge else ''))
    return calls


def step_trd(x):
    calls = x.calls()
    sc = T.trd_scan(calls, int(x.an['trd_permutations']), x.a.seed)
    reg = T.trd_regions(calls, sc, x.an['trd_zmin'])
    tip = int(x.an['terminal_mb'] * 1e6)
    reg['in_tip'] = [bool(r.start < tip or r.end > calls.length(r.chrom) - tip) for r in reg.itertuples()]                      # the first and last terminal_mb Mb: low confidence (assembly ends)
    g95, g99 = np.percentile(sc['G'], [95, 99])
    bins = calls.meta[['chrom', 'LG', 'pos']].assign(YH_share=sc['f'], z=sc['z'])
    x.write('trd_bins.tsv', bins.round(4)); x.write('trd_regions.tsv', reg.round(4))
    ft = x.fish(need=False)
    if ft is not None and len(reg):
        rows = []
        for r in reg.itertuples():
            st = T.region_state(calls, r.chrom, r.start, r.end)
            for col in [c for c in (x.a.group_by.split(',') + ['Sex']) if c in ft.columns]:
                t = T.group_checks(st, D.groups_of(calls, ft, col), col); t.insert(0, 'region', f'{r.LG}:{r.start / 1e6:.1f}-{r.end / 1e6:.1f}'); rows.append(t)
        if rows:
            x.write('trd_groups.tsv', pd.concat(rows, ignore_index=True).round(4))
    x.figure('trd_genomewide.png', P.plot_trd, calls, sc, g95 / (2 * np.sqrt(sc['n'])))
    sig = reg[reg.p_genomewide < 0.05]
    x.out(f"{len(calls.fish)} fish: the genome-wide 5% threshold of |z| is {g95:.2f} (1%: {g99:.2f}), i.e. a share of 0.5 +- {g95 / (2 * np.sqrt(sc['n'])):.3f}; {len(reg)} regions at |z| >= {x.an['trd_zmin']:g}, {len(sig)} significant at genome-wide p < 0.05"
          + ''.join(f"\n  {r.LG} {r.start / 1e6:.1f}-{r.end / 1e6:.1f} Mb: share {r.mean_share:.2f} (peak {r.share_at_peak:.2f}, z {r.z_peak:.1f}, p {r.p_genomewide:.4f})" + ('  [chromosome tip: low confidence]' if r.in_tip else '') for r in sig.itertuples()))
    return calls


def step_linkscan(x):
    a = x.a
    if not a.lg or not a.region:
        raise ValueError('linkscan needs --lg and --region START-END (Mb)')
    calls = x.calls()
    chrom = calls.chrom_of(a.lg)
    lo, hi = mb_pair(a.region)
    best, own = T.linkscan(calls, chrom, lo, hi, int(x.an['linkscan_smooth']))
    x.write('linkscan_best_per_chromosome.tsv', best.round(3)); x.write('linkscan_own_chromosome.tsv', own.round(3))
    top = best[best.chrom != chrom].iloc[0]
    x.out(f"the state of {calls.lg_of(chrom)} {a.region} Mb is inherited together with: its own chromosome up to {own.link.max():.2f}; the strongest link elsewhere is {top.LG} at {top.best_Mb:.1f} Mb ({top.link:.2f}, agreement {top.agree:.2f}); 0 = independent, 1 = always together")
    return calls


def step_publish(x):
    env = x.env()
    if not os.path.abspath(x.root).startswith(os.path.abspath(env.local_master)):
        raise ValueError(f'{x.root} is not inside the local master {env.local_master}, so it has no Dropbox twin')
    if not os.path.isdir(x.root):
        raise ValueError(f'{x.root} does not exist')
    x.out(f'uploading {x.root} (the whole genome folder) to its Dropbox twin and verifying it ...')
    from variant_calling import build as B
    x.out(f'uploaded and verified: {(x.uploader or B.default_uploader(env))(x.root)}')


STEP_FN = {'publish': step_publish, 'marey': step_marey, 'rate': step_rate, 'sexscan': step_sexscan, 'finemap': step_finemap, 'crossovers': step_crossovers, 'trd': step_trd, 'linkscan': step_linkscan}


# ------------------------------------------------------------------ provenance and upload
def provenance(x, argv):
    import scipy
    def git():
        try:
            return subprocess.run(['git', '-C', os.path.dirname(os.path.abspath(__file__)), 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True, timeout=10).stdout.strip() or None
        except Exception:                                                  # noqa
            return None
    inputs = {}
    for k, p in x.inputs.items():
        if os.path.isfile(p):
            inputs[k] = dict(path=p, bytes=os.path.getsize(p), sha1=D.sha1(p), modified=datetime.datetime.fromtimestamp(os.path.getmtime(p)).isoformat(timespec='seconds'))
        else:
            inputs[k] = dict(path=p)
    return dict(step=x.a.step, genome=x.a.genome, run=x.run, argv=argv, written=sorted(os.path.basename(p) for p in x.written), when=datetime.datetime.now().isoformat(timespec='seconds'), host=socket.gethostname(),
                python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__, code_commit=git(), analysis_config=x.an, inputs=inputs)


def upload(x):
    """The output folder to its Dropbox twin, verified (cichlid_common.runner.sync through the lab's local/cloud master folders, as module 05 does)."""
    env = x.env()
    if not os.path.abspath(x.dir).startswith(os.path.abspath(env.local_master)):
        raise ValueError(f'{x.dir} is not inside the local master {env.local_master}, so it has no Dropbox twin: leave out --analysis-dir / --out, or write inside the local master')
    from variant_calling import build as B
    return (x.uploader or B.default_uploader(env))(x.dir)


PATH_ARGS = ('out', 'calls', 'analysis_dir', 'counts', 'counts_before', 'fish_table', 'local_master', 'fasta')


def clean_path(p):
    """A path without doubled slashes, '.' or '..' (a HOME that ends in a slash makes ~/x into //x, and the FileManager maps local to Dropbox paths by string replacement, which a doubled slash defeats).
    A trailing slash is kept."""
    if not p:
        return p
    q = os.path.normpath(p)
    return q + '/' if p.endswith('/') and not q.endswith('/') else q


def main(argv=None, out=print, uploader=None, downloader=None):
    av = sys.argv[1:] if argv is None else list(argv)
    a = build_parser().parse_args(av)
    for name in PATH_ARGS:
        setattr(a, name, clean_path(getattr(a, name)))
    if a.mask and a.mask != 'none':
        a.mask = clean_path(a.mask)
    from bc_analysis.genotyper import sites as S
    cfg = S.load_config(a.config)
    try:
        x = Ctx(a, cfg, out, uploader, downloader)
        STEP_FN[a.step](x)
        if a.step == 'publish':
            return 0
        with open(f'{x.dir}provenance.json', 'w') as fh:
            json.dump(provenance(x, ['analyze'] + av), fh, indent=1, default=str)
        out(f'wrote {len(x.written)} files + provenance.json to {x.dir}')
        if a.upload:
            out(f'uploaded and verified: {upload(x)}')
    except (ValueError, OSError, RuntimeError, KeyError) as e:
        print(f'error: {e}', file=sys.stderr)
        return 2
    return 0
