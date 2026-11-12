"""python -m bc_analysis genotype sites    --vcf <bc_family VCF> [--out DIR]                  informative sites per BC dam, from the parents' cohort VCF
   python -m bc_analysis genotype counts   [--run] [--jobs N] [--samples A,B]                  allele counts of every BC fish at its dam's sites (a plan unless --run); downloads each BAM, counts, deletes it
   python -m bc_analysis genotype hmm      [--jobs N]                                          the HMM on the counts: calls and posteriors on the 100 kb grid, breakpoints, per-fish QC, end support
   python -m bc_analysis genotype compare --b hmm_int2 [--a hmm]                               what a second run (e.g. with --interior-min-mb 2 --tag _int2) changed
   python -m bc_analysis genotype simulate [--scenario S1|S2] [--coverage 0.3 0.5] [--fish N]  the simulated comparison with the original script's model

Outputs go to <local master>/Outputs/BCGenotypes/<genome>/{sites,counts,hmm}/ (local only) unless --out is given. The cross: YH_008_m x MC-1R6R-f -> the F1 sires MCYHF1-002 and -003 -> backcrossed to MC mothers.
Offline (tests, no lab sheets): --sample-table and --alignment-table (exports), --local-master, --cloud-master and --fasta go together."""
import argparse
import concurrent.futures as cf
import json
import os
import sys
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd

from . import counts as C
from . import hmm as H
from . import sites as S


def lg_label(chrom):
    """LG number of a Mzebra_GT3_NCBI chromosome (NC_135167.1 = LG1); other names unchanged."""
    if chrom.startswith('NC_135') and chrom[3:9].isdigit():
        n = int(chrom[3:9]) - 135166
        return f'LG{n}' if 1 <= n <= 23 else chrom
    return chrom


def build_parser():
    ap = argparse.ArgumentParser(prog='python -m bc_analysis genotype', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('step', choices=['sites', 'counts', 'hmm', 'compare', 'simulate'])
    ap.add_argument('--genome', default='Mzebra_GT3_NCBI')
    ap.add_argument('--config', default=None)
    ap.add_argument('--out', default=None, help='the genome folder of the outputs (default <local master>/Outputs/BCGenotypes/<genome>/)')
    ap.add_argument('--vcf', default=None, help='sites: the parents cohort VCF (default <local master>/Outputs/PedigreeVCFs/bc_family.snps.vcf.gz)')
    ap.add_argument('--fai', default=None, help='the reference index (default: <fasta>.fai)')
    ap.add_argument('--run', action='store_true', help='counts: really download and count (default: print the plan)')
    ap.add_argument('--jobs', type=int, default=None)
    ap.add_argument('--samples', default=None, help='counts: only these fish (comma-separated)')
    ap.add_argument('--max-samples', type=int, default=None)
    ap.add_argument('--keep-bams', action='store_true')
    ap.add_argument('--redo', action='store_true', help='counts: count again fish that are current')
    ap.add_argument('--rho', type=float, default=None, help='hmm: crossovers per bp (sensitivity runs)')
    ap.add_argument('--interior-min-mb', type=float, default=None, help='hmm: min_interior_segment_mb')
    ap.add_argument('--terminal-min-mb', type=float, default=None, help='hmm: min_terminal_segment_mb')
    ap.add_argument('--tag', default='', help='hmm: suffix of the output folder (hmm<tag>/), for sensitivity runs')
    ap.add_argument('--a', default='hmm', help='compare: the first run (a folder of the genome folder)')
    ap.add_argument('--b', default=None, help='compare: the second run, e.g. hmm_int2')
    ap.add_argument('--scenario', default='S2', choices=['S1', 'S2'])
    ap.add_argument('--coverage', type=float, nargs='+', default=[0.3, 0.5, 1.0])
    ap.add_argument('--fish', type=int, default=20)
    ap.add_argument('--seed', type=int, default=1)
    o = ap.add_argument_group('offline / tests')
    o.add_argument('--sample-table', default=None)
    o.add_argument('--alignment-table', default=None)
    o.add_argument('--local-master', default=None)
    o.add_argument('--cloud-master', default=None)
    o.add_argument('--fasta', default=None)
    return ap


def get_env(a):
    from variant_calling import build as B
    if a.local_master or a.cloud_master or a.fasta:
        if not (a.local_master and a.cloud_master and a.fasta):
            raise SystemExit('--local-master, --cloud-master and --fasta go together')
        return B.Env(a.local_master, a.cloud_master, a.genome, a.fasta)
    return B.Env.from_file_manager(a.genome)


def genome_dir(a, env=None):
    if a.out:
        return a.out.rstrip('/') + '/'
    env = env or get_env(a)
    return f'{env.local_master}Outputs/BCGenotypes/{a.genome}/'


def read_fai(a, env=None):
    path = a.fai or ((env.genome_fa + '.fai') if env is not None else None)
    if not path or not os.path.exists(path):
        return None
    from variant_calling import build as B
    return B.read_fai(path)


# ------------------------------------------------------------------ which fish, with which dam
def _dam_fish(sample_dt, cfg):
    """Fish whose dam is one of the configured BC dams (so the four 'BC1 hybrid' fish of another cross are not included), with their Sire and Dam."""
    df = sample_dt.copy()
    if 'BroodID' not in df.columns or 'Parents' not in df.columns:
        raise SystemExit('the sample table needs the columns SampleID, BroodID and Parents')
    df = df[df.Parents.astype(str).str.contains(',') & df.BroodID.astype(str).ne('')]
    df = df.assign(Sire=df.Parents.str.split(',').str[0].str.strip(), Dam=df.Parents.str.split(',').str[1].str.strip())
    return df[df.Dam.isin(cfg['sites']['dams'])][['SampleID', 'BroodID', 'Sire', 'Dam']].drop_duplicates('SampleID').reset_index(drop=True)


def _is_backcross(df, cfg):
    pre = tuple(cfg['sites'].get('brood_prefixes') or ())
    return df.BroodID.astype(str).str.startswith(pre) if pre else pd.Series(True, index=df.index)


def bc_fish(sample_dt, cfg):
    """DataFrame SampleID, BroodID, Sire, Dam of the BACKCROSS fish: the dam is one of the configured BC dams and the brood is a backcross brood (BroodID starts with a configured prefix, MCYH-BC1-)."""
    df = _dam_fish(sample_dt, cfg)
    return df[_is_backcross(df, cfg)].reset_index(drop=True)


def left_out(sample_dt, cfg):
    """The fish of a configured dam that are NOT backcross fish (an F1 brood, say), counted by brood, sire and dam."""
    df = _dam_fish(sample_dt, cfg)
    return df[~_is_backcross(df, cfg)].groupby(['BroodID', 'Sire', 'Dam']).size()


# ------------------------------------------------------------------ steps
def step_sites(a, cfg, out=print):
    env = None
    if not (a.vcf and a.out):
        env = get_env(a)
    vcf = a.vcf or f'{env.local_master}Outputs/PedigreeVCFs/bc_family.snps.vcf.gz'
    d = genome_dir(a, env) + 'sites/'
    t0 = time.time()
    res = S.define_sites(vcf, cfg, out)
    S.write_sites(res, d, vcf, cfg, read_fai(a, env))
    out(f'sites ({time.time() - t0:.0f} s), definition key {res["key"]}, written to {d}')
    for r in res['attrition']:
        pass
    att = pd.DataFrame(res['attrition']).pivot(index='step', columns='dam', values='sites').reindex(S.STEPS)
    out(att.to_string())
    return 0


def step_counts(a, cfg, out=print):
    from variant_calling import cohort_cli
    env = get_env(a)
    root = genome_dir(a, env)
    sites_dir = root + 'sites/'
    if not os.path.exists(sites_dir + 'definition.json'):
        raise SystemExit(f'no sites yet: run `genotype sites` first ({sites_dir}definition.json is missing)')
    key = json.load(open(sites_dir + 'definition.json'))['key']
    sample_dt, ledger = cohort_cli.load_tables(SimpleNamespace(sample_table=a.sample_table, alignment_table=a.alignment_table), a.genome)
    fish = bc_fish(sample_dt, cfg)
    if a.samples:
        want = set(a.samples.split(','))
        fish = fish[fish.SampleID.isin(want)]
    led = ledger.dropna(how='all').drop_duplicates('SampleID', keep='last').set_index('SampleID') if ledger is not None and len(ledger) else None
    in_ledger = set(led.index) if led is not None else set()
    missing = fish[~fish.SampleID.isin(in_ledger)].SampleID.tolist()
    fish = fish[fish.SampleID.isin(in_ledger)]

    def size(s):
        try:
            v = float(led.at[s, 'BamSize'])
            return v if v == v else None
        except (KeyError, ValueError, TypeError):
            return None
    fish = fish.assign(bytes=[size(s) for s in fish.SampleID])
    todo = fish if a.redo else fish[[not C.current(root, s, key, b) for s, b in zip(fish.SampleID, fish.bytes)]]
    if a.max_samples:
        todo = todo.head(a.max_samples)
    gb = todo.bytes.sum() / 1e9
    out(f'{len(fish)} BC fish with an alignment ({len(missing)} without one), definition key {key}; {len(todo)} to count ({gb:.0f} GB of BAMs), {len(fish) - len(todo)} current')
    out(fish.groupby(['BroodID', 'Sire', 'Dam']).size().to_string())
    lo = left_out(sample_dt, cfg)
    if len(lo):
        out('left out, not backcross broods (BroodID does not start with ' + ' or '.join(cfg['sites'].get('brood_prefixes') or []) + '): ' + '; '.join(f'{b} (sire {s}, dam {d}): {n}' for (b, s, d), n in lo.items()))
    if missing:
        out('not in the AlignmentDatabase: ' + ', '.join(missing[:12]) + (' ...' if len(missing) > 12 else ''))
    if not a.run:
        out('PLAN ONLY: nothing was downloaded or counted. Add --run to count.')
        return 0
    staged = f'{env.local_master}Cohorts/bc_counts/'
    os.makedirs(staged, exist_ok=True)
    jobs = a.jobs or 4
    done = [0]
    failed = []
    t0 = time.time()

    def one(r):
        m = C.count_sample(env, cfg, r.SampleID, r.Dam, f'{sites_dir}sites_{r.Dam}.tsv.gz', f'{sites_dir}targets_{r.Dam}.tsv', key, root, staged, r.bytes, a.keep_bams, C.RUNNERS['run'], C.RUNNERS['script'])
        done[0] += 1
        if m['status'] == 'ok':
            out(f"[counts] {done[0]}/{len(todo)} {r.SampleID}: {m['n_reads']} reads at {m['n_sites']} sites, {m['yh_fraction']:.1%} YH, {m['seconds']:.0f} s")
        else:
            failed.append((r.SampleID, m['error']))
            out(f"[counts] {done[0]}/{len(todo)} {r.SampleID}: FAILED {m['error']}")
        return m
    with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
        list(ex.map(one, list(todo.itertuples())))
    out(f'[counts] {len(todo) - len(failed)} counted, {len(failed)} failed, {(time.time() - t0) / 60:.1f} min')
    return 0 if not failed else 3


def _fish_job(args):
    sample, counts, lengths, hp = args
    chroms = {ch: (*counts[ch], lengths[ch]) for ch in lengths if ch in counts}
    posts, cs, eps, q = H.genotype_fish(chroms, hp)
    bps, ends, shorts, n_reads = [], [], [], 0.0
    for ch, post in posts.items():
        for sg in H.short_interior_segments(post, hp['bin_bp'], hp['short_segment_mb'], hp['call_threshold']):
            shorts.append(dict(sample=sample, chrom=ch, LG=lg_label(ch), **sg))
        for b in H.breakpoints(post, lengths[ch], hp):
            bps.append(dict(sample=sample, chrom=ch, LG=lg_label(ch), **b))
        ends.append(dict(sample=sample, chrom=ch, LG=lg_label(ch), **H.end_support(cs[ch], post, hp)))
        n_reads += float(cs[ch].sum())
    return sample, posts, eps, q, bps, ends, n_reads, shorts


def step_hmm(a, cfg, out=print):
    env = None if (a.out and a.fai) else get_env(a)
    root = genome_dir(a, env)
    hp = dict(cfg['hmm'])
    for k, v in (('rho', a.rho), ('min_interior_segment_mb', a.interior_min_mb), ('min_terminal_segment_mb', a.terminal_min_mb)):
        if v is not None:
            hp[k] = v
    fai = read_fai(a, env)
    if fai is None:
        raise SystemExit('the reference index is needed for the chromosome lengths: give --fai or a reachable reference')
    from variant_calling import build as B
    chosen = {c for c, _, _ in B.choose_contigs(fai, 'chromosomes', 5_000_000)}
    lengths = {c: L for c, L in fai if c in chosen}                       # whole chromosomes, in the order of the reference
    if not os.path.isdir(root + 'counts'):
        raise SystemExit(f'no counts in {root}counts: run `genotype counts --run` first')
    samples = sorted(f[:-len('.counts.json')] for f in os.listdir(root + 'counts') if f.endswith('.counts.json') and json.load(open(os.path.join(root, 'counts', f))).get('status') == 'ok')
    if not samples:
        raise SystemExit(f'no counts in {root}counts: run `genotype counts --run` first')
    outdir = root + f'hmm{a.tag}/'
    os.makedirs(outdir, exist_ok=True)
    out(f'{len(samples)} fish, {len(lengths)} chromosomes, bin {hp["bin_bp"] // 1000} kb, rho {hp["rho"]:g}, interior minimum {hp["min_interior_segment_mb"]} Mb, terminal minimum {hp["min_terminal_segment_mb"]} Mb')
    t0 = time.time()
    jobs = [(s, C.read_counts(root, s), lengths, hp) for s in samples]
    jn = a.jobs or 4
    results = []
    with cf.ProcessPoolExecutor(max_workers=jn) as ex:
        for k, r in enumerate(ex.map(_fish_job, jobs), 1):
            results.append(r)
            if k % 25 == 0 or k == len(jobs):
                out(f'[hmm] {k}/{len(jobs)} fish, {(time.time() - t0) / 60:.1f} min')
    bp = hp['bin_bp']
    grid = pd.concat([pd.DataFrame({'Chromosome': ch, 'Position': np.arange(int(np.ceil(L / bp))) * bp}) for ch, L in lengths.items()], ignore_index=True)
    post_cols = {s: np.concatenate([posts[ch] for ch in lengths]) for s, posts, *_ in results}
    shorts_all = [x for r in results for x in r[7]]
    pd.DataFrame(shorts_all, columns=['sample', 'chrom', 'LG', 'start', 'end', 'length', 'state', 'mean_posterior']).to_csv(outdir + 'short_segments.tsv', sep='\t', index=False)
    P = pd.DataFrame(post_cols)
    calls = (P > hp['call_threshold']).astype(int)
    pd.concat([grid, calls], axis=1).to_csv(outdir + 'genotypes_calls.csv')                           # the layout of the original OutputGenotypes.csv (an index column, Chromosome, Position, one column per fish)
    pd.concat([grid, P.round(4)], axis=1).to_csv(outdir + 'genotypes_posterior.csv.gz')
    pd.DataFrame([b for r in results for b in r[4]]).to_csv(outdir + 'breakpoints.tsv', sep='\t', index=False)
    pd.DataFrame([e for r in results for e in r[5]]).to_csv(outdir + 'end_support.tsv', sep='\t', index=False)
    lo, hi = hp['uncertain_band']
    qc = []
    for s, posts, eps, q, bps, ends, n_reads, shorts in results:
        nsw = len(bps)
        unc = float(np.mean(np.concatenate([(p > lo) & (p < hi) for p in posts.values()])))
        qc.append(dict(sample=s, reads_used=round(n_reads), eps=round(eps, 4), q_ref_in_het=round(q, 3), n_switches=nsw, switches_per_chromosome=round(nsw / len(lengths), 2), short_interior_segments=len(shorts),
                       terminal_low_confidence=sum(1 for b in bps if b['terminal_low_confidence']), uncertain_fraction=round(unc, 3),
                       flag=bool(eps > hp['flag_max_eps'] or nsw / len(lengths) > hp['flag_max_switches_per_chromosome'])))
    pd.DataFrame(qc).to_csv(outdir + 'sample_qc.tsv', sep='\t', index=False)
    json.dump(dict(parameters=hp, genome=a.genome, n_fish=len(samples), chromosomes=list(lengths), minutes=round((time.time() - t0) / 60, 1)), open(outdir + 'run_info.json', 'w'), indent=1)
    n_short, n_fc = len(shorts_all), len(samples) * len(lengths)
    out(f'{n_short} interior segments shorter than {hp["short_segment_mb"]} Mb in {len({x["sample"] for x in shorts_all})} of {len(samples)} fish ({100 * n_short / n_fc:.1f} per 100 fish-chromosomes): what an interior minimum would remove')
    out(f'wrote {outdir}: genotypes_calls.csv, genotypes_posterior.csv.gz, breakpoints.tsv, short_segments.tsv, end_support.tsv, sample_qc.tsv ({sum(q["flag"] for q in qc)} fish flagged)')
    return 0


def step_compare(a, cfg, out=print):
    """What did the second run change relative to the first? Bins called differently, breakpoints, short interior segments, terminal calls, map length."""
    if not a.b:
        raise SystemExit('compare needs --b: the folder of the second run, e.g. --b hmm_int2')
    env = None if a.out else get_env(a)
    root = genome_dir(a, env)
    da, db = root + a.a + '/', root + a.b + '/'
    for d in (da, db):
        if not os.path.exists(d + 'genotypes_calls.csv'):
            raise SystemExit(f'{d}genotypes_calls.csv is missing: run `genotype hmm` (with --tag) first')
    ca, cb = pd.read_csv(da + 'genotypes_calls.csv', index_col=0), pd.read_csv(db + 'genotypes_calls.csv', index_col=0)
    fish = [c for c in ca.columns[2:] if c in cb.columns]
    if not fish or len(ca) != len(cb):
        raise SystemExit('the two runs do not share fish and grid')
    diff = ca[fish].to_numpy() != cb[fish].to_numpy()

    def table(folder, name):
        p = folder + name
        return pd.read_csv(p, sep='\t') if os.path.exists(p) and os.path.getsize(p) > 1 else pd.DataFrame()
    bpa, bpb = table(da, 'breakpoints.tsv'), table(db, 'breakpoints.tsv')
    sha, shb = table(da, 'short_segments.tsv'), table(db, 'short_segments.tsv')

    def per_fish(df, col=None, cond=None):
        d = df if cond is None or df.empty else df[cond(df)]
        return d.groupby('sample').size().reindex(fish).fillna(0).astype(int) if not d.empty else pd.Series(0, index=fish)
    term = lambda df: df.distance_to_end < 2e6
    tab = pd.DataFrame(dict(sample=fish, bins_differ=diff.sum(0), breakpoints_a=per_fish(bpa).to_numpy(), breakpoints_b=per_fish(bpb).to_numpy(),
                            short_segments_a=per_fish(sha).to_numpy(), short_segments_b=per_fish(shb).to_numpy(),
                            terminal_calls_a=per_fish(bpa, cond=term).to_numpy(), terminal_calls_b=per_fish(bpb, cond=term).to_numpy()))
    tab.to_csv(f'{root}compare_{a.a}_vs_{a.b}.tsv', sep='\t', index=False)
    n_ch = ca.Chromosome.nunique()
    out(f'{a.a} against {a.b}: {len(fish)} fish, {diff.size:,} bins compared')
    out(f'  bins called differently: {int(diff.sum()):,} ({100 * diff.mean():.3f}%), in {int((diff.sum(0) > 0).sum())} fish')
    out(f'  breakpoints: {len(bpa)} -> {len(bpb)}; short interior segments: {len(sha)} -> {len(shb)}; calls within 2 Mb of a chromosome end: {int(tab.terminal_calls_a.sum())} -> {int(tab.terminal_calls_b.sum())}')
    out(f'  mean map length per chromosome: {100 * len(bpa) / (len(fish) * n_ch):.1f} cM -> {100 * len(bpb) / (len(fish) * n_ch):.1f} cM')
    return 0


def step_simulate(a, cfg, out=print):
    from . import simulate as SIM
    hp = dict(cfg['hmm'])
    t0 = time.time()
    rows = SIM.run_comparison(hp, tuple(a.coverage), a.scenario, a.fish, seed=a.seed)
    recall, false = SIM.summarize(rows)
    out(f"scenario {a.scenario}: {SIM.SCENARIOS[a.scenario]['label']}; {a.fish} fish x 6 chromosomes per coverage, {(time.time() - t0) / 60:.1f} min")
    out('\nrecall of the true crossovers by distance to the nearest chromosome end (number of crossovers):')
    out(recall.to_string())
    out('\nfalse calls and accuracy:')
    out(false.to_string())
    if a.out:
        os.makedirs(a.out, exist_ok=True)
        recall.to_csv(os.path.join(a.out, f'simulation_{a.scenario}_recall.tsv'), sep='\t')
        false.to_csv(os.path.join(a.out, f'simulation_{a.scenario}_false_calls.tsv'), sep='\t')
    return 0


def main(argv=None):
    a = build_parser().parse_args(argv)
    cfg = S.load_config(a.config)
    try:
        return {'sites': step_sites, 'counts': step_counts, 'hmm': step_hmm, 'compare': step_compare, 'simulate': step_simulate}[a.step](a, cfg)
    except (ValueError, OSError) as e:
        print(f'error: {e}', file=sys.stderr)
        return 2
