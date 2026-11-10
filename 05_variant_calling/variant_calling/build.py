"""`cohort build`: a cohort VCF with bcftools (module 05; SPEC decisions 59 to 63).

Stages: check the cohort (BLOCK stops everything) -> plan and preflight -> stage the BAMs from Dropbox -> `bcftools mpileup | bcftools call -m` per contig, several at a time
-> concat, rename the columns to the SampleIDs, norm, fill-tags, the biallelic-SNP subset -> a manifest -> upload and verify (`Outputs/PedigreeVCFs/`) -> delete the staged BAMs.
Calling follows Blumer et al. 2025 (their parameters are the defaults of config/cohort_build.yaml); lab and farm stock get no shared Hardy-Weinberg prior (groups per_sample).

Without --run only the plan is printed. A TRIAL (--max-samples) calls the smallest BAMs only, writes under the work folder and is never uploaded, so it cannot reach the folder that
`sexloci run-folder` analyses. Resumable: staged BAMs and finished contigs are kept, and a rerun does only what is missing.

bcftools names the columns of a multi-sample pileup after the read groups, and our BAMs can hold one read group per run; `--ignore-RG` makes it one column per BAM, in the order of
the BAM list, and the columns are renamed to the SampleIDs at the end, after the count is checked. (Built without bcftools at hand: the first real run is the check of these options.)"""
import datetime
import getpass
import json
import os
import shlex
import shutil
import socket
import subprocess
import sys
import time
import types
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import yaml

from cichlid_common.common_light import config_path

DEFAULTS = dict(
    bcftools=dict(min_base_quality=20, min_mapping_quality=13, indel_min_reads=2, indel_min_fraction=0.05, prior=0.001, mpileup_annotate=['FORMAT/AD', 'FORMAT/DP'], call_format_fields=['GQ']),
    fill_tags=['AN', 'AC', 'NS', 'F_MISSING', 'ExcHet'], contigs='chromosomes', min_contig_bp=5_000_000, workers=8, stage_parallel=4, min_free_factor=1.2,
    output_subdir='Outputs/PedigreeVCFs/', work_subdir='Cohorts/')
TOOLS = ('bcftools', 'samtools', 'rclone')


class BuildError(RuntimeError):
    pass


def load_build_config(path=None):
    cfg = {k: (dict(v) if isinstance(v, dict) else v) for k, v in DEFAULTS.items()}
    p = path or config_path('cohort_build.yaml')
    if os.path.exists(p):
        raw = yaml.safe_load(open(p)) or {}
        for k, v in raw.items():
            if isinstance(v, dict) and isinstance(cfg.get(k), dict):
                cfg[k].update(v)
            else:
                cfg[k] = v
    return cfg


# ------------------------------------------------------------------ where things are
class Env:
    """Local and Dropbox locations. From the lab file manager (online), or given explicitly (tests, offline)."""
    def __init__(self, local_master, cloud_master, genome, genome_fa):
        self.local_master = local_master.rstrip('/') + '/'
        self.cloud_master = cloud_master.rstrip('/') + '/'
        self.genome, self.genome_fa = genome, genome_fa

    @classmethod
    def from_file_manager(cls, genome):
        from cichlid_common.file_manager import FileManager
        fm = FileManager()
        fm.setGenome(genome)
        return cls(fm.localMasterDir, fm.cloudMasterDir, genome, fm.localGenomeFile)

    def cloud_sample_dir(self, sample):
        return f'{self.cloud_master}Bamfiles/{self.genome}/{sample}/'

    def project_like(self):
        """What cichlid_common.runner.sync wants: an object with .fm.localMasterDir and .fm.cloudMasterDir."""
        return types.SimpleNamespace(fm=types.SimpleNamespace(localMasterDir=self.local_master, cloudMasterDir=self.cloud_master))


# ------------------------------------------------------------------ pure pieces: regions, contigs, commands
def parse_region(text):
    """'NC_1' -> ('NC_1', None, None); 'NC_1:1-100' or 'NC_1:1,000-100' -> ('NC_1', 1, 100)."""
    if ':' not in text:
        return text, None, None
    contig, rng = text.rsplit(':', 1)
    try:
        a, b = rng.replace(',', '').split('-')
        return contig, int(a), int(b)
    except ValueError:
        raise BuildError(f'a region looks like CONTIG or CONTIG:START-END, not {text!r}')


def region_string(contig, start=None, end=None):
    return contig if start is None else f'{contig}:{start}-{end}'


def read_fai(path):
    return [(l.split('\t')[0], int(l.split('\t')[1])) for l in open(path) if l.strip()]


def choose_contigs(fai, mode, min_bp, region=None):
    """[(contig, start, end)] to call. A region wins; 'chromosomes' = the NC_ contigs of at least min_bp; 'all' = every contig; a list = those named."""
    lengths = dict(fai)
    if region:
        contig, a, b = parse_region(region)
        if contig not in lengths:
            raise BuildError(f'{contig} is not a contig of the reference')
        return [(contig, a, b)]
    if mode == 'chromosomes':
        out = [(c, None, None) for c, n in fai if c.startswith('NC_') and n >= min_bp]
    elif mode == 'all':
        out = [(c, None, None) for c, _ in fai]
    else:
        names = mode if isinstance(mode, (list, tuple)) else str(mode).split(',')
        missing = [c for c in names if c not in lengths]
        if missing:
            raise BuildError(f'not contigs of the reference: {", ".join(missing)}')
        out = [(c, None, None) for c in names]
    if not out:
        raise BuildError(f'no contig is selected (contigs: {mode}, min_contig_bp: {min_bp})')
    return out


def mpileup_cmd(ref, bamlist, region, cfg):
    b = cfg['bcftools']
    cmd = ['bcftools', 'mpileup', '-f', ref, '-b', bamlist, '-r', region, '--ignore-RG', '-q', str(b['min_mapping_quality']), '-Q', str(b['min_base_quality'])]
    if b.get('indel_min_reads'):
        cmd += ['-p', '-m', str(b['indel_min_reads']), '-F', str(b['indel_min_fraction'])]
    if b.get('mpileup_annotate'):
        cmd += ['-a', ','.join(b['mpileup_annotate'])]
    return cmd + ['-Ou']


def call_cmd(out, cfg, groups='per_sample'):
    if groups != 'per_sample':
        raise BuildError(f"groups {groups!r} is not built yet: lab and farm stock use per_sample (no shared Hardy-Weinberg prior, SPEC decision 63)")
    b = cfg['bcftools']
    cmd = ['bcftools', 'call', '-m', '-v', '-P', str(b['prior']), '-G', '-']
    if b.get('call_format_fields'):
        cmd += ['-f', ','.join(b['call_format_fields'])]
    return cmd + ['-Oz', '-o', out]


def contig_script(ref, bamlist, region, out, cfg, groups='per_sample'):
    """One bash script for one contig: the pileup piped into the caller, then the index. pipefail: either tool failing fails the script."""
    return ('set -o pipefail; ' + shlex.join(mpileup_cmd(ref, bamlist, region, cfg)) + ' | ' + shlex.join(call_cmd(out, cfg, groups)) + ' && ' + shlex.join(['bcftools', 'index', '-t', out]))


# ------------------------------------------------------------------ running things
def run_cmd(argv, input_text=None):
    r = subprocess.run(argv, capture_output=True, text=True, input=input_text)
    return r.returncode, r.stdout, r.stderr


def run_script(script):
    r = subprocess.run(['bash', '-c', script], capture_output=True, text=True)
    return r.returncode, r.stdout, r.stderr


def must(rc_out_err, what):
    rc, out, err = rc_out_err
    if rc:
        raise BuildError(f'{what} failed (exit {rc}): {err.strip()[-1500:]}')
    return out


# ------------------------------------------------------------------ the plan
class Plan:
    pass


def _bytes(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def make_plan(name, spec, members, ledger, env, cfg, region=None, contigs=None, max_samples=None):
    p = Plan()
    p.name, p.kind, p.genome, p.caller, p.groups = name, spec.get('kind'), spec.get('genome'), spec.get('caller'), spec.get('groups')
    led = ledger.dropna(how='all').drop_duplicates('SampleID', keep='last').set_index('SampleID') if len(ledger) else pd.DataFrame()
    size = {s: _bytes(led.at[s, 'BamSize']) if 'BamSize' in led.columns and s in led.index else None for s in members}
    p.all_members, p.trial = list(members), bool(max_samples)
    if max_samples:
        keep = set(sorted(members, key=lambda s: (size[s] is None, size[s] or 0))[:max_samples])
        members = [s for s in members if s in keep]
    p.members, p.sizes = list(members), {s: size[s] for s in members}
    p.unknown_sizes = [s for s in members if size[s] is None]
    p.total_bytes = sum(v for v in p.sizes.values() if v)
    p.coverage = {s: _bytes(led.at[s, 'Coverage']) if 'Coverage' in led.columns and s in led.index else None for s in members}
    base = name + ('.trial' if p.trial else '')
    p.work = f"{env.local_master}{cfg['work_subdir']}{name}{'_trial' if p.trial else ''}/"
    p.staged, p.per_contig = p.work + 'staged/', p.work + 'per_contig/'
    out_dir = (p.work + 'trial_out/') if p.trial else f"{env.local_master}{cfg['output_subdir']}"
    p.out_dir = out_dir
    p.files = dict(vcf=f'{out_dir}{base}.vcf.gz', snps=f'{out_dir}{base}.snps.vcf.gz', manifest=f'{out_dir}{base}.manifest.json')
    p.region, p.cfg, p.contig_mode = region, cfg, contigs or cfg['contigs']
    p.contigs = None
    p.notes = []
    fai = env.genome_fa + '.fai'
    if os.path.exists(fai):
        p.contigs = choose_contigs(read_fai(fai), p.contig_mode, cfg['min_contig_bp'], region)
    return p


def contig_label(c):
    contig, a, b = c
    return contig if a is None else f'{contig}_{a}-{b}'


def preflight(env, plan, cfg, free_bytes=None, which=shutil.which, run=run_cmd):
    """Problems that would stop a run, one sentence each; empty when it can start."""
    bad = []
    for t in TOOLS:
        if not which(t):
            bad.append(f'{t} is not on the PATH (activate the environment, then `source tools/env.sh`)')
    if not bad:
        _, out, err = run(['bcftools', 'call'])
        if '--group-samples' not in out + err:
            bad.append("this bcftools has no `call --group-samples` (needed for -G -, no shared Hardy-Weinberg prior): see `bcftools call` usage")
        _, out, err = run(['bcftools', 'mpileup'])
        if '--ignore-RG' not in out + err:
            bad.append("this bcftools has no `mpileup --ignore-RG` (needed for one column per BAM)")
    if not os.path.exists(env.genome_fa):
        bad.append(f'the reference FASTA is not at {env.genome_fa}')
    elif plan.contigs is None:
        bad.append(f'{env.genome_fa}.fai is missing (`samtools faidx {env.genome_fa}`)')
    if plan.unknown_sizes:
        bad.append(f'{len(plan.unknown_sizes)} members have no BamSize in the ledger, so the disk cannot be planned and the download cannot be verified: ' + ', '.join(plan.unknown_sizes[:5]))
    if free_bytes is None:
        probe = plan.work
        while probe and not os.path.exists(probe):
            probe = os.path.dirname(probe.rstrip('/'))
        free_bytes = shutil.disk_usage(probe or '/').free
    plan.free_bytes = free_bytes
    need = plan.total_bytes * cfg['min_free_factor']
    if need > free_bytes:
        bad.append(f'the disk has {free_bytes / 1e9:.0f} GB free but the BAMs need {plan.total_bytes / 1e9:.0f} GB x {cfg["min_free_factor"]:g}')
    return bad


def format_plan(plan, problems, cfg, env):
    out = [f'cohort {plan.name}: {len(plan.members)} samples ({plan.kind}), genome {plan.genome}, caller {plan.caller}, groups {plan.groups}'
           + (f'  TRIAL: the {len(plan.members)} smallest BAMs of {len(plan.all_members)}, local only, never uploaded' if plan.trial else '')]
    out.append(f'BAMs to stage: {len(plan.members)} files, {plan.total_bytes / 1e9:.1f} GB (ledger BamSize); free disk {getattr(plan, "free_bytes", 0) / 1e9:.0f} GB; staged under {plan.staged}')
    if plan.contigs:
        out.append(f'contigs: {len(plan.contigs)} ({plan.contig_mode}' + (f', region {plan.region}' if plan.region else '') + f'), {cfg["workers"]} at a time: ' + ', '.join(contig_label(c) for c in plan.contigs[:4]) + (' ...' if len(plan.contigs) > 4 else ''))
        c = plan.contigs[0]
        bamlist = plan.work + 'bams.list'
        out.append('one contig:  ' + contig_script(env.genome_fa, bamlist, region_string(*c), f'{plan.per_contig}{contig_label(c)}.vcf.gz', cfg, plan.groups or 'per_sample'))
    else:
        out.append('contigs: not listed (the reference .fai is missing)')
    out.append('then:        bcftools concat | reheader (columns renamed to the SampleIDs) | norm | +fill-tags -t ' + ','.join(cfg['fill_tags']) + ' | view -m2 -M2 -v snps')
    out.append('outputs:     ' + ', '.join(os.path.basename(v) for v in plan.files.values()) + (f' in {plan.out_dir} (local only)' if plan.trial else f' in {plan.out_dir}, uploaded and verified on Dropbox'))
    out.append('time:        not measured (the first real build gives the figure)')
    for p in problems:
        out.append(f'[PREFLIGHT] {p}')
    out.append('PLAN ONLY: nothing was staged, called or written. Add --run to build.')
    return '\n'.join(out)


# ------------------------------------------------------------------ staging
def stage_bams(env, plan, cfg, log=print, run=run_cmd):
    os.makedirs(plan.staged, exist_ok=True)

    def one(s):
        d = f'{plan.staged}{s}/'
        bam, bai = f'{d}{s}.all.bam', f'{d}{s}.all.bai'
        want = plan.sizes.get(s)
        ok = os.path.exists(bam) and os.path.exists(bai) and (want is None or os.path.getsize(bam) == int(want)) and run(['samtools', 'quickcheck', '-q', bam])[0] == 0
        if ok:
            return s, 'kept'
        os.makedirs(d, exist_ok=True)
        must(run(['rclone', 'copy', env.cloud_sample_dir(s), d, '--include', f'{s}.all.bam', '--include', f'{s}.all.bai']), f'downloading the BAM of {s}')
        if not (os.path.exists(bam) and os.path.exists(bai)):
            raise BuildError(f'{s}: {s}.all.bam or {s}.all.bai is not in {env.cloud_sample_dir(s)}')
        if want is not None and os.path.getsize(bam) != int(want):
            raise BuildError(f'{s}: the downloaded BAM is {os.path.getsize(bam)} bytes, the ledger says {int(want)}')
        rc, _, err = run(['samtools', 'quickcheck', '-q', bam])
        if rc:
            raise BuildError(f'{s}: samtools quickcheck says the downloaded BAM is damaged: {err.strip()[-300:]}')
        link = bam + '.bai'                                       # htslib finds either name; both are there so that nothing depends on that
        if not os.path.lexists(link):
            os.symlink(os.path.basename(bai), link)
        return s, 'downloaded'

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=cfg['stage_parallel']) as ex:
        results = list(ex.map(one, plan.members))
    kept = sum(1 for _, r in results if r == 'kept')
    log(f'[stage] {len(results) - kept} BAMs downloaded, {kept} already staged ({(time.time() - t0) / 60:.1f} min)')
    with open(plan.work + 'bams.list', 'w') as f:
        f.write('\n'.join(f'{plan.staged}{s}/{s}.all.bam' for s in plan.members) + '\n')
    with open(plan.work + 'samples.txt', 'w') as f:
        f.write('\n'.join(plan.members) + '\n')
    return results


# ------------------------------------------------------------------ calling
def call_contigs(env, plan, cfg, workers=None, log=print):
    os.makedirs(plan.per_contig, exist_ok=True)
    bamlist = plan.work + 'bams.list'
    n = len(plan.members)

    def one(c):
        label = contig_label(c)
        out = f'{plan.per_contig}{label}.vcf.gz'
        done = out + '.done'
        if os.path.exists(done) and os.path.exists(out + '.tbi'):
            return dict(contig=label, file=out, records=int(open(done).read().split()[0]), seconds=0.0, skipped=True)
        t0 = time.time()
        for f in (out, out + '.tbi', done):
            if os.path.exists(f):
                os.remove(f)
        rc, _, err = run_script(contig_script(env.genome_fa, bamlist, region_string(*c), out, cfg, plan.groups or 'per_sample'))
        if rc:
            raise BuildError(f'calling {label} failed (exit {rc}): {err.strip()[-1500:]}')
        cols = must(run_cmd(['bcftools', 'query', '-l', out]), f'reading the samples of {label}').split()
        if len(cols) != n:
            raise BuildError(f'{label}: the call has {len(cols)} sample columns for {n} BAMs: with --ignore-RG there must be one per BAM, so the columns cannot be renamed safely')
        count = must(run_cmd(['bcftools', 'index', '-n', out]), f'counting the records of {label}').strip().split()[0]
        open(done, 'w').write(f'{count} {time.time() - t0:.1f}\n')
        log(f'[call] {label}: {count} variants, {(time.time() - t0) / 60:.1f} min')
        return dict(contig=label, file=out, records=int(count), seconds=time.time() - t0, skipped=False)

    errors, results = [], []
    with ThreadPoolExecutor(max_workers=workers or cfg['workers']) as ex:
        futs = [(c, ex.submit(one, c)) for c in plan.contigs]
        for c, f in futs:
            try:
                results.append(f.result())
            except BuildError as e:
                errors.append(str(e))
    if errors:
        raise BuildError(f'{len(errors)} of {len(plan.contigs)} contigs failed; the finished ones are kept and a rerun does only the rest. First: {errors[0]}')
    skipped = sum(1 for r in results if r['skipped'])
    if skipped:
        log(f'[call] {skipped} contigs were already done')
    return results


# ------------------------------------------------------------------ merging
def merge(env, plan, cfg, results, log=print):
    ref = env.genome_fa
    os.makedirs(plan.out_dir, exist_ok=True)
    w = plan.work
    listing = w + 'contigs.list'
    open(listing, 'w').write('\n'.join(r['file'] for r in results) + '\n')
    raw, renamed, normed = w + 'merged.raw.vcf.gz', w + 'merged.renamed.vcf.gz', w + 'merged.norm.vcf.gz'
    must(run_cmd(['bcftools', 'concat', '-f', listing, '-Oz', '-o', raw]), 'concat')
    must(run_cmd(['bcftools', 'reheader', '-s', w + 'samples.txt', '-o', renamed, raw]), 'renaming the sample columns')
    names = must(run_cmd(['bcftools', 'query', '-l', renamed]), 'reading the renamed samples').split()
    if names != plan.members:
        raise BuildError('the renamed sample columns are not the cohort members in order: ' + ', '.join(names[:4]))
    must(run_cmd(['bcftools', 'norm', '-f', ref, '-Oz', '-o', normed, renamed]), 'norm')
    final = plan.files['vcf']
    must(run_cmd(['bcftools', '+fill-tags', normed, '-Oz', '-o', final, '--', '-t', ','.join(cfg['fill_tags'])]), 'fill-tags')
    must(run_cmd(['bcftools', 'index', '-t', '-f', final]), 'indexing the VCF')
    snps = plan.files['snps']
    must(run_cmd(['bcftools', 'view', '-m2', '-M2', '-v', 'snps', '-Oz', '-o', snps, final]), 'the SNP subset')
    must(run_cmd(['bcftools', 'index', '-t', '-f', snps]), 'indexing the SNP VCF')
    counts = {k: int(must(run_cmd(['bcftools', 'index', '-n', p]), f'counting {k}').strip().split()[0]) for k, p in (('vcf', final), ('snps', snps))}
    for f in (listing, raw, renamed, normed):
        if os.path.exists(f):
            os.remove(f)
    log(f"[merge] {counts['vcf']} variants, {counts['snps']} biallelic SNPs, {len(names)} samples")
    return counts


# ------------------------------------------------------------------ manifest, upload
def tool_line(tool):
    rc, out, err = run_cmd([tool, '--version'])
    return ((out or err).strip().splitlines() or ['n/a'])[0] if not rc else 'n/a'


def write_manifest(plan, spec, report, results, counts, cfg, env, region):
    try:
        from cichlid_common.runner import git_state
        commit, dirty = git_state(None)
    except Exception:
        commit, dirty = 'nogit', None
    man = dict(cohort=plan.name, kind=plan.kind, genome=plan.genome, caller=plan.caller, groups=plan.groups, trial=plan.trial, region=region,
               created=datetime.datetime.now().isoformat(timespec='seconds'), user=getpass.getuser(), host=socket.gethostname(), git_commit=commit, git_dirty=dirty,
               command=' '.join(sys.argv), tools={t: tool_line(t) for t in ('bcftools', 'samtools')},
               parameters=dict(bcftools=cfg['bcftools'], fill_tags=cfg['fill_tags'], contigs=plan.contig_mode, min_contig_bp=cfg['min_contig_bp'], note='--ignore-RG; columns renamed to the SampleIDs'),
               members=[dict(sample=s, bam=f'{s}.all.bam', bam_bytes=plan.sizes.get(s), coverage=plan.coverage.get(s)) for s in plan.members],
               n_members_in_cohort=len(plan.all_members), check=report.to_dict(), contigs=[dict(contig=r['contig'], records=r['records'], seconds=round(r['seconds'], 1), resumed=r['skipped']) for r in results],
               outputs={k: dict(file=os.path.basename(plan.files[k]), bytes=os.path.getsize(plan.files[k]) if os.path.exists(plan.files[k]) else None, records=counts.get(k)) for k in ('vcf', 'snps')},
               spec={k: v for k, v in spec.items() if k != 'name'})
    json.dump(man, open(plan.files['manifest'], 'w'), indent=1, default=str)
    return man


def default_uploader(env):
    from cichlid_common.runner import sync
    proj = env.project_like()
    return lambda path: sync(proj, path)


def default_cloud_exists(env):
    def exists(local_path):
        cloud = local_path.replace(env.local_master, env.cloud_master)
        rc, out, _ = run_cmd(['rclone', 'lsf', cloud])
        return rc == 0 and bool(out.strip())
    return exists


# ------------------------------------------------------------------ the whole build
def build_cohort(name, spec, sample_dt, ledger, env, *, cfg=None, rules=None, region=None, contigs=None, max_samples=None, run_it=False, workers=None, keep_staged=False,
                 keep_work=False, upload=True, force=False, log=print, free_bytes=None, uploader=None, cloud_exists=None, which=shutil.which):
    from variant_calling import rules as R
    cfg = cfg or load_build_config()
    report = R.check_cohort(name, spec, sample_dt, ledger, rules or R.load_rules())
    if report.blocked:
        return dict(status='blocked', report=report)
    plan = make_plan(name, spec, report.members, ledger, env, cfg, region, contigs, max_samples)
    problems = preflight(env, plan, cfg, free_bytes, which)
    if not run_it:
        return dict(status='plan', plan=plan, report=report, problems=problems, text=format_plan(plan, problems, cfg, env))
    if problems:
        raise BuildError('the preflight failed, nothing was started: ' + '; '.join(problems))
    upload = upload and not plan.trial
    existing = [p for p in plan.files.values() if os.path.exists(p)]
    if existing and not force:
        raise BuildError(f'{os.path.basename(existing[0])} already exists in {plan.out_dir}: use --force to rebuild it')
    if upload and not force:
        exists = cloud_exists or default_cloud_exists(env)
        if exists(plan.files['vcf']):
            raise BuildError(f"{os.path.basename(plan.files['vcf'])} is already on Dropbox in {cfg['output_subdir']}: use --force to replace it (module 06 may have analysed it)")
    os.makedirs(plan.work, exist_ok=True)
    t0 = time.time()
    stage_bams(env, plan, cfg, log)
    results = call_contigs(env, plan, cfg, workers, log)
    counts = merge(env, plan, cfg, results, log)
    man = write_manifest(plan, spec, report, results, counts, cfg, env, region)
    if upload:
        up = uploader or default_uploader(env)
        for k in ('vcf', 'snps'):
            up(plan.files[k])
            up(plan.files[k] + '.tbi')
        up(plan.files['manifest'])
        log(f'[upload] {len(plan.files) * 2 - 1} files uploaded and verified -> {cfg["output_subdir"]}')
    if not keep_staged:
        shutil.rmtree(plan.staged, ignore_errors=True)
    if not keep_work:
        shutil.rmtree(plan.per_contig, ignore_errors=True)
    log(f'[done] {plan.name}: {len(plan.members)} samples in {(time.time() - t0) / 60:.1f} min' + (' (TRIAL: local only)' if plan.trial else ''))
    return dict(status='built', plan=plan, report=report, manifest=man, counts=counts)
