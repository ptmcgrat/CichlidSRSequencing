"""python -m align_pipeline run <Genome> [selection and options]      download -> align -> mark duplicates -> chimeric split -> statistics -> slices -> panel -> verified upload -> ledger row, many samples at once
python -m align_pipeline setup-database <Genome>               add the ledger's columns to the header of the genome's AlignmentDatabase worksheet (once, from ONE machine, before the first run)

This is the command line of v6 (alignFastQReads_v2.py) in the module layout, with the decisions of 2026-11-08: only the chimeric split, slices of LG2, LG10 and LG13 (slice_regions/<Genome>.tsv, `END` = the
chromosome's length), the panel stage (preflight, then a measurement per sample; --no-panel to run without), a lean ledger row, and time limits on transfers. The lab's FileManager is used as before for paths, selection
and transfers (cichlid_common.file_manager)."""
import argparse
import hashlib
import os
import re
import resource
import shutil
import subprocess
import sys
import time

from . import __version__
from . import alignment_db as adb
from . import bam_slices as sl
from . import panel as pn
from . import seq_pipeline as sp
from .disk_ledger import DiskLedger
from .procs import configure_rclone_env

HERE = os.path.dirname(os.path.abspath(__file__))
MODULE_DIR = os.path.dirname(HERE)
EXCLUDED_ECOGROUPS = ('Error', 'Exclude')                                      # samples flagged on purpose: not run unless --include-excluded


# ---------------------------------------------------------------- what the ledger row needs to know about a sample (pure, tested)
def _clean(v):
    return '' if v is None or (isinstance(v, float) and v != v) else v


def input_fingerprint(pairs):
    """A short hash of the runs a sample is aligned from: their RunIDs and uBAM sizes, sorted. The same runs and the same files give the same fingerprint; a run added or removed, or a uBAM re-made
    (for example re-downloaded from SRA instead of ENA) changes it. This is the input of the rerun rule (decision 37)."""
    items = []
    for run, size in pairs:
        size = _clean(size)
        items.append((str(run), '' if size == '' else str(int(float(size)))))
    return hashlib.sha1(';'.join(f'{r}:{z}' for r, z in sorted(items)).encode()).hexdigest()[:12]


def sample_meta(sample_dt, reads_dt, samples):
    """{sample: {Species, Instrument, ReadLength, RunIDs (joined with ',,' as in the existing rows), NRuns, InputFingerprint}} from the SampleDatabase and DNAReadsDatabase tables."""
    species = dict(zip(sample_dt.SampleID.astype(str), sample_dt.Species)) if 'Species' in sample_dt.columns else {}
    out = {}
    for s in samples:
        r = reads_dt[reads_dt.SampleID.astype(str) == s]
        runs = [str(x) for x in r.RunID]
        sizes = list(r.FileSize) if 'FileSize' in r.columns else [None] * len(runs)
        instrument = next((x for x in r.Instrument if _clean(x) != ''), '') if 'Instrument' in r.columns else ''
        length = next((x for x in r.ReadLength if _clean(x) != ''), '') if 'ReadLength' in r.columns else ''
        out[s] = dict(Species=_clean(species.get(s)), Instrument=instrument, ReadLength=int(float(length)) if length != '' else '', RunIDs=',,'.join(runs), NRuns=len(runs),
                      InputFingerprint=input_fingerprint(zip(runs, sizes)))
    return out


def short_version(text):
    """'2.28-r1209' -> '2.28', 'The Genome Analysis Toolkit (GATK) v4.6.2.0' -> '4.6.2.0': the same short form as `conda list`, which is what the existing rows hold."""
    m = re.search(r'\d+(?:\.\d+)+', text or '')
    return m.group(0) if m else (text or '').strip()


def tool_versions(run=subprocess.run, pysam_version=None):
    def first(cmd):
        out = run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        return short_version((out.stdout.strip() or out.stderr.strip()).split('\n')[0])
    versions = {'minimap2': first(['minimap2', '--version']), 'samtools': first(['samtools', '--version']), 'gatk': first(['gatk', '--version'])}
    if pysam_version is None:
        import pysam
        pysam_version = pysam.__version__
    versions['pysam'] = pysam_version
    return versions


def read_git_head(path):
    """The short commit of the checkout that contains `path`, read from the .git folder itself (no git program needed: an environment may have none). '' when there is no checkout or it cannot be read.
    Handles a branch checked out normally, a detached HEAD, refs kept only in packed-refs, and linked worktrees (where .git is a file)."""
    d = os.path.abspath(path)
    while True:
        g = os.path.join(d, '.git')
        if os.path.exists(g):
            break
        parent = os.path.dirname(d)
        if parent == d:
            return ''
        d = parent
    try:
        if os.path.isfile(g):                                                       # a linked worktree: ".git" is a file saying "gitdir: <path>"
            line = open(g).read().strip()
            if not line.startswith('gitdir:'):
                return ''
            g = line.split(':', 1)[1].strip()
            g = g if os.path.isabs(g) else os.path.normpath(os.path.join(d, g))
        head = open(os.path.join(g, 'HEAD')).read().strip()
        if not head.startswith('ref:'):                                             # detached HEAD: the commit itself
            return head[:7] if re.fullmatch(r'[0-9a-f]{40}', head) else ''
        ref = head.split(':', 1)[1].strip()
        common = g                                                                  # a worktree keeps its refs in the main .git folder
        if os.path.exists(os.path.join(g, 'commondir')):
            common = os.path.normpath(os.path.join(g, open(os.path.join(g, 'commondir')).read().strip()))
        ref_file = os.path.join(common, ref)
        if os.path.isfile(ref_file):
            sha = open(ref_file).read().strip()
            return sha[:7] if re.fullmatch(r'[0-9a-f]{40}', sha) else ''
        packed = os.path.join(common, 'packed-refs')
        if os.path.isfile(packed):
            for line in open(packed):
                parts = line.split()
                if len(parts) == 2 and parts[1] == ref and re.fullmatch(r'[0-9a-f]{40}', parts[0]):
                    return parts[0][:7]
    except (OSError, UnicodeDecodeError):
        pass
    return ''


def git_commit(path=MODULE_DIR, run=subprocess.run):
    """The short commit of the checkout, recorded in the ledger so BAMs made by an older pipeline can be found. `git` is asked first; when there is no git program or it refuses (an environment
    without git, a checkout owned by another user), the .git folder is read directly."""
    try:
        r = run(['git', '-C', path, 'rev-parse', '--short', 'HEAD'], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout.strip()
    except OSError:
        pass
    return read_git_head(path)


def pipeline_version(run=subprocess.run, path=MODULE_DIR):
    commit = git_commit(path=path, run=run)
    return f'align_pipeline {__version__}' + (f' {commit}' if commit else '')


# ---------------------------------------------------------------- the command line
def build_parser(fm):
    ap = argparse.ArgumentParser(prog='align_pipeline', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    genomes = fm.returnOptions('Genomes')
    sd = sub.add_parser('setup-database', help="add the ledger's columns to the header of the genome's AlignmentDatabase worksheet (once), then exit")
    sd.add_argument('Genome', choices=genomes)
    r = sub.add_parser('run', help='align samples')
    r.add_argument('Genome', choices=genomes, help='Version of the genome to align to')
    r.add_argument('-s', '--SampleIDs', nargs='+', metavar='', choices=fm.returnOptions('Samples'), help='Restrict analysis to the listed sampleIDs')
    r.add_argument('-c', '--Species', nargs='+', metavar='', choices=fm.returnOptions('Species'), help='Restrict analysis to the following species')
    r.add_argument('-p', '--ProjectIDs', nargs='+', metavar='', choices=fm.returnOptions('ProjectIDs'), help='Restrict analysis to specific ProjectIDs')
    r.add_argument('-e', '--Ecogroups', nargs='+', metavar='', choices=fm.returnOptions('Ecogroups'), help='Restrict analysis to a specific Ecogroup')
    r.add_argument('-u', '--Subgroups', nargs='+', metavar='', choices=fm.returnOptions('Subgroups'), help='Restrict analysis to a specific Subgroup')
    r.add_argument('-r', '--Rerun', action='store_true', help='Force realignment of samples that are already in the alignment database (their old row is NOT replaced; remove it first)')
    r.add_argument('--batch-size', type=int, default=None, help='Process only the first N of the selected samples (sorted by SampleID)')
    r.add_argument('--max-cores', type=int, default=None, help='Cores the whole pipeline may use at any time (default: all but 8). Lower it on a shared machine')
    r.add_argument('--max-ram-gb', type=float, default=None, help='RAM the whole pipeline may reserve (default: 70%% of what is available now)')
    r.add_argument('--samples-in-flight', type=int, default=6, help='Samples processed at the same time (default 6)')
    r.add_argument('--align-threads', type=int, default=16, help='minimap2 threads per sample (default 16)')
    r.add_argument('--prefetch', type=int, default=2, help='Samples whose uBAMs may be downloaded beyond those in flight (default 2). Bounds the disk used by uBAMs')
    r.add_argument('-d', '--MaxDownloads', type=int, default=4, help='Files downloaded at the same time (default 4)')
    r.add_argument('--uploads', type=int, default=2, help='Samples uploaded at the same time (default 2)')
    r.add_argument('--markdup-heap-gb', type=int, default=16, help='Java heap of MarkDuplicates (default 16)')
    r.add_argument('--tmp-dir', default=None, help='Base folder (ideally a fast local disk) for the GATK/Picard scratch files: one subfolder per sample')
    r.add_argument('--min-free-gb', type=float, default=200, help='The download queue pauses unless the disk keeps at least this much free after reserving room for the next sample (default 200)')
    r.add_argument('--disk-reserve-factor', type=float, default=3.5, help='Disk reserved per sample, as a multiple of its uBAM size (default 3.5)')
    r.add_argument('--max-pending-uploads', type=int, default=None, help='Finished samples that may wait locally for an upload slot (default: --uploads + 2)')
    r.add_argument('--upload-min-mb-s', type=float, default=5.0, help='The slowest upload speed that is still believable; an upload slower than this for its whole size is killed as hung (default 5)')
    r.add_argument('--status-minutes', type=float, default=10, help='Print a status line this often; 0 = never')
    r.add_argument('--include-excluded', action='store_true', help='Also run samples whose Ecogroup is Error or Exclude (skipped by default)')
    r.add_argument('--keep-local', action='store_true', help='Keep the local BAMs and temp files after the upload is verified (default: delete them)')
    r.add_argument('--slices', default=None, help='Region slice table (default: slice_regions/<Genome>.tsv of this module when it exists)')
    r.add_argument('--no-slices', action='store_true', help='Do not cut region slices')
    r.add_argument('--no-panel', action='store_true', help='Do not run module 04\'s panel measurements (without this the run first checks that the panel is ready and stops if it is not)')
    r.add_argument('--panel-python', default=None, help='The interpreter that can import sd_genotypes (default: the one running this pipeline)')
    r.add_argument('--panel-pythonpath', default=None, help='Where sd_genotypes is checked out, only if `source tools/env.sh` was not run (it puts every module on PYTHONPATH, which the panel inherits)')
    r.add_argument('--panel-modules', default='all', help='Which panel modules to measure (default all)')
    r.add_argument('--panel-expect', nargs='+', metavar='KEY', default=None, help='Definition keys that `panel check` must print, e.g. the pilot\'s frozen keys; the run does not start if one is missing')
    r.add_argument('--no-verify-counts', action='store_true', help='Skip counting the reads in the uBAMs (saves ~2 min per sample, loses the read-count cross-check)')
    r.add_argument('--no-upload', action='store_true', help='Test mode: do everything locally, upload nothing, write no ledger row, delete nothing')
    r.add_argument('--nice', type=int, default=5, help='nice value for every command (default 5)')
    r.add_argument('--list', action='store_true', help='Only list the samples that would be run and the data they need, then exit')
    r.add_argument('--print-commands', metavar='SAMPLE', default=None, help='Print the commands that would be run for this sample, then exit')
    return ap


def main(argv=None, fm=None, which=shutil.which, run=subprocess.run, preflight=pn.preflight, out=print, pysam_version=None):
    if fm is None:
        from cichlid_common.file_manager import FileManager
        fm = FileManager()
    args = build_parser(fm).parse_args(argv)
    fm.setGenome(args.Genome)
    if args.cmd == 'setup-database':
        added = adb.AppendingDatabaseWriter(fm).setup_columns()
        out(('Added to the header of the ' + args.Genome + ' worksheet: ' + ', '.join(added)) if added else 'The header already has every column the pipeline writes. Nothing to do.')
        return 0
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        resource.setrlimit(resource.RLIMIT_NOFILE, (min(4096, hard) if hard != resource.RLIM_INFINITY else 4096, hard))
    except (ValueError, OSError):
        pass
    configure_rclone_env()
    fm.setSamples(projectIDs=args.ProjectIDs, sampleIDs=args.SampleIDs, species=args.Species, ecogroups=args.Ecogroups, subgroups=args.Subgroups, rerun=args.Rerun)
    samples = sorted(fm.samples)
    if not args.include_excluded:
        flagged = set(fm.merged_dt[fm.merged_dt.Ecogroup.isin(EXCLUDED_ECOGROUPS)].SampleID)
        skipped = [s for s in samples if s in flagged]
        if skipped:
            out(f'Skipping {len(skipped)} sample(s) with Ecogroup Error/Exclude (use --include-excluded to run them): ' + ','.join(skipped))
        samples = [s for s in samples if s not in flagged]
    if args.batch_size:
        samples = samples[:args.batch_size]
    if not samples:
        out('No samples to run')
        return 1
    fm.readSampleDatabase()
    fm.readAlignmentDatabase()
    ubam_bytes = fm.merged_dt.groupby('SampleID').FileSize.sum().to_dict() if 'FileSize' in fm.merged_dt else {}          # uBAM bytes per sample, for the disk ledger
    paths = [sp.SamplePaths.from_fm(fm, s, scratch_base=args.tmp_dir, need_bytes=(ubam_bytes.get(s) if ubam_bytes.get(s) == ubam_bytes.get(s) else None)) for s in samples]
    ubam_n = sum(len(p.ubams) for p in paths)
    if args.list:
        out(f'{len(samples)} samples, {ubam_n} uBAM files, {sum(p.need_bytes for p in paths) / 1e9:.0f} GB of uBAMs:')
        out(','.join(samples))
        return 0

    # ---- reference, slices, panel
    ref = fm.localGenomeFile
    index = fm.localMinimapGenomeFile if hasattr(fm, 'localMinimapGenomeFile') and os.path.exists(fm.localMinimapGenomeFile) else ref
    fai = ref + '.fai'
    slice_file = None if args.no_slices else (args.slices or os.path.join(MODULE_DIR, 'slice_regions', args.Genome + '.tsv'))
    slices = []
    if slice_file and os.path.exists(slice_file):
        try:
            slices = sl.load_slices(slice_file, fai=fai if os.path.exists(fai) else None)
        except ValueError as e:
            out(f'The slice table cannot be read: {e}')
            return 1
        problems = sl.check_against_fai(slices, fai) if os.path.exists(fai) else []
        if problems:
            out('Slice table does not fit the reference:\n  ' + '\n  '.join(problems))
            return 1
    panel_cfg = None if args.no_panel else pn.PanelConfig(python=args.panel_python, modules=args.panel_modules, pythonpath=args.panel_pythonpath, expect=args.panel_expect)
    versions_for_plan = {}
    meta = sample_meta(fm.sample_dt, fm.reads_dt, samples)
    cfg = sp.Config(max_cores=args.max_cores, max_ram_gb=args.max_ram_gb, samples_in_flight=args.samples_in_flight, align_threads=args.align_threads, prefetch=args.prefetch,
                    max_downloads=args.MaxDownloads, uploads=args.uploads, keep_local=args.keep_local, slices_file=slice_file, nice=args.nice, verify_counts=not args.no_verify_counts,
                    upload=not args.no_upload, markdup_heap_gb=args.markdup_heap_gb, scratch_base=args.tmp_dir, min_free_gb=args.min_free_gb, disk_reserve_factor=args.disk_reserve_factor,
                    max_pending_uploads=args.max_pending_uploads, status_s=int(args.status_minutes * 60), panel=panel_cfg, pipeline_version=pipeline_version(run=run),
                    upload_min_mb_s=args.upload_min_mb_s)

    def upload_cmds(p):
        copy_cmd = fm.uploadData(p.bam_dir, parallel=True)                       # ['rclone', 'copy', local, cloud]
        exclude = ['--exclude', '*.partial*']
        return copy_cmd + exclude, ['rclone', 'check', copy_cmd[2], copy_cmd[3], '--one-way', '--size-only'] + exclude

    if args.print_commands:                                                      # nothing is executed, nothing is downloaded
        chosen = [p for p in paths if p.sample == args.print_commands]
        if not chosen:
            out(f'{args.print_commands} is not among the selected samples')
            return 1
        pipeline = sp.SeqPipeline(cfg, sp.Runner(sp.ResourceBudget(cfg.max_cores, cfg.max_ram_gb)), ref, index, 0, None, upload_cmds, None, slices=slices, genome_version=fm.genome_version)
        for name, cmd in pipeline.plan(chosen[0]):
            out(f'--- {name}\n{sp.sh_join(cmd)}\n')
        return 0

    # ---- tools, and the panel's preflight: a run that cannot measure does not start (before hours of work)
    for tool in ('gatk', 'minimap2', 'samtools', 'rclone'):
        if which(tool) is None:
            out(f'{tool} not found on PATH (activate the conda environment)')
            return 1
    if panel_cfg:
        ok, text = preflight(panel_cfg)
        if not ok:
            out('The panel is not ready, so nothing was started:\n  ' + text.replace('\n', '\n  ') + '\nFix it, or run with --no-panel.')
            return 1
        out('panel ready; definition keys:\n  ' + text.replace('\n', '\n  '))
    versions = tool_versions(run=run, pysam_version=pysam_version)
    for f in (ref + '.fai', ref.rsplit('.', 1)[0] + '.dict'):
        if not os.path.exists(f):
            out('WARNING: missing reference index file ' + f)
    if os.path.exists(fai):
        ref_size = sum(int(l.split('\t')[1]) for l in open(fai))
    else:
        import pysam
        ref_size = sum(pysam.FastaFile(ref).lengths)

    fm.uploadData(fm.localGenomeDir)                                             # genome data, exactly as before
    fm.downloadData(fm.localGenomeDir)
    os.makedirs(fm.localErrorsDir, exist_ok=True)
    log_path = fm.localErrorsDir + 'seq_pipeline_' + time.strftime('%Y%m%d-%H%M%S') + '.log'

    def log_line(line):
        print(line, flush=True)
        with open(log_path, 'a') as f:
            f.write(line + '\n')

    read_lengths = fm.merged_dt.groupby('SampleID').ReadLength.first().to_dict()
    db_writer = adb.AppendingDatabaseWriter(fm, log=log_line)
    ledger = DiskLedger(fm.localMasterDir, cfg.disk_reserve_factor, cfg.min_free_gb)
    runner = sp.Runner(sp.ResourceBudget(cfg.max_cores, cfg.max_ram_gb), nice=cfg.nice)
    pipeline = sp.SeqPipeline(cfg, runner, ref, index, ref_size, lambda job: fm.downloadData(job.path, parallel=True), upload_cmds, db_writer.add_row, run_ids={s: m['RunIDs'] for s, m in meta.items()},
                              read_lengths=read_lengths, versions=versions, genome_version=fm.genome_version, log=log_line, slices=slices, ledger=ledger, sample_meta=meta)
    if cfg.upload:                                                               # fail BEFORE hours of work if the worksheet cannot take the columns
        lacking = db_writer.check_columns()
        if lacking:
            out(f'The {args.Genome} worksheet of the AlignmentDatabase has no column(s) {lacking}. Run once, from ONE machine:\n  python -m align_pipeline setup-database {args.Genome}')
            return 1
    log_line(f'{len(paths)} samples, {ubam_n} uBAM files; at most {cfg.max_cores} cores and {cfg.max_ram_gb:.0f} GB RAM, {cfg.samples_in_flight} samples in flight, {cfg.align_threads} alignment threads each')
    log_line(f'tools: {versions}; pipeline {cfg.pipeline_version}; panel: ' + ('off' if not panel_cfg else f'{panel_cfg.python} -m sd_genotypes --modules {panel_cfg.modules}'))
    log_line(f'log file: {log_path}')
    if slices:
        log_line(f'{len(slices)} region slices per sample from {slice_file}: ' + ', '.join(s.name for s in slices))
    text, done, failed = pipeline.run_batch(paths)
    log_line(text)
    return 1 if failed else 0
