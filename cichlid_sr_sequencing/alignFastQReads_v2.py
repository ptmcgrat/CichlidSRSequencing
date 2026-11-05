import argparse, os, resource, shutil, subprocess, sys, time
from helper_modules.file_manager import FileManager as FM
from helper_modules import bam_slices as sl
from helper_modules import seq_pipeline as sp
from helper_modules.alignment_db import AppendingDatabaseWriter
from helper_modules.disk_ledger import DiskLedger
from helper_modules.download_queue import configure_rclone_env

HERE = os.path.dirname(os.path.abspath(__file__))
fm_obj = FM()
parser = argparse.ArgumentParser(usage = 'Pipelined download -> alignment -> duplicate marking -> split BAMs -> statistics -> region slices -> verified upload -> database, for many samples at once. No GVCFs.')
parser.add_argument('Genome', type = str, choices = fm_obj.returnOptions('Genomes'), help = 'Version of the genome to align to')
parser.add_argument('-s', '--SampleIDs', nargs = '+', metavar = '', choices = fm_obj.returnOptions('Samples'), help = 'Restrict analysis to the listed sampleIDs')
parser.add_argument('-c', '--Species', nargs = '+', metavar = '', choices = fm_obj.returnOptions('Species'), help = 'Restrict analysis to the following species')
parser.add_argument('-p', '--ProjectIDs', nargs = '+', metavar = '', choices = fm_obj.returnOptions('ProjectIDs'), help = 'Restrict analysis to specific ProjectIDs')
parser.add_argument('-e', '--Ecogroups', nargs = '+', metavar = '', choices = fm_obj.returnOptions('Ecogroups'), help = 'Restrict analysis to a specific Ecogroup')
parser.add_argument('-u', '--Subgroups', nargs = '+', metavar = '', choices = fm_obj.returnOptions('Subgroups'), help = 'Restrict analysis to a specific Subgroup')
parser.add_argument('-r', '--Rerun', action = 'store_true', help = 'Force realignment of samples that are already in the alignment database (their old row is NOT replaced; remove it first)')
parser.add_argument('--batch-size', type = int, default = None, help = 'Process only the first N of the selected samples (sorted by SampleID)')
parser.add_argument('--max-cores', type = int, default = None, help = 'Cores the whole pipeline may use at any time (default: all but 8). Lower it on a shared machine')
parser.add_argument('--max-ram-gb', type = float, default = None, help = 'RAM the whole pipeline may reserve (default: 70%% of what is free now)')
parser.add_argument('--samples-in-flight', type = int, default = 6, help = 'Samples processed at the same time (default 6)')
parser.add_argument('--align-threads', type = int, default = 16, help = 'minimap2 threads per sample (default 16)')
parser.add_argument('--prefetch', type = int, default = 2, help = 'Samples whose uBAMs may be downloaded beyond those in flight (default 2). Bounds the disk used by uBAMs')
parser.add_argument('-d', '--MaxDownloads', type = int, default = 4, help = 'Files downloaded at the same time (default 4)')
parser.add_argument('--uploads', type = int, default = 2, help = 'Samples uploaded at the same time (default 2)')
parser.add_argument('--markdup-heap-gb', type = int, default = 16, help = 'Java heap of MarkDuplicates (default 16). More means fewer scratch files; try 32 on the test run')
parser.add_argument('--tmp-dir', default = None, help = 'Base folder (ideally a fast local disk) for the GATK/Picard scratch files: one subfolder per sample; default: the sample temp folder')
parser.add_argument('--min-free-gb', type = float, default = 200, help = 'The download queue pauses unless the disk keeps at least this much free after reserving room for the next sample (default 200)')
parser.add_argument('--disk-reserve-factor', type = float, default = 3.5, help = 'Disk reserved per sample, as a multiple of its uBAM size (default 3.5; the summary reports the peak actually observed)')
parser.add_argument('--max-pending-uploads', type = int, default = None, help = 'Finished samples that may wait locally for an upload slot (default: --uploads + 2); beyond that the compute side waits')
parser.add_argument('--status-minutes', type = float, default = 10, help = 'Print a status line (what every sample is doing, transfers, disk, cores) this often; 0 = never')
parser.add_argument('--include-excluded', action = 'store_true', help = 'Also run samples whose Ecogroup is Error or Exclude (skipped by default)')
parser.add_argument('--setup-database', action = 'store_true', help = 'Add the new statistics columns to the header of the genome\'s AlignmentDatabase worksheet (once, before the first run), then exit')
parser.add_argument('--keep-local', action = 'store_true', help = 'Keep the local BAMs and temp files after the upload is verified (default: delete them)')
parser.add_argument('--slices', default = None, help = 'Region slice table (default: slice_regions/<Genome>.tsv when it exists)')
parser.add_argument('--no-slices', action = 'store_true', help = 'Do not cut region slices')
parser.add_argument('--no-verify-counts', action = 'store_true', help = 'Skip counting the reads in the uBAMs (saves ~2 min per sample, loses the read-count cross-check)')
parser.add_argument('--no-upload', action = 'store_true', help = 'Test mode: do everything locally, upload nothing, write no database row, delete nothing')
parser.add_argument('--nice', type = int, default = 5, help = 'nice value for every command (default 5)')
parser.add_argument('--list', action = 'store_true', help = 'Only list the samples that would be run and the data they need, then exit')
parser.add_argument('--print-commands', metavar = 'SAMPLE', default = None, help = 'Print the commands that would be run for this sample, then exit')
args = parser.parse_args()

soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
resource.setrlimit(resource.RLIMIT_NOFILE, (4096, hard))
configure_rclone_env()

fm_obj.setGenome(args.Genome)
if args.setup_database:
	added = AppendingDatabaseWriter(fm_obj).setup_columns()
	print(('Added to the header of the ' + args.Genome + ' worksheet: ' + ', '.join(added)) if added else 'The header already has every column the pipeline writes. Nothing to do.')
	sys.exit(0)
fm_obj.setSamples(projectIDs = args.ProjectIDs, sampleIDs = args.SampleIDs, species = args.Species, ecogroups = args.Ecogroups, subgroups = args.Subgroups, rerun = args.Rerun)
samples = sorted(fm_obj.samples)
EXCLUDED_ECOGROUPS = ('Error', 'Exclude')                                       # samples flagged on purpose: not run unless --include-excluded
if not args.include_excluded:
	flagged = set(fm_obj.merged_dt[fm_obj.merged_dt.Ecogroup.isin(EXCLUDED_ECOGROUPS)].SampleID)
	skipped = [s for s in samples if s in flagged]
	if skipped:
		print(f'Skipping {len(skipped)} sample(s) with Ecogroup Error/Exclude (use --include-excluded to run them): ' + ','.join(skipped))
	samples = [s for s in samples if s not in flagged]
if args.batch_size:
	samples = samples[:args.batch_size]
if not samples:
	sys.exit('No samples to run')

# genome data, exactly as before
fm_obj.uploadData(fm_obj.localGenomeDir)
fm_obj.downloadData(fm_obj.localGenomeDir)
ref = fm_obj.localGenomeFile
index = fm_obj.localMinimapGenomeFile if hasattr(fm_obj, 'localMinimapGenomeFile') and os.path.exists(fm_obj.localMinimapGenomeFile) else ref

fm_obj.readSampleDatabase()
fm_obj.readAlignmentDatabase()
ubam_bytes = fm_obj.merged_dt.groupby('SampleID').FileSize.sum().to_dict() if 'FileSize' in fm_obj.merged_dt else {}          # uBAM bytes per sample, for the disk ledger
paths = [sp.SamplePaths.from_fm(fm_obj, s, scratch_base = args.tmp_dir, need_bytes = (ubam_bytes.get(s) if ubam_bytes.get(s) == ubam_bytes.get(s) else None)) for s in samples]
ubam_n = sum(len(p.ubams) for p in paths)
if args.list:
	print(f'{len(samples)} samples, {ubam_n} uBAM files, {sum(p.need_bytes for p in paths) / 1e9:.0f} GB of uBAMs:'); print(','.join(samples)); sys.exit(0)

# tools
for tool in ('gatk', 'minimap2', 'samtools', 'rclone'):
	if shutil.which(tool) is None:
		sys.exit(f'{tool} not found on PATH (activate the conda environment)')
def first_line(cmd):
	out = subprocess.run(cmd, capture_output = True, encoding = 'utf-8')
	return (out.stdout.strip() or out.stderr.strip()).split('\n')[0]
versions = {'minimap2': first_line(['minimap2', '--version']), 'samtools': first_line(['samtools', '--version']), 'gatk': first_line(['gatk', '--version'])}
try:
	import pysam; versions['pysam'] = pysam.__version__
except ImportError:
	sys.exit('pysam is needed by the split step')
for f in (ref + '.fai', ref.rsplit('.', 1)[0] + '.dict'):
	if not os.path.exists(f):
		print('WARNING: missing reference index file ' + f)
ref_size = sum(int(l.split('\t')[1]) for l in open(ref + '.fai')) if os.path.exists(ref + '.fai') else sum(pysam.FastaFile(ref).lengths)

slice_file = None if args.no_slices else (args.slices or os.path.join(HERE, 'slice_regions', args.Genome + '.tsv'))
slices = sl.load_slices(slice_file) if slice_file and os.path.exists(slice_file) else []
if slices and os.path.exists(ref + '.fai'):
	problems = sl.check_against_fai(slices, ref + '.fai')
	if problems:
		sys.exit('Slice table does not fit the reference:\n  ' + '\n  '.join(problems))

merged = fm_obj.merged_dt
read_lengths = merged.groupby('SampleID').ReadLength.first().to_dict()
run_ids = {s: ',,'.join(list(fm_obj.reads_dt[fm_obj.reads_dt.SampleID == s].RunID)) for s in samples}

log_path = fm_obj.localErrorsDir + 'seq_pipeline_' + time.strftime('%Y%m%d-%H%M%S') + '.log'
def log_line(line):
	print(line, flush = True)
	with open(log_path, 'a') as f:
		f.write(line + '\n')

cfg = sp.Config(max_cores = args.max_cores, max_ram_gb = args.max_ram_gb, samples_in_flight = args.samples_in_flight, align_threads = args.align_threads, prefetch = args.prefetch,
                max_downloads = args.MaxDownloads, uploads = args.uploads, keep_local = args.keep_local, slices_file = slice_file, nice = args.nice, verify_counts = not args.no_verify_counts,
                upload = not args.no_upload, markdup_heap_gb = args.markdup_heap_gb, scratch_base = args.tmp_dir, min_free_gb = args.min_free_gb, disk_reserve_factor = args.disk_reserve_factor,
                max_pending_uploads = args.max_pending_uploads, status_s = int(args.status_minutes * 60))
db_writer = AppendingDatabaseWriter(fm_obj, log = log_line)
ledger = DiskLedger(fm_obj.localMasterDir, cfg.disk_reserve_factor, cfg.min_free_gb)
budget = sp.ResourceBudget(cfg.max_cores, cfg.max_ram_gb)
runner = sp.Runner(budget, nice = cfg.nice)

def upload_cmds(p):
	copy_cmd = fm_obj.uploadData(p.bam_dir, parallel = True)                    # ['rclone', 'copy', local, cloud]
	exclude = ['--exclude', '*.partial*']
	return copy_cmd + exclude, ['rclone', 'check', copy_cmd[2], copy_cmd[3], '--one-way', '--size-only'] + exclude

pipeline = sp.SeqPipeline(cfg, runner, ref, index, ref_size, lambda job: fm_obj.downloadData(job.path, parallel = True), upload_cmds, db_writer.add_row,
                          run_ids = run_ids, read_lengths = read_lengths, versions = versions, genome_version = fm_obj.genome_version, log = log_line, slices = slices, ledger = ledger)

if args.print_commands:
	chosen = [p for p in paths if p.sample == args.print_commands]
	if not chosen:
		sys.exit(f'{args.print_commands} is not among the selected samples')
	for name, cmd in pipeline.plan(chosen[0]):
		print(f'--- {name}\n{sp.sh_join(cmd)}\n')
	sys.exit(0)

if cfg.upload:                                                                  # fail BEFORE hours of work if the worksheet cannot take the new columns
	lacking = db_writer.check_columns()
	if lacking:
		sys.exit(f'The {args.Genome} worksheet of the AlignmentDatabase has no column(s) {lacking}. Run once, from ONE machine:\n  python3 alignFastQReads_v2.py {args.Genome} --setup-database')
log_line(f'{len(paths)} samples, {ubam_n} uBAM files; at most {cfg.max_cores} cores and {cfg.max_ram_gb:.0f} GB RAM, {cfg.samples_in_flight} samples in flight, {cfg.align_threads} alignment threads each')
log_line(f'tools: {versions}')
log_line(f'disk: ledger on {fm_obj.localMasterDir}, {ledger.free() / 1e9:.0f} GB free now; a sample is only downloaded when {cfg.disk_reserve_factor} x its uBAM size + {cfg.min_free_gb:.0f} GB margin fit; at most {cfg.disk_slots} samples have uBAMs on disk')
log_line(f'log file: {log_path}')
if slices:
	log_line(f'{len(slices)} region slices per sample from {slice_file}')
text, done, failed = pipeline.run_batch(paths)
log_line(text)
sys.exit(1 if failed else 0)
