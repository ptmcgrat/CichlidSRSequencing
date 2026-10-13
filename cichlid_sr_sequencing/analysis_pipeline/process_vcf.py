"""
process_vcf.py  (Oct 2026 rewrite for the per-analysis folder layout)

Takes the concatenated VCF that callVariants.py wrote for an analysis and produces the PASS-variants VCF used by pca_maker.py.
Concatenation and indexing now happen in callVariants.py (--concat_and_index), so the old --merge / --compress_and_index
options are gone.

Files (all located through nikesh_file_manager.createAnalysisIDFiles):
    input   Outputs/<genome>/<analysis>/vcf_concat_output/<analysis>.vcf.gz
    output  Outputs/<genome>/<analysis>/vcf_filtering/<analysis>_filtered.vcf.gz          every variant, FILTER column filled in
            Outputs/<genome>/<analysis>/vcf_concat_output/<analysis>_pass_variants.vcf.gz  PASS biallelic SNPs (+ .tbi)
            Outputs/<genome>/<analysis>/filtering_stats/                                   vcftools stats, filter summary, per-sample missingness
    logs    logs/<date_time>_process_vcf_<genome>_<analysis>.log / .err   (next to this script, one pair per run)

Usage:
    python3 process_vcf.py Mzebra_GT3_NCBI --analysis_id new_aulonocara --vcf_stats      # optional: stats on the raw VCF
    python3 process_vcf.py Mzebra_GT3_NCBI --analysis_id new_aulonocara --filter
"""
import argparse, os, sys, re, subprocess, datetime, platform, atexit, statistics
from helper_modules.nikesh_file_manager import FileManager as FM

parser = argparse.ArgumentParser(usage = 'python3 process_vcf.py <genome> --analysis_id <name> [--vcf_stats] [--filter]')
parser.add_argument('genome', help = 'genome version used by callVariants.py (e.g. Mzebra_GT3_NCBI)')
parser.add_argument('--analysis_id', required = True, help = 'analysis name used by callVariants.py (e.g. new_aulonocara)')
parser.add_argument('--input_vcf', help = 'use this VCF instead of Outputs/<genome>/<analysis>/vcf_concat_output/<analysis>.vcf.gz (e.g. a _2.vcf.gz re-run)')
parser.add_argument('-v', '--vcf_stats', action = 'store_true', help = 'run vcftools stats (depth, quality, missingness, het, allele freq) on the input VCF')
parser.add_argument('-f', '--filter', action = 'store_true', help = 'run GATK VariantFiltration + SelectVariants to make the PASS-variants VCF')
parser.add_argument('--min_dp', type = int, help = 'site depth (INFO/DP) lower cut-off. Default: half the median site depth of the input VCF')
parser.add_argument('--max_dp', type = int, help = 'site depth (INFO/DP) upper cut-off. Default: twice the median site depth of the input VCF')
parser.add_argument('--overwrite', action = 'store_true', help = 'replace existing filtered / pass-variants files for this analysis')
args = parser.parse_args()
if not (args.vcf_stats or args.filter):
    parser.error('nothing to do: give --vcf_stats and/or --filter')

# Hard filters (GATK JEXL on INFO fields). Same as the lab's previous filtering except the depth cut-offs, which now scale with
# the cohort (the old fixed DP 85-300 was tuned for a much smaller cohort) and the no-call cut-off, which scales with sample number.
MIN_AF = 0.05               # allele frequency across all samples in the VCF
MIN_INBREEDING = -0.6       # strongly negative = excess heterozygosity (collapsed paralogs)
MIN_QD = 2.0                # quality by depth
MAX_FS = 40.0               # Fisher strand bias
MIN_MQ = 50.0               # mapping quality
MAX_NOCALL_FRACTION = 0.10  # NCC (number of no-called samples) may be at most this fraction of samples


def start_run_logs():
    # One .log (stdout) / .err (stderr) pair per run in logs/ next to this script, at the file-descriptor level so gatk,
    # bcftools and vcftools output is captured too. Same approach as callVariants.py.
    log_dir = os.path.dirname(os.path.abspath(__file__)) + '/logs/'
    os.makedirs(log_dir, exist_ok = True)
    start = datetime.datetime.now()
    base = log_dir + f"{start:%Y-%m-%d_%H-%M-%S}_process_vcf_{args.genome}_{re.sub(r'[^A-Za-z0-9._-]+', '_', args.analysis_id)}"
    print(f"process_vcf.py started {start:%Y-%m-%d %H:%M:%S}. Output is being written to:\n  stdout: {base}.log\n  stderr: {base}.err")
    print(f"Follow progress with: tail -f {base}.log")
    sys.stdout.flush(); sys.stderr.flush()
    log_fh, err_fh = open(base + '.log', 'w'), open(base + '.err', 'w')
    os.dup2(log_fh.fileno(), sys.stdout.fileno()); os.dup2(err_fh.fileno(), sys.stderr.fileno())
    sys.stdout.reconfigure(line_buffering = True); sys.stderr.reconfigure(line_buffering = True)
    header = (f"==== process_vcf.py run started {start:%Y-%m-%d %H:%M:%S}\ncommand : python3 {' '.join(sys.argv)}\n"
              f"host    : {platform.node()} (user {os.getenv('USER')})\ncwd     : {os.getcwd()}\n")
    print(header); print(header, file = sys.stderr)
    state = {'done': False}
    def finish():
        end = datetime.datetime.now()
        line = f"\n==== {'RUN COMPLETE' if state['done'] else 'STOPPED BEFORE COMPLETION (see the .err file)'} | finished {end:%Y-%m-%d %H:%M:%S} | elapsed {str(end - start).split('.')[0]}"
        print(line); print(line, file = sys.stderr)
    atexit.register(finish)
    return state


def run(command, **kwargs):
    # Run a command, stop the script if it fails
    print(f"[{datetime.datetime.now():%H:%M:%S}] " + ' '.join(command))
    result = subprocess.run(command, **kwargs)
    if result.returncode != 0:
        sys.exit(f"Command failed (exit {result.returncode}): {' '.join(command)}")
    return result


class VCFProcessor:
    def __init__(self):
        self.fm = FM(args.genome)
        self.fm.createAnalysisIDFiles(args.analysis_id)
        self.input_vcf = args.input_vcf or self.fm.localAnalysisVCFFile
        if not os.path.exists(self.input_vcf) or not os.path.exists(self.input_vcf + '.tbi'):
            sys.exit(f"Input VCF or its .tbi index not found: {self.input_vcf}")
        if not os.path.exists(self.fm.localGenomeFile):
            sys.exit(f"Reference fasta not found: {self.fm.localGenomeFile} (run make_intervals.py {args.genome} on this server first)")
        for directory in [self.fm.localAnalysisFilteringDir, self.fm.localAnalysisStatsDir]:
            os.makedirs(directory, exist_ok = True)
        self.samples = run(['bcftools', 'query', '-l', self.input_vcf], capture_output = True, encoding = 'utf-8').stdout.split()
        self.prefix = self.fm.localAnalysisStatsDir + args.analysis_id
        print(f"Input VCF: {self.input_vcf} ({len(self.samples)} samples, {self.record_count(self.input_vcf):,} variant records)")

    def record_count(self, vcf):
        # number of records, read from the .tbi index (instant)
        out = subprocess.run(['bcftools', 'index', '-n', vcf], capture_output = True, encoding = 'utf-8').stdout.strip()
        return int(out) if out.isdigit() else -1

    def largest_contig(self):
        # contig with the most records (from the index), used for quick per-contig summaries instead of reading the whole
        # genome. Falls back to the first contig if the index has no record counts.
        out = subprocess.run(['bcftools', 'index', '-s', self.input_vcf], capture_output = True, encoding = 'utf-8')
        rows = [l.split('\t') for l in out.stdout.splitlines() if l.strip()]
        if out.returncode == 0 and rows:
            return max(rows, key = lambda f: int(f[2]))[0]
        return run(['tabix', '-l', self.input_vcf], capture_output = True, encoding = 'utf-8').stdout.split()[0]

    def generate_stats(self):
        # vcftools summaries of the input VCF, run in parallel (each reads the whole file)
        print('Generating vcftools stats for the input VCF ...')
        base = ['vcftools', '--gzvcf', self.input_vcf, '--out', self.prefix]
        stats = [['--freq2', '--max-alleles', '2'], ['--depth'], ['--site-mean-depth'], ['--site-quality'], ['--missing-indv'], ['--missing-site'], ['--het']]
        procs = [subprocess.Popen(base + s) for s in stats]
        failed = [' '.join(s) for s, p in zip(stats, procs) if p.wait() != 0]
        if failed:
            sys.exit('vcftools failed for: ' + ', '.join(failed))
        print(f"vcftools stats written to {self.fm.localAnalysisStatsDir}")

    def depth_cutoffs(self):
        # Site depth (INFO/DP) is summed over all samples, so sensible cut-offs depend on cohort size and coverage.
        # Default: half and twice the median site depth, measured on the contig with the most records.
        if args.min_dp and args.max_dp:
            return args.min_dp, args.max_dp
        contig = self.largest_contig()
        values = [int(v) for v in run(['bcftools', 'query', '-r', contig, '-f', '%INFO/DP\n', self.input_vcf], capture_output = True, encoding = 'utf-8').stdout.split() if v.isdigit()]
        median = statistics.median(values)
        min_dp, max_dp = args.min_dp or int(median / 2), args.max_dp or int(median * 2)
        print(f"Median site depth (INFO/DP) on {contig}: {median:,.0f} over {len(values):,} sites -> min_DP {min_dp}, max_DP {max_dp}")
        return min_dp, max_dp

    def filter_variants(self):
        filtered, passed = self.fm.localAnalysisFilteredVCFFile, self.fm.localAnalysisPassVCFFile
        if (os.path.exists(filtered) or os.path.exists(passed)) and not args.overwrite:
            sys.exit(f"Filtered/pass files already exist for {args.analysis_id}. Use --overwrite to replace them:\n  {filtered}\n  {passed}")
        min_dp, max_dp = self.depth_cutoffs()
        max_ncc = round(MAX_NOCALL_FRACTION * len(self.samples))
        filters = [('allele_freq', f'AF < {MIN_AF}'), ('inbreeding_test', f'InbreedingCoeff < {MIN_INBREEDING}'), ('depth_Qual', f'QD < {MIN_QD}'),
                   ('max_DP', f'DP > {max_dp}'), ('min_DP', f'DP < {min_dp}'), ('strand_bias', f'FS > {MAX_FS}'),
                   ('mapping_quality', f'MQ < {MIN_MQ}'), ('no_calls', f'NCC > {max_ncc}')]
        print('Filters:\n' + '\n'.join(f"  {name:16s} {expr}" for name, expr in filters))

        print('RUNNING VARIANT FILTRATION ...')
        command = ['gatk', 'VariantFiltration', '-R', self.fm.localGenomeFile, '-V', self.input_vcf, '-O', filtered, '--verbosity', 'ERROR']
        for name, expr in filters:
            command += ['--filter-name', name, '--filter-expression', expr]
        run(command)

        # Keep PASS sites that are biallelic SNPs (indels and multiallelic sites removed; pca_maker uses biallelic SNPs only)
        print('EXTRACTING PASS BIALLELIC SNPs ...')
        run(['gatk', 'SelectVariants', '-V', filtered, '--exclude-filtered', '--select-type-to-include', 'SNP', '--restrict-alleles-to', 'BIALLELIC', '-O', passed])
        run(['tabix', '-f', '-p', 'vcf', passed]) # re-index with htslib so bcftools can read record counts from the index

        self.summarize(filtered, passed, filters)

    def summarize(self, filtered, passed, filters):
        # Record counts + which filters removed sites (on the largest contig) + per-sample missingness in the PASS file
        n_in, n_pass = self.record_count(self.input_vcf), self.record_count(passed)
        contig = self.largest_contig()
        flags = run(['bcftools', 'query', '-r', contig, '-f', '%FILTER\n', filtered], capture_output = True, encoding = 'utf-8').stdout.split()
        total = len(flags)
        summary = [f"Input records: {n_in:,}", f"PASS biallelic SNPs: {n_pass:,} ({100 * n_pass / max(n_in, 1):.1f}% of input)",
                   f"\nSites failing each filter on {contig} (a site can fail several; {total:,} sites):"]
        for name, expr in filters:
            n = sum(1 for f in flags if name in f.split(';'))
            summary.append(f"  {name:16s} {expr:28s} {n:>12,} ({100 * n / max(total, 1):.1f}%)")
        summary.append(f"  {'PASS':16s} {'':28s} {flags.count('PASS'):>12,} ({100 * flags.count('PASS') / max(total, 1):.1f}%)")

        run(['vcftools', '--gzvcf', passed, '--missing-indv', '--out', self.prefix + '_pass'])
        missing = []
        with open(self.prefix + '_pass.imiss') as fh:
            next(fh)
            for line in fh:
                fields = line.split()
                missing.append((fields[0], float(fields[4])))
        missing.sort(key = lambda m: -m[1])
        summary.append(f"\nPer-sample missing genotypes in the PASS file (highest 15 of {len(missing)}; full table: {self.prefix}_pass.imiss):")
        summary += [f"  {sample:22s} {frac:.3f}" for sample, frac in missing[:15]]

        text = '\n'.join(summary)
        print('\n' + text)
        with open(self.prefix + '_filter_summary.txt', 'w') as fh:
            fh.write(text + '\n')
        print(f"\nPASS variants file: {passed}")


if __name__ == '__main__':
    state = start_run_logs()
    processor = VCFProcessor()
    if args.vcf_stats:
        processor.generate_stats()
    if args.filter:
        processor.filter_variants()
    state['done'] = True