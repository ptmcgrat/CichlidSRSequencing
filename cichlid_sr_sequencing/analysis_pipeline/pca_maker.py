"""
pca_maker.py  (Oct 2026 rewrite for the per-analysis folder layout and the Google Sheets sample database)

Projects every sample in an analysis's PASS-variants VCF onto PCs built from a "core" set of samples, per region:
whole genome, each linkage group, and each inversion region. All samples in the VCF are used (no ecogroup filtering).

Inputs (located through nikesh_file_manager):
    VCF        Outputs/<genome>/<analysis>/vcf_concat_output/<analysis>_pass_variants.vcf.gz   (from process_vcf.py)
    core set   AnalysisSets column given by --core_set (default CORE_PCA); these samples build the PCs
    metadata   SampleDatabase (Ecogroup, Species, Subgroup, Sex, Inversion10), DNAReadsDatabase (ProjectID),
               AnalysisSets column Bionano_data (ground-truth fish)
Outputs:
    Outputs/<genome>/<analysis>/pca/<core_set>/
        core_samples.txt, all_samples.txt, run_info.txt
        vcfs/                         the VCF preprocessing files (step 1 below)
        PCA/<region>/                 plink2 files + <region>_projection_with_metadata.csv
        interactive_PCA_outputs/      <region>_PCA.html
        static_PCA_outputs/           <region>_PCA.pdf (needs the kaleido package)
    logs/<date_time>_pca_maker_<genome>_<analysis>_<core_set>.log / .err   (next to this script, one pair per run)

Method (same tools and steps as the previous pca_maker.py; "core" samples replace the old "subset" samples):
    1. bcftools view all samples (--min-ac=1 --no-update) and the core samples; bcftools +fill-tags AF on the core samples;
       drop core variants with AF < 0.05 or AF > 0.95 (was AF == 1); vcftools --positions keeps those positions for all samples.
    2. Per region: plink2 --vcf (all and core), --set-missing-var-ids @:# on the all-sample pfile, then PCA on the core samples
       (--freq counts --pca allele-wts --max-alleles 2). Whole genome and linkage groups are LD-pruned first
       (--indep-pairwise 50 5 0.1); inversion regions are not.
       New: with fewer than 50 core samples, LD pruning adds --bad-ld (previously the script stopped) and prints a warning.
    3. Remove .acount rows with ALT_CTS == 0 or == max(OBS_CT), then project every sample onto the PCs
       (plink2 --score <weights> 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15).
    4. Plots (plotly): Malinsky ecogroup colours for all samples; Bionano fish as large, translucent black X's;
       non-core samples outlined in black; optional --highlight_set outlined in magenta; Bionano drawn last (on top).

Usage:
    python3 pca_maker.py Mzebra_GT3_NCBI --analysis_id new_aulonocara
    python3 pca_maker.py Mzebra_GT3_NCBI --analysis_id new_aulonocara -r Inversions Whole
    python3 pca_maker.py Mzebra_GT3_NCBI --analysis_id new_aulonocara --plot_only          # re-draw plots from existing results
"""
import argparse, os, sys, re, subprocess, datetime, platform, atexit
import pandas as pd
import plotly.graph_objs as go # pip install plotly (and kaleido for PDFs)
from helper_modules.nikesh_file_manager import FileManager as FM

# Chromosome accessions per genome, from the @SQ lines of Intervals/<genome>/1.interval_list (reference order).
# M. zebra has no LG21, so the last two chromosomes are LG22 and LG23. Unplaced NW_ contigs are left out (they still count in 'Whole').
LINKAGE_GROUPS = {
    'Mzebra_GT3_NCBI': {'LG1': 'NC_135167.1', 'LG2': 'NC_135168.1', 'LG3': 'NC_135169.1', 'LG4': 'NC_135170.1', 'LG5': 'NC_135171.1',
                        'LG6': 'NC_135172.1', 'LG7': 'NC_135173.1', 'LG8': 'NC_135174.1', 'LG9': 'NC_135175.1', 'LG10': 'NC_135176.1',
                        'LG11': 'NC_135177.1', 'LG12': 'NC_135178.1', 'LG13': 'NC_135179.1', 'LG14': 'NC_135180.1', 'LG15': 'NC_135181.1',
                        'LG16': 'NC_135182.1', 'LG17': 'NC_135183.1', 'LG18': 'NC_135184.1', 'LG19': 'NC_135185.1', 'LG20': 'NC_135186.1',
                        'LG22': 'NC_135187.1', 'LG23': 'NC_135188.1'},
}

# Inversion regions from eLife Supplementary File 4: outer bp of the breakpoint intervals, i.e. Breakpoint 1 interval start to
# Breakpoint 2 interval end. Keyed by linkage group number. Mzebra_GT3 coordinates, which are identical on Mzebra_GT3_NCBI.
INVERSIONS = {
    'LG2_inv':   (2,  19751534, 43559056),
    'LG9_inv':   (9,  14912159, 32996230),
    'LG10_inv':  (10, 11459086, 29902559),
    'LG11a_inv': (11,  8281988, 18425215),
    'LG11b_inv': (11, 18425177, 30387602),
    'LG13_inv':  (13,  2398108, 23124401),
    'LG20a_inv': (20, 19616377, 29698815),
    'LG20b_inv': (20, 29682528, 33254382),
}

parser = argparse.ArgumentParser(usage = 'python3 pca_maker.py <genome> --analysis_id <name> [options]')
parser.add_argument('genome', help = 'genome version of the VCF (e.g. Mzebra_GT3_NCBI)')
parser.add_argument('--analysis_id', required = True, help = 'analysis name used by callVariants.py / process_vcf.py (e.g. new_aulonocara)')
parser.add_argument('--core_set', default = 'CORE_PCA', help = 'AnalysisSets column with the samples that build the PCs (default CORE_PCA)')
parser.add_argument('--highlight_set', help = 'optional AnalysisSets column whose samples get a magenta outline')
parser.add_argument('-r', '--regions', nargs = '*', default = ['All', 'Whole', 'Inversions'],
                    help = "regions to run: 'All' (every linkage group), 'Whole' (whole genome), 'Inversions' (all inversion regions), "
                           "or individual names such as LG7 or LG10_inv. Default: All Whole Inversions")
parser.add_argument('--input_vcf', help = 'use this VCF instead of the analysis pass-variants VCF')
parser.add_argument('--plot_only', action = 'store_true', help = 'skip all computation and re-draw plots from existing results')
parser.add_argument('--overwrite', action = 'store_true', help = 'rebuild the preprocessed VCFs even if the input and sample lists have not changed')
parser.add_argument('--threads', type = int, default = 8, help = 'compression threads for bcftools view and bgzip (plink2 uses all available cores)')
args = parser.parse_args()

# Plot styling
MALINSKY_COLORS = {'Mbuna': '#A020F0', 'AC': '#A2CD5A', 'Shallow_Benthic': '#FF6347', 'Deep_Benthic': '#4876FF',
                   'Rhamphochromis': '#8B4513', 'Diplotaxodon': '#FFA54F', 'Utaka': '#006400'}
UNKNOWN_COLOR = '#B0B0B0'     # any ecogroup missing from the palette
HIGHLIGHT_COLOR = '#FF00FF'
SIZE = {'html': {'normal': 5, 'bionano': 7}, 'pdf': {'normal': 3, 'bionano': 4}}
OUTLINE = {'html': {'noncore': 1.0, 'highlight': 2.0}, 'pdf': {'noncore': 0.5, 'highlight': 1.0}}
OPACITY = {'normal': 0.9, 'bionano': 0.6}
MIN_CORE_FOR_LD = 50          # plink2 --indep-pairwise needs this many samples unless --bad-ld


def start_run_logs():
    # One .log / .err pair per run in logs/ next to this script (file-descriptor level, so plink2/bcftools output is captured)
    log_dir = os.path.dirname(os.path.abspath(__file__)) + '/logs/'
    os.makedirs(log_dir, exist_ok = True)
    start = datetime.datetime.now()
    label = '_'.join(re.sub(r'[^A-Za-z0-9._-]+', '_', x) for x in [args.genome, args.analysis_id, args.core_set])
    base = log_dir + f"{start:%Y-%m-%d_%H-%M-%S}_pca_maker_{label}"
    print(f"pca_maker.py started {start:%Y-%m-%d %H:%M:%S}. Output is being written to:\n  stdout: {base}.log\n  stderr: {base}.err")
    print(f"Follow progress with: tail -f {base}.log")
    sys.stdout.flush(); sys.stderr.flush()
    log_fh, err_fh = open(base + '.log', 'w'), open(base + '.err', 'w')
    os.dup2(log_fh.fileno(), sys.stdout.fileno()); os.dup2(err_fh.fileno(), sys.stderr.fileno())
    sys.stdout.reconfigure(line_buffering = True); sys.stderr.reconfigure(line_buffering = True)
    header = (f"==== pca_maker.py run started {start:%Y-%m-%d %H:%M:%S}\ncommand : python3 {' '.join(sys.argv)}\n"
              f"host    : {platform.node()} (user {os.getenv('USER')})\ncwd     : {os.getcwd()}\n")
    print(header); print(header, file = sys.stderr)
    state = {'done': False}
    def finish():
        end = datetime.datetime.now()
        line = f"\n==== {'RUN COMPLETE' if state['done'] else 'STOPPED BEFORE COMPLETION (see the .err file)'} | finished {end:%Y-%m-%d %H:%M:%S} | elapsed {str(end - start).split('.')[0]}"
        print(line); print(line, file = sys.stderr)
    atexit.register(finish)
    return state


def run(command, capture = False):
    # Run a command (list, or a string run through bash with pipefail); stop the script if it fails
    shown = command if isinstance(command, str) else ' '.join(command)
    print(f"[{datetime.datetime.now():%H:%M:%S}] {shown}")
    if isinstance(command, str):
        result = subprocess.run(['bash', '-o', 'pipefail', '-c', command], capture_output = capture, encoding = 'utf-8')
    else:
        result = subprocess.run(command, capture_output = capture, encoding = 'utf-8')
    if result.returncode != 0:
        sys.exit(f"Command failed (exit {result.returncode}): {shown}" + (f"\n{result.stderr}" if capture else ''))
    return result.stdout if capture else None


def run_parallel(commands):
    # Start several commands at once (as the original did with Popen), wait for all of them, stop the script if any failed
    for command in commands:
        print(f"[{datetime.datetime.now():%H:%M:%S}] {' '.join(command)}")
    processes = [subprocess.Popen(command) for command in commands]
    failed = [' '.join(c) for c, p in zip(commands, processes) if p.wait() != 0]
    if failed:
        sys.exit('Command failed: ' + '\n'.join(failed))


class PCAMaker:
    def __init__(self):
        self.fm = FM(args.genome)
        self.fm.createAnalysisIDFiles(args.analysis_id)
        self.in_vcf = args.input_vcf or self.fm.localAnalysisPassVCFFile
        if not os.path.exists(self.in_vcf) or not os.path.exists(self.in_vcf + '.tbi'):
            sys.exit(f"Input VCF or its .tbi not found: {self.in_vcf}\n(run process_vcf.py {args.genome} --analysis_id {args.analysis_id} --filter first)")
        self.out_dir = self.fm.localAnalysisPCADir + re.sub(r'[^A-Za-z0-9._-]+', '_', args.core_set) + '/'
        self.vcf_dir = self.out_dir + 'vcfs/'
        self.html_dir = self.out_dir + 'interactive_PCA_outputs/'
        self.pdf_dir = self.out_dir + 'static_PCA_outputs/'
        for directory in [self.vcf_dir, self.out_dir + 'PCA/', self.html_dir, self.pdf_dir]:
            os.makedirs(directory, exist_ok = True)
        self.core_list = self.out_dir + 'core_samples.txt'
        self.all_list = self.out_dir + 'all_samples.txt'
        # Preprocessing files (old names in brackets)
        self.all_vcf = self.vcf_dir + 'all_samples.vcf.gz'                          # [all_ecogroup_samples.vcf.gz]
        self.core_vcf = self.vcf_dir + 'core_samples.vcf.gz'                        # [subset_samples.vcf.gz]
        self.core_af_vcf = self.vcf_dir + 'af_recalculated_core_samples.vcf.gz'     # [af_recalculated_subset_samples.vcf.gz]
        self.core_filtered_vcf = self.vcf_dir + 'af_filtered_core_samples.vcf.gz'   # [af_filtered_subset_samples.vcf.gz]
        self.positions = self.vcf_dir + 'variants_to_keep.txt'
        self.all_filtered_prefix = self.vcf_dir + 'variants_filtered_all_samples'   # [variants_filtered_ecogroup_samples]
        self.all_filtered_vcf = self.all_filtered_prefix + '.recode.vcf.gz'
        print(f"Input VCF: {self.in_vcf}\nOutput folder: {self.out_dir}")

        self._read_metadata()
        self._choose_samples()
        self._map_linkage_groups()
        self._resolve_regions()

    # ---------------------------------------------------------------- samples and metadata
    def _read_metadata(self):
        # One row per sample: Ecogroup / Species / Subgroup / Sex / Inversion10 from SampleDatabase, ProjectID(s) from DNAReadsDatabase
        self.fm.readSampleDatabase()
        self.fm.readAnalysisSets()
        meta = self.fm.sample_dt[['SampleID', 'Ecogroup', 'Species', 'Subgroup', 'Sex', 'Inversion10']].drop_duplicates('SampleID')
        projects = self.fm.reads_dt.groupby('SampleID').ProjectID.agg(lambda s: ','.join(sorted(set(s.dropna().astype(str)))))
        meta = meta.assign(ProjectID = meta.SampleID.map(projects))
        meta['Inversion10'] = meta.Inversion10.map(lambda v: '' if pd.isna(v) else str(int(v)) if float(v).is_integer() else str(v))
        self.meta = meta.set_index('SampleID')
        self.bionano = self.fm.returnAnalysisSetSamples('Bionano_data')
        self.highlight = self.fm.returnAnalysisSetSamples(args.highlight_set) if args.highlight_set else set()

    def _choose_samples(self):
        # Every sample in the VCF is projected, whatever its coverage
        self.all_samples = run(['bcftools', 'query', '-l', self.in_vcf], capture = True).split()
        core = self.fm.returnAnalysisSetSamples(args.core_set)
        not_in_vcf = sorted(core - set(self.all_samples))
        if not_in_vcf:
            sys.exit(f"{len(not_in_vcf)} {args.core_set} samples are not in the VCF, so the PCs can't be built as intended: {', '.join(not_in_vcf)}\n"
                     f"Either add them to the VCF or remove them from the {args.core_set} column.")
        self.core_samples = [s for s in self.all_samples if s in core]
        if len(self.core_samples) < 11:
            sys.exit(f"Only {len(self.core_samples)} {args.core_set} samples; at least 11 are needed for 10 PCs.")
        self.bad_ld = len(self.core_samples) < MIN_CORE_FOR_LD
        if self.bad_ld:
            print(f"WARNING: {len(self.core_samples)} {args.core_set} samples (< {MIN_CORE_FOR_LD}). LD pruning will run with plink2 --bad-ld; "
                  f"LD estimates from this few samples are noisier, so treat pruned (whole-genome / per-LG) PCAs with some caution.")
        no_meta = sorted(set(self.all_samples) - set(self.meta.index))
        if no_meta:
            print(f"WARNING: {len(no_meta)} VCF samples are not in the SampleDatabase tab (plotted grey, no metadata): {', '.join(no_meta)}")
        print(f"Samples: {len(self.all_samples)} in the PCA, {len(self.core_samples)} build the PCs ({args.core_set}), "
              f"{len(self.bionano & set(self.all_samples))} Bionano, {len(self.highlight & set(self.all_samples))} highlighted")

    def _map_linkage_groups(self):
        # Hard-coded LG -> accession map (LINKAGE_GROUPS); any LG with no variants in the VCF is skipped with a warning
        if args.genome not in LINKAGE_GROUPS:
            sys.exit(f"No linkage-group map for {args.genome}. Add one to LINKAGE_GROUPS at the top of this script "
                     f"(the accessions are in the @SQ lines of Intervals/{args.genome}/1.interval_list).")
        contigs_in_vcf = set(run(['tabix', '-l', self.in_vcf], capture = True).split())
        missing = [lg for lg, contig in LINKAGE_GROUPS[args.genome].items() if contig not in contigs_in_vcf]
        if missing:
            print(f"WARNING: no variants in the VCF for {', '.join(missing)}; skipped")
        self.lg_to_contig = {lg: contig for lg, contig in LINKAGE_GROUPS[args.genome].items() if contig in contigs_in_vcf}

    def _resolve_regions(self):
        # name -> (plot title, bcftools region or None for whole genome, LD prune?)
        regions = {}
        for request in args.regions:
            if request == 'All':
                for lg, contig in self.lg_to_contig.items():
                    regions[lg] = (f"{lg} ({contig})", contig, True)
            elif request == 'Whole':
                regions['Whole'] = ('Whole genome', None, True)
            elif request == 'Inversions':
                for name in INVERSIONS:
                    regions[name] = self._inversion_region(name)
            elif request in self.lg_to_contig:
                regions[request] = (f"{request} ({self.lg_to_contig[request]})", self.lg_to_contig[request], True)
            elif request in INVERSIONS:
                regions[request] = self._inversion_region(request)
            else:
                sys.exit(f"Unknown region '{request}'. Options: All, Whole, Inversions, {', '.join(self.lg_to_contig)}, {', '.join(INVERSIONS)}")
        self.regions = regions
        print('Regions: ' + ', '.join(regions))

    def _inversion_region(self, name):
        lg_number, start, end = INVERSIONS[name]
        contig = self.lg_to_contig.get(f'LG{lg_number}')
        if contig is None:
            sys.exit(f"{name} needs LG{lg_number}, which isn't in the VCF")
        return (f"{name.replace('_inv', '')} inversion ({contig}:{start:,}-{end:,})", f"{contig}:{start}-{end}", False)

    # ---------------------------------------------------------------- VCF preprocessing
    def preprocess(self):
        # Original _create_ecogroup_specific_vcf + _preprocess_vcfs steps. Skipped if the input VCF and sample lists are unchanged.
        stat = os.stat(self.in_vcf)
        run_info = f"input\t{self.in_vcf}\nsize\t{stat.st_size}\nmtime\t{int(stat.st_mtime)}\ncore\t{','.join(self.core_samples)}\nall\t{','.join(self.all_samples)}\n"
        info_file = self.out_dir + 'run_info.txt'
        outputs = [self.core_filtered_vcf, self.core_filtered_vcf + '.tbi', self.all_filtered_vcf, self.all_filtered_vcf + '.tbi']
        if not args.overwrite and os.path.exists(info_file) and open(info_file).read() == run_info and all(os.path.exists(f) for f in outputs):
            print('Preprocessed VCFs already exist for this input and these samples - reusing them (use --overwrite to rebuild).')
            return
        if os.path.exists(info_file):
            os.remove(info_file) # so an interrupted rebuild is never mistaken for a finished one
        with open(self.core_list, 'w') as fh:
            fh.write('\n'.join(self.core_samples) + '\n')
        with open(self.all_list, 'w') as fh:
            fh.write('\n'.join(self.all_samples) + '\n')
        t = str(args.threads)
        print('Writing an all-sample VCF (--min-ac=1) and a core-sample VCF ...')
        run_parallel([['bcftools', 'view', self.in_vcf, '--samples-file', self.all_list, '--min-ac=1', '--no-update', '--threads', t, '-o', self.all_vcf, '-O', 'z'],
                      ['bcftools', 'view', self.in_vcf, '--samples-file', self.core_list, '--threads', t, '-o', self.core_vcf, '-O', 'z']])
        print('Recalculating AF in the core samples ...')
        run(['bcftools', '+fill-tags', self.core_vcf, '-Oz', '-o', self.core_af_vcf, '--', '-t', 'AF'])
        print('Removing variants with core AF < 0.05 or AF > 0.95 ...')
        run(['bcftools', 'view', '--exclude', 'AF<0.05 || AF>0.95', self.core_af_vcf, '-o', self.core_filtered_vcf, '-Oz'])
        run(['bcftools', 'query', '-f', '%CHROM\t%POS\n', self.core_filtered_vcf, '-o', self.positions])
        print('Keeping only those positions in the all-sample VCF (vcftools --positions) ...')
        run(['vcftools', '--positions', self.positions, '--gzvcf', self.all_vcf, '--recode', '--out', self.all_filtered_prefix])
        run(['bgzip', '-f', '-@', t, self.all_filtered_prefix + '.recode.vcf'])
        run_parallel([['tabix', '-f', '-p', 'vcf', self.all_filtered_vcf], ['tabix', '-f', '-p', 'vcf', self.core_filtered_vcf]])
        with open(self.positions) as fh:
            print(f"{sum(1 for _ in fh):,} variants kept after the core-sample AF filter")
        with open(info_file, 'w') as fh:
            fh.write(run_info)

    # ---------------------------------------------------------------- plink2 per region
    def region_prefix(self, name):
        return self.out_dir + 'PCA/' + name + '/' + name

    def run_plink(self, name):
        title, region, ld_prune = self.regions[name]
        prefix = self.region_prefix(name)
        os.makedirs(os.path.dirname(prefix), exist_ok = True)
        print(f"\n===== {title}")
        if region is None: # Whole genome uses the genome-wide files, as before
            core_vcf, all_vcf = self.core_filtered_vcf, self.all_filtered_vcf
        else:
            core_vcf, all_vcf = prefix + '_core.vcf.gz', prefix + '_all.vcf.gz'
            run_parallel([['bcftools', 'view', '-r', region, self.all_filtered_vcf, '-o', all_vcf, '-Oz'],
                          ['bcftools', 'view', '-r', region, self.core_filtered_vcf, '-o', core_vcf, '-Oz']])
        # Steps 1-2: VCF -> pfiles; the all-sample pfile gets chr:pos IDs
        run(['plink2', '--vcf', all_vcf, '--out', prefix + '_all_step_1', '--allow-extra-chr'])
        run(['plink2', '--vcf', core_vcf, '--out', prefix + '_core_step_1', '--allow-extra-chr'])
        run(['plink2', '--pfile', prefix + '_all_step_1', '--set-missing-var-ids', '@:#', '--make-pgen', '--out', prefix + '_all_step_2', '--allow-extra-chr'])

        # Step 3: PCA on the core samples (LD-pruned for whole genome / linkage groups, not for inversions)
        pca_command = ['plink2', '--pfile', prefix + '_core_step_1', '--freq', 'counts', '--pca', 'allele-wts', '--out', prefix + '_pca',
                       '--allow-extra-chr', '--set-missing-var-ids', '@:#', '--max-alleles', '2']
        if ld_prune:
            run(['plink2', '--pfile', prefix + '_core_step_1', '--out', prefix + '_prune', '--allow-extra-chr', '--set-missing-var-ids', '@:#',
                 '--indep-pairwise', '50', '5', '0.1'] + (['--bad-ld'] if self.bad_ld else []))
            pca_command += ['--extract', prefix + '_prune.prune.in']
        run(pca_command)

        # Remove 0-count and fixed alleles from the .acount file (original rule)
        acounts_df = pd.read_csv(prefix + '_pca.acount', sep = '\t')
        no_extremes_df = acounts_df[(acounts_df['ALT_CTS'] != 0) & (acounts_df['ALT_CTS'] != max(acounts_df['OBS_CT']))]
        no_extremes_df.to_csv(prefix + '_pca.acount', sep = '\t', index = False)

        # Step 4: project every sample onto the core PCs
        run(['plink2', '--pfile', prefix + '_all_step_2', '--read-freq', prefix + '_pca.acount', '--score', prefix + '_pca.eigenvec.allele',
             '2', '5', 'header-read', 'no-mean-imputation', 'variance-standardize', '--score-col-nums', '6-15',
             '--out', prefix + '_projection', '--allow-extra-chr'])

    # ---------------------------------------------------------------- plots
    def plot(self, name):
        title, region, ld_prune = self.regions[name]
        prefix = self.region_prefix(name)
        if not os.path.exists(prefix + '_projection.sscore'):
            print(f"WARNING: no results for {name} yet (run without --plot_only) - skipping plot")
            return
        eigenval = pd.read_csv(prefix + '_pca.eigenval', header = None)[0]
        variance = (eigenval / eigenval.sum() * 100).round(2) # % of the variance captured by the computed PCs, as before
        scores = pd.read_csv(prefix + '_projection.sscore', sep = '\t')
        scores = scores.rename(columns = {'#IID': 'SampleID', 'IID': 'SampleID'})
        df = scores[['SampleID', 'PC1_AVG', 'PC2_AVG']].join(self.meta, on = 'SampleID')
        df['Core'] = df.SampleID.isin(self.core_samples)
        df['Bionano'] = df.SampleID.isin(self.bionano)
        df['Highlight'] = df.SampleID.isin(self.highlight) & ~df.Bionano
        df['Role'] = df.Core.map({True: args.core_set, False: 'projected'})
        for column in ['Ecogroup', 'Species', 'Subgroup', 'Sex', 'Inversion10', 'ProjectID']:
            df[column] = df[column].fillna('')
        # Drawing order: ordinary samples, then highlighted, then Bionano on top
        df['order'] = df.Highlight.astype(int) + 2 * df.Bionano.astype(int)
        df = df.sort_values('order', kind = 'stable')
        df.drop(columns = ['order']).to_csv(prefix + '_projection_with_metadata.csv', index = False)

        for kind in ['html', 'pdf']:
            figure = go.Figure(go.Scatter(
                x = df.PC1_AVG, y = df.PC2_AVG, mode = 'markers',
                marker = dict(color = [('black' if b else MALINSKY_COLORS.get(e, UNKNOWN_COLOR)) for b, e in zip(df.Bionano, df.Ecogroup)],
                              symbol = ['x' if b else 'circle' for b in df.Bionano],
                              size = [SIZE[kind]['bionano' if b else 'normal'] for b in df.Bionano],
                              opacity = [OPACITY['bionano' if b else 'normal'] for b in df.Bionano],
                              line = dict(width = [0 if (b or (c and not h)) else OUTLINE[kind]['highlight' if h else 'noncore']
                                                   for b, c, h in zip(df.Bionano, df.Core, df.Highlight)],
                                          color = [HIGHLIGHT_COLOR if h else 'black' for h in df.Highlight])),
                customdata = df[['SampleID', 'Ecogroup', 'Species', 'Subgroup', 'Sex', 'Inversion10', 'ProjectID', 'Role']].values,
                hovertemplate = ('SampleID=%{customdata[0]}<br>Ecogroup=%{customdata[1]}<br>Species=%{customdata[2]}<br>Subgroup=%{customdata[3]}'
                                 '<br>Sex=%{customdata[4]}<br>Inversion10=%{customdata[5]}<br>ProjectID=%{customdata[6]}<br>%{customdata[7]}<extra></extra>')))
            figure.update_layout(title = title, xaxis_title = f"PC1 {variance[0]}%", yaxis_title = f"PC2 {variance[1]}%", showlegend = False)
            if kind == 'html':
                figure.write_html(self.html_dir + name + '_PCA.html')
            else:
                figure.update_layout(width = 330, height = 330, margin = dict(l = 20, r = 20, t = 30, b = 20))
                try:
                    figure.write_image(self.pdf_dir + name + '_PCA.pdf')
                except Exception as error: # usually kaleido not installed
                    print(f"WARNING: PDF not written for {name} ({error}). Install kaleido: pip install kaleido", file = sys.stderr)
        print(f"Plots written for {name}")

    def run(self):
        if not args.plot_only:
            self.preprocess()
            for name in self.regions:
                self.run_plink(name)
        for name in self.regions:
            self.plot(name)


if __name__ == '__main__':
    state = start_run_logs()
    PCAMaker().run()
    state['done'] = True