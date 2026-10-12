import argparse, pdb, os, subprocess, datetime, time, multiprocessing, pathlib, random, re, sys, atexit, platform
import pandas as pd
from pyfaidx import Fasta
from helper_modules.nikesh_file_manager import FileManager as FM

# Shortcut names accepted by -e/--ecogroups, expanded to the Ecogroup values used in the SampleData sheet.
# There is no 'All' shortcut any more: leave -e off for "no ecogroup filter" (keeps samples with a blank Ecogroup, e.g.
# all PM30 samples), or use --whole_alignment_tab for "every sample aligned to this genome".
# 'Shallow_Benthic2' and 'Riverine' were removed: they only existed in the old Ecogroup_PTM column.
ECOGROUP_SHORTCUTS = {'Non_Riverine': ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic', 'Rhamphochromis', 'Diplotaxodon'],
                      'Lake_Malawi': ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic', 'Rhamphochromis', 'Diplotaxodon', 'AC'],
                      'Rock_Sand': ['Mbuna', 'Utaka', 'Shallow_Benthic', 'Deep_Benthic'],
                      'Sand': ['Utaka', 'Shallow_Benthic', 'Deep_Benthic']}
ECOGROUPS = ['Mbuna', 'Shallow_Benthic', 'Deep_Benthic', 'Utaka', 'Diplotaxodon', 'Rhamphochromis', 'AC', 'Hybrid']

parser = argparse.ArgumentParser(usage='This pipeline will take in a set of unaligned bam files generated from Illumina sequencing reads, and call variants using GATK HaplotypeCaller')
parser.add_argument('reference_genome', help = 'genome version = Alignment Database tab name (e.g. Mzebra_GT3_NCBI, Mconophoros_GT1_NCBI)', choices = ['Mzebra_GT3_NCBI', 'Mconophoros_GT1_NCBI', 'Mzebra_GT3', 'Mconophoros_GT1', 'Mzebra_UMD2a'])
parser.add_argument('-e', '--ecogroups', help = 'filter on Ecogroup in the SampleDatabase tab (combined with --project_ids/--species/--subgroups like PI setSamples). Shortcuts: ' + ', '.join(ECOGROUP_SHORTCUTS) + '. Leave off for no ecogroup filter.', choices = ECOGROUPS + list(ECOGROUP_SHORTCUTS), nargs = '*', default = None)
parser.add_argument('-i', '--import_databases', help = 'Call this flag to run GenomicsDBImport on the samples in the database', action = 'store_true')
parser.add_argument('-g', '--genotype', help = 'Call this flag to run GenotypeGVCFs on the samples in the database', action = 'store_true')
parser.add_argument('-d', '--download_GVCF_data', help = 'Use this flag if you need to download the GVCF files from Dropbox to include in the VCF analysis', action = 'store_true')
parser.add_argument('-r', '--regions', help = 'list of linkage groups for which analyses will run', nargs = '*', default = ['All'])
parser.add_argument('-b', '--download_bams', help = 'Download the BAM files from the cloud on which to call HaplotypeCaller', action = 'store_true')
parser.add_argument('-H', '--efficient_haplotypecaller', help = 'use this flag to download BAM files and run HaplotypeCaller on samples', action = 'store_true')
parser.add_argument('-m', '--memory', help = 'How much memory, in GB, to allocate to each child process', default = [4], nargs = 1)
parser.add_argument('-u', '--unmapped', help = 'Use this flag to run -i and -g on the unmapped contigs in the genome', action = 'store_true')
# --- Analysis naming / folders ---
# Each analysis gets its own folders: Outputs/<genome>/<analysis_id>/ and Databases/<genome>/<analysis_id>/.
# Starting a new analysis with a name that already exists creates <name>_2, <name>_3, ... instead of overwriting.
# To run later steps (-g, --concat_and_index, or a re-try) on an analysis that already exists, pass the same --analysis_id with --resume.
parser.add_argument('--analysis_id', '--analysis_name', help = 'name of this analysis = output folder name = final VCF name (ideally = the AnalysisSets column header). Default: the --analysis_set name if given, otherwise VC_<genome>_<YYYY-MM-DD>')
parser.add_argument('--resume', help = 'continue an existing analysis folder (exact --analysis_id) instead of creating a new incremented one. Reuses its frozen sample list.', action = 'store_true')
# --- Sample selection (replaces the old -s/--sampleIDs modes) ---
# Old -s mode     -> new equivalent
#   All            -> -e <ecogroups>  (or --whole_alignment_tab for everything aligned to the genome)
#   alignment_file -> --whole_alignment_tab
#   custom         -> --extra_samples S1 S2 ...   (no interactive prompt, so it works under nohup/screen)
#   v2_column      -> --analysis_set <checkbox column in the AnalysisSets tab>
# All selectors are combined (union), then --exclude_samples are removed, then anything not aligned to <genome> is dropped.
parser.add_argument('--analysis_set', help = 'checkbox column in the AnalysisSets tab; TRUE rows are included')
parser.add_argument('--project_ids', nargs = '*', help = 'select samples by ProjectID, e.g. NewAulonocara (== the PM30 subgroup)')
parser.add_argument('--species', nargs = '*', help = 'select samples by Species')
parser.add_argument('--subgroups', nargs = '*', help = 'select samples by Subgroup')
parser.add_argument('--whole_alignment_tab', action = 'store_true', help = 'select every sample in the Alignment Database tab for this genome')
parser.add_argument('--extra_samples', nargs = '*', help = 'additional SampleIDs to include')
parser.add_argument('--exclude_samples', nargs = '*', help = 'SampleIDs to remove from the selection')
parser.add_argument('--refresh_samples', action = 'store_true', help = 'with --resume: re-select samples from the sheets instead of reusing the frozen sample list')
parser.add_argument('--concat_and_index', help = 'Use this flag to concatenate the interval VCFs output by GenotypeGVCFs into <analysis_name>.vcf.gz in the analysis folder', action = 'store_true')
parser.add_argument('--Output', help = 'Use this flag to specify that the files for use in the pipleine or that need to be downloaded are located at /Output in Utaka', action = 'store_true')
parser.add_argument('--upload', help = 'Use this flag to upload GVCF files from the server to Dropbox', action =  'store_true')
parser.add_argument('--all_sites', help = 'use this flag if you want to gun GenotypeGVCFs in --include-non-varinat-sites mode', action = 'store_true')
parser.add_argument('--concurrent_processes', help = 'specify the number of processes to start concurrently', type = int, default = 4)
parser.add_argument('--local_test', help = 'when this flag is called, variables will be preset to test the code locally', action = 'store_true')
args = parser.parse_args()

# Fail early on argument combinations that can't work
_sample_selectors = [args.analysis_set, args.project_ids, args.species, args.subgroups, args.ecogroups, args.extra_samples]
if args.resume and not (args.analysis_id or args.analysis_set):
    parser.error('--resume needs the exact --analysis_id (or --analysis_set, if that was used as the name) of the analysis folder to continue')
if args.refresh_samples and not args.resume:
    parser.error('--refresh_samples only makes sense with --resume (a new analysis always selects samples fresh)')
if not args.resume and not args.local_test and not args.whole_alignment_tab and not any(_sample_selectors):
    parser.error('no samples selected: give at least one of --analysis_set, --project_ids, --species, --subgroups, -e/--ecogroups, --extra_samples, --whole_alignment_tab')
if not args.resume and args.genotype and not args.import_databases:
    parser.error('-g needs the GenomicsDBs from -i. Run them together, or continue an existing analysis with --resume --analysis_id <name>')
if not args.resume and args.concat_and_index and not args.genotype:
    parser.error('--concat_and_index needs the interval VCFs from -g. Run them together, or continue an existing analysis with --resume --analysis_id <name>')

# Example commands:
#   PM30 run:   python3 callVariants.py Mzebra_GT3_NCBI --analysis_set PM30 -i -g --concat_and_index
#               (analysis/folder name defaults to the column name, here PM30)
#   Cait's run: python3 callVariants.py Mconophoros_GT1_NCBI --analysis_set CaitVCF -i -g --concat_and_index
#               (same samples as --whole_alignment_tab today, but the CaitVCF column stays fixed if more samples get aligned later)
#   Continue Cait's run later, e.g. only concatenating: python3 callVariants.py Mconophoros_GT1_NCBI --analysis_id CaitVCF --resume --concat_and_index

"""
TODO:
1. Paralellize gatk Haplotypecaller to run much more efficiently
    Be sure to include code that will not overwrite files if they exist
2. CODE DOES NOT EXIST ANYMORE TO DOWNLOAD BAM OR GVCF FILES OR ITS BROKEN. NEED TO FIX THIS!!
    GVCF file downloader has been fixed but I do not know how to download it straight to /Output b/c the directory structure is not mirrored on Dropbox and /Output. i think i need to add an option in FIleManager that allows me to specify that the master dir is /Output instead of /Data

PRE-RUN CHECKLIST (Oct 2026 sample-selection update):
- [ ] Server env: pip install gspread gspread-dataframe google-auth
- [ ] If Google returns a 403 "insufficient scopes" error, see the scopes note in nikesh_file_manager._authenticateGS
- [ ] python3 make_intervals.py <genome> has been run ON THE SERVER you're using: it writes Intervals/<genome>/ and also downloads
      the fasta and makes the .fai/.dict that GATK needs there
"""

# Paths of this run's log files (set by start_run_logs) and whether the run reached the end
RUN_LOG = {'log': None, 'err': None, 'start': None, 'completed': False}

def start_run_logs():
    # Every run writes a new pair of files to logs/ in the folder that holds callVariants.py:
    #   logs/<YYYY-MM-DD_HH-MM-SS>_<genome>_<analysis>.log   everything printed to stdout
    #   logs/<YYYY-MM-DD_HH-MM-SS>_<genome>_<analysis>.err   everything printed to stderr (Python errors/tracebacks; GATK and
    #                                                         rclone also write their normal progress messages to stderr)
    # The redirection is done at the file-descriptor level, so it also captures the worker processes and the programs the
    # script calls (gatk, rclone, bcftools), not just Python's print().
    log_dir = os.path.dirname(os.path.abspath(__file__)) + '/logs/'
    os.makedirs(log_dir, exist_ok = True)
    start = datetime.datetime.now()
    label = re.sub(r'[^A-Za-z0-9._-]+', '_', args.analysis_id or args.analysis_set or 'no_analysis_name')
    base = log_dir + f"{start.strftime('%Y-%m-%d_%H-%M-%S')}_{args.reference_genome}_{label}"
    suffix = 1
    while os.path.exists(base + ('' if suffix == 1 else f"_{suffix}") + '.log'): # two runs started in the same second
        suffix += 1
    base += '' if suffix == 1 else f"_{suffix}"
    RUN_LOG.update(log = base + '.log', err = base + '.err', start = start)

    # Last thing shown in the terminal / screen window; everything after this goes to the files
    print(f"callVariants.py started {start:%Y-%m-%d %H:%M:%S}. Output is being written to:\n  stdout: {RUN_LOG['log']}\n  stderr: {RUN_LOG['err']}")
    print(f"Follow progress with: tail -f {RUN_LOG['log']}")
    sys.stdout.flush()
    sys.stderr.flush()
    log_fh, err_fh = open(RUN_LOG['log'], 'w'), open(RUN_LOG['err'], 'w')
    os.dup2(log_fh.fileno(), sys.stdout.fileno())
    os.dup2(err_fh.fileno(), sys.stderr.fileno())
    sys.stdout.reconfigure(line_buffering = True) # keep Python's lines in order with the gatk/rclone output in the same file
    sys.stderr.reconfigure(line_buffering = True)

    header = (f"==== callVariants.py run started {start:%Y-%m-%d %H:%M:%S}\n"
              f"command : python3 {' '.join(sys.argv)}\n"
              f"host    : {platform.node()} (user {os.getenv('USER')})\n"
              f"cwd     : {os.getcwd()}\n")
    print(header + f"stderr  : {RUN_LOG['err']}\n")
    print(header + f"stdout  : {RUN_LOG['log']}\n", file = sys.stderr)
    atexit.register(_finish_run_logs)

def _finish_run_logs():
    # Last line of both files: how the run ended, when, and how long it took
    end = datetime.datetime.now()
    status = 'PIPELINE RUN COMPLETE' if RUN_LOG['completed'] else 'STOPPED BEFORE COMPLETION (see the .err file for the error)'
    line = f"\n==== {status} | finished {end:%Y-%m-%d %H:%M:%S} | elapsed {str(end - RUN_LOG['start']).split('.')[0]}"
    print(line)
    print(line, file = sys.stderr)

class VariantCaller:
    def __init__(self, args):
        self.args = args
        self.genome = args.reference_genome
        self.fm_obj = FM(self.genome)
        self.memory = args.memory
        self.concurrent_processes = args.concurrent_processes
        self.current_time = datetime.datetime.now()

        # Expand ecogroup shortcuts (e.g. Rock_Sand) into the Ecogroup values used in the SampleData sheet.
        # None = no ecogroup filter.
        self.ecogroups = None
        if args.ecogroups:
            expanded = []
            for ecogroup in args.ecogroups:
                expanded += ECOGROUP_SHORTCUTS.get(ecogroup, [ecogroup])
            self.ecogroups = sorted(set(expanded))

        # Interval lists are genome specific (they contain contig names). make_intervals.py writes them to
        # Intervals/<genome>/1.interval_list ... N.interval_list in the folder that holds this script.
        self.intervalDir = os.path.dirname(os.path.abspath(__file__)) + '/Intervals/' + self.genome + '/'
        # TODO: the unmapped-contig interval lists are still the old cwd folder; they're genome specific too (contig names differ
        # between assemblies). make_intervals.py lists every unplaced/mito contig it left out in Intervals/<genome>/excluded_contigs.tsv.
        self.unmappedIntervalDir = os.getcwd() + '/unmapped_contig_intervals/'
        # Interval names used for -i / -g / --concat_and_index: as many N.interval_list files as make_intervals.py wrote (96 by default)
        interval_count = len([f for f in os.listdir(self.intervalDir) if re.fullmatch(r'\d+\.interval_list', f)]) if os.path.isdir(self.intervalDir) else 96
        self.intervals = [str(interval) for interval in range(1, interval_count + 1)]

        # Code block for determining which linkage groups will be processed by the script:
        self.linkage_group_map = {'LG1': 'NC_036780.1', 'LG2':'NC_036781.1', 'LG3':'NC_036782.1', 'LG4':'NC_036783.1', 'LG5':'NC_036784.1', 'LG6':'NC_036785.1',
                             'LG7':'NC_036786.1', 'LG8':'NC_036787.1', 'LG9':'NC_036788.1', 'LG10':'NC_036789.1', 'LG11':'NC_036790.1', 'LG12':'NC_036791.1',
                             'LG13':'NC_036792.1', 'LG14':'NC_036793.1', 'LG15':'NC_036794.1', 'LG16':'NC_036795.1', 'LG17':'NC_036796.1', 'LG18':'NC_036797.1',
                             'LG19':'NC_036798.1', 'LG20':'NC_036799.1', 'LG21':'NC_036800.1', 'LG22':'NC_036801.1', 'mito': 'NC_027944.1'}
        self.linkage_groups = args.regions
        if self.linkage_groups == ['All']:
            self.linkage_groups = list(self.linkage_group_map.values())
        else:
            regions_list = []
            for region in self.linkage_groups:
                if region in self.linkage_group_map.keys():
                    regions_list.append(self.linkage_group_map[region])
                else:
                    raise Exception(region + ' is not a valid option')
            self.linkage_groups = regions_list
            duplicate_set_test = set(self.linkage_groups)
            if len(duplicate_set_test) != len(self.linkage_groups):
                raise Exception('A repeat region has been provided')

        # pre-setting variables for local testing. Pass in the first 3 LGs only since the interval file has been created for only these.
        # (the local test sample list is set in _select_samples)
        if args.local_test:
            self.memory = [1]
            self.linkage_groups = ['NC_036780.1', 'NC_036781.1', 'NC_036782.1']
            self.concurrent_processes = 10

        if not args.local_test and (args.import_databases or args.genotype):
            self._check_interval_files()

        # Per-analysis folders first (the frozen sample list lives inside the analysis folder), then sample selection
        self._setup_analysis_dirs()
        self._select_samples()
        # pdb.set_trace()

    def _setup_analysis_dirs(self):
        # Every analysis gets its own folders so runs never overwrite each other. The folder layout itself is defined in
        # FileManager.createAnalysisIDFiles (so other scripts can find an analysis's files by its ID); this method decides the name:
        # A new analysis never reuses an existing name: if <name> is taken, <name>_2, <name>_3, ... is used.
        # --resume continues the existing analysis with exactly that name.
        if self.args.analysis_id:
            base_name = self.args.analysis_id
        elif self.args.analysis_set:
            base_name = self.args.analysis_set
        else:
            base_name = 'VC_' + self.genome + '_' + self.current_time.strftime('%Y-%m-%d')
        base_name = re.sub(r'[^A-Za-z0-9._-]+', '_', base_name.strip()) # make it folder-safe (sheet column headers may contain spaces)

        if self.args.resume:
            if not self.fm_obj.analysisIDExists(base_name):
                raise Exception(f"--resume: no analysis named {base_name} found for {self.genome}. Existing analyses: {self.fm_obj.returnAnalysisIDs()}")
            analysis_name = base_name
            print(f"Resuming existing analysis {analysis_name}")
        else:
            analysis_name = base_name
            increment = 2
            while self.fm_obj.analysisIDExists(analysis_name):
                analysis_name = f"{base_name}_{increment}"
                increment += 1
            if analysis_name != base_name:
                print(f"An analysis named {base_name} already exists for {self.genome}. Starting a new analysis named {analysis_name} instead.")
                print(f"(To continue the existing one, re-run with --resume --analysis_id {base_name})")

        # Paths come from the file manager (Outputs/<genome>/<name>/..., Databases/<genome>/<name>/); short aliases are kept here
        self.fm_obj.createAnalysisIDFiles(analysis_name)
        self.analysisName = analysis_name
        self.analysisOutputDir = self.fm_obj.localAnalysisOutputDir
        self.intervalVCFDir = self.fm_obj.localAnalysisIntervalVCFDir
        self.analysisConcatDir = self.fm_obj.localAnalysisConcatDir
        self.analysisDatabasesDir = self.fm_obj.localAnalysisDatabasesDir
        self.sampleManifestFile = self.fm_obj.localSampleManifestFile # frozen sample list for this analysis
        self.sampleMapFile = self.fm_obj.localSampleMapFile
        self.compressionListFile = self.fm_obj.localCompressionListFile
        for directory in [self.intervalVCFDir, self.analysisConcatDir, self.analysisDatabasesDir]:
            os.makedirs(directory, exist_ok = True)

        # Record every command run against this analysis folder
        with open(self.fm_obj.localCommandsLogFile, 'a') as fh:
            fh.write(self.current_time.strftime('%Y-%m-%d %H:%M:%S') + '\t' + ' '.join(sys.argv)
                     + (f"\tlog: {RUN_LOG['log']}" if RUN_LOG['log'] else '') + '\n')

        print(f"Analysis name: {self.analysisName}")
        print(f"  Outputs:   {self.analysisOutputDir}")
        print(f"  Databases: {self.analysisDatabasesDir}")

    def _select_samples(self):
        # Decide which samples go into this analysis and freeze that list in the analysis folder.
        # The Google Sheets are live and editable, so every later step (and every --resume) uses the frozen list, not the sheet.
        args = self.args
        previous_samples = None
        if os.path.exists(self.sampleManifestFile):
            with open(self.sampleManifestFile) as fh:
                previous_samples = [x.strip() for x in fh if x.strip()]

        if args.local_test:
            self.sampleIDs = ['CJ_2204_m', 'CV-006-m', 'LA_3006_m', 'MC-008-m', 'OC-001-m', 'small_file1', 'small_file2', 'small_file3', 'small_file4']
            source = 'local test preset'
        elif args.resume and previous_samples is not None and not args.refresh_samples:
            self.sampleIDs = previous_samples
            source = 'frozen sample list ' + self.sampleManifestFile
        else:
            if args.resume and not args.whole_alignment_tab and not any([args.analysis_set, args.project_ids, args.species, args.subgroups, self.ecogroups, args.extra_samples]):
                raise Exception(f"No frozen sample list found at {self.sampleManifestFile} and no sample selectors given. Add the selectors used originally (e.g. --analysis_set).")
            self.sampleIDs = self.fm_obj.selectSamples(analysis_set = args.analysis_set, projectIDs = args.project_ids, species = args.species,
                                                       ecogroups = self.ecogroups, subgroups = args.subgroups, extra_samples = args.extra_samples,
                                                       exclude_samples = args.exclude_samples, whole_alignment_tab = args.whole_alignment_tab)
            source = 'Google Sheets'

        self.sampleIDs = sorted(set(self.sampleIDs))
        if len(self.sampleIDs) == 0:
            raise Exception('No samples selected for this analysis. Check the selectors and that the samples are aligned to ' + self.genome)

        if previous_samples is not None and set(previous_samples) != set(self.sampleIDs):
            # Only reachable with --resume --refresh_samples: keep the old list and say exactly what changed
            backup = self.sampleManifestFile.replace('_samples.txt', '_samples.previous_' + self.current_time.strftime('%Y%m%d_%H%M%S') + '.txt')
            os.rename(self.sampleManifestFile, backup)
            added = sorted(set(self.sampleIDs) - set(previous_samples))
            removed = sorted(set(previous_samples) - set(self.sampleIDs))
            print(f"WARNING: sample list changed from the frozen list (old list saved to {backup}).")
            print(f"  Added ({len(added)}): {','.join(added)}")
            print(f"  Removed ({len(removed)}): {','.join(removed)}")
            print('  GenomicsDBs already built for this analysis still contain the old samples - re-run -i in a new analysis if that matters.')

        with open(self.sampleManifestFile, 'w') as fh:
            fh.write('\n'.join(self.sampleIDs) + '\n')

        print(f"Number of samples for this pipeline run is {len(self.sampleIDs)} (from {source}):")
        print(','.join(self.sampleIDs))

    def _check_interval_files(self):
        # Fail at startup (instead of hours in) if the interval lists for -i / -g are missing
        if self.args.unmapped:
            if not os.path.isdir(self.unmappedIntervalDir):
                raise FileNotFoundError('Unmapped contig interval directory not found: ' + self.unmappedIntervalDir)
            return
        missing = [interval for interval in self.intervals if not os.path.exists(self.intervalDir + interval + '.interval_list')]
        if missing:
            raise FileNotFoundError(f"{len(missing)} of {len(self.intervals)} interval lists missing from {self.intervalDir} (e.g. {', '.join(missing[:5])}). "
                                    f"Generate them with: python3 make_intervals.py {self.genome}")

    def _generate_sample_map(self):
        # Write the GenomicsDBImport sample map into the analysis folder. Stop before GenomicsDBImport starts if any sample is
        # missing its GVCF or the GVCF's .tbi index, so GenomicsDBImport doesn't die partway through the intervals.
        missing = []
        with open(self.sampleMapFile, 'w') as fh:
            for sampleID in self.sampleIDs:
                self.fm_obj.createSampleFiles(sampleID)
                if args.local_test:
                    gvcf_file, gvcf_index = self.fm_obj.localTestGVCFFile, self.fm_obj.localTestGVCFIndex
                elif args.Output:
                    gvcf_file, gvcf_index = self.fm_obj.StorageGVCFFile, self.fm_obj.StorageGVCFIndex
                else:
                    gvcf_file, gvcf_index = self.fm_obj.localGVCFFile, self.fm_obj.localGVCFIndex
                if not os.path.exists(gvcf_file):
                    missing.append(sampleID + '\tGVCF missing\t' + gvcf_file)
                elif not os.path.exists(gvcf_index):
                    missing.append(sampleID + '\t.tbi index missing (GVCF may be partial)\t' + gvcf_index)
                fh.write(sampleID + '\t' + gvcf_file + '\n')
        if missing:
            missing_file = self.fm_obj.localMissingGVCFsFile
            with open(missing_file, 'w') as fh:
                fh.write('\n'.join(missing) + '\n')
            raise FileNotFoundError(f"{len(missing)} of {len(self.sampleIDs)} samples are missing a GVCF or its .tbi index (list written to {missing_file}). "
                                    f"Download them with -d or generate them with -H before running -i.")

    def GVCF_downloader(self, sampleID):
        # Download the GVCF and GVCF.idx files for each sample in self.sampleIDs
        self.fm_obj.createSampleFiles(sampleID)
        if args.local_test: # for downloading files locally and testing if the script works 
            if pathlib.Path(self.fm_obj.localTestGVCFFile).exists():
                print('GCVF file for ' + sampleID + ' exists. Skipping sample...')
            else:
                print(f"Starting GVCF and GVCF index download for sample {sampleID} at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
                self.fm_obj.downloadData(self.fm_obj.localTestGVCFFile)
                self.fm_obj.downloadData(self.fm_obj.localTestGVCFIndex)
                print(f"Download of GVCF and GVCF index for sample {sampleID} complete at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
        elif args.Output: # NOTE: DOES NOT WORK. THE DROPBOX DTRUCTURE MUST BE MIRRORED. SO IF I WANT TO DOWNLOAD INTO /OUTPUTS, IT MUST EXIST IN /OUTPUTS IN DROPBOX. 
            if pathlib.Path(self.fm_obj.StorageGVCFFile).exists():
                print(print('GCVF file for ' + sampleID + ' exists at /Output. Skipping sample...'))
            else:
                print(f"Starting GVCF and GVCF index download for sample {sampleID} at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}. Sample is being stored at /Outputs")
                self.fm_obj.downloadData(self.fm_obj.StorageGVCFFile)
                self.fm_obj.downloadData(self.fm_obj.StorageGVCFIndex)
                print(f"Download of GVCF and GVCF index for sample {sampleID} complete at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
        else:
            # Always hand both files to rclone: it skips files that are already complete locally and re-copies partial ones
            # (an interrupted download used to leave a partial GVCF that was then skipped forever because it "existed").
            print(f"Starting GVCF and GVCF index download for sample {sampleID} at {datetime.datetime.now()} (rclone skips files that are already up to date)")
            self.fm_obj.downloadData(self.fm_obj.localGVCFFile)
            self.fm_obj.downloadData(self.fm_obj.localGVCFIndex)
            print(f"Download of GVCF and GVCF index for sample {sampleID} complete at {datetime.datetime.now()}")

    def uploadGVCFs(self, sampleID):
        # upload GVCF file and index for sampleID in self.sampleIDs if it doesn't already exist on Dropbox. 
        self.fm_obj.createSampleFiles(sampleID)
        if args.local_test: # code for testing using small local files. 
            # check to see if the file has already been uploaded to Dropbox. If it does. Skip it
            subprocess.check_output(['rclone', 'lsf', self.fm_obj.rcloneRemote + self.fm_obj.localBamFile], encoding='utf-8')
            print(f"Checking if {sampleID}'s GVCF file and Index exists on Dropbox...")
            cloud_gvcf_path = self.fm_obj.localTestGVCFFile.replace(self.fm_obj.localMasterDir, self.fm_obj.cloudMasterDir)
            cloud_gvcf_index_path = self.fm_obj.localTestGVCFIndex.replace(self.fm_obj.localMasterDir, self.fm_obj.cloudMasterDir)

            gvcf_code = subprocess.run(['rclone', 'lsf', cloud_gvcf_path]).returncode
            if gvcf_code == 0:
                print(f"Rclone found the GVCF file for {sampleID} on Dropbox. Skipping")
            else:
                print(f"GVCF file for {sampleID} is not on dropbox. Starting upload for GVCF file at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
                self.fm_obj.uploadData(self.fm_obj.localTestGVCFFile)
                print(f"GVCF file upload for {sampleID} complete at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
                
            index_code = subprocess.run(['rclone', 'lsf', cloud_gvcf_index_path]).returncode
            if index_code == 0:
                print(f"Rclone found the GVCF Index for {sampleID} on Dropbox. Skipping...")
            else:
                print(f"Index file for {sampleID} is not on dropbox. Starting upload for GVCF Index at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
                self.fm_obj.uploadData(self.fm_obj.localTestGVCFIndex)
                print(f"GVCF Inxex upload for {sampleID} complete at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
        else: # code for uploading server files
            print(f"Checking if {sampleID}'s GVCF file and Index exists on Dropbox...")
            cloud_gvcf_path = self.fm_obj.localGVCFFile.replace(self.fm_obj.localMasterDir, self.fm_obj.cloudMasterDir)
            cloud_gvcf_index_path = self.fm_obj.localGVCFIndex.replace(self.fm_obj.localMasterDir, self.fm_obj.cloudMasterDir)

            gvcf_code = subprocess.run(['rclone', 'lsf', cloud_gvcf_path]).returncode
            if gvcf_code == 0:
                print(f"Rclone found the GVCF file for {sampleID} on Dropbox. Skipping")
            else:
                print(f"GVCF file for {sampleID} is not on dropbox. Starting upload for GVCF file at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
                self.fm_obj.uploadData(self.fm_obj.localGVCFFile)
                print(f"GVCF file upload for {sampleID} complete at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
                
            index_code = subprocess.run(['rclone', 'lsf', cloud_gvcf_index_path]).returncode
            if index_code == 0:
                print(f"Rclone found the GVCF Index for {sampleID} on Dropbox. Skipping...")
            else:
                print(f"Index file for {sampleID} is not on dropbox. Starting upload for GVCF Index at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
                self.fm_obj.uploadData(self.fm_obj.localGVCFIndex)
                print(f"GVCF Inxex upload for {sampleID} complete at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")

    def EfficientHaplotypeCaller(self, sample):
        """
        NOTE: IMPLEMENT A FILE CHECKER TO MAKE SURE EXISTING FILES DO NOT GET DELETED. INCLUDE CODE FOR UPLOADING GVCF FILES AND INDICIES AS WELL 
        """
        print(f"Processing for sample {sample} started at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
        # to use multiprocess, this function will sequenctially need to perform BAM/BAI download, GVCF generation, and then deleting the BAM/BAI files
        if args.local_test:
            self.fm_obj.createSampleFiles(sample)
            if pathlib.Path(self.fm_obj.localTestBamFile).exists():
                print(f"BAM file for {sample} exists. Skipping download for this sample...")
            else:
                print(f"Downloading BAM and BAM index for {sample}...")
                print(f"BAM Download for sample {sample} started at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
                self.fm_obj.downloadData(self.fm_obj.localTestBamFile)
                self.fm_obj.downloadData(self.fm_obj.localTestBamIndex)
                print(f"BAM Download for sample {sample} ended at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}...")

            if pathlib.Path(self.fm_obj.localTestGVCFFile).exists():
                print('GCVF file for ' + sample + ' exists. Skipping sample... ')
            else:
                print('Generating GCVF file for ' + sample + '...')
                subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'HaplotypeCaller', '--emit-ref-confidence', 'GVCF', '-R', self.fm_obj.localGenomeFile, '-I', self.fm_obj.localTestBamFile, '-O', self.fm_obj.localTestGVCFFile])

            # commenting out bam removal code while GCVF code isn't yet working
            # print('Removing BAM file for ' + saple + '...')
            # os.remove(self.fm_obj.localTestBamFile)
            # os.remove(self.fm_obj.localTestBamIndex)
        else:
            self.fm_obj.createSampleFiles(sample) # create the file structure for the sample
            if pathlib.Path(self.fm_obj.localGVCFIndex).exists(): # If the index for that file exists, then the GVCF file is fine and Haplotypecaller does not need to be rerun 
                print(f"GVCF and Index created for sample {sample}. Skipping HaplotypeCaller run for this sample")
            else: # if index doesn't exist, process the sample
                print(f"No GVCF Index exists for sample {sample}. Starting processing")
                if pathlib.Path(self.fm_obj.localGVCFFile).exists(): # if a GVCF file is found for the sample without the index, that means it will be a partial one. Thus, it needs to be removed. 
                    print(f"Partial GVCF with no Index Detected. Deleting the GVCF file for sample {sample} and regenerating it.")
                    os.remove(self.fm_obj.localGVCFFile)
                new_gcvf_created = False
                if pathlib.Path(self.fm_obj.localBamFile).exists(): # If the BAM file exists for the sample, do not download it.
                    print('BAM file for ' + sample + ' exists. Skipping download for this sample... ')
                else: # download BAM since it does not exist on the server.
                    print('Downloading BAM and BAM index for ' + sample + '...')
                    self.fm_obj.downloadData(self.fm_obj.localBamFile)
                    self.fm_obj.downloadData(self.fm_obj.localBamIndex)
                
                # if the GVCFindex did not exist and a partial file was found, at this stage, the partial file would have been deleted and a new GVCF would need to be generated with the below code:
                print(f"Generating GCVF file for sample {sample}...")
                subprocess.run(['gatk', '--java-options', '-Xmx' +  str(self.memory[0]) + 'G', 'HaplotypeCaller', '--emit-ref-confidence', 'GVCF', '-R', self.fm_obj.localGenomeFile, '-I', self.fm_obj.localBamFile, '-O', self.fm_obj.localGVCFFile])
                new_gcvf_created = True

                if new_gcvf_created:
                    print('Removing BAM file for ' + sample + '...')
                    os.remove(self.fm_obj.localBamFile)
                    os.remove(self.fm_obj.localBamIndex)

        print(f"Sample {sample} finished processing at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")

    def RunGenomicsDBImport(self, interval):
        print(f"Processing for interval {interval} started at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
        if args.local_test: # for testing code locally. Keep or remove --overwrite-existing-genomicsdb-workspace as needed. Default is that it's true
            subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.analysisDatabasesDir + interval + '_database', '--intervals', self.intervalDir + interval + '.interval_list', '--sample-name-map', self.sampleMapFile,'--overwrite-existing-genomicsdb-workspace'])
        elif not args.unmapped: # for running on normal 96 intervals on the servers. Note that --overwrite-existing-genomicsdb-workspace is false, so nothing accidentally gets overwritten on the server. 
            subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.analysisDatabasesDir + interval + '_database', '--intervals', self.intervalDir + interval + '.interval_list', '--sample-name-map', self.sampleMapFile])
          # subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.fm_obj.localDatabasesDir + lg + '_database', '--intervals', os.getcwd() + '/unmapped_contig_intervals/' + lg + '.interval_list', '--sample-name-map', self.sampleMapFile, '--max-num-intervals-to-import-in-parallel', '4']) # this is code used if lgs are passed instaed of intervals. This looks to the all_lg_intervals dir and runs based on the old parallelization methods.
        else: # If running on unmapped contigs, use intervals described in the unmapped intervals dir within analysis_pipeline/
            subprocess.run(['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G', 'GenomicsDBImport', '--genomicsdb-workspace-path', self.analysisDatabasesDir + interval + '_database', '--intervals', self.unmappedIntervalDir + interval + '.interval_list', '--sample-name-map', self.sampleMapFile])
          
        print(f"Interval {interval} finished processing at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")

    def _gendbPath(self, interval):
        # GenotypeGVCFs input for one interval's GenomicsDB workspace. Written as a path relative to the directory the script
        # is launched from (the same idea as the old hard-coded 'gendb://../../../../../../' prefix), but computed, so it's
        # correct from any launch directory on any server.
        return 'gendb://' + os.path.relpath(self.analysisDatabasesDir + interval + '_database')

    def RunGenotypeGVCFs(self, interval):
        # You can ignore this function for now, but I implemented the parallelization code into it
        # check if a vcf index file already exists for the interval you're working on  if it does, skip it and move to another interval
        if pathlib.Path(self.intervalVCFDir + interval + '_output.vcf.gz.tbi').exists():
            print(f"VCF file's index for interval {interval} found. Skipping GenotypeGVCFs for this interval")
        else:
            print(f"Processing for interval {interval} started at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
            if args.local_test: # section is for local_testing code
                if args.all_sites: # This section is for running locally with --include-non-variant-sites true. IMPORTANT: This option is true only in gatk v4.3.0.0 - 24.07.23 NK. This flag is --all-sites on older versions!! Be careful with this option depending on what gatk version you're using!!
                    local_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', self._gendbPath(interval), '-O', self.intervalVCFDir + interval + '_output.vcf.gz', '--heterozygosity', '0.00175', '--include-non-variant-sites', 'true', '--intervals', self.intervalDir + interval + '.interval_list'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
                    local_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
                    local_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
                    local_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
                    local_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
                    subprocess.run(local_command)
                else: # this section is for running locally in normal mode, i.e. WITHOUT --include-non-variant-sites
                    local_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', self._gendbPath(interval), '-O', self.intervalVCFDir + interval + '_output.vcf.gz', '--heterozygosity', '0.00175'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
                    local_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
                    local_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
                    local_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
                    local_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
                    subprocess.run(local_command)
            elif not args.unmapped: # For running the standard 96 intervals (not unmappped and not local_test). gendb path is to the gendb located on the Utaka server. Will need to change if running on Mzebra 
                if args.all_sites: # This section is for running the 96 intervals with --include-non-variant-sites true. IMPORTANT: This option is true only in gatk v4.3.0.0 - 24.07.23 NK. This flag is --all-sites on older versions!! Be careful with this option depending on what gatk version you're using!!
                    genotypegvcfs_unmapped_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', self._gendbPath(interval), '-O', self.intervalVCFDir + interval + '_output.vcf.gz', '--heterozygosity', '0.00175', '--include-non-variant-sites', 'true', '--intervals', self.intervalDir + interval + '.interval_list'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
                    genotypegvcfs_unmapped_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
                    genotypegvcfs_unmapped_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
                    genotypegvcfs_unmapped_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
                    genotypegvcfs_unmapped_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
                    subprocess.run(genotypegvcfs_unmapped_command)
                else: # this section is for running the 96 intervals in normal mode, i.e. WITHOUT --include-non-variant-sites
                    genotypegvcfs_unmapped_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', self._gendbPath(interval), '-O', self.intervalVCFDir + interval + '_output.vcf.gz', '--heterozygosity', '0.00175'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
                    genotypegvcfs_unmapped_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
                    genotypegvcfs_unmapped_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
                    genotypegvcfs_unmapped_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
                    genotypegvcfs_unmapped_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
                    subprocess.run(genotypegvcfs_unmapped_command)
            else: # This section is for running the unmapped contigs.
                if args.all_sites: # This section is for running unmapped contigs with --include-non-variant-sites true. IMPORTANT: This option is true only in gatk v4.3.0.0 - 24.07.23 NK. This flag is --all-sites on older versions!! Be careful with this option depending on what gatk version you're using!!
                    genotypegvcfs_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', self._gendbPath(interval), '-O', self.intervalVCFDir + interval + '_output.vcf.gz', '--heterozygosity', '0.00175', '--include-non-variant-sites', 'true', '--intervals', self.unmappedIntervalDir + interval + '.interval_list'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
                    genotypegvcfs_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
                    genotypegvcfs_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
                    genotypegvcfs_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
                    genotypegvcfs_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
                    subprocess.run(genotypegvcfs_command)
                else: # this section is for running unmapped contigs in normal mode, i.e. WITHOUT --include-non-variant-sites
                    genotypegvcfs_command = ['gatk', '--java-options', '-Xmx' + str(self.memory[0]) + 'G','GenotypeGVCFs', '-R', self.fm_obj.localGenomeFile, '-V', self._gendbPath(interval), '-O', self.intervalVCFDir + interval + '_output.vcf.gz', '--heterozygosity', '0.00175'] # seq divergence estimated to be 0.01 - 0.25% in the Malinksy paper so I've set it at 0.00175 as the average of these values 
                    genotypegvcfs_command += ['-A', 'DepthPerAlleleBySample', '-A', 'Coverage', '-A', 'GenotypeSummaries', '-A', 'TandemRepeat', '-A', 'StrandBiasBySample']
                    genotypegvcfs_command += ['-A', 'ReadPosRankSumTest', '-A', 'AS_ReadPosRankSumTest', '-A', 'AS_QualByDepth', '-A', 'AS_StrandOddsRatio', '-A', 'AS_MappingQualityRankSumTest']
                    genotypegvcfs_command += ['-A', 'FisherStrand',  '-A', 'QualByDepth', '-A', 'RMSMappingQuality', '-A', 'DepthPerSampleHC']
                    genotypegvcfs_command += ['-G', 'StandardAnnotation', '-G', 'AS_StandardAnnotation', '-G', 'StandardHCAnnotation']
                    subprocess.run(genotypegvcfs_command)
            print(f"Interval {interval} finished processing at {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")

    def _merge_vcfs(self):
        # Concatenate this analysis's per-interval (or, with --unmapped, per-contig) VCFs into
        # Outputs/<genome>/<analysis_name>/vcf_concat_output/<analysis_name>.vcf.gz (+ .tbi).
        # An existing final VCF is never overwritten: <analysis_name>_2.vcf.gz, _3, ... is written instead.
        if args.unmapped:
            self.fasta = Fasta(self.fm_obj.localGenomeFile)
            self.unmapped_contigs = [contig for contig in self.fasta.keys() if contig.startswith('NW') or contig == 'NC_027944.1' or contig == 'ptg000146l_obj_unaligned']
            pieces = self.unmapped_contigs
            final_base = self.analysisName + '_unmapped'
        else:
            pieces = self.intervals
            final_base = self.analysisName

        vcf_files = [self.intervalVCFDir + piece + '_output.vcf.gz' for piece in pieces]
        missing = [vcf for vcf in vcf_files if not os.path.exists(vcf) or not os.path.exists(vcf + '.tbi')]
        if missing:
            raise FileNotFoundError(f"{len(missing)} of {len(vcf_files)} VCFs from GenotypeGVCFs (or their .tbi) are missing in {self.intervalVCFDir}, e.g. {missing[:3]}. "
                                    f"Re-run -g with --resume --analysis_id {self.analysisName} before --concat_and_index.")
        with open(self.compressionListFile, 'w') as fh:
            fh.write('\n'.join(vcf_files) + '\n')

        final_vcf = self.analysisConcatDir + final_base + '.vcf.gz'
        increment = 2
        while os.path.exists(final_vcf):
            final_vcf = self.analysisConcatDir + final_base + '_' + str(increment) + '.vcf.gz'
            increment += 1

        print(f"Starting concatenation of {len(vcf_files)} VCF files into {final_vcf} at {datetime.datetime.now()}")
        output = subprocess.run(['bcftools', 'concat', '-f', self.compressionListFile, '-Wtbi', '-O', 'z', '-o', final_vcf, '--threads', str(self.concurrent_processes)])
        if output.returncode != 0:
            raise Exception('bcftools concat failed while writing ' + final_vcf)
        print(f"{final_vcf} finished concatenating and indexing at {datetime.datetime.now()}")

        # Check the final VCF contains exactly the samples in this analysis's frozen sample list
        vcf_samples = set(subprocess.run(['bcftools', 'query', '-l', final_vcf], capture_output = True, encoding = 'utf-8').stdout.split())
        if vcf_samples == set(self.sampleIDs):
            print(f"Sample check passed: {final_vcf} contains exactly the {len(vcf_samples)} samples in {self.sampleManifestFile}")
        else:
            print(f"WARNING: samples in {final_vcf} don't match {self.sampleManifestFile}")
            print(f"  In VCF but not in sample list ({len(vcf_samples - set(self.sampleIDs))}): {','.join(sorted(vcf_samples - set(self.sampleIDs)))}")
            print(f"  In sample list but not in VCF ({len(set(self.sampleIDs) - vcf_samples)}): {','.join(sorted(set(self.sampleIDs) - vcf_samples))}")

    def multiprocess(self, function, sample_type):
        # TODO: find a way to load processes that would take the lonegst time to go first
        # Author: Lauren Sabo; edits made by NK
        if sample_type == 'lg':
            inputs = self.linkage_groups
        elif sample_type == 'sampleID':
            inputs = self.sampleIDs
        elif sample_type == 'interval':
            inputs = self.intervals
        elif sample_type == 'unmapped':
            # define the unmapped contigs
            self.fasta = Fasta(self.fm_obj.localGenomeFile)
            self.unmapped_contigs = [contig for contig in self.fasta.keys() if contig.startswith('NW') or contig == 'NC_027944.1' or contig == 'ptg000146l_obj_unaligned']
            inputs = self.unmapped_contigs
        if args.all_sites and not args.local_test: # if running GenotypeGVCFs in all-sites mode, 96 concurrent processes will deplete all available memory, so change the max concurrent processes to be 48 instead. Makes sure this only applies on server.
            print('All sites run detected. Limiting the number of concurrent processes to 48 to ensure each process has enough RAM to complete.')
            self.concurrent_processes = 48
        elif args.all_sites and args.local_test:
            self.concurrent_processes = 5
        concurrent_processes = min(self.concurrent_processes, len(inputs))
        try:
            with multiprocessing.Pool(processes=concurrent_processes) as pool:
                pool.map(function, inputs)
        except Exception as e:
            print(f"Error occurred during multiprocessing: {e}")

    def run_methods(self):
        if args.download_GVCF_data:
            self.multiprocess(self.GVCF_downloader, 'sampleID')
        if args.upload:
            self.multiprocess(self.uploadGVCFs, 'sampleID')
        if args.efficient_haplotypecaller:
            self.multiprocess(self.EfficientHaplotypeCaller, 'sampleID')
        if args.import_databases:
            # The sample map is only used by GenomicsDBImport. It's built here (after -d / -H) so GVCFs created in this same
            # run are found; it raises if any sample's GVCF or .tbi is missing.
            self._generate_sample_map()
        if args.import_databases and not args.unmapped:
            self.multiprocess(self.RunGenomicsDBImport, 'interval')
        if args.import_databases and args.unmapped:
            self.multiprocess(self.RunGenomicsDBImport, 'unmapped')
        if args.genotype and not args.unmapped:
            self.multiprocess(self.RunGenotypeGVCFs, 'interval')
        if args.genotype and args.unmapped:
            self.multiprocess(self.RunGenotypeGVCFs, 'unmapped')
        if args.concat_and_index:
            self._merge_vcfs()

if __name__ == "__main__":
    start_run_logs() # from here on, stdout/stderr go to logs/<date_time>_<genome>_<analysis>.log / .err
    variant_caller_obj = VariantCaller(args)
    variant_caller_obj.run_methods()
    RUN_LOG['completed'] = True