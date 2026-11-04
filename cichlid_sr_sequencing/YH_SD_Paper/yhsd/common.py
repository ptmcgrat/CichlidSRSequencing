"""Shared helpers: config, project paths (built from FileManager), sample tables, parallel runner."""
import os, sys
import pandas as pd
import yaml
from concurrent.futures import ProcessPoolExecutor, as_completed
from .common_light import sh, PKG_ROOT, config_path   # re-exported


def load_config(path=None):
    with open(path or config_path('paper.yaml')) as f:
        return yaml.safe_load(f)


def parallel(func, items, workers, desc=''):
    """Run func(item) over items in forked processes; returns ({item_index: result}, {item_index: error})."""
    out, errs = {}, {}
    with ProcessPoolExecutor(max_workers=max(1, int(workers))) as ex:
        futs = {ex.submit(func, it): i for i, it in enumerate(items)}
        for n, f in enumerate(as_completed(futs), 1):
            i = futs[f]
            try:
                out[i] = f.result()
            except Exception as e:  # noqa
                errs[i] = f'{items[i] if not isinstance(items[i], tuple) else items[i][1:3]}: {str(e)[-400:]}'
            if n % 10 == 0 or n == len(items):
                print(f'[{desc}] {n}/{len(items)} done, {len(errs)} failed', flush=True)
    if errs:
        print(f'[{desc}] FAILURES:'); [print('  ', v) for v in errs.values()]
    return out, errs


class _LightFM:
    """Picklable stand-in for FileManager inside worker processes (the real one holds a Google API client).
    Borrows FileManager's own rclone methods, which only need localMasterDir / cloudMasterDir."""
    def __init__(self, localMasterDir, cloudMasterDir, localBamfilesDir):
        self.localMasterDir, self.cloudMasterDir, self.localBamfilesDir = localMasterDir, cloudMasterDir, localBamfilesDir

    def downloadData(self, *a, **k):
        from helper_modules.file_manager import FileManager
        return FileManager.downloadData(self, *a, **k)

    def returnCloudFiles(self, *a, **k):
        from helper_modules.file_manager import FileManager
        return FileManager.returnCloudFiles(self, *a, **k)

    def checkCloudFile(self, *a, **k):
        from helper_modules.file_manager import FileManager
        return FileManager.checkCloudFile(self, *a, **k)


class Project:
    def __getstate__(self):
        st = dict(self.__dict__)
        fm = st.pop('fm')
        st['_light'] = (fm.localMasterDir, fm.cloudMasterDir, fm.localBamfilesDir)
        return st

    def __setstate__(self, st):
        light = st.pop('_light')
        self.__dict__.update(st)
        self.fm = _LightFM(*light)

    def __init__(self, cfg):
        d = PKG_ROOT                                           # find the repo folder that contains helper_modules/ by walking up
        for _ in range(6):
            d = os.path.dirname(d)
            if os.path.isdir(os.path.join(d, 'helper_modules')):
                sys.path.insert(0, d); break
        else:
            sys.path.insert(0, os.getcwd())
        from helper_modules.file_manager import FileManager
        self.cfg = cfg
        self.fm = FileManager()
        self.fm.setGenome(cfg['genome_version'])
        self.genome = cfg['genome_version']
        self.cores = cfg.get('cores') or os.cpu_count()
        proj_root = self.fm.localMasterDir + cfg['project_dir_name'] + '/'
        base = proj_root + self.genome + '/'
        self.genomes_dir = proj_root + 'genomes/'                              # every genome used in this project (shared by all analyses)
        self.cloud_genomes_dir = self.fm.cloudMasterDir + cfg['project_dir_name'] + '/genomes/'
        os.makedirs(self.genomes_dir, exist_ok=True)
        self.root = base
        self.cloud_root = self.fm.cloudMasterDir + cfg['project_dir_name'] + '/' + self.genome + '/'
        self.dirs = {k: base + k + '/' for k in ['qc', 'original', 'v2', 'comparison', 'logs', 'tmp']}
        for d in self.dirs.values():
            os.makedirs(d, exist_ok=True)
        self.genome_fa = self.fm.localGenomeFile
        self._chr = None

    def run(self, analysis, step, params=None, tar=()):
        from .runner import Run
        return Run(self, analysis, step, params, tar)

    def workers(self, step):
        return int(self.cfg.get('workers', {}).get(step, max(1, self.cores // 2)))

    # ---- chromosomes
    @property
    def chrom_table(self):
        if self._chr is None:
            p = self.dirs['qc'] + 'chromosomes.tsv'
            if os.path.exists(p):
                self._chr = pd.read_csv(p, sep='\t')
            else:
                from . import genome
                self._chr = genome.build(self.genome_fa, os.path.dirname(self.genome_fa) + '/', self.cfg, p)
        return self._chr

    def _sub(self):
        ct = self.chrom_table
        s = self.cfg.get('chrom_subset')
        return ct[ct.LG.isin(s)] if s else ct

    @property
    def chroms(self):
        return self._sub().chrom.tolist()

    @property
    def lg(self):
        return dict(zip(self.chrom_table.chrom, self.chrom_table.LG))

    def chrom_for_lg(self, name):
        m = {v: k for k, v in self.lg.items()}
        return m[name]

    def chrom_lengths(self):
        sub = self._sub()
        return dict(zip(sub.chrom, sub.length))

    # ---- samples
    def sample_table(self):
        self.fm.readSampleDatabase()
        dt = self.fm.sample_dt.copy()
        par = dt.Parents.fillna('').str.split(',', expand=True)
        dt['Sire'] = par[0].str.strip() if 0 in par else None
        dt['Dam'] = par[1].str.strip() if 1 in par else None
        dt['SireSample'] = dt.Sire.map(self.cfg['f1_sire_names'])
        return dt

    def bc_table(self):
        dt = self.sample_table()
        bc = dt[dt.Species == self.cfg['bc_species']][['SampleID', 'BroodID', 'Sire', 'SireSample', 'Dam', 'Sex']].reset_index(drop=True)
        if self.cfg.get('sample_subset'):                       # test mode: only these BC fish
            bc = bc[bc.SampleID.isin(self.cfg['sample_subset'])].reset_index(drop=True)
        return bc

    # ---- files (built directly, not from FileManager state, so they are safe in worker processes)
    def sample_dir(self, s):
        return f"{self.fm.localBamfilesDir}{self.genome}/{s}/"

    def bam(self, s):
        return self.sample_dir(s) + s + '.all.bam'

    def gvcf(self, s):
        return self.sample_dir(s) + s + '.g.vcf.gz'

    def ensure_local(self, path):
        if not os.path.exists(path):
            os.makedirs(os.path.dirname(path), exist_ok=True)
            self.fm.downloadData(path)
        return path

    def ensure_bam(self, s):
        b = self.ensure_local(self.bam(s))
        if not (os.path.exists(b + '.bai') or os.path.exists(b[:-4] + '.bai')):
            try:
                self.fm.downloadData(b + '.bai')
            except FileNotFoundError:
                sh(['samtools', 'index', '-@', '4', b])
        return b

    def ensure_gvcf(self, s):
        """GVCF from the cloud if present; otherwise call it from the BAM with HaplotypeCaller (deep samples only)."""
        g = self.gvcf(s)
        if not os.path.exists(g):
            try:
                self.fm.downloadData(g)
            except FileNotFoundError:
                mem = self.cfg.get('gatk_mem_gb', 8)
                sh(['gatk', '--java-options', f'-Xmx{mem}g', 'HaplotypeCaller', '-R', self.genome_fa, '-I', self.ensure_bam(s),
                    '-O', g, '-ERC', 'GVCF'])
        if not (os.path.exists(g + '.tbi') or os.path.exists(g + '.idx')):
            try:
                self.fm.downloadData(g + '.tbi')
            except FileNotFoundError:
                sh(['tabix', '-f', '-p', 'vcf', g])
        return g

    def load_pheno(self):
        dt = pd.read_csv(self.cfg['pheno_sheet_csv'])
        p = dt[(dt.species == 'MCYHBC1') & (dt.tagtype == 'PIT')][['samplename', 'sex', 'sex_100', 'DOF', 'implant_date', 'Dead', 'death_date']].copy()
        p['sex_100'] = p.sex_100.fillna('-')
        return p
