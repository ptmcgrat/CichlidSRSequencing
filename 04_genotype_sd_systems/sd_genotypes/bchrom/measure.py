"""B chromosome: raw measurements from ONE coordinate-sorted, indexed BAM aligned to Mzebra_GT3_NCBI (a panel module: measure_bam(bam, cfg, run) and definition_key(cfg)). Only samtools is needed.
Stored per sample: the mean depth of every core block and of a systematic sample of 50 bp bins on the placed chromosomes (their median non-zero value is the genome-wide depth). The call is made later (score.py)."""
import hashlib, json, os, tempfile
import numpy as np
import pandas as pd
from sd_genotypes.gsdf.measure import contigs_in_header

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
MEASUREMENT_KEYS = ('min_mapq', 'blocks_file', 'bins_file', 'meta_file')
MODULE = 'bchrom'


def _paths(cfg):
    d = cfg.get('definition_dir') or DATA
    return [os.path.join(d, cfg[k]) for k in ('blocks_file', 'bins_file', 'meta_file')]


def load_definition(cfg):
    bp, np_, mp = _paths(cfg)
    for p in (bp, np_, mp):
        if not os.path.exists(p):
            raise FileNotFoundError(f'{p} is missing: build the definition with `python -m sd_genotypes.bchrom.definition` (see sd_genotypes/bchrom/data/README.md)')
    return pd.read_csv(bp, sep='\t'), pd.read_csv(np_, sep='\t', header=None, names=['chrom', 'start', 'end']), json.load(open(mp))


def definition_key(cfg):
    """Hash of what is measured: the measurement keys of the config and the CONTENT of the three definition files."""
    load_definition(cfg)                                                   # a missing file is reported with the command that builds it
    h = {k: cfg[k] for k in MEASUREMENT_KEYS}
    h['files'] = [hashlib.md5(open(p, 'rb').read()).hexdigest() for p in _paths(cfg)]
    return hashlib.md5(json.dumps(h, sort_keys=True).encode()).hexdigest()[:12]


def bedcov_command(bam, bed, mapq=0):
    """Total base depth per region (indexed access); default read filters exclude unmapped, secondary, QC-fail and duplicate reads."""
    return ['samtools', 'bedcov'] + (['-Q', str(mapq)] if mapq else []) + [bed, bam]


def regions_depth(run, bam, regions, mapq=0):
    """Mean depth of each region (rows of `regions`: chrom, start, end; 0-based half-open), in order."""
    with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as fh:
        for r in regions.itertuples():
            fh.write(f'{r.chrom}\t{int(r.start)}\t{int(r.end)}\n')
        bed = fh.name
    try:
        out = run(bedcov_command(bam, bed, mapq)).stdout
    finally:
        os.remove(bed)
    sums = {}
    for line in out.splitlines():
        f = line.split('\t')
        if len(f) >= 4:
            sums[(f[0], int(f[1]), int(f[2]))] = float(f[3])
    try:
        return np.array([sums[(r.chrom, int(r.start), int(r.end))] / (int(r.end) - int(r.start)) for r in regions.itertuples()])
    except KeyError as e:
        raise ValueError(f'samtools bedcov returned no row for region {e.args[0]}: wrong BAM or an old samtools?')


def measure_bam(bam, cfg, run):
    """`run(cmd)` returns an object with .stdout (text) and raises on failure (the panel's contract)."""
    blocks, bins, meta = load_definition(cfg)
    contigs = contigs_in_header(run(['samtools', 'view', '-H', bam]).stdout)
    for c, L in meta['lengths'].items():
        if contigs.get(c) != L:
            raise ValueError(f'{bam}: {c} is {contigs.get(c)} bp in the header, the project reference has {L}: not aligned to Mzebra_GT3_NCBI?')
    bd = regions_depth(run, bam, blocks, cfg['min_mapq']); sb = regions_depth(run, bam, bins, cfg['min_mapq'])
    return dict(definition=definition_key(cfg), n_blocks=len(blocks), n_scale_bins=len(bins), block_depth=np.round(bd, 3).tolist(), scale_bins=np.round(sb, 2).tolist())
