"""Depth in the gsdf region from a (region) BAM, and the numbers the call is made from. Pure functions plus command builders; the same code serves a streamed region and a local BAM."""
import numpy as np


def region_string(cfg):
    r = cfg['region']
    return f"{r['chrom']}:{r['start'] + 1}-{r['end']}"                         # samtools regions are 1-based, inclusive


def regions(cfg):
    """The regions of the BAM this module reads (for the retrospective mode, which fetches only these)."""
    return [region_string(cfg)]


def view_command(src, region, out_bam):
    return ['samtools', 'view', '-b', '-o', out_bam, src, region]


def index_command(bam):
    return ['samtools', 'index', bam]


def depth_command(bam, region, mapq=None):
    """Per-position depth over the region, zeros included. samtools depth excludes unmapped, secondary, QC-fail and DUPLICATE reads by default; no MAPQ filter unless mapq is given."""
    return ['samtools', 'depth', '-a', '-r', region] + (['-Q', str(mapq)] if mapq else []) + [bam]


def parse_depth(text, start, end):
    """samtools depth output -> float array over [start, end) (0-based); positions that are missing count as 0."""
    arr = np.zeros(end - start, dtype=np.float64)
    for line in text.splitlines():
        f = line.split('\t')
        if len(f) < 3:
            continue
        i = int(f[1]) - 1 - start
        if 0 <= i < len(arr):
            arr[i] = float(f[2])
    return arr


def bin_means(arr, bin_bp):
    n = len(arr) // bin_bp
    return arr[:n * bin_bp].reshape(n, bin_bp).mean(axis=1)


def window_mean_of_bins(bins, region_start, bin_bp, windows):
    """Mean of the bins whose START lies in any window [w0, w1] (inclusive; genome coordinates): sexloci's rule, the one the calibration was made with."""
    keep = [v for k, v in enumerate(bins) if any(w0 <= region_start + k * bin_bp <= w1 for w0, w1 in windows)]
    return float(np.mean(keep)) if keep else float('nan')


def mean_of_bins_inside(bins, region_start, bin_bp, window):
    """Mean of the bins that lie entirely inside the window [w0, w1) (used for the junction ratio's flank)."""
    w0, w1 = window
    keep = [v for k, v in enumerate(bins) if region_start + k * bin_bp >= w0 and region_start + (k + 1) * bin_bp <= w1]
    return float(np.mean(keep)) if keep else float('nan')


def raw_measurements(depth_all, depth_q, cfg):
    """What `depth` stores per sample: all-read bins, MAPQ >= mapq_unique bins, and the junction window's mean depth (all reads / unique reads)."""
    r = cfg['region']; b = cfg['bin_bp']
    j0, j1 = cfg['junction_window']
    return dict(bins=bin_means(depth_all, b).tolist(), bins_q=bin_means(depth_q, b).tolist(),
                junction_depth=float(depth_all[j0 - r['start']: j1 - r['start']].mean()), junction_depth_q=float(depth_q[j0 - r['start']: j1 - r['start']].mean()))


def metrics(raw, cfg):
    """Unit / flank depth, the depth ratio, and the unique-read junction ratio, from stored raw measurements (bins as a list, or the flat columns of the TSV)."""
    r = cfg['region']; b = cfg['bin_bp']
    unit = window_mean_of_bins(raw['bins'], r['start'], b, [cfg['unit_window']])
    flank = window_mean_of_bins(raw['bins'], r['start'], b, cfg['flank_windows'])
    flank_q = mean_of_bins_inside(raw['bins_q'], r['start'], b, cfg['junction_flank_window'])
    return dict(unit_depth=unit, flank_depth=flank, ratio=unit / flank if flank > 0 else float('nan'), flank_depth_unique=flank_q, junction_depth_unique=raw['junction_depth_q'],
                junction_ratio=raw['junction_depth_q'] / flank_q if flank_q > 0 else float('nan'))


# ---------------------------------------------------------------------------------------------------- measuring a local BAM
MEASUREMENT_KEYS = ('region', 'bin_bp', 'window_rule', 'unit_window', 'flank_windows', 'junction_flank_window', 'junction_window', 'mapq_unique')


def definition_key(cfg):
    """Hash of everything that defines WHAT is measured (not how it is called): stored with every measurement; `gsdf call` refuses to mix measurements made under different definitions."""
    import hashlib, json
    return hashlib.md5(json.dumps({k: cfg[k] for k in MEASUREMENT_KEYS}, sort_keys=True).encode()).hexdigest()[:12]


def contigs_in_header(text):
    out = {}
    for line in text.splitlines():
        if line.startswith('@SQ'):
            kv = dict(f.split(':', 1) for f in line.split('\t')[1:] if ':' in f)
            if 'SN' in kv:
                out[kv['SN']] = int(kv.get('LN', 0))
    return out


def measure_bam(bam, cfg, run):
    """Raw gsdf measurements of one coordinate-sorted, indexed BAM aligned to the project reference. `run(cmd)` returns an object with .stdout (text) and raises on failure.
    Two samtools depth passes over the 130 kb region: all reads (duplicates excluded by samtools' default flags) and MAPQ >= mapq_unique."""
    r = cfg['region']
    contigs = contigs_in_header(run(['samtools', 'view', '-H', bam]).stdout)
    if r['chrom'] not in contigs:
        raise ValueError(f"{bam}: no contig {r['chrom']} in the header: is this BAM aligned to the project reference (Mzebra_GT3_NCBI)?")
    reg = region_string(cfg)
    d_all = parse_depth(run(depth_command(bam, reg)).stdout, r['start'], r['end'])
    d_q = parse_depth(run(depth_command(bam, reg, cfg['mapq_unique'])).stdout, r['start'], r['end'])
    out = raw_measurements(d_all, d_q, cfg)
    out.update(definition=definition_key(cfg), contig_length=contigs[r['chrom']], region=reg)
    return out
