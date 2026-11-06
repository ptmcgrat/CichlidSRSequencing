"""Binned depth tracks. mosdepth (--by <bin>) gives the mean depth of every bin with the same read filters as `samtools depth` / mpileup (unmapped, secondary, QC-fail and duplicate reads
excluded; no MAPQ filter, which matters here: B-derived reads mismap onto their A homologs on purpose). A track is one float32 array per sequence, saved compressed (~50 MB per run)."""
import json, os
import numpy as np
import pandas as pd


def mosdepth_command(prefix, bam, bin_size, threads=4):
    return ['mosdepth', '-n', '-t', str(threads), '--by', str(bin_size), prefix, bam]


def parse_regions(path, bin_size):
    """mosdepth <prefix>.regions.bed.gz -> ({sequence: float32 mean depth per bin}, {sequence: length}). Bins must be contiguous multiples of bin_size from 0."""
    parts = {}
    for df in pd.read_csv(path, sep='\t', header=None, names=['chrom', 'start', 'end', 'depth'],
                          dtype={'chrom': str, 'start': np.int64, 'end': np.int64, 'depth': np.float32}, chunksize=2_000_000):
        for chrom, g in df.groupby('chrom', sort=False):
            parts.setdefault(chrom, []).append((g.start.values, g.end.values, g.depth.values))
    arrays, lengths = {}, {}
    for chrom, p in parts.items():
        start = np.concatenate([x[0] for x in p]); end = np.concatenate([x[1] for x in p]); depth = np.concatenate([x[2] for x in p])
        if not np.array_equal(start, np.arange(len(start), dtype=np.int64) * bin_size):
            raise ValueError(f'{path}: bins of {chrom} are not contiguous multiples of {bin_size}; was mosdepth run with --by {bin_size}?')
        arrays[chrom] = depth.astype(np.float32); lengths[chrom] = int(end[-1])
    if not arrays:
        raise ValueError(f'{path}: no bins')
    return arrays, lengths


def scale_chroms(arrays, bin_size, min_len):
    """Sequences used for genome-wide statistics: the placed chromosomes (length >= min_len), not the unplaced scaffolds that collect the repeats."""
    return [c for c, a in arrays.items() if len(a) * bin_size >= min_len]


def depth_stats(arrays, bin_size, min_len):
    cs = scale_chroms(arrays, bin_size, min_len)
    if not cs:
        raise ValueError('no sequence is long enough to be a chromosome; check bchrom.scale_min_chrom_bp')
    v = np.concatenate([arrays[c] for c in cs]); nz = v[v > 0]
    allv_n = sum(len(a) for a in arrays.values()); allv_s = float(sum(float(a.sum(dtype=np.float64)) for a in arrays.values()))
    return dict(mean=float(v.mean()), median=float(np.median(nz)) if len(nz) else 0.0, mean_all_sequences=allv_s / allv_n,
                scale_chroms=len(cs), n_bins=int(len(v)), zero_fraction=float((v == 0).mean()))


def scale_factor(stats, how='median'):
    if how not in ('median', 'mean'):
        raise ValueError("scale_by must be 'median' or 'mean'")
    sf = stats[how]
    if sf <= 0:
        raise ValueError('scale factor is zero: no coverage')
    return sf


def save_track(path, arrays, lengths, stats, bin_size, extra=None):
    """<path>.npz (the arrays) and <path>.json (lengths, statistics, bin size, provenance)."""
    tmp = path + '.partial.npz'
    np.savez_compressed(tmp, **arrays)
    os.replace(tmp, path + '.npz')
    json.dump(dict(bin_size=bin_size, lengths=lengths, stats=stats, **(extra or {})), open(path + '.json', 'w'), indent=1)


def load_track(path):
    z = np.load(path + '.npz')
    arrays = {k: z[k] for k in z.files}
    meta = json.load(open(path + '.json'))
    return arrays, meta
