"""Build the B-chromosome measurement definition: the core blocks, a systematic sample of 50 bp bins on the placed chromosomes (to estimate the genome-wide depth without a whole-genome pass), and a meta file.
  python -m sd_genotypes.bchrom.definition --core-blocks core_blocks.tsv --block-signal block_signal.tsv --sample-summary sample_summary.tsv --lengths-from track.json|ref.fa.fai --out-dir DIR"""
import argparse, bisect, hashlib, json, os
import numpy as np
import pandas as pd


def read_lengths(path):
    """Sequence lengths from a depth track .json (key `lengths`) or a samtools .fai."""
    if path.endswith('.json'):
        return {k: int(v) for k, v in json.load(open(path))['lengths'].items()}
    out = {}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 2:
            out[f[0]] = int(f[1])
    return out


def _overlaps(iv, starts, s, e):
    i = bisect.bisect_left(starts, e)
    return any(iv[j][1] > s and iv[j][0] < e for j in range(max(0, i - 3), i))


def sample_scale_bins(lengths, blocks, n_bins=12000, bin_bp=50, min_chrom_bp=5_000_000):
    """Evenly spaced full bins on every placed sequence (length >= min_chrom_bp), their number proportional to its length, skipping bins that overlap a core block. Deterministic: no random numbers."""
    placed = {c: L for c, L in lengths.items() if L >= min_chrom_bp}
    if not placed:
        raise ValueError('no sequence is long enough to be a placed chromosome')
    tot = sum(placed.values()); rows = []
    by = {c: sorted(zip(g.start.astype(int), g.end.astype(int))) for c, g in blocks.groupby('chrom')}
    for c, L in sorted(placed.items()):
        nb = L // bin_bp; n = max(1, int(round(n_bins * L / tot))); step = nb / n
        iv = by.get(c, []); starts = [s for s, e in iv]; prev = -1
        for i in range(n):
            b = max(int((i + 0.5) * step), prev + 1)
            while b < nb and _overlaps(iv, starts, b * bin_bp, (b + 1) * bin_bp):
                b += 1
            if b >= nb:
                break
            rows.append((c, b * bin_bp, (b + 1) * bin_bp)); prev = b
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end'])


def clean_flags(core, block_signal, sample_summary, lo=0.7, hi=1.3):
    """A block is `clean` when its median scaled coverage over the B-free samples (roles baseline and control) lies in [lo, hi)."""
    neg = set(sample_summary[sample_summary.role.isin(['baseline', 'control'])].run)
    cols = [c for c in block_signal.columns if c.split('|')[-1] in neg]
    if not cols:
        raise ValueError('no B-free samples found in block_signal.tsv (columns are label|run)')
    med = block_signal[cols].median(axis=1)
    return core.core_id.map(med).between(lo, hi, inclusive='left')


def build(core, lengths, out_dir, block_signal=None, sample_summary=None, n_bins=12000, bin_bp=50, min_chrom_bp=5_000_000, source=None):
    blocks = core[['core_id', 'chrom', 'start', 'end', 'LG', 'tier']].rename(columns={'core_id': 'block_id'}).copy()
    blocks['clean'] = clean_flags(core.set_index('core_id', drop=False), block_signal.set_index('core_id'), sample_summary).values if block_signal is not None else True
    bins = sample_scale_bins(lengths, blocks, n_bins, bin_bp, min_chrom_bp)
    missing = sorted(set(blocks.chrom) - set(lengths))
    if missing:
        raise ValueError(f'core blocks on sequences that are not in the lengths file: {missing}')
    used = {c: lengths[c] for c in sorted(set(bins.chrom) | set(blocks.chrom))}
    os.makedirs(out_dir, exist_ok=True)
    bp, np_, mp = (os.path.join(out_dir, f) for f in ('bchrom_blocks.tsv', 'bchrom_scale_bins.bed', 'bchrom_meta.json'))
    blocks.to_csv(bp, sep='\t', index=False); bins.to_csv(np_, sep='\t', index=False, header=False)
    json.dump(dict(lengths=used, bin_bp=bin_bp, n_scale_bins=len(bins), min_chrom_bp=min_chrom_bp, n_blocks=len(blocks), sampling='evenly spaced full bins per placed sequence, proportional to length, skipping core blocks',
                   source=source or {}), open(mp, 'w'), indent=1)
    return bp, np_, mp, len(blocks), len(bins)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.bchrom.definition')
    ap.add_argument('--core-blocks', required=True); ap.add_argument('--block-signal'); ap.add_argument('--sample-summary'); ap.add_argument('--lengths-from', required=True); ap.add_argument('--out-dir', required=True)
    ap.add_argument('--n-bins', type=int, default=12000)
    a = ap.parse_args(argv)
    md5 = lambda p: hashlib.md5(open(p, 'rb').read()).hexdigest()[:12]
    core = pd.read_csv(a.core_blocks, sep='\t')
    bsig = pd.read_csv(a.block_signal, sep='\t') if a.block_signal else None
    ssum = pd.read_csv(a.sample_summary, sep='\t') if a.sample_summary else None
    src = {'core_blocks_md5': md5(a.core_blocks), 'lengths_from': os.path.basename(a.lengths_from)}
    bp, np_, mp, nb, nbins = build(core, read_lengths(a.lengths_from), a.out_dir, bsig, ssum, a.n_bins, source=src)
    print(f'{nb} blocks ({int((core.tier == "A").sum())} tier A), {nbins} scale bins -> {a.out_dir}; md5 ' + ', '.join(f'{os.path.basename(p)} {md5(p)}' for p in (bp, np_, mp)))


if __name__ == '__main__':
    main()
