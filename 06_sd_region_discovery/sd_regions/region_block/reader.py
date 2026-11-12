"""Reads the VCF ONCE, through sexloci's own reader, variant filter and depth masking, so the genotypes analysed here are the ones sexloci analysed (same QUAL and MQ filters, same `min_dp`, the same call-rate
filter over the samples). One difference, on purpose: the filter of INFO/DP above its 99th percentile is not applied, because it needs a second pass over the whole VCF to find the quantile and removes
about 1% of the sites; the counts therefore differ from sexloci's by that much.

Everything that is needed from the genome is accumulated while streaming (windows, relationship matrices); the sites of the region's chromosome are kept in memory (a few hundred thousand sites x the samples)
for the analyses that need them (PCA, LD, the diagnostic SNPs)."""
import re

import numpy as np

from . import analysis as A


def parse_region(text):
    """'NC_135176.1:10,000,000-30,000,000' -> ('NC_135176.1', 10000000, 30000000)."""
    m = re.fullmatch(r'\s*([^:\s]+):([\d,]+)-([\d,]+)\s*', text or '')
    if not m:
        raise ValueError(f'a region looks like CHROM:START-END, not {text!r}')
    a, b = int(m.group(2).replace(',', '')), int(m.group(3).replace(',', ''))
    if a >= b:
        raise ValueError(f'the region {text!r} ends before it starts')
    return m.group(1), a, b


class Collected:
    """What one pass over the VCF gives."""
    pass


def stream(vcf, names, region, config_path=None, min_dp=3, min_call_rate=0.8, window=250_000, min_maf=0.05, exclude=(), keep_region_chromosome=False, backend=None, log=print):
    """Stream `vcf` for the samples `names` (in that order) and return a Collected with: names, windows (WindowAccumulator), grm_out (GrmAccumulator over the sites outside the region and, unless kept, outside its chromosome),
    chr_pos, chr_raw, chr_dp (the sites of the region's chromosome that passed, with raw genotypes and depths), n_seen, n_used, cfg."""
    from sexloci.config import load_config
    from sexloci.vcf import VariantFilter, iter_chunks, read_header_info, read_samples
    rc, rs, re_ = region
    all_samples = read_samples(vcf)
    missing = [s for s in names if s not in all_samples]
    if missing:
        raise ValueError(f'{len(missing)} samples are not in the VCF: {", ".join(missing[:5])}')
    idx = [all_samples.index(s) for s in names]
    cfg = load_config(config_path)
    cfg.min_dp, cfg.min_call_rate = int(min_dp), float(min_call_rate)
    cfg.variant_depth_quantiles = None                                    # no second pass for the INFO/DP quantile (see the module docstring)
    if backend:
        cfg.backend = backend
    hdr = read_header_info(vcf)
    vf = VariantFilter(cfg, hdr['infos'])
    n = len(names)
    out = Collected()
    out.names, out.cfg, out.window = list(names), cfg, window
    out.windows = A.WindowAccumulator(n, window)
    out.grm_out = A.GrmAccumulator(n, min_maf)
    pos_l, raw_l, dp_l = [], [], []
    out.n_seen = out.n_used = 0
    excl = [(c, a, b) for c, a, b in exclude]
    for k, chunk in enumerate(iter_chunks(vcf, idx, all_samples, cfg, vf)):
        out.n_seen += len(chunk.pos)
        ok = (chunk.dos >= 0).mean(axis=1) >= cfg.min_call_rate               # sexloci's call-rate filter: at least min_call_rate of the samples have a genotype at min_dp
        if not ok.any():
            continue
        c = chunk.subset(ok)
        out.n_used += len(c.pos)
        out.windows.add(c.chrom, c.pos, c.dos)
        on_chr = c.chrom == rc
        in_region = on_chr & (c.pos >= rs) & (c.pos < re_)
        skip = in_region.copy() if keep_region_chromosome else on_chr.copy()
        for ec, ea, eb in excl:
            skip |= (c.chrom == ec) & (c.pos >= ea) & (c.pos < eb)
        out.grm_out.add(c.dos[~skip])
        if on_chr.any():
            pos_l.append(c.pos[on_chr])
            raw_l.append((c.dos_raw if c.dos_raw is not None else c.dos)[on_chr])
            dp_l.append((c.dp if c.dp is not None else np.full(c.dos.shape, np.nan, np.float32))[on_chr])
        if (k + 1) % 20 == 0:
            log(f'  ... {out.n_seen:,} sites read, {out.n_used:,} kept')
    out.chr_pos = np.concatenate(pos_l) if pos_l else np.zeros(0, np.int64)
    out.chr_raw = np.vstack(raw_l) if raw_l else np.zeros((0, n), np.int8)
    out.chr_dp = np.vstack(dp_l) if dp_l else np.zeros((0, n), np.float32)
    return out


def masked(raw, dp, min_dp):
    """Genotypes with depth below min_dp set to missing (the same rule as sexloci's min_dp)."""
    d = raw.copy()
    d[np.nan_to_num(dp, nan=np.inf) < min_dp] = -1
    return d
