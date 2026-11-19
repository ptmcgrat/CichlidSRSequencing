"""Loading what `genotype hmm` wrote, and the table of fish with their sex, sire, dam and brood."""
import hashlib
import os
import re

import numpy as np
import pandas as pd


def lg_number(label, fallback):
    m = re.fullmatch(r'LG(\d+)', str(label))
    return int(m.group(1)) if m else fallback


class Calls:
    """The hard calls of one HMM run: meta (chrom, LG label, lg number, pos), X (bins x fish, int8: 0 = the MC state, 1 = the state that carries the YH allele), fish, and the bin size.
    Rows are sorted by linkage group, then position."""

    def __init__(self, meta, X, fish, bin_bp=100_000, path=None):
        self.meta, self.X, self.fish, self.bin_bp, self.path = meta, X, list(fish), int(bin_bp), path
        self._idx = {c: np.where(meta.chrom.to_numpy() == c)[0] for c in dict.fromkeys(meta.chrom)}
        self.chroms = list(self._idx)

    @classmethod
    def load(cls, path, lg_label=lambda c: c, bin_bp=100_000):
        d = pd.read_csv(path, index_col=0)
        if list(d.columns[:2]) != ['Chromosome', 'Position']:
            raise ValueError(f'{path}: the first columns must be Chromosome and Position (a genotypes_calls.csv of `genotype hmm`)')
        fish = list(d.columns[2:])
        labels = [lg_label(c) for c in d.Chromosome]
        order_of_first_seen = {c: i for i, c in enumerate(dict.fromkeys(d.Chromosome))}
        lg = np.array([lg_number(l, 1000 + order_of_first_seen[c]) for l, c in zip(labels, d.Chromosome)])
        meta = pd.DataFrame({'chrom': d.Chromosome.to_numpy(), 'LG': labels, 'lg': lg, 'pos': d.Position.to_numpy()})
        order = np.lexsort((meta.pos.to_numpy(), meta.lg.to_numpy()))
        vals = d[fish].to_numpy(dtype=float)
        if not np.isin(vals, (0, 1)).all():
            raise ValueError(f'{path}: the calls must be 0 or 1 (the hard calls), not probabilities')
        X = vals.astype(np.int8)[order]
        return cls(meta.iloc[order].reset_index(drop=True), X, fish, bin_bp, path)

    def idx(self, chrom):
        return self._idx[chrom]

    def pos(self, chrom):
        return self.meta.pos.to_numpy()[self._idx[chrom]]

    def length(self, chrom):
        return int(self.pos(chrom).max() + self.bin_bp)

    def lg_of(self, chrom):
        return self.meta.LG.to_numpy()[self._idx[chrom][0]]

    def chrom_of(self, lg):
        """The chromosome of a linkage group given as 10, 'LG10' or an accession."""
        key = f'LG{lg}' if str(lg).isdigit() else str(lg)
        for c in self.chroms:
            if c == key or self.lg_of(c) == key:
                return c
        raise ValueError(f'no chromosome {lg!r} in the calls ({[self.lg_of(c) for c in self.chroms][:6]} ...)')


def fish_table(df):
    """SampleID-indexed table with Sex ('M' or 'F'; anything else, 'M?' and 'U' included, is NaN: not used as a phenotype), Sire, Dam and BroodID where given."""
    d = df.copy()
    if 'SampleID' in d.columns:
        d = d.drop_duplicates('SampleID').set_index('SampleID')
    d.index = d.index.astype(str)
    if 'Sex' in d.columns:
        d['Sex'] = d['Sex'].astype(str).str.strip().where(lambda s: s.isin(['M', 'F']))
    return d


def sex_vector(calls, ft):
    """1 for M, 0 for F, NaN for a fish with no usable sex, in the order of calls.fish."""
    if 'Sex' not in ft.columns:
        raise ValueError('the fish table has no Sex column')
    return ft.reindex(calls.fish)['Sex'].map({'M': 1.0, 'F': 0.0}).to_numpy(dtype=float)


def groups_of(calls, ft, col):
    return ft.reindex(calls.fish)[col] if col in ft.columns else pd.Series(index=calls.fish, dtype=object)


def read_mask(path):
    """{(chrom, bin index)} of the masked 100 kb bins of `genotype mask` (empty when there is no file)."""
    if not path or not os.path.exists(path):
        return set()
    m = pd.read_csv(path, sep='\t')
    return {(r.chrom, int(r.bin)) for r in m.itertuples()}


def sha1(path, chunk=1 << 20):
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h.update(b)
    return h.hexdigest()
