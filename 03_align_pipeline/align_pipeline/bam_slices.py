"""Region slices of the final BAM (a few MB per sample) for later genotyping without going back to the full BAMs.

A slice file is a tab-separated table (see slice_regions/Mzebra_GT3_NCBI.tsv): name, chrom, start, end, note. Coordinates are 1-based and inclusive, `chrom` is the sequence name in the
reference the BAMs are aligned to; `end` may be the word END, which stands for the length of `chrom` taken from the reference index (.fai), so a whole chromosome needs no typed length. Output: <sample dir>/slices/<sample>.slice.<name>.bam (+ .bai)."""
import os
import re

NAME_RE = re.compile(r'^[A-Za-z0-9_.-]+$')


class Slice:
    def __init__(self, name, chrom, start, end, note=''):
        self.name, self.chrom, self.start, self.end, self.note = name, chrom, int(start), int(end), note

    @property
    def region(self):
        return f'{self.chrom}:{self.start}-{self.end}'

    def __repr__(self):
        return f'Slice({self.name} {self.region})'


def fai_lengths(fai_path):
    lengths = {}
    with open(fai_path) as f:
        for line in f:
            p = line.split('\t')
            if len(p) >= 2:
                lengths[p[0]] = int(p[1])
    return lengths


def load_slices(path, fai=None):
    """fai: the path of the reference index or a {sequence: length} dict; needed only when a line says END."""
    lengths = fai_lengths(fai) if isinstance(fai, str) else fai
    slices, seen = [], set()
    with open(path) as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            p = line.split('\t')
            if len(p) < 4:
                raise ValueError(f'{path}:{n}: need name, chrom, start, end (tab separated)')
            end = p[3].strip()
            if end.upper() == 'END':
                if lengths is None:
                    raise ValueError(f'{path}:{n}: END needs the reference index (.fai) to know the length of {p[1]}')
                if p[1] not in lengths:
                    raise ValueError(f'{path}:{n}: sequence {p[1]} is not in the reference index')
                end = lengths[p[1]]
            s = Slice(p[0], p[1], p[2], end, p[4] if len(p) > 4 else '')
            if not NAME_RE.match(s.name):
                raise ValueError(f'{path}:{n}: slice name {s.name!r} may only contain letters, digits, _ . -')
            if s.name in seen:
                raise ValueError(f'{path}:{n}: duplicate slice name {s.name}')
            if s.start < 1 or s.end <= s.start:
                raise ValueError(f'{path}:{n}: {s.name}: need 1 <= start < end')
            seen.add(s.name)
            slices.append(s)
    return slices


def check_against_fai(slices, fai_path):
    """Problems (list of strings) when a slice names a sequence the reference does not have or reaches past its end."""
    lengths = fai_lengths(fai_path)
    problems = []
    for s in slices:
        if s.chrom not in lengths:
            problems.append(f'{s.name}: sequence {s.chrom} is not in the reference')
        elif s.end > lengths[s.chrom]:
            problems.append(f'{s.name}: end {s.end} is past the end of {s.chrom} ({lengths[s.chrom]})')
    return problems


def slice_path(slice_dir, sample, s):
    return os.path.join(slice_dir, f'{sample}.slice.{s.name}.bam')


def view_command(bam, out_bam, s, threads=2):
    return ['samtools', 'view', '-b', '-@', str(threads), '-o', out_bam, bam, s.region]


def index_command(bam, bai, threads=2):
    return ['samtools', 'index', '-@', str(threads), bam, bai]
