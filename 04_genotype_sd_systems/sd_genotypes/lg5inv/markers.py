"""LG5 inversion markers: sites where BOTH inverted haplotypes (LatrZW, MezeOBm) differ from GT3 in the same way, then tested in pools of known composition.

  python -m sd_genotypes.lg5inv.markers candidates --target-fa gt3_lg5.fa --paf LatrZW=LatrZW_vs_GT3.cs.paf --paf MezeOBm=MezeOBm_vs_GT3.cs.paf --out candidates.tsv
  python -m sd_genotypes.lg5inv.markers count --candidates candidates.tsv --ref GT3.fna --slices DIR --runs RUN,RUN,... --out counts.tsv
  python -m sd_genotypes.lg5inv.markers select --counts counts.tsv --negative RUN,.. --positive RUN,.. --out markers.tsv

Why reads decide: GT3 is the only standard-arrangement sequence, so a site where both inverted haplotypes agree against GT3 is as likely a GT3-private allele as an inversion allele. The BB pools tell them apart
(inversion allele: absent in BB, ~the inverted fraction in OB pools; GT3-private allele: fixed in BB). `count` uses NO base-quality filter (SRR18083730 came back from ENA with all qualities Q0) and by default MAPQ >= 20.
Only SNVs are counted from reads in this version; indels are listed in the candidates but not counted.

The alignments are `minimap2 --cs` (short form) of an assembly region against `gt3_lg5.fa` = NC_135171.1:19400001-26980000, renamed LG5_off19400000. The cs string describes the alignment on the forward strand of the
target, so alt alleles are already in GT3 orientation for minus-strand alignments (checked on the reference side by the self-test of every site; the alt side is checked by the pools: an inversion allele is seen in OB reads)."""
import argparse, collections, os, re, subprocess, sys, tempfile
import numpy as np
import pandas as pd

CHROM = 'NC_135171.1'
rc = lambda x: x.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def faidx(ref, a, b):
    """GT3 NC_135171.1 a..b (1-based, inclusive), uppercase, from the reference FASTA."""
    out = subprocess.check_output(['samtools', 'faidx', ref, f'{CHROM}:{a}-{b}']).decode().split('\n')[1:]
    return ''.join(out).upper()
CS = re.compile(r':(\d+)|=([ACGTNacgtn]+)|\*([acgtn])([acgtn])|\+([acgtn]+)|-([acgtn]+)')


def parse_cs(cs):
    """cs (short or long form) -> [('=', n) | ('X', refbase, qbase) | ('I', seq) | ('D', seq)]; raises if any part of the string is not understood."""
    ops, pos = [], 0
    while pos < len(cs):
        m = CS.match(cs, pos)
        if not m:
            raise ValueError(f'cannot parse cs at {pos}: {cs[pos:pos + 20]!r}')
        pos = m.end()
        if m.group(1) is not None:
            ops.append(('=', int(m.group(1))))
        elif m.group(2) is not None:
            ops.append(('=', len(m.group(2))))
        elif m.group(3) is not None:
            ops.append(('X', m.group(3).upper(), m.group(4).upper()))
        elif m.group(5) is not None:
            ops.append(('I', m.group(5).upper()))
        else:
            ops.append(('D', m.group(6).upper()))
    return ops


def normalise(i, s, ref):
    """Left-align an indel. i = 0-based index of the anchor base (the base before the inserted or deleted bases), s = the inserted or deleted bases. Returns (i, s), i = -1 if it cannot be anchored."""
    while i >= 0 and s[-1] == ref[i]:
        s = ref[i] + s[:-1]; i -= 1
    return i, s


def alignment_variants(rec, ref, offset, min_flank=12):
    """One PAF record -> list of (anchor0, (pos, ref, alt, type)). Every SNV/indel needs >= min_flank identical bases on both sides (drops clusters in badly aligned stretches)."""
    f = rec.rstrip('\n').split('\t')
    tstart, tend = int(f[7]), int(f[8])
    cs = next((x[5:] for x in f[12:] if x.startswith('cs:Z:')), None)
    if cs is None:
        raise ValueError('PAF record without a cs tag: align with minimap2 --cs')
    ops = parse_cs(cs)
    t, out = tstart, []
    consumed = sum(o[1] if o[0] == '=' else 1 if o[0] == 'X' else len(o[1]) if o[0] == 'D' else 0 for o in ops)
    if consumed != tend - tstart:
        raise ValueError(f'cs spans {consumed} target bases, the PAF says {tend - tstart} ({f[0]}:{f[2]})')
    for k, o in enumerate(ops):
        before = ops[k - 1][1] if k > 0 and ops[k - 1][0] == '=' else 0
        after = ops[k + 1][1] if k + 1 < len(ops) and ops[k + 1][0] == '=' else 0
        ok = before >= min_flank and after >= min_flank
        if o[0] == '=':
            t += o[1]
        elif o[0] == 'X':
            if ref[t] != o[1]:
                raise ValueError(f'cs says the target base at index {t} is {o[1]} but the target sequence has {ref[t]}: wrong target file or offset')
            if ok and o[2] != 'N' and o[1] != 'N':
                out.append((t, (offset + t + 1, o[1], o[2], 'SNV')))
            t += 1
        elif o[0] == 'I':
            if ok and t - 1 >= 0 and 'N' not in o[1]:
                i, s = normalise(t - 1, o[1], ref)
                if i >= 0:
                    out.append((i, (offset + i + 1, ref[i], ref[i] + s, 'INS')))
        else:
            if ref[t:t + len(o[1])] != o[1]:
                raise ValueError(f'cs deletes {o[1]} at target index {t} but the target sequence has {ref[t:t + len(o[1])]}')
            if ok and t - 1 >= 0 and 'N' not in o[1]:
                i, s = normalise(t - 1, o[1], ref)
                if i >= 0:
                    out.append((i, (offset + i + 1, ref[i] + s, ref[i], 'DEL')))
            t += len(o[1])
    return out


def assembly_variants(paf_lines, ref, offset, min_flank=12, max_div=0.05, min_block=500):
    """-> (variants {(pos, ref, alt): (anchor0, type)}, coverage array, summary). Coverage counts EVERY alignment of >= min_block (a site covered twice is a duplicated or ambiguous stretch of the assembly);
    variants come only from alignments with divergence <= max_div."""
    cov = np.zeros(len(ref) + 1, dtype=np.int32); var = {}; n = n_div = n_short = 0
    for line in paf_lines:
        if not line.strip():
            continue
        f = line.rstrip('\n').split('\t'); n += 1
        if int(f[10]) < min_block:
            n_short += 1; continue
        cov[int(f[7]):int(f[8])] += 1
        de = next((float(x[5:]) for x in f[12:] if x.startswith('de:f:')), 0.0)
        if de > max_div:
            n_div += 1; continue
        for a0, key in alignment_variants(line, ref, offset, min_flank):
            var[key[:3]] = (a0, key[3])
    return var, cov, dict(alignments=n, short=n_short, too_divergent=n_div, variants=len(var))


def read_fasta_single(path):
    name, parts = None, []
    for line in open(path):
        if line.startswith('>'):
            if name is not None:
                raise ValueError(f'{path} has more than one sequence')
            name = line[1:].split()[0]
        else:
            parts.append(line.strip())
    m = re.search(r'_off(\d+)$', name or '')
    if not m:
        raise ValueError(f'target name {name!r} must end in _off<offset>, e.g. LG5_off19400000')
    return ''.join(parts).upper(), int(m.group(1))


def shared_candidates(assemblies, ref, offset, **kw):
    """assemblies = {name: paf lines}. A candidate is a variant seen with the same position, ref and alt in EVERY assembly, at a site covered by exactly one alignment of each. Returns (DataFrame, summary)."""
    per = {n: assembly_variants(lines, ref, offset, **kw) for n, lines in assemblies.items()}
    names = list(per)
    keys = set.intersection(*[set(per[n][0]) for n in names])
    rows = []
    for k in sorted(keys):
        a0 = [per[n][0][k][0] for n in names]
        if all(per[n][1][a] == 1 for n, a in zip(names, a0)):
            rows.append((CHROM, k[0], k[1], k[2], per[names[0]][0][k][1], abs(len(k[1]) - len(k[2]))))
    df = pd.DataFrame(rows, columns=['chrom', 'pos', 'ref', 'alt', 'type', 'indel_len'])
    summary = {n: per[n][2] for n in names}
    summary['shared_before_coverage_filter'] = len(keys); summary['candidates'] = len(df)
    summary['by_type'] = df.type.value_counts().to_dict() if len(df) else {}
    return df, summary


def pileup_counts(bases, refbase):
    """Counts of A/C/G/T/* in a samtools mpileup base string. '.' and ',' are the reference base; ^q (read start + MAPQ char), $ and +N/-N indel sequences are skipped."""
    c = collections.Counter(); i, n = 0, len(bases)
    while i < n:
        ch = bases[i]
        if ch == '^':
            i += 2; continue
        if ch == '$':
            i += 1; continue
        if ch in '+-':
            k = i + 1
            while k < n and bases[k].isdigit():
                k += 1
            i = k + int(bases[i + 1:k]); continue
        if ch in '.,':
            c[refbase] += 1
        elif ch in '*#':
            c['*'] += 1
        elif ch in '<>':
            pass
        else:
            c[ch.upper()] += 1
        i += 1
    return c


def run_mpileup(ref_fna, bam, bed_path, min_mapq):
    cmd = ['samtools', 'mpileup', '-a', '-B', '-Q', '0', '-q', str(min_mapq), '-d', '100000', '-f', ref_fna, '-l', bed_path, bam]
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines()


def count_sites(snvs, bam, ref_fna, min_mapq=20, runner=run_mpileup):
    """snvs: DataFrame with pos, ref, alt (SNVs). Returns (depth array, alt-count array, ref-count array) in the order of snvs. No base-quality filter. A site whose reference base in the FASTA differs from
    the candidate's ref raises (wrong reference)."""
    with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as fh:
        for p in snvs.pos:
            fh.write(f'{CHROM}\t{int(p) - 1}\t{int(p)}\n')
        bed = fh.name
    try:
        lines = runner(ref_fna, bam, bed, min_mapq)
    finally:
        os.remove(bed)
    got = {}
    for line in lines:
        f = line.split('\t')
        if len(f) >= 5:
            got[int(f[1])] = (f[2].upper(), f[4])
    depth, alt, refc = [], [], []
    for p, r, a in zip(snvs.pos, snvs.ref, snvs.alt):
        rb, bases = got.get(int(p), (r, ''))
        if rb != r:
            raise ValueError(f'reference base at {CHROM}:{p} is {rb} in the FASTA but {r} in the candidates: wrong reference')
        c = pileup_counts(bases, r)
        refc.append(c[r]); alt.append(c[a]); depth.append(sum(c[x] for x in 'ACGT'))
    return np.array(depth), np.array(alt), np.array(refc)


def count_pools(cands, slicedir, runs, ref_fna, min_mapq=20, runner=run_mpileup):
    snv = cands[cands.type == 'SNV'].reset_index(drop=True)
    for r in runs:
        d, a, rc_ = count_sites(snv, os.path.join(slicedir, r + '.slice.bam'), ref_fna, min_mapq, runner)
        snv[f'{r}_depth'] = d; snv[f'{r}_alt'] = a
    return snv


def select(table, negative, positive, min_depth=15, neg_max_alt=1, pos_frac=(0.2, 0.8)):
    """Marker = depth >= min_depth in every listed pool, <= neg_max_alt alt reads in every negative pool, alt fraction within pos_frac in every positive pool."""
    t = table.copy(); ok = pd.Series(True, index=t.index)
    for r in negative + positive:
        ok &= t[f'{r}_depth'] >= min_depth
    for r in negative:
        ok &= t[f'{r}_alt'] <= neg_max_alt
    fr = {}
    for r in positive:
        fr[r] = t[f'{r}_alt'] / t[f'{r}_depth'].clip(lower=1)
        ok &= (fr[r] >= pos_frac[0]) & (fr[r] <= pos_frac[1])
    t['mean_positive_fraction'] = sum(fr.values()) / len(positive) if positive else np.nan
    return t[ok].reset_index(drop=True)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.lg5inv.markers', description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('candidates'); c.add_argument('--target-fa', required=True); c.add_argument('--paf', action='append', required=True, help='NAME=file.cs.paf (repeat)'); c.add_argument('--out', required=True)
    c.add_argument('--min-flank', type=int, default=12); c.add_argument('--max-div', type=float, default=0.05); c.add_argument('--min-block', type=int, default=500)
    n = sub.add_parser('count'); n.add_argument('--candidates', required=True); n.add_argument('--ref', required=True); n.add_argument('--slices', required=True); n.add_argument('--runs', required=True)
    n.add_argument('--out', required=True); n.add_argument('--min-mapq', type=int, default=20)
    s = sub.add_parser('select'); s.add_argument('--counts', required=True); s.add_argument('--negative', required=True); s.add_argument('--positive', required=True); s.add_argument('--out', required=True)
    s.add_argument('--min-depth', type=int, default=15); s.add_argument('--neg-max-alt', type=int, default=1); s.add_argument('--pos-min', type=float, default=0.2); s.add_argument('--pos-max', type=float, default=0.8)
    a = ap.parse_args(argv)
    if a.cmd == 'candidates':
        ref, offset = read_fasta_single(a.target_fa)
        asm = {}
        for spec in a.paf:
            name, path = spec.split('=', 1); asm[name] = open(path).read().splitlines()
        df, summary = shared_candidates(asm, ref, offset, min_flank=a.min_flank, max_div=a.max_div, min_block=a.min_block)
        df.to_csv(a.out, sep='\t', index=False)
        for k, v in summary.items():
            print(f'  {k}: {v}')
        print(f'wrote {len(df)} candidates to {a.out}')
    elif a.cmd == 'count':
        cands = pd.read_csv(a.candidates, sep='\t'); runs = a.runs.split(',')
        t = count_pools(cands, a.slices, runs, a.ref, a.min_mapq)
        t.to_csv(a.out, sep='\t', index=False); print(f'counted {len(t)} SNV candidates in {len(runs)} pools -> {a.out}')
    else:
        t = pd.read_csv(a.counts, sep='\t'); neg, pos = a.negative.split(','), a.positive.split(',')
        m = select(t, neg, pos, a.min_depth, a.neg_max_alt, (a.pos_min, a.pos_max)); m.to_csv(a.out, sep='\t', index=False)
        print(f'{len(m)} of {len(t)} SNV candidates pass; positions span {int(m.pos.min()) if len(m) else "-"}-{int(m.pos.max()) if len(m) else "-"}; wrote {a.out}')


if __name__ == '__main__':
    main()
