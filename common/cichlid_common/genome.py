"""Chromosome table for the chosen genome: RefSeq accession <-> LG label <-> length. Written to qc/chromosomes.tsv.
LG labels come from (in order): config chrom_to_lg; the NCBI assembly report; the FASTA headers (RefSeq headers read "... linkage group LG10, ...")."""
import glob, os, re
import pandas as pd
from .common_light import sh


_HDR = re.compile(r'^>(\S+).*?\b(?:linkage group|chromosome)\s+(?:LG|lg)?_?0*(\w+?)\b[,\s]')


def lg_from_headers(header_lines, lengths):
    """{accession: 'LGn'} from FASTA header lines such as '>NC_135176.1 Maylandia zebra isolate ... linkage group LG10, Mzebra_GT3a, whole genome shotgun sequence'.
    Only sequences in `lengths` are used; an unplaced-scaffold header has no LG and is skipped. A label seen twice for different sequences is ambiguous -> dropped."""
    lg, seen = {}, {}
    for line in header_lines:
        m = _HDR.match(line.rstrip('\n') + ' ')
        if not m or m.group(1) not in lengths or 'unplaced' in line:
            continue
        label = 'LG' + m.group(2)
        seen.setdefault(label, []).append(m.group(1)); lg[m.group(1)] = label
    return {a: l for a, l in lg.items() if len(seen[l]) == 1}


def read_headers(genome_fa):
    try:
        return sh(['grep', '-a', '^>', genome_fa]).stdout.splitlines()
    except Exception:                                   # noqa  (no grep / unreadable): stream in Python
        with open(genome_fa) as f:
            return [l for l in f if l.startswith('>')]


def build(genome_fa, genome_dir, cfg, out_tsv):
    fai = genome_fa + '.fai'
    if not os.path.exists(fai):
        sh(['samtools', 'faidx', genome_fa])
    lengths = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(fai)}
    lg = {}
    if cfg.get('chrom_to_lg'):
        lg = dict(cfg['chrom_to_lg'])
    else:
        rep = cfg.get('assembly_report') or next(iter(glob.glob(genome_dir + '*assembly_report.txt')), None)
        if rep and os.path.exists(rep):
            for line in open(rep):
                if line.startswith('#') or not line.strip():
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) >= 7 and f[1] == 'assembled-molecule' and f[6] in lengths:
                    name = f[2] if f[2] not in ('na', '') else f[0]
                    lg[f[6]] = name if name.upper().startswith('LG') else 'LG' + name
    if not lg:
        lg = lg_from_headers(read_headers(genome_fa), lengths)          # NCBI RefSeq headers carry the linkage-group label
    if not lg:
        raise SystemExit('Cannot determine chromosome -> LG mapping for this genome.\n'
                         f'  * put the NCBI *_assembly_report.txt in {genome_dir} (or set assembly_report: in config.yaml), or\n'
                         '  * set chrom_to_lg: {RefSeq accession: LGn} in config/paper.yaml (the FASTA headers had no "linkage group LGn" labels either).\n'
                         'Largest sequences in the .fai for reference:\n' +
                         '\n'.join(f'  {k}\t{v}' for k, v in sorted(lengths.items(), key=lambda x: -x[1])[:30]))
    df = pd.DataFrame([(c, l, lengths[c]) for c, l in lg.items()], columns=['chrom', 'LG', 'length'])
    df['n'] = df.LG.str.extract(r'(\d+)').astype(float)
    df = df.sort_values(['n', 'LG']).drop(columns='n').reset_index(drop=True)
    df.to_csv(out_tsv, sep='\t', index=False)
    return df
