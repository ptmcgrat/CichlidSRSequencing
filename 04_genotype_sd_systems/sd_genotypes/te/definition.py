"""Build the transposon-insertion definition from the built-in sites (sites.py) and the reference: the anchored records plus the reference flanks, so that measuring a BAM needs no FASTA.
  python -m sd_genotypes.te.definition --ref GT3.fna --out-dir DIR [--sets all] [--flanking 250]"""
import argparse, json, os
from . import sites as K

EXPECTED_FINGERPRINTS = {'TE1': 'ALT 1,140 bp, md5 1957b81b6617', 'TE2': 'ALT 9,869 bp, md5 d2b75a3b30e2', 'chr7ins_c1': 'ALT 603 bp, md5 8b7b5cf1d8c7', 'chr7ins_c2': 'ALT 603 bp, md5 8b7b5cf1d8c7'}


def read_lengths_fai(fasta):
    out = {}
    for line in open(fasta + '.fai'):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 2:
            out[f[0]] = int(f[1])
    return out


def build(fetch, lengths, out_dir, sets='all', flanking=250):
    """fetch(chrom, start0, end0) -> reference sequence (0-based half-open). Returns the path of te_definition.json."""
    names = list(K.BUILTIN) if sets == 'all' else [sets] if isinstance(sets, str) else list(sets)
    sites = []
    for s in names:
        for chrom, pos, name, ref, alt, info in K.BUILTIN[s]:
            if chrom not in lengths:
                raise ValueError(f'{chrom} (site {name}) is not in the reference lengths')
            ref, alt = ref.upper(), alt.upper()
            a0 = pos - 1
            found = fetch(chrom, a0, a0 + len(ref)).upper()
            if found != ref:
                raise ValueError(f'site {name}: the anchor/REF {ref!r} does not match the reference at {chrom}:{pos} (found {found!r}): wrong reference, or an off-by-one?')
            left, right = fetch(chrom, a0 - flanking, a0).upper(), fetch(chrom, a0 + len(ref), a0 + len(ref) + flanking).upper()
            if len(left) != flanking or len(right) != flanking:
                raise ValueError(f'site {name}: could not read {flanking} bp of flank on both sides (contig edge?)')
            fp = K.fingerprint(alt)
            if name in EXPECTED_FINGERPRINTS and fp != EXPECTED_FINGERPRINTS[name]:
                raise ValueError(f'site {name}: allele fingerprint {fp} differs from the pinned {EXPECTED_FINGERPRINTS[name]}')
            sites.append(dict(name=name, set=s, chrom=chrom, pos=pos, ref=ref, alt=alt, left=left, right=right, fingerprint=fp, info=info))
    have = {x['name'] for x in sites}
    groups = {g: list(m) for g, m in K.SITE_GROUPS.items() if all(x in have for x in m)}
    d = dict(flanking=flanking, sets=names, lengths={c: lengths[c] for c in sorted({x['chrom'] for x in sites})}, sites=sites, groups=groups,
             group_notes={g: K.SITE_GROUP_NOTES.get(g, '') for g in groups}, source='sexloci 0.1.26 known_sites.py')
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, 'te_definition.json')
    json.dump(d, open(path, 'w'), indent=1)
    return path


def main(argv=None):
    import hashlib, subprocess
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.te.definition'); ap.add_argument('--ref', required=True); ap.add_argument('--out-dir', required=True)
    ap.add_argument('--sets', default='all'); ap.add_argument('--flanking', type=int, default=250)
    a = ap.parse_args(argv)
    if not os.path.exists(a.ref + '.fai'):
        subprocess.run(['samtools', 'faidx', a.ref], check=True)

    def fetch(chrom, s0, e0):
        out = subprocess.run(['samtools', 'faidx', a.ref, f'{chrom}:{max(s0, 0) + 1}-{e0}'], capture_output=True, text=True, check=True).stdout.split('\n')[1:]
        return ''.join(out)
    p = build(fetch, read_lengths_fai(a.ref), a.out_dir, a.sets if a.sets == 'all' else a.sets.split(','), a.flanking)
    print(f'wrote {p} (md5 {hashlib.md5(open(p, "rb").read()).hexdigest()[:12]}); sites: ' + ', '.join(s['name'] for s in json.load(open(p))['sites']))


if __name__ == '__main__':
    main()
