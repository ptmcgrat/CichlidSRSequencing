"""Build the LG5 measurement definition: the marker SNV list and the junction probe sequences (from the reference).
  python -m sd_genotypes.lg5inv.definition --ref GT3.fna --markers markers.tsv --out-dir sd_genotypes/lg5inv/data"""
import argparse, hashlib, json, os
import pandas as pd
from . import markers as K

# From the OB-pool reads clipped at the right breakpoint (junction a: reads starting at 26,726,677; junction b: reads ending at 26,726,685); both flank the same 15 bp of GT3 at 19,657,763-777.
CLIP_L15, CLIP_R15 = 'GAAATCAGATCGTCC', 'GGACGATCTGATTTC'


def junction_probes(fetch):
    """fetch(a, b) -> GT3 NC_135171.1 a..b (1-based, inclusive) uppercase. Returns the 30-mers: A = RC(19,657,763-777)|26,726,677-691 (inverted), B = ..26,726,685|19,657,763-777 (inverted), C = 26,726,671-700 (standard)."""
    if fetch(19657763, 19657777) != CLIP_R15 or K.rc(CLIP_R15) != CLIP_L15:
        raise ValueError('the reference at NC_135171.1:19,657,763-777 does not match the junction reads: wrong reference?')
    ra = fetch(26726677, 26726706); rb = fetch(26726656, 26726685)
    return dict(A=CLIP_L15 + ra[:15], B=rb[-15:] + CLIP_R15, C=rb[-15:] + ra[9:24])


def build(fetch, markers, out_dir):
    m = markers[markers.type == 'SNV'] if 'type' in markers.columns else markers
    m = m[['pos', 'ref', 'alt']].drop_duplicates().sort_values('pos').reset_index(drop=True)
    if m.pos.duplicated().any():
        raise ValueError('two different alleles at one position in the marker list')
    for p, r in zip(m.pos.head(20), m.ref.head(20)):
        if fetch(int(p), int(p)) != r:
            raise ValueError(f'marker ref {r} at {K.CHROM}:{p} does not match the reference: wrong reference or marker file')
    m.insert(0, 'chrom', K.CHROM)
    os.makedirs(out_dir, exist_ok=True)
    mp, jp = os.path.join(out_dir, 'lg5_markers.tsv'), os.path.join(out_dir, 'lg5_junctions.json')
    m.to_csv(mp, sep='\t', index=False)
    json.dump(dict(probes=junction_probes(fetch), note='A, B: inverted junctions; C: standard sequence across the right breakpoint (GT3 26,726,671-700)'), open(jp, 'w'), indent=1)
    return mp, jp, len(m)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.lg5inv.definition'); ap.add_argument('--ref', required=True); ap.add_argument('--markers', required=True); ap.add_argument('--out-dir', required=True)
    a = ap.parse_args(argv)
    mp, jp, n = build(lambda x, y: K.faidx(a.ref, x, y), pd.read_csv(a.markers, sep='\t'), a.out_dir)
    md5 = lambda p: hashlib.md5(open(p, 'rb').read()).hexdigest()[:12]
    print(f'{n} markers -> {mp} (md5 {md5(mp)}); junction probes -> {jp} (md5 {md5(jp)})')


if __name__ == '__main__':
    main()
