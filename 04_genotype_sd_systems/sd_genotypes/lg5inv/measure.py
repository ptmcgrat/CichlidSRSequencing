"""LG5 inversion: raw measurements from ONE coordinate-sorted, indexed BAM aligned to Mzebra_GT3_NCBI (the shape of a panel module: measure_bam(bam, cfg, run) and definition_key(cfg)). Only samtools is needed.
Stored per sample: ref and alt read counts at every marker (samtools mpileup, MAPQ >= min_mapq, no base-quality filter) and the counts of reads containing the junction sequences. The call is made later (call.py)."""
import argparse, hashlib, json, os, tempfile
import pandas as pd
from sd_genotypes.gsdf.measure import contigs_in_header
from . import markers as K

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
MEASUREMENT_KEYS = ('region', 'junction_region', 'min_mapq', 'expected_contig_length', 'markers_file', 'junctions_file')
MODULE = 'lg5inv'


def _paths(cfg):
    d = cfg.get('definition_dir') or DATA
    return os.path.join(d, cfg['markers_file']), os.path.join(d, cfg['junctions_file'])


def load_definition(cfg):
    mp, jp = _paths(cfg)
    for p in (mp, jp):
        if not os.path.exists(p):
            raise FileNotFoundError(f'{p} is missing: build the definition with `python -m sd_genotypes.lg5inv.definition` (see sd_genotypes/lg5inv/data/README.md)')
    return pd.read_csv(mp, sep='\t'), json.load(open(jp))['probes']


def definition_key(cfg):
    """Hash of what is measured: the measurement keys of the config and the CONTENT of the marker and junction files."""
    load_definition(cfg)                                                   # a missing file is reported with the command that builds it
    mp, jp = _paths(cfg)
    h = {k: cfg[k] for k in MEASUREMENT_KEYS}
    h['files'] = [hashlib.md5(open(p, 'rb').read()).hexdigest() for p in (mp, jp)]
    return hashlib.md5(json.dumps(h, sort_keys=True).encode()).hexdigest()[:12]


def mpileup_command(bam, bed, region, mapq):
    return ['samtools', 'mpileup', '-a', '-B', '-Q', '0', '-q', str(mapq), '-d', '100000', '-l', bed, '-r', region, bam]


def junction_counts(sam_lines, probes):
    """Counts of reads (not duplicate/secondary/QC-fail/supplementary) containing probe A, B or C on either strand."""
    res = dict(A=0, B=0, C=0, reads=0)
    for line in sam_lines:
        f = line.split('\t')
        if len(f) < 10:
            continue
        res['reads'] += 1; seq = f[9].upper()
        for k in 'ABC':
            if probes[k] in seq or K.rc(probes[k]) in seq:
                res[k] += 1
    return res


def measure_bam(bam, cfg, run):
    """`run(cmd)` returns an object with .stdout (text) and raises on failure (the panel's contract)."""
    r, j = cfg['region'], cfg['junction_region']
    contigs = contigs_in_header(run(['samtools', 'view', '-H', bam]).stdout)
    if r['chrom'] not in contigs:
        raise ValueError(f"{bam}: no contig {r['chrom']} in the header: is this BAM aligned to the project reference (Mzebra_GT3_NCBI)?")
    if contigs[r['chrom']] != cfg['expected_contig_length']:
        raise ValueError(f"{bam}: {r['chrom']} is {contigs[r['chrom']]} bp, the project reference has {cfg['expected_contig_length']}: not aligned to Mzebra_GT3_NCBI")
    markers, probes = load_definition(cfg)
    with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as fh:
        for p in markers.pos:
            fh.write(f"{r['chrom']}\t{int(p) - 1}\t{int(p)}\n")
        bed = fh.name
    try:
        out = run(mpileup_command(bam, bed, f"{r['chrom']}:{r['start'] + 1}-{r['end']}", cfg['min_mapq'])).stdout
    finally:
        os.remove(bed)
    got = {}
    for line in out.splitlines():
        f = line.split('\t')
        if len(f) >= 5:
            got[int(f[1])] = f[4]
    ref_n, alt_n = [], []
    for p, rb, ab in zip(markers.pos, markers.ref, markers.alt):
        c = K.pileup_counts(got.get(int(p), ''), rb)
        ref_n.append(int(c[rb])); alt_n.append(int(c[ab]))
    jr = f"{j['chrom']}:{j['start'] + 1}-{j['end']}"
    jc = junction_counts(run(['samtools', 'view', '-F', '3840', bam, jr]).stdout.splitlines(), probes)
    return dict(definition=definition_key(cfg), contig_length=contigs[r['chrom']], region=f"{r['chrom']}:{r['start'] + 1}-{r['end']}", n_markers=len(markers), ref=ref_n, alt=alt_n, junction=jc)


def main(argv=None):
    """python -m sd_genotypes.lg5inv.measure --bam S.all.bam --sample S --out-json S.panel.json : the same JSON as `panel measure`, with only this module."""
    import yaml
    from cichlid_common.common_light import config_path
    from sd_genotypes.panel import measure as pm
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.lg5inv.measure'); ap.add_argument('--bam', required=True); ap.add_argument('--sample', required=True); ap.add_argument('--out-json', required=True)
    a = ap.parse_args(argv)
    pm.MODULES.setdefault(MODULE, ('sd_genotypes.lg5inv.measure', 'lg5panel.yaml'))                    # registered at run time: panel/measure.py is untouched
    res = pm.measure(a.bam, a.sample, modules=(MODULE,))
    pm.write_json(a.out_json, res)
    m = res['modules'][MODULE]
    print(f"[lg5inv] {a.sample}: {m['n_markers']} markers, junction a/b/std = {m['junction']['A']}/{m['junction']['B']}/{m['junction']['C']} -> {a.out_json}")


if __name__ == '__main__':
    main()
