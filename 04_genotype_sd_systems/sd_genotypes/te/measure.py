"""Transposon insertions: raw measurements from ONE coordinate-sorted, indexed BAM aligned to Mzebra_GT3_NCBI (a panel module: measure_bam(bam, cfg, run), definition_key(cfg), regions(cfg)).
Stored per site: the numbers of read pairs that support the reference allele, the insertion allele, both equally, or neither (the lab's IndelReadClassifier, vendored in longindel.py). The call is made later (call.py)."""
import hashlib, json, os
from sd_genotypes.gsdf.measure import contigs_in_header
from . import longindel as L

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
MEASUREMENT_KEYS = ('flanking', 'classifier', 'min_mapq', 'include_unmapped', 'definition_file')
MODULE = 'te'


def _path(cfg):
    return os.path.join(cfg.get('definition_dir') or DATA, cfg['definition_file'])


def load_definition(cfg):
    p = _path(cfg)
    if not os.path.exists(p):
        raise FileNotFoundError(f'{p} is missing: build the definition with `python -m sd_genotypes.te.definition --ref <GT3 .fna>` (see sd_genotypes/te/data/README.md)')
    d = json.load(open(p))
    if d['flanking'] != cfg['flanking']:
        raise ValueError(f"the definition was built with {d['flanking']} bp flanks, the config says {cfg['flanking']}")
    return d


def definition_key(cfg):
    """Hash of what is measured: the measurement keys of the config and the CONTENT of the definition file (alleles and reference flanks)."""
    h = {k: cfg[k] for k in MEASUREMENT_KEYS}
    h['file'] = hashlib.md5(open(_path(cfg), 'rb').read()).hexdigest()
    return hashlib.md5(json.dumps(h, sort_keys=True).encode()).hexdigest()[:12]


def classifier_from_site(site, flanking, params, backend):
    """The IndelReadClassifier of a site, built from the STORED flanks. Identical to IndelReadClassifier.from_vcf_record on the reference (tested)."""
    ref, alt, left, right = site['ref'].upper(), site['alt'].upper(), site['left'], site['right']
    anchor = L._common_prefix_len(ref, alt)
    if anchor == 0:
        raise ValueError(f"site {site['name']}: REF and ALT share no anchor base")
    a0 = site['pos'] - 1
    ref_bps = sorted({len(left) + anchor, len(left) + len(ref)}); alt_bps = sorted({len(left) + anchor, len(left) + len(alt)})
    return L.IndelReadClassifier(ref_allele=left + ref + right, alt_allele=left + alt + right, contig=site['chrom'], insertion_pos=a0 + anchor, variant_end_pos=a0 + len(ref),
                                 ref_breakpoints=ref_bps, alt_breakpoints=alt_bps, backend=backend, **params)


def classifiers(cfg, d=None, backend=None):
    d = d or load_definition(cfg)
    be = backend if backend is not None else L.make_backend(cfg.get('backend', 'auto'), cfg['classifier']['gap_open'], cfg['classifier']['gap_extend'])
    return [(s, classifier_from_site(s, d['flanking'], cfg['classifier'], be)) for s in d['sites']]


def merge_regions(wins):
    """(chrom, start0, end0) windows -> merged samtools regions 'chrom:start-end' (1-based, inclusive)."""
    out = []
    for chrom, a, b in sorted(wins):
        if out and out[-1][0] == chrom and a <= out[-1][2]:
            out[-1][2] = max(out[-1][2], b)
        else:
            out.append([chrom, a, b])
    return [f'{c}:{a + 1}-{b}' for c, a, b in out]


def regions(cfg):
    """The regions of the BAM this module reads (for the retrospective mode, which fetches only these)."""
    return merge_regions([(c.contig,) + c.fetch_window() for _, c in classifiers(cfg)])


def measure_bam(bam, cfg, run):
    """`run(cmd)` returns an object with .stdout (text) and raises on failure (the panel's contract). One site failing is recorded for that site (status 'error') and does not stop the others."""
    d = load_definition(cfg)
    contigs = contigs_in_header(run(['samtools', 'view', '-H', bam]).stdout)
    for c, n in d['lengths'].items():
        if contigs.get(c) != n:
            raise ValueError(f'{bam}: {c} is {contigs.get(c)} bp in the header, the project reference has {n}: not aligned to Mzebra_GT3_NCBI?')
    sites = {}; backend = None
    for s, clf in classifiers(cfg, d):
        backend = clf.backend.name
        try:
            w0, w1 = clf.fetch_window()
            lines = run(['samtools', 'view', '-F', '3328', bam, f"{clf.contig}:{w0 + 1}-{w1}"]).stdout.splitlines()
            reads = L.reads_from_samtools_lines(lines, cfg['min_mapq'], cfg['include_unmapped'])
            c = {x: 0 for x in L.Support}
            for pr in clf.classify_read_pairs(reads):
                c[pr.support] += 1
            sites[s['name']] = dict(status='OK', n_ref=c[L.Support.REF], n_alt=c[L.Support.ALT], n_equal=c[L.Support.EQUAL], n_uninf=c[L.Support.UNINFORMATIVE], n_reads=len(reads), message='')
        except Exception as e:
            sites[s['name']] = dict(status='error', n_ref=0, n_alt=0, n_equal=0, n_uninf=0, n_reads=0, message=str(e)[:160])
    return dict(definition=definition_key(cfg), backend=backend, sites=sites)
