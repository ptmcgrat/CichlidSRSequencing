"""Synthetic tests for the genome registry helpers, the inversion-breakpoint pipeline and the centromere module (called from selftest)."""
import json, os, tempfile, types
import numpy as np
import pandas as pd

M = 1_000_000


def _ok(name, cond, extra=''):
    print(('PASS ' if cond else 'FAIL ') + name, extra)
    if not cond:
        raise SystemExit(1)


def line(q, qlen, qs, qe, st, t, tlen, ts, te, cg=None):
    qlen, qs, qe, tlen, ts, te = [int(x) for x in (qlen, qs, qe, tlen, ts, te)]
    n = min(qe - qs, te - ts); f = [q, qlen, qs, qe, st, t, tlen, ts, te, n, max(qe - qs, te - ts), 60]
    return '\t'.join(map(str, f)) + ('\tcg:Z:' + cg if cg else '')


def run():
    from .analyses.inversion_breakpoints import paf, steps
    from .analyses.centromeres import cen, steps as csteps
    from . import genomes_registry as gr
    from .runner import Run
    # ---------------- genome registry helpers
    _ok('guess_accession', gr.guess_accession('ASM4766397v1') == 'GCA_047663971.1' and gr.guess_accession('weird') is None)
    es = {'result': {'uids': ['11'], '11': {'assemblyaccession': 'GCA_047663971.1', 'assemblyname': 'ASM4766397v1', 'speciesname': 'Labeotropheus trewavasae',
                                            'ftppath_genbank': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/047/663/971/GCA_047663971.1_ASM4766397v1'}}}
    rec = gr.parse_esummary(es)[0]
    gz, md5u, fn = gr.download_urls(rec['ftp_path'])
    _ok('NCBI summary -> download URLs', gz.startswith('https://') and gz.endswith('GCA_047663971.1_ASM4766397v1_genomic.fna.gz') and md5u.endswith('md5checksums.txt'))
    _ok('parse esearch / md5checksums', gr.parse_esearch({'esearchresult': {'idlist': ['1', '2']}}) == ['1', '2'] and gr.parse_md5checksums('abc  ./x_genomic.fna.gz\ndef  ./README.txt\n')['x_genomic.fna.gz'] == 'abc')
    d = tempfile.mkdtemp(); fa = d + '/t.fa'; open(fa, 'w').write('>a\nACGTACGT\nAC\n>b\nGGGG\n')
    _ok('fasta_stats', gr.fasta_stats(fa) == (2, 14, 10))

    # ---------------- call step on synthetic alignments: LG5 (two Behrens assemblies) and LG2 (YH, chromosome-level contig with an inversion)
    T5, T2, ref = 'NC_T5', 'NC_T2', 'Mzebra_GT3_NCBI'
    L5, L2 = 35 * M, 40 * M
    adir = tempfile.mkdtemp() + '/'
    latr = [line('ctgA', 13.2 * M, 0, 2.45 * M, '+', T5, L5, 12.0 * M, 14.45 * M), line('ctgA', 13.2 * M, 2.55 * M, 4.5 * M, '-', T5, L5, 19.45 * M, 21.4 * M),
            line('ctgA', 13.2 * M, 4.9 * M, 7.0 * M, '-', T5, L5, 17.35 * M, 19.45 * M), line('ctgA', 13.2 * M, 7.6 * M, 9.45 * M, '-', T5, L5, 15.5 * M, 17.35 * M),
            line('ctgA', 13.2 * M, 9.6 * M, 10.6 * M, '-', T5, L5, 14.55 * M, 15.5 * M), line('ctgA', 13.2 * M, 10.65 * M, 13.2 * M, '+', T5, L5, 21.45 * M, 24.0 * M),
            line('ctgB', 12 * M, 0, 12 * M, '+', T5, L5, 12.0 * M, 24.0 * M)]
    meze = [line('ctgC', 6.9 * M, 0, 6.9 * M, '-', T5, L5, 14.55 * M, 21.4 * M)]
    yh = [line('YH_c2', 41.5 * M, 0, 9.75 * M, '+', T2, L2, 0, 9.75 * M), line('YH_c2', 41.5 * M, 9.85 * M, 15 * M, '-', T2, L2, 27.3 * M, 32.5 * M),
          line('YH_c2', 41.5 * M, 15.6 * M, 24 * M, '-', T2, L2, 18.9 * M, 27.3 * M), line('YH_c2', 41.5 * M, 24.9 * M, 34.8 * M, '-', T2, L2, 9.8 * M, 18.9 * M),
          line('YH_c2', 41.5 * M, 35.0 * M, 41.5 * M, '+', T2, L2, 32.55 * M, 38.0 * M)]
    for name, ls in (('LatrZW', latr), ('MezeOBm', meze), ('YH', yh)):
        open(f'{adir}{name}_vs_{ref}.targets.paf', 'w').write('\n'.join(ls) + '\n')
    cfg = steps.cfg_of(); cfg['inversions'] = {'LG5': {'assemblies': ['LatrZW', 'MezeOBm'], 'reference_is_standard': None, 'label': 'LG5'},
                                               'LG2': {'assemblies': ['YH'], 'reference_is_standard': True, 'label': 'LG2'}}
    steps.cfg_of = lambda: cfg
    loc = tempfile.mkdtemp() + '/'
    proj = types.SimpleNamespace(fm=types.SimpleNamespace(localMasterDir=loc, cloudMasterDir=loc + 'cloud/'), root=loc + 'P/G/', genome=ref, cfg={'runner': {'upload': False}},
                                 chrom_for_lg=lambda lg: {'LG5': T5, 'LG2': T2}[lg], chrom_lengths=lambda: {T5: L5, T2: L2}, lg={T5: 'LG5', T2: 'LG2'}, cores=2, run=None)
    proj.run = lambda a, s, params=None, tar=(): Run(proj, a, s, params, tar)
    os.makedirs(proj.root)
    _real_excl = steps.excluded_regions
    steps.excluded_regions = lambda p_: {}                      # this synthetic data is labelled LG2 but is not the real LG2: do not apply the real exclusion here
    steps.call_step(proj, align_dir=adir)
    steps.excluded_regions = _real_excl
    def last(a, s):
        base = f'{proj.root}{a}/{s}/'; return base + sorted(x for x in os.listdir(base) if x != 'LATEST')[-1] + '/'
    rd = last('inversion_breakpoints', 'call')
    c5 = pd.read_csv(rd + 'LG5_consensus_breakpoints.tsv', sep='\t'); c2 = pd.read_csv(rd + 'LG2_consensus_breakpoints.tsv', sep='\t')
    _ok('LG5: consensus from LatrZW (A+B) and MezeOBm (C)', len(c5) == 1 and abs(c5.left_breakpoint.iloc[0] - 14.55 * M) <= 60_000 and abs(c5.right_breakpoint.iloc[0] - 21.4 * M) <= 60_000 and c5.n_calls.iloc[0] == 3)
    _ok('LG2: YH chromosome-level contig gives the inversion (method A)', len(c2) == 1 and abs(c2.left_breakpoint.iloc[0] - 9.8 * M) <= 100_000 and abs(c2.right_breakpoint.iloc[0] - 32.5 * M) <= 100_000
        and 'YH:A' in c2.sources.iloc[0], str(c2[['left_breakpoint', 'right_breakpoint', 'sources']].values.tolist()))
    p2 = json.load(open(rd + 'LG2_polarity.json')); p5 = json.load(open(rd + 'LG5_polarity.json'))
    _ok('polarity text follows reference_is_standard (true / unknown)', 'INVERTED' in p2['polarity'] and 'not established' in p5['polarity'])
    qcols = pd.read_csv(rd + 'LG2_breakpoints_per_assembly.tsv', sep='\t')
    _ok('query coordinates of the inverted segment are reported', {'inv_q_start', 'inv_q_end'} <= set(qcols.columns) and qcols.inv_q_start.iloc[0] == int(9.85 * M))
    _ok('BED and dotplots written for every chromosome', all(os.path.exists(rd + f) for f in ('LG5_breakpoint_regions.bed', 'LG2_breakpoint_flanks_50kb.bed', 'LG5_LatrZW_dotplot.png', 'LG2_YH_dotplot.svg')))
    _ok('no UMD2a anchor outputs remain', not any('orientation' in f or 'UMD2a' in f for f in os.listdir(rd)))

    # ---------------- centromere module
    rm = ('   SW  perc perc perc  query      position in query           matching  repeat\n score  div. del. ins.  sequence   begin     end    (left)    repeat    class/family\n\n'
          '  1000  5.0  0.0  0.0  chrX  100001  200000  (1000)  +  SAT1  Satellite/centr  1 2 3 1\n  500  8.0  0.0  0.0  chrX  300001  300500  (900)  C  L1  LINE/L1  1 2 3 2\n')
    gff = '##gff-version 3\nchrX\tRepeatMasker\trepeat_region\t100001\t200000\t.\t+\t.\tID=a;Target=Motif:SAT1 1 100;class_family=Satellite/centr\nchrX\tRepeatMasker\trepeat_region\t300001\t300500\t.\t+\t.\tID=b;class_family=LINE/L1\n'
    bed = 'chrX\t100000\t200000\tSAT1\nchrX\t300000\t300500\tSAT2\n'
    ucsc = '#bin\tswScore\tmilliDiv\tmilliDel\tmilliIns\tgenoName\tgenoStart\tgenoEnd\tgenoLeft\tstrand\trepName\trepClass\trepFamily\trepStart\trepEnd\trepLeft\tid\n585\t1000\t50\t0\t0\tchrX\t100000\t200000\t-1\t+\tSAT1\tSatellite\tcentr\t1\t100\t0\t1\n585\t500\t80\t0\t0\tchrX\t300000\t300500\t-1\t+\tL1\tLINE\tL1\t1\t100\t0\t2\n'
    ok = True
    for ext, txt in (('.out', rm), ('.gff3', gff), ('.bed', bed), ('.rmsk.txt', ucsc)):
        f = tempfile.mkdtemp() + '/ann' + ext; open(f, 'w').write(txt)
        df = cen.parse_annotation(f, ['Satellite'])
        ok &= (ext == '.bed' and len(df) == 2) or (ext != '.bed' and len(df) == 1 and df.start.iloc[0] == 100000 and df.end.iloc[0] == 200000)
    _ok('annotation parsing: RepeatMasker .out, GFF3, UCSC rmsk and BED (class filter applied)', ok)
    cov = cen.window_coverage([100, 150, 5000], [400, 600, 5100], 10_000, 1000)
    _ok('window_coverage merges overlaps and splits across windows', abs(cov[0] - 0.5) < 1e-9 and abs(cov[5] - 0.1) < 1e-9 and cov[1:5].sum() == 0)
    d2 = tempfile.mkdtemp(); f2 = d2 + '/g.fa'; open(f2, 'w').write('>c1 desc\nACGTNN\nNNACGT\n>c2\nACGT\nnnnn\n')
    g = cen.fasta_lengths_and_gaps(f2)
    _ok('FASTA N-gap runs (across lines, lower case)', g['c1'] == (12, [(4, 8)]) and g['c2'] == (8, [(4, 8)]))
    conf = {'high': {'min_ratio': 5.0, 'min_dominance': 2.0, 'min_peak': 0.10}, 'medium': {'min_ratio': 3.0, 'min_dominance': 1.5, 'min_peak': 0.05}}
    rng = np.random.default_rng(2)
    base = rng.uniform(0, 0.02, 800)                                # 40 Mb at 50 kb windows
    s1 = base.copy(); s1[360:420] = 0.7                             # satellite array at 18-21 Mb
    r1 = cen.call_centromere(s1, 50_000, 1_000_000, conf)
    two = base.copy(); two[100:160] = 0.6; two[500:560] = 0.55      # two equal peaks: ambiguous
    r2 = cen.call_centromere(two, 50_000, 1_000_000, conf)
    r3 = cen.call_centromere(np.zeros(800), 50_000, 1_000_000, conf)
    _ok('centromere call: dominant array -> high, located; two equal peaks -> low; no satellites -> none',
        r1['confidence'] == 'high' and 17.5 * M <= r1['start'] and r1['end'] <= 21.5 * M and r2['confidence'] == 'low' and r3['confidence'] == 'none', f"{r1['start']}-{r1['end']} ratio {r1['ratio']:.0f}")
    # scan_genome on a small synthetic genome (annotation in RepeatMasker .out format; 3 chromosomes of 6 Mb)
    fa = d2 + '/s.fa'
    with open(fa, 'w') as fh:
        for c in ('chrA', 'chrB', 'tiny'):
            n = 6_000_000 if c != 'tiny' else 100_000
            fh.write(f'>{c}\n'); fh.write(('ACGT' * 25 + '\n') * (n // 100))
    ann = d2 + '/s.out'
    with open(ann, 'w') as fh:
        fh.write('h1\nh2\n\n')
        fh.write('  900  5.0  0.0  0.0  chrA  2000001  3500000  (1)  +  CEN1  Satellite/centr  1 2 3 1\n')
        for k in range(0, 5_000_000, 700_000):
            fh.write(f'  100  5.0  0.0  0.0  chrB  {k + 1}  {k + 3000}  (1)  +  SAT9  Satellite  1 2 3 2\n')
    ccfg = csteps.cfg_of(); ccfg.update(min_chrom_length=1_000_000)
    calls, prof = csteps.scan_genome('toy', fa, ann, ccfg)
    ca = calls[calls.sequence == 'chrA'].iloc[0]
    _ok('scan_genome: chromosome-scale sequences only; chrA centromere found with its family, chrB (scattered) not called high',
        set(calls.sequence) == {'chrA', 'chrB'} and ca.confidence in ('high', 'medium') and 1.8 * M <= ca.start and ca.end <= 3.8 * M and 'CEN1' in ca.top_families
        and calls[calls.sequence == 'chrB'].iloc[0].confidence != 'high', f"{ca.start}-{ca.end} {ca.confidence} {ca.top_families}")
    # lift the called interval through the assembly alignment when the satellite array itself does not align
    yh2 = paf.parse_paf([line('YH_c2', 41.5 * M, 0, 18.0 * M, '+', T2, L2, 0, 18.0 * M, f'{int(18.0 * M)}M'), line('YH_c2', 41.5 * M, 21.0 * M, 41.5 * M, '+', T2, L2, 21.0 * M, 40.0 * M, f'{int(19.0 * M)}M')])
    lf = paf.lift_interval(yh2, 'YH_c2', int(18.5 * M), int(20.5 * M))
    _ok('lift_interval bounds a non-alignable centromere by its alignable flanks', lf['status'] == 'ok' and abs(lf['tstart'] - 18.0 * M) <= 25_000 and abs(lf['tend'] - 21.0 * M) <= 25_000, str(lf))
    # ---------------- overview step reads the call run and a centromere run
    with Run(proj, 'centromeres', 'scan') as r:
        pd.DataFrame([dict(genome=ref, sequence=T5, lg='LG5', start=3 * M, end=4 * M, ref_start=3 * M, ref_end=4 * M, confidence='high', lift_status='reference'),
                      dict(genome='YH', sequence='YH_c2', lg='LG2', start=18 * M, end=20 * M, ref_start=18 * M, ref_end=21 * M, confidence='medium', lift_status='ok')]).to_csv(r.dir + 'centromeres.tsv', sep='\t', index=False)
    steps.figure_step(proj)
    ov = pd.read_csv(last('inversion_breakpoints', 'figure') + 'inversion_overview.tsv', sep='\t').set_index('lg')
    _ok('overview: inversion and centromere per chromosome (centromere lifted from YH where the reference has none)', ov.loc['LG5', 'cen_start'] == 3 * M and ov.loc['LG2', 'cen_end'] == 21 * M
        and abs(ov.loc['LG2', 'inv_start'] - 9.8 * M) <= 100_000 and os.path.exists(last('inversion_breakpoints', 'figure') + 'inversion_overview.svg'))

    # ---------------- annotation registry: names, md5sum.txt, registration, RepeatMasker command
    d3 = tempfile.mkdtemp() + '/'
    open(d3 + 'g.fna', 'w').write('>x\nACGT\n'); open(d3 + 'md5sum.txt', 'w').write(f"{gr.md5_of(d3 + 'g.fna')}  ./g.fna\n")
    _ok('md5sum.txt next to a genome file is checked', gr.verify_against_md5sum_file(d3 + 'g.fna') is True)
    open(d3 + 'md5sum.txt', 'w').write('0' * 32 + '  ./g.fna\n')
    try:
        gr.verify_against_md5sum_file(d3 + 'g.fna'); bad = False
    except RuntimeError:
        bad = True
    _ok('a wrong md5 in md5sum.txt stops the registration', bad)
    out_txt = ('   SW  perc perc perc  query  position in query  matching  repeat\n score  div. del. ins.  sequence begin end (left) repeat class/family\n\n'
               '  900  5.0  0.0  0.0  LG5  100001  200000  (1)  +  SAT1  Satellite  1 2 3 1\n  900  5.0  0.0  0.0  LG5  300001  300500  (1)  +  L1  LINE/L1  1 2 3 2\n  50  5.0  0.0  0.0  chr9  1  500  (1)  +  X  DNA  1 2 3 3\n')
    nm = gr.annotation_names(d3 + 'o.out') if open(d3 + 'o.out', 'w').write(out_txt) else None
    _ok('annotation sequence names read from a RepeatMasker .out', nm == {'LG5': 2, 'chr9': 1})
    v = gr.validate_annotation(nm, {'NC_5': 35 * M, 'NC_9': 36 * M}, {'LG5': 'NC_5', 'LG9': 'NC_9'})
    _ok('names are mapped through LG labels (LG5 -> NC_5, chr9 -> NC_9)', v['n_unmapped'] == 0 and v['mapping'] == {'LG5': 'NC_5', 'chr9': 'NC_9'} and not v['identity'])
    v2 = gr.validate_annotation({'scaffold_77': 5}, {'NC_5': 35 * M}, {'LG5': 'NC_5'})
    _ok('names that cannot be mapped are reported', v2['n_unmapped'] == 1 and v2['unmapped_names'] == ['scaffold_77'])
    # register_annotation with a fake project
    gdir = tempfile.mkdtemp() + '/'; src = tempfile.mkdtemp() + '/'
    os.makedirs(gdir + 'GT3/'); open(gdir + 'GT3/g.fa', 'w').write('>NC_5\nACGT\n>NC_9\nACGT\n'); open(gdir + 'GT3/g.fa.fai', 'w').write(f'NC_5\t{35 * M}\t0\t0\t0\nNC_9\t{36 * M}\t0\t0\t0\n')
    json.dump({'fasta': gdir + 'GT3/g.fa'}, open(gdir + 'GT3/manifest.json', 'w'))
    for fn in ('a.out', 'a.tbl'):
        open(src + fn, 'w').write(out_txt if fn.endswith('.out') else 'tbl')
    fp = types.SimpleNamespace(genomes_dir=gdir, genome='GT3', lg={'NC_5': 'LG5', 'NC_9': 'LG9'}, cfg={'runner': {'upload': False}})
    man = gr.register_annotation(fp, 'GT3', dict(source_dir=src, files=['a.out', 'a.tbl'], main='a.out'), upload=False)
    apath, nmap = gr.annotation_for(fp, 'GT3')
    _ok('register_annotation copies, md5s, validates and exposes the name map', man['usable'] and os.path.exists(gdir + 'GT3/annotation/a.tbl') and len(man['files']['a.out']['md5']) == 32
        and apath.endswith('annotation/a.out') and nmap['LG5'] == 'NC_5')
    open(src + 'b.out', 'w').write(out_txt.replace('LG5', 'weird1'))
    gr.register_annotation(fp, 'GT3', dict(source_dir=src, files=['b.out'], main='b.out'), upload=False)
    try:
        gr.annotation_for(fp, 'GT3'); unusable = False
    except FileNotFoundError:
        unusable = True
    _ok('an annotation whose names do not match the genome is refused', unusable)
    ann2 = d2 + '/lg.out'
    open(ann2, 'w').write('h1\nh2\n\n  900  5.0  0.0  0.0  LGA  2000001  3500000  (1)  +  CEN1  Satellite/centr  1 2 3 1\n')
    calls2, _ = csteps.scan_genome('toy2', fa, ann2, ccfg, name_map={'LGA': 'chrA'})
    _ok('centromere scan applies the annotation -> genome name map', calls2[calls2.sequence == 'chrA'].iloc[0].confidence in ('high', 'medium'))

    # ---------------- strand-independent inversion calling (reverse-complemented contigs, inversion covering half the chromosome)
    def flip(l):
        f = l.split('\t'); qlen, qs, qe = int(f[1]), int(f[2]), int(f[3]); f[2], f[3] = str(qlen - qe), str(qlen - qs); f[4] = '-' if f[4] == '+' else '+'; return '\t'.join(f)
    TC, LC = 'NC_T10', 36 * M
    c10 = [line('c10', 36.2 * M, 0, 11.35 * M, '+', TC, LC, 0, 11.35 * M), line('c10', 36.2 * M, 11.45 * M, 22.5 * M, '-', TC, LC, 18.2 * M, 29.3 * M),
           line('c10', 36.2 * M, 22.9 * M, 29.7 * M, '-', TC, LC, 11.4 * M, 18.2 * M), line('c10', 36.2 * M, 29.8 * M, 36.2 * M, '+', TC, LC, 29.4 * M, 36.0 * M)]
    ra = paf.call_method_AC(paf.parse_paf(c10), TC); rb = paf.call_method_AC(paf.parse_paf([flip(x) for x in c10]), TC)
    same = (len(ra) == len(rb) == 1 and ra.method.iloc[0] == rb.method.iloc[0] == 'A' and ra.left_lo.iloc[0] == rb.left_lo.iloc[0] and ra.right_hi.iloc[0] == rb.right_hi.iloc[0])
    _ok('inversion call is identical for a contig and its reverse complement (even with the inversion at 51% of the chromosome)', same and ra.inv_strand.iloc[0] != rb.inv_strand.iloc[0] and not ra.strand_dependent.iloc[0])
    half = [line('h', 20 * M, 0, 5 * M, '+', TC, LC, 5 * M, 10 * M), line('h', 20 * M, 5.1 * M, 20 * M, '-', TC, LC, 10.1 * M, 25 * M)]
    rh = paf.call_method_AC(paf.parse_paf(half), TC)
    _ok('a contig that ends at a breakpoint gives a half call, flagged strand-dependent', len(rh) == 1 and rh.method.iloc[0] == 'A-half' and bool(rh.strand_dependent.iloc[0]))
    only_plus = paf.call_method_AC(paf.parse_paf([line('p', 8 * M, 0, 7 * M, '+', TC, LC, 10 * M, 17 * M)]), TC)
    _ok('a single forward-strand contig is never called an inversion', len(only_plus) == 0)

    # ---------------- annotation: GT3a-like names (LG*, scaffold*, mitochondrion) against NCBI names, lengths from the masked FASTA
    fai = {f'NC_1351{67 + i}.1': (30 + i) * M for i in range(3)}; fai.update({'NW_027490031.1': 3_000_001, 'NW_027490032.1': 2_000_005, 'NW_027490033.1': 2_000_005})
    names = {'LG1': 500, 'LG2': 500, 'LG3': 500, 'scaffold01': 20, 'scaffold02': 10, 'scaffold03': 10, 'Mzebra_UMD2a_mitochrondrial_genome': 5}
    lens = {'LG1': 30 * M, 'LG2': 31 * M, 'LG3': 32 * M, 'scaffold01': 3_000_001, 'scaffold02': 2_000_005, 'scaffold03': 2_000_005, 'Mzebra_UMD2a_mitochrondrial_genome': 16_900}
    v = gr.validate_annotation(names, fai, {'LG1': 'NC_135167.1', 'LG2': 'NC_135168.1', 'LG3': 'NC_135169.1'}, lens)
    _ok('GT3a-style names: LGs mapped through LG labels, unique-length scaffold mapped, duplicate-length scaffolds and the mito reported but not blocking',
        v['usable'] and v['mapping']['LG1'] == 'NC_135167.1' and v['mapping']['scaffold01'] == 'NW_027490031.1' and 'scaffold02' not in v['mapping'] and
        'Mzebra_UMD2a_mitochrondrial_genome' in v['unmapped_names'] and not v['length_mismatches'], f"unmapped {v['unmapped_names']}")
    bad = dict(lens, LG2=31 * M + 7)
    v2 = gr.validate_annotation(names, fai, {'LG1': 'NC_135167.1', 'LG2': 'NC_135168.1', 'LG3': 'NC_135169.1'}, bad)
    _ok('a chromosome whose length differs between annotation and genome makes the annotation unusable', (not v2['usable']) and v2['length_mismatches'] == ['LG2'])
    v3 = gr.validate_annotation({'LG1': 5, 'LG9': 5}, fai, {'LG1': 'NC_135167.1'}, {'LG1': 30 * M, 'LG9': 36 * M})
    _ok('a missing chromosome (LG9 absent from the genome) makes it unusable', not v3['usable'])
    tbl = tempfile.mkdtemp() + '/x.tbl'
    open(tbl, 'w').write('file name: x.fasta\nsequences:            64\ntotal length:  962191965 bp  (950936190 bp excl N/X-runs)\n\nRepeatMasker version 4.1.2-p1 , default mode\n\nrun with rmblastn version 2.14.1+\nThe query was compared to classified sequences in ".../cichlid_merged_lib.fa"\nFamDB:\n')
    pv = gr.parse_tbl(tbl)
    _ok('.tbl provenance: RepeatMasker, rmblastn, library name, sequences, length', pv['repeatmasker_version'] == '4.1.2-p1' and pv['rmblastn_version'] == '2.14.1+' and pv['library_name'] == 'cichlid_merged_lib.fa'
        and pv['n_sequences'] == 64 and pv['total_length'] == 962191965, str(pv))
    # a second genome annotated with a different RepeatMasker is flagged against the reference
    fp2 = types.SimpleNamespace(genomes_dir=gdir, genome='GT3', lg={'NC_5': 'LG5', 'NC_9': 'LG9'}, cfg={'runner': {'upload': False}})
    os.makedirs(gdir + 'YHX/'); open(gdir + 'YHX/y.fa', 'w').write('>NC_5\nACGT\n'); open(gdir + 'YHX/y.fa.fai', 'w').write(f'NC_5\t{35 * M}\t0\t0\t0\n'); json.dump({'fasta': gdir + 'YHX/y.fa'}, open(gdir + 'YHX/manifest.json', 'w'))
    open(src + 'ref.tbl', 'w').write(open(tbl).read()); open(src + 'y.tbl', 'w').write(open(tbl).read().replace('4.1.2-p1', '4.1.5').replace('cichlid_merged_lib', 'other_lib'))
    open(src + 'y.out', 'w').write('  900  5.0  0.0  0.0  NC_5  100001  200000  (1)  +  SAT1  Satellite  1 2 3 1\n')
    gr.register_annotation(fp2, 'GT3', dict(source_dir=src, files=['a.out', 'ref.tbl'], main='a.out'), upload=False)
    mm = gr.register_annotation(fp2, 'YHX', dict(source_dir=src, files=['y.out', 'y.tbl'], main='y.out', compare_to='GT3'), upload=False)
    _ok('a different RepeatMasker version or library than the reference annotation raises a comparability warning', len(mm['comparability_warnings']) == 2 and mm['usable'], str(mm['comparability_warnings']))
    # FileManager-registered genome by path under the master directory, with an expected md5
    mdir = tempfile.mkdtemp() + '/'; os.makedirs(mdir + '2026_X/Chromosome10s/'); gf = mdir + '2026_X/Chromosome10s/YH_Y.fasta'; open(gf, 'w').write('>h2tg000014l\nACGT\n'); open(gf + '.fai', 'w').write('h2tg000014l\t4\t0\t4\t5\n')
    fp3 = types.SimpleNamespace(genomes_dir=tempfile.mkdtemp() + '/', fm=types.SimpleNamespace(localMasterDir=mdir, localGenomesDir=mdir + 'Genomes/', downloadData=lambda x: None), cfg={'runner': {'upload': False}})
    m3 = gr.register_filemanager(fp3, 'YH_Y_hifiasm', dict(source='filemanager', master_path='2026_X/Chromosome10s/YH_Y.fasta', md5=gr.md5_of(gf)))
    try:
        gr.register_filemanager(fp3, 'YH_Y_bad', dict(source='filemanager', master_path='2026_X/Chromosome10s/YH_Y.fasta', md5='0' * 32)); wrong = False
    except RuntimeError:
        wrong = True
    _ok('FileManager genome registered by path with md5 verified; a wrong expected md5 is refused', m3['md5_verified'] is True and m3['n_contigs'] == 1 and wrong)

    # ---------------- LG labels read from the FASTA headers of the GT3 RefSeq assembly (headers as pasted from utaka)
    from . import genome as gm
    hdr = [f'>NC_{135167 + i}.1 Maylandia zebra isolate NMK-2024a linkage group LG{lgn}, Mzebra_GT3a, whole genome shotgun sequence' for i, lgn in enumerate(list(range(1, 21)) + [22, 23])]
    hdr += ['>NW_027490031.1 Maylandia zebra isolate NMK-2024a unplaced genomic scaffold, Mzebra_GT3a, whole genome shotgun sequence', '>NW_027490032.1 Maylandia zebra isolate NMK-2024a unplaced genomic scaffold, Mzebra_GT3a']
    lens_gt3 = {f'NC_{135167 + i}.1': 30 * M for i in range(22)}; lens_gt3.update({'NW_027490031.1': 2_700_000, 'NW_027490032.1': 4_400_000})
    lgm = gm.lg_from_headers(hdr, lens_gt3)
    _ok('LG labels read from RefSeq FASTA headers (22 chromosomes, no LG21, scaffolds skipped)', len(lgm) == 22 and lgm['NC_135176.1'] == 'LG10' and lgm['NC_135187.1'] == 'LG22' and lgm['NC_135188.1'] == 'LG23'
        and 'LG21' not in lgm.values() and not any(k.startswith('NW_') for k in lgm))
    _ok('headers whose sequences are not in the .fai, or that repeat a label, are ignored', gm.lg_from_headers(hdr + ['>NC_999.1 x linkage group LG5, y'], lens_gt3).get('NC_999.1') is None and
        gm.lg_from_headers(['>a x linkage group LG1, y', '>b x linkage group LG1, y'], {'a': 1, 'b': 1}) == {})
    ftmp = tempfile.mkdtemp() + '/'; fa5 = ftmp + 'g.fna'
    with open(fa5, 'w') as fh:
        for h in hdr[:3]: fh.write(h + '\nACGT\n')
    open(fa5 + '.fai', 'w').write(''.join(f'{k}\t4\t0\t4\t5\n' for k in list(lens_gt3)[:3]))
    tab = gm.build(fa5, ftmp, {}, ftmp + 'chromosomes.tsv')
    _ok('chromosome table builds from the headers when there is no assembly report or config entry', list(tab.LG) == ['LG1', 'LG2', 'LG3'])

    # ---------------- the project folder is the archive: genomes are real, verified copies; links from earlier versions are converted; small run inputs are kept
    mf3 = json.load(open(gr.entry_dir(fp3, 'YH_Y_hifiasm') + 'manifest.json'))
    cp = mf3['fasta']
    _ok('FileManager genome is COPIED into the project folder (real file, same md5, index present), not linked',
        cp.startswith(fp3.genomes_dir) and not os.path.islink(cp) and gr.md5_of(cp) == gr.md5_of(gf) and os.path.exists(cp + '.fai') and mf3['copied_from'] == gf and 'linked_from' not in mf3)
    # a registration made by v16 (symlink + linked_from) is converted by materialize
    old_dir = gr.entry_dir(fp3, 'OldLinked'); os.makedirs(old_dir)
    os.symlink(gf, old_dir + 'YH_Y.fasta')
    json.dump({'key': 'OldLinked', 'fasta': gf, 'linked_from': gf, 'md5_fasta': gr.md5_of(gf)}, open(old_dir + 'manifest.json', 'w'))
    mm2 = gr.materialize(fp3, 'OldLinked')
    _ok('a symlinked registration from an earlier version is converted into a real, verified copy', not os.path.islink(old_dir + 'YH_Y.fasta') and os.path.isfile(old_dir + 'YH_Y.fasta')
        and mm2['fasta'] == old_dir + 'YH_Y.fasta' and 'linked_from' not in mm2 and gr.md5_of(old_dir + 'YH_Y.fasta') == gr.md5_of(gf))
    open(old_dir + 'YH_Y.fasta.part', 'w').write('x')
    try:
        gr.copy_verified(gf, old_dir + 'bad.fa', expect_md5='0' * 32); refused = False
    except RuntimeError:
        refused = True
    _ok('a copy whose md5 differs from the expected one is refused and leaves no half-written file', refused and not os.path.exists(old_dir + 'bad.fa'))
    # runner archives small inputs from outside the project, not files already inside it, not big files
    pr = types.SimpleNamespace(fm=types.SimpleNamespace(localMasterDir=loc, cloudMasterDir=loc + 'cloud/'), root=loc + 'ARCH/G/', genome=ref, cfg={'runner': {'upload': False, 'archive_max_mb': 1}}, run=None)
    ext = tempfile.mkdtemp() + '/'; open(ext + 'small.tsv', 'w').write('a\tb\n1\t2\n'); open(ext + 'big.bin', 'wb').write(b'0' * 1_500_000)
    os.makedirs(loc + 'ARCH/G/', exist_ok=True); open(loc + 'ARCH/G/already_here.txt', 'w').write('x')
    with Run(pr, 'toy', 'archive') as r:
        r.add_input(ext + 'small.tsv'); r.add_input(ext + 'big.bin'); r.add_input(loc + 'ARCH/G/already_here.txt'); r.add_input(ext + 'missing.txt')
        rd2 = r.dir
    recs = {os.path.basename(x['path']): x for x in r.inputs}
    _ok('runner copies small outside inputs into <run>/inputs/, skips big files and files already in the project, flags missing ones',
        recs['small.tsv'].get('archived') == 'inputs/small.tsv' and os.path.exists(rd2 + 'inputs/small.tsv') and 'archived' not in recs['big.bin'] and 'archived' not in recs['already_here.txt'] and recs['missing.txt'].get('missing'))

    # ---------------- failures seen on the real YH / OB assemblies (LG5, LG2, LG13)
    def blk(q, qlen, qs, qe, st, tname, ts, te): return line(q, qlen, qs, qe, st, tname, 50 * M, ts, te)
    # (1) LG5-like: a real inversion called by two assemblies AND repeat-driven calls elsewhere with the same number of supporters: the real one must rank first
    calls = pd.DataFrame([dict(assembly='LatrZW', method='C', left_lo=3.68 * M, left_hi=3.68 * M, right_lo=5.04 * M, right_hi=5.04 * M),
                          dict(assembly='MezeOBm', method='C', left_lo=3.69 * M, left_hi=3.69 * M, right_lo=5.04 * M, right_hi=5.04 * M),
                          dict(assembly='LatrZW', method='A', left_lo=19.657 * M, left_hi=19.658 * M, right_lo=26.726 * M, right_hi=26.727 * M),
                          dict(assembly='MezeOBm', method='A', left_lo=19.657 * M, left_hi=19.658 * M, right_lo=26.726 * M, right_hi=26.727 * M),
                          dict(assembly='LatrZW', method='B', left_lo=3.38 * M, left_hi=3.39 * M, right_lo=4.66 * M, right_hi=4.67 * M)])
    cs = paf.consensus(calls, 150_000)
    _ok('consensus ranks a call with both flanks (A) above repeat-driven C/B calls, even when they have as many supporting calls',
        cs.iloc[0].confidence == 'high' and abs(cs.iloc[0].left_breakpoint - 19.6575 * M) < 5000 and cs.iloc[0].n_calls == 2 and set(cs.confidence) == {'high', 'medium', 'low'}, str(cs[['left_breakpoint', 'confidence']].values.tolist()))
    # (2) LG2-like: reverse-complemented contig [rev end 43.5-48.8][+ 19.8-43.5, speckled with tiny reverse blocks][rev 14.4-20 and 2-4.3][+ hairpin 0-2].
    #     The inversion is ref 19.8-43.5; the old caller reported 2.0-23.8 (a '-' run with '+' flanks in the wrong places)
    TL2 = 'NC_T2L'
    lg2 = [blk('y2', 39.2 * M, 0, 5.3 * M, '-', TL2, 43.5 * M, 48.8 * M),
           blk('y2', 39.2 * M, 5.3 * M, 12.0 * M, '+', TL2, 19.8 * M, 26.5 * M), blk('y2', 39.2 * M, 12.0 * M, 12.06 * M, '-', TL2, 26.5 * M, 26.56 * M),
           blk('y2', 39.2 * M, 12.06 * M, 20.0 * M, '+', TL2, 26.56 * M, 34.5 * M), blk('y2', 39.2 * M, 20.0 * M, 20.05 * M, '-', TL2, 34.5 * M, 34.55 * M),
           blk('y2', 39.2 * M, 20.05 * M, 29.7 * M, '+', TL2, 34.55 * M, 43.5 * M),
           blk('y2', 39.2 * M, 29.8 * M, 35.0 * M, '-', TL2, 14.4 * M, 20.0 * M), blk('y2', 39.2 * M, 35.1 * M, 37.0 * M, '-', TL2, 2.0 * M, 4.3 * M),
           blk('y2', 39.2 * M, 37.1 * M, 39.0 * M, '+', TL2, 0, 2.0 * M)]
    d2l = paf.parse_paf(lg2)
    r2 = paf.call_method_AC(d2l, TL2)
    _ok('LG2-like contig: speckles are absorbed, the hairpin is not called, the inversion is ref 19.8-43.5 Mb (found once, with both flanks)',
        len(r2) == 1 and r2.method.iloc[0] == 'A' and abs(r2.left_hi.iloc[0] - 20.0 * M) < 1000 and abs(r2.right_lo.iloc[0] - 43.5 * M) < 1000, str(r2[['method', 'left_lo', 'left_hi', 'right_lo', 'right_hi']].values.tolist()))
    r2_old = paf.call_method_AC(d2l, TL2, despeckle=0)
    _ok('without despeckling the same contig is NOT resolved (so the despeckling is what fixes it)', not (len(r2_old) == 1 and abs(r2_old.right_lo.iloc[0] - 43.5 * M) < 1000))
    runs = paf.strand_runs(paf.blocks_for_contig(d2l, 'y2', TL2, 20_000), despeckle=300_000)
    _ok('despeckled runs: - + - +  (4 runs)', [r['strand'] for r in runs] == ['-', '+', '-', '+'])
    # (3) flank-run geometry inconsistent with its strand (a '+' run between two '-' runs lying in the wrong order) is not an inversion
    bad = [blk('b', 30 * M, 0, 5 * M, '-', TL2, 10 * M, 15 * M), blk('b', 30 * M, 5.1 * M, 20 * M, '+', TL2, 15.1 * M, 30 * M), blk('b', 30 * M, 20.1 * M, 25 * M, '-', TL2, 30.1 * M, 35 * M)]
    _ok("a '+' run between '-' runs that lie in the wrong order is not an inversion", len(paf.call_method_AC(paf.parse_paf(bad), TL2)) == 0 or paf.call_method_AC(paf.parse_paf(bad), TL2).method.iloc[0] != 'A')
    # (4) LG13-like: a flank adjacency megabases wide is not a breakpoint
    wide = [blk('w', 30 * M, 0, 2.4 * M, '+', TL2, 0, 2.4 * M), blk('w', 30 * M, 2.5 * M, 14 * M, '-', TL2, 12.0 * M, 23.5 * M), blk('w', 30 * M, 14.1 * M, 28 * M, '+', TL2, 23.6 * M, 37.5 * M)]
    rw = paf.call_method_AC(paf.parse_paf(wide), TL2)
    _ok('a breakpoint interval wider than max_interval (flank ends 9.6 Mb before the inverted run) is rejected', len(rw[rw.method == 'A']) == 0 if len(rw) else True)
    # the call step writes the despeckled runs, one dotplot panel per contig, and the polarity text no longer talks about strands
    rd_ = last('inversion_breakpoints', 'call')
    _ok('call step writes the strand-run table next to the dotplots', os.path.exists(rd_ + 'LG5_strand_runs.tsv') and os.path.exists(rd_ + 'LG2_strand_runs.tsv'))
    _ok('polarity text is independent of the strand of the contig', 'strand' not in json.load(open(rd_ + 'LG2_polarity.json'))['polarity'].lower().replace('orientation', ''))

    # ---------------- the reference's known misassembly (LG2 4.0-14.4 Mb) is ignored by the breakpoint caller and the centromere scan
    fakep = types.SimpleNamespace(genome='Mzebra_GT3_NCBI', chrom_for_lg=lambda lg: {'LG2': 'NC_135168.1'}[lg])
    ex = gr.excluded_regions(fakep)
    _ok('config: known misassembly of the reference resolved to the LG2 accession', ex == {'NC_135168.1': [(4_000_000, 14_400_000)]}, str(ex))
    junk = paf.parse_paf([line('j', 20 * M, 0, 7 * M, '-', 'NC_135168.1', 48 * M, 4.3 * M, 11.1 * M), line('k', 20 * M, 0, 6 * M, '+', 'NC_135168.1', 48 * M, 20 * M, 26 * M),
                          line('m', 20 * M, 0, 3 * M, '+', 'NC_135168.1', 48 * M, 13.5 * M, 15.5 * M)])
    kept = paf.drop_blocks_in_regions(junk, ex['NC_135168.1'])
    _ok('blocks inside the misassembled region are dropped, blocks outside or only touching it are kept', list(kept.qname) == ['k', 'm'], str(list(kept.qname)))
    _ok('without the exclusion the repeat junk is called as an inversion (the old LG2 "YH:C 4.3-11.1 Mb")', len(paf.call_method_AC(junk, 'NC_135168.1')) == 1 and len(paf.call_method_AC(kept, 'NC_135168.1')) == 0)
    # centromere scan: a big satellite chunk in the excluded region must not win; the real array outside it is called
    d5 = tempfile.mkdtemp(); fa5 = d5 + '/x.fa'
    with open(fa5, 'w') as fh:
        fh.write('>NC_135168.1\n'); fh.write(('ACGT' * 25 + '\n') * 480000)                       # 48 Mb
    ann5 = d5 + '/x.out'
    with open(ann5, 'w') as fh:
        fh.write('h1\nh2\n\n')
        fh.write('  900  5.0  0.0  0.0  NC_135168.1  6000001  12000000  (1)  +  CHUNK  Satellite  1 2 3 1\n')                 # misplaced repeat chunk
        fh.write('  900  5.0  0.0  0.0  NC_135168.1  26900001  27700000  (1)  +  CEN  Satellite  1 2 3 2\n')                  # the real centromere array
    cfg5 = csteps.cfg_of(); cfg5.update(min_chrom_length=1_000_000)
    c_no, _ = csteps.scan_genome('toy', fa5, ann5, cfg5)
    c_ex, _ = csteps.scan_genome('toy', fa5, ann5, cfg5, exclude=ex)
    _ok('centromere scan: the misplaced chunk wins without the mask, the real array wins with it',
        6_000_000 <= c_no.iloc[0].start < 12_500_000 and 26_000_000 <= c_ex.iloc[0].start and c_ex.iloc[0].end <= 28_500_000 and c_ex.iloc[0].masked_regions == '4000000-14400000', f"{c_no.iloc[0].start} -> {c_ex.iloc[0].start}")

    # ---------------- published centromeres (Table 3) lifted into the reference, including the 10.4 Mb chunk the reference has and the published assembly does not
    import openpyxl
    from .analyses.centromeres import published as pb
    wbp = tempfile.mkdtemp() + '/t3.xlsx'
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = 'Table3 Centromere annotation'
    ws.append(['genome', 'chr', 'start', 'end', 'Size', 'Centroid1', 'Centroid2', 'Distance', 'Gap cen', 'Telo L', 'Telo R', 'Ungapped chr', 'T2T'])
    ws.append(['Maylandia_zebra', 'chr2', 17170767, 17765435, '=D2-C2', 17378589, None, None, 'No', 'Yes', 'Yes', 'Yes', '=AND(J2="Yes",K2="Yes",L2="Yes")'])
    ws.append(['Maylandia_zebra', 'chr10', 7259076, 7998165, '=D3-C3', 7424819, None, None, 'No', 'Yes', 'Yes', 'Yes', '=FALSE'])
    ws.append(['Maylandia_zebra', 'scf00000001', 100, 500, '=D4-C4', None, None, None, 'No', None, None, None, None])
    ws.append(['Astatotilapia_burtoni', 'chr2', 1, 2, '=D5-C5', None, None, None, 'No', None, None, None, None])
    wb.save(wbp)
    tb = pb.read_published(wbp, 'Table3 Centromere annotation', 'Maylandia_zebra')
    _ok('published table: one genome, chrN rows only, formulas not required', list(tb.lg) == ['LG2', 'LG10'] and list(tb['size']) == [594668, 739089] and tb.centroid1.iloc[0] == 17378589)
    NT2, NT10, NEW2, NEW10 = 'NC_135168.1', 'NC_135176.1', 'OZ_chr2', 'OZ_chr10'
    # reference LG2 has the 10.4 Mb chunk at 4.0-14.4 Mb; the published chr2 does not. chr10: the centromere array itself is not alignable (gap 7.3-8.1 in the published assembly)
    pl = [line(NEW2, 38.4 * M, 0, 4.0 * M, '+', NT2, 48.8 * M, 0, 4.0 * M, f'{int(4.0 * M)}M'), line(NEW2, 38.4 * M, 4.0 * M, 38.4 * M, '+', NT2, 48.8 * M, 14.4 * M, 48.8 * M, f'{int(34.4 * M)}M'),
          line(NEW10, 35.0 * M, 0, 7.2 * M, '+', NT10, 35.0 * M, 0, 7.2 * M, f'{int(7.2 * M)}M'), line(NEW10, 35.0 * M, 8.1 * M, 35.0 * M, '+', NT10, 35.0 * M, 8.1 * M, 35.0 * M, f'{int(26.9 * M)}M')]
    pdf = paf.parse_paf(pl)
    cm = pd.DataFrame([(NEW2, NT2, int(38.4 * M), 1.0), (NEW10, NT10, int(35 * M), 1.0)], columns=['contig', 'reference_chrom', 'aligned_bp', 'share'])
    lg2c = {'LG2': NT2, 'LG10': NT10}
    mp = pb.map_published_chr(cm, lg2c, ['LG2', 'LG10', 'LG99'])
    _ok('published chromosome -> assembly sequence through the dominant reference chromosome (and "none" when nothing aligns)', mp['LG2']['seq'] == NEW2 and mp['LG10']['seq'] == NEW10 and mp['LG99'] is None)
    lifted = pb.lift_all(pdf, tb, mp, lg2c, 3_000_000)
    r2 = lifted[lifted.lg == 'LG2'].iloc[0]
    _ok('LG2: published centromere 17.17-17.77 Mb lifts to ~27.57-28.17 Mb in the reference (+10.4 Mb: the chunk), centroid too',
        r2.lift_status == 'ok' and abs(r2.gt3_start - 27_570_767) < 2 and abs(r2.gt3_end - 28_165_435) < 2 and abs(r2.centroid1_ref - 27_778_589) < 2, f'{r2.gt3_start}-{r2.gt3_end}')
    r10 = lifted[lifted.lg == 'LG10'].iloc[0]
    _ok('LG10: an unalignable centromere array is bounded by its alignable flanks (7.2-8.1 Mb)', r10.lift_status == 'ok' and abs(r10.gt3_start - 7.2 * M) < 30_000 and abs(r10.gt3_end - 8.1 * M) < 30_000, f'{r10.lift_status} {r10.gt3_start}-{r10.gt3_end}')
    straddle = paf.parse_paf([line(NEW10, 35.0 * M, 0, 7.2 * M, '+', NT10, 35.0 * M, 0, 7.2 * M, f'{int(7.2 * M)}M'), line(NEW10, 35.0 * M, 8.1 * M, 35.0 * M, '-', NT10, 35.0 * M, 0.5 * M, 27.4 * M, f'{int(26.9 * M)}M')])
    rs = pb.lift_all(straddle, tb, mp, lg2c, 3_000_000).query("lg == 'LG10'").iloc[0]
    _ok('opposite-strand flanks on the expected chromosome still give the region BETWEEN the flank landings (LG13: 32.49-35.76 Mb)', abs(rs.bounded_start - 7.2 * M) < 30_000 and abs(rs.bounded_end - 27.4 * M) < 60_000 and pd.isna(rs.gt3_start))
    _ok('a clean lift has bounded == lifted interval', abs(r10.bounded_start - r10.gt3_start) < 1 and abs(r10.bounded_end - r10.gt3_end) < 1)
    _ok('flanks on opposite strands (an interval across a breakpoint): flagged and NO interval coordinates are reported, but the position and strand each flank lands on are',
        'opposite strands' in rs.lift_status and pd.isna(rs.gt3_start) and rs.left_flank_strand == '+' and rs.right_flank_strand == '-' and abs(rs.left_flank_ref - 7.2 * M) < 30_000 and abs(rs.right_flank_ref - 27.4 * M) < 60_000,
        f'{rs.left_flank_ref}{rs.left_flank_strand} {rs.right_flank_ref}{rs.right_flank_strand}')
    scan_df = pd.DataFrame([dict(lg='LG2', start=26_900_000, end=27_900_000, confidence='high'), dict(lg='LG10', start=0, end=300_000, confidence='low')])
    cmp_ = pb.compare_with_scan(lifted, scan_df, 1_000_000)
    _ok('comparison with the satellite scan: LG2 agrees (after the chunk is accounted for), LG10 (scan at 0-0.3 Mb) disagrees', cmp_[cmp_.lg == 'LG2'].scan_vs_published.iloc[0] == 'agrees' and cmp_[cmp_.lg == 'LG10'].scan_vs_published.iloc[0].startswith('disagrees'))

    # a satellite hit on ANOTHER chromosome right next to the centromere must not be used as a flank (LG3, LG7, LG13, LG20 failed this way on the real data)
    NT3, NT9, NEW3 = 'NC_135169.1', 'NC_135175.1', 'OZ_chr3'
    tb3 = pd.DataFrame([dict(chr='chr3', lg='LG3', start=29_500_000, end=30_000_000, size=500_000, centroid1=np.nan, centroid2=np.nan, gap_in_centromere='No', t2t=False)])
    pl3 = [line(NEW3, 40 * M, 0, 29.0 * M, '+', NT3, 60 * M, 0, 29.0 * M, f'{int(29.0 * M)}M'), line(NEW3, 40 * M, 29.1 * M, 29.4 * M, '+', NT9, 37 * M, 5.0 * M, 5.3 * M, f'{int(0.3 * M)}M'),
           line(NEW3, 40 * M, 30.5 * M, 40 * M, '+', NT3, 60 * M, 30.5 * M, 40 * M, f'{int(9.5 * M)}M')]
    cm3 = pd.DataFrame([(NEW3, NT3, int(39 * M), 0.98)], columns=['contig', 'reference_chrom', 'aligned_bp', 'share'])
    l3 = pb.lift_all(paf.parse_paf(pl3), tb3, pb.map_published_chr(cm3, {'LG3': NT3}, ['LG3']), {'LG3': NT3}, 3_000_000).iloc[0]
    _ok('lift ignores alignments to other chromosomes next to the centromere (cross-chromosome satellite hits no longer break the flanks)',
        l3.lift_status == 'ok' and abs(l3.gt3_start - 29.0 * M) < 30_000 and abs(l3.gt3_end - 30.5 * M) < 30_000, f'{l3.lift_status} {l3.gt3_start}-{l3.gt3_end}')

    # figure: a published centromere that lifts cleanly is drawn from the lifted interval; one with opposite-strand flanks from the region between the flanks; an OLD lift table (no bounded columns) still works
    def lift_run(df):
        with Run(proj, 'centromeres', 'lift') as r_:
            df.to_csv(r_.dir + 'published_centromeres_lifted.tsv', sep='\t', index=False)
    base = dict(lg=['LG2', 'LG5'], gt3_start=[27_129_591, np.nan], gt3_end=[27_702_624, np.nan], lift_status=['ok', 'flanks on opposite strands (interval straddles an inversion breakpoint)'])
    lift_run(pd.DataFrame(dict(base, bounded_start=[27_129_591, 32_490_569], bounded_end=[27_702_624, 35_759_951])))
    steps.figure_step(proj)
    ovb = pd.read_csv(last('inversion_breakpoints', 'figure') + 'inversion_overview.tsv', sep='\t').set_index('lg')
    _ok('overview: clean lift -> published; opposite-strand flanks -> published-bounded interval between the flank landings',
        ovb.loc['LG2', 'cen_start'] == 27_129_591 and ovb.loc['LG2', 'cen_confidence'] == 'published' and ovb.loc['LG5', 'cen_start'] == 32_490_569 and ovb.loc['LG5', 'cen_end'] == 35_759_951 and ovb.loc['LG5', 'cen_confidence'] == 'published-bounded')
    lift_run(pd.DataFrame(base))                                                                  # a table from an older version: no bounded columns
    steps.figure_step(proj)
    _ok('an older lift table without bounded columns does not break the figure step', True)

    # ---------------- array-level centromere definition: largest array of ONE family more than 500 kb from both ends
    ar = cen.merge_arrays([100_000, 130_000, 200_000, 900_000, 905_000, 20_000_000, 20_040_000, 20_120_000, 20_130_000, 31_000_000], [120_000, 160_000, 260_000, 905_000, 910_000, 20_030_000, 20_100_000, 20_125_000, 20_450_000, 31_050_000], 50_000)
    _ok('features of a family are merged into arrays by a gap threshold; covered bp excludes gaps', [(a, b, c, n) for a, b, c, n in ar] == [(100_000, 260_000, 110_000, 3), (900_000, 910_000, 10_000, 2), (20_000_000, 20_450_000, 415_000, 4), (31_000_000, 31_050_000, 50_000, 1)], str(ar))
    L40 = 40_000_000
    res = cen.call_by_arrays([(0, 700_000, 650_000, 9), (18_000_000, 18_450_000, 420_000, 12), (30_000_000, 30_120_000, 110_000, 3), (39_700_000, 39_950_000, 240_000, 5)], L40, 500_000, 100_000)
    _ok('the centromere is the largest INTERNAL array; a bigger terminal array (here 650 kb at the chromosome start) is reported but never chosen',
        res['best']['start'] == 18_000_000 and res['best']['size'] == 420_000 and res['best_overall']['start'] == 0 and res['overall_is_terminal'] and res['n_internal'] == 2 and not res['ambiguous'])
    res2 = cen.call_by_arrays([(10_000_000, 10_400_000, 380_000, 6), (25_000_000, 25_350_000, 330_000, 5)], L40)
    _ok('two comparable internal arrays are flagged ambiguous', res2['ambiguous'] and res2['best']['start'] == 10_000_000)
    _ok('relation to an inversion: pericentric / paracentric / breakpoint inside the centromere',
        cen.relation_to_inversion(27.13e6, 27.60e6, 19.75e6, 43.56e6) == 'pericentric' and cen.relation_to_inversion(7.0e6, 7.42e6, 11.82e6, 29.90e6) == 'paracentric'
        and cen.relation_to_inversion(32.49e6, 32.92e6, 2.40e6, 23.07e6) == 'paracentric' and cen.relation_to_inversion(11.7e6, 11.9e6, 11.82e6, 29.9e6).startswith('breakpoint'))
    # arrays_for_genome on a small synthetic genome: family-47 arrays vs another satellite family and a terminal array
    da = tempfile.mkdtemp(); faa = da + '/a.fa'
    with open(faa, 'w') as fh:
        fh.write('>chrA\n'); fh.write(('ACGT' * 25 + '\n') * 400000)                              # 40 Mb
    anna = da + '/a.out'
    with open(anna, 'w') as fh:
        fh.write('h1\nh2\n\n')
        fh.write('  900  5.0  0.0  0.0  chrA  1  600000  (1)  +  rnd-4_family-47  Satellite  1 2 3 1\n')                                  # terminal family-47 array (600 kb at the start)
        for k in range(0, 420_000, 60_000):
            fh.write(f'  900  5.0  0.0  0.0  chrA  {18_000_001 + k}  {18_000_000 + k + 55_000}  (1)  +  rnd-4_family-47  Satellite  1 2 3 2\n')   # internal array, small gaps
        fh.write('  900  5.0  0.0  0.0  chrA  25000001  26500000  (1)  +  OTHERSAT  Satellite  1 2 3 3\n')                              # bigger array of ANOTHER family: ignored
    cfga = csteps.cfg_of(); cfga.update(min_chrom_length=1_000_000)
    ca, aa = csteps.arrays_for_genome('toy', faa, anna, cfga)
    ra = ca.iloc[0]
    _ok('arrays_for_genome uses only the configured family, skips the terminal array and the other family, reports the terminal one',
        ra.start == 18_000_000 and 18_400_000 <= ra.end <= 18_460_000 and ra.largest_array_is_terminal and ra.n_internal_arrays == 1 and ra.family == 'rnd-4_family-47', f'{ra.start}-{ra.end} size {ra["size"]}')
    # the overview prefers the array call over the published lift and the scan, and adds the relation to the inversion
    with Run(proj, 'centromeres', 'arrays') as r_:
        pd.DataFrame([dict(genome=ref, sequence=T2, lg='LG2', start=27_130_000, end=27_600_000, size=454_000, ambiguous=False),
                      dict(genome=ref, sequence=T5, lg='LG5', start=32_490_000, end=32_920_000, size=419_000, ambiguous=True)]).to_csv(r_.dir + 'centromeres_arrays.tsv', sep='\t', index=False)
    steps.figure_step(proj)
    ova = pd.read_csv(last('inversion_breakpoints', 'figure') + 'inversion_overview.tsv', sep='\t').set_index('lg')
    _ok('overview: array call wins (published lift and scan are only fall-backs), size and pericentric/paracentric are added',
        ova.loc['LG2', 'cen_start'] == 27_130_000 and ova.loc['LG2', 'cen_confidence'] == 'array' and ova.loc['LG5', 'cen_confidence'] == 'array-ambiguous' and ova.loc['LG2', 'cen_size'] == 454_000
        and ova.loc['LG2', 'cen_relation'] == 'pericentric', str(ova[['cen_start', 'cen_confidence', 'cen_relation']].to_dict('index')))

    # ---------------- chr11-like TANDEM (double) inversion: A B' C' D with B and C each inverted in place; chr9-like single inversion; the same contig reverse-complemented
    TD = 'NC_T11'; LD = 40 * M
    dbl = [line('d', 40.6 * M, 0, 6.7 * M, '+', TD, LD, 0, 6.7 * M),
           line('d', 40.6 * M, 6.8 * M, 12.3 * M, '-', TD, LD, 12.3 * M, 17.8 * M), line('d', 40.6 * M, 12.4 * M, 17.9 * M, '-', TD, LD, 6.8 * M, 12.3 * M),
           line('d', 40.6 * M, 18.0 * M, 23.5 * M, '-', TD, LD, 23.2 * M, 28.7 * M), line('d', 40.6 * M, 23.6 * M, 29.1 * M, '-', TD, LD, 17.8 * M, 23.3 * M),
           line('d', 40.6 * M, 29.2 * M, 40.6 * M, '+', TD, LD, 28.8 * M, 40.0 * M)]
    for name, ls in (('as given', dbl), ('contig reverse-complemented', [flip(x) for x in dbl])):
        rd = paf.call_method_AC(paf.parse_paf(ls), TD)
        segs = [tuple(int(v) for v in s_.split('-')) for s_ in rd.segments.iloc[0].split(';')] if len(rd) else []
        _ok(f'double inversion ({name}): ONE inversion 6.8-28.7 Mb made of TWO segments with the internal breakpoint at 17.8 Mb',
            len(rd) == 1 and rd.method.iloc[0] == 'A' and int(rd.n_segments.iloc[0]) == 2 and abs(segs[0][0] - 6.8 * M) < 1000 and abs(segs[0][1] - 17.8 * M) < 1000 and abs(segs[1][1] - 28.7 * M) < 1000
            and abs(int(rd.internal_breakpoints.iloc[0]) - 17.8 * M) < 100_000 and abs(rd.left_hi.iloc[0] - 6.8 * M) < 1000 and abs(rd.right_lo.iloc[0] - 28.7 * M) < 1000, str(rd[['n_segments', 'segments']].values.tolist()))
    sgl = [line('s', 40 * M, 0, 15.9 * M, '+', TD, LD, 0, 15.9 * M), line('s', 40 * M, 16.0 * M, 24.5 * M, '-', TD, LD, 24.5 * M, 33.0 * M), line('s', 40 * M, 24.6 * M, 32.6 * M, '-', TD, LD, 16.5 * M, 24.5 * M),
           line('s', 40 * M, 32.7 * M, 40 * M, '+', TD, LD, 32.7 * M, 40.0 * M)]
    rs1 = paf.call_method_AC(paf.parse_paf(sgl), TD)
    _ok('a single inversion is ONE segment (the blocks of one inverted run fall monotonically in the reference)', len(rs1) == 1 and int(rs1.n_segments.iloc[0]) == 1 and rs1.internal_breakpoints.iloc[0] == '')
    noisy = list(sgl); noisy.insert(2, line('s', 40 * M, 24.5 * M, 24.53 * M, '-', TD, LD, 30.0 * M, 30.03 * M))          # a 30 kb repeat hit that jumps up inside the inversion
    rs2 = paf.call_method_AC(paf.parse_paf(noisy), TD)
    _ok('a small repeat block that jumps inside one segment does not create a second segment', len(rs2) == 1 and int(rs2.n_segments.iloc[0]) == 1, str(rs2[['n_segments', 'segments']].values.tolist()))
    noisy2 = list(sgl); noisy2.insert(2, line('s', 40 * M, 24.5 * M, 24.65 * M, '-', TD, LD, 30.0 * M, 30.15 * M))          # a 150 kb repeat hit (above the small-block cut-off) that jumps up inside the inversion
    rs3 = paf.call_method_AC(paf.parse_paf(noisy2), TD)
    _ok('even a 150 kb out-of-line block does not split one inversion (overlapping segment extents are merged)', len(rs3) == 1 and int(rs3.n_segments.iloc[0]) == 1, str(rs3[['n_segments', 'segments']].values.tolist()))
    cs_ = paf.consensus(pd.concat([paf.call_method_AC(paf.parse_paf(dbl), TD).assign(assembly='YH'), paf.call_method_AC(paf.parse_paf([flip(x) for x in dbl]), TD).assign(assembly='fDipLim2')], ignore_index=True), 150_000)
    _ok('consensus of two assemblies keeps the segments of the best call', len(cs_) == 1 and int(cs_.n_segments.iloc[0]) == 2 and cs_.n_calls.iloc[0] == 2 and cs_.segments.iloc[0].count(';') == 1)
    # the figure script draws the two segments from the overview table (smoke test)
    import subprocess as sp, sys as sy
    from .common_light import PKG_ROOT
    tsv = tempfile.mkdtemp() + '/o.tsv'
    pd.DataFrame([dict(lg='LG11', length=38_663_802, inv_start=6_800_000, inv_end=28_700_000, inv_segments='6800000-17800000;17800000-28700000', cen_start=5_038_192.0, cen_end=5_616_369.0, cen_confidence='array', cen_size=500000.0)]).to_csv(tsv, sep='\t', index=False)
    rc = sp.run([sy.executable, os.path.join(PKG_ROOT, 'figures', 'fig_inversions_to_scale.py'), '--overview', tsv, '--out', os.path.dirname(tsv) + '/fig'], capture_output=True, text=True)
    svg = open(os.path.dirname(tsv) + '/fig.svg', encoding='utf-8').read() if rc.returncode == 0 else ''
    _ok('figure script draws a double inversion as two segments with its label', rc.returncode == 0 and '2 segments' in svg and 'double inversion' in svg, rc.stderr[-200:])

    # ---------------- LG11 as it really looked: a forward flank whose extent is stretched by a stray block, then a double inversion joined by another stray block
    TL, LL = 'NC_T11b', 40 * M
    lg11 = [line('q11', 41 * M, 0.07 * M, 4.0 * M, '+', TL, LL, 0.25 * M, 4.18 * M), line('q11', 41 * M, 4.0 * M, 8.12 * M, '+', TL, LL, 4.18 * M, 8.30 * M),
            line('q11', 41 * M, 8.12 * M, 8.20 * M, '+', TL, LL, 18.40 * M, 18.48 * M),                         # stray forward block: the internal-breakpoint homology, absorbed into the left flank
            line('q11', 41 * M, 8.2 * M, 13.2 * M, '-', TL, LL, 13.4 * M, 18.4 * M), line('q11', 41 * M, 13.2 * M, 18.3 * M, '-', TL, LL, 8.3 * M, 13.4 * M),
            line('q11', 41 * M, 18.3 * M, 18.35 * M, '+', TL, LL, 18.40 * M, 18.45 * M),                         # the dot between the two inverted halves
            line('q11', 41 * M, 18.4 * M, 24.4 * M, '-', TL, LL, 24.4 * M, 30.4 * M), line('q11', 41 * M, 24.4 * M, 30.4 * M, '-', TL, LL, 18.4 * M, 24.4 * M),
            line('q11', 41 * M, 30.5 * M, 40.6 * M, '+', TL, LL, 30.4 * M, 40.5 * M)]
    d11 = paf.parse_paf(lg11)
    runs11 = paf.strand_runs(paf.blocks_for_contig(d11, 'q11', TL, 20_000), 300_000)
    _ok('LG11 geometry: the left forward flank really is stretched to 18.48 Mb by the stray block (the reason the caller dropped the inversion)', runs11[0]['strand'] == '+' and runs11[0]['tmax'] > 18 * M, str([(r['strand'], r['tmin'], r['tmax']) for r in runs11]))
    _core = paf.core_extent(runs11[0])
    _ok('core_extent measures the flank on its colinear core: 0.25-8.30 Mb', abs(_core[0] - 0.25 * M) < 1000 and abs(_core[1] - 8.30 * M) < 1000, str(_core))
    r11 = paf.call_method_AC(d11, TL)
    segs11 = [tuple(int(v) for v in x.split('-')) for x in r11.segments.iloc[0].split(';')] if len(r11) else []
    _ok('LG11 double inversion is now called: ONE inversion 8.30-30.40 Mb, TWO segments, internal breakpoint 18.4 Mb',
        len(r11) == 1 and r11.method.iloc[0] == 'A' and int(r11.n_segments.iloc[0]) == 2 and abs(r11.left_hi.iloc[0] - 8.30 * M) < 20_000 and abs(r11.right_lo.iloc[0] - 30.40 * M) < 20_000
        and abs(int(r11.internal_breakpoints.iloc[0]) - 18.4 * M) < 100_000 and abs(segs11[0][1] - 18.4 * M) < 100_000, str(r11[['method', 'n_segments', 'left_lo', 'left_hi', 'right_lo', 'right_hi', 'internal_breakpoints']].values.tolist()))
    _orig_core = paf.core_extent
    paf.core_extent = lambda run, offset_tol=1_000_000: (run['tmin'], run['tmax'])                      # the previous behaviour: every block counts
    try:
        _ok('without the core extent the same data gives NO call (what happened on the real LG11)', len(paf.call_method_AC(d11, TL)) == 0)
    finally:
        paf.core_extent = _orig_core
    fl2 = [(0.25 * M, 0.25 * M + 4 * M, 0.07 * M, 4.07 * M), (0.25 * M + 4 * M, 0.25 * M + 8 * M, 4.07 * M, 8.07 * M)]
    rev = dict(strand='-', tmin=1, tmax=2, blocks=[(10 * M, 11 * M, 5 * M, 6 * M), (9 * M, 10 * M, 6 * M, 7 * M), (30 * M, 30.05 * M, 6.9 * M, 6.95 * M)])      # '-' flank with one stray block 20 Mb away
    _ce = paf.core_extent(rev)
    _ok("core_extent also works for '-' flanks (anti-diagonal) and drops the stray block", (_ce[0], _ce[1]) == (9 * M, 11 * M), str(_ce))
    _ok('core_extent leaves short runs (fewer than 3 blocks) alone', paf.core_extent(dict(strand='+', tmin=5, tmax=9, blocks=[(5, 6, 0, 1), (8, 9, 3, 4)])) == (5, 9))

    # ---------------- LG9: three assemblies whose breakpoint intervals overlap are ONE inversion even if two midpoints are 181 kb apart
    lg9 = pd.DataFrame([dict(assembly='fDipLim2', method='A', left_lo=14_511_484, left_hi=14_622_365, right_lo=32_763_232, right_hi=32_878_069, n_segments=1, segments='15313086-32161792', internal_breakpoints=''),
                        dict(assembly='YH', method='A', left_lo=14_584_786, left_hi=14_912_159, right_lo=32_763_232, right_hi=32_775_869, n_segments=1, segments='15554102-31239177', internal_breakpoints=''),
                        dict(assembly='fAulStu2', method='A', left_lo=14_585_294, left_hi=14_977_232, right_lo=32_390_715, right_hi=32_763_232, n_segments=1, segments='15291930-31238831', internal_breakpoints=''),
                        dict(assembly='fLetChi2', method='A', left_lo=9_703_206, left_hi=9_713_388, right_lo=12_301_319, right_hi=12_350_695, n_segments=1, segments='', internal_breakpoints='')])
    c9 = paf.consensus(lg9)
    _ok('LG9: the three overlapping calls form one consensus (3 calls, honest interval = their union); the unrelated 2.6 Mb call stays separate',
        len(c9) == 2 and c9.n_calls.iloc[0] == 3 and c9.left_min.iloc[0] == 14_511_484 and c9.left_max.iloc[0] == 14_977_232 and c9.right_min.iloc[0] == 32_390_715 and c9.right_max.iloc[0] == 32_878_069
        and c9.n_calls.iloc[1] == 1, str(c9[['left_breakpoint', 'right_breakpoint', 'n_calls']].values.tolist()))
    _ok('...whereas midpoint-only clustering (the previous behaviour) split them into three rows', len(paf.consensus(lg9, 150_000, gap=-1)) == 4)

    # ---------------- LG9 with the real per-assembly intervals: the shared interval, not the union, is the estimate
    real9 = pd.DataFrame([dict(assembly='fDipLim2', method='A', left_lo=14_511_484, left_hi=14_622_365, right_lo=32_763_232, right_hi=32_878_069),
                          dict(assembly='YH', method='A', left_lo=14_584_786, left_hi=14_912_159, right_lo=32_763_232, right_hi=32_775_869),
                          dict(assembly='fAulStu2', method='A', left_lo=14_585_294, left_hi=14_977_232, right_lo=32_390_715, right_hi=32_763_232)])
    r9 = paf.consensus(real9).iloc[0]
    _ok('LG9 (three real assemblies): left breakpoint = middle of the 37 kb shared interval, right = the single shared position; the union is kept as the outer bound',
        (r9.left_shared_min, r9.left_shared_max) == (14_585_294, 14_622_365) and r9.left_breakpoint == 14_603_829 and (r9.right_shared_min, r9.right_shared_max) == (32_763_232, 32_763_232)
        and r9.right_breakpoint == 32_763_232 and (r9.left_min, r9.left_max, r9.right_min, r9.right_max) == (14_511_484, 14_977_232, 32_390_715, 32_878_069) and r9.n_calls == 3, str(r9.to_dict()))
    split = pd.DataFrame([dict(assembly='A1', method='A', left_lo=1_000_000, left_hi=1_050_000, right_lo=9_000_000, right_hi=9_050_000),
                          dict(assembly='A2', method='A', left_lo=1_080_000, left_hi=1_120_000, right_lo=9_000_000, right_hi=9_050_000)])           # within gap but NOT overlapping on the left
    rs = paf.consensus(split).iloc[0]
    _ok('calls that are close but whose intervals do not intersect fall back to the median midpoint and report no shared interval (the assemblies disagree)',
        rs.n_calls == 2 and np.isnan(rs.left_shared_min) and rs.left_breakpoint == 1_062_500 and rs.right_breakpoint == 9_025_000 and (rs.right_shared_min, rs.right_shared_max) == (9_000_000, 9_050_000), str(rs.to_dict()))
    one = paf.consensus(real9.iloc[[0]]).iloc[0]
    _ok('a single call is unchanged: its own midpoint, no shared interval', one.left_breakpoint == 14_566_924 and np.isnan(one.left_shared_min))
    # the fLetChi2 failure: a stray forward block 22 Mb away stretched the RIGHT flank to 15.0 Mb, so the big inversion was dropped and a small 2.6 Mb piece was called instead
    TQ, LQ = 'NC_T9b', 38 * M
    lg9q = [line('c', 38 * M, 0.35 * M, 9.70 * M, '+', TQ, LQ, 0.35 * M, 9.70 * M), line('c', 38 * M, 9.72 * M, 11.0 * M, '-', TQ, LQ, 11.0 * M, 12.3 * M),
            line('c', 38 * M, 11.0 * M, 13.63 * M, '+', TQ, LQ, 12.35 * M, 14.98 * M), line('c', 38 * M, 13.66 * M, 31.73 * M, '-', TQ, LQ, 14.62 * M, 32.69 * M),
            line('c', 38 * M, 31.76 * M, 31.90 * M, '+', TQ, LQ, 32.8 * M, 32.94 * M), line('c', 38 * M, 31.9 * M, 37.7 * M, '+', TQ, LQ, 32.94 * M, 38.0 * M),
            line('c', 38 * M, 37.70 * M, 37.78 * M, '+', TQ, LQ, 15.0 * M, 15.08 * M)]                           # the stray forward block (ref 15.0 Mb, near the contig end)
    d9 = paf.parse_paf(lg9q)
    r9q = paf.call_method_AC(d9, TQ)
    big = r9q[(r9q.right_hi.astype(float) - r9q.left_lo.astype(float)) > 10 * M] if len(r9q) else r9q
    _ok('fLetChi2-like data: the big 18 Mb inversion is called despite the stray forward block in the right flank', len(big) == 1 and abs(big.left_lo.iloc[0] - 14.62 * M) < 50_000 and abs(big.right_lo.iloc[0] - 32.69 * M) < 50_000, str(r9q[['method', 'left_lo', 'left_hi', 'right_lo', 'right_hi']].values.tolist()))

    # ---------------- the comparison step on a reference WITH an unplaced scaffold (lg = NaN in the calls and in the scan table): this crashed `centromeres arrays` on the real genome
    calls_t = pd.DataFrame([dict(sequence='NC_135168.1', lg='LG2', start=27_130_000, end=27_600_000), dict(sequence='NC_135176.1', lg='LG10', start=7_000_000, end=7_420_000),
                            dict(sequence='NW_027490041.1', lg=np.nan, start=np.nan, end=np.nan), dict(sequence='NW_027490042.1', lg=np.nan, start=1_000_000, end=1_400_000)])
    scan_t = pd.DataFrame([dict(lg='LG2', start=26_900_000, end=27_900_000), dict(lg='LG10', start=0, end=300_000), dict(lg=np.nan, start=np.nan, end=np.nan)])             # the scan has a NaN row too
    pub_t = pd.DataFrame([dict(lg='LG2', gt3_start=27_129_591, gt3_end=27_702_624), dict(lg='LG13', gt3_start=np.nan, gt3_end=np.nan, bounded_start=32_490_569, bounded_end=35_759_951)])
    cmp_ = csteps.compare_calls(calls_t, pub_t, scan_t, {'LG2': (19_800_671, 43_549_813), 'LG10': (11_774_878, 29_891_917)}).set_index('sequence')
    _ok('compare_calls survives unplaced scaffolds (NaN lg) in the calls, the scan table and the published table',
        cmp_.loc['NC_135168.1', 'array_vs_published'] == 'overlaps' and cmp_.loc['NC_135168.1', 'scan_start'] == 26_900_000 and cmp_.loc['NC_135168.1', 'relation_to_inversion'] == 'pericentric'
        and cmp_.loc['NC_135176.1', 'relation_to_inversion'] == 'paracentric' and cmp_.loc['NC_135176.1', 'scan_end'] == 300_000 and np.isnan(cmp_.loc['NC_135176.1', 'published_start'])
        and cmp_.loc['NW_027490041.1', 'relation_to_inversion'] == '' and np.isnan(cmp_.loc['NW_027490042.1', 'scan_start']), str(cmp_[['array_vs_published', 'scan_start', 'relation_to_inversion']].to_dict('index')))
    _ok('...and with no published table and no scan table at all', len(csteps.compare_calls(calls_t, None, None, {})) == 4)
