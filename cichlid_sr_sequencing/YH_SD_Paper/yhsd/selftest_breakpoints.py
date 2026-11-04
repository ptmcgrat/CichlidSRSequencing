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
    steps.call_step(proj, align_dir=adir)
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
    c10 = [line('c10', 36.2 * M, 0, 6.8 * M, '+', TC, LC, 0, 6.8 * M), line('c10', 36.2 * M, 6.9 * M, 18 * M, '-', TC, LC, 18.2 * M, 29.3 * M),
           line('c10', 36.2 * M, 18.4 * M, 25.3 * M, '-', TC, LC, 11.4 * M, 18.2 * M), line('c10', 36.2 * M, 25.4 * M, 36.2 * M, '+', TC, LC, 29.4 * M, 36.0 * M)]
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
