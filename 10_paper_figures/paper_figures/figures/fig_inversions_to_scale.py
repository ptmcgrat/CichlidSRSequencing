"""To-scale chromosome diagram of the four inversions (LG2, LG5, LG10, LG13) with their centromeres, in Mzebra_GT3_NCBI coordinates.
    python fig_inversions_to_scale.py                      # uses the values embedded below (2026-11-04 runs)
    python fig_inversions_to_scale.py --overview <inversion_overview.tsv>     # or the table written by `inversion_breakpoints figure`
    python fig_inversions_to_scale.py --short-arm-left     # mirror a chromosome so its centromere is in the left half (display only)
Outputs SVG (live text, Arial), PDF and PNG. Nothing is hatched or masked: LG2 is drawn as the reference has it."""
import argparse, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle, FancyBboxPatch

plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.family': 'DejaVu Sans', 'font.size': 6.5})

# values of the run of 2026-11-05 (inversion_overview.tsv written by `inversion_breakpoints figure`, YH_SD_Paper v29): replace with --overview when you rerun
EMBEDDED = pd.DataFrame([
    dict(lg='LG2',  length=48_807_247, inv_start=19_800_671, inv_end=43_549_813, inv_segments='20370089-43478612', cen_start=27_129_591, cen_end=27_603_422, cen_size=454_843, cen_confidence='array'),
    dict(lg='LG5',  length=42_954_235, inv_start=19_657_805, inv_end=26_726_680, inv_segments='19657762-26726685', cen_start=3_348_089, cen_end=5_686_970, cen_size=2_275_986, cen_confidence='array'),
    dict(lg='LG9',  length=37_905_347, inv_start=14_621_505, inv_end=32_763_232, inv_segments='15313086-32161792', cen_start=6_069_387, cen_end=6_649_889, cen_size=549_041, cen_confidence='array'),
    dict(lg='LG10', length=34_949_860, inv_start=11_774_878, inv_end=29_891_917, inv_segments='11845191-29902559', cen_start=6_916_061, cen_end=7_422_555, cen_size=403_755, cen_confidence='array'),
    dict(lg='LG11', length=38_663_802, inv_start=8_293_847,  inv_end=30_388_903, inv_segments='8303032-17863061;18502589-30259182', cen_start=5_038_192, cen_end=5_625_078, cen_size=516_656, cen_confidence='array'),
    dict(lg='LG13', length=36_767_363, inv_start=2_405_431,  inv_end=23_099_416, inv_segments='2484781-22835665', cen_start=32_490_570, cen_end=32_918_897, cen_size=419_679, cen_confidence='array-ambiguous')])
# centromeres: largest rnd-4_family-47 satellite array more than 500 kb from either chromosome end; LG13 is "ambiguous" (a second internal array of comparable size)
# colour = what the literature says about the inversion: sex-linked XY (Blumer et al. 2025, Kumar et al. 2025), sex-linked ZW (Behrens et al. 2025), not sex-linked
STYLE = {'LG2': ('#6b7f99', 'LG2', 'YH'), 'LG5': ('#8e44ad', 'LG5  (ZW, orange blotch)', 'L. trewavasae + M. zebra OB'),
         'LG9': ('#e76f51', 'LG9  (XY)', 'YH + paper assemblies'), 'LG10': ('#e76f51', 'LG10  (XY)', 'YH'),
         'LG11': ('#e76f51', 'LG11  (XY)', 'YH + paper assemblies'), 'LG13': ('#6b7f99', 'LG13', 'YH')}
ORDER = ['LG2', 'LG5', 'LG9', 'LG10', 'LG11', 'LG13']      # rows are drawn for the chromosomes present in the table (embedded values: LG2, LG5, LG10, LG13; pass --overview for LG9 and LG11)

ap = argparse.ArgumentParser()
ap.add_argument('--overview', default=None); ap.add_argument('--short-arm-left', action='store_true'); ap.add_argument('--out', default='Fig_inversions_to_scale')
a = ap.parse_args()
ov = pd.read_csv(a.overview, sep='\t') if a.overview else EMBEDDED
ov = ov.set_index('lg'); ov = ov.loc[[g for g in ORDER if g in ov.index]].reset_index()
for c_ in ('inv_segments',):
    if c_ not in ov: ov[c_] = ''
for c_ in ('cen_start', 'cen_end', 'cen_confidence', 'cen_size'):
    if c_ not in ov: ov[c_] = float('nan') if c_ != 'cen_confidence' else ''

def relation(c0, c1, i0, i1):
    if pd.isna(c0) or pd.isna(c1): return ''
    return 'paracentric' if (c1 < i0 or c0 > i1) else ('pericentric' if (c0 >= i0 and c1 <= i1) else 'breakpoint in centromere')


W, X0, SCALE = 7.1, 1.95, 0.088                          # figure width, left edge of the chromosomes (in), inches per Mb
ROW, TOP = 0.78, 0.55
H = TOP + ROW * len(ov) + 0.90
fig = plt.figure(figsize=(W, H)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis('off')
GREY, DK = '#6b6b6b', '#222222'
T = lambda x, y, s, size=6.5, ha='left', va='center', color=DK, weight='normal', style='normal': ax.text(x, y, s, fontsize=size, ha=ha, va=va, color=color, fontweight=weight, fontstyle=style)
mb = lambda v: v / 1e6

for i, r in ov.iterrows():
    yc = H - TOP - ROW * i - 0.12
    col, title, src = STYLE[r.lg]
    L, s, e, c0, c1 = r.length, r.inv_start, r.inv_end, r.cen_start, r.cen_end
    rel = r.get('cen_relation', '') if isinstance(r.get('cen_relation', ''), str) and r.get('cen_relation', '') else relation(c0, c1, s, e)
    has_cen = pd.notna(c0) and pd.notna(c1)
    if a.short_arm_left and has_cen and (c0 + c1) / 2 > L / 2:           # mirror for display only
        s, e, c0, c1 = L - e, L - s, L - c1, L - c0
    x = lambda v: X0 + mb(v) * SCALE
    bh = 0.15
    ax.add_patch(FancyBboxPatch((x(0), yc - bh / 2), mb(L) * SCALE, bh, boxstyle='round,pad=0,rounding_size=0.07', fc='#efefef', ec='#7a7a7a', lw=0.7, zorder=1))
    segs = [tuple(int(v) for v in sg.split('-')) for sg in str(r.inv_segments).split(';') if '-' in str(sg)] if isinstance(r.inv_segments, str) and r.inv_segments else []
    if a.short_arm_left and has_cen and (c0 + c1) / 2 > L / 2:
        segs = [(L - b_, L - a_) for a_, b_ in segs]
    if len(segs) < 2:
        segs = [(s, e)]
    segs = sorted(segs)
    cuts = [s] + [(segs[k][1] + segs[k + 1][0]) / 2 for k in range(len(segs) - 1)] + [e]       # the two halves of a double inversion meet at the internal breakpoint
    ib = cuts[1:-1]
    for k, (sa, sb) in enumerate(zip(cuts[:-1], cuts[1:])):
        ax.add_patch(Rectangle((x(sa), yc - bh / 2), mb(sb - sa) * SCALE, bh, fc=col, ec='white' if len(cuts) > 2 else 'none', lw=0.8, alpha=0.92 if k % 2 == 0 else 0.62, zorder=2))
    for v in ib:
        ax.plot([x(v), x(v)], [yc - bh / 2 - 0.035, yc - bh / 2], color=DK, lw=0.7, zorder=6); T(x(v), yc - bh / 2 - 0.085, f'{mb(v):.2f}', 5.0, ha='center', color=GREY)

    # centromere: dark bar across the chromosome (minimum visible width); lighter dashed bar when it is only bounded by flanks
    cw = max(mb(c1 - c0) * SCALE, 0.055) if has_cen else 0
    if not has_cen:
        pass
    elif r.cen_confidence == 'array-ambiguous':
        ax.add_patch(Rectangle(((x(c0) + x(c1)) / 2 - cw / 2, yc - 0.115), cw, 0.23, fc='#222222', ec='none', alpha=0.45, zorder=4))
    elif r.cen_confidence == 'published-bounded':
        ax.add_patch(Rectangle((x(c0), yc - 0.115), mb(c1 - c0) * SCALE, 0.23, fc='#333333', ec='#333333', lw=0.8, ls=(0, (2, 1.5)), alpha=0.35, zorder=4))
    else:
        ax.add_patch(Rectangle(((x(c0) + x(c1)) / 2 - cw / 2, yc - 0.115), cw, 0.23, fc='#222222', ec='none', zorder=4))
    # inversion bracket above the chromosome
    yb = yc + 0.215
    ax.plot([x(s), x(s), x(e), x(e)], [yb - 0.05, yb, yb, yb - 0.05], color=col, lw=1.0, solid_capstyle='butt', zorder=5)
    T((x(s) + x(e)) / 2, yb + 0.09, f'{mb(e - s):.1f} Mb' + (f'  ({len(segs)} segments)' if len(segs) > 1 else ''), 6.8, ha='center', weight='bold', color=col)
    T(x(s) - 0.04, yb + 0.005, f'{mb(s):.2f}', 5.3, ha='right', color=GREY); T(x(e) + 0.04, yb + 0.005, f'{mb(e):.2f}', 5.3, ha='left', color=GREY)
    # row label
    T(0.12, yc + 0.09, title, 7.2, weight='bold')
    T(0.12, yc - 0.075, (rel + ' ' if rel else '') + ('double inversion ' if len(segs) == 2 else 'inversion ') + f'{mb(s):.2f}-{mb(e):.2f} Mb', 5.3, color=GREY)
    if len(ib) == 1:
        T(0.12, yc - 0.305, f'internal breakpoint {mb(ib[0]):.2f} Mb', 5.3, color=GREY)
    csz = f', {r.cen_size / 1e3:.0f} kb' + (', ambiguous' if r.cen_confidence == 'array-ambiguous' else ' array') if pd.notna(r.cen_size) else ''
    T(0.12, yc - 0.19, (f'centromere {mb(c0):.2f}-{mb(c1):.2f} Mb' + csz + (' (bounded)' if r.cen_confidence == 'published-bounded' else '')) if has_cen else 'centromere not located', 5.3, color=GREY)
    T(x(L) + 0.07, yc, f'{mb(L):.1f} Mb', 5.3, color=GREY)
# scale bar
ys = 1.05
ax.plot([X0, X0 + 10 * SCALE], [ys, ys], color=DK, lw=1.2, solid_capstyle='butt'); T(X0 + 5 * SCALE, ys + 0.09, '10 Mb', 6.0, ha='center')
# legend
ly = 0.80
def key(xx, fc, label, dashed=False, w=0.16, alpha=1.0, ec='none'):
    ax.add_patch(Rectangle((xx, ly - 0.06), w, 0.12, fc=fc, ec=('#333333' if dashed else ec), lw=0.8, ls=(0, (2, 1.5)) if dashed else '-', alpha=alpha)); T(xx + w + 0.07, ly, label, 5.8)
key(0.12, '#e76f51', 'sex-linked, XY'); key(1.30, '#8e44ad', 'sex-linked, ZW'); key(2.46, '#6b7f99', 'not sex-linked')
key(3.62, '#222222', 'centromere', w=0.06)
if (ov.cen_confidence == 'published-bounded').any(): key(4.62, '#333333', 'centromere bounded by flanks', dashed=True, alpha=0.35)
if (ov.cen_confidence == 'array-ambiguous').any(): key(4.62, '#222222', 'centromere, ambiguous (second comparable array)', alpha=0.45, w=0.06)
NOTES = ['Coordinates: Mzebra_GT3_NCBI (GCF_041146795.1); all chromosomes share one scale. Bracket: inversion with its size and breakpoints (Mb); a double inversion is drawn as two segments.',
         'Sex-linked systems (literature): XY on LG9 (eukambuzi), LG10 (deep benthics) and LG11 (other shallow benthics, utaka); ZW on LG5 (mbuna, orange blotch).',
         'Breakpoints: LG2, LG10, LG13 from the YH assembly only; LG9 from four and LG11 from three assemblies (shared breakpoint intervals); LG5 from two assemblies.',
         'Centromeres: largest rnd-4_family-47 satellite array more than 500 kb from either chromosome end.']
for k, line in enumerate(NOTES):
    T(0.12, 0.54 - 0.105 * k, line, 5.0, color=GREY)
for ext in ('svg', 'pdf', 'png'):
    fig.savefig(f'{a.out}.{ext}', dpi=300)
s = open(f'{a.out}.svg', encoding='utf-8').read()
s = re.sub(r"font:\s*([^;\"']*?)DejaVu Sans", r"font: \1Arial, Helvetica, sans-serif", s).replace('DejaVu Sans', 'Arial, Helvetica, sans-serif')
open(f'{a.out}.svg', 'w', encoding='utf-8').write(s)
print(f'wrote {a.out}.svg/.pdf/.png')
