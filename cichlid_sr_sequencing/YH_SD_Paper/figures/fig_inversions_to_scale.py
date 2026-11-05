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

# values from the real runs (call 20261104-193337, lift 20261104-200249, flank landings of LG13 from the v22 lift): replace with --overview when you rerun
EMBEDDED = pd.DataFrame([
    dict(lg='LG2',  length=48_807_247, inv_start=19_800_671, inv_end=43_549_813, cen_start=27_130_000, cen_end=27_600_000, cen_size=454_000, cen_confidence='array'),
    dict(lg='LG5',  length=42_954_235, inv_start=19_657_805, inv_end=26_726_680, cen_start=3_372_816,  cen_end=5_686_969,  cen_size=float('nan'), cen_confidence='published'),
    dict(lg='LG10', length=34_949_860, inv_start=11_774_878, inv_end=29_891_917, cen_start=7_000_000,  cen_end=7_420_000,  cen_size=394_000, cen_confidence='array'),
    dict(lg='LG13', length=36_767_363, inv_start=2_405_431,  inv_end=23_099_416, cen_start=32_490_000, cen_end=32_920_000, cen_size=419_000, cen_confidence='array')])
# centromeres of LG2, LG10, LG13: largest rnd-4_family-47 satellite array > 500 kb from either chromosome end (array-level call, GT3 coordinates);
# LG5: published M. zebra assembly lifted to GT3 (no array-level call yet)
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
H = TOP + ROW * len(ov) + 0.86
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
    for k, (sa, sb) in enumerate(sorted(segs)):                                   # one block per inverted segment; a tandem (double) inversion is drawn as two adjacent shades
        ax.add_patch(Rectangle((x(sa), yc - bh / 2), mb(sb - sa) * SCALE, bh, fc=col, ec='white' if len(segs) > 1 else 'none', lw=0.8, alpha=0.92 if k % 2 == 0 else 0.62, zorder=2))
    # centromere: dark bar across the chromosome (minimum visible width); lighter dashed bar when it is only bounded by flanks
    cw = max(mb(c1 - c0) * SCALE, 0.055) if has_cen else 0
    if not has_cen:
        pass
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
    T(0.12, yc - 0.075, (rel + ' ' if rel else '') + (f'double inversion ' if len(segs) == 2 else 'inversion ') + f'{mb(s):.2f}-{mb(e):.2f} Mb', 5.3, color=GREY)
    csz = f', {r.cen_size / 1e3:.0f} kb array' if pd.notna(r.cen_size) else ''
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
NOTES = ['Coordinates: Mzebra_GT3_NCBI (GCF_041146795.1); all chromosomes share one scale. Bracket: inverted segment with its size and breakpoints (Mb).',
         'Sex-linked systems (literature): XY on LG9 (eukambuzi), LG10 (deep benthics) and LG11 (other shallow benthics, utaka); ZW on LG5 (mbuna, orange blotch).',
         'Breakpoints: LG2, LG9, LG10, LG11 and LG13 from the YH assembly (plus paper assemblies for LG9 and LG11); LG5 from two orange-blotch assemblies.',
         'Centromeres: largest rnd-4_family-47 satellite array more than 500 kb from either chromosome end; LG5: published M. zebra assembly lifted to GT3.']
for k, line in enumerate(NOTES):
    T(0.12, 0.54 - 0.105 * k, line, 5.0, color=GREY)
for ext in ('svg', 'pdf', 'png'):
    fig.savefig(f'{a.out}.{ext}', dpi=300)
s = open(f'{a.out}.svg', encoding='utf-8').read()
s = re.sub(r"font:\s*([^;\"']*?)DejaVu Sans", r"font: \1Arial, Helvetica, sans-serif", s).replace('DejaVu Sans', 'Arial, Helvetica, sans-serif')
open(f'{a.out}.svg', 'w', encoding='utf-8').write(s)
print(f'wrote {a.out}.svg/.pdf/.png')
