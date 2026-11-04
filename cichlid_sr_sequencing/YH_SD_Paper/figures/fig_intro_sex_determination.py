"""Intro figure: sex-determination systems across the three radiations seeded from an A. calliptera-like ancestor in Lake Malawi.  python fig_intro_sex_determination.py
Edit the DATA block (rows, columns, cell states) to change the content; everything else is layout. Outputs SVG (editable text) + PDF + PNG.
Sources: Munby et al. 2021 preprint (gsdf-dup, chr7-ins, chr19-ins), Behrens et al. 2025 (LG5 ZW/OB inversion), Clark et al. 2019 (B chromosome),
Blumer et al. 2025 and Kumar et al. 2025 (chr9/10/11 inversions). Backbone: Malinsky et al. 2018 model (A. calliptera lineage seeds pelagic, then benthic+utaka, then mbuna)
with the hybrid benthic origin of Blumer et al. 2025 (arrow A); benthic subgroups from Blumer Data S1 and the SampleDatabase Subgroup labels."""
import re, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch, Rectangle, Polygon

plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.family': 'DejaVu Sans', 'font.size': 6.5})

# ---------------------------------------------------------------- DATA (edit here)
# default = intro figure (literature only). `--with-this-study` adds the 'this study' stars and the preliminary rocky-Aulonocara marker.
SHOW_THIS_STUDY = '--with-this-study' in sys.argv
OUT = 'Fig_intro_sex_determination_Malawi' + ('_with_this_study' if SHOW_THIS_STUDY else '')
W, H = 7.1, 6.81                                          # figure size in inches (journal two-column width)
ROWS = ['A. calliptera', 'Mbuna', 'Utaka*', 'True deeps*', 'Rocky Aulonocara', 'Shallow Lethrinops', 'Other shallow benthics', 'Eukambuzi', 'Diplotaxodon', 'Rhamphochromis']
RADIATIONS = [('A. calliptera', 0, 0, '#9a9a9a', 'riverine'), ('mbuna', 1, 1, '#6b8f71', 'mbuna'), ('benthic + utaka', 2, 7, '#b08968', 'benthic + utaka'), ('pelagic', 8, 9, '#7a8fa6', 'pelagic')]
ROW_COLOR = {'Utaka*': '#2e8b57', 'True deeps*': '#2f6db5', 'Rocky Aulonocara': '#7fb0e0', 'Shallow Lethrinops': '#e8664a', 'Other shallow benthics': '#e8664a', 'Eukambuzi': '#e8664a'}
# (key, header, group, color, heterogamety label)
COLS = [('gsdf', 'gsdf\nduplication', 0, '#2a9d8f', 'XY'),
        ('c7ins', 'chr7 TE\ninsertion', 0, '#2a9d8f', 'XY'),
        ('c19ins', 'chr19 TE\ninsertion', 0, '#2a9d8f', 'XY'),
        ('lg5', 'LG5\ninversion', 1, '#8e44ad', 'ZW'),
        ('bchr', 'B chromo-\nsome', 1, '#b07cc6', 'B'),
        ('c9', 'chr9\ninversion', 2, '#e76f51', 'XY'),
        ('c11', 'chr11\ninversion', 2, '#e76f51', 'XY'),
        ('c10', 'chr10\ninversion', 2, '#e76f51', 'XY')]
GROUPS = ['gsdf region and id3\nY alleles (XY)', 'Dominant female\n(W) systems', 'Benthic XY: introgressed Y-like\nhaplotype of an inversion']
# cell states: 'S' = documented sex-determining / sex-linked, 'P' = variant present, sex role not documented; optional text; '*' adds a 'this study' star
CELLS = {
    ('A. calliptera', 'gsdf'): ('S', '85%'), ('A. calliptera', 'c7ins'): ('S', '4%'), ('A. calliptera', 'c19ins'): ('S', '11%'),
    ('Mbuna', 'gsdf'): ('P', '~75%'), ('Mbuna', 'c7ins'): ('P', '~13%'), ('Mbuna', 'lg5'): ('S', ''), ('Mbuna', 'bchr'): ('S', ''),
    ('Utaka*', 'gsdf'): ('P', ''), ('Utaka*', 'c11'): ('S', 'some spp.'),
    ('True deeps*', 'gsdf'): ('P', ''), ('True deeps*', 'c10'): ('S*', ''),
    ('Rocky Aulonocara', 'gsdf'): ('P', ''), ('Rocky Aulonocara', 'c10'): ('Q*', ''),
    ('Shallow Lethrinops', 'gsdf'): ('P', ''),
    ('Other shallow benthics', 'gsdf'): ('P*', ''), ('Other shallow benthics', 'c11'): ('S', ''),
    ('Eukambuzi', 'gsdf'): ('P', ''), ('Eukambuzi', 'c9'): ('S', ''),
}
if not SHOW_THIS_STUDY:
    CELLS = {k: (v[0].rstrip('*'), v[1]) for k, v in CELLS.items() if not v[0].startswith('Q')}
# ---------------------------------------------------------------- layout (inches)
X_TREE0, X_TREE1 = 0.12, 1.20
X_LAB = 1.52
X_COL0, COL_W = 2.88, 0.52
Y_TOP_ROW, ROW_H = 5.36, 0.255
def ry(i): return Y_TOP_ROW - i * ROW_H
def cx(j): return X_COL0 + COL_W * (j + 0.5)

fig = plt.figure(figsize=(W, H)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis('off')
GREY, DK = '#6b6b6b', '#222222'

def T(x, y, s, size=6.5, ha='left', va='center', color=DK, weight='normal', style='normal', rot=0):
    return ax.text(x, y, s, fontsize=size, ha=ha, va=va, color=color, fontweight=weight, fontstyle=style, rotation=rot, linespacing=1.05)

# panel letters
T(0.02, H - 0.07, 'A', 9, weight='bold'); T(0.02, 2.88, 'B', 9, weight='bold')

# ---------------------------------------------------------------- tree: the A. calliptera lineage is the trunk; three radiations split off in the order pelagic -> benthic+utaka -> mbuna
y = {n: ry(i) for i, n in enumerate(ROWS)}
def mid(*ys): return sum(ys) / len(ys)
def hline(x0, x1, yy, color=DK, lw=0.9, ls='-'): ax.plot([x0, x1], [yy, yy], color=color, lw=lw, ls=ls, solid_capstyle='butt')
def vline(xx, y0, y1, color=DK, lw=0.9): ax.plot([xx, xx], [y0, y1], color=color, lw=lw, solid_capstyle='butt')
L = X_TREE1
X_PEL, X_PELN, X_BEN, X_MB = 0.22, 0.50, 0.66, 0.92                  # split positions along the trunk (schematic, not to scale)
y0 = y['A. calliptera']
yPEL = mid(y['Diplotaxodon'], y['Rhamphochromis'])
yN8 = mid(y['Other shallow benthics'], y['Eukambuzi']); yN7 = mid(y['Shallow Lethrinops'], yN8); yDeep = mid(y['True deeps*'], y['Rocky Aulonocara']); yN6 = mid(yDeep, yN7); yN5 = mid(y['Utaka*'], yN6)
xN5, xN6, xN7, xN8, xDeep = 0.86, 0.96, 1.04, 1.12, 1.08
# trunk (thick) from the riverine-generalist ancestor to extant A. calliptera
hline(X_TREE0, L, y0, lw=2.0, color='#555555')
T(X_TREE0, y0 + 0.17, 'Three radiations seeded from an\nA. calliptera-like riverine generalist', 5.3, color=GREY, va='bottom')
# 1: pelagic radiation (seeded first)
vline(X_PEL, y0, yPEL); hline(X_PEL, X_PELN, yPEL); vline(X_PELN, y['Diplotaxodon'], y['Rhamphochromis'])
hline(X_PELN, L, y['Diplotaxodon']); hline(X_PELN, L, y['Rhamphochromis'])
# 2: benthic + utaka radiation (hybrid origin)
vline(X_BEN, y0, yN5); hline(X_BEN, xN5, yN5); vline(xN5, y['Utaka*'], yN6)
hline(xN5, L, y['Utaka*'], ls=(0, (2, 1.4))); hline(xN5, xN6, yN6)
vline(xN6, yDeep, yN7); hline(xN6, xDeep, yDeep); hline(xN6, xN7, yN7)
vline(xDeep, y['True deeps*'], y['Rocky Aulonocara']); hline(xDeep, L, y['True deeps*'], ls=(0, (2, 1.4))); hline(xDeep, L, y['Rocky Aulonocara'])
vline(xN7, y['Shallow Lethrinops'], yN8); hline(xN7, L, y['Shallow Lethrinops']); hline(xN7, xN8, yN8)
vline(xN8, y['Other shallow benthics'], y['Eukambuzi']); hline(xN8, L, y['Other shallow benthics']); hline(xN8, L, y['Eukambuzi'])
# 3: mbuna radiation (seeded last)
vline(X_MB, y0, y['Mbuna']); hline(X_MB, L, y['Mbuna'])
# seeding order markers 1-3
for n, xx in enumerate([X_PEL, X_BEN, X_MB], start=1):
    ax.add_patch(Circle((xx, y0), 0.052, fc='white', ec='#555555', lw=0.8, zorder=6)); T(xx, y0 - 0.002, str(n), 5.0, ha='center', weight='bold', color='#333333').set_zorder(7)
# A: founding admixture of the benthics (Diplotaxodon-like lineage -> benthic stem); brought the chr9, 11, 20 inversions
ax.annotate('', xy=(X_BEN - 0.012, yN5), xytext=(0.575, y['Diplotaxodon']),
            arrowprops=dict(arrowstyle='-|>', connectionstyle='angle,angleA=90,angleB=0,rad=0', ls=(0, (3, 1.6)), color='#8c4a2f', lw=0.9, shrinkA=0, shrinkB=0, mutation_scale=7))
T(0.54, mid(y['Other shallow benthics'], y['Eukambuzi']) + 0.02, 'A', 6.2, ha='center', weight='bold', color='#8c4a2f')
# radiation bars + labels (between tree and taxon names)
for name, i0, i1, col, lab in RADIATIONS:
    ya, yb = ry(i0) + ROW_H * 0.42, ry(i1) - ROW_H * 0.42
    ax.add_patch(Rectangle((1.255, yb), 0.035, ya - yb, fc=col, ec='none'))
    if i0 != 0:
        T(1.355, (ya + yb) / 2, lab, 4.7 if i0 == i1 else 5.3, ha='center', color=col, weight='bold', rot=90)
# taxon labels + ecogroup colour chips (Blumer palette: shallow red, deep blue, utaka green)
for i, r in enumerate(ROWS):
    if r in ROW_COLOR: ax.add_patch(Rectangle((X_LAB - 0.085, ry(i) - 0.045), 0.04, 0.09, fc=ROW_COLOR[r], ec='none'))
    italic = r in ('Rhamphochromis', 'Diplotaxodon', 'A. calliptera')
    T(X_LAB, ry(i), r, 6.6, style='italic' if italic else 'normal')
# row guides
for i in range(len(ROWS)):
    ax.add_patch(Rectangle((X_COL0 - 0.04, ry(i) - ROW_H / 2), COL_W * len(COLS) + 0.08, ROW_H, fc='#f4f4f4' if i % 2 == 0 else 'white', ec='none', zorder=0))

# ---------------------------------------------------------------- column headers + group bars
gx = {}
for j, (k, hdr, g, col, het) in enumerate(COLS):
    gx.setdefault(g, []).append(j)
    T(cx(j), Y_TOP_ROW + 0.42, hdr, 5.9, ha='center', color=DK)
    # chromosome pictogram above the header
    px, py = cx(j), Y_TOP_ROW + 0.82
    ax.add_patch(FancyBboxPatch((px - 0.045, py - 0.17), 0.09, 0.34, boxstyle='round,pad=0,rounding_size=0.045', fc='white', ec=GREY, lw=0.6))
    ax.add_patch(Circle((px, py + 0.01), 0.018, fc=GREY, ec='none'))
    if k in ('lg5', 'c9', 'c11', 'c10'):
        ax.add_patch(Rectangle((px - 0.045, py - 0.14), 0.09, 0.13 if k == 'lg5' else 0.20, fc=col, ec='none', alpha=0.9))
    elif k == 'bchr':
        ax.add_patch(FancyBboxPatch((px - 0.025, py - 0.10), 0.05, 0.20, boxstyle='round,pad=0,rounding_size=0.025', fc=col, ec='none'))
    else:
        ax.add_patch(Rectangle((px - 0.045, py + 0.07), 0.09, 0.05, fc=col, ec='none'))
    T(px + 0.075, py - 0.10, {'c9': 'B', 'c11': 'C', 'c10': 'E'}.get(k, ''), 5.2, color=GREY, ha='left')
for g, js in gx.items():
    x0, x1 = cx(js[0]) - COL_W * 0.42, cx(js[-1]) + COL_W * 0.42
    ax.plot([x0, x1], [Y_TOP_ROW + 1.14, Y_TOP_ROW + 1.14], color=DK, lw=0.8)
    T((x0 + x1) / 2, Y_TOP_ROW + 1.25, GROUPS[g], 5.8, ha='center', weight='bold')

# ---------------------------------------------------------------- cells
for (r, k), (st, txt) in CELLS.items():
    i, j = ROWS.index(r), [c[0] for c in COLS].index(k); col, het = COLS[j][3], COLS[j][4]
    x0, y0 = cx(j), ry(i)
    solid = st.startswith('S'); prelim = st.startswith('Q')
    ax.add_patch(Circle((x0, y0), 0.083, fc=col, ec=col, lw=0.9, alpha=1.0 if solid else (0.45 if prelim else 0.22), zorder=3))
    if not solid: ax.add_patch(Circle((x0, y0), 0.083, fc='none', ec=col, lw=0.9, ls=(0, (1.6, 1.2)) if prelim else '-', zorder=4))
    if solid: T(x0, y0 - 0.002, het, 4.6, ha='center', color='white', weight='bold').set_zorder(5)
    if prelim: T(x0, y0 - 0.002, het + '?', 4.6, ha='center', color='white', weight='bold').set_zorder(5)
    if txt: T(x0 + 0.10, y0, txt, 4.9, color=GREY)
    if st.endswith('*'): T(x0 - 0.115, y0 + 0.02, '\u2605', 7, ha='center', color='#d4a017').set_zorder(6)

# introgression labels (Blumer events B, C, E) next to the benthic inversion columns
# ---------------------------------------------------------------- panel B: haplotype polarity of the four inversion-linked systems
panels = [('chr9 \u00b7 eukambuzi', 'XY', 'X', 'Y-like', '#e76f51', 'Y-like haplotype introgressed from\noutside the radiation (event B)', False),
          ('chr11 \u00b7 other shallow, utaka', 'XY', 'X', 'Y-like', '#e76f51', 'Y-like haplotype re-introgressed from an\nA. calliptera-like lineage (event C)', False),
          ('chr10 \u00b7 true deeps', 'XY', 'X', 'Y-like', '#e76f51', 'Y-like haplotype introgressed (event E);\nXY with inverted X (Kumar et al. 2025)' if not SHOW_THIS_STUDY else 'Y-like haplotype introgressed (event E);\ncausal mutation identified in this study', SHOW_THIS_STUDY),
          ('LG5 \u00b7 mbuna (orange blotch)', 'ZW', 'Z', 'W', '#8e44ad', 'W haplotype carries the inversion, a\ndominant female locus and OB (pax7a)', False)]
PB_X0, PB_W = 0.12, 1.74
for n, (title, het, a, b, col, note, star) in enumerate(panels):
    x0 = PB_X0 + n * (PB_W + 0.0)
    T(x0 + 0.02, 2.80, title, 6.2, weight='bold')
    if star: T(x0 + PB_W - 0.12, 2.80, '\u2605', 8, ha='center', color='#d4a017')
    inv_first = (het == 'XY')                                   # XY: X inverted, Y-like non-inverted; ZW: Z non-inverted, W inverted
    for m, (lab, inverted) in enumerate([(a, het == 'XY'), (b, het == 'ZW')]):
        px = x0 + 0.38 + m * 0.62; ybot, ytop = 1.62, 2.50
        ax.add_patch(FancyBboxPatch((px - 0.07, ybot), 0.14, ytop - ybot, boxstyle='round,pad=0,rounding_size=0.07', fc='white', ec=GREY, lw=0.7))
        ax.add_patch(Circle((px, ytop - 0.30), 0.028, fc=GREY, ec='none'))
        blk = Rectangle((px - 0.07, ybot + 0.18), 0.14, 0.40, fc=col if inverted else '#cfd8e3', ec='none', alpha=0.95 if inverted else 1)
        ax.add_patch(blk)
        if inverted:
            ax.add_patch(Polygon([[px - 0.04, ybot + 0.30], [px + 0.04, ybot + 0.30], [px, ybot + 0.46]], closed=True, fc='white', ec='none'))
        else:
            ax.add_patch(Polygon([[px - 0.04, ybot + 0.46], [px + 0.04, ybot + 0.46], [px, ybot + 0.30]], closed=True, fc='#7a8ba3', ec='none'))
        T(px, ybot - 0.10, lab, 6.4, ha='center', weight='bold')
        T(px, ybot - 0.23, 'inverted' if inverted else 'non-inverted', 5.2, ha='center', color=GREY)
    T(x0 + 0.38 + 0.31, 2.56, ('\u2642' if het == 'XY' else '\u2640') + ' ' + het, 6.0, ha='center', color=GREY)
    T(x0 + 0.02, 1.14, note, 5.3, va='top', color=DK)
# ---------------------------------------------------------------- legend + sources
ly = 0.72
ax.add_patch(Circle((0.20, ly), 0.07, fc='#555555', ec='#555555')); T(0.32, ly, 'documented sex-determining system', 5.6)
ax.add_patch(Circle((1.95, ly), 0.07, fc='#555555', ec='#555555', alpha=0.22)); ax.add_patch(Circle((1.95, ly), 0.07, fc='none', ec='#555555', lw=0.9))
T(2.07, ly, 'variant present, sex role not documented', 5.6)
if SHOW_THIS_STUDY:
    ax.add_patch(Circle((4.10, ly), 0.07, fc='#555555', ec='#555555', alpha=0.45)); ax.add_patch(Circle((4.10, ly), 0.07, fc='none', ec='#555555', lw=0.9, ls=(0, (1.6, 1.2))))
    T(4.22, ly, 'preliminary (this study)', 5.6)
    T(5.42, ly, '\u2605', 8, ha='center', color='#d4a017'); T(5.52, ly, 'analysed in this study', 5.6)
FOOT = ['1-3: order in which the three radiations were seeded from the A. calliptera lineage (Malinsky et al. 2018). A: the benthic clade arose by admixture of',
        'Diplotaxodon- and A. calliptera-like lineages, which brought the chr9, 11 and 20 inversions (Blumer et al. 2025). *Utaka and true deeps are not monophyletic (rocky Aulonocara are nested within).',
        'Percentages: A. calliptera = share of Masoko males carrying the allele; Mbuna = share of species (Munby et al. 2021). Schematic, not to scale.',
        'Sources: Munby et al. 2021; Behrens et al. 2025; Clark et al. 2019; Blumer et al. 2025 (events B, C, E in panel B); Kumar et al. 2025; Malinsky et al. 2018.']
for k, line in enumerate(FOOT):
    T(0.12, 0.50 - 0.105 * k, line, 5.0, color=GREY)

for ext in ('svg', 'pdf', 'png'):
    fig.savefig(f'{OUT}.{ext}', dpi=300)
# make the SVG text editable in Illustrator/Inkscape with a common font stack
s = open(f'{OUT}.svg', encoding='utf-8').read()
s = re.sub(r"font:\s*([^;\"']*?)DejaVu Sans", r"font: \1Arial, Helvetica, sans-serif", s).replace('DejaVu Sans', 'Arial, Helvetica, sans-serif')
open(f'{OUT}.svg', 'w', encoding='utf-8').write(s)
print(f'wrote {OUT}.svg/.pdf/.png')
