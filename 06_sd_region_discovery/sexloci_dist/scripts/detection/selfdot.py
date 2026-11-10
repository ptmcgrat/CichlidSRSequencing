#!/usr/bin/env python3
"""selfdot.py - self dot plot of one genome region, to find tandem / segmental duplications.

Needs only numpy and matplotlib. Example (the gsdf region of the M. zebra GT3a reference):

  python3 selfdot.py GENOME.fna --chrom NC_135173.1 --start 50700000 --end 50830000 \
          --mark 50754902-50776300 --out gsdf_selfdot

GENOME.fna may be plain or .gz. Writes <out>.png (dot plot + identity along the duplicated diagonal) and <out>_segments.tsv,
and prints a short summary: the repeat period, and the start/end of both copies.
How it works: every k-mer that occurs 2..max_occ times in the region is a dot (blue = same strand, red = inverted). A tandem
duplication is a line parallel to the main diagonal, shifted by the repeat period. The period is the most common offset between
matching k-mers; the copies are the stretch along that offset where k-mers still match (identity is estimated from the fraction
of matching k-mers: identity ~ fraction ** (1/k)).
"""
import argparse
import gzip
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_region(path, chrom, start, end):
    """Stream the FASTA and keep only chrom:start-end (1-based, inclusive)."""
    op = gzip.open if path.endswith(".gz") else open
    parts, on, pos = [], False, 0
    with op(path, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                if on:
                    break
                on = line[1:].split()[0] == chrom
                continue
            if on:
                line = line.rstrip("\n")
                n = len(line)
                lo, hi = max(start - 1 - pos, 0), min(end - pos, n)
                if hi > lo:
                    parts.append(line[lo:hi])
                pos += n
                if pos >= end:
                    break
    if not parts:
        sys.exit(f"{chrom}:{start}-{end} not found in {path} (check the chromosome name: grep '>' {path} | head)")
    return "".join(parts).upper()


def kmer_codes(seq, k):
    lut = np.full(256, 255, np.uint8)
    for c, v in zip("ACGT", range(4)):
        lut[ord(c)] = v
    a = lut[np.frombuffer(seq.encode(), np.uint8)]
    n = len(a) - k + 1
    ok = np.ones(n, bool)
    fw = np.zeros(n, np.uint64)
    rc = np.zeros(n, np.uint64)
    for t in range(k):
        x = a[t:t + n]
        ok &= x < 4
        xv = np.where(x < 4, x, 0).astype(np.uint64)
        fw = (fw << np.uint64(2)) | xv
        rc |= (np.uint64(3) - xv) << np.uint64(2 * t)
    return fw, rc, ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True, help="1-based")
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--k", type=int, default=17)
    ap.add_argument("--max-occ", type=int, default=4, help="ignore k-mers seen more often than this (repeats / TEs)")
    ap.add_argument("--min-offset", type=int, default=2000, help="smallest repeat period looked for (ignores short tandem repeats)")
    ap.add_argument("--window", type=int, default=500)
    ap.add_argument("--min-frac", type=float, default=0.15, help="fraction of matching k-mers for a window to count as duplicated")
    ap.add_argument("--gap-windows", type=int, default=4, help="bridge gaps of up to this many windows (inserted repeats, indels)")
    ap.add_argument("--mark", help="start-end of a feature to shade, e.g. the gene: 50754902-50776300")
    ap.add_argument("--out", default="selfdot")
    a = ap.parse_args()

    seq = read_region(a.fasta, a.chrom, a.start, a.end)
    k, L = a.k, len(seq)
    print(f"{a.chrom}:{a.start}-{a.end}: {L:,} bp read; N bases {seq.count('N'):,}")
    fw, rc, ok = kmer_codes(seq, k)
    idx = np.flatnonzero(ok)
    order = np.argsort(fw[idx], kind="stable")
    sc, si = fw[idx][order], idx[order]
    groups = np.split(si, np.flatnonzero(np.diff(sc)) + 1)
    table = {}
    direct = []
    for g in groups:
        if 2 <= len(g) <= a.max_occ:
            for x in range(len(g)):
                for y in range(x + 1, len(g)):
                    direct.append((g[x], g[y]))
        if 1 <= len(g) <= a.max_occ:
            table[int(fw[g[0]])] = g
    direct = np.array(direct, dtype=np.int64).reshape(-1, 2)
    inverted = []
    for j in idx:
        g = table.get(int(rc[j]))
        if g is not None:
            for i in g:
                if i < j:
                    inverted.append((i, j))
    inverted = np.array(inverted, dtype=np.int64).reshape(-1, 2)
    print(f"dots: {len(direct):,} same-strand, {len(inverted):,} inverted (k={k}, k-mers seen 2..{a.max_occ} times)")

    # ---- repeat period: the most common offset between same-strand matches
    off = direct[:, 1] - direct[:, 0]
    off = off[off >= a.min_offset]
    if len(off) < 50:
        print("no tandem/segmental duplication found at this k and --min-offset (fewer than 50 matching k-mers)")
        L0, segs = None, []
    else:
        h, e = np.histogram(off, bins=np.arange(a.min_offset, off.max() + 21, 20))
        top = []
        for b in np.argsort(h)[::-1]:
            c = (e[b] + e[b + 1]) / 2
            if h[b] >= 20 and all(abs(c - t[0]) > 300 for t in top):
                top.append((c, int(h[b])))
            if len(top) == 5:
                break
        print("most common offsets between matching k-mers (bp: matching k-mers):", ", ".join(f"{int(c):,}: {n:,}" for c, n in top))
        L0 = int(np.median(off[np.abs(off - top[0][0]) <= 150]))
        # ---- where along the sequence do k-mers still match at that offset (allowing small indels, +-D)?
        D, w = 150, a.window
        nw = L // w + 1
        best = np.zeros(nw)
        bestd = np.zeros(nw, int)
        npos = len(fw)
        for d in range(-D, D + 1):
            o = L0 + d
            if o <= 0 or o >= npos:
                continue
            m = (fw[:npos - o] == fw[o:]) & ok[:npos - o] & ok[o:]
            cnt = np.bincount(np.arange(npos - o) // w, weights=m, minlength=nw)[:nw]
            upd = cnt > best
            best[upd], bestd[upd] = cnt[upd], d
        frac = best / w
        hit = frac >= a.min_frac
        segs, wi = [], 0
        while wi < nw:
            if hit[wi]:
                wj = wi
                while True:
                    nxt = next((x for x in range(wj + 1, min(wj + 2 + a.gap_windows, nw)) if hit[x]), None)
                    if nxt is None:
                        break
                    wj = nxt
                lo, hi = wi * w, min((wj + 1) * w, L)
                off_med = L0 + int(np.median(bestd[wi:wj + 1]))
                ident = float(np.mean(frac[wi:wj + 1][hit[wi:wj + 1]]) ** (1 / k))
                segs.append((lo, hi, off_med, ident))
                wi = wj + 1
            else:
                wi += 1
        segs = [s for s in segs if s[1] - s[0] >= 2000]
        print(f"\nrepeat period ~ {L0:,} bp. Duplicated segments (copy 1 -> copy 2 = copy 1 + offset):")
        with open(a.out + "_segments.tsv", "w") as fh:
            fh.write("chrom\tcopy1_start\tcopy1_end\tcopy2_start\tcopy2_end\tlength\toffset\test_identity\n")
            for lo, hi, o, idt in segs:
                row = (a.chrom, a.start + lo, a.start + hi - 1, a.start + lo + o, a.start + hi - 1 + o, hi - lo, o, round(idt, 3))
                fh.write("\t".join(map(str, row)) + "\n")
                print(f"  copy 1 {row[1]:,}-{row[2]:,}   copy 2 {row[3]:,}-{row[4]:,}   {row[5]:,} bp   offset {o:,}   identity ~{idt:.1%}")
        if not segs:
            print("  none long enough (>= 2 kb); try a smaller --min-frac or k")

    # ---- figure
    fig, ax = plt.subplots(1, 2, figsize=(14, 6.4), gridspec_kw=dict(width_ratios=[1, 1.15]))
    X = lambda v: (a.start + v) / 1e6
    if len(direct):
        s = direct[np.random.default_rng(1).permutation(len(direct))[:150000]]
        ax[0].scatter(X(s[:, 0]), X(s[:, 1]), s=0.6, c="#2166ac", lw=0)
    if len(inverted):
        s = inverted[:60000]
        ax[0].scatter(X(s[:, 0]), X(s[:, 1]), s=0.6, c="#b2182b", lw=0)
    ax[0].plot([X(0), X(L)], [X(0), X(L)], c="0.8", lw=0.6)
    if a.mark:
        m0, m1 = (int(t) for t in a.mark.split("-"))
        for axis in ("x", "y"):
            (ax[0].axvspan if axis == "x" else ax[0].axhspan)(m0 / 1e6, m1 / 1e6, color="gold", alpha=0.3, lw=0)
        ax[1].axvspan(m0 / 1e6, m1 / 1e6, color="gold", alpha=0.3, lw=0, label="marked feature")
    ax[0].set_xlabel(f"{a.chrom} position (Mb)")
    ax[0].set_ylabel("position of the matching k-mer (Mb)")
    ax[0].set_title("Self dot plot (blue: same strand, red: inverted)", fontsize=10)
    ax[0].set_xlim(X(0), X(L))
    ax[0].set_ylim(X(0), X(L))
    if L0 is not None:
        xs = (np.arange(nw) * w + w / 2)
        ax[1].plot(X(xs), np.where(frac > 0, frac, np.nan) ** (1 / k), c="#2166ac", lw=1.2, label=f"identity at offset ~{L0:,} bp")
        ax[1].set_ylim(0, 1.02)
        ax[1].axhline(a.min_frac ** (1 / k), c="0.6", ls=":", lw=0.8)
        for lo, hi, o, idt in segs:
            ax[1].axvspan(X(lo), X(hi), ymax=0.12, color="#2166ac", alpha=0.5, lw=0)
            ax[1].axvspan(X(lo + o), X(hi + o), ymax=0.12, color="#d6604d", alpha=0.5, lw=0)
        ax[1].set_title("Estimated identity between a position and the position one repeat period downstream\n"
                        "(blue bar: copy 1, red bar: copy 2)", fontsize=9)
        ax[1].legend(fontsize=8, frameon=False, loc="upper right")
    ax[1].set_xlabel(f"{a.chrom} position (Mb)")
    ax[1].set_ylabel("estimated identity")
    fig.tight_layout()
    fig.savefig(a.out + ".png", dpi=150)
    print(f"\nwrote {a.out}.png" + (f" and {a.out}_segments.tsv" if L0 is not None else ""))


if __name__ == "__main__":
    main()
