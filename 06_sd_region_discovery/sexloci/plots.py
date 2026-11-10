"""Genome-wide and per-region diagnostic plots."""
from __future__ import annotations

import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .knownloci import place_known_loci


def _natkey(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def label_of(chrom, cfg):
    return (getattr(cfg, "contig_labels", None) or {}).get(chrom, chrom)


def chrom_order(chroms, cfg):
    """Linkage groups in numeric order (LG1, LG2, ...), anything unlabelled after them in natural order."""
    def key(c):
        lab = label_of(c, cfg)
        return (0, int(lab[2:])) if lab.startswith("LG") and lab[2:].isdigit() else (1, _natkey(c))
    return sorted(chroms, key=key)


def _offsets(w, window_size, cfg=None):
    chroms = chrom_order(w["chrom"].unique(), cfg) if cfg is not None else sorted(w["chrom"].unique(), key=_natkey)
    off, tot = {}, 0
    ends = w.groupby("chrom")["end"].max()
    for c in chroms:
        off[c] = tot
        tot += int(ends[c]) + window_size
    return chroms, off, tot


def _panels(w):
    p = []
    if "assoc_max_chi2" in w:
        p.append(("assoc_max_chi2", "Association\nmax chi2 / window", "assoc"))
    if "ped_sire_T" in w or "ped_dam_T" in w:
        p.append(("ped", "Pedigree: mean |r|\n(inherited haplotype\nvs sex)", "ped"))
    if "fst" in w:
        p.append(("fst", "Fst\n(male vs female)", "plain"))
        p.append(("dhet", "dHet\n(male - female)", "plain"))
    if "cov_log2_mf" in w:
        p.append(("cov_log2_mf", "log2 depth\nmale / female", "plain"))
    return p


def overlay_known_loci(axes, chroms, off, w, cfg):
    """Dashed vertical lines at the known sex-determination loci on every panel, labelled above the first panel (loci that would overlap
    at genome scale, like the two insertions 216 bp apart, share one label). Returns [(locus, x)] for what was drawn."""
    if not getattr(cfg, "show_known_loci", True):
        return []
    ends = w.groupby("chrom")["end"].max()
    marks = []
    for loc, c in place_known_loci(chroms, lambda k: label_of(k, cfg)):
        if c in off and loc.mid <= ends[c] + cfg.window_size:                 # skip a locus beyond the windows that were analysed
            marks.append((loc, off[c] + loc.mid))
    if not marks:
        return []
    for ax in axes:
        for loc, x in marks:
            ax.axvline(x, color=loc.color, ls="--" if loc.kind == "gene" else ":", lw=1.1, alpha=0.9, zorder=1)
    tot = max(off.values()) + float(ends.max()) + cfg.window_size
    groups, cur = [], []
    for loc, x in sorted(marks, key=lambda m: m[1]):
        if cur and x - cur[-1][1] > 0.02 * tot:
            groups.append(cur)
            cur = []
        cur.append((loc, x))
    groups.append(cur)
    for g in groups:
        axes[0].text(np.mean([x for _, x in g]), 1.02, "/".join(l.name for l, _ in g), transform=axes[0].get_xaxis_transform(),
                     ha="center", va="bottom", fontsize=7.5, fontweight="bold", color=g[0][0].color, clip_on=False)
    return marks


def genome_plot(w, regions, thr_chi2, cfg, path):
    panels = _panels(w)
    if not panels or w.empty:
        return
    chroms, off, tot = _offsets(w, cfg.window_size, cfg)
    fig, axes = plt.subplots(len(panels), 1, figsize=(14, 2.2 * len(panels) + 0.8), sharex=True)
    axes = np.atleast_1d(axes)
    x = (w["chrom"].map(off) + (w["start"] + w["end"]) / 2).to_numpy()
    cc = w["chrom"].map({c: i % 2 for i, c in enumerate(chroms)}).to_numpy()
    cols = np.array(["#4a4a4a", "#9a9a9a"])
    for ax, (col, label, kind) in zip(axes, panels):
        if kind == "ped":
            for role, color in (("sire", "#1f77b4"), ("dam", "#d62728")):
                c = f"ped_{role}_absr"
                if c in w:
                    y = w[c].to_numpy(float)
                    sig = w[f"ped_{role}_padj"].to_numpy(float) <= cfg.alpha
                    ax.scatter(x, y, s=8, c=color, alpha=0.45, label=role)
                    ax.scatter(x[sig], y[sig], s=26, facecolors="none", edgecolors=color, linewidth=0.9)
            ax.set_ylim(0, 1.05)
            ax.legend(loc="lower right", fontsize=7, ncol=2, frameon=False)
            ax.text(0.005, 0.97, "circled = significant after correction", transform=ax.transAxes, fontsize=6.5,
                    va="top", color="0.3")
        else:
            y = w[col].to_numpy(float)
            ax.scatter(x, y, s=8, c=cols[cc])
            if kind == "assoc":
                sig = w["assoc_padj"].to_numpy(float) <= cfg.alpha
                ax.scatter(x[sig], y[sig], s=26, c="#d62728", edgecolor="k", linewidth=0.5)
                if np.isfinite(thr_chi2):
                    ax.axhline(thr_chi2, ls="--", c="#d62728", lw=0.8)
        ax.set_ylabel(label, fontsize=8)
        for _, r in regions.iterrows():
            ax.axvspan(off[r.chrom] + r.start, off[r.chrom] + r.end, color="gold", alpha=0.35, lw=0)
    mids = [off[c] + (w[w.chrom == c]["end"].max()) / 2 for c in chroms]
    axes[-1].set_xticks(mids)
    axes[-1].set_xticklabels([label_of(c, cfg) for c in chroms], rotation=90, fontsize=7)
    marked = overlay_known_loci(axes, chroms, off, w, cfg)
    fig.suptitle("Sex-linkage evidence by window, whole genome (gold = called regions"
                 + ("; dashed/dotted lines = known sex-determination loci)" if marked else ")"), fontsize=10, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    fig.savefig(path, dpi=130)
    fig.savefig(os.path.splitext(path)[0] + ".pdf")
    plt.close(fig)


def region_plot(w, r, thr_chi2, cfg, path, pad_windows=4):
    pad = pad_windows * cfg.window_size
    sub = w[(w.chrom == r.chrom) & (w.end >= r.start - pad) & (w.start <= r.end + pad)]
    panels = _panels(sub)
    if sub.empty or not panels:
        return
    fig, axes = plt.subplots(len(panels), 1, figsize=(9, 1.9 * len(panels) + 0.8), sharex=True)
    axes = np.atleast_1d(axes)
    x = ((sub["start"] + sub["end"]) / 2e6).to_numpy()
    for ax, (col, label, kind) in zip(axes, panels):
        if kind == "ped":
            for role, color in (("sire", "#1f77b4"), ("dam", "#d62728")):
                c = f"ped_{role}_absr"
                if c in sub:
                    ax.plot(x, sub[c].to_numpy(float), "o-", ms=3, c=color, label=role)
            ax.set_ylim(0, 1.05)
            ax.legend(fontsize=7, frameon=False)
        else:
            ax.plot(x, sub[col].to_numpy(float), "o-", ms=3, c="#333")
            if kind == "assoc" and np.isfinite(thr_chi2):
                ax.axhline(thr_chi2, ls="--", c="#d62728", lw=0.8)
        ax.axvspan(r.start / 1e6, r.end / 1e6, color="gold", alpha=0.35, lw=0)
        ax.set_ylabel(label, fontsize=8)
    axes[-1].set_xlabel(f"{label_of(r.chrom, cfg)} ({r.chrom}) position (Mb)")
    fig.suptitle(f"{label_of(r.chrom, cfg)}:{r.start:,}-{r.end:,}  tier {r.tier}  model {r.model}", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------------------
# Inheritance along the genome (the "haplotype blocks" view)
# ----------------------------------------------------------------------------------------------------------------
def align_states(U: np.ndarray) -> np.ndarray:
    """Turn per-window consensus scores (windows x offspring) into inherited-haplotype states (0/1).

    The sign of each window's principal component is arbitrary, so every window is flipped to agree with the
    previous one; the state is then the sign of the aligned score. Orientation is arbitrary per chromosome but
    constant along it, so recombinants appear as a single offspring switching state.
    """
    out = np.empty(U.shape)
    prev = None
    for j in range(len(U)):
        u = U[j]
        if prev is not None:
            sgn = np.sign(np.nansum(u * prev))
            u = u * (1.0 if sgn == 0 else sgn)
        out[j] = u
        prev = u
    return (out > 0).astype(float)


def haplotype_genome_plot(tracks, family_meta, win, regions, cfg, path):
    """One panel per family and parent: rows = offspring (females then males), columns = genome windows, colour =
    which of the parent's two haplotypes was inherited. A sex-linked region shows as a block where the colours split
    by sex."""
    from matplotlib.colors import ListedColormap
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Rectangle
    if not tracks or win.empty:
        return
    chroms, off, tot = _offsets(win, cfg.window_size, cfg)
    ncol = int(tot / cfg.window_size) + 1
    fams = sorted(family_meta)
    rows = [(f, role) for f in fams for role in ("sire", "dam") if (f, role) in tracks]
    if not rows:
        return
    heights = [max(len(family_meta[f]["sex"]), 4) for f, _ in rows]
    fig = plt.figure(figsize=(18, 1.0 + 0.14 * sum(heights) + 0.35 * len(rows)))
    gs = GridSpec(len(rows), 2, width_ratios=[1, 0.008], height_ratios=heights, hspace=0.25, wspace=0.005,
                  top=0.955, bottom=0.07, left=0.075, right=0.985)
    cmap = ListedColormap(["#2166ac", "#f4a582"])
    cmap.set_bad("white")
    sexcol = ListedColormap(["#d6604d", "#4393c3"])           # female red, male blue
    first = None
    for i, (f, role) in enumerate(rows):
        meta = family_meta[f]
        order = np.argsort(meta["sex"], kind="stable")          # females (0) first, then males (1)
        M = np.full((len(order), ncol), np.nan)
        for c, (widx, U) in tracks[(f, role)].items():
            states = align_states(U)
            base = int(off[c] / cfg.window_size)
            M[:, base + widx] = states[:, order].T if states.shape[1] == len(order) else np.nan
        ax = fig.add_subplot(gs[i, 0], sharex=first)
        first = first or ax
        ax.imshow(np.ma.masked_invalid(M), aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
        axs = fig.add_subplot(gs[i, 1])
        axs.imshow(meta["sex"][order][:, None], aspect="auto", cmap=sexcol, vmin=0, vmax=1)
        axs.set_xticks([]); axs.set_yticks([])
        if len(order) <= 40:
            ax.set_yticks(range(len(order)))
            ax.set_yticklabels([meta["names"][k] for k in order], fontsize=5)
        else:
            ax.set_yticks([])
        ax.set_ylabel(f"{f.replace('fam:', '')}\n{role} haplotype", fontsize=6)
        for c in chroms:
            ax.axvline(off[c] / cfg.window_size - 0.5, color="k", lw=0.4)
        for _, r in regions.iterrows():
            x0 = off[r.chrom] / cfg.window_size + (r.start - 1) / cfg.window_size
            ax.add_patch(Rectangle((x0 - 0.5, -0.5), (r.end - r.start + 1) / cfg.window_size, len(order),
                                   fill=False, ec="gold", lw=1.6))
        ax.tick_params(axis="x", labelbottom=False)
    mids = [off[c] / cfg.window_size + win[win.chrom == c].widx.max() / 2 for c in chroms]
    ax.set_xticks(mids)
    ax.set_xticklabels([label_of(c, cfg) for c in chroms], rotation=90, fontsize=7)
    ax.tick_params(axis="x", labelbottom=True)
    fig.suptitle("Inherited parental haplotype along the genome (rows: offspring, females red / males blue in the bar at right; "
                 "the two colours are the parent's two haplotypes; gold = called regions)", fontsize=9)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    fig.savefig(os.path.splitext(path)[0] + ".pdf", bbox_inches="tight")
    plt.close(fig)


def chromosome_summary_plot(fam_tables, family_meta, cfg, path):
    """Per chromosome and family: the strongest association between inherited parental haplotype and offspring sex
    (|r| over windows), sire vs dam. A sex chromosome stands out for one parent in every family."""
    if not fam_tables:
        return
    rows = []
    for f, t in fam_tables.items():
        n = len(family_meta[f]["sex"])
        for (role, chrom), g in t.groupby(["role", "chrom"]):
            rows.append(dict(family=f, role=role, chrom=chrom, r=float(np.sqrt(g["T"].max() / n))))
    df = pd.DataFrame(rows)
    chroms = chrom_order(df.chrom.unique(), cfg)
    x = {c: i for i, c in enumerate(chroms)}
    fams = sorted(df.family.unique())
    marks = ["o", "s", "^", "D", "v", "P", "X"]
    fig, ax = plt.subplots(figsize=(14, 4))
    for k, f in enumerate(fams):
        for role, color in (("sire", "#1f77b4"), ("dam", "#d62728")):
            d = df[(df.family == f) & (df.role == role)]
            ax.scatter([x[c] + (-0.15 if role == "sire" else 0.15) for c in d.chrom], d.r, c=color,
                       marker=marks[k % len(marks)], s=34, alpha=0.8, edgecolor="k", linewidth=0.4,
                       label=f"{f.replace('fam:', '')} {role}")
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels([label_of(c, cfg) for c in chroms], rotation=90, fontsize=8)
    ax.set_ylabel("max |r| (inherited haplotype vs sex)")
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=6, ncol=2, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.set_title("Strongest haplotype-sex association per chromosome (blue = sire, red = dam; marker shape = family)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    fig.savefig(os.path.splitext(path)[0] + ".pdf")
    plt.close(fig)


def retention_plot(by_window, cfg, path, window_size=None):
    """Whole-genome view of what an existing filtered VCF kept: per-window retention, depleted windows highlighted."""
    if by_window is None or by_window.empty:
        return
    ws = window_size or cfg.window_size
    w = by_window.assign(start=by_window["window"] * ws, end=(by_window["window"] + 1) * ws)
    chroms, off, tot = _offsets(w.rename(columns={"window": "widx"}), ws, cfg)
    x = (w.chrom.map(off) + (w.start + w.end) / 2).to_numpy()
    fig, ax = plt.subplots(figsize=(14, 3.6))
    ax.scatter(x, w.retained_filtered_file, s=6, c="#4a4a4a", label="kept by the filtered file")
    dep = w.depleted_in_filtered_file.to_numpy()
    ax.scatter(x[dep], w.retained_filtered_file.to_numpy()[dep], s=22, c="#d62728", label="depleted window")
    ax.plot(x, w.retained_ours, lw=0.5, c="#1f77b4", alpha=0.5, label="kept by sexloci filter")
    mids = [off[c] + w[w.chrom == c].end.max() / 2 for c in chroms]
    ax.set_xticks(mids)
    ax.set_xticklabels([label_of(c, cfg) for c in chroms], rotation=90, fontsize=7)
    ax.set_ylabel("fraction of raw sites kept")
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=7, frameon=False, ncol=3, loc="lower right")
    ax.set_title("Sites retained per window after filtering (whole genome)", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    fig.savefig(os.path.splitext(path)[0] + ".pdf")
    plt.close(fig)


def filter_audit_plot(df, cfg, path):
    """Whole-genome view of what the site filter removes: retention per window and, below, the share of sites each reason flags."""
    if df is None or df.empty:
        return
    chroms, off, tot = _offsets(df, cfg.window_size, cfg)
    x = (df.chrom.map(off) + (df.start + df.end) / 2).to_numpy()
    reasons = [c for c in df.columns if c not in ("chrom", "widx", "sites", "removed", "start", "end", "retained")]
    reasons = sorted(reasons, key=lambda r: -df[r].sum())[:5]
    fig, axes = plt.subplots(2, 1, figsize=(14, 6), sharex=True)
    axes[0].scatter(x, df.retained, s=6, c="#4a4a4a")
    axes[0].set_ylabel("fraction of sites kept")
    axes[0].set_ylim(0, 1.02)
    for r, col in zip(reasons, ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e"]):
        axes[1].plot(x, df[r] / df.sites, lw=0.7, c=col, label=r)
    axes[1].set_ylabel("share of sites flagged")
    axes[1].set_ylim(0, 1.02)
    axes[1].legend(fontsize=7, frameon=False, ncol=5, loc="upper right")
    mids = [off[c] + df[df.chrom == c].end.max() / 2 for c in chroms]
    axes[1].set_xticks(mids)
    axes[1].set_xticklabels([label_of(c, cfg) for c in chroms], rotation=90, fontsize=7)
    fig.suptitle("What the site filter removes, by window (whole genome)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    fig.savefig(os.path.splitext(path)[0] + ".pdf")
    plt.close(fig)
