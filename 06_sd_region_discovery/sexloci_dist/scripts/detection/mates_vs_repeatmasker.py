#!/usr/bin/env python3
"""mates_vs_repeatmasker.py - identify the transposon family of an insertion from the reference's own RepeatMasker annotation.

The idea: an insertion that is absent from the reference but is a copy of a multi-copy transposon leaves reads whose mates map to OTHER copies of
that transposon elsewhere in the genome. Those other copies are already annotated by RepeatMasker, so the repeat names at the mates' positions say
which family the inserted element belongs to. Mates that are unmapped, on the same chromosome nearby or in unannotated sequence say nothing.

For each BAM and region this collects the mapped reads whose mate is on another chromosome (or more than --far bp away), looks up the mate's
position in the RepeatMasker annotation and reports, per repeat family, how many mates land there compared with how many would land there if mates
fell at random (the family's share of the genome). A real signal is one or a few families with large enrichment, in the carriers and not in controls.

    python3 mates_vs_repeatmasker.py --repeatmasker GCF_041146795.1_Mzebra_GT3a_rm.out.gz \\
        --bams SAMEA117806260.region.bam SAMEA117806261.region.bam AG_552_1.region.bam \\
        --regions NC_135173.1:50764700-50765300 NC_135173.1:50786430-50787030

With no arguments it fetches the lab's annotation for --genome-version from Dropbox (the .out file, which carries the repeat class, plus name_map.tsv and
the manifest) into ~/Temp/SexDeterminationAnalysis/annotation/<genome>/, translates the annotation's sequence names to the BAM's, and scans the two
carriers and one control in ~/Temp/SexDeterminationAnalysis/insertion_scan/bams at the two reference copies of candidate #1:

    python3 mates_vs_repeatmasker.py

--repeatmasker uses a local file instead (a RepeatMasker .out(.gz), or a BED/GFF with the repeat name in column 4 / in the Name=, Target= or ID=
attribute; a BED may carry the class in column 5). BAMs are indexed automatically if they are not. Needs samtools, and rclone for the download.
"""
import argparse
import collections
import gzip
import json
import os
import re
import subprocess
import sys
from bisect import bisect_right

GENOME_SIZE_DEFAULT = 962_000_000          # M_zebra_GT3a: 962 Mb (Kumar et al. 2025)
REMOTE_TEMPLATE = "ptm_dropbox:/COS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/YH_SD_Paper/genomes/{genome}/annotation/"
LOCAL_DIR_TEMPLATE = "~/Temp/SexDeterminationAnalysis/annotation/{genome}"
ANNOTATION_FILES = ["Mzebra_GT3a.fasta.out", "Mzebra_GT3a.fasta.out.gff", "name_map.tsv", "annotation_manifest.json"]
BAM_DIR = "~/Temp/SexDeterminationAnalysis/insertion_scan/bams"
DEFAULT_BAMS = ["SAMEA117806260.region.bam", "SAMEA117806261.region.bam", "AG_552_1.region.bam"]     # two carriers, then a control
DEFAULT_REGIONS = ["NC_135173.1:50764700-50765300", "NC_135173.1:50786500-50787300"]                  # candidate #1 in both reference copies
DEFAULT_SITES = [50765052, 50786786]                                                                  # its junction in each copy (used to split reads left / right of the site)


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


# ------------------------------------------------------------------------------------------------ annotation
def read_annotation(path):
    """{chrom: [(start, end, name, class)]} with 1-based inclusive coordinates."""
    ivs = collections.defaultdict(list)
    low = path.lower().replace(".gz", "")
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            if low.endswith(".out"):
                f = line.split()
                if len(f) < 11 or not f[0].isdigit():
                    continue
                try:
                    extra = None
                    if len(f) >= 14:                         # position in the library consensus: '+' = begin end (left); 'C' = (left) end begin
                        num = lambda x: int(x.strip("()"))
                        extra = ("+", num(f[11]), num(f[12]), num(f[13])) if f[8] == "+" else ("C", num(f[13]), num(f[12]), num(f[11]))
                    ivs[f[4]].append((int(f[5]), int(f[6]), f[9], f[10], extra))
                except ValueError:
                    continue
            elif low.endswith((".gff", ".gff3", ".gtf")):
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                attrs = {}
                for part in f[8].replace("; ", ";").split(";"):
                    if "=" in part:
                        k, v = part.split("=", 1)
                        attrs[k.strip()] = v.strip()
                name = attrs.get("Name") or attrs.get("Target", "").split(" ")[0].replace("Motif:", "") or attrs.get("ID", f[2])
                ivs[f[0]].append((int(f[3]), int(f[4]), name, attrs.get("repeat_class") or attrs.get("class") or f[2], None))
            else:                                                       # BED: chrom start(0-based) end name [class]
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                ivs[f[0]].append((int(f[1]) + 1, int(f[2]), f[3], f[4] if len(f) > 4 else "", None))
    return ivs


class Index:
    def __init__(self, ivs):
        self.chrom = {}
        for c, lst in ivs.items():
            lst = sorted(lst)
            starts = [x[0] for x in lst]
            pm, m = [], 0
            for x in lst:
                m = max(m, x[1])
                pm.append(m)
            self.chrom[c] = (lst, starts, pm)

    def overlapping(self, c, a, b):
        if c not in self.chrom:
            return []
        lst, starts, pm = self.chrom[c]
        j = bisect_right(starts, b) - 1
        out = []
        while j >= 0 and pm[j] >= a:
            if lst[j][1] >= a and lst[j][0] <= b:
                out.append(lst[j])
            j -= 1
        return out


def genome_shares(ivs, genome_size):
    """Share of the genome covered by each repeat name and class (nested hits overcount slightly; good enough for a background)."""
    by_name, by_class = collections.Counter(), collections.Counter()
    for lst in ivs.values():
        for s, e, name, cls, _ in lst:
            by_name[name] += e - s + 1
            by_class[cls] += e - s + 1
    return ({k: v / genome_size for k, v in by_name.items()}, {k: v / genome_size for k, v in by_class.items()})


def consensus_lengths(ivs):
    """Length of each library consensus, from the .out columns: end + (left) on '+' hits, else the longest end seen."""
    plus, other = {}, {}
    for lst in ivs.values():
        for _, _, name, _, ex in lst:
            if ex:
                if ex[0] == "+":
                    plus[name] = max(plus.get(name, 0), ex[2] + ex[3])
                else:
                    other[name] = max(other.get(name, 0), ex[2])
    return {n: plus.get(n) or other.get(n, 0) for n in set(plus) | set(other)}


def consensus_position(hit, mate_start, mate_end):
    """Where on the library consensus a mate lies, assuming it sits at the same offset inside the annotated copy (ignores small indels)."""
    s, e, _, _, ex = hit
    if not ex:
        return None
    mid = (mate_start + mate_end) // 2
    off = min(max(mid - s, 0), e - s)
    cp = ex[1] + off if ex[0] == "+" else ex[2] - off
    return max(cp, 1)


# ------------------------------------------------------------------------------------------------ fetching the lab's annotation
def fetch_annotation(remote: str, local_dir: str, refresh=False):
    """Copy the annotation files from Dropbox with rclone (only those that exist there and are not already here). Returns {file name: local path}."""
    local_dir = os.path.expanduser(local_dir)
    os.makedirs(local_dir, exist_ok=True)
    have = {f: os.path.join(local_dir, f) for f in ANNOTATION_FILES if os.path.exists(os.path.join(local_dir, f)) and os.path.getsize(os.path.join(local_dir, f)) > 0}
    need_out = not any(f in have for f in ("Mzebra_GT3a.fasta.out", "Mzebra_GT3a.fasta.out.gff"))
    if refresh or need_out or "name_map.tsv" not in have:
        cmd = ["rclone", "copy", remote, local_dir] + [x for f in ANNOTATION_FILES for x in ("--include", f)]
        print(f"[annotation] rclone copy {remote} -> {local_dir}", flush=True)
        try:
            p = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError:
            sys.exit("rclone is not on the PATH; install it or pass --repeatmasker with a local file")
        if p.returncode != 0:
            sys.exit(f"rclone copy failed: {p.stderr.strip()[:300]}\nCheck the remote with: rclone lsf {remote}")
        have = {f: os.path.join(local_dir, f) for f in ANNOTATION_FILES if os.path.exists(os.path.join(local_dir, f)) and os.path.getsize(os.path.join(local_dir, f)) > 0}
    else:
        print(f"[annotation] using the copy in {local_dir} (--refresh to download again)")
    return have


def load_name_map(path, annotation_names):
    """{annotation sequence name: NCBI name}. The file's layout is not assumed: the column that matches the annotation's sequence names and the column
    of NC_/NW_ names are found by overlap. Returns ({}, message) when no translation is needed or none can be found."""
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t") if "\t" in line else line.split()
            if len(f) >= 2 and f[0] and not line.startswith("#"):
                rows.append(f)
    if not rows:
        return {}, f"{os.path.basename(path)} is empty"
    ncols = max(len(r) for r in rows)
    cols = [[r[i] for r in rows if len(r) > i] for i in range(ncols)]
    ann = set(annotation_names)
    over = [len(ann & set(c)) for c in cols]
    nc = [sum(1 for v in c if re.match(r"^N[CWZ]_\d+\.\d+$", v)) for c in cols]
    ca, cn = max(range(ncols), key=lambda i: over[i]), max(range(ncols), key=lambda i: nc[i])
    if over[ca] == 0:
        return {}, f"no column of {os.path.basename(path)} matches the annotation's sequence names (e.g. {sorted(ann)[:2]})"
    if nc[cn] == 0:
        return {}, f"{os.path.basename(path)} has no NC_/NW_ names"
    if ca == cn:
        return {}, "the annotation already uses the NCBI sequence names"
    return {r[ca]: r[cn] for r in rows if len(r) > max(ca, cn)}, f"translated {over[ca]} annotation sequence names to NCBI names using {os.path.basename(path)}"


def ensure_indexed(bam):
    if not (os.path.exists(bam + ".bai") or os.path.exists(bam + ".csi") or os.path.exists(re.sub(r"\.bam$", ".bai", bam))):
        print(f"[index] samtools index {bam}", flush=True)
        p = subprocess.run(["samtools", "index", bam], capture_output=True, text=True)
        if p.returncode != 0:
            sys.exit(f"samtools index {bam} failed: {p.stderr.strip()[:200]} (is the BAM coordinate-sorted?)")


# ------------------------------------------------------------------------------------------------ reads
def _aligned_len(cigar, default):
    """Reference span of the aligned part of a read (soft clips excluded), so a clipped read's position reflects where it really aligns."""
    n = [(int(x), o) for x, o in re.findall(r"(\d+)([MIDNSHP=X])", cigar)] if cigar != "*" else []
    return sum(x for x, o in n if o in "MDN=X") if n else default


def anomalous_mates(sam_lines, far=3000):
    """(mate chrom, mate start, mate end, kind) for mapped reads whose mate is elsewhere; plus counts of unmapped mates."""
    out, unmapped_mate = [], 0
    for line in sam_lines:
        if not line or line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        flag, rname, rnext, pnext, tlen, seq = int(f[1]), f[2], f[6], int(f[7]), int(f[8]), f[9]
        if flag & 0x4 or flag & (0x100 | 0x400 | 0x800):
            continue
        if flag & 0x8:
            unmapped_mate += 1
        elif rnext not in ("=", "*") and rnext != rname:
            out.append((rnext, pnext, pnext + len(seq) - 1, "other chromosome", int(f[3]) + _aligned_len(f[5], len(seq)) // 2))
        elif rnext in ("=", rname) and abs(tlen) > far:
            out.append((rname, pnext, pnext + len(seq) - 1, f"same chromosome, >{far} bp", int(f[3]) + _aligned_len(f[5], len(seq)) // 2))
    return out, unmapped_mate


def side_table(label, mates, sides, index, cons_len, top=6):
    """Per family: how many mates come from reads left vs right of the site, and where on the consensus they fall. If two families are the two
    ends of one element, each should take mates from one side."""
    rows = collections.defaultdict(lambda: {"left": [], "right": []})
    for m_, side in zip(mates, sides):
        seen = set()
        for h in index.overlapping(m_[0], m_[1], m_[2]):
            if h[2] in seen:
                continue
            seen.add(h[2])
            rows[h[2]][side].append(consensus_position(h, m_[1], m_[2]))
    if not rows:
        return []
    out = [f"{label}: mates by the side of the site their read lies on"]
    for name, d in sorted(rows.items(), key=lambda kv: -(len(kv[1]['left']) + len(kv[1]['right'])))[:top]:
        def span(v):
            v = sorted(x for x in v if x is not None)
            return f"{len(v):2d}" + (f" (consensus {v[0]}-{v[-1]})" if v else "")
        out.append(f"      {name:20s} {cons_len.get(name, 0):>6,} bp   left of site: {span(d['left']):28s} right of site: {span(d['right'])}")
    return out


def view(bam, region):
    p = subprocess.run(["samtools", "view", bam, region], capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"samtools view {bam} {region} failed: {p.stderr.strip()[:200]} (is the BAM indexed? samtools index {bam})")
    return p.stdout.splitlines()


# ------------------------------------------------------------------------------------------------ report
def summarise(label, mates, unmapped_mate, index, share_name, share_class, top=10, verbose=False, cons_len=None, end_zone=600):
    lines = [f"{label}: {len(mates)} mates mapped elsewhere ({sum(1 for m in mates if m[3] == 'other chromosome')} on other chromosomes), "
             f"{unmapped_mate} reads with an unmapped mate"]
    if not mates:
        return lines, collections.Counter(), collections.Counter()
    names, classes, unannot = collections.Counter(), collections.Counter(), 0
    detail, cpos = [], collections.defaultdict(list)
    for m_ in mates:
        c, a, b, kind = m_[:4]
        hits = index.overlapping(c, a, b)
        if not hits:
            unannot += 1
            detail.append(f"{c}:{a}-{b}  ({kind})  not in a repeat annotation")
            continue
        # one vote per distinct (name, class) so nested annotations do not count a mate twice for the same family
        uniq = {(h[2], h[3]) for h in hits}
        for name, cls in uniq:
            names[name] += 1
            classes[cls] += 1
        seen = set()
        for h in hits:
            cp = consensus_position(h, a, b)
            if cp is not None and h[2] not in seen:
                seen.add(h[2])
                cpos[h[2]].append(cp)
        detail.append(f"{c}:{a}-{b}  ({kind})  " + "; ".join(f"{n} [{k}]" for n, k in sorted(uniq)))
    n = len(mates)
    if unannot == n and not any(m_[0] in index.chrom for m_ in mates):
        lines.append(f"   WARNING: none of the mates' chromosomes (e.g. {mates[0][0]}) appear in the annotation (which has e.g. "
                     f"{next(iter(index.chrom))}): the sequence names differ between the BAM and the annotation file")
    lines.append(f"   in an annotated repeat: {n - unannot}/{n} ({(n - unannot) / n:.0%}); unannotated: {unannot}")
    lines.append("   family / class (mates there, share of genome, enrichment):")
    for cls, k in classes.most_common(top):
        exp = share_class.get(cls, 0) * n
        lines.append(f"      {cls:28s} {k:3d} mates   genome share {share_class.get(cls, 0):6.2%}   {k / exp:6.1f}x" if exp else
                     f"      {cls:28s} {k:3d} mates   genome share    n/a")
    lines.append("   repeat names:")
    for name, k in names.most_common(top):
        exp = share_name.get(name, 0) * n
        lines.append(f"      {name:28s} {k:3d} mates   genome share {share_name.get(name, 0):6.3%}   {k / exp:6.1f}x" if exp else
                     f"      {name:28s} {k:3d} mates   genome share    n/a")
    cons_len = cons_len or {}
    for name, k in names.most_common(top):
        ps = sorted(cpos.get(name, []))
        if len(ps) >= 3:
            L = cons_len.get(name, 0)
            near = sum(1 for x in ps if x <= end_zone or (L and x >= L - end_zone))
            lines.append(f"   {name}: library consensus {L:,} bp; {len(ps)} mates fall at consensus positions "
                         f"{', '.join(str(x) for x in ps[:40])}{' ...' if len(ps) > 40 else ''}"
                         + (f"  ->  {near}/{len(ps)} within {end_zone} bp of an end" if L > 3 * end_zone else
                            f"  (consensus under {3 * end_zone:,} bp: the end zones would cover most of it, so look at the positions themselves)"))
    if verbose:
        lines.append("   each mate:")
        lines += ["      " + d for d in detail]
    return lines, names, classes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repeatmasker", help="a local RepeatMasker .out(.gz), BED or GFF; default: fetch the lab's annotation from Dropbox")
    ap.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    ap.add_argument("--annotation-remote", help=f"rclone folder holding the annotation (default {REMOTE_TEMPLATE})")
    ap.add_argument("--annotation-dir", help=f"where the downloaded annotation is kept (default {LOCAL_DIR_TEMPLATE})")
    ap.add_argument("--refresh", action="store_true", help="download the annotation again")
    ap.add_argument("--bams", nargs="+", help=f"region BAMs, carriers first, controls after (default: {', '.join(DEFAULT_BAMS)} in {BAM_DIR})")
    ap.add_argument("--regions", nargs="+", default=DEFAULT_REGIONS, help="chrom:start-end windows around the site (both reference copies if the region is duplicated)")
    ap.add_argument("--sites", nargs="+", type=int, help="the junction position in each region (same order), to split reads left / right of it; "
                    f"default {DEFAULT_SITES} for the default regions")
    ap.add_argument("--far", type=int, default=3000, help="same-chromosome mates farther than this count as anomalous (default 3000)")
    ap.add_argument("--genome-size", type=float, default=GENOME_SIZE_DEFAULT, help="for the background share (default 962 Mb, GT3a)")
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--verbose", action="store_true", help="list every mate with its annotation")
    a = ap.parse_args()
    name_map_path = None
    if a.repeatmasker:
        annot = a.repeatmasker
    else:
        files = fetch_annotation(a.annotation_remote or REMOTE_TEMPLATE.format(genome=a.genome_version),
                                 a.annotation_dir or LOCAL_DIR_TEMPLATE.format(genome=a.genome_version), a.refresh)
        annot = files.get("Mzebra_GT3a.fasta.out") or files.get("Mzebra_GT3a.fasta.out.gff")
        if not annot:
            sys.exit("no .out or .out.gff file was found in the annotation folder")
        if annot.endswith(".gff"):
            print("[annotation] only the .out.gff is available: it has repeat names but no repeat classes, so the class table will not be informative")
        name_map_path = files.get("name_map.tsv")
        if files.get("annotation_manifest.json"):
            try:
                print("[annotation] manifest:", json.dumps(json.load(open(files["annotation_manifest.json"])))[:500])
            except Exception:
                pass
    ivs = read_annotation(annot)
    if not ivs:
        sys.exit(f"no annotation records read from {annot}")
    if name_map_path:
        mapping, msg = load_name_map(name_map_path, ivs.keys())
        print("[annotation]", msg)
        if mapping:
            moved = collections.defaultdict(list)
            for c, lst in ivs.items():
                moved[mapping.get(c, c)] += lst
            ivs = moved
    idx = Index(ivs)
    cons_len = consensus_lengths(ivs)
    sn, sc = genome_shares(ivs, a.genome_size)
    print(f"annotation: {sum(len(v) for v in ivs.values()):,} records on {len(ivs)} sequences; "
          f"{sum(sc.values()):.0%} of a {a.genome_size / 1e6:.0f} Mb genome (before de-duplicating nested hits)\n")
    bams = a.bams or [os.path.join(os.path.expanduser(BAM_DIR), b) for b in DEFAULT_BAMS]
    for bam in bams:
        if not os.path.exists(bam):
            sys.exit(f"{bam} not found; give --bams, or run find_insertion.py first so the region BAMs exist")
        ensure_indexed(bam)
    for bam in bams:
        label = os.path.basename(bam).replace(".region.bam", "").replace(".bam", "")
        mates, um, sides = [], 0, []
        sites = a.sites or (DEFAULT_SITES if a.regions == DEFAULT_REGIONS else None)
        for i, r in enumerate(a.regions):
            m, u = anomalous_mates(view(bam, r), a.far)
            mates += m
            um += u
            if sites and i < len(sites):
                sides += ["left" if x[4] < sites[i] else "right" for x in m]
        lines, names, classes = summarise(label, mates, um, idx, sn, sc, a.top, a.verbose, cons_len)
        print("\n".join(lines))
        if mates and len(sides) == len(mates):
            print("\n".join(side_table(label, mates, sides, idx, cons_len)))
        print()
    print("A transposon family shows as a name or class with many mates, a large enrichment (>10x) and nothing similar in the controls. "
          "A handful of mates is only suggestive.")


if __name__ == "__main__":
    main()
