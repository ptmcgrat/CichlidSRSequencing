"""Build the candidate-QTN browser: one self-contained HTML file.

Pulls everything from cloud storage, assembles a single HTML file with the
genotype matrix embedded, and pushes it back to the WebServer folder. Nothing
needs to be exported by hand and nothing needs to be running to view the result --
the output is a static file that opens by double-clicking it out of Dropbox.

    python3 buildBrowser.py                      # download, build, upload
    python3 buildBrowser.py --no-download        # rebuild from local copies
    python3 buildBrowser.py --no-upload          # build locally only

Inputs (all fetched via FileManager):
    <NikeshDir>/QTG_Candidates/          per-sample VCFs + manifests + sites files
    <NikeshDir>/QTG_Candidates/candidateQTNs_all.tsv
    sample database, alignment database  (Google Sheets)
    <NikeshDir>/Annotation/<gtf>         gene models for the region

Output:
    <NikeshDir>/WebServer/qtn_browser.html
"""

import argparse
import glob
import gzip
import json
import os
import re
import sys
from collections import defaultdict

import pandas as pd

from helper_modules.file_manager import FileManager as FM
from helper_modules import pipeline_checks as pc
from helper_modules.pipeline_checks import log, warn, PipelineError

TEMPLATE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "helper_modules", "browser_template.html")


def parse_args():
    p = argparse.ArgumentParser(description="Build the candidate-QTN browser.")
    p.add_argument("--genome-version", default="Mzebra_GT3_NCBI")
    p.add_argument("--gtf", default=None,
                   help="Override the GTF path. Defaults to FileManager's "
                        "localGTFFile, downloaded from cloud storage. Pass a local "
                        "path to use a different annotation without uploading it.")
    p.add_argument("--candidates", default="QTG_Candidates/candidateQTNs_all.tsv",
                   help="Master candidate table, relative to the Nikesh directory.")
    p.add_argument("--out-name", default="qtn_browser.html")
    p.add_argument("--min-depth", type=int, default=6,
                   help="Initial Min DP setting in the UI (adjustable in-page)")
    p.add_argument("--no-download", action="store_true",
                   help="Use whatever is already local. Faster for iterating on layout.")
    p.add_argument("--no-upload", action="store_true",
                   help="Build locally without pushing to the WebServer folder.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------

def fetch(fm_obj, args):
    cand_dir = fm_obj.localNikeshDir + "QTG_Candidates/"
    tsv_path = fm_obj.localNikeshDir + args.candidates

    # An explicit --gtf is treated as a local file and never downloaded, so a
    # different annotation can be tried without uploading it first.
    gtf_path = args.gtf or getattr(fm_obj, "localGTFFile", None)
    if gtf_path is None:
        warn("FileManager has no localGTFFile attribute and no --gtf was given; "
             "the gene track will be empty")

    if args.no_download:
        log("--no-download: using local copies")
    else:
        log("downloading QTG_Candidates/ (VCFs, manifests, sites files)")
        fm_obj.downloadData(cand_dir.rstrip("/"))

        try:
            fm_obj.downloadData(tsv_path)
        except FileNotFoundError:
            warn(f"candidate table not found in cloud storage at {tsv_path}")

        if gtf_path and not args.gtf:
            # If the annotation gets bgzipped later, localGTFFile may still name the
            # plain file (or the reverse). Try the stated path, then the other form.
            alt = gtf_path[:-3] if gtf_path.endswith(".gz") else gtf_path + ".gz"
            for cand in (gtf_path, alt):
                try:
                    fm_obj.downloadData(cand)
                    gtf_path = cand
                    break
                except FileNotFoundError:
                    continue
            else:
                warn(f"GTF not found in cloud storage as {gtf_path} or {alt}")

    pc.require_file(tsv_path, "candidate table")
    if gtf_path and not os.path.exists(gtf_path):
        warn(f"no GTF at {gtf_path}; the gene track will be empty")
        gtf_path = None
    if gtf_path:
        log(f"gene models from {gtf_path}")
    return cand_dir, tsv_path, gtf_path


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def read_variants(tsv_path):
    dt = pd.read_csv(tsv_path, sep="\t")
    out = []
    for _, r in dt.iterrows():
        info = str(r.Info)
        m = re.search(r"TYPE=([A-Za-z]+)", info)
        ref, alt = str(r.Reference).upper(), str(r.Alt).upper()
        out.append({"pos": int(r.Position), "chrom": str(r.Chromosome),
                    "ref": ref, "alt": alt,
                    "type": m.group(1) if m else "?",
                    "big": len(ref) > 10 or len(alt) > 10,
                    "len": abs(len(alt) - len(ref)),
                    "q": round(float(r.Q), 1) if pd.notna(r.Q) else None})
    out.sort(key=lambda v: v["pos"])
    return out


def read_samples(fm_obj):
    """Sample metadata joined to per-sample coverage.

    OtherYHs is derived here rather than hardcoded in the UI: if the sample sheet
    grows a real OtherYHs category it is used as-is, and the fallback quietly
    stops firing.
    """
    fm_obj.readSampleDatabase()
    fm_obj.readAlignmentDatabase()
    sdt, adt = fm_obj.sample_dt, fm_obj.alignment_dt
    cov = dict(zip(adt.SampleID, pd.to_numeric(adt.get("Coverage"), errors="coerce")))

    meta = {}
    for _, r in sdt.iterrows():
        cat = r.Category if pd.notna(r.get("Category")) else ""
        if not cat and isinstance(r.Species, str) and "yellow head" in r.Species.lower():
            cat = "OtherYHs"
        meta[r.SampleID] = {
            "id": r.SampleID,
            "sp": r.Species if pd.notna(r.Species) else "",
            "eco": r.Ecogroup if pd.notna(r.Ecogroup) else "Unassigned",
            "sub": r.Subgroup if pd.notna(r.get("Subgroup")) else "Unassigned",
            "cat": cat,
            "sex": r.Sex if pd.notna(r.Sex) else "",
            "inv": int(r.Inversion10) if pd.notna(r.get("Inversion10")) else -1,
            "lab": int(r.LabReared) if pd.notna(r.get("LabReared")) else -1,
            "cov": round(float(cov.get(r.SampleID) or 0), 1),
        }
    return meta


def read_genotypes(cand_dir, variants, meta):
    """Read every per-sample VCF into flat genotype arrays.

    Records are keyed on chrom+pos, NOT on alleles: bcftools re-anchors indels, so
    an output record may carry different REF/ALT than the candidate table asked
    for. No two candidates share a position, so position alone is unique.
    """
    pos_index = {(v["chrom"], v["pos"]): i for i, v in enumerate(variants)}
    if len(pos_index) != len(variants):
        raise PipelineError("two candidates share a position; the position key "
                            "is no longer unique and the loader needs alleles")

    vcfs = sorted(glob.glob(cand_dir + "*_candidate_QTNs.vcf.gz"))
    log(f"found {len(vcfs)} per-sample VCFs")
    if not vcfs:
        raise PipelineError(f"no VCFs in {cand_dir}")

    samples, skipped = [], []
    for path in vcfs:
        sid = os.path.basename(path).replace("_candidate_QTNs.vcf.gz", "")
        if sid not in meta:
            skipped.append(sid)
            continue
        samples.append((sid, path))
    if skipped:
        warn(f"{len(skipped)} VCFs have no sample-database row and were skipped: "
             f"{skipped[:5]}")

    NV, NS = len(variants), len(samples)
    gt = [-1] * (NV * NS)
    dp = [0] * (NV * NS)
    ada = [0] * (NV * NS)
    gq = [0] * (NV * NS)
    unc = [1] * (NV * NS)          # assume uncallable until a record says otherwise
    rerep = defaultdict(int)

    for si, (sid, path) in enumerate(samples):
        n = 0
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                key = (f[0], int(f[1]))
                vi = pos_index.get(key)
                if vi is None:
                    continue
                v = variants[vi]
                if f[3].upper() != v["ref"] or f[4].upper() != v["alt"]:
                    rerep[vi] += 1
                fmt = f[8].split(":")
                val = f[9].split(":")
                d = dict(zip(fmt, val))
                k = vi * NS + si
                g = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1,
                     "1/1": 2, "1|1": 2}.get(d.get("GT", "./."), -1)
                ad = d.get("AD", "")
                a0 = a1 = 0
                if "," in ad:
                    try:
                        parts = [int(x) for x in ad.split(",")[:2]]
                        a0, a1 = parts[0], parts[1]
                    except ValueError:
                        pass
                depth = int(d["DP"]) if d.get("DP", "").isdigit() else a0 + a1
                gt[k] = g
                dp[k] = min(255, depth)
                ada[k] = min(255, a1)
                gq[k] = min(99, int(d["GQ"])) if d.get("GQ", "").isdigit() else 0
                unc[k] = 0                      # a record exists -> it was callable
                n += 1
        if n == 0:
            warn(f"{sid}: VCF contained no matching records")

    if rerep:
        log(f"{len(rerep)} site(s) came back with re-anchored alleles in at least "
            f"one sample; matched on position as intended")

    return [meta[s] for s, _ in samples], gt, dp, ada, gq, unc


def apply_manifests(cand_dir, variants, samples, unc):
    """Mark sites the pipeline recorded as missing for a given sample.

    read_genotypes already marks anything without a record as uncallable; the
    manifests confirm it and let us report WHY (no reads vs no usable reads).
    """
    pos_index = {v["pos"]: i for i, v in enumerate(variants)}
    sid_index = {s["id"]: i for i, s in enumerate(samples)}
    NS = len(samples)
    reasons = defaultdict(list)
    n = 0
    for path in sorted(glob.glob(cand_dir + "*.manifest.json")):
        try:
            m = json.load(open(path))
        except Exception:
            continue
        si = sid_index.get(m.get("sample_id"))
        if si is None:
            continue
        for e in m.get("missing_sites", []):
            vi = pos_index.get(e.get("pos"))
            if vi is None:
                continue
            unc[vi * NS + si] = 1
            n += 1
            reasons[vi].append(0 if not e.get("depth_raw") else
                               1 if not e.get("depth_mq20") else 2)
    log(f"manifests marked {n} sample-site pairs as uncallable")
    # 0 = no reads, 1 = reads but none usable, 2 = usable reads yet no record
    summary = {vi: {"nodata": r.count(0), "unmappable": r.count(1),
                    "unexplained": r.count(2)} for vi, r in reasons.items()}
    return summary


def annotate_mappability(variants, dp, unc, samples):
    """Flag sites that are systematically starved of reads across all samples.

    Each site's median depth is compared to each sample's own coverage, so a site
    that reads shallow everywhere is distinguished from a shallow sample. This is
    the same signal as the MQ0 fraction and needs no BAM access.
    """
    NS = len(samples)
    covs = [s["cov"] or 1 for s in samples]
    for vi, v in enumerate(variants):
        ratios = []
        for si in range(NS):
            k = vi * NS + si
            if unc[k]:
                continue
            ratios.append(dp[k] / max(1.0, covs[si]))
        if not ratios:
            v["lowmap"] = True
            v["depth_ratio"] = 0.0
            continue
        ratios.sort()
        med = ratios[len(ratios) // 2]
        v["depth_ratio"] = round(med, 2)
        v["lowmap"] = med < 0.35
    n = sum(1 for v in variants if v["lowmap"])
    log(f"{n}/{len(variants)} sites flagged low-mappability (median depth < 35% "
        f"of sample coverage)")


# ---------------------------------------------------------------------------
# Gene models
# ---------------------------------------------------------------------------

def open_maybe_gz(path):
    """Open a text file whether or not it is gzipped.

    Sniffs the two-byte gzip magic number rather than trusting the extension, so
    a compressed file named .gtf or a plain file named .gz both read correctly.
    """
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def read_gtf(path, chrom, lo, hi):
    """Collapse each gene to the union of its exons across all transcripts."""
    if not path:
        return []
    genes = {}
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != chrom:
                continue
            start, end = int(f[3]), int(f[4])
            if end < lo or start > hi:
                continue
            attr = f[8]
            gid = re.search(r'gene_id "([^"]+)"', attr)
            if not gid:
                continue
            gid = gid.group(1)
            name = re.search(r'gene "([^"]+)"', attr) or re.search(r'gene_name "([^"]+)"', attr)
            bt = re.search(r'gene_biotype "([^"]+)"', attr)
            g = genes.setdefault(gid, {"name": name.group(1) if name else gid,
                                       "start": start, "end": end, "strand": f[6],
                                       "biotype": bt.group(1) if bt else "", "exons": []})
            if f[2] in ("gene", "transcript"):
                g["start"] = min(g["start"], start)
                g["end"] = max(g["end"], end)
                if bt and not g["biotype"]:
                    g["biotype"] = bt.group(1)
            elif f[2] == "exon":
                g["exons"].append([start, end])

    out = []
    for g in genes.values():
        ex = sorted(g["exons"])
        merged = []
        for a, b in ex:
            if merged and a <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        g["exons"] = merged
        out.append(g)
    out.sort(key=lambda g: g["start"])
    log(f"{len(out)} genes in the region")
    return out


def annotate_context(variants, genes):
    """Describe each candidate's position relative to the annotation."""
    for v in variants:
        p = v["pos"]
        inside = [g for g in genes if g["start"] <= p <= g["end"]]
        if inside:
            g = inside[0]
            exonic = any(a <= p <= b for a, b in g["exons"])
            v["ctx"] = f"{'exon' if exonic else 'intron'} of {g['name']}"
            v["gene"] = g["name"]
        elif genes:
            near = min(genes, key=lambda g: min(abs(p - g["start"]), abs(p - g["end"])))
            d = min(abs(p - near["start"]), abs(p - near["end"]))
            side = "upstream" if p < near["start"] else "downstream"
            v["ctx"] = f"{d/1000:.1f} kb {side} of {near['name']}"
            v["gene"] = ""
        else:
            v["ctx"] = ""
            v["gene"] = ""


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    pc.require_tools(("bcftools",))
    fm_obj = FM(genome_version=args.genome_version)

    cand_dir, tsv_path, gtf_path = fetch(fm_obj, args)

    variants = read_variants(tsv_path)
    log(f"{len(variants)} candidates, "
        f"{variants[0]['pos']:,}-{variants[-1]['pos']:,}")

    meta = read_samples(fm_obj)
    samples, gt, dp, ada, gq, unc = read_genotypes(cand_dir, variants, meta)
    log(f"{len(samples)} samples loaded")

    miss = apply_manifests(cand_dir, variants, samples, unc)
    annotate_mappability(variants, dp, unc, samples)

    chrom = variants[0]["chrom"]
    lo, hi = variants[0]["pos"] - 5000, variants[-1]["pos"] + 5000
    genes = read_gtf(gtf_path, chrom, lo, hi)
    annotate_context(variants, genes)
    for vi, v in enumerate(variants):
        v["miss"] = miss.get(vi, {})

    payload = {"samples": samples, "variants": variants, "genes": genes,
               "gt": gt, "dp": dp, "ada": ada, "gq": gq, "unc": unc,
               "minDP": args.min_depth,
               "built": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M")}

    pc.require_file(TEMPLATE, "browser template")
    html = open(TEMPLATE).read().replace("__DATA__",
                                         json.dumps(payload, separators=(",", ":")))

    web_dir = fm_obj.localNikeshDir + "WebServer/"
    os.makedirs(web_dir, exist_ok=True)
    out_path = web_dir + args.out_name
    with open(out_path, "w") as fh:
        fh.write(html)
    log(f"wrote {out_path} ({os.path.getsize(out_path)/1e6:.1f} MB)")

    if not args.no_upload:
        fm_obj.uploadData(out_path)
        log(f"uploaded to WebServer/{args.out_name}")
    else:
        log("--no-upload: left in place locally")


if __name__ == "__main__":
    main()