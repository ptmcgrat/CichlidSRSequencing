"""sv-genotype: genotype the gsdf duplication AND long insertions (transposons) in every sample of the database, streaming only the
parts of the BAM files that are needed, and write the results back into the Sample Database.

Per sample ONE `samtools view` request through the read-only rclone web server fetches all regions together:
  * the 130 kb around gsdf (copy number from read depth, as `gsdf-genotype`);
  * a ~1 kb window around each insertion site (reads, including unmapped mates stored at their mate's position).
The insertion genotype is the lab's IndelReadClassifier (longindel.py): each read pair is aligned to a local reference allele and a local
alternative allele (flank + inserted sequence + flank) and counted as supporting REF, ALT, equal or uninformative when it crosses a
breakpoint; a binomial model gives 0/0, 0/1, 1/1 or ./. with a phred quality.

Each locus has its own cache, so a new insertion site only fetches its own window for samples already done, and gsdf results from an
earlier `gsdf-genotype` run are imported instead of measured again.
"""
from __future__ import annotations

import gzip
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from fnmatch import fnmatch
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from .gsdf import (CHROM as GSDF_CHROM, DEFAULT_R2, DEFAULT_STEP, MEAS_COLS, REGION as GSDF_REGION, HttpServer, Runner, _classify_failure,
                   call_genotypes, default_runner, fit_anchors, measure_bam, pedigree_roles, resolve_extra_columns)
from .longindel import IndelReadClassifier, call_from_counts, make_backend, reads_from_samtools_lines

TE_COLS = ["SampleID", "Site", "Status", "GT", "n_ref", "n_alt", "n_equal", "n_uninf", "GQ", "n_reads", "Message"]


# ------------------------------------------------------------------------------------------------ the insertion sites
@dataclass
class IndelSite:
    name: str
    chrom: str
    pos: int                       # 1-based, VCF anchor base
    ref: str
    alt: str
    info: Dict[str, str] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return len(self.alt) - len(self.ref)


def load_indel_sites(path: str, names: Optional[List[str]] = None) -> List[IndelSite]:
    """Read a VCF-style file of anchored indels (CHROM POS ID REF ALT [QUAL FILTER INFO ...]); '#' lines are skipped, there need not be a
    header, extra trailing columns are ignored. Site names: --indel-names, else the ID column, else TE_<pos>."""
    sites: List[IndelSite] = []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                raise ValueError(f"{path}: expected at least 5 tab-separated columns (CHROM POS ID REF ALT), got {len(f)}: {line[:80]!r}")
            if "," in f[4]:
                raise ValueError(f"{path}: multi-allelic ALT at {f[0]}:{f[1]} is not supported; split the record")
            if not (set(f[3].upper()) <= set("ACGTN") and set(f[4].upper()) <= set("ACGTN")):
                raise ValueError(f"{path}: REF/ALT at {f[0]}:{f[1]} must contain only A, C, G, T, N (symbolic alleles are not supported)")
            info = dict(kv.split("=", 1) for kv in re.split(r"[;,]", f[7].strip('"') if len(f) > 7 else "") if "=" in kv)
            sites.append(IndelSite(name="", chrom=f[0], pos=int(f[1]), ref=f[3].upper(), alt=f[4].upper(), info=info))
            sites[-1].name = f[2] if f[2] not in (".", "") else f"TE_{f[1]}"
    if not sites:
        raise ValueError(f"no records found in {path}")
    if names:
        if len(names) != len(sites):
            raise ValueError(f"--indel-names gives {len(names)} names for {len(sites)} sites in {path}")
        for s, n in zip(sites, names):
            s.name = n
    if len({s.name for s in sites}) != len(sites):
        raise ValueError("site names must be unique: " + ", ".join(s.name for s in sites))
    return sites


def load_builtin_sites(name: Union[str, Sequence[str]] = "", names: Optional[List[str]] = None) -> List[IndelSite]:
    """The insertion sites that ship with sexloci (known_sites.py). `name` is one set, a list of sets, or 'all'; empty means the default set
    (deep_benthic_TE: TE1, TE2). Sites keep the order of the sets given. Names: `names` (one per site), else the ID of each record."""
    from .known_sites import BUILTIN, DEFAULT_SET
    want = [name] if isinstance(name, str) else list(name)
    want = [w for w in want if w] or [DEFAULT_SET]
    if "all" in want:
        want = list(BUILTIN)
    for w in want:
        if w not in BUILTIN:
            raise ValueError(f"unknown built-in site set {w!r}; available: all, {', '.join(BUILTIN)}")
    want = list(dict.fromkeys(want))
    sites = [IndelSite(name=i, chrom=c, pos=pos, ref=r.upper(), alt=a.upper(), info=dict(kv.split("=", 1) for kv in re.split(r"[;,]", info) if "=" in kv))
             for w in want for c, pos, i, r, a, info in BUILTIN[w]]
    if names:
        if len(names) != len(sites):
            raise ValueError(f"--indel-names gives {len(names)} names for the {len(sites)} sites of the built-in set(s) {', '.join(want)}; "
                             f"to rename the sites of one set, choose it with --builtin-sites")
        for st, n in zip(sites, names):
            st.name = n
    return sites


def site_groups(sites: List[IndelSite]) -> Dict[str, List[str]]:
    """The combined calls that apply: known_sites.SITE_GROUPS entries whose member sites are ALL being genotyped (and whose name is not a site name)."""
    from .known_sites import SITE_GROUPS
    have = {s.name for s in sites}
    return {g: list(m) for g, m in SITE_GROUPS.items() if g not in have and all(x in have for x in m)}


def add_group_columns(full: pd.DataFrame, groups: Dict[str, List[str]], min_informative_pairs: int = 5) -> pd.DataFrame:
    """One combined genotype per sample for each group of sites that are the same locus at several places in the reference: the informative pairs
    of the members are summed and called with the sexloci model (call_from_counts), so a combined call equals what a single site with those counts
    would give. Adds <group>_GT, _ref, _alt, _GQ, _reads and _status. A sample with no usable member gets an empty call and NaN counts."""
    for g, members in groups.items():
        def num(col):
            return pd.to_numeric(full[col], errors="coerce").fillna(0) if col in full.columns else pd.Series(0.0, index=full.index)
        n_ref, n_alt, n_reads = (sum(num(f"{m}_{k}") for m in members) for k in ("ref", "alt", "reads"))
        sts = pd.concat([full[f"{m}_status"].fillna("") if f"{m}_status" in full.columns else pd.Series("", index=full.index) for m in members], axis=1)
        ok = (sts == "OK").any(axis=1)
        first_bad = sts.apply(lambda r: next((v for v in r.tolist() if v and v != "OK"), ""), axis=1)
        calls = [call_from_counts(int(r_), int(a_), min_informative_pairs=min_informative_pairs) if o else ("", float("nan"), None)
                 for r_, a_, o in zip(n_ref, n_alt, ok)]
        full[f"{g}_GT"] = [c[0] for c in calls]
        full[f"{g}_ref"], full[f"{g}_alt"], full[f"{g}_reads"] = n_ref.where(ok), n_alt.where(ok), n_reads.where(ok)
        full[f"{g}_GQ"] = [round(c[1], 1) if c[1] == c[1] else np.nan for c in calls]
        full[f"{g}_status"] = np.where(ok, "OK", first_bad)
    return full


# ------------------------------------------------------------------------------------------------ reference windows
class _Window:
    def __init__(self, seq: str, start: int):
        self.seq, self.start = seq, start

    def __getitem__(self, sl: slice) -> str:
        a, b = sl.start, sl.stop
        a = self.start if a is None else a
        b = self.start + len(self.seq) if b is None else b
        a, b = max(a, self.start), min(b, self.start + len(self.seq))      # clipped: callers check the length
        return self.seq[a - self.start:b - self.start] if b > a else ""


class RefWindows:
    """Just the reference windows that are needed, indexed like pysam: ref[chrom][a:b] with absolute 0-based coordinates."""

    def __init__(self):
        self._w: Dict[str, _Window] = {}

    def __getitem__(self, chrom: str) -> _Window:
        return self._w[chrom]


def load_reference_windows(fasta_path: str, wins: Dict[str, Tuple[int, int]]) -> RefWindows:
    """Stream the FASTA once and keep wins[chrom] = (start0, end0) of each wanted chromosome."""
    op = gzip.open if fasta_path.endswith(".gz") else open
    out, parts, cur, pos = RefWindows(), {}, None, 0
    remaining = set(wins)
    with op(fasta_path, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                if cur is not None and cur in wins:
                    remaining.discard(cur)
                    if not remaining:
                        break
                cur, pos = line[1:].split()[0], 0
                continue
            if cur in wins:
                seq = line.rstrip("\n")
                lo, hi = wins[cur]
                a, b = max(lo - pos, 0), min(hi - pos, len(seq))
                if b > a:
                    parts.setdefault(cur, []).append(seq[a:b])
                pos += len(seq)
    for chrom, (lo, hi) in wins.items():
        if chrom not in parts:
            raise ValueError(f"chromosome {chrom} (window {lo}-{hi}) not found in {fasta_path}")
        out._w[chrom] = _Window("".join(parts[chrom]).upper(), lo)
    return out


def find_reference_fasta(genome_dir: str, explicit: Optional[str] = None) -> str:
    """The genome FASTA. An NCBI download folder also holds rna.fna, cds_from_genomic.fna, protein files and so on; those are ignored, and when
    more than one whole-genome FASTA remains the one named *_genomic.fna wins. Otherwise ask for --reference-fasta."""
    if explicit:
        if not os.path.exists(explicit):
            raise FileNotFoundError(f"reference FASTA not found: {explicit}")
        return explicit
    import glob
    allc = sorted({p for ext in ("*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz", "*.fasta.gz") for p in glob.glob(os.path.join(genome_dir, ext))})
    not_genome = ("rna", "cds", "protein", "translated", "pseudo")
    cands = [p for p in allc if not any(w in os.path.basename(p).lower() for w in not_genome)]
    if len(cands) > 1:
        g = [p for p in cands if "_genomic." in os.path.basename(p).lower()]
        cands = g if len(g) == 1 else cands
    if len(cands) == 1:
        return cands[0]
    names = ", ".join(map(os.path.basename, allc))
    raise FileNotFoundError(("no genome FASTA found" if not allc else f"cannot tell which of these is the genome: {names}") + f" in {genome_dir}; give --reference-fasta")


def build_classifiers(sites: List[IndelSite], ref: RefWindows, flanking: int = 250, backend="auto", **kw) -> List[IndelReadClassifier]:
    """One classifier per site. Checks that the REF allele of the VCF really is the sequence of this reference at that position."""
    be = make_backend(backend, kw.get("gap_open", 5), kw.get("gap_extend", 2)) if isinstance(backend, str) else backend
    out = []
    for s in sites:
        found = ref[s.chrom][s.pos - 1: s.pos - 1 + len(s.ref)].upper()
        if found != s.ref.upper():
            hint = ""
            for d in (-1, 1):
                if ref[s.chrom][s.pos - 1 + d: s.pos - 1 + d + len(s.ref)].upper() == s.ref.upper():
                    hint = f" (it matches at {s.chrom}:{s.pos + d}: off by one?)"
            raise ValueError(f"site {s.name}: the VCF REF {s.ref!r} does not match the reference at {s.chrom}:{s.pos} (found {found!r}){hint}. "
                             f"Was the VCF made against this assembly?")
        out.append(IndelReadClassifier.from_vcf_record(ref, s.chrom, s.pos, s.ref, s.alt, flanking=flanking, backend=be, **kw))
    return out


def merge_regions(wins: List[Tuple[str, int, int]]) -> List[str]:
    """(chrom, start0, end0) windows -> merged samtools regions 'chrom:start-end' (1-based, inclusive)."""
    out: List[List[Any]] = []
    for chrom, a, b in sorted(wins):
        if out and out[-1][0] == chrom and a <= out[-1][2]:
            out[-1][2] = max(out[-1][2], b)
        else:
            out.append([chrom, a, b])
    return [f"{c}:{a + 1}-{b}" for c, a, b in out]


# ------------------------------------------------------------------------------------------------ one sample
def reads_from_region(run: Runner, bam: str, region: str, min_mapq: int = 0, include_unmapped: bool = True):
    rc, out, err = run(["samtools", "view", "-F", "3328", bam, region])
    if rc != 0:
        raise RuntimeError(f"samtools view {region} failed: {err.strip()[:150]}")
    return reads_from_samtools_lines(out.splitlines(), min_mapq, include_unmapped)


def genotype_site(clf: IndelReadClassifier, reads, min_informative_pairs: int = 5) -> Dict[str, Any]:
    pairs = clf.classify_read_pairs(reads)
    call = clf.call_genotype(pairs, min_informative_pairs=min_informative_pairs)
    return dict(GT=call.genotype, n_ref=call.n_ref, n_alt=call.n_alt, n_equal=call.n_equal, n_uninf=call.n_uninformative,
                GQ=round(call.quality, 1), n_reads=len(reads))


def process_sample_sv(sample: str, base_url: str, workdir: str, keep_bam: bool, run: Runner, need_gsdf: bool,
                      need_sites: List[Tuple[IndelSite, IndelReadClassifier]], retries: int = 2, min_informative_pairs: int = 5) -> Dict:
    res: Dict[str, Any] = dict(SampleID=sample, fetch="error", message="", gsdf=None, te={})
    wins = [(c.contig,) + c.fetch_window() for _, c in need_sites]
    regions = ([f"{GSDF_CHROM}:{GSDF_REGION[0]}-{GSDF_REGION[1]}"] if need_gsdf else []) + merge_regions(wins)
    raw, bam = os.path.join(workdir, f"{sample}_sv.raw.bam"), os.path.join(workdir, f"{sample}_sv.bam")
    url = f"{base_url}/{sample}/{sample}.all.bam"
    err = ""
    for attempt in range(retries + 1):
        rc, _, err = run(["samtools", "view", "-b", "-o", raw, url] + regions)
        if rc == 0:
            break
        kind = _classify_failure(err)
        if kind != "error" or attempt == retries:
            res.update(fetch=kind, message=err.strip().splitlines()[-1][:160] if err.strip() else "")
            return res
        time.sleep(2 * (attempt + 1))
    try:
        for cmd in (["samtools", "sort", "-o", bam, raw], ["samtools", "index", bam]):
            rc, _, err = run(cmd)
            if rc != 0:
                res.update(fetch="error", message=f"{cmd[1]} failed: {err.strip()[:120]}")
                return res
        res["fetch"] = "OK"
        if need_gsdf:
            try:
                m = measure_bam(bam, run)
                m["Status"] = "OK" if np.isfinite(m["DepthRatio"]) else "no_coverage"
                m["Message"] = ""
                res["gsdf"] = m
            except Exception as e:
                res["gsdf"] = dict(Status="error", Message=str(e)[:160], FlankDepth=np.nan, UnitDepth=np.nan, DepthRatio=np.nan, JunctionJ=np.nan)
        for site, clf in need_sites:
            try:
                reads = reads_from_region(run, bam, f"{site.chrom}:{clf.fetch_window()[0] + 1}-{clf.fetch_window()[1]}")
                r = genotype_site(clf, reads, min_informative_pairs)
                r.update(Status="OK", Message="")
            except Exception as e:
                r = dict(Status="error", Message=str(e)[:160], GT="", n_ref=0, n_alt=0, n_equal=0, n_uninf=0, GQ=np.nan, n_reads=0)
            res["te"][site.name] = r
    finally:
        for p in (raw, raw + ".bai", bam + ".bai") + (() if keep_bam else (bam,)):
            if os.path.exists(p):
                os.remove(p)
    return res


# ------------------------------------------------------------------------------------------------ the whole database
def _fmt(v) -> str:
    return "" if (isinstance(v, float) and np.isnan(v)) else str(v)


def run_sv_genotype(tables, cfg, outdir: str, lab, genome_version: str, sites: List[IndelSite], classifiers: List[IndelReadClassifier],
                    do_gsdf: bool = True, jobs: int = 4, species: Optional[str] = None, ecogroup: Optional[str] = None,
                    samples: Optional[List[str]] = None, limit: Optional[int] = None, min_coverage: float = 2.0, refresh: bool = False,
                    keep_bams: bool = True, anchors: Optional[str] = None, theory: bool = False, port: Optional[int] = None,
                    min_informative_pairs: int = 5, min_gq: float = 10.0, gsdf_cache: Optional[str] = None,
                    runner: Optional[Runner] = None, server_factory: Optional[Callable] = None, log: Callable = print,
                    extra_columns: Optional[List[str]] = None) -> pd.DataFrame:
    os.makedirs(outdir, exist_ok=True)
    bamdir = os.path.join(outdir, "region_bams")
    os.makedirs(bamdir, exist_ok=True)
    run = runner or default_runner
    gcache, tcache = os.path.join(outdir, "gsdf_measurements.tsv"), os.path.join(outdir, "te_measurements.tsv")
    by_name = {s.name: (s, c) for s, c in zip(sites, classifiers)}

    # ---- who is genotyped
    extra = resolve_extra_columns(tables.sample, extra_columns, log)
    meta = tables.sample.copy()
    for c in ("Species", "Ecogroup", "Sex"):
        if c not in meta.columns:
            meta[c] = ""
    meta = meta.drop_duplicates("SampleID").set_index("SampleID")[["Species", "Ecogroup", "Sex"] + list(extra.values())].fillna("")
    meta = meta.rename(columns={v: k for k, v in extra.items()})
    meta["Pedigree"] = pedigree_roles(tables.sample).reindex(meta.index).fillna("")
    cov, aligned = pd.Series(dtype=float), set()
    if tables.alignment is not None:
        al = tables.alignment.drop_duplicates("SampleID").set_index("SampleID")
        aligned = set(al.index)
        if "Coverage" in al.columns:
            cov = pd.to_numeric(al["Coverage"], errors="coerce")
    pick = meta
    if species:
        pick = pick[pick.Species.str.contains(species, case=False, regex=False)]
    if ecogroup:
        pick = pick[pick.Ecogroup.str.lower() == ecogroup.lower()]
    if samples:
        pick = pick[[any(fnmatch(s, p) for p in samples) for s in pick.index]]
    todo = [s for s in pick.index if s in aligned and (min_coverage <= 0 or not (s in cov.index and cov[s] < min_coverage))]
    if limit:
        todo = todo[:limit]

    # ---- what is already measured (per locus). --refresh measures the SELECTED samples again and keeps everyone else's results.
    if do_gsdf and not os.path.exists(gcache) and gsdf_cache and os.path.exists(gsdf_cache):
        old = pd.read_csv(gsdf_cache, sep="\t")
        old = old[old.Status == "OK"]
        old.to_csv(gcache, sep="\t", index=False)
        log(f"[sv] imported {len(old)} finished gsdf measurements from {gsdf_cache}")
    gdone = pd.read_csv(gcache, sep="\t") if (do_gsdf and os.path.exists(gcache)) else pd.DataFrame(columns=MEAS_COLS)
    tdone = pd.read_csv(tcache, sep="\t", dtype={"GT": str, "Site": str}).fillna({"GT": ""}) if os.path.exists(tcache) else pd.DataFrame(columns=TE_COLS)
    if refresh:
        sel = set(todo)
        gdone, tdone = gdone[~gdone.SampleID.isin(sel)], tdone[~tdone.SampleID.isin(sel)]
        log(f"[sv] --refresh: the {len(sel)} selected samples are measured again; results of the other samples are kept")
    g_have = set(gdone[gdone.Status == "OK"].SampleID) if len(gdone) else set()
    t_have = {(r.SampleID, r.Site) for r in tdone[tdone.Status == "OK"].itertuples()} if len(tdone) else set()

    jobs_list = []
    for s in todo:
        need_g = do_gsdf and s not in g_have
        need_s = [by_name[n] for n in by_name if (s, n) not in t_have]
        if need_g or need_s:
            jobs_list.append((s, need_g, need_s))
    log(f"[sv] {len(meta)} samples in the database; {len(pick)} selected; {len(todo)} have a {genome_version} alignment "
        f"(coverage >= {min_coverage}x); {len(jobs_list)} need at least one locus fetched "
        f"(gsdf: {sum(j[1] for j in jobs_list)}, insertion sites: {sum(len(j[2]) for j in jobs_list)} sample-sites)")

    g_rows, t_rows = gdone.to_dict("records"), tdone.to_dict("records")
    if jobs_list:
        remote = lab.cloud_master_dir + f"Bamfiles/{genome_version}/"
        factory = server_factory or (lambda: HttpServer(remote, port))
        lock = threading.Lock()
        if do_gsdf:                                                    # rewrite the caches without the samples about to be measured again
            gdone[MEAS_COLS].to_csv(gcache, sep="\t", index=False) if len(gdone) else open(gcache, "w").write("\t".join(MEAS_COLS) + "\n")
        tdone[TE_COLS].to_csv(tcache, sep="\t", index=False) if len(tdone) else open(tcache, "w").write("\t".join(TE_COLS) + "\n")
        with factory() as base:
            with ThreadPoolExecutor(max_workers=max(1, jobs)) as ex:
                futs = {ex.submit(process_sample_sv, s, base, bamdir, keep_bams, run, ng, ns, 2, min_informative_pairs): s for s, ng, ns in jobs_list}
                need_by = {j[0]: (j[1], j[2]) for j in jobs_list}
                for i, f in enumerate(as_completed(futs), 1):
                    r = f.result()
                    sid, tag = r["SampleID"], []
                    ng, ns = need_by[sid]
                    with lock:
                        if r["fetch"] != "OK":                      # record the failure for every locus that was needed, so it is retried
                            tag.append(r["fetch"] + (f" ({r['message']})" if r["message"] else ""))
                            if ng:
                                g = dict(SampleID=sid, Status=r["fetch"], FlankDepth=np.nan, UnitDepth=np.nan, DepthRatio=np.nan, JunctionJ=np.nan,
                                         Message=r["message"])
                                g_rows.append(g)
                                with open(gcache, "a") as fh:
                                    fh.write("\t".join(_fmt(g.get(c, "")) for c in MEAS_COLS) + "\n")
                            for site, _c in ns:
                                row = dict(SampleID=sid, Site=site.name, Status=r["fetch"], GT="", n_ref=0, n_alt=0, n_equal=0, n_uninf=0, GQ=np.nan,
                                           n_reads=0, Message=r["message"])
                                t_rows.append(row)
                                with open(tcache, "a") as fh:
                                    fh.write("\t".join(_fmt(row.get(c, "")) for c in TE_COLS) + "\n")
                        else:
                            if r["gsdf"] is not None:
                                g = dict(SampleID=sid, **r["gsdf"])
                                g_rows.append(g)
                                with open(gcache, "a") as fh:
                                    fh.write("\t".join(_fmt(g.get(c, "")) for c in MEAS_COLS) + "\n")
                                tag.append(f"gsdf ratio {g['DepthRatio']:.2f}" if g["Status"] == "OK" else f"gsdf {g['Status']}")
                            for name, t in r["te"].items():
                                row = dict(SampleID=sid, Site=name, **t)
                                t_rows.append(row)
                                with open(tcache, "a") as fh:
                                    fh.write("\t".join(_fmt(row.get(c, "")) for c in TE_COLS) + "\n")
                                tag.append(f"{name} {t['GT'] or t['Status']}")
                    log(f"[sv] [{i}/{len(jobs_list)}] {sid}: " + "; ".join(tag))

    # ---- gsdf calls (as gsdf-genotype)
    meas = pd.DataFrame(g_rows, columns=MEAS_COLS).drop_duplicates("SampleID", keep="last").set_index("SampleID") if do_gsdf else None
    gcols = ["Status", "FlankDepth", "DepthRatio", "JunctionJ", "EstCopies", "CopyNumber", "Genotype", "Flag"]
    if do_gsdf and len(meas):
        ok = meas[meas.Status == "OK"]
        if theory:
            r2, step, how = 0.50, 0.25, "theory (copies / 4)"
        elif anchors:
            adf = pd.read_csv(anchors, sep="\t")
            a = adf.set_index("SampleID" if "SampleID" in adf.columns else "sample").copies
            r2, step = fit_anchors(ok, a)
            how = f"fitted on {int(a.index.isin(ok.index).sum())} anchor samples"
        else:
            r2, step, how = getattr(cfg, "gsdf_r2", DEFAULT_R2), getattr(cfg, "gsdf_step", DEFAULT_STEP), "MC family (daughters 2 copies, sons 3)"
        log(f"[sv] gsdf calibration: {how} (ratio at 2 copies {r2:.3f}, +{step:.3f} per extra copy)")
        called = call_genotypes(meas, r2, step)[gcols].rename(columns={"Status": "gsdf_status"})
    else:
        called = pd.DataFrame(columns=[c if c != "Status" else "gsdf_status" for c in gcols])
    full = meta.join(called, how="left")
    full["Genotype"] = full["Genotype"].fillna("") if "Genotype" in full else ""
    full["Flag"] = full["Flag"].fillna("") if "Flag" in full else ""

    # ---- insertion calls
    tt = pd.DataFrame(t_rows, columns=TE_COLS).drop_duplicates(["SampleID", "Site"], keep="last")
    for s in sites:
        x = tt[tt.Site == s.name].set_index("SampleID")
        for col, src in ((f"{s.name}_GT", "GT"), (f"{s.name}_ref", "n_ref"), (f"{s.name}_alt", "n_alt"), (f"{s.name}_GQ", "GQ"),
                         (f"{s.name}_reads", "n_reads"), (f"{s.name}_status", "Status")):
            full[col] = x[src].reindex(full.index) if len(x) else np.nan
        full[f"{s.name}_GT"] = full[f"{s.name}_GT"].fillna("")
    groups = site_groups(sites)
    add_group_columns(full, groups, min_informative_pairs)
    full["Coverage"] = cov.reindex(full.index)
    st_cols = (["gsdf_status"] if do_gsdf and "gsdf_status" in full.columns else []) + [f"{x.name}_status" for x in sites]
    overall = []
    for sid, row in full[st_cols].iterrows():
        sts = [v for v in row.tolist() if isinstance(v, str) and v]
        bad = [v for v in sts if v != "OK"]
        overall.append(None if not sts else (bad[0] if bad else "OK"))
    full["Status"] = pd.Series(overall, index=full.index)
    full["Status"] = full.Status.where(full.Status.notna(), np.where(full.index.isin(aligned), "not selected / skipped", f"no {genome_version} alignment"))
    full.loc[full.index.isin(aligned) & (full.Status == "not selected / skipped") & (full.Coverage < min_coverage), "Status"] = f"skipped: coverage < {min_coverage}x"
    full = full.reset_index().rename(columns={"index": "SampleID"})
    site_cols = [c for s in sites for c in (f"{s.name}_GT", f"{s.name}_ref", f"{s.name}_alt", f"{s.name}_GQ", f"{s.name}_reads", f"{s.name}_status")]
    site_cols += [c for g in groups for c in (f"{g}_GT", f"{g}_ref", f"{g}_alt", f"{g}_GQ", f"{g}_reads", f"{g}_status")]
    cols = (["SampleID", "Species", "Ecogroup"] + list(extra) + ["Sex", "Pedigree"] + ["Genotype", "CopyNumber", "EstCopies", "DepthRatio", "FlankDepth", "JunctionJ", "Flag"]
            + site_cols + ["Coverage", "Status"])
    for c in cols:
        if c not in full.columns:
            full[c] = ""
    full = full[cols].sort_values(["Ecogroup"] + list(extra) + ["Species", "SampleID"]).reset_index(drop=True)
    path = os.path.join(outdir, "sv_genotypes.tsv")
    full.round({"EstCopies": 2, "DepthRatio": 3, "FlankDepth": 1, "JunctionJ": 3, "Coverage": 1}).to_csv(path, sep="\t", index=False)

    # ---- summary
    for s in sites:
        called_s = full[full[f"{s.name}_GT"].isin(["0/0", "0/1", "1/1"]) & (pd.to_numeric(full[f"{s.name}_GQ"], errors="coerce") >= min_gq)]
        log(f"\n[sv] {s.name} ({s.chrom}:{s.pos}, +{s.length} bp): {len(called_s)} confident calls (GQ >= {min_gq:g}) of "
            f"{int((full[f'{s.name}_status'] == 'OK').sum())} genotyped")
        if len(called_s):
            log(pd.crosstab(called_s.Ecogroup.replace("", "(none)"), called_s[f"{s.name}_GT"]).to_string())
    for g, members in groups.items():
        called_g = full[full[f"{g}_GT"].isin(["0/0", "0/1", "1/1"]) & (pd.to_numeric(full[f"{g}_GQ"], errors="coerce") >= min_gq)]
        log(f"\n[sv] {g} (combined: {' + '.join(members)}, informative pairs summed): {len(called_g)} confident calls (GQ >= {min_gq:g}) of "
            f"{int((full[f'{g}_status'] == 'OK').sum())} genotyped")
        if len(called_g):
            log(pd.crosstab(called_g.Ecogroup.replace("", "(none)"), called_g[f"{g}_GT"]).to_string())
    log(f"[sv] wrote {path}")
    return full


def sheet_specs(full: pd.DataFrame, sites: List[IndelSite], gsdf_column: str = "gsdf number", do_gsdf: bool = True, include_flagged: bool = False,
                min_gq: float = 10.0, include_lowqual: bool = False, version: str = "") -> List[Dict[str, Any]]:
    """What to write into the Sample Database: gsdf copy number (unflagged calls) and, per insertion site, the genotype (GQ >= min_gq)."""
    specs = []
    if do_gsdf:
        g = full[full.Genotype.fillna("").astype(str) != ""]
        g = g if include_flagged else g[g.Flag.fillna("") == ""]
        specs.append(dict(column=gsdf_column, values={str(r.SampleID).strip(): int(r.CopyNumber) for r in g.itertuples()},
                          note=(f"Total gsdf copies in the diploid genome (2 = no duplicated chromosome, 3 = one, 4 = two), from short-read depth over "
                                f"the 21.7 kb gsdf unit relative to its flanks; sexloci sv-genotype {version}. Blank = no call or a flagged call.")))
    groups = site_groups(sites)
    in_group = {m for ms in groups.values() for m in ms}                 # the per-copy columns stay in the TSV; the sheet gets the combined call
    for s in sites:
        if s.name in in_group:
            continue
        gt = full[full[f"{s.name}_GT"].isin(["0/0", "0/1", "1/1"])]
        if not include_lowqual:
            gt = gt[pd.to_numeric(gt[f"{s.name}_GQ"], errors="coerce") >= min_gq]
        specs.append(dict(column=s.name, values={str(r.SampleID).strip(): getattr(r, f"{s.name}_GT") for r in gt.itertuples()},
                          note=(f"Genotype of the {s.length} bp insertion at {s.chrom}:{s.pos} from short reads (0/0 absent, 0/1 heterozygous, 1/1 "
                                f"homozygous); read pairs crossing the junctions are scored against the reference and insertion alleles; "
                                f"sexloci sv-genotype {version}. Blank = fewer than the minimum informative pairs or GQ < {min_gq:g}.")))
    from .known_sites import SITE_GROUP_NOTES
    for g, members in groups.items():
        gt = full[full[f"{g}_GT"].isin(["0/0", "0/1", "1/1"])]
        if not include_lowqual:
            gt = gt[pd.to_numeric(gt[f"{g}_GQ"], errors="coerce") >= min_gq]
        base = SITE_GROUP_NOTES.get(g, f"Combined genotype of {', '.join(members)} (informative read pairs summed).")
        specs.append(dict(column=g, values={str(r.SampleID).strip(): getattr(r, f"{g}_GT") for r in gt.itertuples()},
                          note=f"{base} sexloci sv-genotype {version}. Blank = fewer than the minimum informative pairs or GQ < {min_gq:g}."))
    return specs
