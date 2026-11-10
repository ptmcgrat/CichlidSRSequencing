"""Genotype the gsdf-region tandem duplication from short reads, for every sample with an alignment to the GT3a reference.

The M. zebra GT3a reference carries two tandem copies of a 21.7 kb unit (NC_135173.1:50,751,322-50,794,777). Reads from every genomic
copy map to both reference copies, so a sample's depth over the unit relative to its own flanks is proportional to its total copy number:

    ratio = (copies in the diploid genome) / 4        (2 copies -> ~0.5, 3 -> ~0.75, 4 -> ~1.0, the reference's own state)

For each sample only the 130 kb around the locus is read from the cloud BAM (a few hundred kB), through a read-only rclone web
server on the Bamfiles folder; nothing is written to the cloud and the whole BAM is never downloaded.

Genotype notation (assumes each chromosome carries one or two copies of the unit):
    -/-      two copies in total (both chromosomes single-copy)
    dup/-    three copies (one chromosome carries the tandem duplication, e.g. a male whose Y does)
    dup/dup  four copies (both chromosomes do, like the reference)
    CN1, CN5, CN6   anything else (check the Flag column)
Depth alone cannot say WHICH chromosome carries the duplication.
"""
from __future__ import annotations

import io
import os
import socket
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from fnmatch import fnmatch
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

CHROM = "NC_135173.1"
REGION = (50_700_000, 50_830_000)
UNIT = (50_755_000, 50_785_000)                         # 5 kb bins lying fully inside the two copies
FLANK = ((0, 50_735_000), (50_800_000, 50_815_000))     # outside the block and before the repeat-rich bump at 50.825 Mb
JUNCTION = (50_772_650, 50_773_150)                     # where the two copies join (unique reads only)
BIN = 5000
DEFAULT_R2, DEFAULT_STEP = 0.530, 0.268                 # calibrated on the MC family: 2 copies -> 0.530, +0.268 per extra copy
MIN_FLANK_DEPTH = 3.0

MEAS_COLS = ["SampleID", "Status", "FlankDepth", "UnitDepth", "DepthRatio", "JunctionJ", "Message"]
OUT_COLS = ["SampleID", "Species", "Ecogroup", "Sex", "Pedigree", "Genotype", "CopyNumber", "EstCopies", "DepthRatio", "FlankDepth", "JunctionJ",
            "Flag", "Coverage", "Status"]

Runner = Callable[[List[str]], Tuple[int, str, str]]


def default_runner(args: List[str], timeout: int = 600) -> Tuple[int, str, str]:
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError:
        raise RuntimeError(f"'{args[0]}' is not installed or not on PATH") from None
    except subprocess.TimeoutExpired:
        return 124, "", f"timeout after {timeout}s"
    return p.returncode, p.stdout, p.stderr


# ------------------------------------------------------------------------------------------------ reading the region
class HttpServer:
    """A read-only rclone web server on the cloud Bamfiles folder (the only network step)."""

    def __init__(self, remote_path: str, port: Optional[int] = None):
        self.remote_path, self.port, self.proc = remote_path, port, None

    def __enter__(self):
        if self.port is None:
            with socket.socket() as s:
                s.bind(("127.0.0.1", 0))
                self.port = s.getsockname()[1]
        try:
            self.proc = subprocess.Popen(["rclone", "serve", "http", "--read-only", "--addr", f"127.0.0.1:{self.port}", self.remote_path],
                                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except FileNotFoundError:
            raise RuntimeError("rclone is not installed or not on PATH") from None
        for _ in range(100):
            if self.proc.poll() is not None:
                raise RuntimeError(f"rclone serve exited immediately; check that {self.remote_path} exists (rclone lsf '{self.remote_path}')")
            with socket.socket() as s:
                if s.connect_ex(("127.0.0.1", self.port)) == 0:
                    return f"http://127.0.0.1:{self.port}"
            time.sleep(0.3)
        raise RuntimeError("rclone serve did not start within 30 s")

    def __exit__(self, *exc):
        if self.proc is not None:
            self.proc.terminate()
            try:
                self.proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.proc.kill()


def _classify_failure(err: str) -> str:
    e = err.lower()
    if "404" in e or "not found" in e or "no such file" in e:
        return "no_bam"
    if "index" in e or ".bai" in e or ".csi" in e:
        return "no_index"
    return "error"


def measure_bam(bam: str, run: Runner) -> Dict[str, float]:
    """Depth over the unit relative to the flanks (5 kb bins, duplicates excluded as samtools depth does by default) and the unique
    (MAPQ>=30) coverage at the junction between the copies, relative to the unique flank depth."""
    reg = f"{CHROM}:{REGION[0]}-{REGION[1]}"
    rc, out, err = run(["samtools", "depth", "-a", "-r", reg, bam])
    if rc != 0:
        raise RuntimeError(f"samtools depth failed: {err.strip()[:150]}")
    d = pd.read_csv(io.StringIO(out), sep="\t", header=None, names=["c", "pos", "d"], usecols=[1, 2])
    bins = d.groupby(d.pos // BIN * BIN).d.mean()
    unit = bins[(bins.index >= UNIT[0]) & (bins.index <= UNIT[1])].mean()
    flank = bins[[any(lo <= b <= hi for lo, hi in FLANK) for b in bins.index]].mean()
    res = dict(FlankDepth=flank, UnitDepth=unit, DepthRatio=unit / flank if flank and flank > 0 else np.nan, JunctionJ=np.nan)
    try:
        def q30(lo, hi):
            rc2, o2, _ = run(["samtools", "depth", "-a", "-Q", "30", "-r", f"{CHROM}:{lo}-{hi}", bam])
            if rc2 != 0 or not o2.strip():
                return np.nan
            return float(pd.read_csv(io.StringIO(o2), sep="\t", header=None, usecols=[2])[2].mean())
        fq, jq = q30(REGION[0], FLANK[0][1]), q30(*JUNCTION)
        if fq and fq > 0:
            res["JunctionJ"] = jq / fq
    except Exception:
        pass
    return res


def process_sample(sample: str, base_url: str, workdir: str, keep_bam: bool, run: Runner, retries: int = 2) -> Dict:
    row = dict(SampleID=sample, Status="error", FlankDepth=np.nan, UnitDepth=np.nan, DepthRatio=np.nan, JunctionJ=np.nan, Message="")
    bam = os.path.join(workdir, f"{sample}_gsdf.bam")
    reg = f"{CHROM}:{REGION[0]}-{REGION[1]}"
    url = f"{base_url}/{sample}/{sample}.all.bam"
    err = ""
    for attempt in range(retries + 1):
        rc, _, err = run(["samtools", "view", "-b", "-o", bam, url, reg])
        if rc == 0:
            break
        kind = _classify_failure(err)
        if kind != "error" or attempt == retries:
            row.update(Status=kind, Message=err.strip().splitlines()[-1][:160] if err.strip() else "")
            return row
        time.sleep(2 * (attempt + 1))
    try:
        rc, _, err = run(["samtools", "index", bam])
        if rc != 0:
            row.update(Status="error", Message="index failed: " + err.strip()[:120])
            return row
        row.update(measure_bam(bam, run))
        row["Status"] = "OK" if np.isfinite(row["DepthRatio"]) else "no_coverage"
    except Exception as e:
        row.update(Status="error", Message=str(e)[:160])
    finally:
        if not keep_bam:
            for p in (bam, bam + ".bai"):
                if os.path.exists(p):
                    os.remove(p)
    return row


# ------------------------------------------------------------------------------------------------ calling
def genotype_label(n) -> str:
    if n is None or pd.isna(n):
        return ""
    n = int(n)
    return {2: "-/-", 3: "dup/-", 4: "dup/dup"}.get(n, f"CN{n}")


def fit_anchors(meas: pd.DataFrame, anchors: pd.Series) -> Tuple[float, float]:
    a = anchors[anchors.index.isin(meas.index)]
    if a.nunique() < 2:
        raise ValueError("--anchors needs samples with at least two different copy numbers that were measured")
    step, icpt = np.polyfit(a.values.astype(float), meas.loc[a.index, "DepthRatio"].values.astype(float), 1)
    return float(icpt + 2 * step), float(step)


def call_genotypes(m: pd.DataFrame, r2: float, step: float, min_flank_depth: float = MIN_FLANK_DEPTH) -> pd.DataFrame:
    m = m.copy()
    m["EstCopies"] = 2 + (m.DepthRatio - r2) / step
    n = m.EstCopies.round().clip(1, 6)
    m["CopyNumber"] = n.astype("Int64")
    m["Genotype"] = [genotype_label(x) for x in m.CopyNumber]
    flags = []
    for _, r in m.iterrows():
        f = []
        if pd.notna(r.FlankDepth) and r.FlankDepth < min_flank_depth:
            f.append("low depth")
        if pd.notna(r.EstCopies):
            if abs(r.EstCopies - round(r.EstCopies)) > 0.35:
                f.append("uncertain")
            if r.EstCopies < 0.6 or r.EstCopies > 5.4:
                f.append("outside calibrated range")
        flags.append("; ".join(f))
    m["Flag"] = flags
    return m


# ------------------------------------------------------------------------------------------------ extra sample-sheet columns
def resolve_extra_columns(sample: pd.DataFrame, requested: Optional[List[str]], log: Callable = print) -> Dict[str, str]:
    """output column name -> sample-sheet column. 'Subgroup' is picked up automatically (any spelling: SubGroup, sub_group, ...)."""
    norm = lambda x: "".join(ch for ch in str(x).lower() if ch.isalnum())
    cols = {norm(c): c for c in sample.columns}
    out: Dict[str, str] = {}
    if "subgroup" in cols:
        out["Subgroup"] = cols["subgroup"]
    for want in requested or []:
        if norm(want) in cols:
            out[cols[norm(want)]] = cols[norm(want)]
        else:
            log(f"[gsdf] WARNING: the sample sheet has no column '{want}' (columns: {', '.join(map(str, sample.columns))})")
    return {k: v for k, v in out.items() if k not in ("SampleID", "Species", "Ecogroup", "Sex")}


_UNKNOWN_PARENT = {"", "?", "na", "nan", "none", "unknown", "unk", "-"}


def pedigree_roles(sample: pd.DataFrame) -> pd.Series:
    """SampleID -> 'offspring' (has a Parents entry) or 'parent' (named in someone's Parents entry); other samples are absent.
    Fish from the same family share haplotypes, so they must not be counted as independent replicates."""
    roles: Dict[str, str] = {}
    if "Parents" not in sample.columns:
        return pd.Series(dtype=object)
    for sid, par in zip(sample.SampleID.astype(str), sample.Parents):
        toks = [t.strip() for t in str(par).split(",")] if pd.notna(par) else []
        toks = [t for t in toks if t.lower() not in _UNKNOWN_PARENT]
        if toks:
            roles[sid] = "offspring"
            for t in toks:
                roles.setdefault(t, "parent")
    return pd.Series(roles, dtype=object)


# ------------------------------------------------------------------------------------------------ the whole database
def run_gsdf_genotype(tables, cfg, outdir: str, lab, genome_version: str, jobs: int = 4, species: Optional[str] = None,
                      ecogroup: Optional[str] = None, samples: Optional[List[str]] = None, limit: Optional[int] = None,
                      min_coverage: float = 2.0, refresh: bool = False, keep_bams: bool = True, anchors: Optional[str] = None,
                      theory: bool = False, port: Optional[int] = None, runner: Optional[Runner] = None,
                      server_factory: Optional[Callable] = None, log: Callable = print,
                      extra_columns: Optional[List[str]] = None) -> pd.DataFrame:
    os.makedirs(outdir, exist_ok=True)
    bamdir = os.path.join(outdir, "region_bams")
    os.makedirs(bamdir, exist_ok=True)
    run = runner or default_runner
    cache = os.path.join(outdir, "gsdf_measurements.tsv")

    extra = resolve_extra_columns(tables.sample, extra_columns, log)
    meta = tables.sample.copy()
    for c in ("Species", "Ecogroup", "Sex"):
        if c not in meta.columns:
            meta[c] = ""
    meta = meta.drop_duplicates("SampleID").set_index("SampleID")[["Species", "Ecogroup", "Sex"] + list(extra.values())].fillna("")
    meta = meta.rename(columns={v: k for k, v in extra.items()})
    meta["Pedigree"] = pedigree_roles(tables.sample).reindex(meta.index).fillna("")
    cov = pd.Series(dtype=float)
    aligned = set()
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

    done = pd.DataFrame(columns=MEAS_COLS)
    if os.path.exists(cache) and not refresh:
        done = pd.read_csv(cache, sep="\t")
    have = set(done[done.Status == "OK"].SampleID) if len(done) else set()
    todo_fetch = [s for s in todo if s not in have]
    log(f"[gsdf] {len(meta)} samples in the database; {len(pick)} selected; {len(todo)} have a {genome_version} alignment "
        f"(coverage >= {min_coverage}x); {len(have & set(todo))} already measured; {len(todo_fetch)} to fetch")

    results = [] if refresh else [r for r in done.to_dict("records") if r["SampleID"] not in set(todo_fetch)]
    if todo_fetch:
        remote = lab.cloud_master_dir + f"Bamfiles/{genome_version}/"
        factory = server_factory or (lambda: HttpServer(remote, port))
        lock = threading.Lock()
        if not os.path.exists(cache) or refresh:
            with open(cache, "w") as fh:
                fh.write("\t".join(MEAS_COLS) + "\n")
        with factory() as base:
            with ThreadPoolExecutor(max_workers=max(1, jobs)) as ex:
                futs = {ex.submit(process_sample, s, base, bamdir, keep_bams, run): s for s in todo_fetch}
                for i, f in enumerate(as_completed(futs), 1):
                    r = f.result()
                    with lock:
                        results.append(r)
                        with open(cache, "a") as fh:
                            fh.write("\t".join("" if (isinstance(r.get(c), float) and np.isnan(r.get(c))) else str(r.get(c, "")) for c in MEAS_COLS) + "\n")
                    tag = f"ratio {r['DepthRatio']:.2f}" if r["Status"] == "OK" else r["Status"] + (f" ({r['Message']})" if r["Message"] else "")
                    log(f"[gsdf] [{i}/{len(todo_fetch)}] {r['SampleID']}: {tag}")
    meas = pd.DataFrame(results, columns=MEAS_COLS).drop_duplicates("SampleID", keep="last").set_index("SampleID")
    nb = int((meas.Status == "no_bam").sum())
    if nb and nb >= 0.5 * max(len(meas), 1):
        log(f"[gsdf] WARNING: {nb} of {len(meas)} BAMs were not found; check the cloud layout: rclone lsf '{lab.cloud_master_dir}Bamfiles/{genome_version}/' | head")

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
    log(f"[gsdf] calibration: {how} (ratio at 2 copies {r2:.3f}, +{step:.3f} per extra copy)")

    called = call_genotypes(meas, r2, step)
    full = meta.join(called[["Status", "FlankDepth", "DepthRatio", "JunctionJ", "EstCopies", "CopyNumber", "Genotype", "Flag"]], how="left")
    full["Genotype"] = full.Genotype.fillna("")
    full["Flag"] = full.Flag.fillna("")
    full["Coverage"] = cov.reindex(full.index)
    full["Status"] = full.Status.where(full.Status.notna(), np.where(full.index.isin(aligned), "not selected / skipped", f"no {genome_version} alignment"))
    full.loc[full.index.isin(aligned) & full.Status.eq("not selected / skipped") & (full.Coverage < min_coverage), "Status"] = f"skipped: coverage < {min_coverage}x"
    full = full.reset_index().rename(columns={"index": "SampleID"})
    cols = ["SampleID", "Species", "Ecogroup"] + list(extra) + OUT_COLS[3:]
    full["Pedigree"] = full.Pedigree.fillna("")
    full = full[cols].sort_values(["Ecogroup"] + list(extra) + ["Species", "SampleID"]).reset_index(drop=True)
    path = os.path.join(outdir, "gsdf_genotypes.tsv")
    full.round({"EstCopies": 2, "DepthRatio": 3, "FlankDepth": 1, "JunctionJ": 3, "Coverage": 1}).to_csv(path, sep="\t", index=False)
    scored = full[full.Genotype.astype(str) != ""]
    if len(scored):
        log("\n[gsdf] genotypes by ecogroup (samples with a call):")
        log(pd.crosstab(scored.Ecogroup.replace("", "(none)"), scored.Genotype).to_string())
        for name in extra:
            if 1 < scored[name].nunique() <= 60:
                log(f"\n[gsdf] genotypes by {name}:")
                log(pd.crosstab(scored[name].replace("", "(none)"), scored.Genotype).to_string())
        log(f"[gsdf] flagged: {(scored.Flag != '').sum()} of {len(scored)}  |  statuses: " + ", ".join(f"{k}: {v}" for k, v in full.Status.value_counts().items()))
    log(f"[gsdf] wrote {path}")
    return full


# ------------------------------------------------------------------------------------------------ summaries by any column
def _fisher(a: int, b: int, c: int, d: int) -> float:
    from scipy.stats import fisher_exact
    return float(fisher_exact([[a, b], [c, d]])[1])


def summarize_by(table: pd.DataFrame, by: str, min_n: int = 1, independent: bool = True) -> pd.DataFrame:
    """Carrier frequency (>= 3 copies, i.e. at least one duplicated chromosome) per group, split by sex.

    pattern: absent (no carrier among >= 5 samples) | male-limited (>= 80% of >= 5 males carry it, <= 10% of >= 3 females) |
    male-biased (Fisher p < 0.05, >= 3 of each sex, more males than females carry it) | no sex information (< 3 males or < 3 females of
    known sex) | present, not sex-associated.
    independent=True drops fish that belong to a pedigree (Pedigree = offspring or parent): relatives share haplotypes, so one family
    can make a result look far stronger than it is. Needs the Pedigree column written by gsdf-genotype."""
    if by not in table.columns:
        raise ValueError(f"no column '{by}' in the table (columns: {', '.join(table.columns)})")
    t = table[table.Genotype.fillna("").astype(str) != ""].copy()
    if independent and "Pedigree" in t.columns:
        t = t[t.Pedigree.fillna("") == ""]
    t["carrier"] = t.CopyNumber >= 3
    rows = []
    for grp, x in t.groupby(t[by].fillna("").replace("", "(none)")):
        if len(x) < min_n:
            continue
        m, f = x[x.Sex == "M"], x[x.Sex == "F"]
        mc, fc = int(m.carrier.sum()), int(f.carrier.sum())
        p = _fisher(mc, len(m) - mc, fc, len(f) - fc) if len(m) >= 3 and len(f) >= 3 else np.nan
        if x.carrier.sum() == 0 and len(x) >= 5:
            pat = "absent"
        elif len(m) < 3 or len(f) < 3:
            pat = "no sex information"
        elif len(m) >= 5 and mc / len(m) >= 0.8 and fc / len(f) <= 0.1:
            pat = "male-limited"
        elif p < 0.05 and mc / len(m) > fc / len(f):
            pat = "male-biased"
        else:
            pat = "present, not sex-associated" if x.carrier.sum() else "absent"
        rows.append({by: grp, "n": len(x), "carriers": int(x.carrier.sum()), "pct_carriers": round(100 * x.carrier.mean(), 1),
                     "males": len(m), "male_carriers": mc, "females": len(f), "female_carriers": fc,
                     "unknown_sex": int((~x.Sex.isin(["M", "F"])).sum()), "unknown_sex_carriers": int(x[~x.Sex.isin(["M", "F"])].carrier.sum()),
                     "dup/dup": int((x.CopyNumber == 4).sum()), "5+ copies": int((x.CopyNumber >= 5).sum()),
                     "flagged": int((x.Flag.fillna("") != "").sum()), "male_vs_female_p": p, "pattern": pat})
    out = pd.DataFrame(rows)
    return out.sort_values(["carriers", "n"], ascending=False).reset_index(drop=True) if len(out) else out
