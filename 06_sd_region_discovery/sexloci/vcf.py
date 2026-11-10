"""Streaming VCF reader producing genotype 'chunks' (alt-allele dosage, depth).

Two backends:
  * builtin - dependency-free text parser (fine for small/medium data; slow for whole-genome VCFs)
  * allel   - scikit-allel's compiled parser (recommended for real data)

For large cohorts pre-filter with bcftools (e.g. ``bcftools view -m2 -M2 -v snps``) and/or run
per-chromosome with ``--contigs``.
"""
from __future__ import annotations

import gzip
import operator
import re
from collections import Counter
from dataclasses import dataclass
from typing import Iterator, List, Optional

import numpy as np

_GT = {"0/0": 0, "0|0": 0, "0/1": 1, "1/0": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}


@dataclass
class Chunk:
    chrom: np.ndarray          # (m,) str
    pos: np.ndarray            # (m,) int64
    dos: np.ndarray            # (m, n) int8: alt dosage 0/1/2, -1 = missing (after the strict min_dp/min_gq mask)
    dp: Optional[np.ndarray]   # (m, n) float32 raw depth (NaN if unknown) or None
    gq: Optional[np.ndarray] = None        # (m, n) float32 raw GQ or None
    dos_raw: Optional[np.ndarray] = None   # genotypes with no depth/GQ masking (only the GT field)
    dos_role: Optional[np.ndarray] = None  # genotypes masked with each sample's role-specific thresholds

    def subset(self, keep: np.ndarray) -> "Chunk":
        f = lambda a: None if a is None else a[keep]
        return Chunk(self.chrom[keep], self.pos[keep], self.dos[keep], f(self.dp), f(self.gq), f(self.dos_raw),
                     f(self.dos_role))


def mask_by_sample(dos_raw, dp, gq, min_dp, min_gq):
    """Mask genotypes with a per-sample depth / GQ threshold (vectors of length n_samples)."""
    out = dos_raw.copy()
    if dp is not None:
        out[np.nan_to_num(dp, nan=np.inf) < np.asarray(min_dp, float)[None, :]] = -1
    if gq is not None:
        out[np.nan_to_num(gq, nan=np.inf) < np.asarray(min_gq, float)[None, :]] = -1
    return out


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_samples(path: str) -> List[str]:
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise ValueError(f"No #CHROM header line found in {path}")


_OPS = {"<": operator.lt, ">": operator.gt, "<=": operator.le, ">=": operator.ge}


def read_header_info(path: str) -> dict:
    """Samples plus the INFO / FILTER ids declared in the header (cheap: reads only the header)."""
    samples, infos, filters = [], set(), set()
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID="):
                infos.add(line.split("ID=", 1)[1].split(",", 1)[0])
            elif line.startswith("##FILTER=<ID="):
                filters.add(line.split("ID=", 1)[1].split(",", 1)[0])
            elif line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
            elif not line.startswith("#"):
                break
    return dict(samples=samples, infos=infos, filters=filters)


class VariantFilter:
    """Site-level filter: FILTER labels, hard INFO thresholds, optional INFO/DP quantile bounds."""

    def __init__(self, cfg, infos=None):
        self.enabled = bool(cfg.apply_variant_filter)
        self.hard = [(f, op, float(t)) for f, op, t in cfg.variant_filters]
        self.infos = None if infos is None else set(infos)
        self.dp_lo = self.dp_hi = None
        self.stats = Counter()
        self.dp_values = []
        self._wanted = {f for f, _, _ in self.hard if f != "QUAL"} | {"DP"}
        self.ignore = set(getattr(cfg, "ignore_filter_labels", []) or [])

    def clone(self):
        c = VariantFilter.__new__(VariantFilter)
        c.__dict__.update(self.__dict__)
        c.stats, c.dp_values = Counter(), []
        return c

    def applicable(self):
        return [h for h in self.hard if h[0] == "QUAL" or self.infos is None or h[0] in self.infos]

    def not_applicable(self):
        a = set(self.applicable())
        return [h for h in self.hard if h not in a]

    def set_depth_bounds(self, values, quantiles):
        if not quantiles or len(values) == 0:
            return
        lo, hi = (list(quantiles) + [None, None])[:2]
        v = np.asarray(values, float)
        self.dp_lo = float(np.quantile(v, lo)) if lo is not None else None
        self.dp_hi = float(np.quantile(v, hi)) if hi is not None else None

    @staticmethod
    def label(h):
        return f"{h[0]}{h[1]}{h[2]:g}"

    def _info(self, info):
        out = {}
        for kv in info.split(";"):
            k, _, v = kv.partition("=")
            if k in self._wanted and v:
                try:
                    out[k] = float(v.split(",")[0])
                except ValueError:
                    pass
        return out

    def reasons(self, filt, qual, info):
        r = []
        if filt not in ("PASS", ".", ""):
            r.extend("FILTER=" + x for x in filt.split(";") if x not in self.ignore)
        vals = self._info(info)
        try:
            vals["QUAL"] = float(qual)
        except ValueError:
            pass
        for h in self.hard:
            v = vals.get(h[0])
            if v is not None and _OPS[h[1]](v, h[2]):
                r.append(self.label(h))
        dp = vals.get("DP")
        if dp is not None:
            if self.dp_hi is not None and dp > self.dp_hi:
                r.append("DP>q_hi")
            if self.dp_lo is not None and dp < self.dp_lo:
                r.append("DP<q_lo")
        return r, dp

    def accept(self, filt, qual, info, record_dp=False):
        self.stats["seen"] += 1
        if not self.enabled:
            self.stats["kept"] += 1
            if record_dp:
                dp = self._info(info).get("DP")
                if dp is not None:
                    self.dp_values.append(dp)
            return True
        r, dp = self.reasons(filt, qual, info)
        if record_dp and dp is not None:
            self.dp_values.append(dp)
        for x in r:
            self.stats[x] += 1
        if r:
            self.stats["removed"] += 1
            return False
        self.stats["kept"] += 1
        return True

    def mask_arrays(self, qual, info_arrays, filter_flags):
        """Vectorised version for the scikit-allel backend. info_arrays: name -> float array;
        filter_flags: label -> bool array (True = site carries that FILTER label)."""
        m = len(qual)
        fail = np.zeros(m, bool)
        self.stats["seen"] += m
        if self.enabled:
            for lab, flag in filter_flags.items():
                if lab in self.ignore:
                    continue
                self.stats["FILTER=" + lab] += int(flag.sum())
                fail |= flag
            vals = dict(info_arrays)
            vals["QUAL"] = qual
            for h in self.hard:
                v = vals.get(h[0])
                if v is None:
                    continue
                with np.errstate(invalid="ignore"):
                    bad = _OPS[h[1]](v, h[2]) & np.isfinite(v)
                self.stats[self.label(h)] += int(bad.sum())
                fail |= bad
            dp = vals.get("DP")
            if dp is not None:
                with np.errstate(invalid="ignore"):
                    for bound, op, nm in ((self.dp_hi, operator.gt, "DP>q_hi"), (self.dp_lo, operator.lt, "DP<q_lo")):
                        if bound is not None:
                            bad = op(dp, bound) & np.isfinite(dp)
                            self.stats[nm] += int(bad.sum())
                            fail |= bad
        self.stats["removed"] += int(fail.sum())
        self.stats["kept"] += int((~fail).sum())
        return ~fail

    def summary(self):
        seen = self.stats["seen"]
        return dict(seen=seen, kept=self.stats["kept"], removed=self.stats["removed"],
                    removed_fraction=(self.stats["removed"] / seen) if seen else 0.0,
                    reasons={k: v for k, v in self.stats.items() if k not in ("seen", "kept", "removed")})


def scan_sites(path: str, cfg, vf: VariantFilter, record_dp: bool = True) -> Counter:
    """Read only the first 8 columns of every record (fast): applies `vf`, records INFO/DP and
    counts FILTER labels. Used for the depth quantiles and by `describe`."""
    keep = _contig_filter(cfg)
    labels = Counter()
    with _open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if not keep(f[0]):
                continue
            ref, alt = f[3], f[4]
            if "," in alt or alt in (".", "*") or (cfg.snps_only and (len(ref) != 1 or len(alt) != 1)):
                continue
            labels[f[6]] += 1
            vf.accept(f[6], f[5], f[7], record_dp=record_dp)
    return labels


def _contig_filter(cfg):
    if cfg.contigs:
        allowed = set(cfg.contigs)
        return lambda c: c in allowed
    if cfg.contig_regex:
        rx = re.compile(cfg.contig_regex)
        return lambda c: rx.search(c) is not None
    return lambda c: True


def _mask(dos, dp, gq, cfg):
    dos = dos.copy()
    if dp is not None and cfg.min_dp:
        dos[np.nan_to_num(dp, nan=np.inf) < cfg.min_dp] = -1
    if gq is not None and cfg.min_gq:
        dos[np.nan_to_num(gq, nan=np.inf) < cfg.min_gq] = -1
    return dos


def _iter_builtin(path, sample_idx, cfg, keep_contig, chunk_size, vf=None) -> Iterator[Chunk]:
    chroms, poss, rows, dps, gqs = [], [], [], [], []
    have_dp = False
    cache = {}

    def flush():
        nonlocal chroms, poss, rows, dps, gqs
        if not chroms:
            return None
        dos = np.array(rows, dtype=np.int8)
        dp = np.array(dps, dtype=np.float32) if have_dp else None
        gq = np.array(gqs, dtype=np.float32) if have_dp else None
        if gq is not None and not np.isfinite(gq).any():
            gq = None
        ch = Chunk(np.array(chroms), np.array(poss, dtype=np.int64), _mask(dos, dp, gq, cfg), dp, gq, dos)
        chroms, poss, rows, dps, gqs = [], [], [], [], []
        return ch

    nan = float("nan")
    with _open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if not keep_contig(f[0]):
                continue
            ref, alt = f[3], f[4]
            if "," in alt or alt in (".", "*"):
                continue
            if cfg.snps_only and (len(ref) != 1 or len(alt) != 1):
                continue
            if vf is not None and not vf.accept(f[6], f[5], f[7]):
                continue
            fmt = f[8]
            if fmt not in cache:
                keys = fmt.split(":")
                cache[fmt] = (
                    keys.index("GT") if "GT" in keys else None,
                    keys.index("DP") if "DP" in keys else None,
                    keys.index("GQ") if "GQ" in keys else None,
                )
            gi, di, qi = cache[fmt]
            if gi is None:
                continue
            row, drow, qrow = [], [], []
            for si in sample_idx:
                p = f[9 + si].split(":")
                row.append(_GT.get(p[gi], -1))
                if di is not None:
                    try:
                        drow.append(float(p[di]))
                    except (ValueError, IndexError):
                        drow.append(nan)
                else:
                    drow.append(nan)
                if qi is not None:
                    try:
                        qrow.append(float(p[qi]))
                    except (ValueError, IndexError):
                        qrow.append(nan)
                else:
                    qrow.append(nan)
            if di is not None:
                have_dp = True
            chroms.append(f[0])
            poss.append(int(f[1]))
            rows.append(row)
            dps.append(drow)
            gqs.append(qrow)
            if len(chroms) >= chunk_size:
                yield flush()
    last = flush()
    if last is not None:
        yield last


def _iter_allel(path, sample_idx, all_samples, cfg, keep_contig, chunk_size, vf=None) -> Iterator[Chunk]:
    import allel  # noqa: optional dependency

    wanted = [all_samples[i] for i in sample_idx]
    fields = ["variants/CHROM", "variants/POS", "variants/REF", "variants/ALT",
              "variants/numalt", "calldata/GT", "calldata/DP"]
    try:
        hdr = allel.read_vcf_headers(path)
        has_gq = "GQ" in hdr.formats
        has_dp = "DP" in hdr.formats
    except Exception:  # pragma: no cover
        has_gq, has_dp = True, True
    if not has_dp:
        fields.remove("calldata/DP")
    if has_gq:
        fields.append("calldata/GQ")
    info_fields, filt_ids = [], []
    if vf is not None and vf.enabled:
        try:
            info_fields = [k for k in sorted(vf._wanted) if k in set(hdr.infos)]
            filt_ids = [x for x in hdr.filters if x != "PASS"]
        except Exception:  # pragma: no cover
            pass
        fields += ["variants/QUAL"] + [f"variants/{k}" for k in info_fields] + [f"variants/FILTER_{x}" for x in filt_ids]
    _, file_samples, _, it = allel.iter_vcf_chunks(
        path, fields=fields, samples=wanted, alt_number=1, chunk_length=chunk_size,
        types={"calldata/DP": "i4", "calldata/GQ": "i4"})
    order = np.array([list(file_samples).index(s) for s in wanted])
    for chunk, _, _, _ in it:
        chrom = np.asarray(chunk["variants/CHROM"]).astype(str)
        ok = chunk["variants/numalt"] == 1
        if cfg.snps_only:
            ref = np.char.str_len(np.asarray(chunk["variants/REF"]).astype(str)) == 1
            alt = np.char.str_len(np.asarray(chunk["variants/ALT"]).astype(str)) == 1
            ok &= ref & alt
        ok &= np.array([keep_contig(c) for c in chrom]) if len(chrom) else ok
        if not ok.any():
            continue
        if vf is not None:
            qual = np.asarray(chunk["variants/QUAL"], float)[ok] if "variants/QUAL" in chunk else np.full(ok.sum(), np.nan)
            infos = {k: np.asarray(chunk[f"variants/{k}"], float)[ok] for k in info_fields}
            flags = {x: np.asarray(chunk[f"variants/FILTER_{x}"], bool)[ok] for x in filt_ids}
            keepv = vf.mask_arrays(qual, infos, flags)
            idx_ok = np.flatnonzero(ok)[keepv]
            ok = np.zeros_like(ok)
            ok[idx_ok] = True
            if not ok.any():
                continue
        gt = chunk["calldata/GT"][ok][:, order]
        valid = (gt >= 0).all(axis=-1)
        dos = np.where(valid, (gt > 0).sum(axis=-1), -1).astype(np.int8)
        dp = gq = None
        if "calldata/DP" in chunk:
            dp = chunk["calldata/DP"][ok][:, order].astype(np.float32)
            dp[dp < 0] = np.nan
        if "calldata/GQ" in chunk:
            gq = chunk["calldata/GQ"][ok][:, order].astype(np.float32)
            gq[gq < 0] = np.nan
        yield Chunk(chrom[ok], np.asarray(chunk["variants/POS"])[ok].astype(np.int64),
                    _mask(dos, dp, gq, cfg), dp, gq, dos)


def iter_chunks(path: str, sample_idx: List[int], all_samples: List[str], cfg, vf=None) -> Iterator[Chunk]:
    keep = _contig_filter(cfg)
    backend = cfg.backend
    if backend == "auto":
        try:
            import allel  # noqa: F401
            backend = "allel"
        except ImportError:
            backend = "builtin"
    if backend == "allel":
        return _iter_allel(path, sample_idx, all_samples, cfg, keep, cfg.chunk_size, vf)
    return _iter_builtin(path, sample_idx, cfg, keep, cfg.chunk_size, vf)
