"""Download every VCF in a cloud folder (through the built-in read-only LabData / rclone) and analyse each
one on its own. Every VCF is a separate analysis unit; a blocked or failing VCF is recorded and the batch continues.
"""
from __future__ import annotations

import datetime
import fnmatch
import os
import traceback
from dataclasses import dataclass
from typing import Dict, List, Optional

import pandas as pd

from .config import Config
from .pipeline import PreflightBlocked, preflight, run_pipeline
from .report import findings_frame
from .tables import load_tables_labdata

@dataclass
class Entry:
    vcf: str
    genome_version: Optional[str] = None
    skip: bool = False
    note: str = ""


def read_manifest(path: Optional[str]) -> Dict[str, Entry]:
    """Optional TSV/CSV with columns vcf, genome_version, skip, note (per-file overrides)."""
    if not path:
        return {}
    with open(path) as fh:
        sep = "\t" if "\t" in fh.readline() else ","
    df = pd.read_csv(path, sep=sep, dtype=str).fillna("")
    out = {}
    for _, r in df.iterrows():
        out[r["vcf"].strip()] = Entry(r["vcf"].strip(), r.get("genome_version", "").strip() or None,
                                      str(r.get("skip", "")).strip().lower() in ("1", "true", "yes", "y"),
                                      r.get("note", ""))
    return out


class FolderRunner:
    def __init__(self, lab, default_genome: str, cfg: Config, cloud_subdir: Optional[str] = None,
                 local_dir: Optional[str] = None, log=print):
        self.lab, self.default_genome, self.cfg, self.log = lab, default_genome, cfg, log
        sub = (cloud_subdir or cfg.cloud_vcf_subdir).strip("/") + "/"
        self.local_dir = (local_dir.rstrip("/") + "/") if local_dir else lab.local_master_dir + sub
        self._tables: Dict[str, object] = {}

    def tables(self, genome):
        if genome not in self._tables:
            self.log(f"reading sample/alignment databases for {genome} ...")
            self._tables[genome] = load_tables_labdata(self.lab, genome)
        return self._tables[genome]

    def list_vcfs(self, include: Optional[str] = None) -> List[str]:
        names = self.lab.list_files(self.local_dir)
        vcfs = sorted(n for n in names if n.endswith((".vcf.gz", ".vcf")))
        return [n for n in vcfs if not include or fnmatch.fnmatch(n, include)]

    def fetch(self, name: str, refresh: bool = False) -> str:
        """Bring the VCF (and its index) up to date with the cloud copy; returns the local path. The action taken for
        the VCF is kept in `self.last_action`."""
        for fname in (name, name + ".tbi"):
            local = self.local_dir + fname
            try:
                action = self.lab.sync(local, fast=fname.endswith(".gz"), force=refresh)
            except FileNotFoundError:
                if fname.endswith(".tbi"):
                    continue                    # the index is optional for us
                raise
            if fname == name:
                self.last_action = action
            self.log(f"{fname}: {action}" + (" (the cloud copy differs from the local one)" if action == "updated" else ""))
        return self.local_dir + name

    def run(self, outdir, manifest=None, include=None, check_only=False, refresh=False, quiet=True) -> pd.DataFrame:
        manifest = manifest or {}
        rows = []
        self.last_findings = {}
        for name in self.list_vcfs(include):
            ent = manifest.get(name, Entry(name))
            stem = name[:-7] if name.endswith(".vcf.gz") else name[:-4]
            row = dict(vcf=name, genome_version=ent.genome_version or self.default_genome, status="", n_block=0,
                       download="", vcf_modified="", n_warn=0, block_codes="", tier_cap="", n_regions=0, top_tier="", top_region="", model="",
                       variants_kept=None, removed_fraction=None, message="")
            if ent.skip:
                row["status"] = "SKIPPED"
                rows.append(row)
                continue
            out_i = os.path.join(outdir, stem)
            try:
                path = self.fetch(name, refresh)
                row["download"] = self.last_action
                row["vcf_modified"] = datetime.datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M")
                tables = self.tables(row["genome_version"])
                if check_only:
                    unit, F, _ = preflight(path, tables, self.cfg, out_i)
                    row["status"] = "BLOCKED" if F.blocked else "CHECK_OK"
                    self._fill(row, F)
                    self.last_findings[name] = F
                    row["message"] = row["block_codes"]
                else:
                    out = run_pipeline(path, tables, out_i, self.cfg, quiet=quiet)
                    row["status"] = "OK"
                    self._fill(row, out["findings_obj"])
                    self.last_findings[name] = out["findings_obj"]
                    reg = out["regions"]
                    row["n_regions"] = len(reg)
                    vs = out.get("variant_filter", {})
                    row["variants_kept"], row["removed_fraction"] = out["n_variants_used"], vs.get("removed_fraction")
                    if len(reg):
                        r = reg.iloc[0]
                        row.update(top_tier=r.tier, top_region=f"{r.chrom}:{r.start}-{r.end}", model=r.model)
            except PreflightBlocked as e:
                row["status"] = "BLOCKED"
                self._fill(row, e.findings)
                self.last_findings[name] = e.findings
                row["message"] = row["block_codes"]
            except Exception as e:                          # keep the batch going
                row["status"], row["message"] = "ERROR", f"{type(e).__name__}: {e}"[:300]
                os.makedirs(out_i, exist_ok=True)
                with open(os.path.join(out_i, "error.txt"), "w") as fh:
                    fh.write(traceback.format_exc())
            self.log(f"{name}: {row['status']} {row['message']}")
            rows.append(row)
        df = pd.DataFrame(rows)
        os.makedirs(outdir, exist_ok=True)
        df.to_csv(os.path.join(outdir, "folder_summary.tsv"), sep="\t", index=False)
        frames = [findings_frame(F).assign(vcf=n) for n, F in self.last_findings.items() if len(F)]
        if frames:
            pd.concat(frames)[["vcf", "severity", "code", "message", "why", "fix", "n_samples", "samples", "tier_cap"]].to_csv(
                os.path.join(outdir, "folder_findings.tsv"), sep="\t", index=False)
        return df

    @staticmethod
    def _fill(row, F):
        row["n_block"] = sum(f.severity == "BLOCK" for f in F)
        row["block_codes"] = ";".join(dict.fromkeys(f.code for f in F if f.severity == "BLOCK"))
        row["n_warn"] = sum(f.severity == "WARN" for f in F)
        row["tier_cap"] = F.tier_cap() or "" if hasattr(F, "tier_cap") else ""
