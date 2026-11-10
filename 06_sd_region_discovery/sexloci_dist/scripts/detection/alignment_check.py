#!/usr/bin/env python3
"""Which of the chr7-insertion carriers have alignments, and against which genome?

Read-only. For each genome version the lab tracks, this asks the alignment sheet (through sexloci's LabData, the same way sexloci does) whether
each carrier is listed and at what coverage, and (unless --no-cloud) lists the Dropbox folder Bamfiles/<genome>/ with rclone to see whether a
BAM folder exists for it. Every ID we have for a fish is tried: its biosample ID, Munby's sample ID and the provider ID.

run in the environment where the `sexloci` command works:
    conda activate ncbi
    python3 alignment_check.py --carriers chr7ins_carriers_cross_check_Blumer_Malinsky_tables.tsv
writes carrier_alignments.tsv next to where you run it.
"""
import argparse
import re
import subprocess
import sys

import pandas as pd

GENOMES = ["Mzebra_GT3_NCBI", "Mzebra_GT3", "Mzebra_UMD2a", "Mconophoros_GT1"]
norm = lambda x: re.sub(r"[^A-Za-z0-9]", "", str(x)).upper()


def rclone_dirs(path, runner=None):
    """Names of the folders directly under an rclone path, normalised; None if the listing failed."""
    run = runner or (lambda a: subprocess.run(a, capture_output=True, text=True))
    p = run(["rclone", "lsf", "--dirs-only", path])
    if p.returncode != 0:
        return None
    return {norm(x.strip("/")) for x in p.stdout.split() if x.strip("/")}


def check(carriers, tables, bamdirs, genomes=GENOMES):
    """carriers: DataFrame (BioSample, Munby_sample_ID, provider_ID, ...); tables: {genome: DataFrame or None}; bamdirs: {genome: set or None}."""
    index = {}
    for g, a in tables.items():
        if a is None:
            continue
        cov = pd.to_numeric(a["Coverage"], errors="coerce") if "Coverage" in a.columns else pd.Series(float("nan"), index=a.index)
        index[g] = {norm(s): (s, c) for s, c in zip(a["SampleID"], cov)}
    rows = []
    for r in carriers.to_dict("records"):
        ids = [v for v in (r.get("BioSample", ""), r.get("Munby_sample_ID", ""), r.get("provider_ID", "")) if v and str(v) != "nan"]
        out = {k: r.get(k, "") for k in ("BioSample", "Munby_sample_ID", "provider_ID", "location", "sex", "chr7-ins GT", "chr7-ins AC", "in_SampleDatabase")}
        for g in genomes:
            hit = next(((i, index[g][norm(i)]) for i in ids if g in index and norm(i) in index[g]), None)
            out[f"{g}_sheet"] = "no table" if g not in index else (f"yes ({hit[1][1]:.1f}x)" if hit and pd.notna(hit[1][1]) else "yes" if hit else "no")
            out[f"{g}_matched_id"] = hit[0] if hit else ""
            bd = bamdirs.get(g)
            out[f"{g}_cloud_folder"] = "not checked" if bd is None else "yes" if any(norm(i) in bd for i in ids) else "no"
        rows.append(out)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--carriers", default="chr7ins_carriers_cross_check_Blumer_Malinsky_tables.tsv")
    ap.add_argument("--no-cloud", action="store_true", help="skip the rclone listing of Bamfiles/<genome>/")
    ap.add_argument("--out", default="carrier_alignments.tsv")
    a = ap.parse_args()
    from sexloci.config import Config
    from sexloci.labdata import LabData
    lab = LabData.from_config(Config())
    carriers = pd.read_csv(a.carriers, sep="\t", dtype=str).fillna("")
    tables, bamdirs = {}, {}
    for g in GENOMES:
        try:
            t = lab.read_alignment_table(g)
            t = t.loc[:, [c for c in t.columns if not str(c).startswith("Unnamed")]].dropna(subset=["SampleID"])
            t["SampleID"] = t["SampleID"].astype(str).str.strip()
            tables[g] = t
            print(f"alignment sheet {g}: {len(t)} rows")
        except Exception as e:
            tables[g] = None
            print(f"alignment sheet {g}: not read ({str(e)[:90]})")
        if not a.no_cloud:
            bamdirs[g] = rclone_dirs(lab.cloud_master_dir + f"Bamfiles/{g}/")
            print(f"cloud folder Bamfiles/{g}/:", "listing failed" if bamdirs[g] is None else f"{len(bamdirs[g])} sample folders")
    res = check(carriers, tables, bamdirs)
    res.to_csv(a.out, sep="\t", index=False)
    print(f"\n{len(res)} carriers checked; wrote {a.out}\n")
    for g in GENOMES:
        s = res[f"{g}_sheet"].str.startswith("yes").sum()
        c = (res[f"{g}_cloud_folder"] == "yes").sum()
        print(f"{g:16s} in the alignment sheet: {s:3d}   cloud BAM folder: {c:3d}")
    anyal = res[[f"{g}_sheet" for g in GENOMES]].apply(lambda r: any(str(x).startswith("yes") for x in r), axis=1)
    print(f"\ncarriers with an alignment to at least one genome in the sheet: {anyal.sum()} of {len(res)}")


if __name__ == "__main__":
    main()
