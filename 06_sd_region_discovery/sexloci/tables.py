"""Access to the SampleDatabase / DNAReadsDatabase / alignment sheets.

Two interchangeable sources:
  * the lab sheets, read through `labdata.LabData` (the default; read-only)
  * TSV/CSV exports of the sheets (works offline)

The sample/reads tables are never merged: that merge yields one row per sequencing run and would duplicate samples.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd


@dataclass
class Tables:
    sample: pd.DataFrame                    # SampleDatabase
    alignment: Optional[pd.DataFrame]       # alignment sheet for the genome version
    reads: Optional[pd.DataFrame]           # DNAReadsDatabase (batch information)
    genome_version: str
    gff: Optional[str] = None
    source: str = "tsv"


def _clean(df: Optional[pd.DataFrame]) -> Optional[pd.DataFrame]:
    if df is None:
        return None
    df = df.loc[:, [c for c in df.columns if not str(c).startswith("Unnamed")]].copy()
    if "SampleID" not in df.columns:
        raise ValueError("table has no 'SampleID' column")
    df = df.dropna(subset=["SampleID"])
    df["SampleID"] = df["SampleID"].astype(str).str.strip()
    return df[df["SampleID"] != ""].reset_index(drop=True)


def _read_table(path: Optional[str]) -> Optional[pd.DataFrame]:
    if not path:
        return None
    with open(path) as fh:
        first = fh.readline()
    sep = "\t" if "\t" in first else ","
    return pd.read_csv(path, sep=sep, dtype=object)


def load_tables_tsv(sample_table, genome_version, alignment_table=None, reads_table=None) -> Tables:
    sample = _clean(_read_table(sample_table))
    for col in ("Sex", "Species"):
        if col not in sample.columns:
            raise ValueError(f"sample table is missing required column '{col}'")
    align = _clean(_read_table(alignment_table))
    if align is not None and "Coverage" in align:
        align["Coverage"] = pd.to_numeric(align["Coverage"], errors="coerce")
    return Tables(sample, align, _clean(_read_table(reads_table)), genome_version, source="tsv")


def load_tables_labdata(lab, genome_version, log=None) -> Tables:
    """Read SampleDatabase, DNAReadsDatabase and the alignment sheet through a LabData object."""
    sample, reads = lab.read_sample_tables()
    align = _clean(lab.read_alignment_table(genome_version))
    if "Coverage" in align:
        align["Coverage"] = pd.to_numeric(align["Coverage"], errors="coerce")
    return Tables(_clean(sample), align, _clean(reads), genome_version, lab.ensure_gff(genome_version), "labdata")
