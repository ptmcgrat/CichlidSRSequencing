"""Configuration defaults. Every field can be overridden from a YAML file or the CLI."""
from __future__ import annotations

from dataclasses import dataclass, asdict, field, fields
from typing import Dict, List, Optional

import yaml


DEFAULT_GENOME = "Mzebra_GT3_NCBI"


@dataclass
class Config:
    # ---- windows / runtime
    window_size: int = 250_000
    contigs: Optional[List[str]] = None
    contig_regex: Optional[str] = None
    backend: str = "auto"               # 'auto' | 'builtin' | 'allel'
    chunk_size: int = 20_000
    seed: int = 1
    n_perm: int = 1000
    alpha: float = 0.05

    # ---- genotype / marker QC
    snps_only: bool = True
    min_dp: int = 5
    min_gq: int = 0                     # genotypes below this GQ -> missing (0 = off; the lab pipeline filters on DP only)
    min_call_rate: float = 0.8
    min_mac: int = 2

    # ---- variant-level filtering (bad sites removed before any statistic is computed)
    # A site is removed if its FILTER column holds a label other than PASS/'.', or if any hard filter below
    # fires. Hard filters are skipped when the INFO field is absent from the VCF. Thresholds follow the lab's
    # published recipe (QD<2, FS>40, MQ<50) plus the GATK-recommended SOR/MQRankSum/ReadPosRankSum.
    # DELIBERATELY NOT applied: allele-frequency, InbreedingCoeff/ExcessHet and a lower depth cutoff.
    # They are population-genomic filters that would remove sex-linked sites (XY: all males heterozygous;
    # hemizygous X/Y regions have low depth in one sex) and rare parent-heterozygous markers.
    ignore_filter_labels: List[str] = field(default_factory=list)   # FILTER labels NOT to honour, e.g. [allele_freq, min_DP]
    apply_variant_filter: bool = True
    variant_filters: List[List] = field(default_factory=lambda: [
        ["QUAL", "<", 30.0], ["QD", "<", 2.0], ["FS", ">", 40.0], ["MQ", "<", 50.0],
        ["SOR", ">", 3.0], ["MQRankSum", "<", -12.5], ["ReadPosRankSum", "<", -8.0]])
    # [lower, upper] quantiles of INFO/DP; null = off. Upper tail removes collapsed repeats / paralogs.
    variant_depth_quantiles: Optional[List[Optional[float]]] = field(default_factory=lambda: [None, 0.99])
    max_filtered_fraction_warn: float = 0.6

    # ---- preflight thresholds (sample number)
    min_per_sex_block: int = 3          # fewer of either sex -> BLOCK
    min_per_sex_warn: int = 6           # fewer -> WARN (low power)
    min_per_sex_good: int = 10
    family_min_offspring_warn: int = 12
    min_independent_families: int = 3   # fewer -> WARN (no replication)
    # 'pair': every distinct sire x dam pair is an independent brood (one male crossed to several females is normal);
    # 'parent': broods that share a sire or dam count as ONE family (the strict view)
    family_independence: str = "pair"
    # gsdf-region copy-number genotyping (sexloci gsdf-genotype): depth ratio at 2 copies and the step per extra copy
    gsdf_r2: float = 0.530
    gsdf_step: float = 0.268
    min_stratum_size: int = 4           # strata smaller than this are skipped in association
    min_offspring_per_sex: int = 2      # per family, for the transmission test

    # ---- preflight thresholds (coverage / QC)
    cov_warn: float = 5.0
    cov_bad: float = 3.0
    parent_cov_warn: float = 10.0
    cov_spread_fold: float = 5.0
    sex_depth_ratio_lo: float = 0.6
    sex_depth_ratio_hi: float = 1.67
    sex_depth_p: float = 0.01
    het_depth_rho: float = 0.5
    trio_error_warn: float = 0.05
    trio_min_sites: int = 200
    confound_p: float = 0.01
    qc_outlier_z: float = 6.0
    sample_missing_warn: float = 0.5
    max_call_rate_drop_warn: float = 0.3      # warn when the call-rate filter removed more than this fraction of sites
    # genome version -> ecogroup of the reference species (for a reference-divergence warning)
    genome_ecogroup: Dict[str, str] = field(default_factory=lambda: {
        "Mzebra_UMD2a": "Mbuna", "Mzebra_GT3": "Mbuna", "Mzebra_GT3_NCBI": "Mbuna",
        "Mconophoros_GT1": "Shallow_Benthic"})

    # ---- association / pedigree
    store_min_chi2: float = 10.83
    ped_min_markers: int = 10
    ped_max_markers_per_window: int = 400
    ped_min_valid_frac: float = 0.5     # drop markers with more missing offspring than this
    # The transmission test is deliberately tolerant of low coverage: it reads which parental haplotype each offspring
    # inherited from hundreds of informative markers per window (recombination is rare), and low depth only turns
    # heterozygotes into homozygotes, so single genotypes can be poor. Parents must still be called well.
    ped_min_dp: int = 2                 # offspring genotypes below this depth -> missing (transmission test only)
    ped_min_gq: int = 0                 # offspring GQ filter for the transmission test (0 = off)
    ped_parent_min_dp: int = 6          # parent genotypes below this depth -> missing
    ped_parent_min_gq: int = 20
    ped_cov_warn: float = 2.0           # family offspring below this mean depth: WARN
    ped_cov_bad: float = 1.0            # ... below this: unusable
    # sex contrasts (association, Fst, heterozygosity) need confident genotypes: samples that keep fewer than
    # (1 - this) of their genotypes at min_dp/min_gq are left out of those contrasts (they still count for transmission)
    max_sample_missing_contrast: float = 0.5
    prescan_variants: int = 20000       # sites read first to measure each sample's genotype retention

    # ---- popgen
    min_per_sex_called: int = 3
    diag_tol: float = 0.1
    min_window_markers: int = 5
    min_evidence_markers: int = 30      # a window needs this many markers before Fst / dHet / depth outliers count as evidence
    narrow_spike_bp: int = 20_000       # a region whose significant SNPs span less than this (and are few) is flagged
    narrow_spike_max_markers: int = 10

    # ---- region calling
    z_thr: float = 5.0
    min_fst: float = 0.10
    min_dhet: float = 0.05
    min_diag: int = 3
    min_cov_log2: float = 0.25
    merge_gap: int = 1
    peak_frac: float = 0.9
    min_support_units: int = 2          # independent units that must each be significant for tier A
    block_flag_mb: float = 5.0

    gff: Optional[str] = None
    # mark the known sex-determination loci (gsdf, the two LG10 transposon insertions) on genome_overview.png
    show_known_loci: bool = True
    # YH-style VCFs use the UMD2a chromosome names for GT3 sequence: match the annotation's chromosomes by linkage group
    annotation_alias_by_linkage_group: bool = True
    # display names for plots (the NC_ -> linkage-group table from the lab's Assignment3 script; check it fits your assembly)
    contig_labels: Dict[str, str] = field(default_factory=lambda: {
        f"NC_0{36780 + i}.1": f"LG{n}" for i, n in enumerate(list(range(1, 21)) + [22, 23])})

    # ---- lab data access (read-only; replaces the external FileManager)
    output_root: str = "~/Temp/SexDeterminationAnalysis/"   # default output location
    rclone_remote: str = "ptm_dropbox:/"
    cloud_master_dir: str = "COS/BioSci/BioSci-McGrath/Apps/CichlidSequencingData/"
    local_master_dir: Optional[str] = None                  # default: ~/Temp/CichlidSequencingData/
    cloud_vcf_subdir: str = "Outputs/PedigreeVCFs/"
    sample_sheet_id: str = "1NmgB_TWoO01Qz2ufvECuZFkxXayhUsyu8wQGStVB_8k"      # SampleDatabase + DNAReadsDatabase
    alignment_sheet_id: str = "1vlA_eJ09RSKCaM0foFMc0n5tL95XK4deB3dr7kAFx3E"   # one worksheet per genome version
    sheets_read_only: bool = True                           # request read-only Google scopes

    # WARN codes you have reviewed and accept (shown as INFO, prefixed "[acknowledged]"). BLOCK findings and
    # tier caps are never removed by acknowledging, e.g. [reference_divergence]
    acknowledged_codes: List[str] = field(default_factory=list)


def load_config(path: Optional[str] = None, **overrides) -> Config:
    cfg = Config()
    valid = {f.name for f in fields(Config)}
    data = {}
    if path:
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
    data.update({k: v for k, v in overrides.items() if v is not None})
    unknown = set(data) - valid
    if unknown:
        raise ValueError(f"Unknown config keys: {sorted(unknown)}")
    for k, v in data.items():
        setattr(cfg, k, v)
    return cfg


def config_dict(cfg: Config) -> dict:
    return asdict(cfg)
