"""Simulate a cohort with a known sex locus (XY, ZW or none) plus matching sample/alignment/reads tables.

Pedigree families come from meiosis with crossovers (suppressed inside the sex-linked block in the
heterogametic parent). Unrelated animals are drawn from a population where the heterogametic haplotype
carries fixed differences inside the block. The block may be a single narrow locus or a wide segment;
the pipeline treats both identically. `problems` injects data-quality issues for testing the checks.
"""
from __future__ import annotations

import gzip
import json
import os

import numpy as np
import pandas as pd


def _gamete(rng, h0, h1, pos, chrom_len, no_recomb, sex_pos):
    cuts = np.sort(rng.uniform(0, chrom_len, rng.poisson(1.0)))
    if no_recomb is not None:
        cuts = cuts[(cuts < no_recomb[0]) | (cuts > no_recomb[1])]
    start = rng.integers(2)
    which = (start + np.searchsorted(cuts, pos)) % 2
    origin = (start + np.searchsorted(cuts, sex_pos)) % 2
    return np.where(which == 0, h0, h1), origin


def simulate(outdir, system="XY", mode="both", n_chrom=4, chrom_len=3_000_000, snps_per_chrom=2500,
             sex_chrom=2, sex_pos=1_500_000, block_half_width=250_000, n_families=3, n_off=16,
             n_unrelated=15, depth=18, err=0.01, missing=0.02, seed=1, gz=False,
             species="Mchenga conophoros", genome_version="Mconophoros_GT1", problems=None, prefix=""):
    problems = problems or {}
    rng = np.random.default_rng(seed)
    os.makedirs(outdir, exist_ok=True)
    chroms = [f"chr{i + 1}" for i in range(n_chrom)]
    pos = {c: np.unique(rng.integers(1, chrom_len, snps_per_chrom)) for c in chroms}
    freqs = {c: rng.uniform(0.05, 0.5, len(pos[c])) for c in chroms}
    sc = f"chr{sex_chrom}"
    block = None if system == "none" else (sex_pos - block_half_width, sex_pos + block_half_width)
    diag = {c: rng.random(len(pos[c])) < 0.6 for c in chroms}

    def hap(c, kind):
        h = rng.random(len(freqs[c])) < freqs[c]
        if c == sc and block is not None:
            inb = (pos[c] >= block[0]) & (pos[c] <= block[1]) & diag[c]
            h[inb] = (kind == "het")
        return h

    def parent_haps(het_sex):
        return {c: ((hap(c, "het"), hap(c, "hom")) if (c == sc and het_sex and system != "none")
                    else (hap(c, "hom"), hap(c, "hom"))) for c in chroms}

    samples, sex, parents, brood, geno = [], [], [], [], []

    def add(name, s, par="", br="", g=None):
        samples.append(name); sex.append(s); parents.append(par); brood.append(br); geno.append(g)

    gen = lambda h: {c: h[c][0].astype(int) + h[c][1].astype(int) for c in chroms}
    if mode in ("pedigree", "both"):
        male_het = system == "XY"
        for fi in range(n_families):
            fam = f"F{fi + 1}"
            shared = bool(problems.get("shared_sire")) and fi > 0
            if not shared:
                sh = parent_haps(male_het)
            dh = parent_haps(system == "ZW")
            sire, dam = ("F1_sire_m" if problems.get("shared_sire") else f"{fam}_sire_m"), f"{fam}_dam_f"
            if not shared:
                add(sire, "M", g=gen(sh))
            add(dam, "F", g=gen(dh))
            for j in range(n_off):
                gam, origin = {}, None
                for role, ph in (("sire", sh), ("dam", dh)):
                    g = {}
                    for c in chroms:
                        nr = block if (c == sc and system != "none" and (role == "sire") == male_het) else None
                        g[c], o = _gamete(rng, ph[c][0], ph[c][1], pos[c], chrom_len, nr, sex_pos)
                        if c == sc and nr is not None:
                            origin = o
                    gam[role] = g
                if system == "none":
                    is_male = rng.random() < 0.5
                else:
                    is_male = (origin == 0) if system == "XY" else (origin != 0)
                g = {c: gam["sire"][c].astype(int) + gam["dam"][c].astype(int) for c in chroms}
                add(f"{fam}_o{j + 1}", "M" if is_male else "F", f"{sire},{dam}", f"{fam}-brood", g)
    if mode in ("unrelated", "both") and n_unrelated > 0:
        for lab in ("M", "F"):
            for j in range(n_unrelated):
                het = (system == "XY" and lab == "M") or (system == "ZW" and lab == "F")
                add(f"U_{lab}{j + 1}", lab, g=gen(parent_haps(het)))

    if prefix:
        samples = [prefix + x for x in samples]
        parents = [",".join(prefix + t for t in p.split(",")) if p else "" for p in parents]
    n = len(samples)
    depths = np.full(n, float(depth))
    for i in range(int(problems.get("low_depth_samples", 0))):
        depths[-(i + 1)] = 2.0
    if problems.get("low_depth_offspring"):
        import re as _re
        for i, nm in enumerate(samples):
            if _re.search(r"_o\d+$", nm):
                depths[i] = float(problems["low_depth_offspring"])
    if problems.get("sex_depth_bias"):
        depths[[i for i in range(n) if sex[i] == "M"]] *= 0.5

    path = os.path.join(outdir, "sim.vcf.gz" if gz else "sim.vcf")
    op = gzip.open if gz else open
    gts = ["0/0", "0/1", "1/1"]
    omit = set(problems.get("omit_samples", []))                 # in the sample table but missing from the VCF
    cols = [i for i, nm in enumerate(samples) if nm not in omit]
    n_bad_total = [0]
    with op(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c in chroms:
            fh.write(f"##contig=<ID={c},length={chrom_len}>\n")
        for t in ("GT", "DP", "GQ"):
            fh.write(f'##FORMAT=<ID={t},Number=1,Type=String,Description="{t}">\n')
        if not problems.get("no_info"):
            for k in ("QD", "FS", "MQ", "SOR", "DP"):
                fh.write(f'##INFO=<ID={k},Number=1,Type=Float,Description="{k}">\n')
        if problems.get("annotate_filter"):
            fh.write('##FILTER=<ID=PASS,Description="All filters passed">\n##FILTER=<ID=BadSite,Description="bad">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples[i] for i in cols) + "\n")
        for c in chroms:
            G = np.array([[g[c][k] for g in geno] for k in range(len(pos[c]))])
            G = np.where(rng.random(G.shape) < err, rng.integers(0, 3, G.shape), G)
            bad = rng.random(len(pos[c])) < problems.get("bad_variant_frac", 0.0)
            G[bad] = rng.integers(0, 3, G[bad].shape)               # bad sites carry noise genotypes
            D = rng.poisson(depths[None, :], G.shape)
            # low depth: a heterozygote looks homozygous when every read shows the same allele, p = 0.5**(depth-1)
            under = (G == 1) & (D >= 1) & (rng.random(G.shape) < 0.5 ** (np.maximum(D, 1) - 1))
            G = np.where(under, rng.choice([0, 2], G.shape), G)
            miss = rng.random(G.shape) < missing
            n_bad_total[0] += int(bad.sum())
            for k in range(len(pos[c])):
                if problems.get("no_info"):
                    qual, info, filt = "50", ".", "PASS"
                elif bad[k]:
                    qual, info = "20", f"QD=0.5;FS=80.0;MQ=30.0;SOR=5.0;DP={int(D[k].sum())}"
                    filt = "BadSite" if problems.get("annotate_filter") else "."
                else:
                    qual = "900"
                    info = f"QD={rng.normal(20, 4):.2f};FS={abs(rng.normal(2, 1)):.2f};MQ=60.0;SOR=0.8;DP={int(D[k].sum())}"
                    filt = "PASS" if problems.get("annotate_filter") else "."
                cells = ["./.:0:0" if miss[k, i] else f"{gts[G[k, i]]}:{D[k, i]}:{min(99, 3 * D[k, i])}"
                         for i in cols]
                fh.write(f"{c}\t{pos[c][k]}\t.\tA\tG\t{qual}\t{filt}\t{info}\tGT:DP:GQ\t" + "\t".join(cells) + "\n")

    sx = list(sex)
    for i in range(int(problems.get("unsexed", 0))):
        sx[-(i + 1)] = "U"
    sp = [species] * n
    if problems.get("mixed_species"):
        sp[-1] = "Copadichromis virginalis"
    if problems.get("sex_batch_confound"):
        project = ["PM22" if s == "F" else "PM30" for s in sex]
        for i in rng.choice(n, max(1, n // 12), replace=False):
            project[i] = "PM30" if project[i] == "PM22" else "PM22"
    else:
        project, seen = [], {"M": 0, "F": 0}
        for sx_ in sex:                      # alternate within each sex -> batches are balanced by sex
            project.append("PM22" if seen[sx_] % 2 == 0 else "PM30")
            seen[sx_] += 1
    eco = {"Mchenga conophoros": "Shallow_Benthic"}.get(species, "Shallow_Benthic")
    st = pd.DataFrame({"SampleID": samples, "Sex": sx, "Species": sp, "DoB": "", "BroodID": brood,
                       "Parents": parents, "Ecogroup": eco, "LabReared": 1,
                       "Inversion10": 0, "Category": "ignored"})
    if problems.get("drop_from_table"):
        st = st.iloc[: -int(problems["drop_from_table"])]
    st.to_csv(os.path.join(outdir, "sample_table.tsv"), sep="\t", index=False)
    al = pd.DataFrame({"SampleID": samples, "GenomeVersion": genome_version, "RunIDs": [f"{p}_{s}" for p, s in zip(project, samples)],
                       "Coverage": depths * 0.95, "TotalReads": 1e8, "UnmappedReads": 0.005, "DiscordantReads": 0.03,
                       "InversionReads": 0.01, "DuplicationReads": 0.002, "ClippedReads": 0.03, "ChimericReads": 0.04,
                       "minimap2_version": 2.31, "gatk_version": "4.6.2.0", "pysam_version": "0.24.1", "BamSize": 1e10})
    if problems.get("qc_outlier"):
        al.loc[3, "InversionReads"] = 0.25
    if problems.get("wrong_genome"):
        al.loc[0, "GenomeVersion"] = "Mzebra_GT3"
    al.to_csv(os.path.join(outdir, "alignment_table.tsv"), sep="\t", index=False)
    rd = pd.DataFrame({"SampleID": samples, "ProjectID": project, "RunID": [f"{s}_L001" for s in samples],
                       "Instrument": "NovaSeq", "LibraryID": "lib"})
    rd.to_csv(os.path.join(outdir, "reads_table.tsv"), sep="\t", index=False)
    truth = dict(system=system, chrom=sc, block=block, sex_pos=sex_pos, n_bad_variants=n_bad_total[0])
    with open(os.path.join(outdir, "truth.json"), "w") as fh:
        json.dump(truth, fh, indent=2)
    return dict(vcf=path, sample_table=os.path.join(outdir, "sample_table.tsv"),
                alignment_table=os.path.join(outdir, "alignment_table.tsv"),
                reads_table=os.path.join(outdir, "reads_table.tsv"), truth=truth, genome_version=genome_version)
