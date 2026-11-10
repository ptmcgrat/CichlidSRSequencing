"""Command line interface:  sexloci run | check | run-folder | list-vcfs | describe | filter-audit | filter-compare | gsdf-genotype | sv-genotype | write-builtin-vcf | gsdf-summary | simulate

Output goes to ~/Temp/SexDeterminationAnalysis/ by default (config key output_root)."""
from __future__ import annotations

import argparse
import os
import sys

from .config import DEFAULT_GENOME, load_config


def _cfg_args(p):
    p.add_argument("--config", help="YAML overriding defaults (config/example_config.yaml)")
    p.add_argument("--window-size", type=int)
    p.add_argument("--n-perm", type=int)
    p.add_argument("--contigs", nargs="+")
    p.add_argument("--contig-regex")
    p.add_argument("--backend", choices=["auto", "builtin", "allel"])
    p.add_argument("--gff")
    p.add_argument("--seed", type=int)
    p.add_argument("--family-independence", choices=["pair", "parent"],
                   help="pair (default): each sire x dam pair is an independent brood; parent: broods sharing a parent count as one")


def _lab_args(p):
    p.add_argument("--genome-version", default=DEFAULT_GENOME, help="selects the alignment sheet (default: %(default)s)")
    p.add_argument("--local-master-dir", help="local mirror of the cloud data (default ~/Temp/CichlidSequencingData/)")


def _add_common(p):
    p.add_argument("--vcf", required=True, help="multi-sample VCF(.gz) containing ONE species")
    _lab_args(p)
    off = p.add_argument_group("offline alternative to the lab sheets (TSV/CSV exports)")
    off.add_argument("--sample-table", help="SampleDatabase export (needs SampleID, Sex, Species)")
    off.add_argument("--alignment-table", help="alignment-sheet export")
    off.add_argument("--reads-table", help="DNAReadsDatabase export (used for batch checks)")
    p.add_argument("--outdir", help="default: ~/Temp/SexDeterminationAnalysis/<vcf name>/")
    _cfg_args(p)


def _config(a):
    return load_config(getattr(a, "config", None), window_size=getattr(a, "window_size", None),
                       n_perm=getattr(a, "n_perm", None), contigs=getattr(a, "contigs", None),
                       contig_regex=getattr(a, "contig_regex", None), backend=getattr(a, "backend", None),
                       gff=getattr(a, "gff", None), seed=getattr(a, "seed", None),
                       local_master_dir=getattr(a, "local_master_dir", None),
                       family_independence=getattr(a, "family_independence", None))


def output_root(cfg) -> str:
    return os.path.expanduser(cfg.output_root).rstrip("/") + "/"


def default_outdir(cfg, vcf) -> str:
    base = os.path.basename(vcf)
    stem = base[:-7] if base.endswith(".vcf.gz") else base[:-4] if base.endswith(".vcf") else base
    return output_root(cfg) + stem


def _load(a):
    from .labdata import LabData
    from .tables import load_tables_labdata, load_tables_tsv
    cfg = _config(a)
    if a.sample_table:
        t = load_tables_tsv(a.sample_table, a.genome_version, a.alignment_table, a.reads_table)
    else:
        t = load_tables_labdata(LabData.from_config(cfg), a.genome_version)
    return cfg, t


def _print_findings(F):
    for f in F.sorted():
        print(f"[{f.severity}] {f.code}: {f.message}")
        if f.why:
            print(f"    why: {f.why}")
        if f.fix:
            print(f"    fix: {f.fix}")


def build_parser():
    ap = argparse.ArgumentParser(prog="sexloci", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    _add_common(sub.add_parser("run", help="preflight checks, then map sex-linked loci"))
    _add_common(sub.add_parser("check", help="preflight checks only (no genotype pass)"))
    # ---- folder workflow (download through the built-in LabData, analyse each VCF)
    rf = sub.add_parser("run-folder", help="download every VCF in a cloud folder and analyse each one")
    lv = sub.add_parser("list-vcfs", help="list the VCFs in the cloud folder")
    for p_ in (rf, lv):
        _lab_args(p_)
        p_.add_argument("--cloud-subdir", help="folder under the cloud master dir (default Outputs/PedigreeVCFs/)")
        p_.add_argument("--local-dir", help="where to download (default: <local master dir>/<cloud-subdir>)")
        p_.add_argument("--include", help="glob on file names, e.g. 'YH*'")
    lv.add_argument("--config", help="YAML overriding defaults")
    _cfg_args(rf)
    rf.add_argument("--outdir", help="default: ~/Temp/SexDeterminationAnalysis/ (one sub-folder per VCF)")
    rf.add_argument("--manifest", help="TSV: vcf, genome_version, skip, note (per-file overrides)")
    rf.add_argument("--check-only", action="store_true", help="download + preflight only")
    rf.add_argument("--refresh", action="store_true", help="re-download even if present locally")
    gg = sub.add_parser("gsdf-genotype", help="genotype the gsdf-region duplication for every sample with a GT3a alignment (reads only the region)")
    _lab_args(gg)
    gg.add_argument("--config", help="YAML overriding defaults (gsdf_r2 / gsdf_step set the calibration)")
    gg.add_argument("--outdir", help="default: ~/Temp/SexDeterminationAnalysis/gsdf_genotypes/")
    gg.add_argument("--jobs", type=int, default=4, help="samples fetched in parallel (default 4; Dropbox may throttle more)")
    gg.add_argument("--species", help="only samples whose Species contains this text")
    gg.add_argument("--ecogroup", help="only this Ecogroup, e.g. Mbuna")
    gg.add_argument("--samples", nargs="+", help="only these SampleIDs (globs allowed, e.g. 'MC_*')")
    gg.add_argument("--limit", type=int, help="process at most this many samples (for a trial run)")
    gg.add_argument("--min-coverage", type=float, default=2.0, help="skip samples with genome-wide coverage below this (default 2)")
    gg.add_argument("--refresh", action="store_true", help="measure everything again (default: resume, reusing finished samples)")
    gg.add_argument("--no-keep-bams", action="store_true", help="delete each region BAM after measuring it")
    gg.add_argument("--anchors", help="TSV (SampleID, copies) of samples of known copy number, to calibrate on instead of the MC family")
    gg.add_argument("--theory", action="store_true", help="no calibration: ratio = copies / 4")
    gg.add_argument("--extra-columns", nargs="+", help="more sample-sheet columns to carry into the table (Subgroup is added automatically)")
    wr = gg.add_argument_group("writing the result into the Sample Database sheet (opt-in; dry run unless --yes)")
    wr.add_argument("--write-sample-database", action="store_true",
                    help="add a column 'gsdf number' (total gsdf copies) to the SampleDatabase worksheet. Shows a plan only unless --yes")
    wr.add_argument("--yes", action="store_true", help="really write (without it --write-sample-database is a dry run)")
    wr.add_argument("--overwrite", action="store_true", help="replace values already in the column (default: fill blanks only)")
    wr.add_argument("--include-flagged", action="store_true", help="also write flagged calls (low depth, uncertain, out of range); default: leave them blank")
    wr.add_argument("--column-name", default="gsdf number", help="name of the column to create or fill (default: 'gsdf number')")
    gg.add_argument("--port", type=int, help="port for the local rclone web server (default: a free one)")
    off = gg.add_argument_group("offline alternative to the lab sheets (TSV/CSV exports)")
    off.add_argument("--sample-table", help="SampleDatabase export (SampleID, Species, Ecogroup, Sex)")
    off.add_argument("--alignment-table", help="alignment-sheet export (SampleID, Coverage)")
    sv = sub.add_parser("sv-genotype", help="genotype the gsdf duplication AND long insertions (transposons) for every sample, streaming only the needed BAM regions")
    _lab_args(sv)
    sv.add_argument("--config", help="YAML overriding defaults (gsdf_r2 / gsdf_step set the gsdf calibration)")
    sv.add_argument("--outdir", help="default: ~/Temp/SexDeterminationAnalysis/sv_genotypes/")
    sv.add_argument("--indel-vcf", help="VCF-style file of anchored insertions/deletions (CHROM POS ID REF ALT ...). Default: the insertions built into sexloci "
                                        "(--builtin-sites)")
    sv.add_argument("--builtin-sites", nargs="+", default=["all"], metavar="SET",
                    help="built-in site set(s) used when --indel-vcf is not given: deep_benthic_TE (TE1, TE2 on NC_135176.1), chr7_gsdf_ins (chr7ins_c1, chr7ins_c2: the "
                         "transposon insertion upstream of gsdf, at both GT3a copies; also written as one combined `chr7ins` call), or all (default)")
    sv.add_argument("--no-indel", action="store_true", help="skip the insertions (gsdf copy number only)")
    sv.add_argument("--indel-names", nargs="+", help="names for the sites, in order (also the Sample Database column names). Default: the ID of each record "
                                                    "(TE1 TE2 chr7ins_c1 chr7ins_c2 for the built-in sets; for a VCF the ID column, else TE_<pos>). Renaming "
                                                    "the members of a combined call (chr7ins_c1, chr7ins_c2) switches the combined call off")
    sv.add_argument("--reference-fasta", help="the genome FASTA the VCF was made against (default: the only FASTA in the genome folder)")
    sv.add_argument("--flanking", type=int, default=250, help="reference bases kept on each side of an insertion in the local alleles (default 250)")
    sv.add_argument("--aligner", choices=["auto", "parasail", "numpy"], default="auto",
                    help="read aligner: parasail (fast; pip install parasail) or numpy (slower, no extra install). Default: parasail if installed")
    sv.add_argument("--min-informative-pairs", type=int, default=5, help="read pairs that cross a junction needed for a call (default 5)")
    sv.add_argument("--min-gq", type=float, default=10.0, help="phred-scaled genotype quality needed to write a call into the sheet (default 10)")
    sv.add_argument("--no-gsdf", action="store_true", help="skip the gsdf copy number (insertions only)")
    sv.add_argument("--gsdf-cache", help="finished gsdf measurements to reuse (default: the gsdf-genotype output folder, if present)")
    sv.add_argument("--jobs", type=int, default=4, help="samples fetched in parallel (default 4; Dropbox may throttle more)")
    sv.add_argument("--species", help="only samples whose Species contains this text")
    sv.add_argument("--ecogroup", help="only this Ecogroup, e.g. Deep_Benthic")
    sv.add_argument("--samples", nargs="+", help="only these SampleIDs (globs allowed, e.g. 'AU_*')")
    sv.add_argument("--limit", type=int, help="process at most this many samples (for a trial run)")
    sv.add_argument("--min-coverage", type=float, default=2.0, help="skip samples with genome-wide coverage below this (default 2)")
    sv.add_argument("--refresh", action="store_true",
                    help="measure the SELECTED samples again, e.g. after adding sequencing to them (use with --samples/--species); everyone else's results are kept. "
                         "Default: resume, per locus")
    sv.add_argument("--no-keep-bams", action="store_true", help="delete each region BAM after use")
    sv.add_argument("--anchors", help="TSV (SampleID, copies) of samples of known gsdf copy number, to calibrate on instead of the MC family")
    sv.add_argument("--theory", action="store_true", help="gsdf: no calibration, ratio = copies / 4")
    sv.add_argument("--extra-columns", nargs="+", help="more sample-sheet columns to carry into the table (Subgroup is added automatically)")
    wv = sv.add_argument_group("writing the results into the Sample Database sheet (opt-in; dry run unless --yes)")
    wv.add_argument("--write-sample-database", action="store_true",
                    help="add/fill one column per locus: 'gsdf number' and one genotype column per insertion site. Shows a plan only unless --yes")
    wv.add_argument("--yes", action="store_true", help="really write (without it --write-sample-database is a dry run)")
    wv.add_argument("--overwrite", action="store_true", help="replace values already in the columns (default: fill blanks only)")
    wv.add_argument("--include-flagged", action="store_true", help="gsdf: also write flagged calls (low depth, uncertain)")
    wv.add_argument("--include-lowqual", action="store_true", help="insertions: also write calls below --min-gq (never writes no-calls)")
    wv.add_argument("--gsdf-column-name", default="gsdf number", help="name of the gsdf column (default: 'gsdf number')")
    sv.add_argument("--port", type=int, help="port for the local rclone web server (default: a free one)")
    offv = sv.add_argument_group("offline alternative to the lab sheets (TSV/CSV exports)")
    offv.add_argument("--sample-table", help="SampleDatabase export (SampleID, Species, Ecogroup, Sex)")
    offv.add_argument("--alignment-table", help="alignment-sheet export (SampleID, Coverage)")
    bv = sub.add_parser("write-builtin-vcf", help="write the built-in insertion sites (TE1/TE2 and the chr7 insertion upstream of gsdf) as a VCF file")
    bv.add_argument("out", help="output path, e.g. transposons.vcf")
    bv.add_argument("--set", default="deep_benthic_TE", help="built-in site set: deep_benthic_TE (default), chr7_gsdf_ins, or all")
    gs = sub.add_parser("gsdf-summary", help="summarise gsdf_genotypes.tsv by any column (Ecogroup, Subgroup, Species, ...)")
    gs.add_argument("--table", help="default: ~/Temp/SexDeterminationAnalysis/gsdf_genotypes/gsdf_genotypes.tsv")
    gs.add_argument("--by", nargs="+", default=["Ecogroup"], help="column(s) to group by (default Ecogroup)")
    gs.add_argument("--min-n", type=int, default=1, help="hide groups with fewer genotyped samples")
    gs.add_argument("--include-pedigree", action="store_true",
                    help="also count fish that belong to a pedigree (default: independent fish only, because relatives share haplotypes)")
    ds = sub.add_parser("describe", help="inspect a local VCF: FILTER annotations, INFO fields, what the filter removes")
    ds.add_argument("--vcf", nargs="+", required=True)
    ds.add_argument("--config")
    fa = sub.add_parser("filter-audit", help="where along the genome does the site filter remove sites? (table + whole-genome plot)")
    fa.add_argument("--vcf", required=True)
    fa.add_argument("--outdir", help="default: ~/Temp/SexDeterminationAnalysis/filter_audit/<vcf name>/")
    fa.add_argument("--config")
    fa.add_argument("--window-size", type=int)
    fa.add_argument("--contig-regex")
    fc = sub.add_parser("filter-compare", help="compare a RAW VCF with an already-filtered VCF of the same data")
    fc.add_argument("--raw", required=True)
    fc.add_argument("--filtered", required=True)
    fc.add_argument("--outdir", help="default: ~/Temp/SexDeterminationAnalysis/filter_compare/")
    fc.add_argument("--config")
    s = sub.add_parser("simulate", help="write a simulated VCF + sample/alignment/reads tables")
    s.add_argument("--outdir", help="default: ~/Temp/SexDeterminationAnalysis/simulated/")
    s.add_argument("--system", choices=["XY", "ZW", "none"], default="XY")
    s.add_argument("--mode", choices=["pedigree", "unrelated", "both"], default="both")
    s.add_argument("--seed", type=int, default=1)
    return ap


def main(argv=None):
    a = build_parser().parse_args(argv)

    if a.cmd == "simulate":
        from .simulate import simulate
        a.outdir = a.outdir or output_root(load_config()) + "simulated"
        print(simulate(a.outdir, system=a.system, mode=a.mode, seed=a.seed))
        return 0
    if a.cmd == "describe":
        import json
        from .describe import describe_vcf
        cfg = load_config(a.config)
        for v in a.vcf:
            d = describe_vcf(v, cfg)
            print(f"\n== {v}\n{json.dumps(d, indent=2, default=str)}")
        return 0
    if a.cmd == "filter-audit":
        from .describe import audit_filter_by_window
        from .plots import filter_audit_plot
        cfg = load_config(a.config, window_size=a.window_size, contig_regex=a.contig_regex)
        out = a.outdir or output_root(cfg) + "filter_audit/" + os.path.basename(default_outdir(cfg, a.vcf))
        os.makedirs(out, exist_ok=True)
        df = audit_filter_by_window(a.vcf, cfg)
        df.to_csv(os.path.join(out, "filter_audit_by_window.tsv"), sep="\t", index=False)
        filter_audit_plot(df, cfg, os.path.join(out, "filter_audit_genome.png"))
        tot, rem = int(df.sites.sum()), int(df.removed.sum())
        print(f"{tot:,} SNPs, {rem:,} flagged ({rem / tot:.1%}); results in {out}")
        low = df[df.sites >= 20].nsmallest(8, "retained")[["chrom", "start", "sites", "retained"]]
        print("windows keeping the fewest sites:\n" + low.round(3).to_string(index=False))
        return 0
    if a.cmd == "filter-compare":
        import json
        from .describe import compare_filters
        cfg = load_config(a.config)
        os.makedirs(a.outdir, exist_ok=True)
        a.outdir = a.outdir or output_root(cfg) + "filter_compare"
        os.makedirs(a.outdir, exist_ok=True)
        r = compare_filters(a.raw, a.filtered, cfg)
        from .plots import retention_plot
        retention_plot(r["by_window"], cfg, os.path.join(a.outdir, "retention_genome.png"))
        r["by_contig"].to_csv(os.path.join(a.outdir, "retention_by_contig.tsv"), sep="\t", index=False)
        r["by_window"].to_csv(os.path.join(a.outdir, "retention_by_window.tsv"), sep="\t", index=False)
        with open(os.path.join(a.outdir, "filter_compare_summary.json"), "w") as fh:
            json.dump(r["summary"], fh, indent=2)
        print(json.dumps(r["summary"], indent=2))
        return 0
    if a.cmd == "gsdf-genotype":
        from .gsdf import run_gsdf_genotype
        from .labdata import LabData
        from .tables import load_tables_labdata, load_tables_tsv
        cfg = _config(a)
        lab = LabData.from_config(cfg, log=lambda m: print(m, flush=True))
        t = (load_tables_tsv(a.sample_table, a.genome_version, a.alignment_table) if a.sample_table
             else load_tables_labdata(lab, a.genome_version))
        out = a.outdir or output_root(cfg) + "gsdf_genotypes"
        if a.write_sample_database and a.sample_table:
            sys.exit("--write-sample-database needs the live Sample Database, not an exported --sample-table")
        full = run_gsdf_genotype(t, cfg, out, lab, a.genome_version, a.jobs, a.species, a.ecogroup, a.samples, a.limit, a.min_coverage,
                                 a.refresh, not a.no_keep_bams, a.anchors, a.theory, a.port, log=lambda m: print(m, flush=True),
                                 extra_columns=a.extra_columns)
        if a.write_sample_database:
            import dataclasses
            from . import __version__
            from .samplesheet import update_sample_database
            lab_w = LabData.from_config(dataclasses.replace(cfg, sheets_read_only=False), log=lambda m: print(m, flush=True))
            update_sample_database(lab_w.open_sample_worksheet(), full, out, a.column_name, a.include_flagged, a.overwrite, a.yes,
                                   __version__, log=lambda m: print(m, flush=True))
        return 0
    if a.cmd == "write-builtin-vcf":
        from .known_sites import fingerprint, write_builtin_vcf, BUILTIN
        write_builtin_vcf(a.out, a.set)
        for c_, p_, i_, r_, al_, _ in (BUILTIN[a.set] if a.set != "all" else [x for v in BUILTIN.values() for x in v]):
            print(f"{i_}: {c_}:{p_}  {fingerprint(al_)}")
        print(f"wrote {a.out}")
        return 0
    if a.cmd == "sv-genotype":
        from .labdata import LabData
        from .svgeno import (build_classifiers, find_reference_fasta, load_builtin_sites, load_indel_sites, load_reference_windows, run_sv_genotype,
                             sheet_specs, site_groups)
        from .tables import load_tables_labdata, load_tables_tsv
        cfg = _config(a)
        log = lambda m: print(m, flush=True)
        if a.write_sample_database and a.sample_table:
            sys.exit("--write-sample-database needs the live Sample Database, not an exported --sample-table")
        if a.no_gsdf and a.no_indel:
            sys.exit("nothing to genotype: --no-gsdf and --no-indel together")
        lab = LabData.from_config(cfg, log=log)
        t = (load_tables_tsv(a.sample_table, a.genome_version, a.alignment_table) if a.sample_table
             else load_tables_labdata(lab, a.genome_version))
        out = a.outdir or output_root(cfg) + "sv_genotypes"
        sites, clfs = [], []
        if not a.no_indel:
            from .known_sites import fingerprint
            if a.indel_vcf:
                sites = load_indel_sites(a.indel_vcf, a.indel_names)
                src = a.indel_vcf
            else:
                sites = load_builtin_sites(a.builtin_sites, a.indel_names)
                src = "built-in set(s) " + ", ".join(f"'{x}'" for x in a.builtin_sites)
            fasta = find_reference_fasta(lab.genome_dir(a.genome_version), a.reference_fasta)
            wins = {}
            for s_ in sites:
                lo, hi = max(0, s_.pos - 1 - a.flanking - 10), s_.pos - 1 + len(s_.ref) + a.flanking + 10
                old = wins.get(s_.chrom, (lo, hi))
                wins[s_.chrom] = (min(old[0], lo), max(old[1], hi))
            log(f"[sv] reading the reference windows from {fasta}")
            clfs = build_classifiers(sites, load_reference_windows(fasta, wins), a.flanking, a.aligner)
            be = clfs[0].backend.name
            log(f"[sv] {len(sites)} insertion sites from {src}; aligner: {be}")
            for s_ in sites:
                log(f"[sv]   {s_.name}: {s_.chrom}:{s_.pos}, +{s_.length:,} bp inserted; {fingerprint(s_.alt)}")
            for g_, m_ in site_groups(sites).items():
                log(f"[sv] combined call {g_} = {' + '.join(m_)} (informative pairs summed)")
            if be == "numpy":
                log("[sv] note: parasail is not installed, so the slower numpy aligner is used (pip install parasail for the fast one)")
        gcache = a.gsdf_cache or output_root(cfg) + "gsdf_genotypes/gsdf_measurements.tsv"
        full = run_sv_genotype(t, cfg, out, lab, a.genome_version, sites, clfs, not a.no_gsdf, a.jobs, a.species, a.ecogroup, a.samples, a.limit,
                               a.min_coverage, a.refresh, not a.no_keep_bams, a.anchors, a.theory, a.port, a.min_informative_pairs, a.min_gq,
                               gcache, log=log, extra_columns=a.extra_columns)
        if a.write_sample_database:
            import dataclasses
            from . import __version__
            from .samplesheet import update_sample_columns
            lab_w = LabData.from_config(dataclasses.replace(cfg, sheets_read_only=False), log=log)
            specs = sheet_specs(full, sites, a.gsdf_column_name, not a.no_gsdf, a.include_flagged, a.min_gq, a.include_lowqual, __version__)
            update_sample_columns(lab_w.open_sample_worksheet(), specs, out, a.overwrite, a.yes, __version__, log=log)
        return 0
    if a.cmd == "gsdf-summary":
        import pandas as pd
        from .gsdf import summarize_by
        cfg = _config(a)
        path = a.table or output_root(cfg) + "gsdf_genotypes/gsdf_genotypes.tsv"
        tab = pd.read_csv(path, sep="\t", dtype={"Genotype": str})
        if not a.include_pedigree:
            if "Pedigree" in tab.columns:
                n_ped = int(((tab.Pedigree.fillna("") != "") & (tab.Genotype.fillna("").astype(str) != "")).sum())
                print(f"[gsdf] independent fish only: {n_ped} genotyped fish from pedigrees left out (--include-pedigree keeps them)")
            else:
                print("[gsdf] WARNING: this table has no Pedigree column, so related fish are still counted; rerun gsdf-genotype "
                      "(cached, no downloads) to add it")
        for col in a.by:
            res = summarize_by(tab, col, a.min_n, independent=not a.include_pedigree)
            dest = os.path.join(os.path.dirname(path), f"gsdf_summary_by_{col}.tsv")
            res.round(4).to_csv(dest, sep="\t", index=False)
            with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_rows", 500):
                print(f"\n== by {col} ({len(res)} groups) -> {dest}")
                print(res.drop(columns=["unknown_sex", "unknown_sex_carriers"]).round({"male_vs_female_p": 4}).to_string(index=False))
        return 0
    if a.cmd in ("run-folder", "list-vcfs"):
        from .folder import FolderRunner, read_manifest
        from .labdata import LabData
        cfg = _config(a)
        lab = LabData.from_config(cfg, log=lambda m: print(m, flush=True))
        fr = FolderRunner(lab, a.genome_version, cfg, a.cloud_subdir, a.local_dir,
                          log=lambda m: print("[sexloci]", m, flush=True))
        if a.cmd == "list-vcfs":
            print("\n".join(fr.list_vcfs(a.include)))
            return 0
        outdir = a.outdir or output_root(cfg)
        df = fr.run(outdir, read_manifest(a.manifest), a.include, a.check_only, a.refresh)
        cols = ["vcf", "status", "n_block", "n_warn", "block_codes", "n_regions", "top_tier", "top_region", "model"]
        print("\n" + df[cols].to_string(index=False))
        for name, F in fr.last_findings.items():
            shown = [f for f in F.sorted() if f.severity in ("BLOCK", "WARN")]
            if shown:
                print(f"\n== {name}")
                for f in shown:
                    print(f"  [{f.severity}] {f.code}: {f.message[:700]}")
                    if f.why:
                        print(f"      why: {f.why}")
                    if f.fix:
                        print(f"      fix: {f.fix}")
            if any(f.rows for f in F):
                print(f"  -> per-sample list: {outdir.rstrip('/')}/{name[:-7] if name.endswith('.vcf.gz') else name[:-4]}/fix_list.tsv")
        print(f"\nresults in {outdir}   (all findings: {outdir.rstrip('/')}/folder_findings.tsv)")
        return 0 if (df.status.isin(["OK", "CHECK_OK", "SKIPPED"])).all() else 2
    from .pipeline import PreflightBlocked, preflight, run_pipeline
    cfg, tables = _load(a)
    a.outdir = a.outdir or default_outdir(cfg, a.vcf)
    if a.cmd == "check":
        _, F, _ = preflight(a.vcf, tables, cfg, a.outdir)
        _print_findings(F)
        return 2 if F.blocked else 0
    try:
        out = run_pipeline(a.vcf, tables, a.outdir, cfg)
    except PreflightBlocked as e:
        print("\nRUN BLOCKED - fix these and rerun:")
        _print_findings(e.findings)
        print(f"\n(details: {a.outdir}/preflight_report.md)")
        return 2
    _print_findings(out["findings_obj"])
    reg = out["regions"]
    print(f"\n{len(reg)} candidate region(s); see {a.outdir}/report.md")
    if len(reg):
        print(reg[["chrom", "start", "end", "tier", "model", "units_supporting", "units_total"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
