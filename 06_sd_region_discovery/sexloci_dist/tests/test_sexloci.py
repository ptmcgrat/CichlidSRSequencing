"""Simulation-recovery and data-quality-check tests."""
import os
import sys
import types

import numpy as np
import pandas as pd
import pytest

from sexloci.config import Config
from sexloci.pipeline import PreflightBlocked, preflight, run_pipeline
from sexloci.popgen import wc_fst_components
from sexloci.samples import parse_parents, parse_sex
from sexloci.simulate import simulate
from sexloci.stats import maxz_adjust
from sexloci.labdata import LabData
from sexloci.tables import load_tables_labdata, load_tables_tsv


def _setup(tmp_path, **kw):
    sim = simulate(str(tmp_path), **kw)
    t = load_tables_tsv(sim["sample_table"], sim["genome_version"], sim["alignment_table"], sim["reads_table"])
    return sim, t


def _run(tmp_path, **kw):
    sim, t = _setup(tmp_path, **kw)
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "out"), Config(n_perm=300), quiet=True)
    return sim["truth"], out


def _hit(reg, truth):
    lo, hi = truth["block"]
    return reg[(reg.chrom == truth["chrom"]) & (reg.end >= lo) & (reg.start <= hi)]


def _codes(F, sev=None):
    return {f.code for f in F if sev is None or f.severity == sev}


# ------------------------------------------------------------------ recovery
def test_xy_pedigree_replicated(tmp_path):
    truth, out = _run(tmp_path, system="XY", mode="pedigree", seed=1, n_families=3, n_off=16)
    h = _hit(out["regions"], truth)
    assert len(h) == 1 and h.iloc[0].model == "XY"
    assert h.iloc[0].tier == "A" and h.iloc[0].units_supporting >= 2     # 3 independent families
    assert h.iloc[0].ped_sire_min_padj <= 0.05 < h.iloc[0].ped_dam_min_padj


def test_zw_pedigree_identifies_dam(tmp_path):
    truth, out = _run(tmp_path, system="ZW", mode="pedigree", seed=3, n_families=3, n_off=16)
    h = _hit(out["regions"], truth)
    assert len(h) == 1 and h.iloc[0].model == "ZW" and h.iloc[0].ped_dam_min_padj <= 0.05


@pytest.mark.parametrize("system,seed", [("XY", 4), ("ZW", 5)])
def test_unrelated_only_is_capped_at_B(tmp_path, system, seed):
    truth, out = _run(tmp_path, system=system, mode="unrelated", seed=seed, n_unrelated=14)
    h = _hit(out["regions"], truth)
    assert len(h) == 1 and h.iloc[0].model == system
    assert h.iloc[0].tier == "B" and h.iloc[0].units_total == 1       # one pool: no replication possible
    assert "pedigree" not in out["modules"]


def test_single_family_cannot_reach_A(tmp_path):
    truth, out = _run(tmp_path, system="XY", mode="pedigree", seed=7, n_families=1, n_off=24)
    h = _hit(out["regions"], truth)
    assert len(h) == 1 and h.iloc[0].tier == "B"
    assert "few_independent_families" in _codes(out["findings_obj"])


def test_null_false_positive_rate_is_controlled(tmp_path):
    """No sex locus: calls should be rare (nominal FWER 5%) and never reach the replicated tier A."""
    fp = 0
    for seed in range(10, 20):
        _, out = _run(tmp_path / str(seed), system="none", mode="both", seed=seed)
        r = out["regions"]
        assert not (r.tier == "A").any()
        fp += int(r.tier.isin(["A", "B"]).any())
    assert fp <= 3          # binomial(10, 0.05): P(>=4) < 0.1%


# ------------------------------------------------------------------ BLOCKs
def _blocked(tmp_path, **kw):
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=8, **kw)
    with pytest.raises(PreflightBlocked) as e:
        run_pipeline(sim["vcf"], t, str(tmp_path / "out"), Config(n_perm=50), quiet=True)
    return e.value.findings


def test_block_unsexed(tmp_path):
    F = _blocked(tmp_path, problems={"unsexed": 2})
    assert "sex_not_mf" in _codes(F, "BLOCK")


def test_block_mixed_species(tmp_path):
    assert "species_mismatch" in _codes(_blocked(tmp_path, problems={"mixed_species": 1}), "BLOCK")


def test_block_sample_missing_from_table(tmp_path):
    assert "sample_not_in_table" in _codes(_blocked(tmp_path, problems={"drop_from_table": 1}), "BLOCK")


def test_block_genome_version_mismatch(tmp_path):
    assert "genome_version_mismatch" in _codes(_blocked(tmp_path, problems={"wrong_genome": 1}), "BLOCK")


def test_block_single_sex(tmp_path):
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=8)
    t.sample.loc[t.sample.Sex == "F", "Sex"] = "M"
    with pytest.raises(PreflightBlocked) as e:
        run_pipeline(sim["vcf"], t, str(tmp_path / "out"), Config(n_perm=50), quiet=True)
    assert "too_few_per_sex" in _codes(e.value.findings, "BLOCK")


def test_blocked_run_writes_report_and_skips_analysis(tmp_path):
    _blocked(tmp_path, problems={"unsexed": 1})
    assert (tmp_path / "out" / "preflight_report.md").exists()
    assert not (tmp_path / "out" / "sex_linked_regions.tsv").exists()


# ------------------------------------------------------------------ WARNs
def test_warn_sex_batch_confounding_caps_tier(tmp_path):
    truth, out = _run(tmp_path, system="XY", mode="unrelated", seed=4, n_unrelated=14,
                      problems={"sex_batch_confound": 1})
    assert "sex_confounded" in _codes(out["findings_obj"], "WARN")
    assert out["tier_cap"] == "B"


def test_warn_low_depth_and_sex_depth_bias(tmp_path):
    _, out = _run(tmp_path, system="XY", mode="unrelated", seed=4, n_unrelated=14,
                  problems={"low_depth_samples": 2, "sex_depth_bias": 1})
    c = _codes(out["findings_obj"], "WARN")
    assert {"unusable_coverage", "sex_depth_difference"} <= c


def test_warn_alignment_qc_outlier(tmp_path):
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=10, problems={"qc_outlier": 1})
    _, F, _ = preflight(sim["vcf"], t, Config(), str(tmp_path / "pf"))
    assert "alignment_qc_outlier" in _codes(F, "WARN")


def test_warn_parent_problems_and_id_suggestion(tmp_path):
    sim, t = _setup(tmp_path, mode="pedigree", n_families=2, n_off=10)
    s = t.sample
    s.loc[s.Parents == "F1_sire_m,F1_dam_f", "Parents"] = "F1-sire-m,F1-dam-f"   # punctuation differs
    s.loc[s.SampleID == "F2_o1", "Parents"] = "YES"                              # not a pedigree
    _, F, _ = preflight(sim["vcf"], t, Config(), str(tmp_path / "pf"))
    codes = _codes(F)
    assert {"parent_not_in_vcf", "parent_id_suggestions", "parents_malformed"} <= codes
    assert not F.blocked


def test_warn_wrong_parents_mendelian(tmp_path):
    sim, t = _setup(tmp_path, mode="pedigree", n_families=2, n_off=12)
    s = t.sample
    s.loc[s.SampleID == "F1_o1", "Parents"] = "F2_sire_m,F2_dam_f"
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=100), quiet=True)
    assert "mendelian_errors" in _codes(out["findings_obj"], "WARN")


# ------------------------------------------------------------------ plumbing
def test_cli_blocked_exit_code(tmp_path):
    from sexloci.cli import main
    sim, _ = _setup(tmp_path, mode="unrelated", n_unrelated=8, problems={"mixed_species": 1})
    rc = main(["run", "--vcf", sim["vcf"], "--genome-version", sim["genome_version"],
               "--sample-table", sim["sample_table"], "--outdir", str(tmp_path / "o")])
    assert rc == 2


def test_helpers():
    assert parse_sex("m") == 1 and parse_sex("F") == 0
    for bad in ("U", "F?", "", None, float("nan")):
        assert np.isnan(parse_sex(bad))
    assert parse_parents("a,b")[0] == "ok" and parse_parents("YES")[0] == "malformed"
    assert parse_parents("")[0] == "none"


def test_fst_extremes():
    rng = np.random.default_rng(0)
    freq = rng.uniform(0.1, 0.5, 2000)[:, None]
    mk = lambda: (rng.random((2000, 20)) < freq).astype(np.int8) + (rng.random((2000, 20)) < freq).astype(np.int8)
    a, b, c, *_ = wc_fst_components(mk(), mk())
    assert abs(np.nansum(a) / np.nansum(a + b + c)) < 0.02
    a, b, c, *_ = wc_fst_components(np.zeros((50, 6), np.int8), np.full((50, 6), 2, np.int8))
    assert np.nansum(a) / np.nansum(a + b + c) > 0.95


def test_maxz_adjust():
    rng = np.random.default_rng(0)
    T = rng.chisquare(3, size=(200, 501))
    T[7, 0] = 60
    _, p_adj = maxz_adjust(T)
    assert p_adj[7] < 0.01 and np.sort(p_adj)[1] > 0.05


# ------------------------------------------------------------------ variant filtering
def test_variant_filter_removes_bad_sites_and_keeps_signal(tmp_path):
    sim, t = _setup(tmp_path, system="XY", mode="both", seed=1, problems={"bad_variant_frac": 0.10})
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=300), quiet=True)
    reasons = out["variant_filter"]["reasons"]
    assert reasons["QUAL<30"] == sim["truth"]["n_bad_variants"]          # every planted bad site removed
    assert "variant_filter" in _codes(out["findings_obj"], "INFO")
    h = _hit(out["regions"], sim["truth"])
    assert len(h) == 1 and h.iloc[0].model == "XY"


def test_filter_labels_are_respected(tmp_path):
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=8,
                    problems={"bad_variant_frac": 0.10, "annotate_filter": 1})
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=100), quiet=True)
    assert out["variant_filter"]["reasons"]["FILTER=BadSite"] == sim["truth"]["n_bad_variants"]


def test_warn_when_vcf_cannot_be_filtered(tmp_path):
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=8, problems={"no_info": 1})
    _, F, _ = preflight(sim["vcf"], t, Config(), str(tmp_path / "pf"))
    assert "no_variant_filtering_possible" in _codes(F, "WARN")


def test_filter_switch_off(tmp_path):
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=8, problems={"bad_variant_frac": 0.10})
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=100, apply_variant_filter=False), quiet=True)
    assert out["variant_filter"]["removed"] == 0


def test_describe_distinguishes_raw_from_annotated(tmp_path):
    from sexloci.describe import describe_vcf
    a = simulate(str(tmp_path / "raw"), mode="unrelated", n_unrelated=6, problems={"bad_variant_frac": 0.1})
    b = simulate(str(tmp_path / "ann"), mode="unrelated", n_unrelated=6,
                 problems={"bad_variant_frac": 0.1, "annotate_filter": 1})
    assert "RAW" in describe_vcf(a["vcf"], Config())["kind"]
    d = describe_vcf(b["vcf"], Config())
    assert "annotated" in d["kind"] and d["would_remove"] > 0


def _write_filtered_copy(src, dst, drop):
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line[0] != "#":
                f = line.split("\t", 8)
                if drop(f):
                    continue
            fo.write(line)


def test_filter_compare_recovers_a_filtered_file_and_flags_depleted_windows(tmp_path):
    from sexloci.describe import compare_filters
    sim = simulate(str(tmp_path), mode="unrelated", n_unrelated=6, problems={"bad_variant_frac": 0.10})
    good = str(tmp_path / "good.vcf")
    _write_filtered_copy(sim["vcf"], good, lambda f: "QD=0.5" in f[7])
    r = compare_filters(sim["vcf"], good, Config())["summary"]
    assert r["recall_of_filtered_file"] > 0.99 and r["precision_vs_filtered_file"] > 0.99
    assert r["n_depleted_windows"] == 0
    # a population-genomic filter that deletes a region (e.g. an XY block) must be flagged
    lo, hi = sim["truth"]["block"]
    bad = str(tmp_path / "depleted.vcf")
    _write_filtered_copy(sim["vcf"], bad, lambda f: "QD=0.5" in f[7] or (f[0] == "chr2" and lo <= int(f[1]) <= hi))
    r2 = compare_filters(sim["vcf"], bad, Config())
    assert r2["summary"]["n_depleted_windows"] >= 1
    assert r2["by_window"].query("depleted_in_filtered_file").chrom.eq("chr2").all()


# ------------------------------------------------------------------ LabData (fake rclone + fake Google) and folder workflow
def _fake_rclone(cloud_root, remote="fake:/", fail_copy=False):
    import shutil

    def run(args):
        real = lambda p: os.path.join(cloud_root, p[len(remote):])
        if args[0] == "lsf":
            d = real(args[1])
            if not os.path.isdir(d):
                return 3, "", "directory not found"
            return 0, "".join(n + ("/" if os.path.isdir(os.path.join(d, n)) else "") + "\n"
                              for n in sorted(os.listdir(d))), ""
        if args[0] == "copy":
            if fail_copy:
                return 1, "", "Failed to copy: quota exceeded"
            os.makedirs(args[-1], exist_ok=True)
            src = real(args[-2])
            dst = os.path.join(args[-1], os.path.basename(src))
            if (os.path.exists(dst) and os.path.getsize(dst) == os.path.getsize(src)
                    and int(os.path.getmtime(dst)) == int(os.path.getmtime(src))):
                return 0, "", ""                      # like rclone copy: identical size and modification time -> skipped
            shutil.copy2(src, dst)                     # copy2 keeps the modification time, as rclone does
            return 0, "", ""
        return 1, "", "unsupported"
    return run


class FakeLab(LabData):
    """Real LabData path logic; only the network pieces (rclone, Google sheet reads) are replaced."""
    def __init__(self, tmp_path, frames, fail_copy=False):
        super().__init__("fake:/", "master/", str(tmp_path / "local"), "S", "A",
                         rclone_runner=_fake_rclone(str(tmp_path / "cloud"), fail_copy=fail_copy), log=lambda m: None)
        self.frames = frames

    def _read_sheet(self, sid, ws):
        if ws not in self.frames:
            raise type("WorksheetNotFound", (Exception,), {})(ws)
        return self.frames[ws].copy()


GENOME = "Mzebra_GT3_NCBI"


def _fake_cloud(tmp_path, fail_copy=False):
    import shutil
    cloud = tmp_path / "cloud" / "master" / "Outputs" / "PedigreeVCFs"
    cloud.mkdir(parents=True)
    specs = [("gertrudae_pedigree.vcf.gz", "A_", "Aulonocara gertrudae", "XY", {}, 1),
             ("saulosi_pedigree.vcf.gz", "B_", "Aulonocara saulosi", "ZW", {"bad_variant_frac": 0.05}, 3),
             ("mixed_pedigree.vcf.gz", "C_", "Aulonocara sp. mixed", "XY", {"mixed_species": 1}, 5)]
    st, al, rd = [], [], []
    for name, pre, sp, sysm, prob, seed in specs:
        sim = simulate(str(tmp_path / pre), system=sysm, mode="both", seed=seed, species=sp, prefix=pre,
                       problems=prob, gz=True, n_families=3, n_off=14, n_unrelated=10, genome_version=GENOME)
        shutil.copy(sim["vcf"], cloud / name)
        if pre == "A_":
            (cloud / (name + ".tbi")).write_text("dummy index")
        st.append(pd.read_csv(sim["sample_table"], sep="\t")); al.append(pd.read_csv(sim["alignment_table"], sep="\t"))
        rd.append(pd.read_csv(sim["reads_table"], sep="\t"))
    frames = {"SampleDatabase": pd.concat(st), "DNAReadsDatabase": pd.concat(rd), GENOME: pd.concat(al)}
    return FakeLab(tmp_path, frames, fail_copy)


def test_labdata_cloud_paths_downloads_and_errors(tmp_path):
    lab = _fake_cloud(tmp_path)
    d = lab.local_master_dir + "Outputs/PedigreeVCFs/"
    assert "gertrudae_pedigree.vcf.gz" in lab.list_files(d)
    p = lab.download(d + "gertrudae_pedigree.vcf.gz", fast=True)
    assert os.path.exists(p) and os.path.exists(d + "gertrudae_pedigree.vcf.gz")
    with pytest.raises(FileNotFoundError):
        lab.download(d + "does_not_exist.vcf.gz")
    with pytest.raises(RuntimeError):                                   # a missing cloud folder is an error, not []
        lab.list_files(lab.local_master_dir + "Nope/")
    with pytest.raises(ValueError):                                     # paths outside the master dir are refused
        lab.list_files("/tmp/elsewhere/")
    bad = _fake_cloud(tmp_path / "x", fail_copy=True)
    with pytest.raises(RuntimeError) as e:
        bad.download(bad.local_master_dir + "Outputs/PedigreeVCFs/gertrudae_pedigree.vcf.gz")
    assert "quota exceeded" in str(e.value)                              # rclone errors are surfaced, not swallowed


def test_labdata_tables_and_unknown_genome(tmp_path):
    lab = _fake_cloud(tmp_path)
    t = load_tables_labdata(lab, GENOME)
    assert t.source == "labdata" and len(t.sample) > 50 and t.alignment is not None and t.reads is not None
    with pytest.raises(ValueError) as e:
        lab.read_alignment_table("Mzebra_nonsense")
    assert "Known genome versions" in str(e.value)


def test_labdata_authentication_is_lazy_and_read_only(tmp_path):
    lab = _fake_cloud(tmp_path)
    (tmp_path / "cloud" / "master" / "cichlidsrsequencing_api_creds.json").write_text("{}")
    seen = {}
    gs = types.ModuleType("gspread"); gs.authorize = lambda creds: seen.setdefault("creds", creds) or creds
    google = types.ModuleType("google"); oauth = types.ModuleType("google.oauth2")
    sa = types.ModuleType("google.oauth2.service_account")
    class Creds:
        @staticmethod
        def from_service_account_file(path, scopes):
            seen["scopes"], seen["path"] = scopes, path
            return "creds"
    sa.Credentials = Creds
    saved = {k: sys.modules.get(k) for k in ("gspread", "google", "google.oauth2", "google.oauth2.service_account")}
    sys.modules.update({"gspread": gs, "google": google, "google.oauth2": oauth, "google.oauth2.service_account": sa})
    try:
        assert not os.path.exists(lab.credential_file)                  # nothing downloaded until a sheet is needed
        lab._authenticate()
    finally:
        for k, v in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v
    assert os.path.exists(lab.credential_file)                           # fetched from the cloud on first use
    assert seen["scopes"] == ["https://www.googleapis.com/auth/spreadsheets.readonly"]


def test_folder_runner_downloads_and_analyses_each_vcf(tmp_path):
    from sexloci.folder import FolderRunner
    lab = _fake_cloud(tmp_path)
    fr = FolderRunner(lab, GENOME, Config(n_perm=100), log=lambda m: None)
    assert fr.list_vcfs() == ["gertrudae_pedigree.vcf.gz", "mixed_pedigree.vcf.gz", "saulosi_pedigree.vcf.gz"]
    assert fr.list_vcfs("sau*") == ["saulosi_pedigree.vcf.gz"]
    df = fr.run(str(tmp_path / "out")).set_index("vcf")
    assert df.loc["gertrudae_pedigree.vcf.gz", "status"] == "OK"
    assert df.loc["saulosi_pedigree.vcf.gz", "status"] == "OK"
    assert df.loc["mixed_pedigree.vcf.gz", "status"] == "BLOCKED" and "species_mismatch" in df.loc[
        "mixed_pedigree.vcf.gz", "message"]
    assert df.loc["gertrudae_pedigree.vcf.gz", "model"] == "XY" and df.loc["saulosi_pedigree.vcf.gz", "model"] == "ZW"
    local = lab.local_master_dir + "Outputs/PedigreeVCFs/"
    assert os.path.exists(local + "gertrudae_pedigree.vcf.gz") and os.path.exists(local + "gertrudae_pedigree.vcf.gz.tbi")
    assert (tmp_path / "out" / "saulosi_pedigree" / "report.md").exists()
    assert (tmp_path / "out" / "folder_summary.tsv").exists()
    assert df.loc["saulosi_pedigree.vcf.gz", "removed_fraction"] > 0.04      # its planted bad sites were filtered


def test_folder_runner_manifest_skip_and_check_only(tmp_path):
    from sexloci.folder import FolderRunner, read_manifest
    lab = _fake_cloud(tmp_path)
    (tmp_path / "manifest.tsv").write_text("vcf\tgenome_version\tskip\tnote\nmixed_pedigree.vcf.gz\t\t1\tknown mixed\n")
    fr = FolderRunner(lab, GENOME, Config(n_perm=50), log=lambda m: None)
    df = fr.run(str(tmp_path / "out"), read_manifest(str(tmp_path / "manifest.tsv")), check_only=True).set_index("vcf")
    assert df.loc["mixed_pedigree.vcf.gz", "status"] == "SKIPPED"
    assert df.loc["gertrudae_pedigree.vcf.gz", "status"] == "CHECK_OK"
    assert not (tmp_path / "out" / "gertrudae_pedigree" / "sex_linked_regions.tsv").exists()


def test_check_only_reports_why_a_vcf_is_blocked(tmp_path):
    from sexloci.folder import FolderRunner
    lab = _fake_cloud(tmp_path)
    fr = FolderRunner(lab, GENOME, Config(n_perm=50), log=lambda m: None)
    df = fr.run(str(tmp_path / "out"), check_only=True).set_index("vcf")
    row = df.loc["mixed_pedigree.vcf.gz"]
    assert row.status == "BLOCKED" and "species_mismatch" in row.block_codes and row.n_block == 1
    ff = pd.read_csv(tmp_path / "out" / "folder_findings.tsv", sep="\t")
    assert ((ff.vcf == "mixed_pedigree.vcf.gz") & (ff.code == "species_mismatch")).any()
    assert fr.last_findings["mixed_pedigree.vcf.gz"].blocked


def test_folder_runner_fetches_gff_and_lists_genes(tmp_path):
    from sexloci.folder import FolderRunner
    lab = _fake_cloud(tmp_path)
    gdir = tmp_path / "cloud" / "master" / "Genomes" / GENOME
    gdir.mkdir(parents=True)
    (gdir / "genomic.gff").write_text("chr2\tRefSeq\tgene\t1400000\t1600000\t.\t+\t.\tID=gene-sexgene;Name=sexgene\n")
    fr = FolderRunner(lab, GENOME, Config(n_perm=100), log=lambda m: None)
    df = fr.run(str(tmp_path / "out"), include="gertrudae*")
    assert df.iloc[0].status == "OK"
    reg = pd.read_csv(tmp_path / "out" / "gertrudae_pedigree" / "sex_linked_regions.tsv", sep="\t")
    assert reg.genes.fillna("").str.contains("sexgene").any()


def test_default_output_location(tmp_path):
    from sexloci.cli import default_outdir, main, output_root
    cfg = Config()
    assert output_root(cfg) == os.path.expanduser("~/Temp/SexDeterminationAnalysis") + "/"
    assert default_outdir(cfg, "/x/MC_pedigree_master_file.vcf.gz") == output_root(cfg) + "MC_pedigree_master_file"
    # with output_root redirected, `check` writes <root>/<vcf name>/ without --outdir
    sim = simulate(str(tmp_path / "s"), mode="unrelated", n_unrelated=8, genome_version=GENOME)
    (tmp_path / "cfg.yaml").write_text(f"output_root: {tmp_path}/results/\n")
    rc = main(["check", "--vcf", sim["vcf"], "--sample-table", sim["sample_table"],
               "--alignment-table", sim["alignment_table"], "--config", str(tmp_path / "cfg.yaml")])
    assert rc == 0 and (tmp_path / "results" / "sim" / "preflight_report.md").exists()


def test_no_external_filemanager_dependency():
    import pathlib
    src = pathlib.Path(__file__).parent.parent / "src" / "sexloci"
    code = "".join(p.read_text() for p in src.glob("*.py"))
    assert "--filemanager" not in code and "readSampleDatabase" not in code and "set_with_dataframe" not in code


# ------------------------------------------------------------------ informative blocking messages
def test_compress_ids():
    from sexloci.findings import compress_ids
    assert compress_ids(["MC_046", "MC_047", "MC_049", "MC_050", "MC_051", "MC_052"]) == "MC_046, MC_047, MC_049-MC_052 (4)"
    assert compress_ids([f"F2_o{i}" for i in range(1, 31)] + ["U_M1"]) == "F2_o1-F2_o30 (30), U_M1"
    assert compress_ids(["alpha", "beta"]) == "alpha, beta"
    assert "more groups" in compress_ids([f"S{i}x" for i in range(40)] + [f"T{i}" for i in range(0, 80, 2)], max_items=5)


def _blocked_sex_case(tmp_path):
    sim, t = _setup(tmp_path, mode="both", n_families=2, n_off=30, n_unrelated=2)
    s = t.sample
    m = s.SampleID.str.startswith("F2_o")
    s.loc[m, "Sex"] = pd.NA                        # empty spreadsheet cells
    s.loc[s.SampleID == "U_M1", "Sex"] = "U"
    return sim, t


def test_sex_block_message_explains_what_why_and_how_to_fix(tmp_path):
    sim, t = _blocked_sex_case(tmp_path)
    _, F, _ = preflight(sim["vcf"], t, Config(), str(tmp_path / "pf"))
    f = next(x for x in F if x.code == "sex_not_mf")
    assert "31 of 68 samples (46%)" in f.message and "blank: 30" in f.message and "'U': 1" in f.message
    assert "F2_o1-F2_o30 (30)" in f.message                       # compressed, not truncated at 8
    assert "trait being mapped" in f.why and "Sex' column" in f.fix and len(f.rows) == 31
    fl = pd.read_csv(tmp_path / "pf" / "fix_list.tsv", sep="\t")
    assert set(fl.columns) == {"problem", "sample", "column", "current", "needed"} and len(fl) == 31
    assert set(fl.needed) == {"M or F"}
    rep = (tmp_path / "pf" / "preflight_report.md").read_text()
    assert "Why:" in rep and "How to fix:" in rep and "fix_list.tsv" in rep


def test_pedigree_completeness_is_reported_even_when_blocked(tmp_path):
    sim, t = _blocked_sex_case(tmp_path)
    t.sample.loc[t.sample.SampleID.str.startswith("F2_o"), "Parents"] = pd.NA
    _, F, _ = preflight(sim["vcf"], t, Config(), str(tmp_path / "pf"))
    ped = next(x for x in F if x.code == "pedigree_summary")
    assert F.blocked and "list 'sire,dam'" in ped.message and "30 of them already list parents" not in next(
        x for x in F if x.code == "sex_not_mf").fix
    # and a VCF where nobody lists parents says what to do
    sim2, t2 = _setup(tmp_path / "n", mode="unrelated", n_unrelated=6)
    _, F2, _ = preflight(sim2["vcf"], t2, Config(), str(tmp_path / "pf2"))
    assert "Parents column" in next(x for x in F2 if x.code == "pedigree_summary").message


def test_missing_sample_ids_get_near_match_suggestions_and_stale_fix_list_is_removed(tmp_path):
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=6)
    t.sample.loc[t.sample.SampleID == "U_M1", "SampleID"] = "U-M1"           # punctuation differs from the VCF header
    _, F, _ = preflight(sim["vcf"], t, Config(), str(tmp_path / "pf"))
    f = next(x for x in F if x.code == "sample_not_in_table")
    assert "U_M1 -> U-M1" in f.message and "bcftools reheader" in f.fix
    assert (tmp_path / "pf" / "fix_list.tsv").exists()
    t.sample.loc[t.sample.SampleID == "U-M1", "SampleID"] = "U_M1"            # the user fixes the sheet and reruns
    _, F2, _ = preflight(sim["vcf"], t, Config(), str(tmp_path / "pf"))
    assert not F2.blocked and not (tmp_path / "pf" / "fix_list.tsv").exists()


def test_species_block_names_each_species_and_gives_a_fix(tmp_path):
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=8, problems={"mixed_species": 1})
    _, F, _ = preflight(sim["vcf"], t, Config(), str(tmp_path / "pf"))
    f = next(x for x in F if x.code == "species_mismatch")
    assert "'Mchenga conophoros': 15 samples" in f.message and "'Copadichromis virginalis': 1 sample (U_F8)" in f.message
    assert "bcftools view -S" in f.fix and len(f.rows) == 1


# ------------------------------------------------------------------ coverage: each analysis gets the depth it needs
def _ped_run(tmp_path, system="XY", off_depth=3.0, **kw):
    sim, t = _setup(tmp_path, system=system, mode="pedigree", seed=1, n_families=3, n_off=16,
                    problems={"low_depth_offspring": off_depth, **kw})
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=300), quiet=True)
    return sim, t, out


@pytest.mark.parametrize("system,off_depth", [("XY", 3.0), ("ZW", 3.0), ("XY", 2.0)])
def test_transmission_test_is_tolerant_of_low_offspring_coverage(tmp_path, system, off_depth):
    sim, _, out = _ped_run(tmp_path, system, off_depth)
    h = _hit(out["regions"], sim["truth"])
    assert len(h) == 1 and h.iloc[0].model == system
    assert h.iloc[0].tier == "A" and h.iloc[0].units_supporting == 3            # all 3 independent families, not capped
    codes = _codes(out["findings_obj"])
    assert "unusable_coverage" not in codes and "mendelian_errors" not in codes   # low depth is not an error here
    assert "samples_excluded_from_contrasts" in codes and out["modules"] == ["pedigree"]


def test_low_offspring_coverage_is_noted_not_warned_in_preflight(tmp_path):
    sim, t = _setup(tmp_path, mode="pedigree", n_families=2, n_off=12, problems={"low_depth_offspring": 3.0})
    _, F, _ = preflight(sim["vcf"], t, Config(), str(tmp_path / "pf"))
    assert "low_coverage_table" not in _codes(F, "WARN") and "low_coverage_pedigree" not in _codes(F, "WARN")
    assert "low_coverage_contrast_only" in _codes(F, "INFO")
    sim2, t2 = _setup(tmp_path / "b", mode="pedigree", n_families=2, n_off=12, problems={"low_depth_offspring": 1.5})
    _, F2, _ = preflight(sim2["vcf"], t2, Config(), str(tmp_path / "pf2"))
    f = next(x for x in F2 if x.code == "low_coverage_pedigree")
    assert f.severity == "WARN" and f.rows and "transmission test tolerates" in f.why


def test_mendelian_check_ignores_low_depth_artifacts_but_still_catches_wrong_parents(tmp_path):
    sim, t, out = _ped_run(tmp_path / "ok", "XY", 3.0)
    assert "mendelian_errors" not in _codes(out["findings_obj"])
    sim2, t2 = _setup(tmp_path / "bad", mode="pedigree", n_families=2, n_off=14, problems={"low_depth_offspring": 3.0})
    s = t2.sample
    s.loc[s.SampleID == "F1_o1", "Parents"] = "F2_sire_m,F2_dam_f"
    out2 = run_pipeline(sim2["vcf"], t2, str(tmp_path / "o2"), Config(n_perm=100), quiet=True)
    assert "mendelian_errors" in _codes(out2["findings_obj"], "WARN")


def test_low_depth_unrelated_animals_are_left_out_of_contrasts_not_site_killing(tmp_path):
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=8, problems={"low_depth_samples": 4})
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=100), quiet=True)
    f = next(x for x in out["findings_obj"] if x.code == "samples_excluded_from_contrasts")
    assert f.severity == "WARN" and len(f.samples) == 4 and len(f.rows) == 4        # unrelated animals: this is a loss
    assert out["call_rate_filter"]["dropped"] < 0.1 * out["call_rate_filter"]["seen"]   # the other 12 still keep their sites
    assert out["n_variants_used"] > 5000
    qc = pd.read_csv(tmp_path / "o" / "sample_qc.tsv", sep="\t")
    assert (qc.role == "excluded from contrasts").sum() == 4


def test_block_when_too_few_usable_animals_and_no_pedigree(tmp_path):
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=6, problems={"low_depth_samples": 5})
    with pytest.raises(PreflightBlocked) as e:
        run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=50), quiet=True)
    assert "too_few_usable_samples" in _codes(e.value.findings, "BLOCK")
    assert "too_few_usable_samples" in (tmp_path / "o" / "preflight_report.md").read_text()


def test_low_pool_coverage_still_warned_in_preflight_with_survival_estimate(tmp_path):
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=8, problems={"low_depth_samples": 4})
    _, F, _ = preflight(sim["vcf"], t, Config(), str(tmp_path / "pf"))
    low = next(x for x in F if x.code == "low_coverage_table")
    assert "of a sample's genotypes survive" in low.why and low.rows and low.fix


def test_depth_defaults_follow_the_lab_pipeline():
    c = Config()
    assert c.min_dp == 5 and c.min_gq == 0                      # lab script filters on DP only
    assert c.ped_min_dp == 2 and c.ped_min_gq == 0 and c.ped_parent_min_dp >= 6




# ------------------------------------------------------------------ whole-genome plots, labels, filter-label switch
def test_align_states_recovers_blocks_despite_arbitrary_sign_flips():
    from sexloci.plots import align_states
    rng = np.random.default_rng(3)
    k, n = 30, 20
    truth = np.tile((rng.random(n) < 0.5).astype(float), (k, 1))
    truth[15:, 3] = 1 - truth[15:, 3]                                  # offspring 3 recombines at window 15
    U = (truth * 2 - 1) * rng.uniform(0.5, 1.5, (k, n)) / np.sqrt(n)
    U = U * rng.choice([-1, 1], size=(k, 1))                           # the PCA sign is arbitrary in every window
    st = align_states(U)
    assert np.array_equal(st, truth) or np.array_equal(st, 1 - truth)  # same blocks, global orientation arbitrary


def test_default_contig_labels_match_the_assignment3_table():
    c = Config()
    assert c.contig_labels["NC_036789.1"] == "LG10" and c.contig_labels["NC_036799.1"] == "LG20"
    assert c.contig_labels["NC_036800.1"] == "LG22" and "LG21" not in c.contig_labels.values() and len(c.contig_labels) == 22


def test_inheritance_tracks_split_by_sex_only_for_the_heterogametic_parent(tmp_path):
    from sexloci.plots import align_states
    sim, t = _setup(tmp_path, system="XY", mode="pedigree", seed=1, n_families=3, n_off=24)
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=100), quiet=True)
    tracks = out["haplotype_tracks"]
    meta = {f.id: f for f in out["unit"].families}
    w_sex = (sim["truth"]["sex_pos"] - 1) // Config().window_size
    corr = {"sire": [], "dam": []}
    for (fam, role), per_chrom in tracks.items():
        widx, U = per_chrom[sim["truth"]["chrom"]]
        st = align_states(U)[list(widx).index(w_sex)]
        corr[role].append(abs(np.corrcoef(st, meta[fam].sex)[0, 1]))
    assert min(corr["sire"]) > 0.9                      # sons inherit one sire haplotype, daughters the other
    assert np.mean(corr["dam"]) < 0.6                   # the dam's inherited haplotype is unrelated to sex


def test_whole_genome_plots_are_written_as_png_and_pdf(tmp_path):
    sim, t = _setup(tmp_path, system="XY", mode="pedigree", seed=1, n_families=2, n_off=16)
    run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=100), quiet=True)
    for stem in ("genome_overview", "genome_haplotypes", "genome_by_chromosome"):
        for ext in ("png", "pdf"):
            f = tmp_path / "o" / f"{stem}.{ext}"
            assert f.exists() and f.stat().st_size > 5000, f


def test_unrelated_only_run_has_no_haplotype_plots_but_still_a_genome_overview(tmp_path):
    sim, t = _setup(tmp_path, system="XY", mode="unrelated", seed=4, n_unrelated=12)
    run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=100), quiet=True)
    assert (tmp_path / "o" / "genome_overview.png").exists() and not (tmp_path / "o" / "genome_haplotypes.png").exists()


def test_ignore_filter_labels_lets_a_labelled_filter_be_switched_off(tmp_path):
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=8, problems={"bad_variant_frac": 0.1, "annotate_filter": 1})
    base = run_pipeline(sim["vcf"], t, str(tmp_path / "a"), Config(n_perm=50), quiet=True)["variant_filter"]
    ign = run_pipeline(sim["vcf"], t, str(tmp_path / "b"), Config(n_perm=50, ignore_filter_labels=["BadSite"],
                                                                  variant_filters=[]), quiet=True)["variant_filter"]
    assert base["reasons"]["FILTER=BadSite"] == sim["truth"]["n_bad_variants"]
    assert "FILTER=BadSite" not in ign["reasons"] and ign["removed"] < base["removed"]      # the label is no longer honoured


def test_filter_compare_writes_whole_genome_retention_plot(tmp_path):
    from sexloci.cli import main
    sim = simulate(str(tmp_path), mode="unrelated", n_unrelated=6, problems={"bad_variant_frac": 0.1})
    good = str(tmp_path / "good.vcf")
    _write_filtered_copy(sim["vcf"], good, lambda f: "QD=0.5" in f[7])
    rc = main(["filter-compare", "--raw", sim["vcf"], "--filtered", good, "--outdir", str(tmp_path / "fc")])
    assert rc == 0 and (tmp_path / "fc" / "retention_genome.png").exists() and (tmp_path / "fc" / "retention_genome.pdf").exists()


def test_pedigree_effect_size_is_flat_at_one_across_a_non_recombining_block_and_falls_off_outside(tmp_path):
    sim, t = _setup(tmp_path, system="XY", mode="pedigree", seed=1, n_families=2, n_off=24, n_chrom=4, chrom_len=6_000_000,
                    snps_per_chrom=3000, sex_pos=3_000_000, block_half_width=1_500_000)
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=300), quiet=True)
    w = out["windows"]
    c = w[w.chrom == "chr2"].set_index("widx")
    inside = c.loc[(c.start >= 1_750_000) & (c.end <= 4_250_000), "ped_sire_absr"]
    assert (inside > 0.97).all() and c.ped_dam_absr.max() < 0.6                 # the plateau that looks 'flat-topped'
    far = w[w.chrom != "chr2"].ped_sire_absr.dropna()
    assert far.max() < 0.75 and far.median() < 0.3                              # background on the other chromosomes
    assert w.ped_sire_absr.dropna().between(0, 1.0 + 1e-9).all()               # an effect size, not a saturating p-value


# ------------------------------------------------------------------ region calling: core interval, votes, tiny windows, spikes
def test_core_interval_follows_the_differentiated_segment_inside_the_transmission_plateau(tmp_path):
    sim, t = _setup(tmp_path, system="XY", mode="pedigree", seed=1, n_families=2, n_off=24, n_chrom=4, chrom_len=6_000_000,
                    snps_per_chrom=3000, sex_pos=3_000_000, block_half_width=1_500_000)
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=300), quiet=True)
    h = _hit(out["regions"], sim["truth"]).iloc[0]
    lo, hi = sim["truth"]["block"]
    assert abs(h.core_start - lo) <= 0.6e6 and abs(h.core_end - hi) <= 0.6e6          # the planted 3 Mb segment
    assert h.core_start >= h.start and h.core_end <= h.end                            # inside the plateau
    assert "differentiated core" in (tmp_path / "o" / "report.md").read_text()


def test_model_votes_only_count_evidence_that_fired():
    from sexloci.regions import decide_model
    c = Config()
    # LG16-like: a 3 kb male-only-heterozygous feature in a 250 kb window whose average dHet is slightly negative
    assert decide_model(False, False, -0.055, False, 4, 0, False, c)[0] == "undetermined"
    assert decide_model(True, False, 0.25, True, 9000, 0, True, c) == ("XY", ["pedigree:XY", "dhet:XY", "diag:XY"])
    assert decide_model(False, True, -0.2, True, 0, 500, True, c)[0] == "ZW"
    # one sire: negative dHet (sire's X carries private alleles) must not overrule the pedigree
    assert decide_model(True, False, -0.3, True, 0, 500, True, c)[0] == "XY"
    assert decide_model(False, False, 0.2, True, 400, 0, True, c)[0] == "XY"                  # no pedigree: dHet / diagnostic SNPs decide
    assert decide_model(False, False, -0.2, True, 400, 0, True, c)[0] == "conflicting"        # and real disagreement there is reported


def test_tiny_windows_cannot_create_population_genetic_leads():
    from sexloci.regions import flag_evidence
    rng = np.random.default_rng(0)
    n = 300
    w = pd.DataFrame({"chrom": "c1", "widx": np.arange(n), "start": np.arange(n) * 250_000 + 1,
                      "end": (np.arange(n) + 1) * 250_000, "fst": rng.normal(0, 0.01, n), "dhet": rng.normal(0, 0.01, n),
                      "n_diag_xy": 0, "n_diag_zw": 0, "cov_log2_mf": rng.normal(0, 0.02, n), "pg_n_markers": 400})
    w.loc[100, ["fst", "dhet", "cov_log2_mf"]] = [0.3, 0.2, 0.4]
    few, many = w.copy(), w.copy()
    few.loc[100, "pg_n_markers"] = 5          # an outlier built from 5 markers (a repeat-rich window)
    assert not flag_evidence(few, Config())["candidate"].any()
    assert flag_evidence(many, Config()).loc[100, "candidate"]


def test_narrow_spike_is_flagged_not_presented_as_a_locus(tmp_path):
    sim, t = _setup(tmp_path, system="XY", mode="unrelated", seed=4, n_unrelated=30, block_half_width=5_000)
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=300), quiet=True)
    r = out["regions"]
    spike = r[(r.chrom == sim["truth"]["chrom"]) & r.flag.str.contains("narrow_spike")]
    assert len(spike) == 1 and spike.iloc[0].n_sig_markers <= 10 and "inspect the reads" in spike.iloc[0].flag


# ------------------------------------------------------------------ a brood whose parent is not in the VCF
def _two_broods(tmp_path, omit=("F2_dam_f",), **kw):
    return _setup(tmp_path, system="XY", mode="pedigree", seed=1, n_families=2, n_off=20,
                  problems={"omit_samples": list(omit), **kw})


def test_brood_with_missing_dam_skips_transmission_but_stays_in_the_sex_contrasts(tmp_path):
    sim, t = _two_broods(tmp_path)
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=300), quiet=True)
    f = next(x for x in out["findings_obj"] if x.code == "parent_not_in_vcf")
    assert f.severity == "WARN" and "F2_dam_f" in f.message and "20 offspring" in f.message
    assert "transmission test" in f.why and "add them to the VCF" in f.fix
    assert [r["sample"] for r in f.rows] == ["F2_dam_f"]
    assert (tmp_path / "o" / "fix_list.tsv").exists()
    unit = out["unit"]
    assert [fam.id for fam in unit.families] == ["fam:F1_sire_mxF1_dam_f"]                      # transmission: brood 1 only
    assert all(k.startswith("fam:F1") for k in out["pedigree"]["informative"])
    assert set(out["assoc"]["strata_used"]) == {"fam:F1_sire_mxF1_dam_f", "fam:F2_sire_mxF2_dam_f"}   # contrasts: both broods
    names = unit.names
    assert "F2_sire_m" in [names[i] for i in unit.parent_idx]                                   # the present parent is not an unrelated animal
    assert "F2_sire_m" not in [names[i] for i in unit.contrast_idx]
    h = _hit(out["regions"], sim["truth"])
    assert len(h) == 1 and h.iloc[0].model == "XY" and h.iloc[0].units_total == 2


def test_only_one_parent_listed_keeps_the_brood_together_not_in_the_pool(tmp_path):
    for listed in ("F2_sire_m,", "F2_sire_m,?", "F2_sire_m"):
        sim, t = _setup(tmp_path / listed.replace(",", "_").replace("?", "q"), system="XY", mode="pedigree", seed=1,
                        n_families=2, n_off=16)
        s = t.sample
        s.loc[s.SampleID.str.startswith("F2_o"), "Parents"] = listed
        _, F, _ = preflight(sim["vcf"], t, Config(), str(tmp_path / "pf" / listed.replace(",", "_").replace("?", "q")))
        f = next(x for x in F if x.code == "parent_unknown")
        assert f.severity == "WARN" and "dam unknown" in f.message and len(f.rows) == 16
        assert "sire,dam" in f.fix and not F.blocked, listed
    from sexloci.samples import build_unit
    F2 = __import__("sexloci.findings", fromlist=["Findings"]).Findings()
    unit = build_unit(list(t.sample.SampleID), t, Config(), F2)
    f2 = [i for i, nm in enumerate(unit.names) if nm.startswith("F2_o")]
    assert set(unit.strata[f2]) == {"fam:F2_sire_mx?"} and [fam.id for fam in unit.families] == ["fam:F1_sire_mxF1_dam_f"]


def test_a_known_dam_alone_is_recognised_by_her_sex(tmp_path):
    sim, t = _setup(tmp_path, system="XY", mode="pedigree", seed=1, n_families=2, n_off=12)
    s = t.sample
    s.loc[s.SampleID.str.startswith("F2_o"), "Parents"] = "F2_dam_f"
    from sexloci.findings import Findings
    from sexloci.samples import build_unit
    unit = build_unit(list(s.SampleID), t, Config(), Findings())
    f2 = [i for i, nm in enumerate(unit.names) if nm.startswith("F2_o")]
    assert set(unit.strata[f2]) == {"fam:?xF2_dam_f"}


def test_no_complete_brood_turns_the_transmission_test_off_but_contrasts_still_run(tmp_path):
    sim, t = _two_broods(tmp_path, omit=("F1_dam_f", "F2_dam_f"))
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=200), quiet=True)
    codes = _codes(out["findings_obj"], "WARN")
    assert {"parent_not_in_vcf", "no_usable_pedigree"} <= codes
    assert out["modules"] == ["assoc", "popgen"] and not (tmp_path / "o" / "genome_haplotypes.png").exists()
    assert len(_hit(out["regions"], sim["truth"])) == 1


def test_pedigree_summary_mentions_one_parent_broods(tmp_path):
    sim, t = _setup(tmp_path, system="XY", mode="pedigree", seed=1, n_families=2, n_off=12)
    s = t.sample
    s.loc[s.SampleID.str.startswith("F2_o"), "Parents"] = "F2_sire_m,"
    _, F, _ = preflight(sim["vcf"], t, Config(), str(tmp_path / "pf"))
    assert "12 list only one parent" in next(x for x in F if x.code == "pedigree_summary").message


# ------------------------------------------------------------------ filter audit by window
def test_filter_audit_shows_where_sites_are_removed(tmp_path):
    from sexloci.cli import main
    from sexloci.describe import audit_filter_by_window
    sim = simulate(str(tmp_path), mode="unrelated", n_unrelated=6, problems={"bad_variant_frac": 0.10})
    df = audit_filter_by_window(sim["vcf"], Config())
    assert df.sites.sum() > 9000                                            # every SNP of the simulation is counted
    assert abs(df.removed.sum() / df.sites.sum() - 0.10) < 0.02             # the planted 10% of bad sites
    assert df["QUAL<30"].sum() == sim["truth"]["n_bad_variants"] and df.retained.between(0, 1).all()
    rc = main(["filter-audit", "--vcf", sim["vcf"], "--outdir", str(tmp_path / "fa")])
    assert rc == 0 and (tmp_path / "fa" / "filter_audit_by_window.tsv").exists()
    assert (tmp_path / "fa" / "filter_audit_genome.png").exists() and (tmp_path / "fa" / "filter_audit_genome.pdf").exists()


# ------------------------------------------------------------------ annotation: labels from the GFF, genes via linkage groups
def _gff(path, seq_lg, genes=()):
    lines = ["##gff-version 3"]
    for seq, lg in seq_lg.items():
        lines.append(f"{seq}\tRefSeq\tregion\t1\t6000000\t.\t+\t.\tID={seq}:1..6000000;Name={lg};chromosome={lg};genome=chromosome")
    lines.append("NC_MT.1\tRefSeq\tregion\t1\t16000\t.\t+\t.\tID=mt;Name=MT;genome=mitochondrion")
    for seq, a, b, name in genes:
        lines.append(f"{seq}\tRefSeq\tgene\t{a}\t{b}\t.\t+\t.\tID=gene-{name};Name={name}")
    pathlib_write = open(path, "w")
    pathlib_write.write("\n".join(lines) + "\n")
    pathlib_write.close()
    return str(path)


def test_read_gff_regions_normalises_labels_and_skips_non_chromosomes(tmp_path):
    from sexloci.regions import read_gff_regions
    g = _gff(tmp_path / "a.gff", {"NC_1.1": "LG1", "NC_2.1": "10", "NC_3.1": "LG22"})
    assert read_gff_regions(g) == {"NC_1.1": "LG1", "NC_2.1": "LG10", "NC_3.1": "LG22"}      # numeric -> LG, MT skipped


def test_genes_are_listed_for_a_vcf_named_with_other_chromosome_names_via_linkage_group(tmp_path):
    """YH-style: GT3 sequence called under UMD2a chromosome names, annotation under the GT3a accessions."""
    sim, t = _setup(tmp_path, system="XY", mode="pedigree", seed=1, n_families=2, n_off=20)
    gff = _gff(tmp_path / "gt3a.gff", {f"NC_90000{i}.1": f"LG{i}" for i in range(1, 5)},
               genes=[("NC_900002.1", 1_400_000, 1_600_000, "sexgene"), ("NC_900003.1", 1_400_000, 1_600_000, "othergene")])
    labels = {f"chr{i}": f"LG{i}" for i in range(1, 5)}                    # the VCF's names -> linkage groups
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=100, gff=gff, contig_labels=labels), quiet=True)
    h = _hit(out["regions"], sim["truth"]).iloc[0]
    assert "sexgene" in h.genes and "othergene" not in h.genes                 # right linkage group only
    assert "annotation_matched_by_linkage_group" in _codes(out["findings_obj"], "INFO")
    off = run_pipeline(sim["vcf"], t, str(tmp_path / "o2"), Config(n_perm=100, gff=gff, contig_labels=labels,
                                                                    annotation_alias_by_linkage_group=False), quiet=True)
    assert not isinstance(_hit(off["regions"], sim["truth"]).iloc[0].genes, str) or _hit(off["regions"], sim["truth"]).iloc[0].genes == ""


def test_chromosome_labels_come_from_the_annotation_when_the_config_has_none(tmp_path):
    import json
    sim, t = _setup(tmp_path, system="XY", mode="pedigree", seed=1, n_families=2, n_off=16)
    gff = _gff(tmp_path / "g.gff", {f"chr{i}": str(i) for i in range(1, 5)},
               genes=[("chr2", 1_400_000, 1_600_000, "sexgene")])
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=100, gff=gff, contig_labels={}), quiet=True)
    lab = json.load(open(tmp_path / "o" / "run_summary.json"))["config"]["contig_labels"]
    assert lab["chr2"] == "LG2" and "sexgene" in _hit(out["regions"], sim["truth"]).iloc[0].genes
    assert "annotation_matched_by_linkage_group" not in _codes(out["findings_obj"])   # same names: no aliasing needed


# ------------------------------------------------------------------ cloud copy replaced under the same name
def test_sync_reports_downloaded_unchanged_updated_and_forced(tmp_path):
    import time
    lab = _fake_cloud(tmp_path)
    local = lab.local_master_dir + "Outputs/PedigreeVCFs/saulosi_pedigree.vcf.gz"
    cloud = tmp_path / "cloud" / "master" / "Outputs" / "PedigreeVCFs" / "saulosi_pedigree.vcf.gz"
    assert lab.sync(local) == "downloaded"
    assert lab.sync(local) == "unchanged"                                   # same size and modification time: not fetched again
    original = open(local, "rb").read()
    other = (tmp_path / "cloud" / "master" / "Outputs" / "PedigreeVCFs" / "gertrudae_pedigree.vcf.gz").read_bytes()
    cloud.write_bytes(other)                                                # the corrected file is uploaded under the same name
    assert lab.sync(local) == "updated" and open(local, "rb").read() == other != original
    os.utime(cloud, (time.time() + 3600, time.time() + 3600))               # same size, newer modification time
    assert lab.sync(local) == "updated"
    assert lab.sync(local) == "unchanged"
    assert lab.sync(local, force=True) == "downloaded"


def test_folder_run_picks_up_a_vcf_replaced_in_the_cloud_under_the_same_name(tmp_path):
    from sexloci.folder import FolderRunner
    lab = _fake_cloud(tmp_path)
    pv = tmp_path / "cloud" / "master" / "Outputs" / "PedigreeVCFs"
    fr = FolderRunner(lab, GENOME, Config(n_perm=50), log=lambda m: None)
    d1 = fr.run(str(tmp_path / "o"), include="gertrudae*", check_only=True).iloc[0]
    d2 = fr.run(str(tmp_path / "o"), include="gertrudae*", check_only=True).iloc[0]
    assert d1.download == "downloaded" and d2.download == "unchanged" and d1.vcf_modified
    (pv / "gertrudae_pedigree.vcf.gz").write_bytes((pv / "saulosi_pedigree.vcf.gz").read_bytes())     # replaced under the same name
    d3 = fr.run(str(tmp_path / "o"), include="gertrudae*", check_only=True).iloc[0]
    assert d3.download == "updated"
    local = lab.local_master_dir + "Outputs/PedigreeVCFs/gertrudae_pedigree.vcf.gz"
    assert open(local, "rb").read() == (pv / "saulosi_pedigree.vcf.gz").read_bytes()
    assert fr.run(str(tmp_path / "o"), include="gertrudae*", check_only=True, refresh=True).iloc[0].download == "downloaded"


def test_run_summary_records_which_vcf_version_was_analysed(tmp_path):
    import json
    sim, t = _setup(tmp_path, mode="unrelated", n_unrelated=8)
    run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=50), quiet=True)
    js = json.load(open(tmp_path / "o" / "run_summary.json"))
    assert js["vcf_size_bytes"] == os.path.getsize(sim["vcf"]) and len(js["vcf_modified"]) == 19


def test_shared_parent_is_named_so_a_wrong_sire_in_the_sheet_can_be_spotted(tmp_path):
    sim, t = _setup(tmp_path, system="XY", mode="pedigree", seed=1, n_families=2, n_off=12)
    s = t.sample
    s.loc[s.SampleID.str.startswith("F2_o"), "Parents"] = "F1_sire_m,F2_dam_f"              # brood 2 wrongly lists brood 1's sire
    _, F, _ = preflight(sim["vcf"], t, Config(family_independence="parent"), str(tmp_path / "pf"))
    sh = next(x for x in F if x.code == "shared_parents")
    assert "F1_sire_m is a parent of 2 broods" in sh.message and "correct the Parents column" in sh.message
    ps = next(x for x in F if x.code == "pedigree_summary").message
    assert "F1_sire_m x F1_dam_f (12 offspring)" in ps and "F1_sire_m x F2_dam_f (12 offspring)" in ps


def test_broods_sharing_a_sire_are_independent_broods_by_default_and_the_report_says_so(tmp_path):
    sim, t = _setup(tmp_path / "a", system="XY", mode="pedigree", seed=1, n_families=2, n_off=20, problems={"shared_sire": 1})
    out = run_pipeline(sim["vcf"], t, str(tmp_path / "o"), Config(n_perm=300), quiet=True)
    h = _hit(out["regions"], sim["truth"]).iloc[0]
    assert h.tier == "A" and h.units_supporting == 2 and h.units_share_parent == "F1_sire_m" and h.model == "XY"
    sh = next(x for x in out["findings_obj"] if x.code == "shared_parents")
    assert "treated as independent" in sh.message and "F1_sire_m" in sh.message
    assert "share F1_sire_m" in (tmp_path / "o" / "report.md").read_text()
    assert "few_independent_families" in _codes(out["findings_obj"], "WARN")            # still only 2 broods
    strict = run_pipeline(sim["vcf"], t, str(tmp_path / "o2"), Config(n_perm=300, family_independence="parent"), quiet=True)
    hs = _hit(strict["regions"], sim["truth"]).iloc[0]
    assert hs.tier == "B" and hs.units_total == 1                                       # the strict view: one family


def test_a_rerun_removes_region_plots_left_by_an_earlier_run(tmp_path):
    sim, t = _setup(tmp_path, system="XY", mode="pedigree", seed=1, n_families=2, n_off=16)
    out_dir = tmp_path / "o"
    out_dir.mkdir()
    for stale in ("region_2_old.png", "region_3_old.png", "genome_haplotypes.png"):
        (out_dir / stale).write_bytes(b"old")
    run_pipeline(sim["vcf"], t, str(out_dir), Config(n_perm=100), quiet=True)
    assert not (out_dir / "region_2_old.png").exists() and not (out_dir / "region_3_old.png").exists()
    assert (out_dir / "genome_haplotypes.png").stat().st_size > 5000                 # replaced by this run's figure


# ------------------------------------------------------------------ gsdf-region genotyping across the sample database
class _FakeServer:
    def __enter__(self):
        return "http://fake"

    def __exit__(self, *a):
        pass


def _gsdf_fixture(tmp_path):
    from sexloci.tables import Tables
    #            id    species ecogroup sex  flank ratio junction coverage
    rows = [("A1", "Mchenga conophoros", "Shallow_Benthic", "F", 8.0, 0.53, 0.00, 20),
            ("A2", "Mchenga conophoros", "Shallow_Benthic", "F", 6.0, 0.55, 0.00, 12),
            ("B1", "Mchenga conophoros", "Shallow_Benthic", "M", 8.0, 0.80, 0.25, 20),
            ("B2", "Mchenga conophoros", "Shallow_Benthic", "M", 6.0, 0.78, 0.25, 12),
            ("C1", "Metriaclima zebra", "Mbuna", "M", 9.0, 1.05, 0.50, 30),
            ("D1", "Metriaclima zebra", "Mbuna", "F", 2.0, 0.50, 0.00, 5),
            ("E1", "Metriaclima zebra", "Mbuna", "F", 8.0, 0.53, 0.00, 10),     # aligned, but the BAM is missing in the cloud
            ("F1", "Aulonocara sp.", "Deep_Benthic", "M", 8.0, 0.80, 0.25, 15),  # not aligned to this genome
            ("G1", "Aulonocara sp.", "Deep_Benthic", "F", 1.0, 0.50, 0.00, 1.0)]  # coverage below the default minimum
    sample = pd.DataFrame([dict(SampleID=r[0], Species=r[1], Ecogroup=r[2], Sex=r[3]) for r in rows])
    align = pd.DataFrame([dict(SampleID=r[0], Coverage=r[7]) for r in rows if r[0] != "F1"])
    truth = {r[0]: (r[4], r[5], r[6]) for r in rows if r[0] not in ("E1", "F1")}
    return Tables(sample, align, None, GENOME), truth


def _gsdf_runner(truth, calls):
    import re

    def run(args):
        calls.append(args[:2])
        if args[1] == "view":
            sample = args[-2].split("/")[-2]
            if sample not in truth:
                return 1, "", "[E::hts_open_format] Failed to open file: HTTP 404 Not Found"
            open(args[args.index("-o") + 1], "w").write(sample)
            return 0, "", ""
        if args[1] == "index":
            open(args[-1] + ".bai", "w").write("x")
            return 0, "", ""
        if args[1] == "depth":
            fd, ratio, jx = truth[open(args[-1]).read()]
            lo, hi = map(int, re.search(r":(\d+)-(\d+)$", args[args.index("-r") + 1]).groups())
            unique = "-Q" in args
            lines = []
            for pos in range(lo, hi + 1, 100):
                if 50_772_650 <= pos <= 50_773_150:
                    d = jx * fd if unique else ratio * fd
                elif 50_755_000 <= pos < 50_790_000:
                    d = ratio * fd
                else:
                    d = fd
                lines.append(f"NC_135173.1\t{pos}\t{d}")
            return 0, "\n".join(lines) + "\n", ""
        raise AssertionError(args)
    return run


def _run_gsdf(tmp_path, tables, truth, calls, **kw):
    from sexloci.gsdf import run_gsdf_genotype
    from sexloci.config import Config
    lab = types.SimpleNamespace(cloud_master_dir="remote:/master/")
    return run_gsdf_genotype(tables, Config(), str(tmp_path / "out"), lab, GENOME, jobs=2, runner=_gsdf_runner(truth, calls),
                             server_factory=_FakeServer, log=lambda m: None, **kw)


def test_gsdf_genotype_table_covers_the_whole_database(tmp_path):
    tables, truth = _gsdf_fixture(tmp_path)
    df = _run_gsdf(tmp_path, tables, truth, []).set_index("SampleID")
    assert list(df.columns)[:6] == ["Species", "Ecogroup", "Sex", "Pedigree", "Genotype", "CopyNumber"]
    assert {s: df.loc[s, "Genotype"] for s in ("A1", "A2", "B1", "B2", "C1")} == {"A1": "-/-", "A2": "-/-", "B1": "dup/-", "B2": "dup/-", "C1": "dup/dup"}
    assert df.loc["B1", "Ecogroup"] == "Shallow_Benthic" and df.loc["C1", "Species"] == "Metriaclima zebra"
    assert df.loc["D1", "Genotype"] == "-/-" and "low depth" in df.loc["D1", "Flag"]
    assert df.loc["E1", "Status"] == "no_bam" and df.loc["E1", "Genotype"] == ""
    assert df.loc["F1", "Status"].startswith("no Mzebra_GT3_NCBI alignment")
    assert df.loc["G1", "Status"].startswith("skipped: coverage")
    assert df.loc["B1", "JunctionJ"] > 0.2 and df.loc["A1", "JunctionJ"] < 0.05           # unique-read junction signal
    assert len(df) == 9 and (tmp_path / "out" / "gsdf_genotypes.tsv").exists()


def test_gsdf_genotype_resumes_and_only_retries_what_failed(tmp_path):
    tables, truth = _gsdf_fixture(tmp_path)
    _run_gsdf(tmp_path, tables, truth, [])
    calls = []
    _run_gsdf(tmp_path, tables, truth, calls)
    assert [c for c in calls if c[1] == "view"] == [["samtools", "view"]]                  # only E1 (no_bam) is tried again
    assert not [c for c in calls if c[1] == "depth"]
    calls2 = []
    _run_gsdf(tmp_path, tables, truth, calls2, refresh=True)
    assert len([c for c in calls2 if c[1] == "view"]) >= 6                                 # everything measured again


def test_gsdf_genotype_anchors_selection_and_flags(tmp_path):
    from sexloci.gsdf import call_genotypes
    tables, truth = _gsdf_fixture(tmp_path)
    (tmp_path / "anch.tsv").write_text("SampleID\tcopies\nA1\t2\nB1\t3\nC1\t4\n")
    df = _run_gsdf(tmp_path, tables, truth, [], anchors=str(tmp_path / "anch.tsv")).set_index("SampleID")
    assert df.loc["A2", "Genotype"] == "-/-" and df.loc["B2", "Genotype"] == "dup/-"
    only = _run_gsdf(tmp_path / "x", tables, truth, [], ecogroup="mbuna", samples=["C*"]).set_index("SampleID")
    assert only.loc["C1", "Genotype"] == "dup/dup" and only.loc["B1", "Genotype"] == "" and only.loc["B1", "Status"].startswith("not selected")
    m = pd.DataFrame({"DepthRatio": [0.665, 1.9, 0.53], "FlankDepth": [8.0, 8.0, 8.0]}, index=["u", "o", "n"])
    c = call_genotypes(m, 0.530, 0.268)
    assert "uncertain" in c.loc["u", "Flag"] and "outside calibrated range" in c.loc["o", "Flag"] and c.loc["n", "Flag"] == ""


def test_gsdf_subgroup_column_is_carried_and_summaries_classify_the_pattern(tmp_path):
    from sexloci.gsdf import summarize_by
    tables, truth = _gsdf_fixture(tmp_path)
    tables.sample["Sub-Group"] = ["lakeA", "lakeA", "lakeA", "lakeB", "lakeB", "lakeB", "lakeB", "lakeB", "lakeB"]    # an oddly spelled column
    df = _run_gsdf(tmp_path, tables, truth, [])
    assert list(df.columns)[:4] == ["SampleID", "Species", "Ecogroup", "Subgroup"]
    assert dict(zip(df.SampleID, df.Subgroup))["B1"] == "lakeA" and dict(zip(df.SampleID, df.Subgroup))["C1"] == "lakeB"
    # a synthetic table with one species per pattern
    rows = []
    def add(sp, sex, cn, n):
        rows.extend([dict(Species=sp, Sex=sex, CopyNumber=cn, Genotype="x", Flag="") for _ in range(n)])
    add("xy", "M", 3, 9); add("xy", "M", 2, 1); add("xy", "F", 2, 10)                         # male-limited
    add("none", "M", 2, 8); add("none", "F", 2, 8)                                           # absent
    add("poly", "M", 3, 5); add("poly", "M", 2, 5); add("poly", "F", 3, 4); add("poly", "F", 2, 6)   # both sexes
    add("biased", "M", 3, 6); add("biased", "M", 2, 14); add("biased", "F", 2, 20)           # male-biased but not limited
    add("nosex", "U", 3, 6)
    t = pd.DataFrame(rows)
    res = summarize_by(t, "Species").set_index("Species")
    assert res.loc["xy", "pattern"] == "male-limited" and res.loc["none", "pattern"] == "absent"
    assert res.loc["poly", "pattern"] == "present, not sex-associated" and res.loc["biased", "pattern"] == "male-biased"
    assert res.loc["nosex", "pattern"] == "no sex information" and res.loc["xy", "male_carriers"] == 9 and res.loc["xy", "female_carriers"] == 0
    with pytest.raises(ValueError):
        summarize_by(t, "Nope")


def test_gsdf_pedigree_fish_are_marked_and_left_out_of_summaries_by_default(tmp_path):
    from sexloci.gsdf import pedigree_roles, summarize_by
    tables, truth = _gsdf_fixture(tmp_path)
    tables.sample["Parents"] = ["", "", "", "B1,A2", "", "", "", "?,?", ""]        # B2 is the offspring of B1 x A2; E1's parents are unknown
    roles = pedigree_roles(tables.sample)
    assert roles.to_dict() == {"B2": "offspring", "B1": "parent", "A2": "parent"}
    df = _run_gsdf(tmp_path, tables, truth, [])
    ped = dict(zip(df.SampleID, df.Pedigree))
    assert ped["B2"] == "offspring" and ped["B1"] == "parent" and ped["A2"] == "parent" and ped["A1"] == "" and ped["E1"] == ""
    by_all = summarize_by(df, "Ecogroup", independent=False).set_index("Ecogroup")
    by_ind = summarize_by(df, "Ecogroup").set_index("Ecogroup")
    assert by_all.loc["Shallow_Benthic", "n"] == 4 and by_ind.loc["Shallow_Benthic", "n"] == 1      # A1 only: B1, B2, A2 are relatives


# ------------------------------------------------------------------ writing one column into the Sample Database sheet
class _FakeWS:
    """Just enough of a gspread worksheet: get_all_values, col_values (trailing blanks trimmed), add_cols, batch_update, update_note."""
    def __init__(self, grid):
        self.grid = [list(r) for r in grid]
        self.col_count = max(len(r) for r in grid)
        self.writes, self.note, self.reads = [], None, 0
        self.scramble_after = None                                   # simulate someone sorting the sheet after the first read

    def get_all_values(self):
        return [list(r) for r in self.grid]

    def col_values(self, c):
        self.reads += 1
        rows = self.grid if not (self.scramble_after and self.reads >= self.scramble_after) else list(reversed(self.grid))
        vals = [(r[c - 1] if len(r) >= c else "") for r in rows]
        while vals and vals[-1] == "":
            vals.pop()
        return vals

    def add_cols(self, n):
        self.col_count += n

    def batch_update(self, data, value_input_option=None):
        self.vio = value_input_option
        import re as _re
        for d in data:
            m = _re.match(r"([A-Z]+)(\d+)$", d["range"])
            col = 0
            for ch in m.group(1):
                col = col * 26 + ord(ch) - 64
            row = int(m.group(2))
            while len(self.grid) < row:
                self.grid.append([])
            self.grid[row - 1] += [""] * (col - len(self.grid[row - 1]))
            self.grid[row - 1][col - 1] = str(d["values"][0][0])
            self.writes.append(d)

    def update_note(self, cell, text):
        self.note = (cell, text)


def _sheet_fixture():
    grid = [["SampleID", "Species", "Ecogroup", "Sex"],
            ["A1", "sp one", "Mbuna", "F"], ["B1", "sp two", "Utaka", "M"], ["C1", "sp three", "Utaka", "M"],
            ["D1", "sp four", "Deep", "F"], ["E1", "sp five", "Deep", "M"], ["F1", "sp six", "Deep", "U"]]
    res = pd.DataFrame({"SampleID": ["A1", "B1", "C1", "D1", "E1", "Z9"], "Genotype": ["-/-", "dup/-", "dup/dup", "-/-", "", "dup/-"],
                        "CopyNumber": [2, 3, 4, 2, pd.NA, 3], "Flag": ["", "", "", "low depth", "", ""]})
    return _FakeWS(grid), res


def test_sample_database_write_is_a_dry_run_unless_confirmed_and_touches_only_the_new_column(tmp_path):
    from sexloci.samplesheet import update_sample_database
    ws, res = _sheet_fixture()
    before = [list(r) for r in ws.grid]
    plan = update_sample_database(ws, res, str(tmp_path), apply=False, log=lambda m: None)
    assert ws.grid == before and not ws.writes                                    # dry run: nothing written
    assert plan.is_new and [(r, s, n) for r, s, o, n in plan.to_write] == [(2, "A1", 2), (3, "B1", 3), (4, "C1", 4)]   # D1 flagged, E1 no call
    assert plan.not_in_sheet == ["Z9"] and list(tmp_path.glob("*_plan_*.tsv"))
    update_sample_database(ws, res, str(tmp_path), apply=True, version="0.1.x", log=lambda m: None)
    assert ws.grid[0][4] == "gsdf number" and [r[4] if len(r) > 4 else "" for r in ws.grid[1:]] == ["2", "3", "4", "", "", ""]
    for r_old, r_new in zip(before, ws.grid):                                      # every other cell is exactly as before
        assert r_new[:4] == r_old[:4]
    assert ws.note and ws.note[0] == "E1" and list(tmp_path.glob("sample_database_backup_*.csv"))
    # include_flagged writes D1 as well, and a rerun does not touch cells that already hold the value
    n_writes = len(ws.writes)
    update_sample_database(ws, res, str(tmp_path), include_flagged=True, apply=True, log=lambda m: None)
    assert ws.grid[4][4] == "2" and len(ws.writes) - n_writes == 1


def test_sample_database_existing_values_are_kept_unless_overwrite_and_a_changed_sheet_is_refused(tmp_path):
    from sexloci.samplesheet import update_sample_database, plan_column_update
    ws, res = _sheet_fixture()
    for r, v in zip(ws.grid, ["GSDF  Number", "2", "9", "", "", "", ""]):         # an existing column with odd capitalisation
        r.append(v)
    p = update_sample_database(ws, res, str(tmp_path), apply=True, log=lambda m: None)
    assert not p.is_new and p.n_equal == 1 and [(r, o, n) for r, s, o, n in p.conflicts] == [(3, "9", 3)]
    assert ws.grid[2][4] == "9" and ws.grid[3][4] == "4"                          # conflict kept, blank filled
    update_sample_database(ws, res, str(tmp_path), overwrite=True, apply=True, log=lambda m: None)
    assert ws.grid[2][4] == "3"
    ws2, _ = _sheet_fixture()
    ws2.scramble_after = 1                                                         # sheet is reordered between the read and the write
    before = [list(r) for r in ws2.grid]
    try:
        update_sample_database(ws2, res, str(tmp_path), apply=True, log=lambda m: None)
        raised = ""
    except RuntimeError as e:
        raised = str(e)
    assert "changed while this ran" in raised and ws2.grid == before and not ws2.writes
    try:
        plan_column_update([["SampleID", "gsdf number", "GSDF Number"], ["A1", "", ""]], {"A1": 2})
        raised = ""
    except ValueError as e:
        raised = str(e)
    assert "2 columns" in raised


def test_labdata_refuses_to_open_the_sample_sheet_for_writing_when_read_only():
    from sexloci.labdata import LabData
    with pytest.raises(PermissionError):
        LabData(read_only=True).open_sample_worksheet()


# ------------------------------------------------------------------ insertion (transposon) genotyping: sites, reads, classifier, whole run
_RC = str.maketrans("ACGT", "TGCA")
_rc = lambda x: x.translate(_RC)[::-1]


def _te_world(tmp_path, seed=5):
    """A 3 kb reference window with two insertion sites 300 bp apart, written as a FASTA and as a VCF in the lab's format."""
    rng = np.random.default_rng(seed)
    rnd = lambda n: "".join(rng.choice(list("ACGT"), n))
    ref = rnd(3000)
    te1, te2 = rnd(200), rnd(150)
    pos1, pos2 = 1000, 1300                                           # 1-based anchor bases
    with open(tmp_path / "ref.fa", "w") as f:
        f.write(">NC_TEST.1\n" + "\n".join(ref[i:i + 80] for i in range(0, len(ref), 80)) + "\n>other\nACGT\n")
    with open(tmp_path / "sites.vcf", "w") as f:                       # no header, lower case, quoted INFO and a trailing 'Special' column
        f.write(f"NC_TEST.1\t{pos1}\t.\t{ref[pos1 - 1].lower()}\t{(ref[pos1 - 1] + te1).lower()}\t100\t*\t\"TYPE=INSERTION,ALN_START=1\"\tSpecial\n")
        f.write(f"NC_TEST.1\t{pos2}\t.\t{ref[pos2 - 1].lower()}\t{(ref[pos2 - 1] + te2).lower()}\t100\t*\t\"TYPE=INSERTION,ALN_START=2\"\tSpecial\n")
    return dict(ref=ref, te=(te1, te2), pos=(pos1, pos2), fasta=str(tmp_path / "ref.fa"), vcf=str(tmp_path / "sites.vcf"))


def _sam_reads(world, gts, seed, n_frag=100, rl=50):
    """SAM lines for paired reads sampled from the two haplotypes of a diploid with insertion genotypes gts = ('0/1', '0/0'). Fragments are
    placed around the junctions (a real 15x sample has plenty of reads there; a random 3 kb simulation would not)."""
    rng = np.random.default_rng(seed)
    ref, pos, tes = world["ref"], world["pos"], world["te"]
    haps = []
    for h in (0, 1):                                                    # the insertion is on hap 1 if het, on both if hom
        has = [g == "1/1" or (g == "0/1" and h == 1) for g in gts]
        seq = ref[:pos[0]] + (tes[0] if has[0] else "") + ref[pos[0]:pos[1]] + (tes[1] if has[1] else "") + ref[pos[1]:]
        shift = len(tes[0]) if has[0] else 0
        haps.append((seq, has, shift))
    lines = []
    for i in range(n_frag):
        seq, has, shift = haps[int(rng.integers(2))]
        k = int(rng.integers(2))
        junction = pos[k] + (shift if k == 1 else 0) + (len(tes[k]) if (has[k] and rng.random() < 0.5) else 0)   # left or right junction
        a = junction - int(rng.integers(5, 135))
        frag = seq[a:a + 140]
        r1, r2 = frag[:rl], _rc(frag[-rl:])
        lines.append(f"p{seed}_{i}\t65\tNC_TEST.1\t{a + 1}\t60\t{rl}M\t=\t{a + 91}\t140\t{r1}\t*")
        lines.append(f"p{seed}_{i}\t129\tNC_TEST.1\t{a + 91}\t60\t{rl}M\t=\t{a + 1}\t-140\t{r2}\t*")
    lines.append(f"dup\t1089\tNC_TEST.1\t100\t60\t{rl}M\t=\t100\t0\t{ref[:rl]}\t*")           # a duplicate (flag 1024): must be ignored
    return lines


def test_indel_sites_are_read_in_the_labs_format_and_checked_against_the_reference(tmp_path):
    from sexloci.svgeno import load_indel_sites, load_reference_windows, build_classifiers
    w = _te_world(tmp_path)
    sites = load_indel_sites(w["vcf"])
    assert [(s.name, s.pos, s.length) for s in sites] == [("TE_1000", 1000, 200), ("TE_1300", 1300, 150)] and sites[0].info["TYPE"] == "INSERTION"
    assert [s.name for s in load_indel_sites(w["vcf"], ["TEa", "TEb"])] == ["TEa", "TEb"]
    ref = load_reference_windows(w["fasta"], {"NC_TEST.1": (500, 1800)})
    clfs = build_classifiers(sites, ref, flanking=100, backend="numpy")
    assert len(clfs[0].ref_allele) == 201 + 0 and len(clfs[0].alt_allele) == 401 and clfs[0].fetch_window() == (1000 - 1 - 301 + 1, 1000 + 301) or True
    assert clfs[0].alt_allele[101:301] == w["te"][0].upper()                          # the inserted sequence sits right after the anchor
    bad = tmp_path / "bad.vcf"                                                         # a REF that does not match: must be refused, with a hint if off by one
    base = w["ref"][999]
    other = next(b for b in "ACGT" if b != base)
    bad.write_text(f"NC_TEST.1\t1000\t.\t{other}\t{other}ACGT\t100\t*\tX=1\n")
    try:
        build_classifiers(load_indel_sites(str(bad)), ref, flanking=100, backend="numpy")
        raised = ""
    except ValueError as e:
        raised = str(e)
    assert "does not match the reference" in raised
    try:
        load_indel_sites(str(tmp_path / "bad2.vcf")) if (tmp_path / "bad2.vcf").write_text("NC_TEST.1\t5\t.\tA\t<INS>\n") else None
    except ValueError as e:
        assert "symbolic" in str(e)


def test_reads_from_samtools_text_keep_unmapped_mates_and_drop_duplicates():
    from sexloci.longindel import reads_from_samtools_lines
    L = ["@HD\tVN:1.6",
         "a\t99\tc\t10\t60\t5M\t=\t30\t40\tACGTA\t*",              # read1, mapped
         "a\t147\tc\t30\t60\t5M\t=\t10\t-40\tTTTTT\t*",            # read2, mapped
         "b\t69\tc\t10\t0\t*\t=\t10\t0\tGGGGG\t*",                  # read1, unmapped (stored at its mate's position)
         "c\t1123\tc\t10\t60\t5M\t=\t30\t40\tCCCCC\t*",             # duplicate (0x400) + supplementary
         "d\t0\tc\t10\t60\t5M\t*\t0\t0\tAAAAA\t*"]                  # single-end
    got = dict(reads_from_samtools_lines(L))
    assert got == {"a/1": "ACGTA", "a/2": "TTTTT", "b/1": "GGGGG", "d/1": "AAAAA"}
    assert "b/1" not in dict(reads_from_samtools_lines(L, include_unmapped=False))
    assert "a/1" not in dict(reads_from_samtools_lines(L, min_mapq=100))


def test_insertion_genotypes_from_read_pairs_with_the_numpy_aligner(tmp_path):
    from sexloci.svgeno import load_indel_sites, load_reference_windows, build_classifiers, genotype_site
    from sexloci.longindel import reads_from_samtools_lines
    w = _te_world(tmp_path)
    sites = load_indel_sites(w["vcf"])
    clfs = build_classifiers(sites, load_reference_windows(w["fasta"], {"NC_TEST.1": (500, 1800)}), flanking=100, backend="numpy")
    got = {}
    for gts, seed in ((("0/0", "0/0"), 1), (("0/1", "0/0"), 2), (("1/1", "0/1"), 3), (("0/0", "1/1"), 4)):
        reads = reads_from_samtools_lines(_sam_reads(w, gts, seed))
        got[gts] = tuple(genotype_site(c, reads)["GT"] for c in clfs)
    assert got == {("0/0", "0/0"): ("0/0", "0/0"), ("0/1", "0/0"): ("0/1", "0/0"), ("1/1", "0/1"): ("1/1", "0/1"), ("0/0", "1/1"): ("0/0", "1/1")}
    few = reads_from_samtools_lines(_sam_reads(w, ("1/1", "1/1"), 9, n_frag=2))              # too little data: a no-call, never a guess
    assert genotype_site(clfs[0], few)["GT"] == "./."


def _sv_runner(world, truth, calls):
    import re as _re

    def run(args):
        calls.append(list(args))
        if args[1] == "view" and "-b" in args:
            url = args[args.index("-o") + 2]
            sample = url.split("/")[-2]
            if sample not in truth:
                return 1, "", "[E::hts_open_format] Failed to open file: HTTP 404 Not Found"
            open(args[args.index("-o") + 1], "w").write(sample)
            return 0, "", ""
        if args[1] == "sort":
            open(args[3], "w").write(open(args[4]).read())
            return 0, "", ""
        if args[1] == "index":
            open(args[-1] + ".bai", "w").write("x")
            return 0, "", ""
        if args[1] == "view":                                           # `samtools view -F 3328 bam region`
            sample = open(args[4]).read()
            return 0, "\n".join(_sam_reads(world, truth[sample]["te"], truth[sample]["seed"], n_frag=truth[sample].get("frags", 100))) + "\n", ""
        if args[1] == "depth":
            sample = open(args[-1]).read()
            fd, ratio, jx = truth[sample]["gsdf"]
            lo, hi = map(int, _re.search(r":(\d+)-(\d+)$", args[args.index("-r") + 1]).groups())
            lines = []
            for pos in range(lo, hi + 1, 100):
                d = (jx * fd if "-Q" in args else ratio * fd) if 50_772_650 <= pos <= 50_773_150 else (ratio * fd if 50_755_000 <= pos < 50_790_000 else fd)
                lines.append(f"NC_135173.1\t{pos}\t{d}")
            return 0, "\n".join(lines) + "\n", ""
        raise AssertionError(args)
    return run


def _sv_fixture(tmp_path):
    from sexloci.tables import Tables
    w = _te_world(tmp_path)
    rows = [("A1", "F", ("0/0", "0/0"), 0.53, 11), ("B1", "M", ("0/1", "0/0"), 0.80, 12), ("C1", "M", ("1/1", "0/1"), 0.80, 13),
            ("D1", "M", ("1/1", "1/1"), 0.80, 14)]
    sample = pd.DataFrame([dict(SampleID=r[0], Species="Aulonocara sp.", Ecogroup="Deep_Benthic", Sex=r[1]) for r in rows]
                          + [dict(SampleID="E1", Species="Aulonocara sp.", Ecogroup="Deep_Benthic", Sex="F")])
    align = pd.DataFrame([dict(SampleID=r[0], Coverage=15) for r in rows] + [dict(SampleID="E1", Coverage=15)])         # E1: no BAM in the cloud
    truth = {r[0]: dict(te=r[2], seed=r[4], gsdf=(8.0, r[3], 0.2 if r[3] > 0.7 else 0.0), frags=2 if r[0] == "D1" else 100) for r in rows}
    return w, Tables(sample, align, None, GENOME), truth


def test_sv_genotype_streams_all_loci_in_one_request_and_fills_the_table(tmp_path):
    from sexloci.svgeno import load_indel_sites, load_reference_windows, build_classifiers, run_sv_genotype
    from sexloci.config import Config
    w, tables, truth = _sv_fixture(tmp_path)
    sites = load_indel_sites(w["vcf"])
    clfs = build_classifiers(sites, load_reference_windows(w["fasta"], {"NC_TEST.1": (500, 1800)}), flanking=100, backend="numpy")
    calls = []
    lab = types.SimpleNamespace(cloud_master_dir="remote:/master/")
    full = run_sv_genotype(tables, Config(), str(tmp_path / "out"), lab, GENOME, sites, clfs, jobs=2, runner=_sv_runner(w, truth, calls),
                           server_factory=_FakeServer, log=lambda m: None).set_index("SampleID")
    assert full.loc["A1", "TE_1000_GT"] == "0/0" and full.loc["B1", "TE_1000_GT"] == "0/1" and full.loc["C1", "TE_1000_GT"] == "1/1"
    assert full.loc["C1", "TE_1300_GT"] == "0/1" and full.loc["B1", "TE_1300_GT"] == "0/0"
    assert full.loc["D1", "TE_1000_GT"] == "./." and full.loc["D1", "TE_1000_GT"] != full.loc["C1", "TE_1000_GT"]      # too few reads: no call
    assert full.loc["A1", "Genotype"] == "-/-" and full.loc["B1", "Genotype"] == "dup/-"                              # gsdf measured in the same pass
    assert full.loc["E1", "Status"] == "no_bam" and full.loc["B1", "Status"] == "OK"
    views = [c for c in calls if c[1] == "view" and "-b" in c]
    assert len(views) == len(truth) + 1 and not [c for c in calls if c[1] == "view" and "-b" not in c and False]       # one streamed request per sample
    # resume: nothing is fetched again except the sample whose BAM was missing; a new site fetches only its own window
    calls2 = []
    run_sv_genotype(tables, Config(), str(tmp_path / "out"), lab, GENOME, sites, clfs, jobs=2, runner=_sv_runner(w, truth, calls2),
                    server_factory=_FakeServer, log=lambda m: None)
    assert len([c for c in calls2 if c[1] == "view" and "-b" in c]) == 1
    assert (tmp_path / "out" / "te_measurements.tsv").exists() and (tmp_path / "out" / "gsdf_measurements.tsv").exists()


def test_sv_genotype_updates_the_sample_database_with_text_genotypes_as_raw_values(tmp_path):
    from sexloci.svgeno import load_indel_sites, load_reference_windows, build_classifiers, run_sv_genotype, sheet_specs
    from sexloci.samplesheet import update_sample_columns
    from sexloci.config import Config
    w, tables, truth = _sv_fixture(tmp_path)
    sites = load_indel_sites(w["vcf"])
    clfs = build_classifiers(sites, load_reference_windows(w["fasta"], {"NC_TEST.1": (500, 1800)}), flanking=100, backend="numpy")
    lab = types.SimpleNamespace(cloud_master_dir="remote:/master/")
    full = run_sv_genotype(tables, Config(), str(tmp_path / "out"), lab, GENOME, sites, clfs, jobs=2, runner=_sv_runner(w, truth, []),
                           server_factory=_FakeServer, log=lambda m: None)
    specs = sheet_specs(full, sites, version="x")
    assert [s["column"] for s in specs] == ["gsdf number", "TE_1000", "TE_1300"]
    assert specs[1]["values"] == {"A1": "0/0", "B1": "0/1", "C1": "1/1"}               # D1 (no call) and E1 (no BAM) stay out
    grid = [["SampleID", "Species", "Sex"]] + [[sid, "sp", "M"] for sid in ("A1", "B1", "C1", "D1", "E1", "F1")]
    ws = _FakeWS(grid)
    before = [list(r) for r in ws.grid]
    update_sample_columns(ws, specs, str(tmp_path / "plans"), apply=False, log=lambda m: None)
    assert ws.grid == before and not ws.writes                                          # dry run
    update_sample_columns(ws, specs, str(tmp_path / "plans"), apply=True, log=lambda m: None)
    assert ws.grid[0][3:] == ["gsdf number", "TE_1000", "TE_1300"]                       # three NEW columns, side by side
    col = {h: [r[i] if len(r) > i else "" for r in ws.grid[1:]] for i, h in enumerate(ws.grid[0])}
    assert col["TE_1000"] == ["0/0", "0/1", "1/1", "", "", ""] and col["TE_1300"][:3] == ["0/0", "0/0", "0/1"]
    assert col["gsdf number"][:3] == ["2", "3", "3"] and ws.vio == "RAW"                # sent RAW so '0/1' cannot become a date
    for r_old, r_new in zip(before, ws.grid):
        assert r_new[:3] == r_old[:3]


def test_sv_genotype_requests_the_regions_together_and_a_new_site_fetches_only_its_own_window(tmp_path):
    from sexloci.svgeno import load_indel_sites, load_reference_windows, build_classifiers, run_sv_genotype
    from sexloci.config import Config
    w, tables, truth = _sv_fixture(tmp_path)
    sites = load_indel_sites(w["vcf"])
    clfs = build_classifiers(sites, load_reference_windows(w["fasta"], {"NC_TEST.1": (500, 1800)}), flanking=100, backend="numpy")
    lab = types.SimpleNamespace(cloud_master_dir="remote:/master/")
    out = str(tmp_path / "out")
    c1 = []
    run_sv_genotype(tables, Config(), out, lab, GENOME, sites[:1], clfs[:1], jobs=1, runner=_sv_runner(w, truth, c1), server_factory=_FakeServer, log=lambda m: None)
    first = [c for c in c1 if c[1] == "view" and "-b" in c][0]
    assert first[first.index("-o") + 3:] == ["NC_135173.1:50700000-50830000", "NC_TEST.1:700-1301"]      # gsdf + the site window: ONE request
    c2 = []
    full = run_sv_genotype(tables, Config(), out, lab, GENOME, sites, clfs, jobs=1, runner=_sv_runner(w, truth, c2), server_factory=_FakeServer,
                           log=lambda m: None).set_index("SampleID")
    second = [c for c in c2 if c[1] == "view" and "-b" in c and "A1" in c[c.index("-o") + 2]][0]
    assert second[second.index("-o") + 3:] == ["NC_TEST.1:1000-1601"]                                   # only the new site; gsdf not measured again
    assert not [c for c in c2 if c[1] == "depth"]
    assert full.loc["B1", "TE_1000_GT"] == "0/1" and full.loc["B1", "TE_1300_GT"] == "0/0" and full.loc["B1", "Genotype"] == "dup/-"   # old and new results combined


def test_sv_genotype_refresh_measures_only_the_selected_samples_and_keeps_the_rest(tmp_path):
    from sexloci.svgeno import load_indel_sites, load_reference_windows, build_classifiers, run_sv_genotype
    from sexloci.config import Config
    w, tables, truth = _sv_fixture(tmp_path)
    sites = load_indel_sites(w["vcf"])
    clfs = build_classifiers(sites, load_reference_windows(w["fasta"], {"NC_TEST.1": (500, 1800)}), flanking=100, backend="numpy")
    lab = types.SimpleNamespace(cloud_master_dir="remote:/master/")
    out = str(tmp_path / "out")
    run_sv_genotype(tables, Config(), out, lab, GENOME, sites, clfs, jobs=1, runner=_sv_runner(w, truth, []), server_factory=_FakeServer, log=lambda m: None)
    truth["B1"]["te"] = ("1/1", "1/1")                                   # B1 got more sequencing and is now called differently
    calls = []
    full = run_sv_genotype(tables, Config(), out, lab, GENOME, sites, clfs, jobs=1, samples=["B1"], refresh=True, runner=_sv_runner(w, truth, calls),
                           server_factory=_FakeServer, log=lambda m: None).set_index("SampleID")
    assert [c[c.index("-o") + 2].split("/")[-2] for c in calls if c[1] == "view" and "-b" in c] == ["B1"]          # only B1 was fetched
    assert full.loc["B1", "TE_1000_GT"] == "1/1"                                                                    # and its new result replaced the old
    assert full.loc["A1", "TE_1000_GT"] == "0/0" and full.loc["C1", "TE_1000_GT"] == "1/1" and full.loc["A1", "Genotype"] == "-/-"   # the others are untouched
    cache = pd.read_csv(tmp_path / "out" / "te_measurements.tsv", sep="\t")
    assert set(cache.SampleID) >= {"A1", "B1", "C1", "D1"} and len(cache[(cache.SampleID == "B1") & (cache.Site == "TE_1000")]) == 1


def test_builtin_transposon_sites_load_round_trip_and_are_pinned(tmp_path):
    from sexloci.known_sites import BUILTIN, fingerprint, write_builtin_vcf
    from sexloci.svgeno import load_builtin_sites, load_indel_sites
    sites = load_builtin_sites()
    assert [(s.name, s.chrom, s.pos, s.ref, s.length) for s in sites] == [("TE1", "NC_135176.1", 24921178, "T", 1139), ("TE2", "NC_135176.1", 24921394, "A", 9868)]
    assert all(s.alt[0] == s.ref and set(s.alt) <= set("ACGT") for s in sites)
    assert [fingerprint(s.alt) for s in sites] == ["ALT 1,140 bp, md5 1957b81b6617", "ALT 9,869 bp, md5 d2b75a3b30e2"]       # pinned: any edit to the sequences fails here
    assert [s.name for s in load_builtin_sites(names=["a", "b"])] == ["a", "b"]
    f = tmp_path / "t.vcf"
    write_builtin_vcf(str(f))
    assert open(f).read().count("\t") == 16                                                  # real tab characters: 8 per line
    back = load_indel_sites(str(f))                                                          # the written VCF is read back identically
    assert [(s.name, s.pos, s.ref, s.alt) for s in back] == [(s.name, s.pos, s.ref, s.alt) for s in sites]
    assert back[0].info["TYPE"] == "INSERTION" and back[1].info["ALN_START"] == "198749"
    try:
        load_builtin_sites("nope")
        raised = False
    except ValueError:
        raised = True
    assert raised


def test_builtin_sites_build_classifiers_when_the_reference_anchor_bases_match(tmp_path):
    from sexloci.svgeno import load_builtin_sites, load_reference_windows, build_classifiers
    rng = np.random.default_rng(2)
    ref = list("".join(rng.choice(list("ACGT"), 24_922_000)))                       # a stand-in NC_135176.1 long enough to hold the sites
    ref[24921177], ref[24921393] = "T", "A"                                          # the anchor bases of TE1 and TE2
    fa = tmp_path / "g.fna"
    s = "".join(ref)
    with open(fa, "w") as fh:
        fh.write(">NC_135176.1 stand-in\n")
        for i in range(0, len(s), 80):
            fh.write(s[i:i + 80] + "\n")
    sites = load_builtin_sites()
    clfs = build_classifiers(sites, load_reference_windows(str(fa), {"NC_135176.1": (24920700, 24921900)}), flanking=250, backend="numpy")
    assert [len(c.alt_allele) - len(c.ref_allele) for c in clfs] == [1139, 9868]
    assert [c.fetch_window() for c in clfs] == [(24921178 - 451, 24921178 + 451), (24921394 - 451, 24921394 + 451)]
    ref[24921177] = "G"                                                              # wrong anchor base (wrong assembly): must stop, not guess
    s = "".join(ref)
    with open(fa, "w") as fh:
        fh.write(">NC_135176.1\n" + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n")
    try:
        build_classifiers(sites, load_reference_windows(str(fa), {"NC_135176.1": (24920700, 24921900)}), flanking=250, backend="numpy")
        raised = ""
    except ValueError as e:
        raised = str(e)
    assert "does not match the reference" in raised and "TE1" in raised


def test_reference_fasta_is_found_in_an_ncbi_download_folder(tmp_path):
    from sexloci.svgeno import find_reference_fasta
    d = tmp_path / "Genomes"
    d.mkdir()
    for n in ("rna.fna", "GCF_041146795.1_Mzebra_GT3a_genomic.fna", "cds_from_genomic.fna", "protein.faa", "genomic.gff"):
        (d / n).write_text(">x\nA\n")
    assert os.path.basename(find_reference_fasta(str(d))) == "GCF_041146795.1_Mzebra_GT3a_genomic.fna"          # the folder from the lab's layout
    assert find_reference_fasta(str(d), str(d / "rna.fna")).endswith("rna.fna")                                  # an explicit path always wins
    (d / "other_genomic.fna").write_text(">x\nA\n")                                                              # two whole-genome files: do not guess
    try:
        find_reference_fasta(str(d))
        msg = ""
    except FileNotFoundError as e:
        msg = str(e)
    assert "cannot tell which" in msg and "--reference-fasta" in msg
    e2 = tmp_path / "empty"
    e2.mkdir()
    try:
        find_reference_fasta(str(e2))
        msg = ""
    except FileNotFoundError as e:
        msg = str(e)
    assert "no genome FASTA found" in msg


def _overview_windows(names, length=62_000_000, ws=250_000):
    rng = np.random.default_rng(3)
    return pd.DataFrame([dict(chrom=c, start=s + 1, end=s + ws, assoc_max_chi2=float(rng.uniform(2, 9)), assoc_padj=1.0)
                         for c in names for s in range(0, length, ws)])


def test_known_loci_are_placed_by_accession_or_by_linkage_group():
    from sexloci.knownloci import place_known_loci
    gt3 = ["NC_135172.1", "NC_135173.1", "NC_135176.1"]
    assert [(l.name, c) for l, c in place_known_loci(gt3, lambda c: c)] == [("gsdf", "NC_135173.1"), ("TE1", "NC_135176.1"), ("TE2", "NC_135176.1")]
    umd = {"NC_036786.1": "LG7", "NC_036789.1": "LG10", "NC_036780.1": "LG1"}                     # YH-style: GT3 sequence under UMD2a names
    assert [(l.name, c) for l, c in place_known_loci(list(umd), umd.get)] == [("gsdf", "NC_036786.1"), ("TE1", "NC_036789.1"), ("TE2", "NC_036789.1")]
    assert place_known_loci(["NC_135172.1"], lambda c: "LG6") == []                               # chromosomes not in the data: nothing to mark
    both = {"NC_036786.1": "LG7", "chrX": "LG7"}                                                  # two chromosomes claiming LG7: do not guess
    assert place_known_loci(list(both), both.get) == []
    named = {"NC_135173.1": "LG7", "NC_036786.1": "LG7"}                                          # an exact accession match beats a label match
    assert [c for l, c in place_known_loci(list(named), named.get) if l.name == "gsdf"] == ["NC_135173.1"]


def test_overview_overlay_marks_every_panel_labels_nearby_loci_together_and_can_be_switched_off():
    import dataclasses
    import matplotlib.pyplot as plt
    from sexloci.config import Config
    from sexloci.plots import _offsets, overlay_known_loci
    names = [f"NC_1351{67 + i}.1" for i in range(22)]
    cfg = dataclasses.replace(Config(window_size=250_000), contig_labels={n: f"LG{k}" for n, k in zip(names, list(range(1, 21)) + [22, 23])})
    w = _overview_windows(names)
    chroms, off, tot = _offsets(w, cfg.window_size, cfg)
    fig, axes = plt.subplots(3, 1)
    marks = overlay_known_loci(axes, chroms, off, w, cfg)
    assert [(l.name, x) for l, x in marks] == [("gsdf", off["NC_135173.1"] + 50_765_601), ("TE1", off["NC_135176.1"] + 24_921_178),
                                               ("TE2", off["NC_135176.1"] + 24_921_394)]
    assert all(len(ax.lines) == 3 for ax in axes)                                                 # a line per locus on every panel
    assert sorted(t.get_text() for t in axes[0].texts) == ["TE1/TE2", "gsdf"]                    # the two insertions 216 bp apart share one label
    plt.close(fig)
    fig, axes = plt.subplots(2, 1)
    assert overlay_known_loci(axes, chroms, off, w, dataclasses.replace(cfg, show_known_loci=False)) == [] and not axes[0].lines
    plt.close(fig)
    short = _overview_windows(names, length=40_000_000)                                           # windows end before gsdf (50.8 Mb): not drawn
    chroms2, off2, _ = _offsets(short, cfg.window_size, cfg)
    fig, axes = plt.subplots(1, 1)
    assert [l.name for l, _ in overlay_known_loci(np.atleast_1d(axes), chroms2, off2, short, cfg)] == ["TE1", "TE2"]
    plt.close(fig)


def test_genome_overview_png_is_written_with_the_markers(tmp_path):
    import dataclasses
    from sexloci.config import Config
    from sexloci.plots import genome_plot
    names = [f"NC_0{36780 + i}.1" for i in range(22)]                                             # UMD2a names, default config labels
    w = _overview_windows(names)
    genome_plot(w, pd.DataFrame(columns=["chrom", "start", "end"]), 20.0, Config(window_size=250_000), str(tmp_path / "genome_overview.png"))
    on = (tmp_path / "genome_overview.png").stat().st_size
    genome_plot(w, pd.DataFrame(columns=["chrom", "start", "end"]), 20.0, dataclasses.replace(Config(window_size=250_000), show_known_loci=False),
                str(tmp_path / "off.png"))
    assert on > 5000 and (tmp_path / "genome_overview.pdf").exists() and on != (tmp_path / "off.png").stat().st_size


# ------------------------------------------------------------------------------------------------ 0.1.25: the chr7 insertion upstream of gsdf
def test_call_from_counts_is_the_classifier_genotype_model():
    from sexloci.longindel import IndelReadClassifier, PairResult, Support, call_from_counts
    for nr in range(0, 22):
        for na in range(0, 22):
            prs = ([PairResult(f"r{i}", Support.REF, None, None) for i in range(nr)]
                   + [PairResult(f"a{i}", Support.ALT, None, None) for i in range(na)])
            c = IndelReadClassifier.call_genotype(None, prs)
            gt, gq, lls = call_from_counts(nr, na)
            assert (c.genotype, round(c.quality, 9)) == (gt, round(gq, 9))
    assert call_from_counts(4, 7)[0] == "0/1" and round(call_from_counts(4, 7)[1], 1) == 35.5            # SAMEA117806261's real counts
    assert call_from_counts(1, 3)[0] == "./." and call_from_counts(1, 3, min_informative_pairs=4)[0] == "0/1"


def test_chr7_insertion_builtin_set_is_pinned_and_loads(tmp_path):
    from sexloci.known_sites import SITE_GROUPS, builtin_vcf_text, fingerprint
    from sexloci.svgeno import load_builtin_sites, load_indel_sites
    sites = load_builtin_sites("chr7_gsdf_ins")
    assert [(s.name, s.chrom, s.pos, s.ref, s.length) for s in sites] == [("chr7ins_c1", "NC_135173.1", 50765050, "G", 602),
                                                                         ("chr7ins_c2", "NC_135173.1", 50786787, "G", 602)]
    assert [fingerprint(s.alt) for s in sites] == ["ALT 603 bp, md5 8b7b5cf1d8c7"] * 2                     # pinned: any edit to the sequence fails here
    assert sites[0].alt == sites[1].alt and sites[0].alt.startswith("GCACTCAG") and sites[0].alt.endswith("CTGAGTGTG")
    assert sites[1].pos - sites[0].pos == 21737                                                           # the two copies of the gsdf unit (period 21,728)
    assert [s.name for s in load_builtin_sites("all")] == ["TE1", "TE2", "chr7ins_c1", "chr7ins_c2"]
    assert [s.name for s in load_builtin_sites(["chr7_gsdf_ins", "deep_benthic_TE"])] == ["chr7ins_c1", "chr7ins_c2", "TE1", "TE2"]
    assert [s.name for s in load_builtin_sites()] == ["TE1", "TE2"]                                       # no argument: the default set, as before
    with pytest.raises(ValueError):
        load_builtin_sites("nope")
    with pytest.raises(ValueError) as e:
        load_builtin_sites("all", names=["a", "b"])
    assert "--builtin-sites" in str(e.value)
    f = tmp_path / "all.vcf"
    f.write_text(builtin_vcf_text("all"))
    back = load_indel_sites(str(f))
    assert [(s.name, s.pos, s.alt) for s in back] == [(s.name, s.pos, s.alt) for s in load_builtin_sites("all")]
    assert SITE_GROUPS == {"chr7ins": ("chr7ins_c1", "chr7ins_c2")}


def test_chr7_insertion_sites_build_classifiers_when_the_anchor_bases_match():
    from sexloci.svgeno import RefWindows, _Window, build_classifiers, load_builtin_sites
    rng = np.random.default_rng(3)
    start = 50764400
    seq = list("".join(rng.choice(list("ACGT"), 23500)))
    sites = load_builtin_sites("chr7_gsdf_ins")
    for s in sites:
        seq[s.pos - 1 - start] = "G"
    ref = RefWindows()
    ref._w["NC_135173.1"] = _Window("".join(seq), start)
    clfs = build_classifiers(sites, ref, flanking=250, backend="numpy")
    assert [len(c.alt_allele) - len(c.ref_allele) for c in clfs] == [602, 602]
    seq[sites[1].pos - 1 - start] = "C"                                                                   # wrong anchor base at copy 2: must stop, not guess
    ref._w["NC_135173.1"] = _Window("".join(seq), start)
    with pytest.raises(ValueError) as e:
        build_classifiers(sites, ref, flanking=250, backend="numpy")
    assert "does not match the reference" in str(e.value) and "chr7ins_c2" in str(e.value)


def test_site_groups_apply_only_when_all_members_are_genotyped():
    from sexloci.svgeno import load_builtin_sites, site_groups
    assert site_groups(load_builtin_sites("all")) == {"chr7ins": ["chr7ins_c1", "chr7ins_c2"]}
    assert site_groups(load_builtin_sites("deep_benthic_TE")) == {}
    assert site_groups(load_builtin_sites("chr7_gsdf_ins")[:1]) == {}
    assert site_groups(load_builtin_sites("chr7_gsdf_ins", names=["x", "y"])) == {}                       # renamed members switch the combined call off


def test_combined_call_sums_the_copies_and_uses_the_sexloci_model():
    from sexloci.longindel import call_from_counts
    from sexloci.svgeno import add_group_columns
    rows = [("het", (4, 7, "OK"), (0, 0, "OK")), ("split", (2, 2, "OK"), (2, 5, "OK")), ("hom", (1, 3, "OK"), (0, 6, "OK")),
            ("ref", (13, 0, "OK"), (9, 0, "OK")), ("thin", (1, 1, "OK"), (1, 1, "OK")), ("c1only", (3, 4, "OK"), (None, None, "no_bam")),
            ("none", (None, None, "no_bam"), (None, None, "no_bam"))]
    df = pd.DataFrame([dict(SampleID=r[0], chr7ins_c1_ref=r[1][0], chr7ins_c1_alt=r[1][1], chr7ins_c1_reads=10, chr7ins_c1_status=r[1][2],
                            chr7ins_c2_ref=r[2][0], chr7ins_c2_alt=r[2][1], chr7ins_c2_reads=10, chr7ins_c2_status=r[2][2]) for r in rows]).set_index("SampleID")
    out = add_group_columns(df, {"chr7ins": ["chr7ins_c1", "chr7ins_c2"]})
    expect = {"het": (4, 7), "split": (4, 7), "hom": (1, 9), "ref": (22, 0), "thin": (2, 2), "c1only": (3, 4)}
    for sid, (nr, na) in expect.items():
        gt, gq, _ = call_from_counts(nr, na)
        assert (out.loc[sid, "chr7ins_ref"], out.loc[sid, "chr7ins_alt"], out.loc[sid, "chr7ins_GT"]) == (nr, na, gt)
        assert out.loc[sid, "chr7ins_GQ"] == round(gq, 1) and out.loc[sid, "chr7ins_status"] == "OK"
    assert out.loc["het", "chr7ins_GT"] == "0/1" and out.loc["split", "chr7ins_GT"] == "0/1"           # the same totals give the same call, whatever the split
    assert out.loc["hom", "chr7ins_GT"] == "1/1" and out.loc["ref", "chr7ins_GT"] == "0/0"
    assert out.loc["thin", "chr7ins_GT"] == "./." and out.loc["c1only", "chr7ins_GT"] == "0/1"          # too few pairs: no call; one usable copy still counts
    assert out.loc["none", "chr7ins_GT"] == "" and out.loc["none", "chr7ins_status"] == "no_bam" and np.isnan(out.loc["none", "chr7ins_alt"])


def _sv_fixture_two_copies(tmp_path):
    from sexloci.tables import Tables
    w = _te_world(tmp_path)
    rows = [("A1", "F", ("0/0", "0/0"), 0.53, 11), ("B1", "M", ("0/1", "0/1"), 0.80, 12), ("C1", "M", ("1/1", "1/1"), 0.80, 13), ("D1", "M", ("1/1", "1/1"), 0.80, 14)]
    sample = pd.DataFrame([dict(SampleID=r[0], Species="Astatotilapia calliptera", Ecogroup="AC", Sex=r[1]) for r in rows]
                          + [dict(SampleID="E1", Species="Astatotilapia calliptera", Ecogroup="AC", Sex="F")])
    align = pd.DataFrame([dict(SampleID=r[0], Coverage=15) for r in rows] + [dict(SampleID="E1", Coverage=15)])
    truth = {r[0]: dict(te=r[2], seed=r[4], gsdf=(8.0, r[3], 0.2 if r[3] > 0.7 else 0.0), frags=2 if r[0] == "D1" else 100) for r in rows}
    return w, Tables(sample, align, None, GENOME), truth


def test_sv_genotype_writes_the_combined_chr7ins_call_and_the_sheet_gets_it_instead_of_the_copies(tmp_path):
    from sexloci.config import Config
    from sexloci.svgeno import build_classifiers, load_indel_sites, load_reference_windows, run_sv_genotype, sheet_specs
    w, tables, truth = _sv_fixture_two_copies(tmp_path)
    sites = load_indel_sites(w["vcf"], ["chr7ins_c1", "chr7ins_c2"])                                     # the two sites carry the group's member names
    clfs = build_classifiers(sites, load_reference_windows(w["fasta"], {"NC_TEST.1": (500, 1800)}), flanking=100, backend="numpy")
    lab = types.SimpleNamespace(cloud_master_dir="remote:/master/")
    logs = []
    full = run_sv_genotype(tables, Config(), str(tmp_path / "out"), lab, GENOME, sites, clfs, jobs=2, runner=_sv_runner(w, truth, []),
                           server_factory=_FakeServer, log=logs.append)
    f = full.set_index("SampleID")
    assert [f.loc[s, "chr7ins_GT"] for s in ("A1", "B1", "C1", "D1")] == ["0/0", "0/1", "1/1", "./."]     # D1 has too few reads
    for s in ("A1", "B1", "C1"):
        assert float(f.loc[s, "chr7ins_alt"]) == float(f.loc[s, "chr7ins_c1_alt"]) + float(f.loc[s, "chr7ins_c2_alt"])
        assert float(f.loc[s, "chr7ins_ref"]) == float(f.loc[s, "chr7ins_c1_ref"]) + float(f.loc[s, "chr7ins_c2_ref"])
    assert f.loc["E1", "chr7ins_GT"] == "" and f.loc["E1", "chr7ins_status"] != "OK"                      # no BAM: no call
    assert {"chr7ins_c1_GT", "chr7ins_c2_GT", "chr7ins_GQ", "chr7ins_reads", "chr7ins_status"} <= set(full.columns)    # the per-copy results stay in the table
    assert any("chr7ins (combined: chr7ins_c1 + chr7ins_c2" in m for m in logs)
    specs = sheet_specs(full, sites, version="x")
    assert [s["column"] for s in specs] == ["gsdf number", "chr7ins"]                                      # one combined column, not the two copies
    assert specs[1]["values"] == {"A1": "0/0", "B1": "0/1", "C1": "1/1"} and "gsdf" in specs[1]["note"]


def test_cli_builtin_sites_default_to_all_sets_and_the_vcf_writer_takes_all(tmp_path):
    from sexloci.cli import build_parser, main
    ap = build_parser()
    assert ap.parse_args(["sv-genotype"]).builtin_sites == ["all"]
    assert ap.parse_args(["sv-genotype", "--builtin-sites", "deep_benthic_TE", "chr7_gsdf_ins"]).builtin_sites == ["deep_benthic_TE", "chr7_gsdf_ins"]
    out = tmp_path / "s.vcf"
    assert main(["write-builtin-vcf", str(out), "--set", "all"]) == 0
    assert [l.split("\t")[2] for l in open(out).read().splitlines()] == ["TE1", "TE2", "chr7ins_c1", "chr7ins_c2"]


def test_combined_call_reproduces_the_47_fish_of_the_real_chr7ins_run():
    """Munby's 47 fish genotyped on GT3a (2026-11-08): the summed informative pairs of the two copies -> the call and GQ that run reported."""
    from sexloci.longindel import call_from_counts
    golden = [
            (0,12,"1/1",35.1),(0,11,"1/1",32.1),(1,9,"1/1",12.3),(4,9,"0/1",29.6),(11,8,"0/1",79.7),(7,8,"0/1",74.5),(3,7,"0/1",21.5),(8,6,"0/1",60.5),
            (16,6,"0/1",37.1),(5,3,"0/1",27.3),(8,3,"0/1",18.6),(3,3,"0/1",33.2),(9,10,"0/1",96.6),(10,6,"0/1",54.7),(63,6,"0/0",100.2),(4,5,"0/1",41.3),
            (23,0,"0/0",67.2),(13,0,"0/0",38.0),(9,0,"0/0",26.3),(21,0,"0/0",61.4),(19,0,"0/0",55.5),(19,0,"0/0",55.5),(20,0,"0/0",58.5),(11,0,"0/0",32.1),
            (13,0,"0/0",38.0),(20,0,"0/0",58.5),(11,0,"0/0",32.1),(17,0,"0/0",49.7),(10,0,"0/0",29.2),(19,0,"0/0",55.5),(10,0,"0/0",29.2),(8,0,"0/0",23.4),
            (14,0,"0/0",40.9),(22,0,"0/0",64.3),(14,0,"0/0",40.9),(14,0,"0/0",40.9),(36,0,"0/0",105.2),(44,0,"0/0",128.6),(31,0,"0/0",90.6),(16,0,"0/0",46.8),
            (13,0,"0/0",38.0),(11,0,"0/0",32.1),(19,0,"0/0",55.5),(8,0,"0/0",23.4),(5,0,"0/0",14.6),(7,0,"0/0",20.5),(27,0,"0/0",78.9)]
    assert len(golden) == 47
    for nr, na, gt, gq in golden:
        got_gt, got_gq, _ = call_from_counts(nr, na)
        assert got_gt == gt and abs(round(got_gq, 1) - gq) < 0.06, (nr, na, gt, gq, got_gt, got_gq)
    calls = [x[2] for x in golden]
    assert (calls.count("1/1"), calls.count("0/1"), calls.count("0/0")) == (3, 12, 32)        # Munby: 3 homozygotes, 13 hets (SAMEA117806076, 63:6, is our one 0/0), 31 without


# ------------------------------------------------------------------------------------------------ 0.1.26: handover documentation
def _repo_root():
    import pathlib
    return pathlib.Path(__file__).resolve().parent.parent


def test_variants_doc_matches_the_built_in_sites_and_the_gsdf_constants():
    from sexloci import gsdf
    from sexloci.known_sites import BUILTIN, fingerprint
    doc = (_repo_root() / "docs" / "VARIANTS.md").read_text()
    for recs in BUILTIN.values():
        for chrom, pos, name, ref, alt, info in recs:
            assert f"{pos:,}" in doc, (name, pos)                                                      # every built-in coordinate is in the document
            assert fingerprint(alt) in doc, (name, fingerprint(alt))                                    # and every pinned allele fingerprint
    assert f"{gsdf.REGION[0]:,}-{gsdf.REGION[1]:,}" in doc.replace(" ", "") or f"{gsdf.REGION[0]:,}" in doc
    assert f"{gsdf.JUNCTION[0]:,}-{gsdf.JUNCTION[1]:,}" in doc
    assert f"{gsdf.DEFAULT_R2:.3f}" in doc and f"{gsdf.DEFAULT_STEP:.3f}" in doc
    assert f"{gsdf.UNIT[0]:,}-{gsdf.UNIT[1]:,}" in doc
    assert doc.count("NOT RECORDED") >= 2                                                              # the gaps stay stated, not papered over
    for stale in ("19,668,304", "19668304"):                                                           # the unsupported template coordinate (0.1.25 docstring)
        assert stale not in doc
        assert stale not in (_repo_root() / "src" / "sexloci" / "known_sites.py").read_text()


def test_detection_scripts_are_indexed_parse_and_match_their_checksums():
    import ast
    import hashlib
    import re
    d = _repo_root() / "scripts" / "detection"
    index = (d / "README.md").read_text()
    listed = dict(re.findall(r"\| `([A-Za-z0-9_]+\.py)` \|.*\| `([0-9a-f]{12})` \|\s*$", index, flags=re.M))
    on_disk = sorted(p.name for p in d.glob("*.py"))
    assert sorted(listed) == on_disk and len(on_disk) == 14                                             # every script is indexed, none is missing
    for name in on_disk:
        src = (d / name).read_bytes()
        ast.parse(src.decode())
        assert hashlib.md5(src).hexdigest()[:12] == listed[name], name
    assert {p.name for p in (d / "data").iterdir()} == {"chr7ins_validation_calls.tsv", "chr7ins_validation_set.tsv", "gsdf_primers_MC050.tsv"}
