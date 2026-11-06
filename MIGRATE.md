# Moving to the layout of SPEC.md

This tree is the code of `YH_SD_Paper` v34 and the download tool v9d reorganised into the numbered module folders of `SPEC.md`. Nothing in it changes what the code does; the tests of every module pass in the new layout
(`python tools/run_tests.py`), and the one new feature (B-chromosome reads through the SRA source, `04`) has its own tests.

## What moved where

| Before | Now |
|---|---|
| `cichlid_sr_sequencing/YH_SD_Paper/yhsd/` shared files (`common`, `runner`, `genomes_registry`, `genome`), `analyses/inversion_breakpoints/paf.py`, `ensure_paf` | `common/cichlid_common/` |
| `analyses/centromeres`, `analyses/inversion_breakpoints` | `01_structural_landmarks/landmarks/` |
| `helper_modules/ena_fetch, sra_fetch, ena_scheduler, grab_ena_v2, download_manifest`, `downloadShortReadData_v2.py`, `downloadManifest.py`, `unit_scripts/ena_bench.py`, `patch_file_manager.py` | `02_download_reads/download_reads/` (`sources/ena.py`, `sources/sra.py`, `scheduler.py`, `worker.py`, `manifest.py`, `run.py`, `manifest_tool.py`, `bench.py`, `tools/`) |
| `analyses/panel, gsdf, bchrom, sample_table` | `04_genotype_sd_systems/sd_genotypes/` |
| `callset.py` | `05_variant_calling/variant_calling/` |
| `analyses/parental_sex_region, sex_association`, `popgen.py`, `inheritance.py` | `06_sd_region_discovery/sd_regions/` |
| `analyses/bc_genotyping`, `bc_ancestry_hmm.py` | `07_bc_analysis/bc_analysis/` |
| `figures/` | `10_paper_figures/paper_figures/figures/` |
| `selftest*.py` | `<module>/tests/` (the cross-module one in `common/tests/selftest_integration.py`) |
| `config/*.yaml` | `<module>/config/` (`common/config/` for `paper.yaml` and `genomes.yaml`) |

Commands: `python -m yhsd.cli centromeres arrays` becomes `python -m landmarks centromeres arrays`; every option is unchanged. `python -m download_reads run|manifest|bench`; the old script names still work through the wrappers in `02_download_reads/`.
Set the path once with `source tools/env.sh`, or use `tools/cichlid <package> <command>`.

## What is NOT moved: module 03 (alignment)

`alignFastQReads_v2.py`, `seq_pipeline.py`, `bam_split.py`, `bam_slices.py`, `bam_stats.py`, `alignment_db.py`, `disk_ledger.py`, `download_manager.py` and the old `file_manager.py` stay in `cichlid_sr_sequencing/` for now:
an alignment run on utaka is still using the old file manager, which can overwrite another machine's sheet changes. The new code finds them there (`cichlid_common.file_manager` is a shim that imports `helper_modules.file_manager`).
Move them as one change after that run has finished and every machine has pulled: files to `03_align_pipeline/align_pipeline/`, `file_manager.py` to `common/cichlid_common/file_manager_impl.py` (the shim then picks it up), imports
rewritten as in the other modules, tests moved with them. That patch needs the files I have not seen (`alignment_worker.py`, `download_queue.py`, ...), so it is a separate step.

## Order of operations

1. **Now, on the laptop, on a branch** (the machines that are running jobs stay on `YH_SD_Paper`):
   `git checkout -b reorg_spec YH_SD_Paper`, copy this tree over the repo root (`rsync -a --exclude .git <this>/ ~/CichlidSRSequencing/`), `python tools/run_tests.py`, commit, `git push -u origin reorg_spec`.
2. **Try it on one machine** that is not running jobs: `git fetch && git checkout reorg_spec`, `python tools/run_tests.py`, then one real command, e.g. `tools/cichlid download_reads manifest status Blumer Malinsky`.
   The FileManager methods of the sheet are found through the old copies until step 4, so this works before anything is retired.
3. **When no download or alignment job runs anywhere**: merge `reorg_spec` into `YH_SD_Paper`, and on every machine `git pull`.
4. **Retire the old copies** (each is a dry run until `--execute`; it refuses to run while jobs are running or the tree is dirty):
   `tools/migrate_lab_files.sh check`, `retire-old-download`, `upgrade-file-manager`, `retire-old-scripts`, `retire-yhsd`. Commit after each.
5. **After the alignment run on utaka**: the module 03 move described above.

## Things that behave differently

* `config_path(name)` searches `common/config/` and then every `<module>/config/`; config file names must stay unique across modules.
* The worker processes of the downloader are started with `python -m download_reads.worker` and `PYTHONPATH` set from `cichlid_common.common_light.subprocess_env()`; nothing is launched by file path any more.
* `PKG_ROOT` is kept as an alias of `REPO_ROOT`. Run folders record the commit of this repository (the check for tracked files looks at `common/config/paper.yaml`).
* A module may import the package of another module it depends on, never its internals: 04 uses `download_reads.sources` (public), 06 and 07 use `variant_calling.callset`; everything may use `cichlid_common`.
