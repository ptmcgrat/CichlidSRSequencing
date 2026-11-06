# CichlidSRSequencing

Short-read data of Lake Malawi cichlids, from the read archive to sex-determination results. The goals, the data contracts between the parts and the state of each part are in [`SPEC.md`](SPEC.md); how to get from the
old layout to this one is in [`MIGRATE.md`](MIGRATE.md); decisions and the measurements behind them are in [`docs/decisions.md`](docs/decisions.md).

| Folder | Package | What it does |
|---|---|---|
| `common/` | `cichlid_common` | shared library |
| `01_structural_landmarks/` | `landmarks` | centromeres and inversion breakpoints |
| `02_download_reads/` | `download_reads` | ENA / SRA download, DownloadDatabase sheet |
| `03_align_pipeline/` | (in `cichlid_sr_sequencing/` until moved) | alignment, duplicate marking, splits, slices, upload |
| `04_genotype_sd_systems/` | `sd_genotypes` | gsdf duplication, transposon insertions, inversions, B chromosome |
| `05_variant_calling/` | `variant_calling` | VCFs for families and species |
| `06_sd_region_discovery/` | `sd_regions` | SD regions in many species |
| `07_bc_analysis/` | `bc_analysis` | backcross genotyping and analysis |
| `08_ase_variants/`, `09_candidate_qtn_webserver/` | `ase_variants`, `candidate_qtn` | not started |
| `10_paper_figures/` | `paper_figures` | figures and tables |

```
source tools/env.sh                       # puts every module on PYTHONPATH (or: tools/cichlid <package> <command> ...)
python tools/run_tests.py                 # the tests of every module, each in its own process
python -m download_reads manifest status Blumer Malinsky
python -m sd_genotypes bchrom status
```

Environments: `envs/` (`download.yml` for modules 02 and 04's downloads, `analysis.yml` for the rest). Machines and data locations: SPEC.md, sections 4 and 7.
