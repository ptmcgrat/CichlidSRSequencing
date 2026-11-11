# common: the shared library (`cichlid_common`)

`common_light` (paths, config lookup, subprocess helper, no heavy imports), `common` (the Project: paths built from the FileManager, sample tables, parallel runner), `runner` (immutable run folders on Dropbox with
manifests), `genomes_registry` and `genome` (all genomes of the project), `paf` and `genome_align` (lifting between assemblies), `file_manager` (a shim that finds the lab's FileManager), `cli_base` (the shared options and
dispatch of the module CLIs), `gspread_retry` (retry of every Google call: connection resets, 429, 5xx), `data_manifest` (the index of every run folder and file the steps have put on Dropbox: `python -m cichlid_common manifest build | find`). `python -m cichlid_common chroms | genomes fetch|annotate|list`. Tests: `common/tests/` (the cross-module integration self-test).

`cichlid_common/sheet_writer.py`: guarded writing of computed values into columns of an existing Google worksheet (a plan by default, only the named cells, nothing overwritten unless asked, a read-back check and a journal); `fake_sheet.py` is its test double.
