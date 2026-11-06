# common: the shared library (`cichlid_common`)

`common_light` (paths, config lookup, subprocess helper, no heavy imports), `common` (the Project: paths built from the FileManager, sample tables, parallel runner), `runner` (immutable run folders on Dropbox with
manifests), `genomes_registry` and `genome` (all genomes of the project), `paf` and `genome_align` (lifting between assemblies), `file_manager` (a shim that finds the lab's FileManager), `cli_base` (the shared options and
dispatch of the module CLIs). `python -m cichlid_common chroms | genomes fetch|annotate|list`. Tests: `common/tests/` (the cross-module integration self-test).
