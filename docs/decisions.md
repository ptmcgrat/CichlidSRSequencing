# Decision log

Newest last. Each entry: what was decided, why, and the measurement or fact behind it.

1. **Reads come from NCBI's copy on AWS when ENA is slow** (2026-11-06). ENA moved about 1 Mbit/s over Aspera and 0.3 MB/s over https from three unrelated networks (utaka, D07, a home connection) and four files, at the same time as the AWS copy
   ran at 17-40 MB/s on one connection. Not a Georgia Tech firewall problem; no request to IT. `--source auto` probes Aspera and orders the sources; the other is the fallback.
2. **Every SRA run is checked against ENA.** vdb-validate; fasterq-dump's spot count against ENA's `read_count`; the first 1,000 reads (sequence and quality) against ENA's FASTQ. NCBI also distributes SRA Lite files with simplified qualities.
3. **ENA's `read_count` has no single meaning** (2026-11-06). ERR13389718: 109,512,218 = 2 x 54,756,109 spots (reads). ERR14734229: 40,613,306 against 40,650,256 spots (about pairs). All 8 runs measured have 0.04-0.40% MORE spots in SRA than ENA's count.
   The check tries both readings and accepts up to 1% more, rejects more than 0.1% fewer. ReadLength comes from the Run Info's AvgSpotLen, not from a guess.
4. **Binned qualities are not a sign of SRA Lite.** ERR13389718 (HiSeq 2500) has 6 distinct quality symbols at ENA and at NCBI. Only one or two distinct symbols is treated as Lite-like. No instrument-based rule.
5. **The sheet is written cell by cell, never by whole worksheet** (setSampleDatabase rewrites worksheets and flattens formulas; two machines would overwrite each other). Runs are claimed before work; a claim expires after 6 h; an unreachable sheet
   does not stop a download (journal, replayed later).
6. **`DownloadDatabase` is a spreadsheet of its own**, one tab per Run Info file, so nothing that rewrites the SampleDatabase or the alignment sheets can touch it.
7. **Speeds measured:** D07 11 MB/s (8 connections); mzebra 104 MB/s for one download, about 150 MB/s for three; utaka 40 MB/s (single connection, curl). A run's .sra file is 70% of its FASTQ size; fasterq-dump 38 s; the uBAM conversion 12 min (D07).
   The upload of the uBAMs to Dropbox (about 3.7 TB for the 328 runs) is the likely limit: stage timings are in the sheet since v9d.
8. **The B chromosome is found by replicating Clark et al.'s coverage ratio on the project reference, not by lifting his blocks** (56% of his core sequence is on unplaced contigs of the old assembly).
9. **B-chromosome reads are downloaded through the SRA source** (`bchrom` download.source: auto) with ENA as the fallback; disk is reserved at 9x the ENA FASTQ size for the SRA route (.sra + uncompressed FASTQ + scratch).
10. **One module, one folder, contracts instead of imports** (SPEC.md). A module may import the package of an upstream module it depends on, never its internals.
11. **Module 03 is not moved with the others** because the alignment run on utaka uses the old file manager. The shim `cichlid_common.file_manager` finds the old one; the move is one later change.
12. **The new layout lives on a branch** until no job runs anywhere; the old script names keep working through wrappers for one release.
