"""Analyses of the BC genotypes (module 07), as functions and as `python -m bc_analysis analyze <step>`:
marey (genetic against physical position), rate (is a region's recombination rate unusual), sexscan and finemap (where the sex determiner is, and its limits from the recombinant fish), crossovers (the crossovers
that bound it, localized from the READS), trd (transmission ratio distortion) and linkscan (which part of the genome a segment is inherited with). Every step writes its tables, figures and a provenance.json
under <genome outputs>/analysis/<step>/ and can copy them to Dropbox (--upload)."""
