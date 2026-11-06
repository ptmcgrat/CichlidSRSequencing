"""bchrom: B-chromosome blocks of Lake Malawi cichlids, re-discovered in the coordinates of the project reference (Mzebra_GT3_NCBI).

Re-implements the coverage-ratio analysis of Clark, Conte & Kocher (2018, Genes 9:610) on public Illumina reads (PRJNA552724: 12 B+ females; PRJNA808791: male and female
pools), aligned with the SAME recipe as the mbuna samples (lab pipeline: minimap2 + GATK MergeBamAlignment + MarkDuplicates), so that the blocks, and later the per-sample B
genotype, are defined in one coordinate system with one mapping behaviour. No liftover from the old assembly is needed to define the blocks.

  python -m yhsd.cli bchrom prepare    download -> align -> binned depth, one run at a time, resumable; the reads and BAMs are deleted as soon as the depth track is saved
  python -m yhsd.cli bchrom call       B blocks per carrier, the core set, signal tables, leave-one-out evaluation, control false-positive rates
  python -m yhsd.cli bchrom all        prepare + call
  python -m yhsd.cli bchrom status     what is done, what is running, disk space
  python -m yhsd.cli bchrom clean      delete the depth tracks and the work folder (after `call` has been uploaded)
  python -m yhsd.cli bchrom compare    compare the core set with Clark's blocks lifted to the reference (needs a PAF of his assembly against the reference)
"""
