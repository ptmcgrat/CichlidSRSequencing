"""Single-pass replacement for unit_scripts/split_bamfile_by_contig.py.

The old approach started one python process per (sample, contig), each re-opening the BAM, and then ran `gatk MergeSamFiles` six times per sample to stitch ~64 small BAMs per type
back together. This version reads the coordinate-sorted BAM ONCE, writes the six outputs directly (already sorted, no merge step), and decides every read with `route()`, which is a
line-by-line transcription of the old decision logic (see tests/test_bam_split.py: it is compared with the old logic on tens of thousands of randomised reads). Differences, on purpose:
  * a read without base qualities is skipped by the quality-dependent rules instead of crashing the whole contig job;
  * the mean base quality is only computed when a rule needs it (the old script summed it for every read, duplicates and mapped reads included).
Unchanged on purpose: reads with no position (both mates unmapped) are NOT written anywhere here, exactly as before (the old per-contig fetch never returned them); they stay in all.bam
and are counted separately in the statistics. Single-end reads would all be routed to `discordant` by the old rules (mate_is_unmapped is False for unpaired reads); that is kept."""
import os

OUTPUTS = ('unmapped', 'discordant', 'inversion', 'duplication', 'clipped', 'chimeric')
MIN_MEAN_QUALITY = 25
MIN_CLIP = 5
MIN_INSERT_FOR_DUPLICATION = 102


def route(read):
    """Which outputs a read belongs to: a tuple (possibly empty) of names from OUTPUTS. `read` needs: flag, query_qualities, has_tag(), cigartuples, mapping_quality,
    reference_id, next_reference_id, reference_start, next_reference_start, template_length (all pysam AlignedSegment attributes)."""
    flag = read.flag
    if flag & 0x400:                                         # duplicate
        return ()
    quals = read.query_qualities
    if flag & 0x4:                                           # unmapped (placed at its mate's position)
        if quals is not None and len(quals) and sum(quals) / len(quals) > MIN_MEAN_QUALITY:
            return ('unmapped',)
        return ()
    out = []
    if read.has_tag('SA'):                                   # chimeric / supplementary evidence
        out.append('chimeric')
    else:
        cig = read.cigartuples
        if cig is not None and sum(1 for op, _ in cig if op == 4) == 1:          # exactly one soft clip
            if not read.has_tag('XM') and not (flag & 0x100) and read.mapping_quality != 0:
                if cig[0][0] == 4 and cig[0][1] > MIN_CLIP:
                    n = cig[0][1]
                    if quals is not None and sum(quals[:n]) / n > MIN_MEAN_QUALITY:
                        out.append('clipped')
                elif cig[-1][0] == 4 and cig[-1][1] > MIN_CLIP:
                    n = cig[-1][1]
                    if quals is not None and sum(quals[-n:]) / n > MIN_MEAN_QUALITY:
                        out.append('clipped')
    if not (flag & 0x8):                                     # mate is mapped
        if read.reference_id != read.next_reference_id:      # chromosome fusion
            out.append('discordant')
        elif bool(flag & 0x10) == bool(flag & 0x20):         # both reads on the same strand: inversion
            out.append('inversion')
        else:
            pos, mpos = read.reference_start, read.next_reference_start
            if ((pos < mpos and flag & 0x10) or (pos > mpos and flag & 0x20)) and abs(read.template_length) > MIN_INSERT_FOR_DUPLICATION:
                out.append('duplication')                    # outward-facing pair, large insert: tandem duplication
    return tuple(out)


def split_bam(bam_path, prefix, threads=2):
    """Read `bam_path` once and write {prefix}.{type}.bam for the six types (written as .partial and renamed when complete). Returns a dict of counts for the QC summary."""
    import pysam
    counts = {'reads_seen': 0, 'duplicates_skipped': 0, 'no_position_skipped': 0, 'no_quality': 0}
    counts.update({k: 0 for k in OUTPUTS})
    reader = pysam.AlignmentFile(bam_path, threads=threads)
    writers, partial = {}, {}
    try:
        for name in OUTPUTS:
            partial[name] = f'{prefix}.{name}.bam.partial'
            writers[name] = pysam.AlignmentFile(partial[name], mode='wb', template=reader, threads=1)
        for read in reader:
            counts['reads_seen'] += 1
            if read.reference_id < 0:                        # no position (both mates unmapped): never written, as before
                counts['no_position_skipped'] += 1
                continue
            if read.flag & 0x400:
                counts['duplicates_skipped'] += 1
                continue
            if read.query_qualities is None:
                counts['no_quality'] += 1
            for name in route(read):
                writers[name].write(read)
                counts[name] += 1
    finally:
        reader.close()
        for w in writers.values():
            w.close()
    for name in OUTPUTS:
        os.replace(partial[name], f'{prefix}.{name}.bam')
    return counts
