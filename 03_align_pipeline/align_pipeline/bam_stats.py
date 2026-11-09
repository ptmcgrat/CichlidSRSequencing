"""Alignment statistics for the sequencing pipeline, computed from tool output that is fast to produce and unambiguous.

Why this module exists (what the original code did):
  * `UnmappedReads` was the number of records in unmapped.bam divided by ALL records. unmapped.bam only holds unmapped reads whose MATE is mapped (orphans), not duplicates, mean
    quality > 25: read pairs with both reads unmapped are never in it. It is NOT the unmapped fraction.
  * `DuplicationReads` was the fraction of records in duplication.bam, i.e. read pairs that are evidence of a tandem DUPLICATION (an orientation signal). It is NOT the PCR/optical
    duplicate rate. The duplicate rate was never recorded: the Picard metrics file was written and never read.
  * every split-BAM fraction divided by a record count that includes duplicates, while the split BAMs exclude duplicates.
The six split-BAM fractions (UnmappedReads, DiscordantReads, InversionReads, DuplicationReads, ClippedReads, ChimericReads) are no longer computed: the ledger drops them and only the chimeric split is
written (decisions 37 and 41). TotalReads and Coverage keep their original meaning, and new, explicitly named columns are added:
  PrimaryReads, MappedFraction, UnmappedFraction, UnmappedNoPositionFraction, UnmappedOrphanFraction, DuplicateFraction (Picard definition), DuplicateFractionFromBam,
  OpticalDuplicateFraction, MeanReadLength, MappedBases, CoverageMapped, CoverageNonDuplicate, ProperlyPairedFraction, MQ0Fraction, InsertSizeMean, InsertSizeSD, ErrorRate, QCWarnings.
Everything is parsed from `samtools stats`, `samtools idxstats` and the Picard MarkDuplicates metrics file; the parsers are pure functions that are unit tested."""
import math
import re

LEGACY_BAM_TYPES = ('all',)
LEGACY_COLUMNS = ('TotalReads', 'Coverage')
NEW_COLUMNS = ('PrimaryReads', 'MappedFraction', 'UnmappedFraction', 'UnmappedNoPositionFraction', 'UnmappedOrphanFraction', 'DuplicateFraction', 'DuplicateFractionFromBam',
               'OpticalDuplicateFraction', 'MeanReadLength', 'MappedBases', 'CoverageMapped', 'CoverageNonDuplicate', 'ProperlyPairedFraction', 'MQ0Fraction', 'InsertSizeMean',
               'InsertSizeSD', 'ErrorRate', 'UBamReads', 'QCWarnings', 'PipelineVersion')
NAN = float('nan')


def _num(raw):
    raw = raw.strip()
    return int(raw) if re.fullmatch(r'-?\d+', raw) else float(raw)


def parse_samtools_stats(text):
    """The `SN` summary lines of `samtools stats` -> {'raw total sequences': 123, 'error rate': 0.002, ...}. 'raw total sequences' excludes secondary and supplementary alignments."""
    sn = {}
    for line in text.splitlines():
        if not line.startswith('SN\t'):
            continue
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        try:
            sn[parts[1].strip().rstrip(':').strip()] = _num(parts[2])
        except ValueError:
            continue
    return sn


def parse_idxstats(text):
    """`samtools idxstats` -> list of dicts (name, length, mapped, unmapped). The '*' row holds reads without a position (both mates unmapped)."""
    rows = []
    for line in text.splitlines():
        p = line.split('\t')
        if len(p) >= 4:
            rows.append(dict(name=p[0], length=int(p[1]), mapped=int(p[2]), unmapped=int(p[3])))
    return rows


def idxstats_totals(rows):
    """records = every record in the file (the number `gatk CountReads` reported); unplaced = reads with no position at all."""
    unplaced = sum(r['unmapped'] + r['mapped'] for r in rows if r['name'] == '*')
    return dict(records=sum(r['mapped'] + r['unmapped'] for r in rows), unplaced_unmapped=unplaced, mapped=sum(r['mapped'] for r in rows if r['name'] != '*'))


def parse_picard_metrics(text):
    """The DuplicationMetrics table of a Picard/GATK MarkDuplicates metrics file -> one dict per LIBRARY (numbers converted; empty cells None)."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith('## METRICS CLASS'):
            header = lines[i + 1].split('\t')
            rows, j = [], i + 2
            while j < len(lines) and lines[j].strip() != '':
                cells = lines[j].split('\t')
                row = {}
                for k, h in enumerate(header):
                    v = cells[k] if k < len(cells) else ''
                    try:
                        row[h] = _num(v) if v.strip() != '' else None
                    except ValueError:
                        row[h] = v
                rows.append(row)
                j += 1
            return rows
    return []


def duplication_summary(libraries):
    """Picard's definition, summed over libraries: PERCENT_DUPLICATION = (UNPAIRED_READ_DUPLICATES + 2 * READ_PAIR_DUPLICATES) / (UNPAIRED_READS_EXAMINED + 2 * READ_PAIRS_EXAMINED).
    Only mapped reads are examined; unmapped reads are in neither the numerator nor the denominator. `unmapped_reads` is Picard's UNMAPPED_READS: EVERY unmapped primary read, including those whose mate is
    mapped (they sit at the mate's position), so it equals samtools' `reads unmapped`, not the number of reads without a position (the `*` line of idxstats: both mates unmapped). Measured 2026-11-09: 1,896,658 in both,
    against 582,860 reads without a position."""
    g = lambda r, k: r.get(k) or 0
    examined = sum(g(r, 'UNPAIRED_READS_EXAMINED') + 2 * g(r, 'READ_PAIRS_EXAMINED') for r in libraries)
    dups = sum(g(r, 'UNPAIRED_READ_DUPLICATES') + 2 * g(r, 'READ_PAIR_DUPLICATES') for r in libraries)
    optical = sum(2 * g(r, 'READ_PAIR_OPTICAL_DUPLICATES') for r in libraries)
    return dict(examined_reads=examined, duplicate_reads=dups, optical_reads=optical, fraction=(dups / examined) if examined else NAN,
                optical_fraction=(optical / examined) if examined else NAN, unmapped_reads=sum(g(r, 'UNMAPPED_READS') for r in libraries), n_libraries=len(libraries))


def _ratio(a, b):
    return (a / b) if b else NAN


def legacy_columns(record_counts, read_length=None, reference_size=None):
    """The two columns the database always had that are kept, with their ORIGINAL meaning. record_counts: {'all': n} = number of records in all.bam (idxstats), the number `gatk CountReads` gave."""
    total = record_counts['all']
    return {'TotalReads': total, 'Coverage': (total * (read_length / 2) / reference_size) if (read_length and reference_size) else NAN}      # the original formula, kept as it was


def derive_columns(sn, idx_all, picard_libraries, reference_size, ubam_reads=None, dup_tolerance=0.01, unmapped_tolerance=0.005, pipeline_version='align_pipeline'):
    """The new, explicitly named columns plus a list of warnings when independent measurements disagree."""
    primary = sn.get('raw total sequences', 0)
    mapped, unmapped = sn.get('reads mapped', 0), sn.get('reads unmapped', 0)
    dup = duplication_summary(picard_libraries) if picard_libraries else None
    dup_bam = _ratio(sn.get('reads duplicated', 0), mapped)
    cols = dict(PrimaryReads=primary, MappedFraction=_ratio(mapped, primary), UnmappedFraction=_ratio(unmapped, primary),
                UnmappedNoPositionFraction=_ratio(idx_all['unplaced_unmapped'], primary), UnmappedOrphanFraction=_ratio(unmapped - idx_all['unplaced_unmapped'], primary),
                DuplicateFraction=dup['fraction'] if dup else dup_bam, DuplicateFractionFromBam=dup_bam, OpticalDuplicateFraction=dup['optical_fraction'] if dup else NAN,
                MeanReadLength=sn.get('average length', NAN), MappedBases=sn.get('bases mapped (cigar)', 0),
                CoverageMapped=_ratio(sn.get('bases mapped (cigar)', 0), reference_size),
                CoverageNonDuplicate=_ratio(sn.get('bases mapped (cigar)', 0) - sn.get('bases duplicated', 0), reference_size),
                ProperlyPairedFraction=_ratio(sn.get('reads properly paired', 0), primary), MQ0Fraction=_ratio(sn.get('reads MQ0', 0), primary),
                InsertSizeMean=sn.get('insert size average', NAN), InsertSizeSD=sn.get('insert size standard deviation', NAN), ErrorRate=sn.get('error rate', NAN),
                UBamReads=ubam_reads if ubam_reads is not None else NAN)
    warnings = []
    if primary == 0:
        warnings.append('no primary reads counted')
    if mapped + unmapped != primary:
        warnings.append(f'mapped ({mapped}) + unmapped ({unmapped}) != primary reads ({primary})')
    if dup is None:
        warnings.append('no Picard duplication metrics: DuplicateFraction is the flag-based value')
    elif not math.isnan(dup_bam) and not math.isnan(dup['fraction']) and abs(dup['fraction'] - dup_bam) > dup_tolerance:
        warnings.append(f"duplicate rate differs: Picard {dup['fraction']:.4f} vs flags in the BAM {dup_bam:.4f}")
    if dup is not None and primary and abs(dup['unmapped_reads'] - unmapped) > unmapped_tolerance * primary:             # both count every unmapped primary read (see duplication_summary)
        warnings.append(f"unmapped reads differ: Picard {dup['unmapped_reads']} vs samtools stats {unmapped} (check)")
    if ubam_reads is not None and ubam_reads != primary:
        warnings.append(f'reads in the unmapped BAM(s) ({ubam_reads}) != primary reads in the final BAM ({primary})')
    cols['QCWarnings'] = '; '.join(warnings)
    cols['PipelineVersion'] = pipeline_version
    return cols
