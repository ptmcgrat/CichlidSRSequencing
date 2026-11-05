"""Compare the BAMs of the NEW pipeline with the BAMs of the OLD pipeline for one sample, record by record.

    python unit_scripts/validate_against_old.py <old sample folder> <new sample folder> <SampleID>

For each of all / unmapped / discordant / inversion / duplication / clipped / chimeric it streams `samtools view` of both files and compares the number of records and an
order-independent fingerprint of every record's alignment fields (QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN; sequence, qualities and tags are not compared). Identical
fingerprints mean the same reads with the same alignments and flags. Copy the old folder first: a new run writes into the same folder name."""
import subprocess, sys, zlib

TYPES = ('all', 'unmapped', 'discordant', 'inversion', 'duplication', 'clipped', 'chimeric')


def fingerprint(lines):
    """(number of records, 64-bit order-independent checksum) of SAM text lines."""
    n, acc = 0, 0
    for line in lines:
        key = '\t'.join(line.split('\t', 9)[:9]).encode()
        acc = (acc + (zlib.crc32(key) | (zlib.crc32(key[::-1]) << 32))) & 0xFFFFFFFFFFFFFFFF
        n += 1
    return n, acc


def bam_fingerprint(path):
    p = subprocess.Popen(['samtools', 'view', path], stdout=subprocess.PIPE, text=True)
    result = fingerprint(p.stdout)
    if p.wait() != 0:
        raise RuntimeError(f'samtools view failed for {path}')
    return result


if __name__ == '__main__':
    old, new, sample = (a.rstrip('/') for a in sys.argv[1:4])
    bad = 0
    print(f'{"BAM":12s} {"old records":>14s} {"new records":>14s}  verdict')
    for t in TYPES:
        name = f'{sample}.{t}.bam'
        a, b = bam_fingerprint(f'{old}/{name}'), bam_fingerprint(f'{new}/{name}')
        verdict = 'IDENTICAL' if a == b else ('same number of records, DIFFERENT content' if a[0] == b[0] else 'DIFFERENT number of records')
        bad += a != b
        print(f'{t:12s} {a[0]:14d} {b[0]:14d}  {verdict}')
    sys.exit(1 if bad else 0)
