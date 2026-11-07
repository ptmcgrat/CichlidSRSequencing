"""Add the AlignmentDatabase rows of samples whose BAMs are already uploaded (audit class `needs_row`), WITHOUT re-aligning or re-uploading anything. Python 3.7 compatible (utaka's environment).

It follows the code of alignment_worker.uploadAndUpdateDatabase for one sample: read counts of the seven BAMs (gatk CountReads, as there), coverage = reads * read_length / reference size with read_length = ReadLength / 2
from the DNAReadsDatabase, every other count as a fraction of the total, RunIDs joined with ',,', tool versions from `conda list`, BamSize. The counting runs in parallel; ALL rows are then added with ONE read and ONE write
of the sheet (after a CSV backup, with the grid enlarged if needed and a read-back check).

    python3 register_alignments.py --genome G --from-audit audit.tsv --dry-run --limit 3        # count three samples, write nothing
    python3 register_alignments.py --genome G --check-existing 3                                # recount 3 registered samples and compare with the sheet (validates the counting)
    python3 register_alignments.py --genome G --from-audit audit.tsv --jobs 16                  # the real thing

Needs the BAMs on THIS machine (Bamfiles/<genome>/<sample>/): samples without them are reported and skipped. Only samples without a row are added; nothing already in the sheet is changed."""
import argparse, concurrent.futures, csv, os, subprocess, sys, time

ENDINGS = ['.all.bam', '.unmapped.bam', '.discordant.bam', '.inversion.bam', '.duplication.bam', '.clipped.bam', '.chimeric.bam']
KEYS = ['all', 'unmapped', 'discordant', 'inversion', 'duplication', 'clipped', 'chimeric']


def count_gatk(path, java_options='-Xmx4g'):
    """The count of alignment_worker: line 2 of the stdout of `gatk CountReads -I file`."""
    r = subprocess.run(['gatk', '--java-options', java_options, 'CountReads', '-I', path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    for line in r.stdout.split('\n'):
        if line.strip().isdigit():
            return int(line.strip())
    raise RuntimeError('gatk CountReads gave no count for %s: %s' % (path, (r.stderr or '')[-200:]))


def count_samtools(path, threads=4):
    r = subprocess.run(['samtools', 'view', '-c', '-@', str(threads), path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if r.returncode != 0 or not r.stdout.strip().isdigit():
        raise RuntimeError('samtools view -c failed for %s: %s' % (path, (r.stderr or '')[-200:]))
    return int(r.stdout.strip())


def count_pysam(path, threads=4):
    """samtools view -c with the samtools that is built into pysam (nothing to install; several threads)."""
    import pysam
    out = pysam.view('-c', '-@', str(threads), path)
    if not str(out).strip().isdigit():
        raise RuntimeError('pysam view -c gave no count for %s: %r' % (path, out))
    return int(str(out).strip())


COUNTERS = {'gatk': count_gatk, 'samtools': count_samtools, 'pysam': count_pysam}


def sample_files(bam_dir, sample):
    return {k: os.path.join(bam_dir, sample, sample + e) for k, e in zip(KEYS, ENDINGS)}


def count_sample(bam_dir, sample, counter):
    """{key: reads} of the seven BAMs of one sample; raises if a file is missing."""
    paths = sample_files(bam_dir, sample)
    missing = [k for k, p in paths.items() if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError('%s: missing locally: %s' % (sample, ','.join(missing)))
    return {k: counter(p) for k, p in paths.items()}


def build_row(sample, counts, genome, run_ids, read_length, reference_size, versions, bam_size):
    """The row exactly as alignment_worker builds it."""
    total = counts['all']
    if total <= 0:
        raise ValueError('%s: no reads in the all.bam' % sample)
    frac = {k: (v / float(total) if k != 'all' else v) for k, v in counts.items()}
    row = {'SampleID': sample, 'GenomeVersion': genome, 'RunIDs': ',,'.join(run_ids), 'Coverage': total * read_length / float(reference_size), 'TotalReads': total, 'UnmappedReads': frac['unmapped'],
           'DiscordantReads': frac['discordant'], 'InversionReads': frac['inversion'], 'DuplicationReads': frac['duplication'], 'ClippedReads': frac['clipped'], 'ChimericReads': frac['chimeric']}
    row.update(versions); row['BamSize'] = bam_size
    return row


def tool_versions(run=subprocess.run):
    try:
        out = run(['conda', 'list'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.decode('utf-8').split('\n')
    except Exception:
        out = []
    def v(prefix):
        hit = [x.split()[1] for x in out if x.startswith(prefix) and len(x.split()) > 1]
        return hit[0] if hit else 'unknown'
    return {'minimap2_version': v('minimap2'), 'gatk_version': v('gatk4'), 'pysam_version': v('pysam')}


def _ensure_rows(fm, tab, need, log):
    try:
        ws = fm.gc.open_by_key(fm.a_ID).worksheet(tab)
        if ws.row_count < need:
            old = ws.row_count; ws.resize(rows=need + 50); log('  %s: grid enlarged from %d to %d rows' % (tab, old, need + 50))
    except Exception as e:
        log('  %s: could not check the grid size (%s: %s); continuing' % (tab, type(e).__name__, str(e)[:100]))


def samples_from_audit(path, classes):
    with open(path) as f:
        return [r['SampleID'] for r in csv.DictReader(f, delimiter='\t') if r['class'] in classes]


def main(argv=None, fm=None, counter=None, versions=None, reference_size=None, out=print):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--genome', required=True); ap.add_argument('--from-audit', default=None, help='the TSV of audit_alignments.py'); ap.add_argument('--classes', nargs='+', default=['needs_row'])
    ap.add_argument('--samples', nargs='+', default=None); ap.add_argument('--jobs', type=int, default=8); ap.add_argument('--counter', choices=sorted(COUNTERS), default='gatk', help='gatk = what alignment_worker used (slow: one thread per BAM); samtools / pysam = samtools view -c, many times faster, same count for well-formed BAMs: check it with --check-existing first')
    ap.add_argument('--allow-unlisted', action='store_true', help='also register samples that are not in the SampleDatabase (ReadLength is taken from the DNAReadsDatabase)')
    ap.add_argument('--limit', type=int, default=None); ap.add_argument('--dry-run', action='store_true'); ap.add_argument('--check-existing', type=int, default=0, metavar='N')
    a = ap.parse_args(argv)
    if fm is None:
        sys.path.insert(0, os.getcwd())
        from helper_modules.file_manager import FileManager
        fm = FileManager()
    fm.setGenome(a.genome); fm.readAlignmentDatabase(); fm.readSampleDatabase()
    counter = counter or COUNTERS[a.counter]
    old = fm.alignment_dt.dropna(how='all'); have = set(str(x) for x in old.SampleID.dropna())
    bam_dir = fm.localBamRefDir
    if reference_size is None:
        import pysam
        reference_size = sum(pysam.FastaFile(fm.localGenomeFile).lengths)
    versions = versions or tool_versions()
    pool = concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs)

    def run_counts(samples):
        """Every BAM of every sample is its own task (a few samples do not wait for their slowest file); a line is printed as each sample completes."""
        t0 = time.time(); res, err, left, futs = {s: {} for s in samples}, {}, {}, {}
        for s in samples:
            paths = sample_files(bam_dir, s); miss = [k for k, p in paths.items() if not os.path.exists(p)]
            if miss:
                err[s] = '%s: missing locally: %s' % (s, ','.join(miss)); continue
            left[s] = len(KEYS)
            for k in KEYS:
                futs[pool.submit(counter, paths[k])] = (s, k)
        done, total = 0, len(left)
        for f in concurrent.futures.as_completed(futs):
            s, k = futs[f]
            try:
                res[s][k] = f.result()
            except Exception as e:
                err.setdefault(s, '%s: %s' % (s, str(e)[:200]))
            left[s] -= 1
            if left[s] == 0:
                done += 1; sec = int(time.time() - t0); out('  [%dm%02ds] %d/%d samples counted (%s)' % (sec // 60, sec % 60, done, total, s))
        return {s: r for s, r in res.items() if s not in err and len(r) == len(KEYS)}, err

    if a.check_existing:                                                                         # validate the counting against rows that are already in the sheet
        cand = [s for s in sorted(have) if all(os.path.exists(p) for p in sample_files(bam_dir, s).values())][:a.check_existing]
        out('checking %d registered samples that have their BAMs here: %s' % (len(cand), ' '.join(cand)))
        res, err = run_counts(cand); bad = 0
        for s in cand:
            if s in err:
                out('  %s: %s' % (s, err[s])); bad += 1; continue
            stored = old[old.SampleID.astype(str) == s].iloc[0]
            rd = fm.merged_dt[fm.merged_dt['SampleID'] == s]['ReadLength'].values
            if not len(rd):
                rd = fm.reads_dt[fm.reads_dt['SampleID'] == s]['ReadLength'].values
            run_ids = list(fm.reads_dt[fm.reads_dt['SampleID'] == s].RunID)
            row = build_row(s, res[s], a.genome, run_ids, (rd[0] / 2) if len(rd) else float('nan'), reference_size, versions, os.path.getsize(sample_files(bam_dir, s)['all']))
            diffs, n = [], 0
            for col in ['TotalReads', 'UnmappedReads', 'DiscordantReads', 'InversionReads', 'DuplicationReads', 'ClippedReads', 'ChimericReads', 'Coverage', 'RunIDs']:
                if col not in stored.index or row[col] != row[col]:
                    continue                                                                           # no such column, or no ReadLength to compute the coverage from
                n += 1
                try:
                    same = (str(stored[col]) == str(row[col])) if col == 'RunIDs' else abs(float(stored[col]) - float(row[col])) <= 1e-6 * max(1.0, abs(float(row[col])))
                except (TypeError, ValueError):
                    same = False
                if not same:
                    diffs.append('%s: sheet %s, computed %s' % (col, stored[col], row[col]))
            out('  %s: %d values compared with the sheet (reads %d, coverage %.2f) -> %s' % (s, n, res[s]['all'], row['Coverage'], 'same' if not diffs else 'DIFFERENT: ' + '; '.join(diffs)))
            bad += 1 if diffs else 0
        out('the counting %s the sheet' % ('MATCHES' if not bad and cand else 'DOES NOT MATCH (or could not be checked against)'))
        return bad == 0 and bool(cand)
    todo = a.samples or (samples_from_audit(a.from_audit, a.classes) if a.from_audit else [])
    if not todo:
        raise SystemExit('give --from-audit <tsv> or --samples ...')
    already = [s for s in todo if s in have]; todo = [s for s in todo if s not in have]
    if a.limit:
        todo = todo[:a.limit]
    out('%d samples to register (%d already have a row and are left alone)' % (len(todo), len(already)))
    t0 = time.time(); counts, errors = run_counts(todo)
    out('counted %d samples in %.0f s; %d failed' % (len(counts), time.time() - t0, len(errors)))
    for s, e in sorted(errors.items())[:20]:
        out('  not registered: %s' % e)
    rows = []
    for s in todo:
        if s not in counts:
            continue
        rd = fm.merged_dt[fm.merged_dt['SampleID'] == s]['ReadLength'].values                          # as alignment_worker: SampleDatabase x DNAReadsDatabase
        if not len(rd):
            own = fm.reads_dt[fm.reads_dt['SampleID'] == s]['ReadLength'].values                       # the DNAReadsDatabase has the ReadLength itself
            if not len(own):
                out('  not registered: %s is not in the DNAReadsDatabase' % s); continue
            if not a.allow_unlisted:
                out('  not registered: %s has reads in the DNAReadsDatabase but no row in the SampleDatabase (alignment_worker skips such samples too; --allow-unlisted registers them)' % s); continue
            rd = own
        run_ids = list(fm.reads_dt[fm.reads_dt['SampleID'] == s].RunID)
        try:
            rows.append(build_row(s, counts[s], a.genome, run_ids, rd[0] / 2, reference_size, versions, os.path.getsize(sample_files(bam_dir, s)['all'])))
        except (ValueError, OSError) as e:
            out('  not registered: %s' % e)
    for r in rows[:5]:
        out('  %s: coverage %.1f, %d reads, %.3f unmapped, runs %s' % (r['SampleID'], r['Coverage'], r['TotalReads'], r['UnmappedReads'], r['RunIDs']))
    if a.dry_run or not rows:
        out('dry run: nothing written' if a.dry_run else 'nothing to add'); return rows
    cols = set(old.columns)
    if not set(rows[0]) <= cols:
        raise SystemExit('the sheet has no column for %s; nothing written' % sorted(set(rows[0]) - cols))
    import pandas as pd
    stamp = time.strftime('%Y%m%d-%H%M%S'); bdir = '%sdatabase_backup_%s/' % (fm.localTempDir, stamp); os.makedirs(bdir, exist_ok=True)
    fm.alignment_dt.to_csv(bdir + 'AlignmentDatabase_%s.csv' % a.genome, index=False)
    out('  the sheet was saved to %s (Google Sheets also keeps its version history)' % bdir)
    fm.alignment_dt = pd.concat([old, pd.DataFrame(rows)], ignore_index=True)
    assert fm.alignment_dt.SampleID.is_unique
    _ensure_rows(fm, a.genome, len(fm.alignment_dt) + 1, out)
    if not fm.setAlignmentDatabase():
        raise SystemExit('the write failed after retries; the sheet may be partly written: restore from %s, then run again' % bdir)
    want = list(fm.alignment_dt.SampleID.astype(str)); fm.readAlignmentDatabase()
    if list(fm.alignment_dt.dropna(how='all').SampleID.astype(str)) != want:
        raise SystemExit('the sheet does not hold what was written (expected %d rows, found %d): see %s' % (len(want), len(fm.alignment_dt), bdir))
    out('  verified: AlignmentDatabase %s %d -> %d rows' % (a.genome, len(old), len(fm.alignment_dt)))
    return rows


if __name__ == '__main__':
    main()
