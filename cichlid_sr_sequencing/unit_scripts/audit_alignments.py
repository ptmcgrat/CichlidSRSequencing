"""Read-only audit of the alignments of ONE genome: AlignmentDatabase rows against what is on Dropbox and on this machine. Changes nothing anywhere. Python 3.7 compatible (utaka's environment).

    python3 audit_alignments.py --genome Mzebra_GT3_NCBI [--out audit.tsv] [--min-share 0.8]

For every sample it finds in the sheet, on Dropbox (Bamfiles/<genome>/<sample>/) or on this machine, it says which of these it is:
  ok              row in the AlignmentDatabase and every usual file on Dropbox
  needs_row       no row, every usual file on Dropbox (the upload finished, the row was not added)
  partial_upload  no row, some files on Dropbox but not all, or an empty file, or a size different from the local copy
  needs_upload    no row, nothing on Dropbox, a complete set of files on this machine
  local_incomplete  no row, nothing on Dropbox, an incomplete set on this machine (the alignment did not finish)
  row_only        a row, but files missing or empty on Dropbox
  not_started     reads in DNAReadsDatabase, but no row, no files here, no files on Dropbox
"The usual files" are not assumed: they are the file-name endings (.all.bam, .all.bam.bai, ...) that at least --min-share of the samples that have a row and files on Dropbox have."""
import argparse, collections, json, os, subprocess, sys


def suffix_of(sample, name):
    return name[len(sample):] if name.startswith(sample) else name


def typical_suffixes(rowed, cloud, min_share=0.8):
    """Endings present for at least min_share of the samples that have a row and some files on Dropbox."""
    base = [s for s in rowed if cloud.get(s)]
    if not base:
        return set()
    n = collections.Counter(suf for s in base for suf in {suffix_of(s, f) for f in cloud[s]})
    return {suf for suf, c in n.items() if c / float(len(base)) >= min_share}


def classify(rowed, cloud, local, reads_samples, typical):
    """-> {sample: (class, details)}"""
    out = {}
    for s in sorted(set(rowed) | set(cloud) | set(local) | set(reads_samples)):
        c = {suffix_of(s, f): sz for f, sz in cloud.get(s, {}).items()}; l = {suffix_of(s, f): sz for f, sz in local.get(s, {}).items()}
        miss_c, miss_l = sorted(typical - set(c)), sorted(typical - set(l))
        empty = sorted(k for k, v in c.items() if v == 0)
        differ = sorted(k for k in c if k in l and c[k] != l[k])
        bad_cloud = bool(miss_c or empty)
        if s in rowed:
            cls = 'ok' if (c and not bad_cloud) else 'row_only'
        elif c:
            cls = 'needs_row' if (not bad_cloud and not differ) else 'partial_upload'
        elif l:
            cls = 'needs_upload' if not miss_l else 'local_incomplete'
        else:
            cls = 'not_started'
        det = []
        if c and miss_c: det.append('missing on Dropbox: ' + ','.join(miss_c))
        if empty: det.append('empty on Dropbox: ' + ','.join(empty))
        if differ: det.append('size differs from the local copy: ' + ','.join(differ))
        if l and miss_l and cls in ('local_incomplete',): det.append('missing locally: ' + ','.join(miss_l))
        out[s] = (cls, '; '.join(det))
    return out


def cloud_listing(cloud_dir, run=subprocess.run):
    """{sample: {file name: size}} of Bamfiles/<genome>/ with ONE recursive rclone lsjson."""
    r = run(['rclone', 'lsjson', '-R', '--files-only', cloud_dir], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if r.returncode == 3:
        return {}
    if r.returncode != 0:
        raise RuntimeError('rclone lsjson %s failed (exit %d): %s' % (cloud_dir, r.returncode, (r.stderr or '')[-200:]))
    out = {}
    for x in json.loads(r.stdout or '[]'):
        parts = x['Path'].split('/')
        if len(parts) >= 2:
            out.setdefault(parts[0], {})['/'.join(parts[1:])] = int(x.get('Size', 0))
    return out


def local_listing(root):
    out = {}
    if not os.path.isdir(root):
        return out
    for s in os.listdir(root):
        d = os.path.join(root, s)
        if os.path.isdir(d):
            for dp, _, fs in os.walk(d):
                for f in fs:
                    out.setdefault(s, {})[os.path.relpath(os.path.join(dp, f), d)] = os.path.getsize(os.path.join(dp, f))
    return out


def main(argv=None, fm=None, cloud_fn=cloud_listing, local_fn=local_listing, out=print):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--genome', required=True); ap.add_argument('--out', default=None); ap.add_argument('--min-share', type=float, default=0.8)
    a = ap.parse_args(argv)
    if fm is None:
        sys.path.insert(0, os.getcwd())
        from helper_modules.file_manager import FileManager
        fm = FileManager()
    fm.setGenome(a.genome); fm.readAlignmentDatabase(); fm.readSampleDatabase()
    rowed = set(str(x) for x in fm.alignment_dt.SampleID.dropna())
    reads_samples = set(str(x) for x in fm.reads_dt.SampleID.dropna())
    cloud = cloud_fn(fm.localBamRefDir.replace(fm.localMasterDir, fm.cloudMasterDir)); local = local_fn(fm.localBamRefDir)
    typical = typical_suffixes(rowed, cloud, a.min_share)
    res = classify(rowed, cloud, local, reads_samples, typical)
    out('genome %s: %d rows in the AlignmentDatabase, %d sample folders on Dropbox, %d on this machine (%s), %d samples with reads' % (a.genome, len(rowed), len(cloud), len(local), fm.localBamRefDir, len(reads_samples)))
    out('usual files: ' + (', '.join(sorted(typical)) or '(none: no sample has a row and files on Dropbox)'))
    by = collections.defaultdict(list)
    for s, (cls, det) in res.items():
        by[cls].append(s)
    for cls in ('ok', 'needs_row', 'partial_upload', 'needs_upload', 'local_incomplete', 'row_only', 'not_started'):
        out('%-17s %5d  %s' % (cls, len(by[cls]), ' '.join(by[cls][:12]) + (' ...' if len(by[cls]) > 12 else '') if cls != 'ok' else ''))
    eco = {}
    if 'Ecogroup' in fm.sample_dt.columns:
        eco = dict(zip(fm.sample_dt.SampleID.astype(str), fm.sample_dt.Ecogroup.astype(str)))
    ns = collections.Counter(eco.get(s, '?') for s in by['not_started'])
    if ns:
        out('not_started by Ecogroup: ' + ', '.join('%s %d' % kv for kv in sorted(ns.items())))
    path = a.out or os.path.join(getattr(fm, 'localTempDir', '.'), 'audit_alignments_%s.tsv' % a.genome)
    with open(path, 'w') as f:
        f.write('SampleID\tclass\tdetails\tEcogroup\tfiles_dropbox\tfiles_local\n')
        for s, (cls, det) in sorted(res.items()):
            f.write('%s\t%s\t%s\t%s\t%d\t%d\n' % (s, cls, det, eco.get(s, ''), len(cloud.get(s, {})), len(local.get(s, {}))))
    out('written: ' + path)
    return res


if __name__ == '__main__':
    main()
