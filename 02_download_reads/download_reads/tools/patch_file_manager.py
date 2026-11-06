"""Add the download-database methods to helper_modules/file_manager.py WITHOUT replacing the file: two insertions, idempotent, with a diff preview and a backup.

  python -m download_reads patch-file-manager --dry-run        # show the diff, change nothing
  python -m download_reads patch-file-manager                  # apply (a backup file_manager.py.bak-<time> is written)

Inserted: `self.dl_ID` (the DownloadDatabase spreadsheet) after the line that sets `self.a_ID`, and three methods before `_authenticateGS`: downloadManifest(tab, who) -> a Manifest for one tab,
downloadTabs(), readDownloadTab(tab). They only DELEGATE to helper_modules/download_manifest.py (cell-level writes, claims, journal); nothing existing is changed, and the old file managers that are
running are not affected (Python has them loaded already)."""
import argparse, difflib, os, re, sys, time

DL_ID = '1gNFAd_B8tWxQZRJeEmuhVl1qsqJf1KJtKYpehoEMfrQ'
ID_LINE = f"\t\tself.dl_ID = '{DL_ID}' # Download Database (one tab per Run Info file; see 02_download_reads/download_reads/manifest.py)\n"
IMPORT_BLOCK = '''\t\timport sys, glob
\t\there = os.path.dirname(os.path.abspath(__file__))
\t\tfor up in (2, 3):                                              # <repo>/cichlid_sr_sequencing/helper_modules (original layout) or <repo>/common/cichlid_common (after the move)
\t\t\troot = here
\t\t\tfor _ in range(up): root = os.path.dirname(root)
\t\t\tfor d in glob.glob(os.path.join(root, '02_download_reads')):
\t\t\t\tif d not in sys.path: sys.path.insert(0, d)
\t\ttry:
\t\t\tfrom download_reads import manifest as dm
\t\texcept ImportError:
\t\t\tfrom helper_modules import download_manifest as dm        # the original place, until the move
'''
OLD_IMPORT = '''\t\ttry:
\t\t\tfrom helper_modules import download_manifest as dm
\t\texcept ImportError:
\t\t\timport download_manifest as dm
'''
METHODS = '''\t# ---- Download database: one tab per Run Info file. The logic is in 02_download_reads/download_reads/manifest.py: only the cells of changed rows are written (batched), runs are claimed before they
\t# are downloaded, and updates that cannot be written are journalled locally. This is NOT like setSampleDatabase (which rewrites whole worksheets).
\tdef downloadManifest(self, tab, who = None, flush_seconds = 15.0):
''' + IMPORT_BLOCK + '''\t\twho = who or (platform.node().split('.')[0] + '-' + str(os.getpid()))
\t\treturn dm.Manifest(dm.GspreadBackend(self.gc, self.dl_ID), tab, who, journal = self.localTempDir + 'download_manifest_' + tab.replace(' ', '_') + '.journal.jsonl', flush_s = flush_seconds)

\tdef downloadTabs(self):
\t\treturn self.downloadManifest('_').b.tabs()

\tdef readDownloadTab(self, tab):
\t\treturn self.downloadManifest(tab).load()

'''


def patch(text):
    """Returns (new text, list of what was done). Raises ValueError if an anchor is missing."""
    done = []
    if 'self.dl_ID' not in text:
        m = re.search(r'^\t\tself\.a_ID = .*\n', text, re.M)
        if not m:
            raise ValueError("anchor not found: the line that sets self.a_ID")
        text = text[:m.end()] + ID_LINE + text[m.end():]; done.append('dl_ID')
    if 'def downloadManifest' in text and OLD_IMPORT in text:                       # patched by an earlier version: upgrade the import of the manifest module
        text = text.replace(OLD_IMPORT, IMPORT_BLOCK); done.append('methods upgraded')
    if 'def downloadManifest' not in text:
        i = text.find('\tdef _authenticateGS(self):')
        if i < 0:
            raise ValueError("anchor not found: def _authenticateGS")
        text = text[:i] + METHODS + text[i:]; done.append('methods')
    return text, done


def default_file():
    """The FileManager to patch: common/cichlid_common/file_manager_impl.py once it has been moved there, else the original helper_modules/file_manager.py."""
    from cichlid_common.common_light import REPO_ROOT
    import glob
    for pat in ('common/cichlid_common/file_manager_impl.py', 'cichlid_sr_sequencing/helper_modules/file_manager.py', '*/helper_modules/file_manager.py'):
        for f in glob.glob(os.path.join(REPO_ROOT, pat)):
            return f
    return os.path.join(REPO_ROOT, 'cichlid_sr_sequencing', 'helper_modules', 'file_manager.py')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--file', default=default_file()); ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args(argv)
    old = open(a.file).read()
    new, done = patch(old)
    if not done:
        print(f'{a.file}: already patched, nothing to do'); return 0
    print(''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True), 'file_manager.py', 'file_manager.py (patched)', n=1)))
    if a.dry_run:
        print('dry run: nothing written'); return 0
    bak = a.file + '.bak-' + time.strftime('%Y%m%d-%H%M%S')
    open(bak, 'w').write(old); open(a.file, 'w').write(new)
    compile(new, a.file, 'exec')
    print(f'patched ({", ".join(done)}); backup {bak}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
