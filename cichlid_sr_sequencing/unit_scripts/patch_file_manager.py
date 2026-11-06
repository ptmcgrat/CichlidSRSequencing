"""Add the download-database methods to helper_modules/file_manager.py WITHOUT replacing the file: two insertions, idempotent, with a diff preview and a backup.

  python3 unit_scripts/patch_file_manager.py --dry-run        # show the diff, change nothing
  python3 unit_scripts/patch_file_manager.py                  # apply (a backup file_manager.py.bak-<time> is written)

Inserted: `self.dl_ID` (the DownloadDatabase spreadsheet) after the line that sets `self.a_ID`, and three methods before `_authenticateGS`: downloadManifest(tab, who) -> a Manifest for one tab,
downloadTabs(), readDownloadTab(tab). They only DELEGATE to helper_modules/download_manifest.py (cell-level writes, claims, journal); nothing existing is changed, and the old file managers that are
running are not affected (Python has them loaded already)."""
import argparse, difflib, os, re, sys, time

DL_ID = '1gNFAd_B8tWxQZRJeEmuhVl1qsqJf1KJtKYpehoEMfrQ'
ID_LINE = f"\t\tself.dl_ID = '{DL_ID}' # Download Database (one tab per Run Info file; see helper_modules/download_manifest.py)\n"
METHODS = '''\t# ---- Download database: one tab per Run Info file. The logic is in helper_modules/download_manifest.py: only the cells of changed rows are written (batched), runs are claimed before they
\t# are downloaded, and updates that cannot be written are journalled locally. This is NOT like setSampleDatabase (which rewrites whole worksheets).
\tdef downloadManifest(self, tab, who = None, flush_seconds = 15.0):
\t\ttry:
\t\t\tfrom helper_modules import download_manifest as dm
\t\texcept ImportError:
\t\t\timport download_manifest as dm
\t\twho = who or (platform.node().split('.')[0] + '-' + str(os.getpid()))
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
    if 'def downloadManifest' not in text:
        i = text.find('\tdef _authenticateGS(self):')
        if i < 0:
            raise ValueError("anchor not found: def _authenticateGS")
        text = text[:i] + METHODS + text[i:]; done.append('methods')
    return text, done


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--file', default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'helper_modules', 'file_manager.py')); ap.add_argument('--dry-run', action='store_true')
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
