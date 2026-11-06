"""The download manifest: one Google-Sheets TAB per Run Info file (Blumer, Malinsky, ...), one ROW per run, with the run's metadata plus what happened to it: status, when, from which site, by which
machine, how many attempts, why it failed, integrity, size, speed.  Lives in its OWN spreadsheet (FileManager.dl_ID), so nothing that rewrites the SampleDatabase / alignment sheets can touch it.

How it writes, and why NOT like FileManager.setSampleDatabase (which rewrites whole worksheets from a DataFrame and flattens formulas, and would let two machines overwrite each other):
  * only the CELLS of the rows that changed are written (one batched call), and the row of a run is looked up by its RunID at every flush, so sorting or inserting rows in the sheet is harmless;
  * a run is CLAIMED (ClaimedBy / ClaimedAt, re-read to confirm) before it is downloaded, so two machines do not take the same run; a claim expires after a few hours;
  * when the sheet cannot be reached, the updates go to a local journal and are replayed at the next flush: a Google outage never stops a download.
The storage is behind a small backend interface (GspreadBackend for Google, MemoryBackend for the tests). GspreadBackend is the only part that talks to Google: it uses worksheet.get_all_values /
col_values / row_values / batch_update / append_rows / add_worksheet / add_cols, which have the same signature in gspread 5 and 6."""
import calendar, json, os, random, time
import pandas as pd

RUNINFO_COLUMNS = ['Run', 'AvgSpotLen', 'Bases', 'BioProject', 'BioSample', 'Experiment', 'Instrument', 'Library Name', 'LibraryLayout', 'LibrarySelection', 'LibrarySource', 'Organism', 'Platform', 'SRA Study']
TRACK_COLUMNS = ['DownloadStatus', 'InDatabase', 'DownloadedOn', 'Source', 'Machine', 'ClaimedBy', 'ClaimedAt', 'UpdatedAt', 'Attempts', 'FailureKind', 'LastError', 'Verified', 'QualityCheck', 'ReadCountCheck', 'ReadsSource',
                 'ReadsInUBAM', 'UBamBytes', 'UBamPath', 'DownloadSeconds', 'MBPerSec', 'ConvertSeconds', 'UploadSeconds', 'PipelineVersion', 'Notes', 'Skip']
STATUSES = ['pending', 'claimed', 'downloading', 'converting', 'uploaded', 'in_database', 'failed', 'skipped']
ACTIVE = ('claimed', 'downloading', 'converting')
DONE = ('uploaded', 'in_database', 'skipped')
RETRYABLE = ('', 'transient', 'throttle', 'source_unavailable')           # every source was tried and none delivered: worth another try after a fix or when a source changes


def col_letters(n):
    s = ''
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def a1(r, c):
    return f'{col_letters(c)}{r}'


def utc_iso(t=None):
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(t if t is not None else time.time()))


def parse_iso(s):
    try:
        return calendar.timegm(time.strptime(s, '%Y-%m-%dT%H:%M:%SZ'))
    except (ValueError, TypeError):
        return None


def cell(v):
    if v is None:
        return ''
    if isinstance(v, bool):
        return 'TRUE' if v else 'FALSE'
    if isinstance(v, float):
        return f'{v:.4g}' if abs(v) < 1e4 else str(int(v))
    return str(v)


# ------------------------------------------------------------------------------------------------------------------ backends
class MemoryBackend:
    """An in-memory spreadsheet with the interface of GspreadBackend. `fail(name)` may raise to simulate Google being unreachable; `calls` records every call."""
    def __init__(self):
        self.data, self.calls, self.fail = {}, [], None

    def _c(self, name):
        self.calls.append(name)
        if self.fail:
            self.fail(name)

    def tabs(self):
        self._c('tabs'); return list(self.data)

    def ensure_tab(self, title, header):
        self._c('ensure_tab')
        if title not in self.data:
            self.data[title] = [list(header)]

    def values(self, title):
        self._c('values'); return [list(r) for r in self.data[title]]

    def header(self, title):
        self._c('header'); return list(self.data[title][0])

    def column(self, title, c):
        self._c('column'); return [r[c - 1] if len(r) >= c else '' for r in self.data[title]]

    def write_cells(self, title, cells):
        self._c('write_cells')
        rows = self.data[title]
        for r, c, v in cells:
            while len(rows) < r:
                rows.append([])
            while len(rows[r - 1]) < c:
                rows[r - 1].append('')
            rows[r - 1][c - 1] = v

    def append_rows(self, title, new):
        self._c('append_rows'); self.data[title].extend([list(r) for r in new])

    def add_columns(self, title, names):
        self._c('add_columns'); self.data[title][0].extend(names)


class GspreadBackend:
    """Google Sheets through gspread. `gc` is FileManager.gc (an authorised gspread client); `spreadsheet_id` is FileManager.dl_ID."""
    def __init__(self, gc, spreadsheet_id):
        self.gc, self.id, self._sh = gc, spreadsheet_id, None

    @property
    def sh(self):
        if self._sh is None:
            self._sh = self.gc.open_by_key(self.id)
        return self._sh

    def _ws(self, title):
        return self.sh.worksheet(title)

    def tabs(self):
        return [w.title for w in self.sh.worksheets()]

    def ensure_tab(self, title, header):
        if title in self.tabs():
            return
        ws = self.sh.add_worksheet(title=title, rows=1000, cols=max(26, len(header) + 2))
        ws.batch_update([{'range': f'{a1(1, 1)}:{a1(1, len(header))}', 'values': [list(header)]}], value_input_option='RAW')
        try:
            ws.freeze(rows=1)
        except Exception:
            pass

    def values(self, title):
        return self._ws(title).get_all_values()

    def header(self, title):
        return self._ws(title).row_values(1)

    def column(self, title, c):
        return self._ws(title).col_values(c)

    def write_cells(self, title, cells):
        """cells: [(row, col, value)], 1-based. Adjacent cells of a row go into one range; sent in chunks of 200 ranges (RAW: no formula interpretation)."""
        by_row = {}
        for r, c, v in cells:
            by_row.setdefault(r, {})[c] = v
        data = []
        for r, cs in sorted(by_row.items()):
            run, prev = [], None
            for c in sorted(cs):
                if prev is not None and c != prev + 1:
                    data.append({'range': f'{a1(r, run[0][0])}:{a1(r, run[-1][0])}', 'values': [[x[1] for x in run]]}); run = []
                run.append((c, cs[c])); prev = c
            data.append({'range': f'{a1(r, run[0][0])}:{a1(r, run[-1][0])}', 'values': [[x[1] for x in run]]})
        ws = self._ws(title)
        for i in range(0, len(data), 200):
            ws.batch_update(data[i:i + 200], value_input_option='RAW')

    def append_rows(self, title, rows):
        self._ws(title).append_rows(rows, value_input_option='RAW')

    def add_columns(self, title, names):
        ws = self._ws(title); header = ws.row_values(1)
        need = len(header) + len(names)
        if need > ws.col_count:
            ws.add_cols(need - ws.col_count + 2)
        ws.batch_update([{'range': f'{a1(1, len(header) + 1)}:{a1(1, need)}', 'values': [list(names)]}], value_input_option='RAW')


# ------------------------------------------------------------------------------------------------------------------ the manifest
class Manifest:
    def __init__(self, backend, tab, who, now=time.time, sleep=time.sleep, journal=None, flush_s=15.0, version='', log=print):
        self.b, self.tab, self.who, self.now, self.sleep, self.journal = backend, tab, who, now, sleep, journal
        self.flush_s, self.version, self.log = flush_s, version, log
        self.pending, self._last_flush = {}, -1e9

    # ---- structure
    def setup(self):
        """Create the tab (header: the Run Info columns, then the tracking columns) or add whatever columns are missing at the right end. Nothing is ever reordered or deleted."""
        want = RUNINFO_COLUMNS + TRACK_COLUMNS
        self.b.ensure_tab(self.tab, want)
        have = self.b.header(self.tab)
        missing = [c for c in want if c not in have]
        if missing:
            self.b.add_columns(self.tab, missing)
        return missing

    def load(self):
        """The whole tab as a DataFrame of strings, with `_row` = the row number in the sheet."""
        v = self.b.values(self.tab)
        if not v or 'Run' not in v[0]:                                      # an empty tab (a new spreadsheet has an empty 'Sheet1'), or one that is not a manifest tab
            return pd.DataFrame(columns=RUNINFO_COLUMNS + TRACK_COLUMNS + ['_row'])
        hdr = v[0]
        rows = [r + [''] * (len(hdr) - len(r)) for r in v[1:]]
        df = pd.DataFrame(rows, columns=hdr, dtype=str)
        df['_row'] = range(2, len(df) + 2)
        return df[df['Run'] != ''] if 'Run' in df else df

    # ---- import
    def import_runinfo(self, ri, in_db=(), on_cloud=()):
        """Append the runs of a Run Info DataFrame that are not in the tab yet (existing rows are never touched). A run already in the database starts as in_database (Source legacy), one whose uBAM is
        on the cloud as uploaded, every other as pending. Returns dict(added, existing)."""
        self.setup()
        hdr = self.b.header(self.tab); have = set(self.load()['Run']) if len(self.b.values(self.tab)) > 1 else set()
        new = ri[~ri['Run'].isin(have)]
        rows = []
        for _, r in new.iterrows():
            d = {c: cell(r[c]) if c in ri.columns else '' for c in RUNINFO_COLUMNS}
            run = d['Run']
            d['DownloadStatus'] = 'in_database' if run in in_db else ('uploaded' if run in on_cloud else 'pending')
            d['InDatabase'] = cell(run in in_db)
            if run in in_db or run in on_cloud:
                d['Source'] = 'legacy'
            d['UpdatedAt'] = utc_iso(self.now()); d['Notes'] = ''
            rows.append([d.get(c, '') for c in hdr])
        if rows:
            self.b.append_rows(self.tab, rows)
        return dict(added=len(rows), existing=len(ri) - len(rows))

    def reconcile(self, in_db=(), on_cloud=()):
        """Rows whose run is now in the database / on the cloud: mark them (in_database also for `uploaded` rows). Returns the number of rows changed."""
        df = self.load(); n = 0
        for r in df.itertuples():
            if r.DownloadStatus in ('in_database', 'skipped'):
                continue
            if r.Run in in_db:                                              # also upgrades `uploaded` (a --no-db run whose row has been added since)
                self.update(r.Run, DownloadStatus='in_database', InDatabase=True, ClaimedBy='', ClaimedAt=''); n += 1
            elif r.DownloadStatus == 'uploaded':
                continue
            elif r.Run in on_cloud and r.DownloadStatus not in ACTIVE:
                self.update(r.Run, DownloadStatus='uploaded', ClaimedBy='', ClaimedAt=''); n += 1
        return n

    # ---- which runs
    def runs_to_download(self, retry_failed=True, ttl_s=6 * 3600):
        """Rows to work on: pending, failed in a retryable way (or any failure with retry_failed='all'), and abandoned claims (older than ttl_s). Skip = TRUE rows and final rows are left out."""
        df = self.load(); t = self.now(); idx = []
        for i, r in enumerate(df.itertuples()):
            if str(r.Skip).upper() == 'TRUE' or r.DownloadStatus in DONE:
                continue
            ca = parse_iso(r.ClaimedAt)
            if r.DownloadStatus in ('', 'pending'):
                ok = True
            elif r.DownloadStatus == 'failed':
                ok = bool(retry_failed) and (retry_failed == 'all' or r.FailureKind in RETRYABLE)
            elif r.DownloadStatus in ACTIVE:
                ok = r.ClaimedBy == self.who or ca is None or t - ca >= ttl_s     # abandoned (or ours from an earlier crash)
            else:
                ok = False
            if ok:
                idx.append(i)
        return df.iloc[idx]

    def claim(self, run, ttl_s=6 * 3600):
        """Take a run: refuses if it is final, skipped, or claimed by someone else less than ttl_s ago. Writes the claim, then re-reads it: True only if the claim is still ours.
        If the sheet cannot be reached the run is taken UNCLAIMED (True, with a warning): a Google outage must not stop the downloads; the worst case is two machines doing the same run, which
        the database refuses (a run can be added only once) and which costs only the download."""
        try:
            df = self.load(); row = df[df['Run'] == run]
            if row.empty:
                return False
            r = row.iloc[0]; ca = parse_iso(r.ClaimedAt)
            if str(r.Skip).upper() == 'TRUE' or r.DownloadStatus in DONE:
                return False
            if r.DownloadStatus in ACTIVE and r.ClaimedBy not in ('', self.who) and ca is not None and self.now() - ca < ttl_s:
                return False
            self.update(run, DownloadStatus='claimed', ClaimedBy=self.who, ClaimedAt=utc_iso(self.now()))
            if not self.flush(force=True):
                self.log(f'  manifest {self.tab}: could not record the claim of {run}; going ahead without it'); return True
            self.sleep(random.uniform(0.2, 1.0))
            back = self.load(); b = back[back['Run'] == run]
            return (not b.empty) and b.iloc[0].ClaimedBy == self.who
        except Exception as e:
            self.log(f'  manifest {self.tab}: sheet unreachable while claiming {run} ({type(e).__name__}); going ahead without a claim')
            return True

    # ---- updates
    def update(self, run, **fields):
        d = self.pending.setdefault(run, {})
        d.update({k: cell(v) for k, v in fields.items()}); d['UpdatedAt'] = utc_iso(self.now())
        if self.version and 'PipelineVersion' not in d:
            d['PipelineVersion'] = self.version

    def _absorb_journal(self):
        if not (self.journal and os.path.exists(self.journal)):
            return
        for line in open(self.journal):
            if line.strip():
                d = json.loads(line)
                self.pending[d['run']] = {**d['fields'], **self.pending.get(d['run'], {})}       # the journal is older than what is pending now
        os.remove(self.journal)

    def flush(self, force=False):
        """Write the queued updates (one read of the RunID column and the header, one batched write). On any failure they go to the journal and False is returned; nothing is lost."""
        if not force and self.now() - self._last_flush < self.flush_s:
            return True
        self._absorb_journal()
        if not self.pending:
            self._last_flush = self.now(); return True
        try:
            hdr = self.b.header(self.tab)
            col = self.b.column(self.tab, hdr.index('Run') + 1)
            rowof = {v: i + 1 for i, v in enumerate(col) if v and i > 0}
            cells, unknown = [], set()
            for run, fields in self.pending.items():
                r = rowof.get(run)
                if r is None:
                    unknown.add(run); continue
                for k, v in fields.items():
                    if k not in hdr:
                        self.b.add_columns(self.tab, [k]); hdr = self.b.header(self.tab)
                    cells.append((r, hdr.index(k) + 1, v))
            if cells:
                self.b.write_cells(self.tab, cells)
            if unknown:
                self.log(f'  manifest {self.tab}: {len(unknown)} runs are not in the sheet (deleted?): {sorted(unknown)[:5]}')
            self.pending.clear(); self._last_flush = self.now()
            return True
        except Exception as e:
            self._to_journal()
            self.log(f'  manifest {self.tab}: could not write ({type(e).__name__}: {str(e)[:150]}); {len(self._journalled)} run updates kept in {self.journal}')
            return False

    def _to_journal(self):
        self._journalled = dict(self.pending)
        if not self.journal:                                  # no journal file: keep the updates queued and retry at the next flush
            self._last_flush = self.now(); return
        if self.journal:
            os.makedirs(os.path.dirname(os.path.abspath(self.journal)), exist_ok=True)
            with open(self.journal, 'a') as f:
                for run, fields in self.pending.items():
                    f.write(json.dumps(dict(tab=self.tab, run=run, fields=fields)) + '\n')
        self.pending.clear(); self._last_flush = self.now()
