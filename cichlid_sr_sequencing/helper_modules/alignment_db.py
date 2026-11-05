"""Writing one sample's row to the AlignmentDatabase (a worksheet per genome version of the Google spreadsheet).

Differences from FileManager.addAlignmentRow, which this replaces for the new pipeline:
  * extra (new) columns are allowed; the legacy columns must all be present;
  * the worksheet is RESIZED to the data (set_with_dataframe(resize=True)): the new pipeline adds about 20 columns and a default 26-column worksheet would reject the write;
  * a failed write raises (FileManager.setAlignmentDatabase only printed the exception and returned False, and the caller ignored it);
  * a sample that is already in the database is left alone (returns 'exists') instead of raising, so a rerun after a crash between the write and the cleanup is harmless;
  * one lock: samples finish on different threads and every write is read-modify-write of the whole sheet."""
import threading
import time

import pandas as pd

LEGACY_KEYS = ('SampleID', 'GenomeVersion', 'RunIDs', 'Coverage', 'TotalReads', 'UnmappedReads', 'DiscordantReads', 'InversionReads', 'DuplicationReads', 'ClippedReads',
               'ChimericReads', 'minimap2_version', 'gatk_version', 'pysam_version', 'BamSize')


class DatabaseWriter:
    def __init__(self, fm, attempts=4, backoff=10.0, sleep=time.sleep):
        self.fm, self.attempts, self.backoff, self.sleep, self.lock = fm, attempts, backoff, sleep, threading.Lock()

    def _write(self):
        from gspread_dataframe import set_with_dataframe
        last = None
        for i in range(self.attempts):
            try:
                ws = self.fm.gc.open_by_key(self.fm.a_ID).worksheet(self.fm.genome_version)
                set_with_dataframe(ws, self.fm.alignment_dt, resize=True)
                return
            except Exception as e:                                           # rate limits and network errors are common: retry
                last = e
                self.sleep(self.backoff * (i + 1))
        raise RuntimeError(f'could not write the alignment database after {self.attempts} attempts: {type(last).__name__}: {last}')

    def add_row(self, row):
        missing = [k for k in LEGACY_KEYS if k not in row]
        if missing:
            raise ValueError(f'alignment row is missing {missing}')
        with self.lock:
            self.fm.readAlignmentDatabase()
            if row['SampleID'] in set(self.fm.alignment_dt.SampleID.dropna()):
                return 'exists'
            self.fm.alignment_dt = pd.concat([self.fm.alignment_dt, pd.DataFrame([row])], ignore_index=True)
            self._write()
        return 'added'


class AppendingDatabaseWriter:
    """Writes one sample's row by APPENDING it to the worksheet with a single atomic Sheets call, then reads it back to verify.

    Why not DatabaseWriter (read the whole sheet, add a row, write the whole sheet back): with two machines finishing samples every few minutes, two whole-sheet writes that overlap
    in time lose a row (the later write overwrites the earlier one with a copy that does not contain it), and a lost row means the sample is aligned a second time. An append never
    touches other rows, so concurrent appends from two machines both survive. Rewriting ~1,000 rows for every sample is also slow and counts against the API quota.

    The row is laid out by the worksheet's HEADER ROW (column order of the sheet, not of our dict). New columns must be added to the header once, before the runs start:
    `alignFastQReads_v2.py <Genome> --setup-database`. A run refuses to start writing if a column it needs is missing. After every append the SampleID column is read again:
    the sample must be there exactly once (a missing row is appended again, a duplicate is reported)."""
    NEEDED = None                                    # set below: legacy columns + the new statistics columns

    def __init__(self, fm, attempts=4, backoff=10.0, sleep=time.sleep, worksheet_fn=None, log=print):
        self.fm, self.attempts, self.backoff, self.sleep, self.log = fm, attempts, backoff, sleep, log
        self.worksheet_fn = worksheet_fn or (lambda: fm.gc.open_by_key(fm.a_ID).worksheet(fm.genome_version))
        self.lock = threading.Lock()

    @staticmethod
    def needed_columns():
        from . import bam_stats as bs
        return list(LEGACY_KEYS) + [c for c in bs.NEW_COLUMNS if c not in LEGACY_KEYS]

    @staticmethod
    def _cell(v):
        """A value Google Sheets accepts: plain python numbers/strings; empty for NaN/None."""
        try:
            if v is None or (isinstance(v, float) and v != v):
                return ''
            if hasattr(v, 'item'):
                v = v.item()                         # numpy scalar -> python
            if isinstance(v, float) and v != v:
                return ''
            return v
        except Exception:
            return str(v)

    def check_columns(self):
        """Columns this pipeline writes that the worksheet's header lacks (empty list = ready)."""
        header = self.worksheet_fn().row_values(1)
        return [c for c in self.needed_columns() if c not in header]

    def setup_columns(self, columns=None):
        """Add the columns the new pipeline writes to the header row (idempotent). Returns the list of columns that were added."""
        columns = columns or self.needed_columns()
        with self.lock:
            ws = self.worksheet_fn()
            header = ws.row_values(1)
            missing = [c for c in columns if c not in header]
            if not missing:
                return []
            if len(header) + len(missing) > ws.col_count:
                ws.add_cols(len(header) + len(missing) - ws.col_count)
            for i, name in enumerate(missing):
                ws.update_cell(1, len(header) + 1 + i, name)
            return missing

    def add_row(self, row):
        missing_legacy = [k for k in LEGACY_KEYS if k not in row]
        if missing_legacy:
            raise ValueError(f'alignment row is missing {missing_legacy}')
        sample = row['SampleID']
        with self.lock:
            last = None
            for attempt in range(self.attempts):
                try:
                    ws = self.worksheet_fn()
                    header = ws.row_values(1)
                    lacking = [c for c in row if c not in header]
                    if lacking:
                        raise RuntimeError(f'the worksheet {self.fm.genome_version!r} has no column {lacking}: run `alignFastQReads_v2.py {self.fm.genome_version} --setup-database` once first')
                    col = header.index('SampleID') + 1
                    ids = ws.col_values(col)
                    if sample in ids:                                                       # already there: a rerun, or another machine got there first
                        return 'exists' if attempt == 0 else 'added'
                    ws.append_rows([[self._cell(row.get(h, '')) for h in header]], value_input_option='RAW', table_range='A1')
                    ids = ws.col_values(col)                                                # read it back: present, and exactly once
                    n = ids.count(sample)
                    if n == 0:
                        raise RuntimeError('the row was not found after appending')
                    if n > 1:
                        self.log(f'WARNING: {sample} now appears {n} times in the alignment database (another machine added it at the same moment); remove the extra rows by hand')
                    return 'added'
                except RuntimeError as e:
                    if 'run `alignFastQReads_v2.py' in str(e):
                        raise                                                               # a setup problem, not a transient one
                    last = e
                except Exception as e:                                                      # rate limits and network errors are common: retry
                    last = e
                self.sleep(self.backoff * (attempt + 1))
            raise RuntimeError(f'could not write {sample} to the alignment database after {self.attempts} attempts: {type(last).__name__}: {last}')
