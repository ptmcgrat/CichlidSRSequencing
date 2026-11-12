"""The ledger: one row per aligned sample in the AlignmentDatabase (a worksheet per genome version of the Google spreadsheet).

What the sheet carries (decisions 37 and 41): the columns kept from before (SampleID, GenomeVersion, RunIDs, Coverage, TotalReads, the tool versions, BamSize), DuplicationRate and Species, and the
columns agreed with them (usable coverage, mapped fraction, QC warnings, panel status, provenance, the runs' count and fingerprint, instrument and read length). The six split-read fractions are gone.
Everything else the pipeline measures (about twenty statistics) goes to <sample>.qc.json next to the BAM, not to the sheet. Whether the sheet should carry more of it is an open decision, so the
column list is ONE constant, LEDGER_COLUMNS, and `ledger_row()` is the one place that maps the statistics to it.

Writing a row (AppendingDatabaseWriter): a single atomic append, then a read-back. A failed write raises. A sample that is already in the sheet is left alone (a rerun after a crash between the write
and the cleanup is harmless). Old rows keep their old columns (the six fractions and the rest stay in the sheet, blank for new rows)."""
import threading
import time

REQUIRED_KEYS = ('SampleID', 'GenomeVersion', 'RunIDs', 'Coverage', 'TotalReads', 'minimap2_version', 'gatk_version', 'pysam_version', 'BamSize')
LEDGER_COLUMNS = ('SampleID', 'GenomeVersion', 'Species', 'RunIDs', 'InputFingerprint', 'Coverage', 'TotalReads', 'MappedFraction', 'DuplicationRate', 'BamSize',
                  'QCWarnings', 'PanelStatus', 'AlignedOn', 'Machine', 'PipelineVersion', 'minimap2_version', 'gatk_version', 'pysam_version')
# Columns the pipeline wrote until 2026-11 and no longer does (SPEC decision 75): UsableCoverage is replaced by EffectiveCoverage (module 05, `effcov`, the depth the variant caller sees), NRuns by RunIDs, and Instrument and
# ReadLength live in the ReadsDatabase. The pipeline still READS the instrument and the read length (the coverage of the legacy columns needs the read length); it just does not write them to the sheet. A sheet that still has
# these columns is fine (their cells stay blank in new rows); a sheet without them is fine too. They are dropped from the sheet by hand, after this code is on every machine.
RETIRED_COLUMNS = ('NRuns', 'UsableCoverage', 'Instrument', 'ReadLength')
FROM_QC = {'DuplicationRate': 'DuplicateFraction'}              # ledger name <- name in qc.json (v6's statistics keep their names there)


def ledger_row(qc, meta):
    """The row of the sheet: every LEDGER_COLUMNS key (empty when unknown), taken from `meta` (what is known about the sample and the run) and `qc` (the statistics of qc.json)."""
    row = {}
    for col in LEDGER_COLUMNS:
        if col in meta:
            row[col] = meta[col]
        elif FROM_QC.get(col, col) in qc:
            row[col] = qc[FROM_QC.get(col, col)]
        else:
            row[col] = ''
    return row


class AppendingDatabaseWriter:
    """Writes one sample's row by APPENDING it to the worksheet with a single atomic Sheets call, then reads it back to verify.

    Why not read the whole sheet, add a row and write the whole sheet back: with two machines finishing samples every few minutes, two whole-sheet writes that overlap in time lose a row (the
    later write overwrites the earlier one with a copy that does not contain it), and a lost row means the sample is aligned a second time. An append never touches other rows, so concurrent
    appends from two machines both survive. Rewriting ~1,000 rows for every sample is also slow and counts against the API quota.

    The row is laid out by the worksheet's HEADER ROW (column order of the sheet, not of our dict). New columns must be added to the header once, before the runs start:
    `python -m align_pipeline setup-database <Genome>`. A run refuses to start writing if a column it needs is missing. After every append the SampleID column is read again:
    the sample must be there exactly once (a missing row is appended again, a duplicate is reported)."""

    def __init__(self, fm, attempts=4, backoff=10.0, sleep=time.sleep, worksheet_fn=None, log=print):
        self.fm, self.attempts, self.backoff, self.sleep, self.log = fm, attempts, backoff, sleep, log
        self.worksheet_fn = worksheet_fn or (lambda: fm.gc.open_by_key(fm.a_ID).worksheet(fm.genome_version))
        self.lock = threading.Lock()

    @staticmethod
    def needed_columns():
        return list(LEDGER_COLUMNS)

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
        """Add the columns the pipeline writes to the header row (idempotent). Returns the list of columns that were added."""
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
        missing_required = [k for k in REQUIRED_KEYS if k not in row]
        if missing_required:
            raise ValueError(f'alignment row is missing {missing_required}')
        sample = row['SampleID']
        with self.lock:
            last = None
            for attempt in range(self.attempts):
                try:
                    ws = self.worksheet_fn()
                    header = ws.row_values(1)
                    lacking = [c for c in row if c not in header]
                    if lacking:
                        raise RuntimeError(f'the worksheet {self.fm.genome_version!r} has no column {lacking}: run `python -m align_pipeline setup-database {self.fm.genome_version}` once first')
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
                    if 'setup-database' in str(e):
                        raise                                                               # a setup problem, not a transient one
                    last = e
                except Exception as e:                                                      # rate limits and network errors are common: retry
                    last = e
                self.sleep(self.backoff * (attempt + 1))
            raise RuntimeError(f'could not write {sample} to the alignment database after {self.attempts} attempts: {type(last).__name__}: {last}')
