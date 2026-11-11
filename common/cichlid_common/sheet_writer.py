"""Guarded writing of computed values into COLUMNS of an existing Google worksheet (the AlignmentDatabase now, the SampleDatabase next).

The sheets are shared and partly hand-curated, and pipelines append to them while we work, so this writer does the least that can do the job and checks what it did:

  * it touches ONLY the cells you name: `{SampleID: {column: value}}`, found by the SampleID in the key column and the column's name in the header row; never a whole row, never another cell;
  * the default is a PLAN that writes nothing: for every cell it says `write` (the cell is empty), `same` (already that value), `conflict` (it holds something else: left alone unless overwrite=True),
    `overwrite`, `no_row` (the sample is not in the sheet), `duplicate_row` (the sample is in more than one row: ambiguous, left alone) or `no_column`;
  * applying re-reads the sheet and re-plans immediately before writing (the row numbers come from THAT read, so rows appended or inserted since the plan do not matter), writes with
    `value_input_option='RAW'` in a few batched calls (idempotent: gspread_retry repeats them after a network error or HTTP 429), then reads the sheet again and checks that every written cell
    sits in the row of the right sample and holds the value; a mismatch is reported loudly with the cell addresses;
  * missing columns are added at the END of the header (the same calls as the alignment pipeline's `setup-database`), which running pipelines do not mind: they lay a row out by the header they read at
    that moment, so a column they do not know stays blank in their rows;
  * a journal (TSV) lists every planned and every verified change, with the old value, so a mistake can be undone by hand.

The one thing it cannot do is write a cell only if its row has not moved: Sheets offers no conditional write. The window is the few milliseconds between the final read and the write (it writes
through the worksheet it has just read), and pipelines only APPEND, which moves nothing. If a person inserted or deleted a row in exactly that moment, a value could land in the neighbouring row;
the check afterwards catches it and names the cells, but the value that cell held before is not in the journal. That is harmless for columns that are ours (new, empty except for our values) and is
why the SampleDatabase, whose columns hold hand-entered data, will be written with --limit trials first and outside working hours if it matters.

Numbers are written as numbers. A cell that a column format displays rounded ('4.01' for 4.0123) still counts as the same value: the comparison tolerates the rounding the display shows."""
import datetime
import os
import re


class SheetWriteError(RuntimeError):
    pass


WRITABLE = ('write', 'overwrite')


def col_letter(n):
    """1 -> A, 26 -> Z, 27 -> AA."""
    if n < 1:
        raise ValueError('a column number starts at 1')
    s = ''
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def a1(row, col):
    return f'{col_letter(col)}{row}'


def is_blank(text):
    return text is None or str(text).strip() == ''


def same_value(text, value):
    """Does the cell text `text` hold `value`? A number matches a text that parses to it up to the rounding that text shows; anything else must match as stripped text."""
    t = '' if text is None else str(text).strip()
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        pct = t.endswith('%')
        body = t.rstrip('%').replace(',', '')
        try:
            x = float(body) / (100.0 if pct else 1.0)
        except ValueError:
            return False
        decimals = (len(body.split('.')[1]) if '.' in body else 0) + (2 if pct else 0)
        return abs(x - float(value)) <= 0.5 * 10 ** -decimals + 1e-12
    return t == str(value).strip()


class Plan:
    """The cells of an update and what would happen to each: `changes` is a list of dicts with sample, column, row, col, old, new, status."""

    def __init__(self, changes, header_missing, header_used, ws=None):
        self.changes, self.header_missing, self.header_used, self.ws = changes, list(header_missing), header_used, ws

    @property
    def writable(self):
        return [c for c in self.changes if c['status'] in WRITABLE]

    def count(self, status=None, column=None):
        return sum(1 for c in self.changes if (status is None or c['status'] == status) and (column is None or c['column'] == column))

    def format(self, preview=8):
        cols = list(dict.fromkeys(c['column'] for c in self.changes))
        lines = [f'{len(self.changes)} cells in {len({c["sample"] for c in self.changes})} samples; ' + ', '.join(f'{s} {self.count(s)}' for s in
                 ('write', 'overwrite', 'same', 'conflict', 'no_row', 'duplicate_row', 'no_column') if self.count(s))]
        if self.header_missing:
            lines.append(f'columns to ADD at the end of the header: {", ".join(self.header_missing)}')
        for col in cols:
            lines.append(f'  {col}: ' + ', '.join(f'{s} {self.count(s, col)}' for s in ('write', 'overwrite', 'same', 'conflict') if self.count(s, col)))
        for status, what in (('conflict', 'a different value is already there (left alone)'), ('no_row', 'the sample is not in the sheet'), ('duplicate_row', 'the sample is in more than one row (left alone)')):
            bad = [c for c in self.changes if c['status'] == status]
            if bad:
                names = ', '.join(sorted({c['sample'] for c in bad})[:preview])
                lines.append(f'  {what}: {len({c["sample"] for c in bad})} samples' + (f' ({names}{", ..." if len({c["sample"] for c in bad}) > preview else ""})' if names else ''))
        w = self.writable
        if w:
            lines.append('  first cells to write:')
            for c in w[:preview]:
                lines.append(f"    {a1(c['row'], c['col']):>7}  {c['sample']:<18}{c['column']:<22}{c['old'] or '(empty)'} -> {c['new']}")
        return '\n'.join(lines)


class Report:
    def __init__(self, plan, written, failures, added_columns, journal):
        self.plan, self.written, self.failures, self.added_columns, self.journal = plan, written, failures, added_columns, journal

    @property
    def ok(self):
        return not self.failures


class ColumnUpdater:
    """`worksheet_fn()` returns a gspread Worksheet (or a test double with get_all_values, col_count, add_cols, update_cell and batch_update). `key` is the header of the column that identifies a row."""

    def __init__(self, worksheet_fn, key='SampleID', chunk=200, log=print):
        self.worksheet_fn, self.key, self.chunk, self.log = worksheet_fn, key, chunk, log

    # ---------------------------------------------------------------- reading
    def _read(self):
        ws = self.worksheet_fn()
        values = ws.get_all_values()
        if not values:
            raise SheetWriteError('the worksheet is empty: there is no header row')
        header = [str(h).strip() for h in values[0]]
        used = len(header)
        while used and not header[used - 1]:
            used -= 1
        return ws, header, used, values[1:]

    # ---------------------------------------------------------------- the plan
    def plan(self, updates, overwrite=False):
        """What `apply` would do, without writing anything. `updates` is {sample: {column: value}}."""
        ws, header, used, rows = self._read()
        wanted = list(dict.fromkeys(c for cols in updates.values() for c in cols))
        for name in [self.key] + wanted:
            if header.count(name) > 1:
                raise SheetWriteError(f'the header has the column {name!r} more than once: the right one cannot be told')
        if self.key not in header:
            raise SheetWriteError(f'the worksheet has no key column {self.key!r}')
        kc = header.index(self.key)
        missing = [c for c in wanted if c not in header]
        col_of = {c: header.index(c) + 1 for c in wanted if c in header}
        col_of.update({c: used + 1 + i for i, c in enumerate(missing)})                  # where a missing column will be
        by_key = {}
        for i, r in enumerate(rows):
            k = (r[kc] if kc < len(r) else '').strip()
            if k:
                by_key.setdefault(k, []).append(i + 2)
        changes = []
        for sample, cols in updates.items():
            sample = str(sample).strip()
            for column, value in cols.items():
                c = dict(sample=sample, column=column, row=None, col=col_of.get(column), old='', new=value, status='')
                nums = by_key.get(sample, [])
                if not nums:
                    c['status'] = 'no_row'
                elif len(nums) > 1:
                    c['status'] = 'duplicate_row'
                else:
                    c['row'] = nums[0]
                    r = rows[nums[0] - 2]
                    c['old'] = '' if column in missing else (r[col_of[column] - 1] if col_of[column] - 1 < len(r) else '')
                    if is_blank(c['old']):
                        c['status'] = 'write'
                    elif same_value(c['old'], value):
                        c['status'] = 'same'
                    else:
                        c['status'] = 'overwrite' if overwrite else 'conflict'
                changes.append(c)
        return Plan(changes, missing, used, ws)

    # ---------------------------------------------------------------- columns
    def add_columns(self, names):
        """Append the header cells that are missing, at the end of the header (the alignment pipeline's setup-database does the same). Returns what was added."""
        ws, header, used, _ = self._read()
        missing = [n for n in dict.fromkeys(names) if n not in header]
        if not missing:
            return []
        if used + len(missing) > ws.col_count:
            ws.add_cols(used + len(missing) - ws.col_count)
        for i, name in enumerate(missing):
            ws.update_cell(1, used + 1 + i, name)
        return missing

    # ---------------------------------------------------------------- writing
    def apply(self, updates, overwrite=False, limit=None, journal=None):
        """Add missing columns, re-plan on a fresh read, write the writable cells, verify them. `limit` keeps only the first N samples (alphabetically): a first trial on the live sheet."""
        first = self.plan(updates, overwrite)
        added = self.add_columns(first.header_missing) if first.header_missing else []
        plan = self.plan(updates, overwrite)
        todo = plan.writable
        if limit is not None:
            keep = set(sorted({c['sample'] for c in todo})[:limit])
            todo = [c for c in todo if c['sample'] in keep]
        self._journal(journal, 'planned', todo)
        ws = plan.ws                                                  # the worksheet that was just read: no further round trip between that read and the write
        for i in range(0, len(todo), self.chunk):
            batch = todo[i:i + self.chunk]
            ws.batch_update([{'range': a1(c['row'], c['col']), 'values': [[c['new']]]} for c in batch], value_input_option='RAW')
        failures = self._verify(todo)
        self._journal(journal, 'FAILED' if failures else 'verified', todo, failures)
        return Report(plan, len(todo) - len(failures), failures, added, journal)

    def _verify(self, todo):
        if not todo:
            return []
        _, header, _, rows = self._read()
        kc = header.index(self.key)
        failures = []
        for c in todo:
            r = rows[c['row'] - 2] if c['row'] - 2 < len(rows) else []
            owner = (r[kc] if kc < len(r) else '').strip()
            cell = r[c['col'] - 1] if c['col'] - 1 < len(r) else ''
            if owner != c['sample']:
                failures.append(dict(c, problem=f"row {c['row']} now belongs to {owner or '(nobody)'}, not {c['sample']}: the sheet changed while writing"))
            elif not same_value(cell, c['new']):
                failures.append(dict(c, problem=f"{a1(c['row'], c['col'])} reads back {cell!r}, not {c['new']!r}"))
        return failures

    def _journal(self, path, phase, changes, failures=()):
        if not path:
            return
        os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
        new = not os.path.exists(path)
        bad = {(f['sample'], f['column']): f['problem'] for f in failures}
        with open(path, 'a') as f:
            if new:
                f.write('time\tphase\tsample\tcolumn\tcell\told\tnew\tnote\n')
            now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            for c in changes:
                note = bad.get((c['sample'], c['column']), c['status'] if phase == 'planned' else '')
                f.write('\t'.join(str(x).replace('\t', ' ').replace('\n', ' ') for x in (now, phase, c['sample'], c['column'], a1(c['row'], c['col']), c['old'], c['new'], note)) + '\n')
