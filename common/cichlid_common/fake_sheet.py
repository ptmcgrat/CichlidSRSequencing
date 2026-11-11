"""A test double for a gspread Worksheet, used by the tests of the modules that write sheets (nothing in production imports it).

It keeps the cells as a list of rows, answers the calls ColumnUpdater makes (get_all_values, col_count, add_cols, update_cell, batch_update) and records every write call, so a test can show that a plan
wrote nothing, that only the named cells changed, and how many API calls a write took. Hooks simulate what the real sheets do: `before_write` runs before the first batch is applied (someone edits the
sheet between the plan and the write), `after_write` runs after the batches and before the check (rows shift), `show_decimals` shows numbers rounded the way a formatted column does."""
import re


class FakeWorksheet:
    title = 'Mzebra_GT3_NCBI'

    def __init__(self, header, rows, extra_columns=0, show_decimals=None):
        self.grid = [list(header) + [''] * extra_columns] + [list(r) + [''] * (len(header) + extra_columns - len(r)) for r in rows]
        self.calls = []
        self.show_decimals = show_decimals or {}
        self.before_write = None
        self.after_write = None
        self._fired = False

    @property
    def col_count(self):
        return max(len(r) for r in self.grid)

    def get_all_values(self):
        w = self.col_count
        out = []
        for ri, r in enumerate(self.grid):
            row = []
            for ci in range(w):
                v = r[ci] if ci < len(r) else ''
                if ri and isinstance(v, float) and self.grid[0][ci] in self.show_decimals:
                    v = f'{v:.{self.show_decimals[self.grid[0][ci]]}f}'
                row.append('' if v is None else str(v))
            out.append(row)
        return out

    def row_values(self, i):
        r = self.get_all_values()[i - 1]
        while r and not r[-1]:
            r.pop()
        return r

    def add_cols(self, n):
        self.calls.append(('add_cols', n))
        for r in self.grid:
            r.extend([''] * n)

    def update_cell(self, row, col, value):
        self.calls.append(('update_cell', row, col, value))
        while len(self.grid) < row:
            self.grid.append([''] * self.col_count)
        self.grid[row - 1][col - 1] = value

    def batch_update(self, data, value_input_option=None):
        if self.before_write and not self._fired:
            self._fired = True
            self.before_write(self)
        self.calls.append(('batch_update', len(data), value_input_option))
        for d in data:
            m = re.fullmatch(r'([A-Z]+)(\d+)', d['range'])
            col = 0
            for ch in m.group(1):
                col = col * 26 + ord(ch) - 64
            self.grid[int(m.group(2)) - 1][col - 1] = d['values'][0][0]
        if self.after_write:
            self.after_write(self)

    def writes(self):
        return [c for c in self.calls if c[0] in ('batch_update', 'update_cell', 'add_cols')]

    def column(self, name):
        i = self.grid[0].index(name)
        return [r[i] for r in self.grid[1:]]

    def insert_row_at_top(self, row):
        """A person inserts a row under the header."""
        self.grid.insert(1, list(row) + [''] * (self.col_count - len(row)))
