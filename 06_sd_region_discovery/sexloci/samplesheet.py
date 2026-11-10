"""Write ONE new column (default 'gsdf number') into the SampleDatabase worksheet, safely.

This is the only place in sexloci that can modify a lab sheet, and it only runs when `gsdf-genotype --write-sample-database` is given.
What it guarantees:
  * dry run by default: without `--yes` nothing is written, the plan is only printed and saved;
  * only the cells of that one column are touched (changed cells only, in a single batch request); no other column or row is edited;
  * existing values in the column are kept unless `--overwrite`; conflicts are reported;
  * rows are matched by SampleID, never by position, and the write is refused if the SampleID column changed since it was read
    (someone sorted or inserted rows while the script ran);
  * a full snapshot of the sheet is saved before any write, and every change is logged (old and new value);
  * the column is read back after writing and compared with what was intended;
  * values are sent as RAW: numbers stay numbers and text stays text. (USER_ENTERED would let Google Sheets turn a genotype like
    '1/1' or '0/1' into a date.)
Several columns can be written in one guarded batch (`update_sample_columns`). The service account in the lab credential file needs EDIT
access to the spreadsheet.
"""
from __future__ import annotations

import csv
import datetime
import os
import re
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple

import pandas as pd

DEFAULT_COLUMN = "gsdf number"


def _norm(x) -> str:
    return re.sub(r"\s+", " ", str(x).strip().lower())


def a1(row: int, col: int) -> str:
    letters = ""
    while col:
        col, r = divmod(col - 1, 26)
        letters = chr(65 + r) + letters
    return f"{letters}{row}"


@dataclass
class ColumnPlan:
    column: str
    col_index: int                                        # 1-based column that will hold the values
    is_new: bool
    id_col: int                                           # 1-based SampleID column
    id_snapshot: List[str] = field(default_factory=list)  # SampleID column as read, row 1 = header
    to_write: List[Tuple[int, str, str, Any]] = field(default_factory=list)       # (row, SampleID, old, new)
    conflicts: List[Tuple[int, str, str, Any]] = field(default_factory=list)      # existing value differs, not overwritten
    n_equal: int = 0
    not_in_sheet: List[str] = field(default_factory=list)


def _clean(v):
    """int stays int (numpy ints too), everything else becomes text."""
    try:
        import numpy as np
        if isinstance(v, np.integer):
            return int(v)
    except ImportError:
        pass
    return v if isinstance(v, int) and not isinstance(v, bool) else str(v)


def plan_column_update(grid: List[List[str]], values: Dict[str, Any], column: str = DEFAULT_COLUMN, overwrite: bool = False,
                       id_header: str = "SampleID", first_free_col: Optional[int] = None) -> ColumnPlan:
    """Decide which cells to write from a snapshot of the sheet (list of rows, as returned by get_all_values). values: SampleID -> int or
    text. first_free_col (1-based): where a NEW column may go at the earliest, so that several new columns get consecutive positions."""
    if not grid:
        raise ValueError("the SampleDatabase worksheet is empty")
    header = [str(h) for h in grid[0]]
    idx = {}
    for i, h in enumerate(header):
        idx.setdefault(_norm(h), []).append(i)
    if len(idx.get(_norm(id_header), [])) != 1:
        raise ValueError(f"the worksheet needs exactly one '{id_header}' column (found {len(idx.get(_norm(id_header), []))})")
    id_c = idx[_norm(id_header)][0]
    found = idx.get(_norm(column), [])
    if len(found) > 1:
        raise ValueError(f"the worksheet has {len(found)} columns named '{column}'; refusing to guess which to update")
    last = max((i for i, h in enumerate(header) if str(h).strip()), default=-1)
    col_c = found[0] if found else max(last + 1, (first_free_col or 1) - 1)
    plan = ColumnPlan(column=column, col_index=col_c + 1, is_new=not found, id_col=id_c + 1)
    plan.id_snapshot = [(r[id_c].strip() if len(r) > id_c else "") for r in grid]
    seen = set()
    for rno, row in enumerate(grid[1:], start=2):
        sid = plan.id_snapshot[rno - 1]
        if sid not in values:
            continue
        seen.add(sid)
        new = _clean(values[sid])
        old = (row[col_c] if len(row) > col_c else "").strip()
        if old == "":
            plan.to_write.append((rno, sid, old, new))
        elif old == str(new):
            plan.n_equal += 1
        elif overwrite:
            plan.to_write.append((rno, sid, old, new))
        else:
            plan.conflicts.append((rno, sid, old, new))
    plan.not_in_sheet = sorted(set(values) - seen)
    return plan


def _retry(fn, what: str, attempts: int = 3):
    last = None
    for i in range(attempts):
        try:
            return fn()
        except Exception as e:                                    # transient API / rate-limit errors
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"{what} failed after {attempts} attempts: {last}") from last


def apply_column_plans(ws, plans: List[ColumnPlan], notes: Dict[str, str], log: Callable = print) -> None:
    """Write the cells of every plan in ONE batch request, after checking that the SampleID column is unchanged."""
    if not plans:
        return
    id_col, snap = plans[0].id_col, plans[0].id_snapshot
    now = [c.strip() for c in _retry(lambda: ws.col_values(id_col), "re-reading the SampleID column")]
    n = max(len(now), len(snap))
    if now + [""] * (n - len(now)) != snap + [""] * (n - len(snap)):
        raise RuntimeError("the SampleID column of the worksheet changed while this ran (rows were sorted, inserted or deleted); "
                           "nothing was written. Rerun the command.")
    data = []
    need_cols = max(p.col_index for p in plans)
    if any(p.is_new for p in plans) and need_cols > ws.col_count:
        _retry(lambda: ws.add_cols(need_cols - ws.col_count), "adding columns")
    for plan in plans:
        if plan.is_new:
            data.append({"range": a1(1, plan.col_index), "values": [[plan.column]]})
        for rno, _sid, _old, new in plan.to_write:
            data.append({"range": a1(rno, plan.col_index), "values": [[new]]})
    if not data:
        log("[sheet] nothing to write")
        return
    _retry(lambda: ws.batch_update(data, value_input_option="RAW"), "writing the columns")
    for plan in plans:
        if plan.is_new and notes.get(plan.column):
            try:
                ws.update_note(a1(1, plan.col_index), notes[plan.column])
            except Exception:
                pass                                               # notes are optional (older gspread)
        back = _retry(lambda: ws.col_values(plan.col_index), f"reading column '{plan.column}' back")
        bad = [(r, new, back[r - 1] if len(back) >= r else "") for r, _s, _o, new in plan.to_write
               if (back[r - 1] if len(back) >= r else "").strip() != str(new)]
        if bad:
            raise RuntimeError(f"read-back check failed for column '{plan.column}': {len(bad)} cells, e.g. row {bad[0][0]}: "
                               f"wrote {bad[0][1]}, found '{bad[0][2]}'")
        log(f"[sheet] wrote {len(plan.to_write)} cells in column '{plan.column}'" + (" (new column created)" if plan.is_new else "")
            + "; read-back check passed")


def apply_column_plan(ws, plan: ColumnPlan, outdir: str, note: Optional[str] = None, log: Callable = print) -> None:
    apply_column_plans(ws, [plan], {plan.column: note or ""}, log)


def update_sample_columns(ws, specs: List[Dict[str, Any]], outdir: str, overwrite: bool = False, apply: bool = False, version: str = "",
                          log: Callable = print) -> List[ColumnPlan]:
    """Several columns at once. specs: [{'column': name, 'values': {SampleID: int|text}, 'note': text}, ...]. Dry run unless apply=True;
    one snapshot, one guarded batch write, one plan file per column."""
    os.makedirs(outdir, exist_ok=True)
    grid = _retry(lambda: ws.get_all_values(), "reading the SampleDatabase worksheet")
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    plans, free = [], None
    for sp in specs:
        plan = plan_column_update(grid, sp["values"], sp["column"], overwrite, first_free_col=free)
        if plan.is_new:
            free = plan.col_index + 1
        plans.append(plan)
        plan_path = os.path.join(outdir, f"sample_database_{sp['column'].replace(' ', '_')}_plan_{stamp}.tsv")
        with open(plan_path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["row", "SampleID", "old_value", "new_value", "action"])
            for r, sid, o, nw in plan.to_write:
                w.writerow([r, sid, o, nw, "write" if o == "" else "overwrite"])
            for r, sid, o, nw in plan.conflicts:
                w.writerow([r, sid, o, nw, "kept (existing value differs; use --overwrite)"])
            for sid in plan.not_in_sheet:
                w.writerow(["", sid, "", sp["values"][sid], "not in the sheet"])
        log(f"[sheet] column '{sp['column']}': {'will be CREATED at the end of the header' if plan.is_new else 'exists (column %d)' % plan.col_index}"
            f" | {len(sp['values'])} values | to write {len(plan.to_write)}, already equal {plan.n_equal}, "
            f"existing value differs and is kept {len(plan.conflicts)}, not found in the sheet {len(plan.not_in_sheet)}  (plan: {plan_path})")
    if not apply:
        log("[sheet] DRY RUN: nothing was written. Review the plan files, then rerun with --yes to write.")
        return plans
    backup = os.path.join(outdir, f"sample_database_backup_{stamp}.csv")
    with open(backup, "w", newline="") as fh:
        csv.writer(fh).writerows(grid)
    log(f"[sheet] snapshot of the whole sheet saved to {backup}")
    notes = {sp["column"]: (sp.get("note") or f"sexloci {version}, {datetime.date.today().isoformat()}") for sp in specs}
    apply_column_plans(ws, plans, notes, log)
    return plans


def update_sample_database(ws, results: pd.DataFrame, outdir: str, column: str = DEFAULT_COLUMN, include_flagged: bool = False,
                           overwrite: bool = False, apply: bool = False, version: str = "", log: Callable = print) -> ColumnPlan:
    """results: the gsdf-genotype table (SampleID, CopyNumber, Genotype, Flag). Dry run unless apply=True."""
    os.makedirs(outdir, exist_ok=True)
    called = results[results.Genotype.fillna("").astype(str) != ""].copy()
    called["Flag"] = called.Flag.fillna("")
    n_flagged = int((called.Flag != "").sum())
    use = called if include_flagged else called[called.Flag == ""]
    values = {str(r.SampleID).strip(): int(r.CopyNumber) for r in use.itertuples()}
    grid = _retry(lambda: ws.get_all_values(), "reading the SampleDatabase worksheet")
    plan = plan_column_update(grid, values, column, overwrite)
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    plan_path = os.path.join(outdir, f"sample_database_{column.replace(' ', '_')}_plan_{stamp}.tsv")
    with open(plan_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["row", "SampleID", "old_value", "new_value", "action"])
        for r, s, o, nw in plan.to_write:
            w.writerow([r, s, o, nw, "write" if o == "" else "overwrite"])
        for r, s, o, nw in plan.conflicts:
            w.writerow([r, s, o, nw, "kept (existing value differs; use --overwrite)"])
        for s in plan.not_in_sheet:
            w.writerow(["", s, "", values[s], "not in the sheet"])
    log(f"[sheet] column '{column}': {'will be CREATED at the end of the header' if plan.is_new else 'exists (column %d)' % plan.col_index}")
    log(f"[sheet] {len(values)} samples have a value to write ({n_flagged} flagged calls "
        f"{'included' if include_flagged else 'left blank; --include-flagged writes them'}); {len(results) - len(called)} samples have no call")
    log(f"[sheet] cells to write: {len(plan.to_write)} | already equal: {plan.n_equal} | existing value differs, kept: {len(plan.conflicts)} "
        f"| samples not found in the sheet: {len(plan.not_in_sheet)}")
    log(f"[sheet] plan saved to {plan_path}")
    if not apply:
        log("[sheet] DRY RUN: nothing was written. Review the plan, then rerun with --yes to write.")
        return plan
    backup = os.path.join(outdir, f"sample_database_backup_{stamp}.csv")
    with open(backup, "w", newline="") as fh:
        csv.writer(fh).writerows(grid)
    log(f"[sheet] snapshot of the whole sheet saved to {backup}")
    note = (f"Total gsdf copies in the diploid genome (2 = no duplicated chromosome, 3 = one, 4 = two), from short-read depth over the "
            f"21.7 kb gsdf unit relative to its flanks; sexloci gsdf-genotype {version}, {datetime.date.today().isoformat()}. "
            f"Blank = no call (no GT3a alignment, low coverage or a flagged call).")
    apply_column_plan(ws, plan, outdir, note, log)
    return plan
