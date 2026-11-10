"""Structured warnings. BLOCK stops the run, WARN is reported (and may cap confidence), INFO is logged."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

BLOCK, WARN, INFO = "BLOCK", "WARN", "INFO"
_ORDER = {BLOCK: 0, WARN: 1, INFO: 2}


@dataclass
class Finding:
    severity: str
    code: str
    message: str
    samples: List[str] = field(default_factory=list)
    cap: Optional[str] = None      # best tier still allowed ('B' or 'C') when this finding applies
    why: str = ""                  # why this matters / why it blocks
    fix: str = ""                  # what to do about it
    rows: List[dict] = field(default_factory=list)   # per-sample detail for fix_list.tsv


class Findings(list):
    def add(self, severity, code, message, samples=None, cap=None, why="", fix="", rows=None):
        self.append(Finding(severity, code, message, list(samples or []), cap, why, fix, list(rows or [])))

    def block(self, code, message, samples=None, **kw):
        self.add(BLOCK, code, message, samples, **kw)

    def warn(self, code, message, samples=None, cap=None, **kw):
        self.add(WARN, code, message, samples, cap, **kw)

    def info(self, code, message, samples=None, **kw):
        self.add(INFO, code, message, samples, **kw)

    @property
    def blocked(self) -> bool:
        return any(f.severity == BLOCK for f in self)

    def acknowledge(self, codes):
        """Downgrade reviewed WARNs to INFO. Never touches BLOCKs; a finding's tier cap is kept."""
        for f in self:
            if f.severity == WARN and f.code in set(codes):
                f.severity = INFO
                f.message = "[acknowledged] " + f.message

    def sorted(self):
        return sorted(self, key=lambda f: (_ORDER[f.severity], f.code))

    def tier_cap(self) -> Optional[str]:
        caps = [f.cap for f in self if f.cap]
        return max(caps) if caps else None


def short_list(items, n=8) -> str:
    items = list(items)
    s = ", ".join(map(str, items[:n]))
    return s + (f", ... (+{len(items) - n} more)" if len(items) > n else "")


_ID = __import__("re").compile(r"^(.*?)(\d+)$")


def compress_ids(ids, max_items: int = 14) -> str:
    """'MC_046, MC_047, MC_049 ... MC_094' -> 'MC_046, MC_047, MC_049-MC_094 (46)'. Runs of >= 3 consecutive numbers
    sharing a prefix and zero-padding are collapsed; IDs without a numeric tail are listed as they are."""
    ids = list(dict.fromkeys(map(str, ids)))
    groups, loose = {}, []
    for x in ids:
        m = _ID.match(x)
        if m:
            groups.setdefault(m.group(1), []).append((int(m.group(2)), x))
        else:
            loose.append(x)
    parts = []
    for _, items in sorted(groups.items(), key=lambda kv: kv[0]):
        items.sort()
        run = [items[0]]
        for it in items[1:] + [None]:
            if it is not None and it[0] == run[-1][0] + 1:
                run.append(it)
                continue
            parts.append(f"{run[0][1]}-{run[-1][1]} ({len(run)})" if len(run) >= 3 else ", ".join(r[1] for r in run))
            run = [it] if it is not None else []
    parts += sorted(loose)
    if len(parts) > max_items:
        return ", ".join(parts[:max_items]) + f", ... (+{len(parts) - max_items} more groups; {len(ids)} samples in total)"
    return ", ".join(parts)
