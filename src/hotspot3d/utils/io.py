"""Canonical TSV / JSON writers — Output Contract IX.4 (FROZEN conventions).

All TSVs: UTF-8, LF, tab-separated, header row mandatory, ``NA`` for missing
(never empty, never NaN), booleans ``TRUE``/``FALSE``, floats at 6 significant
digits unless stated, frozen column order, ``schema_version`` column.
All JSON: ``schema_version`` key.

Windows-safe naming (P7) is enforced at write time, not discovered at download
time: filenames match ``[A-Za-z0-9._-]+``, relative paths stay <= 150 chars.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Iterable, Sequence

from .errors import BlockedError

SCHEMA_VERSION = "1.0.0"
NA = "NA"
_FILENAME_RE = re.compile(r"^[A-Za-z0-9._-]+$")
MAX_REL_PATH = 150


# --- naming / safety --------------------------------------------------------

def check_windows_safe(path: str | Path, run_root: str | Path | None = None) -> None:
    p = Path(path)
    for part in p.parts:
        if part in ("/", "\\"):
            continue
        if not _FILENAME_RE.match(part):
            raise BlockedError(
                f"Output Contract P7 violation: path component {part!r} in {p} is not "
                f"Windows-safe (must match [A-Za-z0-9._-]+)."
            )
    if run_root is not None:
        rel = str(p.relative_to(run_root)) if str(p).startswith(str(run_root)) else str(p)
        if len(rel) > MAX_REL_PATH:
            raise BlockedError(
                f"Output Contract P7 violation: relative path is {len(rel)} chars "
                f"(limit {MAX_REL_PATH}): {rel}"
            )


# --- value formatting -------------------------------------------------------

def fmt_value(value: Any, sig: int = 6) -> str:
    """Format one cell under the frozen TSV conventions."""
    if value is None:
        return NA
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, str):
        if value == "":
            return NA
        return value.replace("\t", " ").replace("\n", " ").replace("\r", " ")
    if isinstance(value, (int,)):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return NA
        if value == 0:
            return "0"
        return f"{value:.{sig}g}"
    if isinstance(value, (list, tuple, set)):
        if not value:
            return NA
        return ";".join(fmt_value(v, sig) for v in value)
    return str(value)


def parse_value(text: str) -> Any:
    if text == NA:
        return None
    if text == "TRUE":
        return True
    if text == "FALSE":
        return False
    try:
        if re.fullmatch(r"-?\d+", text):
            return int(text)
        return float(text)
    except ValueError:
        return text


# --- TSV --------------------------------------------------------------------

def write_tsv(path: str | Path, rows: Iterable[dict], columns: Sequence[str],
              schema_version: str = SCHEMA_VERSION, sig: int = 6,
              run_root: str | Path | None = None) -> Path:
    """Write a canonical TSV with a frozen column order.

    ``schema_version`` is prepended automatically if not already in ``columns``.
    An empty ``rows`` still writes the header — an empty table is explicit, never
    a missing file (P3).
    """
    path = Path(path)
    check_windows_safe(path, run_root)
    path.parent.mkdir(parents=True, exist_ok=True)

    cols = list(columns)
    if "schema_version" not in cols:
        cols = ["schema_version"] + cols

    lines = ["\t".join(cols)]
    for row in rows:
        record = dict(row)
        record.setdefault("schema_version", schema_version)
        unknown = set(record) - set(cols)
        if unknown:
            raise BlockedError(
                f"Schema violation writing {path.name}: unexpected column(s) "
                f"{sorted(unknown)}; frozen schema is {cols}"
            )
        lines.append("\t".join(fmt_value(record.get(c), sig) for c in cols))

    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def read_tsv(path: str | Path, allowed_columns: Sequence[str] | None = None,
             forbidden_columns: Sequence[str] | None = None) -> list[dict]:
    """Read a canonical TSV.

    ``allowed_columns`` implements the Stage B **column allowlist**: loading any
    column outside the allowlist aborts as a leakage event. ``forbidden_columns``
    is the explicit complement used by handoff enforcement.
    """
    from .errors import LeakageError

    path = Path(path)
    if not path.is_file():
        raise BlockedError(f"BLOCKED — required input not found: {path}")
    raw = path.read_text(encoding="utf-8").splitlines()
    if not raw:
        raise BlockedError(f"BLOCKED — empty file (no header): {path}")
    header = raw[0].split("\t")

    if forbidden_columns:
        hit = sorted(set(header) & set(forbidden_columns))
        if hit:
            raise LeakageError(
                f"LEAKAGE EVENT reading {path.name}: forbidden column(s) {hit} present. "
                f"Stage aborted."
            )
    if allowed_columns is not None:
        extra = sorted(set(header) - set(allowed_columns) - {"schema_version"})
        if extra:
            raise LeakageError(
                f"LEAKAGE EVENT reading {path.name}: column(s) {extra} are outside the "
                f"permitted allowlist {sorted(allowed_columns)}. Stage aborted."
            )

    out = []
    for line in raw[1:]:
        if not line:
            continue
        values = line.split("\t")
        out.append({k: parse_value(v) for k, v in zip(header, values)})
    return out


def read_tsv_columns(path: str | Path, columns: Sequence[str]) -> list[dict]:
    """Project a TSV down to an explicit column subset (allowlist-safe read)."""
    rows = read_tsv(path)
    return [{c: r.get(c) for c in columns} for r in rows]


# --- JSON -------------------------------------------------------------------

def write_json(path: str | Path, payload: dict, schema_version: str = SCHEMA_VERSION,
               run_root: str | Path | None = None) -> Path:
    path = Path(path)
    check_windows_safe(path, run_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    body = {"schema_version": schema_version, **payload}
    text = json.dumps(body, indent=2, sort_keys=False, default=_json_default,
                      ensure_ascii=False)
    path.write_text(text + "\n", encoding="utf-8", newline="\n")
    return path


def read_json(path: str | Path) -> dict:
    path = Path(path)
    if not path.is_file():
        raise BlockedError(f"BLOCKED — required input not found: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _json_default(obj: Any) -> Any:
    import numpy as np
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        value = float(obj)
        return None if (math.isnan(value) or math.isinf(value)) else value
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (set, frozenset)):
        return sorted(obj)
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"not JSON serializable: {type(obj)}")


def write_text(path: str | Path, text: str, crlf: bool = False,
               run_root: str | Path | None = None) -> Path:
    """Plain-text writer. ``crlf=True`` only for README.txt / transfer_info.txt (IX.0)."""
    path = Path(path)
    check_windows_safe(path, run_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    newline = "\r\n" if crlf else "\n"
    path.write_text(text if text.endswith("\n") else text + "\n",
                    encoding="utf-8", newline=newline)
    return path
