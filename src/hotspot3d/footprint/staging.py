"""Per-stage bookkeeping shared by 07_*, 08_* and 09_* (Output Contract IX.3).

Absence is explicit (P3): every expected file appears in ``stage_manifest.tsv``,
created or not, so a missing artifact is a visible ``NOT_CREATED`` row rather than
something a reader has to notice is absent.
"""
from __future__ import annotations

from pathlib import Path

from ..utils.hashing import sha256_file
from ..utils.io import write_tsv

MANIFEST_NAME = "stage_manifest.tsv"

STAGE_MANIFEST_COLUMNS = [
    "file", "expected", "created", "status", "sha256", "size_bytes", "description",
]


def write_stage_manifest(stage_dir: Path, expected: list[tuple[str, str]],
                         extra_created: list[Path] | None = None) -> Path:
    """``expected`` is a list of ``(relative_path, description)`` pairs."""
    rows = []
    declared = set()
    for rel, description in expected:
        declared.add(rel)
        path = stage_dir / rel
        if rel == MANIFEST_NAME:
            # The manifest cannot hash itself; the row is still present so the
            # expected-file list is complete rather than quietly one short.
            rows.append({
                "file": rel, "expected": True, "created": True,
                "status": "CREATED", "sha256": None, "size_bytes": None,
                "description": f"{description} (hash omitted: self-reference)",
            })
            continue
        created = path.is_file()
        rows.append({
            "file": rel,
            "expected": True,
            "created": created,
            "status": "CREATED" if created else "NOT_CREATED",
            "sha256": sha256_file(path) if created else None,
            "size_bytes": path.stat().st_size if created else None,
            "description": description,
        })
    for path in sorted(extra_created or []):
        rel = str(path.relative_to(stage_dir))
        if rel in declared:
            continue
        rows.append({
            "file": rel, "expected": False, "created": True, "status": "CREATED",
            "sha256": sha256_file(path), "size_bytes": path.stat().st_size,
            "description": "additional retained intermediate",
        })
    return write_tsv(stage_dir / MANIFEST_NAME, rows, STAGE_MANIFEST_COLUMNS)


def stage_files(stage_dir: Path) -> list[Path]:
    return sorted(p for p in stage_dir.rglob("*")
                  if p.is_file() and "__pycache__" not in p.parts)
