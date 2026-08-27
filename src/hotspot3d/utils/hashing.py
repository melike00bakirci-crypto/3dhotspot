"""SHA-256 hashing and code versioning.

Every handoff carries a manifest of hashes; every consumer recomputes them
(trust is never transitive, Part IV).
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

_CHUNK = 1 << 20


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_text(text: str) -> str:
    return sha256_bytes(text.encode("utf-8"))


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(_CHUNK):
            h.update(chunk)
    return h.hexdigest()


def manifest_for(paths: list[str | Path], root: str | Path | None = None) -> dict[str, str]:
    """Map relative path -> sha256 for every existing file in ``paths``."""
    root = Path(root) if root is not None else None
    out: dict[str, str] = {}
    for p in sorted(str(x) for x in paths):
        pp = Path(p)
        if not pp.is_file():
            continue
        key = str(pp.relative_to(root)) if root is not None else str(pp)
        out[key] = sha256_file(pp)
    return out


def verify_manifest(manifest: dict[str, str], root: str | Path) -> list[str]:
    """Return the list of mismatching / missing entries. Empty list == verified."""
    root = Path(root)
    bad: list[str] = []
    for rel, expected in sorted(manifest.items()):
        target = root / rel
        if not target.is_file():
            bad.append(f"{rel}: MISSING")
        elif sha256_file(target) != expected:
            bad.append(f"{rel}: HASH_MISMATCH")
    return bad


def code_version(package_root: str | Path | None = None) -> str:
    """Deterministic content hash of the pipeline source tree.

    Used to pin the frozen APIs that later stages re-execute, and as a component
    of RUN_ID. Depends only on file contents and relative paths, so it is stable
    across machines and checkouts.
    """
    root = Path(package_root) if package_root else Path(__file__).resolve().parents[1]
    h = hashlib.sha256()
    for path in sorted(root.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        h.update(str(path.relative_to(root)).encode("utf-8"))
        h.update(path.read_bytes())
    return h.hexdigest()


def dir_files(root: str | Path) -> list[Path]:
    """All files under ``root``, sorted, excluding caches."""
    root = Path(root)
    if not root.exists():
        return []
    return sorted(
        p for p in root.rglob("*")
        if p.is_file() and "__pycache__" not in p.parts and not p.name.endswith(".pyc")
    )


def freeze_directory(root: str | Path, relative_to: str | Path | None = None) -> dict[str, str]:
    """Hash every file under ``root`` — the FREEZE GATE primitive (agent 3 §7.2)."""
    base = Path(relative_to) if relative_to is not None else Path(root)
    return {str(p.relative_to(base)): sha256_file(p) for p in dir_files(root)}


def env_fingerprint() -> dict[str, str]:
    import platform
    import sys
    return {
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor() or "unknown",
        "cpu_count": str(os.cpu_count()),
    }
