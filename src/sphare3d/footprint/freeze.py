"""The FREEZE GATE between Phase C and Phase D (agent §7.2).

From the moment ``08_FINAL_FOOTPRINT/footprint_freeze.json`` is written,
``FP_original``, ``r_fp``, the descriptor definitions and the QC thresholds are
immutable for the remainder of the run. Phase C is never re-entered: if Phase D
exposes a genuine defect in the Phase C solution, that is an escalation to the
Lead for a re-run under a **new RUN_ID**, never an in-place patch.

The gate is a hash set, not a promise. Phase D verifies it at entry and again at
stage end, and any change is a BLOCKING baseline-integrity failure — because a
robustness number computed against a moved baseline measures nothing.
"""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

from ..utils.errors import BlockedError
from ..utils.hashing import freeze_directory, sha256_file
from ..utils.io import read_json, write_json
from ..utils.runctx import RunContext

FROZEN_STAGES = ("07_FOOTPRINT_RADIUS", "08_FINAL_FOOTPRINT")

#: Written *after* the gate by construction, therefore outside the frozen set.
#: Listed explicitly so the exclusion is auditable rather than implicit.
POST_GATE_FILES = (
    "footprint_freeze.json", "handoff_03.json", "stage_status.json",
    "stage_manifest.tsv", "stage_c_report.md", "NOT_RUN.txt",
    *(f"warnings_{stage}.tsv" for stage in FROZEN_STAGES),
)


def _frozen_map(ctx: RunContext) -> dict[str, str]:
    out: dict[str, str] = {}
    for stage in FROZEN_STAGES:
        stage_dir = ctx.full_results / stage
        for rel, digest in freeze_directory(stage_dir, ctx.run_root).items():
            if Path(rel).name in POST_GATE_FILES:
                continue
            out[rel] = digest
    return dict(sorted(out.items()))


def write_freeze(ctx: RunContext, r_fp: float | None, code_version: str,
                 api_version: str, extra: dict | None = None) -> dict:
    """Seal Phase C. Returns the payload that Phase D re-verifies."""
    payload = {
        "run_id": ctx.run_id,
        "config_sha256": ctx.config.sha256,
        "frozen_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "footprint_radius": r_fp,
        "footprint_code_version": code_version,
        "footprint_api_version": api_version,
        "frozen_stages": list(FROZEN_STAGES),
        "excluded_written_after_gate": list(POST_GATE_FILES),
        "immutability_statement": (
            "FP_original, r_fp, the descriptor definitions and the QC thresholds are "
            "immutable for the remainder of this run. Phase C is never re-entered; "
            "no robustness observation may cause a re-selection, re-tuning or "
            "re-export of r_fp or FP_original."),
        "files_sha256": _frozen_map(ctx),
    }
    if extra:
        payload.update(extra)
    return {"payload": payload,
            "path": write_json(
                ctx.full_results / "08_FINAL_FOOTPRINT" / "footprint_freeze.json",
                payload)}


def load_freeze(ctx: RunContext) -> dict:
    path = ctx.full_results / "08_FINAL_FOOTPRINT" / "footprint_freeze.json"
    if not path.is_file():
        raise BlockedError(
            "BLOCKED — footprint_freeze.json is absent. Phase D may not start "
            "without the Phase C freeze record (agent §6.8)."
        )
    return read_json(path)


def verify_freeze(ctx: RunContext, payload: dict, when: str) -> None:
    """Re-verify every frozen hash. Any change is BLOCKING."""
    recorded = payload.get("files_sha256", {})
    if not recorded:
        raise BlockedError("BLOCKED — the freeze record contains no file hashes.")
    problems = []
    for rel, expected in sorted(recorded.items()):
        target = ctx.run_root / rel
        if not target.is_file():
            problems.append(f"{rel}: MISSING")
        elif sha256_file(target) != expected:
            problems.append(f"{rel}: HASH_MISMATCH")
    if problems:
        raise BlockedError(
            f"BASELINE INTEGRITY FAILURE ({when}) — frozen Phase C artifacts changed: "
            f"{problems[:5]}. FP_original and r_fp are hash-frozen before Phase D "
            f"starts; a robustness result computed against a moved baseline is void. "
            f"Stopping and reporting to the Lead."
        )


def verify_upstream_unchanged(hashes: dict[str, str], paths: dict[str, Path],
                              when: str) -> None:
    """§9.20 — upstream artifacts must be bit-identical at stage end."""
    problems = []
    for name, expected in sorted(hashes.items()):
        path = paths[name]
        if not path.is_file():
            problems.append(f"{name}: MISSING")
        elif sha256_file(path) != expected:
            problems.append(f"{name}: HASH_MISMATCH")
    if problems:
        raise BlockedError(
            f"BASELINE INTEGRITY FAILURE ({when}) — upstream artifacts changed "
            f"during the stage: {problems}. Stopping immediately and reporting."
        )
