"""The deterministic stage-contract validator — NOT an LLM agent.

Runs after any stage completes, mechanically, in milliseconds, every time. It
checks the SHAPE of what a stage produced against the Output Contract and the
handoff schemas already declared in :mod:`hotspot3d.orchestration.contracts` and
:mod:`hotspot3d.reporting.manifest` — it does not recompute, re-derive, or judge
any scientific value. A radius of 21 Å and zero significant centers both pass
this validator if the artifacts that report them are well-formed; whether either
number is the CORRECT answer is not this module's question.

This is the always-on half of the QA architecture. The failure-triggered half —
:mod:`hotspot3d.utils.qa_taxonomy` and the ``pipeline-qa-debugger`` agent it
backs — is invoked only when this validator (or a scientific agent's own anomaly
flag) reports something wrong. Nothing here spawns an agent, makes a judgment
call, or is slow enough to matter next to the stage it validates.

What is checked (Section 6 of the QA architecture spec):
  * the stage produced its declared PRIMARY outputs, or is explicitly NOT_RUN
    with a NOT_RUN.txt explaining why (P3 — absence is explicit, never silent);
  * every artifact that exists is well-formed: readable, parses as declared
    (JSON/TSV), and is not a zero-byte or truncated file masquerading as content;
  * the handoff (if the stage produces one) carries every required payload key
    (:data:`HANDOFF_REQUIRED_KEYS`) and its manifest hashes verify against disk;
  * ``run_id`` and ``config_sha256`` are identical across every handoff a stage
    both reads and writes — a run must never silently mix two configurations;
  * ``status``/``outcome_type`` are members of the frozen vocabulary, and a
    non-COMPLETED stage carries a NOT_RUN.txt;
  * declared counts that MUST reconcile do (e.g. a written manifest's row count
    against the files actually on disk).

What is deliberately NOT checked here: whether a p-value is small enough,
whether a radius is "the right" radius, whether a hotspot is biologically
plausible — anything requiring scientific judgment routes to a human or, for
implementation-shaped anomalies, to ``pipeline-qa-debugger``.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .hashing import sha256_file
from .runctx import RunContext
from .status import OutcomeType, Status

VALID_STATUS_VALUES = {s.value for s in Status}
VALID_OUTCOME_VALUES = {o.value for o in OutcomeType}


@dataclass
class ContractViolation:
    """One mechanical failure. ``rule`` is a stable short id (grep-able in a bug
    brief's EVIDENCE_ARTIFACT); ``detail`` is the specific, concrete finding."""

    rule: str
    detail: str
    path: str | None = None

    def __str__(self) -> str:
        where = f" ({self.path})" if self.path else ""
        return f"[{self.rule}]{where} {self.detail}"


@dataclass
class ContractResult:
    stage: str
    violations: list[ContractViolation] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not self.violations

    def add(self, rule: str, detail: str, path: str | None = None) -> None:
        self.violations.append(ContractViolation(rule=rule, detail=detail, path=path))

    def summary(self) -> str:
        if self.passed:
            return f"{self.stage}: PASS (stage-contract validator)"
        lines = [f"{self.stage}: FAIL ({len(self.violations)} violation(s))"]
        lines += [f"  - {v}" for v in self.violations]
        return "\n".join(lines)

    def raise_if_failed(self) -> None:
        if not self.passed:
            from .errors import BlockedError
            raise BlockedError(
                f"BLOCKED — stage-contract validation failed for {self.stage}:\n"
                + self.summary())


def _try_read_json(path: Path, result: ContractResult, rule: str) -> dict | None:
    try:
        raw = path.read_text(encoding="utf-8")
    except OSError as exc:
        result.add(rule, f"could not read: {exc}", str(path))
        return None
    if not raw.strip():
        result.add(rule, "file exists but is empty (0 meaningful bytes)", str(path))
        return None
    try:
        return json.loads(raw)
    except json.JSONDecodeError as exc:
        result.add(rule, f"malformed JSON, not silently accepted: {exc}", str(path))
        return None


def _expected_for(stage: str) -> list[tuple[str, str, str]]:
    """PRIMARY-tier expected files for ``stage``, from the ONE canonical source
    (:data:`hotspot3d.reporting.manifest.EXPECTED`) rather than a second list
    this module would have to keep in sync by hand."""
    from ..reporting.manifest import EXPECTED
    return [row for row in EXPECTED.get(stage, []) if row[1] == "PRIMARY"]


def validate_stage_contract(ctx: RunContext, stage: str) -> ContractResult:
    """Run every mechanical check for one stage directory. Always returns a
    :class:`ContractResult` — never raises on its own; call
    :meth:`ContractResult.raise_if_failed` if the caller wants BLOCKED semantics.
    """
    result = ContractResult(stage=stage)
    stage_dir = ctx.full_results / stage

    if not stage_dir.is_dir():
        result.add("STAGE_DIR_MISSING", "stage directory does not exist at all")
        return result

    # -- stage_status.json: must exist, parse, and carry a valid vocabulary -----
    status_path = stage_dir / "stage_status.json"
    status_payload = None
    if not status_path.is_file():
        result.add("STAGE_STATUS_MISSING", "stage_status.json is required for "
                   "every stage directory without exception (P3)", str(status_path))
    else:
        status_payload = _try_read_json(status_path, result, "STAGE_STATUS_MALFORMED")

    status_value = None
    outcome_value = None
    if status_payload is not None:
        status_value = status_payload.get("status")
        outcome_value = status_payload.get("outcome_type")
        if status_value not in VALID_STATUS_VALUES:
            result.add("STAGE_STATUS_INVALID_VOCABULARY",
                      f"status={status_value!r} is not one of {sorted(VALID_STATUS_VALUES)}",
                      str(status_path))
        if outcome_value not in VALID_OUTCOME_VALUES:
            result.add("STAGE_STATUS_INVALID_VOCABULARY",
                      f"outcome_type={outcome_value!r} is not one of "
                      f"{sorted(VALID_OUTCOME_VALUES)}", str(status_path))

        # -- terminal-state propagation: non-COMPLETED must explain itself -------
        if status_value not in (None, Status.COMPLETED.value):
            not_run = stage_dir / "NOT_RUN.txt"
            if not not_run.is_file():
                result.add("MISSING_NOT_RUN_EXPLANATION",
                          f"status={status_value} but no NOT_RUN.txt exists — P3 "
                          f"requires absence to be explicit, never silent",
                          str(not_run))

    # -- PRIMARY outputs: present when COMPLETED-like, absent+explained otherwise
    expected = _expected_for(stage)
    is_completed_like = status_value in (Status.COMPLETED.value,
                                         Status.COMPLETED_NEGATIVE.value)
    for name, _tier, _desc in expected:
        path = stage_dir / name
        exists = path.is_file()
        if is_completed_like and not exists:
            result.add("EXPECTED_OUTPUT_MISSING",
                      f"stage status is {status_value} but the declared PRIMARY "
                      f"output {name!r} does not exist", str(path))
        if exists and path.stat().st_size == 0:
            result.add("EXPECTED_OUTPUT_EMPTY",
                      f"{name!r} exists but is zero bytes — a malformed/partial "
                      f"artifact must never be silently accepted as valid", str(path))
        if exists and name.endswith(".json"):
            _try_read_json(path, result, "EXPECTED_OUTPUT_MALFORMED")

    return result


def validate_handoff_contract(ctx: RunContext, handoff, *, upstream_handoffs:
                              list[Any] | None = None) -> ContractResult:
    """Handoff-specific checks: required keys, manifest integrity, and
    run_id/config_sha256 consistency against any upstream handoffs this stage
    itself consumed (a run must never silently mix two configurations).
    """
    from ..orchestration.contracts import HANDOFF_REQUIRED_KEYS
    from .hashing import verify_manifest

    result = ContractResult(stage=handoff.name)

    missing = [k for k in HANDOFF_REQUIRED_KEYS.get(handoff.name, [])
              if k not in handoff.payload]
    if missing:
        result.add("HANDOFF_SCHEMA_INCOMPLETE",
                  f"missing required payload keys: {missing}")

    if handoff.manifest:
        bad = verify_manifest(handoff.manifest, ctx.run_root)
        if bad:
            result.add("HANDOFF_MANIFEST_HASH_MISMATCH",
                      f"{len(bad)} file(s) failed hash verification: {bad[:5]}"
                      + (" ..." if len(bad) > 5 else ""))
    else:
        result.add("HANDOFF_MANIFEST_EMPTY",
                  "handoff carries no manifest entries — nothing to verify against")

    for upstream in (upstream_handoffs or []):
        if upstream.run_id != handoff.run_id:
            result.add("RUN_ID_MISMATCH",
                      f"{handoff.name} carries run_id={handoff.run_id!r} but its "
                      f"upstream {upstream.name} carries {upstream.run_id!r}")
        if upstream.config_sha256 != handoff.config_sha256:
            result.add("CONFIG_SHA256_MISMATCH",
                      f"{handoff.name} carries config_sha256="
                      f"{handoff.config_sha256[:12]}... but its upstream "
                      f"{upstream.name} carries {upstream.config_sha256[:12]}...")

    return result


def validate_manifest_row_count(ctx: RunContext, manifest_rows: list[dict]) -> ContractResult:
    """The run manifest's row count for CREATED files must equal the number of
    files it claims exist — catching a manifest that silently drifted from disk."""
    result = ContractResult(stage="00_RUN_SUMMARY")
    created = [r for r in manifest_rows if r.get("status") == "CREATED"]
    missing_on_disk = [r["relative_path"] for r in created
                       if not (ctx.run_root / r["relative_path"]).is_file()]
    if missing_on_disk:
        result.add("MANIFEST_CLAIMS_FILE_NOT_ON_DISK",
                  f"{len(missing_on_disk)} row(s) marked CREATED but absent from "
                  f"disk: {missing_on_disk[:5]}" + (" ..." if len(missing_on_disk) > 5 else ""))
    return result
