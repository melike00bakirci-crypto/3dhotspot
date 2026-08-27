"""Stage status, warnings and the frozen warning-code enum — Output Contract IX.3/IX.7.

Principles enforced:
  P3  Absence is explicit — every stage ALWAYS emits stage_status.json.
  P4  Scientific negative != technical failure — ``outcome_type`` separates them.
  IX.7 Severity names the required action: INFO / ADVISORY / MAJOR / BLOCKING.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Iterable

from .io import write_json, write_tsv


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class Severity(str, Enum):
    INFO = "INFO"            # recorded for provenance; no action
    ADVISORY = "ADVISORY"    # worth knowing; does not affect validity
    MAJOR = "MAJOR"          # may materially change interpretation; must be read
    BLOCKING = "BLOCKING"    # a stage or the run stopped; result unusable as-is


class Status(str, Enum):
    COMPLETED = "COMPLETED"
    COMPLETED_NEGATIVE = "COMPLETED_NEGATIVE"
    #: Workflow v2 §0/§5.4 — the stage ran and rejected nothing, but the design
    #: could not have rejected anything. Distinct from COMPLETED_NEGATIVE, which
    #: asserts the opposite about the design, and from a technical failure.
    UNDERPOWERED = "UNDERPOWERED"
    NOT_RUN = "NOT_RUN"
    BLOCKED = "BLOCKED"
    FAILED = "FAILED"


class OutcomeType(str, Enum):
    COMPLETED = "COMPLETED"
    SCIENTIFIC_NEGATIVE = "SCIENTIFIC_NEGATIVE"
    #: Workflow v2 Appendix A — the run is UNINFORMATIVE about the gene. It is
    #: neither a finding nor a fault, and must never be reported as either.
    UNINFORMATIVE = "UNINFORMATIVE"
    TECHNICAL_FAILURE = "TECHNICAL_FAILURE"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class RunStatus(str, Enum):
    COMPLETED = "COMPLETED"
    COMPLETED_WITH_WARNINGS = "COMPLETED_WITH_WARNINGS"
    COMPLETED_NEGATIVE = "COMPLETED_NEGATIVE"
    UNDERPOWERED = "UNDERPOWERED"
    PARTIAL = "PARTIAL"
    FAILED = "FAILED"
    BLOCKED = "BLOCKED"


BANNERS = {
    RunStatus.COMPLETED: "ANALYSIS COMPLETE: {gene}",
    RunStatus.COMPLETED_WITH_WARNINGS: "ANALYSIS COMPLETE (WITH MAJOR WARNINGS): {gene}",
    RunStatus.COMPLETED_NEGATIVE: "ANALYSIS COMPLETE — NEGATIVE RESULT: {gene}",
    # Deliberately not worded as a result. Workflow v2 Appendix A forbids
    # reporting this as an absence of hotspots.
    RunStatus.UNDERPOWERED: "ANALYSIS UNINFORMATIVE — TEST CANNOT REJECT: {gene}",
    RunStatus.PARTIAL: "ANALYSIS PARTIAL: {gene}",
    RunStatus.FAILED: "ANALYSIS FAILED: {gene}",
    RunStatus.BLOCKED: "RUN BLOCKED: {gene}",
}

EXIT_CODES = {
    RunStatus.COMPLETED: 0,
    RunStatus.COMPLETED_WITH_WARNINGS: 0,
    RunStatus.COMPLETED_NEGATIVE: 0,
    # Non-zero: an underpowered run produced no interpretable result, so a caller
    # that checks only the exit code must not mistake it for a completed analysis.
    RunStatus.UNDERPOWERED: 3,
    RunStatus.PARTIAL: 0,
    RunStatus.FAILED: 1,
    RunStatus.BLOCKED: 2,
}


# --- FROZEN warning-code enum (IX.7) ---------------------------------------
# Extensible only by pre-registered revision. Default severity is fixed here so
# an agent cannot quietly downgrade a MAJOR finding.
WARNING_CODES: dict[str, Severity] = {
    "PERMUTATION_RESOLUTION_LIMITED": Severity.MAJOR,
    "FALLBACK_RADIUS_DOMAIN": Severity.MAJOR,
    "BOUNDARY_OPTIMUM_WARNING": Severity.MAJOR,
    "EXCESSIVE_FOOTPRINT_COVERAGE": Severity.ADVISORY,
    "FOOTPRINT_COVERAGE_CONSTRAINT_FAILURE": Severity.BLOCKING,
    "GLOBAL_CLUSTERING_NON_SIGNIFICANT": Severity.MAJOR,
    "NO_SIGNIFICANT_HOTSPOTS": Severity.MAJOR,
    "RESIDUE_CLASS_CONFLICT": Severity.ADVISORY,
    "STRUCTURAL_CONFIDENCE_SENSITIVE": Severity.MAJOR,
    "REVIEW_STATUS_SENSITIVE": Severity.MAJOR,
    "LOO_SPARSE_EVIDENCE": Severity.ADVISORY,
    "INSUFFICIENT_CLASSIFIED_RESIDUES": Severity.BLOCKING,
    "STRUCTURE_MAPPING_FAILURE": Severity.BLOCKING,
    "MULTI_FRAGMENT_AFDB_ENTRY": Severity.BLOCKING,
    "RENDERING_UNAVAILABLE": Severity.INFO,
    "PDF_UNAVAILABLE": Severity.INFO,
    "REVIEW_PACK_SIZE_EXCEEDED": Severity.ADVISORY,
    "NEAR_TIE_RADIUS_SELECTION": Severity.ADVISORY,
    "OBJECTIVE_REDUNDANCY": Severity.ADVISORY,
    # additional pre-registered operational codes
    # Instability of the geometric footprint under center perturbation. A dedicated
    # code is required because agent 3 §15 obliges the stage to surface instability
    # EXPLICITLY, and reusing an unrelated code (e.g. the Stage B LOO diagnostic)
    # would mislabel the finding rather than report it.
    "LIMITED_GEOMETRIC_ROBUSTNESS": Severity.MAJOR,
    # Agent 3 §15 lists "iterations failing en masse" as an escalation condition, so it
    # must reach the run-level warnings and the terminal summary — a JSON block alone
    # would leave an escalation invisible at the two places a reader actually looks.
    "ROBUSTNESS_ITERATIONS_FAILED_EN_MASSE": Severity.MAJOR,
    # n_S = 2 gives K_MAX = 1 and exactly two iterations: the profile is valid but rests
    # on two observations. II.11 requires that be "reported honestly as minimal".
    "ROBUSTNESS_EVIDENCE_MINIMAL": Severity.ADVISORY,
    "ROBUSTNESS_NOT_EVALUABLE": Severity.MAJOR,
    "UNDERPOWERED_COHORT": Severity.MAJOR,
    # --- Workflow v2 Step 1 (§4, §5.4) -------------------------------------
    # The per-center test could not have rejected any center regardless of the
    # data: p_floor > q/m. BLOCKING because Appendix A forbids reporting such a
    # run as a negative, and a lesser severity would let it pass as one.
    "TEST_CANNOT_REJECT": Severity.BLOCKING,
    # §4 — the Pareto set had one member, or every objective was degenerate, so
    # the radius was not selected by optimization; it was the only survivor.
    "VACUOUS_PARETO_SELECTION": Severity.MAJOR,
    # §4 — more than half the scanned radii were ruled inadmissible, so
    # admissibility, not the objectives, decided the outcome.
    "ADMISSIBILITY_DOMINATES_SELECTION": Severity.ADVISORY,
    # §4 — an objective identical at every admissible radius contributes nothing
    # to distance-to-utopia and is excluded from it; reported, never silent.
    "DEGENERATE_OBJECTIVE": Severity.ADVISORY,
    # --- Workflow v2 remaining migration -----------------------------------
    # §2 — the biological assembly is an obligate oligomer and no multimer model
    # could be used; the run is stopped rather than analysing inter-subunit
    # neighbourhoods no monomer model contains.
    "OLIGOMERIC_ASSEMBLY_UNAVAILABLE": Severity.BLOCKING,
    # §2 — the exception path config/pipeline.yaml documents: proceeding on the
    # monomer with the limitation declared and hotspots labelled intra-subunit
    # only. Not a fault; a required disclosure whenever it is taken.
    "OLIGOMERIC_ASSEMBLY_MONOMER_DECLARED": Severity.ADVISORY,
    # §4 — a selected radius whose sphere covers most of U_struct localizes
    # nothing; see config biological_plausibility.structural_coverage_*.
    "EXCESSIVE_STRUCTURAL_COVERAGE": Severity.MAJOR,
    # §4 — the single largest sphere at the selected radius absorbs most of the
    # classified cohort; a soft, ADVISORY-only penalty (no MAJOR bound in v2).
    "HIGH_COHORT_ABSORPTION": Severity.ADVISORY,
    # §4 — the selected radius is large relative to the modelled structure's own
    # radius of gyration.
    "HIGH_R_TO_DOMAIN_RATIO": Severity.ADVISORY,
    # §8 — a robustness iteration is recorded as underpowered rather than as a
    # footprint that failed to reproduce, so pooling the two never happens.
    "ROBUSTNESS_ITERATION_UNDERPOWERED": Severity.ADVISORY,
    "BH_BY_DISAGREEMENT": Severity.MAJOR,
    "VOXEL_RESOLUTION_FALLBACK": Severity.ADVISORY,
    "ANNOTATION_SOURCE_UNVERSIONED": Severity.BLOCKING,
    "POSTHOC_GATE_FAILED": Severity.INFO,
    "SYNTHETIC_INPUT_MODE": Severity.INFO,
    # LiveAnnotationSource (annotation/live_sources.py) implements five of the
    # six AnnotationSource methods; literature_mechanisms is explicitly out of
    # scope and must never be indistinguishable from a real, searched, empty
    # result (agent §9.7's "found nothing" is a different, positive claim).
    "LITERATURE_MECHANISMS_NOT_RUN": Severity.MAJOR,
}

WARNING_COLUMNS = [
    "severity", "warning_code", "stage", "agent_owner", "message",
    "potential_consequence", "recommended_action", "affected_output",
    "status", "detected_utc",
]

STAGE_STATUS_SCHEMA = "1.0.0"


@dataclass
class Warning_:
    """One row of warnings_<stage>.tsv."""
    warning_code: str
    stage: str
    agent_owner: str
    message: str
    potential_consequence: str = ""
    recommended_action: str = ""
    affected_output: str = ""
    status: str = "UNRESOLVED"        # UNRESOLVED | RESOLVED | ACCEPTED_BY_DESIGN
    severity: Severity | None = None
    detected_utc: str = field(default_factory=utc_now)

    def __post_init__(self):
        if self.warning_code not in WARNING_CODES:
            raise ValueError(
                f"'{self.warning_code}' is not in the FROZEN warning-code enum. "
                f"Codes are extensible only by pre-registered revision."
            )
        if self.severity is None:
            self.severity = WARNING_CODES[self.warning_code]

    def as_row(self) -> dict:
        return {
            "severity": self.severity.value, "warning_code": self.warning_code,
            "stage": self.stage, "agent_owner": self.agent_owner,
            "message": self.message,
            "potential_consequence": self.potential_consequence or "NA",
            "recommended_action": self.recommended_action or "NA",
            "affected_output": self.affected_output or "NA",
            "status": self.status, "detected_utc": self.detected_utc,
        }


_SEVERITY_ORDER = {Severity.BLOCKING: 0, Severity.MAJOR: 1,
                   Severity.ADVISORY: 2, Severity.INFO: 3}


class WarningCollector:
    """Per-stage warning accumulator. The Lead concatenates these (P2)."""

    def __init__(self, stage: str, agent_owner: str):
        self.stage = stage
        self.agent_owner = agent_owner
        self._items: list[Warning_] = []

    def add(self, warning_code: str, message: str, **kwargs) -> Warning_:
        w = Warning_(warning_code=warning_code, stage=self.stage,
                     agent_owner=self.agent_owner, message=message, **kwargs)
        self._items.append(w)
        return w

    def extend(self, warnings: Iterable[Warning_]) -> None:
        self._items.extend(warnings)

    @property
    def items(self) -> list[Warning_]:
        return sorted(self._items, key=lambda w: (_SEVERITY_ORDER[w.severity], w.warning_code))

    def counts(self) -> dict[str, int]:
        out = {s.value: 0 for s in Severity}
        for w in self._items:
            out[w.severity.value] += 1
        return out

    def has(self, code: str) -> bool:
        return any(w.warning_code == code for w in self._items)

    def write(self, stage_dir: str | Path) -> Path:
        return write_tsv(Path(stage_dir) / f"warnings_{self.stage}.tsv",
                         [w.as_row() for w in self.items], WARNING_COLUMNS)


@dataclass
class StageStatus:
    """stage_status.json — written for EVERY stage, always (P3)."""

    stage: str
    agent_owner: str
    status: Status = Status.NOT_RUN
    outcome_type: OutcomeType = OutcomeType.NOT_APPLICABLE
    reason: str = ""
    started_utc: str = ""
    ended_utc: str = ""
    wall_seconds: float = 0.0
    upstream_handoff: str = "NA"
    upstream_qc_status: str = "NA"
    n_outputs_expected: int = 0
    n_outputs_created: int = 0
    warnings: dict = field(default_factory=lambda: {s.value: 0 for s in Severity})
    recommended_action: str = ""
    negative_result: dict | None = None

    def as_dict(self) -> dict:
        payload = {
            "stage": self.stage, "agent_owner": self.agent_owner,
            "status": self.status.value, "outcome_type": self.outcome_type.value,
            "reason": self.reason, "started_utc": self.started_utc,
            "ended_utc": self.ended_utc, "wall_seconds": round(self.wall_seconds, 3),
            "upstream_handoff": self.upstream_handoff,
            "upstream_qc_status": self.upstream_qc_status,
            "n_outputs_expected": self.n_outputs_expected,
            "n_outputs_created": self.n_outputs_created,
            "warnings": self.warnings,
            "recommended_action": self.recommended_action,
        }
        if self.negative_result is not None:
            payload["negative_result"] = self.negative_result
        return payload

    def write(self, stage_dir: str | Path) -> Path:
        stage_dir = Path(stage_dir)
        stage_dir.mkdir(parents=True, exist_ok=True)
        path = write_json(stage_dir / "stage_status.json", self.as_dict(),
                          schema_version=STAGE_STATUS_SCHEMA)
        if self.status != Status.COMPLETED:
            _write_not_run(stage_dir, self)
        return path


def _write_not_run(stage_dir: Path, st: StageStatus) -> None:
    """NOT_RUN.txt — the human paragraph required whenever status != COMPLETED (IX.3)."""
    # A scientific negative has TWO distinct forms and they must never share wording:
    #   (A) the analysis RAN and found nothing  -> status COMPLETED_NEGATIVE
    #   (B) the input could not SUPPORT the analysis -> status BLOCKED
    # Printing (A)'s wording for a (B) outcome would claim "we looked and found
    # nothing" when in fact nothing was looked for — an overclaim, and precisely the
    # conflation outcome_type exists to prevent.
    if st.outcome_type is OutcomeType.SCIENTIFIC_NEGATIVE:
        if st.status is Status.BLOCKED:
            kind = (
                "This is a SCIENTIFIC NEGATIVE OF THE SECOND KIND: the input could not "
                "support the analysis. NOTHING WAS TESTED, so this run makes no "
                "statement at all about whether the effect is present — it is not a "
                "finding of 'no effect'. The cause lies in the evidence base, not in "
                "the code, and downstream stages were correctly prevented from running.")
        else:
            kind = (
                "This is a TRUE SCIENTIFIC NEGATIVE, not a technical problem. The "
                "analysis ran correctly and the result is that the effect was not "
                "present. This is a complete, valid outcome and it licenses no "
                "modification of the method.")
    else:
        kind = {
            OutcomeType.TECHNICAL_FAILURE: (
                "This is a TECHNICAL FAILURE, not a scientific result. Nothing can be "
                "concluded about the biology from this stage."),
            OutcomeType.NOT_APPLICABLE: (
                "This stage did not run because an upstream stage terminated the chain. "
                "See the upstream stage's NOT_RUN.txt for the originating cause."),
            # Workflow v2 §5.4 / Appendix A. The third kind, and the one this
            # pipeline previously had no way to say: the analysis RAN, rejected
            # nothing, and COULD NOT HAVE rejected anything. It is neither a
            # finding nor a fault, and the wording must not let it be read as
            # either — Appendix A prohibits "no hotspot is detectable in this
            # gene" unless the §5.4 certificate passed.
            OutcomeType.UNINFORMATIVE: (
                "This run is UNINFORMATIVE, and that is neither a scientific result "
                "nor a technical failure. Under the configuration executed, the test "
                "could not have rejected any center regardless of the data: the "
                "smallest p-value attainable (p_floor) exceeded the rank-1 "
                "Benjamini-Hochberg critical value (q/m). NOTHING may be concluded "
                "about the presence or absence of hotspots in this gene. See "
                "power_certificate.json for p_floor, q/m, their ratio, and which "
                "floor is binding — a permutation-resolution floor is remediable by "
                "increasing B, a combinatorial floor is not."),
            OutcomeType.COMPLETED: "This stage completed.",
        }[st.outcome_type]
    text = (
        f"STAGE NOT COMPLETED: {st.stage}\n"
        f"{'=' * (21 + len(st.stage))}\n\n"
        f"Status       : {st.status.value}\n"
        f"Outcome type : {st.outcome_type.value}\n"
        f"Owner        : {st.agent_owner}\n"
        f"Detected     : {st.ended_utc or st.started_utc or 'NA'}\n\n"
        f"What happened\n-------------\n{st.reason or 'No reason recorded.'}\n\n"
        f"How to read this\n----------------\n{kind}\n\n"
        f"Recommended action\n------------------\n"
        f"{st.recommended_action or 'None recorded.'}\n"
    )
    (stage_dir / "NOT_RUN.txt").write_text(text, encoding="utf-8", newline="\n")
