"""The pipeline-qa-debugger's fixed vocabulary — Lead-owned, structural, not advisory.

Three things live here, each because leaving it to per-incident judgment is exactly
what this module exists to prevent:

1. :class:`FailureClass` / :class:`ConfigSubclass` — the closed taxonomy every QA
   classification must use. Adding a category is a Lead decision, not something a
   bug brief invents in passing.
2. :func:`classify_config_key` — the INFRA_CONFIG / SCIENTIFIC_CONFIG split, keyed
   by pattern over the real ``config/pipeline.yaml`` namespace, defaulting an
   UNRECOGNIZED key to SCIENTIFIC_CONFIG. That default is deliberate: this file is
   named "FROZEN PIPELINE CONFIGURATION" for a reason, and a config key nobody has
   classified yet must fail closed (routed to a module owner) rather than open
   (silently editable by QA).
3. The QA write-access manifest (:data:`QA_DIRECT_WRITE_GLOBS`,
   :data:`QA_SYMBOL_CARVEOUTS`, :func:`qa_may_edit_directly`,
   :func:`owner_of_path`) — which files QA may repair itself versus which must be
   routed to the owning scientific agent, down to individual symbols where a pure
   glue concern (e.g. HTTP retry/backoff) sits inside an otherwise scientific file.

:class:`BugBrief` is the structured artifact every QA finding must produce, and
:class:`RepairTrail` is the bounded-repair-loop tracker that stops a fix from being
attempted forever against the same root cause.
"""
from __future__ import annotations

import fnmatch
import json
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

from .runctx import STAGE_OWNERS

# --------------------------------------------------------------------------- #
# 1. the fixed failure taxonomy                                               #
# --------------------------------------------------------------------------- #


class FailureClass(str, Enum):
    """Every QA classification is EXACTLY one of these six. Stable by design —
    see the module docstring. Do not add a seventh without a Lead decision."""

    IMPLEMENTATION_BUG = "IMPLEMENTATION_BUG"
    SCIENTIFIC_NEGATIVE = "SCIENTIFIC_NEGATIVE"
    UNDERPOWERED = "UNDERPOWERED"
    BLOCKED = "BLOCKED"
    EXTERNAL_DATA_FAILURE = "EXTERNAL_DATA_FAILURE"
    CONFIGURATION_PROBLEM = "CONFIGURATION_PROBLEM"


class ConfigSubclass(str, Enum):
    """CONFIGURATION_PROBLEM has exactly two subclasses, and every
    CONFIGURATION_PROBLEM brief must state which one applies."""

    INFRA_CONFIG = "INFRA_CONFIG"
    SCIENTIFIC_CONFIG = "SCIENTIFIC_CONFIG"


#: What kind of evidence artifact each failure class requires. Enforced by
#: :meth:`BugBrief.validate` — a classification without the right evidence TYPE is
#: invalid, not merely under-documented.
REQUIRED_EVIDENCE_KIND: dict[FailureClass, str] = {
    FailureClass.IMPLEMENTATION_BUG: (
        "log excerpt, traceback, or expected-vs-actual schema/output diff"),
    FailureClass.SCIENTIFIC_NEGATIVE: (
        "computed values establishing the terminal state was statistically correct "
        "(e.g. the power certificate that PASSED, the realized center set)"),
    FailureClass.UNDERPOWERED: (
        "power-certificate fields (p_floor, c_1, ratio, binding_floor) establishing "
        "the design could not have rejected, not merely that it didn't"),
    FailureClass.BLOCKED: (
        "the specific precondition that failed and the log line that raised it"),
    FailureClass.EXTERNAL_DATA_FAILURE: (
        "network/HTTP/retry log evidence showing the failure originated outside "
        "the pipeline's own logic"),
    FailureClass.CONFIGURATION_PROBLEM: (
        "the specific config key, its current value, and its ConfigSubclass"),
}


# --------------------------------------------------------------------------- #
# 2. INFRA_CONFIG vs SCIENTIFIC_CONFIG — pattern-based, fail-closed            #
# --------------------------------------------------------------------------- #

#: Prefixes over the config/pipeline.yaml dotted-key namespace that carry NO
#: scientific semantics: paths, execution mode, thread/worker counts, output
#: encoding/formatting/retention/archival mechanics. Enumerated deliberately as an
#: ALLOWLIST rather than trying to enumerate the ~150 scientific keys — anything
#: not matched here is SCIENTIFIC_CONFIG by default (see module docstring).
#:
#: `execution.allow_network` / `execution.offline_mode` gate WHETHER a run touches
#: the network, not what counts as significant, included, or admissible — an
#: execution-mode toggle, not a scientific decision. `execution.n_jobs` is a
#: worker count, the textbook INFRA_CONFIG example. `output.*` is entirely
#: serialization/retention/archival mechanics (encoding, TSV formatting, DPI,
#: exit-code mapping, what gets tarred). `schema_version`/`config_version` are
#: file-format metadata.
INFRA_CONFIG_PREFIXES: tuple[str, ...] = (
    "execution.",
    "output.",
    "schema_version",
    "config_version",
)

#: Explicit, audited traps: keys that LOOK infra-shaped (a data-source name, a
#: boolean toggle, a numeric threshold with an operational-sounding name) but are
#: FROZEN scientific-methodology decisions in this codebase. Listed even though
#: they would already fall through to the SCIENTIFIC_CONFIG default, because the
#: whole point of Section 2's split is that this must never be "per-incident
#: judgment" — a reviewer should be able to find the reasoning here, not re-derive
#: it. This is CASE 9's worked example: `structure.source` reads like "which data
#: source to hit," but it is II.13's FROZEN acquisition-methodology pin, not an
#: infrastructure choice.
SCIENTIFIC_CONFIG_TRAPS: dict[str, str] = {
    "structure.source": (
        "which structural database to use is a II.13 FROZEN acquisition-"
        "methodology decision, not a data-source convenience setting"),
    "structure.model_template": (
        "which AlphaFold model VERSION is analysed; a pre-registered revision "
        "under a new RUN_ID, never a routine config edit (see the Lead ruling "
        "recorded beside this key in config/pipeline.yaml)"),
    "seeding.MASTER_SEED": (
        "FROZEN (II.12); changing it changes every derived permutation stream — "
        "exactly the 'change seeds to select a favorable result' the firewall "
        "prohibits"),
    "seeding.derivation": "FROZEN (II.12), same reasoning as MASTER_SEED",
    "robustness.N_CAP": (
        "a config knob with a stated cost model, but explicitly named "
        "'robustness budgets' in the firewall (Section 3) — the Lead sets it, "
        "never QA, regardless of how expensive an escalation looks"),
    "robustness.max_wall_seconds": "same reasoning as robustness.N_CAP",
    "footprint_geometry.voxel_h_fallback_A": (
        "a computational-feasibility fallback whose VALUE changes the geometric "
        "resolution actually used — it changes the computed footprint, not just "
        "whether the run finishes"),
    "output.tsv.float_significant_digits": (
        "borderline: a DISPLAY/serialization precision, not a computation "
        "precision — kept as INFRA_CONFIG via the output.* prefix, noted here so "
        "a future reviewer sees the call was made deliberately, not missed"),
}


def classify_config_key(dotted_key: str) -> ConfigSubclass:
    """INFRA_CONFIG if (and only if) the key matches an allowlisted infra prefix;
    SCIENTIFIC_CONFIG otherwise, including every key nobody has classified yet.

    This is CASE 9's mechanism: ``structure.source`` does not match any
    ``INFRA_CONFIG_PREFIXES`` entry, so it falls through to SCIENTIFIC_CONFIG
    regardless of how it reads on its face.
    """
    for prefix in INFRA_CONFIG_PREFIXES:
        if dotted_key == prefix or dotted_key.startswith(prefix):
            return ConfigSubclass.INFRA_CONFIG
    return ConfigSubclass.SCIENTIFIC_CONFIG


# --------------------------------------------------------------------------- #
# 3. QA write-access / ownership manifest                                     #
# --------------------------------------------------------------------------- #

#: File-glob patterns QA may edit DIRECTLY — pure engineering/glue territory.
#: Everything else routes to the owning scientific agent (see
#: :func:`owner_of_path`), with the narrow exception of
#: :data:`QA_SYMBOL_CARVEOUTS`.
QA_DIRECT_WRITE_GLOBS: tuple[str, ...] = (
    # cross-stage plumbing, provenance, reporting, archival — Lead-owned
    # infrastructure that QA may repair without a module-owner round trip.
    "src/hotspot3d/orchestration/pipeline.py",
    "src/hotspot3d/orchestration/cli.py",
    "src/hotspot3d/reporting/*.py",
    "src/hotspot3d/utils/io.py",
    "src/hotspot3d/utils/hashing.py",
    "src/hotspot3d/utils/errors.py",
    "src/hotspot3d/utils/status.py",
    "src/hotspot3d/utils/runctx.py",
    "src/hotspot3d/utils/stage_contract.py",
    "src/hotspot3d/utils/qa_taxonomy.py",
    # QA's own tests and regression fixtures.
    "tests/integration/test_qa_*.py",
    "tests/regression/**",
)

#: Symbol-level carve-outs: a FILE that is otherwise a scientific module's
#: territory, but specific NAMES inside it carry no scientific semantics — the
#: retry/backoff/chunking constants of the ClinVar downloader are the canonical
#: example (Section 2's "regardless of whether it lives in a config file," applied
#: to source instead of YAML). QA may edit these symbols directly; touching
#: anything else in the same file still routes to the owner.
QA_SYMBOL_CARVEOUTS: dict[str, tuple[str, ...]] = {
    "src/hotspot3d/data/sources.py": (
        "CLINVAR_CHUNK_BYTES", "CLINVAR_MAX_ATTEMPTS", "CLINVAR_BACKOFF_BASE_S",
        "CLINVAR_STREAM_BYTES", "ClinVarClient.timeout", "UniProtClient.timeout",
    ),
    "src/hotspot3d/structure/sources.py": (
        "AlphaFoldClient.timeout",
    ),
}


def owner_of_path(path: str | Path) -> str:
    """The scientific agent (or 'lead') that owns a repo-relative path.

    Reuses :data:`STAGE_OWNERS` for anything under ``FULL_RESULTS/<NN_STAGE>/``.
    For source files, walks ``src/hotspot3d/<package>/`` against the same
    ownership partition :data:`AGENT_WRITE_SCOPE`'s stage list implies (a package
    belongs to whichever agent owns the stages it implements) — kept as an
    explicit table here rather than re-derived, so ownership is one lookup, not
    an inference.
    """
    rel = str(path).replace("\\", "/")
    for stage, owner in STAGE_OWNERS.items():
        if f"FULL_RESULTS/{stage}/" in rel or rel.startswith(f"{stage}/"):
            return owner

    package_owners = {
        "src/hotspot3d/data/": "data-structure",
        "src/hotspot3d/structure/": "data-structure",
        "src/hotspot3d/spatial/": "hotspot-statistics",
        "src/hotspot3d/hotspot/": "hotspot-statistics",
        "src/hotspot3d/footprint/": "footprint-robustness",
        "src/hotspot3d/robustness/": "footprint-robustness",
        "src/hotspot3d/annotation/": "biological-annotation",
        "src/hotspot3d/orchestration/": "lead",
        "src/hotspot3d/reporting/": "lead",
        "src/hotspot3d/utils/": "lead",
    }
    for prefix, owner in package_owners.items():
        if rel.startswith(prefix):
            return owner
    return "lead"


def qa_may_edit_directly(path: str | Path, symbol: str | None = None) -> bool:
    """True iff QA may repair ``path`` (optionally scoped to ``symbol``) itself.

    A glob match on :data:`QA_DIRECT_WRITE_GLOBS` is sufficient on its own. A
    match ONLY in :data:`QA_SYMBOL_CARVEOUTS` requires ``symbol`` to be given and
    to be one of the carved-out names — omitting ``symbol`` against a carve-out
    file returns False, so a caller cannot claim "the whole file is fine" by
    forgetting to name what it's actually touching.
    """
    rel = str(path).replace("\\", "/")
    if any(fnmatch.fnmatch(rel, pat) for pat in QA_DIRECT_WRITE_GLOBS):
        return True
    carved = QA_SYMBOL_CARVEOUTS.get(rel)
    if carved and symbol is not None:
        return symbol in carved
    return False


# --------------------------------------------------------------------------- #
# 4. the bug brief — Section 5's protocol, as an enforced structure            #
# --------------------------------------------------------------------------- #


@dataclass
class BugBrief:
    """One QA finding, in the exact shape Section 5 requires.

    ``validate()`` is what makes the protocol structural rather than advisory:
    Lead must reject a brief that fails it (see
    :func:`assert_brief_acceptable_to_lead`), so a brief missing its required
    evidence type, or one that both claims a methodology change is required AND
    that QA applied a fix, cannot silently pass as "documented."
    """

    bug_id: str
    failure_class: FailureClass
    stage: str
    owner: str
    observed_behavior: str
    expected_behavior: str
    minimal_reproduction: str
    root_cause: str
    evidence_artifact: str
    why_this_is_an_implementation_defect: str
    scientific_methodology_change_required: bool
    files_or_components_affected: list[str]
    targeted_test_required: str
    config_subclass: ConfigSubclass | None = None
    #: The literal dotted config key, REQUIRED whenever failure_class is
    #: CONFIGURATION_PROBLEM. Having the key as its own field (not just prose
    #: inside EVIDENCE_ARTIFACT) is what lets :meth:`validate` re-derive the
    #: correct subclass via :func:`classify_config_key` and catch a
    #: misclassification structurally (CASE 9) rather than trusting whatever
    #: subclass was typed in.
    config_key: str | None = None
    #: Set True only once a fix has actually been applied and re-validated —
    #: never at brief-creation time. Checked against
    #: ``scientific_methodology_change_required`` in :meth:`validate`.
    fix_applied_by_qa: bool = False
    repair_attempt_number: int = 1
    root_cause_signature: str | None = None

    def validate(self) -> list[str]:
        """Structural problems with this brief. Empty list == acceptable to Lead."""
        problems: list[str] = []

        if self.failure_class is FailureClass.CONFIGURATION_PROBLEM:
            if self.config_subclass is None:
                problems.append(
                    "FAILURE_CLASS=CONFIGURATION_PROBLEM requires CONFIG_SUBCLASS "
                    "(INFRA_CONFIG or SCIENTIFIC_CONFIG); none given")
            if not self.config_key:
                problems.append(
                    "FAILURE_CLASS=CONFIGURATION_PROBLEM requires CONFIG_KEY (the "
                    "literal dotted key) so the declared subclass can be checked "
                    "against the manifest, not merely asserted")
            elif self.config_subclass is not None:
                # CASE 9: re-derive from the ONE source of truth and compare —
                # a declared subclass is never trusted on its own.
                actual = classify_config_key(self.config_key)
                if actual is not self.config_subclass:
                    problems.append(
                        f"CONFIG_SUBCLASS={self.config_subclass.value} was declared "
                        f"for {self.config_key!r}, but classify_config_key reports "
                        f"{actual.value} — a SCIENTIFIC_CONFIG value must never be "
                        f"treated as INFRA_CONFIG regardless of how it was labelled")
        elif self.config_subclass is not None:
            problems.append(
                f"CONFIG_SUBCLASS is only valid when FAILURE_CLASS="
                f"CONFIGURATION_PROBLEM, got {self.failure_class.value}")

        if not self.evidence_artifact.strip():
            problems.append(
                f"EVIDENCE_ARTIFACT is empty; {self.failure_class.value} requires "
                f"{REQUIRED_EVIDENCE_KIND[self.failure_class]}")

        if self.scientific_methodology_change_required and self.fix_applied_by_qa:
            problems.append(
                "SCIENTIFIC_METHODOLOGY_CHANGE_REQUIRED=TRUE but "
                "fix_applied_by_qa=TRUE — QA must never perform a methodology "
                "change; this brief should have escalated instead")

        if not self.owner.strip():
            problems.append("OWNER is empty")

        if self.fix_applied_by_qa and not qa_may_edit_directly_hint(self):
            problems.append(
                f"fix_applied_by_qa=TRUE but owner={self.owner!r} is a scientific "
                f"agent, not 'lead' or an infra path — a QA-applied fix outside "
                f"the write manifest is invalid regardless of severity")

        return problems

    def as_dict(self) -> dict[str, Any]:
        d = {
            "BUG_ID": self.bug_id,
            "FAILURE_CLASS": self.failure_class.value,
            "CONFIG_SUBCLASS": self.config_subclass.value if self.config_subclass else None,
            "CONFIG_KEY": self.config_key,
            "STAGE": self.stage,
            "OWNER": self.owner,
            "OBSERVED_BEHAVIOR": self.observed_behavior,
            "EXPECTED_BEHAVIOR": self.expected_behavior,
            "MINIMAL_REPRODUCTION": self.minimal_reproduction,
            "ROOT_CAUSE": self.root_cause,
            "EVIDENCE_ARTIFACT": self.evidence_artifact,
            "WHY_THIS_IS_AN_IMPLEMENTATION_DEFECT": self.why_this_is_an_implementation_defect,
            "SCIENTIFIC_METHODOLOGY_CHANGE_REQUIRED": self.scientific_methodology_change_required,
            "FILES_OR_COMPONENTS_AFFECTED": list(self.files_or_components_affected),
            "TARGETED_TEST_REQUIRED": self.targeted_test_required,
            "fix_applied_by_qa": self.fix_applied_by_qa,
            "repair_attempt_number": self.repair_attempt_number,
            "root_cause_signature": self.root_cause_signature,
        }
        return d


def qa_may_edit_directly_hint(brief: BugBrief) -> bool:
    """Loose structural check backing :meth:`BugBrief.validate`'s last rule:
    a brief that claims QA applied the fix must name an owner QA is actually
    allowed to touch (itself, or 'lead' for the shared infra it may edit)."""
    if brief.owner in ("lead", "pipeline-qa-debugger"):
        return True
    if brief.failure_class is FailureClass.EXTERNAL_DATA_FAILURE:
        return True
    if (brief.failure_class is FailureClass.CONFIGURATION_PROBLEM
            and brief.config_subclass is ConfigSubclass.INFRA_CONFIG):
        return True
    return False


def assert_brief_acceptable_to_lead(brief: BugBrief) -> None:
    """Section 11 — Lead must reject a brief lacking its required evidence
    artifact, or any structurally invalid brief. Raises, does not warn."""
    from .errors import BlockedError

    problems = brief.validate()
    if problems:
        raise BlockedError(
            f"BLOCKED — bug brief {brief.bug_id} rejected by Lead:\n"
            + "\n".join(f"  - {p}" for p in problems))


# --------------------------------------------------------------------------- #
# 5. bounded repair loop                                                      #
# --------------------------------------------------------------------------- #

#: Section 8's default: at most this many owner-fix cycles for the SAME
#: (stage, root_cause_signature) before Lead must stop and report, not retry again.
DEFAULT_MAX_REPAIR_ATTEMPTS = 2


@dataclass
class RepairAttempt:
    attempt_number: int
    root_cause_signature: str
    bug_id: str
    owner: str
    outcome: str  # "pending" | "fixed" | "reappeared" | "failed"


@dataclass
class RepairTrail:
    """Tracks repair attempts for one (run_id, stage) pair, JSON-persisted so the
    trail survives across separate agent invocations within a run.

    ``exceeded()`` is Section 8's stop condition. ``QA_REPAIR_EXHAUSTED`` in the
    module-level report is exactly ``exceeded()`` evaluated at the point Lead
    checks it.
    """

    run_id: str
    stage: str
    max_attempts: int = DEFAULT_MAX_REPAIR_ATTEMPTS
    attempts: list[RepairAttempt] = field(default_factory=list)

    def record_attempt(self, *, root_cause_signature: str, bug_id: str,
                       owner: str, outcome: str = "pending") -> RepairAttempt:
        attempt = RepairAttempt(
            attempt_number=len(self.attempts) + 1,
            root_cause_signature=root_cause_signature, bug_id=bug_id,
            owner=owner, outcome=outcome)
        self.attempts.append(attempt)
        return attempt

    def signature_reappeared(self, signature: str) -> bool:
        """True if this EXACT root-cause signature was already marked 'fixed'
        and has now resurfaced — Section 8: escalate immediately, don't re-patch."""
        fixed_before = any(a.root_cause_signature == signature and a.outcome == "fixed"
                           for a in self.attempts[:-1])
        return fixed_before and bool(self.attempts) and \
            self.attempts[-1].root_cause_signature == signature

    def exceeded(self) -> bool:
        same_stage_attempts = [a for a in self.attempts if a.outcome != "fixed"]
        return len(same_stage_attempts) > self.max_attempts

    def as_dict(self) -> dict:
        return {
            "run_id": self.run_id, "stage": self.stage,
            "max_attempts": self.max_attempts,
            "attempts": [
                {"attempt_number": a.attempt_number,
                 "root_cause_signature": a.root_cause_signature,
                 "bug_id": a.bug_id, "owner": a.owner, "outcome": a.outcome}
                for a in self.attempts
            ],
            "QA_REPAIR_EXHAUSTED": self.exceeded(),
        }

    def write(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.as_dict(), indent=2), encoding="utf-8")
        return path

    @classmethod
    def read(cls, path: Path, *, run_id: str, stage: str,
             max_attempts: int = DEFAULT_MAX_REPAIR_ATTEMPTS) -> "RepairTrail":
        if not path.is_file():
            return cls(run_id=run_id, stage=stage, max_attempts=max_attempts)
        raw = json.loads(path.read_text(encoding="utf-8"))
        trail = cls(run_id=raw["run_id"], stage=raw["stage"],
                   max_attempts=raw.get("max_attempts", max_attempts))
        trail.attempts = [
            RepairAttempt(attempt_number=a["attempt_number"],
                          root_cause_signature=a["root_cause_signature"],
                          bug_id=a["bug_id"], owner=a["owner"], outcome=a["outcome"])
            for a in raw.get("attempts", [])
        ]
        return trail


def root_cause_signature(*, stage: str, failure_class: FailureClass,
                         detail: str) -> str:
    """A stable, short identifier for "the same underlying defect" — used to
    detect reappearance (Section 8) without relying on free-text equality."""
    import hashlib

    digest = hashlib.sha256(f"{stage}|{failure_class.value}|{detail}".encode()).hexdigest()
    return f"{stage}:{failure_class.value}:{digest[:12]}"
