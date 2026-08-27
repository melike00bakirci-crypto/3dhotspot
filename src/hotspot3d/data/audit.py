"""F12 enforcement — a mechanical proof that no Stage A branch keys on a star.

The Plan forbids review stars from ever acting as an inclusion criterion. A
prose promise is unfalsifiable, so this module parses the AST of every module in
the Stage A code path and rejects any *decision position* whose test expression
mentions review-status or star vocabulary.

Decision positions audited:
  ``if`` / ``elif``, conditional expressions, ``while``, comprehension ``if``
  clauses, ``assert`` tests, and calls to ``filter`` / ``itertools.filterfalse``.

Reporting positions are deliberately NOT audited: ``max_star`` and
``review_star_distribution.tsv`` are required outputs, and the >=1*/>=2*
sensitivity channel of METHOD_SPEC II.7 must be able to group by star. The
prohibition is on *selection*, not on *description*.

The complementary runtime guarantee — that permuting every star level leaves the
cohort bit-identical — is asserted in ``tests/data_structure``.
"""
from __future__ import annotations

import ast
from dataclasses import dataclass
from importlib import import_module
from pathlib import Path

from ..utils.errors import LeakageError

#: Vocabulary that may never appear in a Stage A decision position.
STAR_TOKENS: frozenset[str] = frozenset({
    "star", "stars", "star_level", "star_levels", "starlevel", "max_star",
    "min_star", "review_status", "reviewstatus", "review_stars", "review_star",
    "review_status_stars", "star_filter", "review_star_filter",
})

#: Every module that can influence which record or residue enters the cohort.
STAGE_A_AUDITED_MODULES: tuple[str, ...] = (
    "hotspot3d.data.clinvar",
    "hotspot3d.data.classify",
    "hotspot3d.data.cohort",
    "hotspot3d.data.review_status",
    "hotspot3d.data.hgvs",
    "hotspot3d.data.sources",
    "hotspot3d.data.stage",
    "hotspot3d.structure.sources",
    "hotspot3d.structure.cif",
    "hotspot3d.structure.mapping",
    "hotspot3d.structure.plddt",
    "hotspot3d.structure.qc",
)

_FILTER_CALLS = frozenset({"filter", "filterfalse", "takewhile", "dropwhile", "query"})


@dataclass(frozen=True)
class AuditFinding:
    module: str
    path: str
    lineno: int
    position: str
    token: str

    def as_text(self) -> str:
        return f"{self.path}:{self.lineno} [{self.position}] mentions {self.token!r}"


def _tokens_in(node: ast.AST) -> set[str]:
    """Every identifier-like token appearing anywhere in a subtree."""
    found: set[str] = set()
    for sub in ast.walk(node):
        if isinstance(sub, ast.Name):
            found.add(sub.id.lower())
        elif isinstance(sub, ast.Attribute):
            found.add(sub.attr.lower())
        elif isinstance(sub, ast.arg):
            found.add(sub.arg.lower())
        elif isinstance(sub, ast.keyword) and sub.arg:
            found.add(sub.arg.lower())
        elif isinstance(sub, ast.Constant) and isinstance(sub.value, str):
            found.add(sub.value.strip().lower())
    return found


def _decision_positions(tree: ast.AST):
    """Yield ``(position_name, test_node)`` for every decision in the module."""
    for node in ast.walk(tree):
        if isinstance(node, ast.If):
            yield "if", node.test
        elif isinstance(node, ast.IfExp):
            yield "conditional_expression", node.test
        elif isinstance(node, ast.While):
            yield "while", node.test
        elif isinstance(node, ast.Assert):
            yield "assert", node.test
        elif isinstance(node, ast.comprehension):
            for test in node.ifs:
                yield "comprehension_if", test
        elif isinstance(node, ast.Call):
            name = getattr(node.func, "id", None) or getattr(node.func, "attr", None)
            if str(name or "").lower() in _FILTER_CALLS:
                yield f"call:{name}", node


def audit_module_source(source: str, module_name: str, path: str) -> list[AuditFinding]:
    """Return every star-keyed decision found in one module's source."""
    tree = ast.parse(source, filename=path)
    findings: list[AuditFinding] = []
    for position, test in _decision_positions(tree):
        hits = _tokens_in(test) & STAR_TOKENS
        for token in sorted(hits):
            findings.append(AuditFinding(module=module_name, path=path,
                                         lineno=getattr(test, "lineno", 0),
                                         position=position, token=token))
    return findings


def assert_no_star_based_inclusion(
    modules: tuple[str, ...] = STAGE_A_AUDITED_MODULES,
) -> dict:
    """Assert F12 mechanically. Raises :class:`LeakageError` on any violation.

    Returns the audit record written into ``cohort_summary.json`` and provenance,
    so the check is not merely performed but *evidenced*.
    """
    findings: list[AuditFinding] = []
    scanned: list[str] = []
    for name in modules:
        module = import_module(name)
        path = Path(getattr(module, "__file__", "") or "")
        if not path.is_file():
            raise LeakageError(
                f"F12 audit cannot verify {name!r}: module source not on disk. "
                f"An unverifiable code path is treated as a failed audit."
            )
        scanned.append(name)
        findings.extend(audit_module_source(path.read_text(encoding="utf-8"), name, str(path)))

    if findings:
        raise LeakageError(
            "F12 VIOLATION — a Stage A decision keys on ClinVar review status or "
            "review star. Review stars are metadata and may never be an inclusion "
            "criterion, under any name, in any code path:\n  "
            + "\n  ".join(f.as_text() for f in findings)
        )

    return {
        "check": "no_star_based_inclusion",
        "frozen_decision": "F12",
        "passed": True,
        "modules_scanned": scanned,
        "n_modules_scanned": len(scanned),
        "tokens_prohibited_in_decisions": sorted(STAR_TOKENS),
        "decision_positions_audited": [
            "if", "conditional_expression", "while", "assert", "comprehension_if",
            "filter/filterfalse/takewhile/dropwhile/query calls",
        ],
        "note": (
            "Reporting positions (max_star, review_star_distribution.tsv, the "
            "separate >=1*/>=2* sensitivity channel) are intentionally outside the "
            "audit: the prohibition is on selection, not description."
        ),
    }
