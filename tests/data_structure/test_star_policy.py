"""F12 — review stars are metadata, never an inclusion criterion.

Two independent guarantees are tested here, because either alone is weak:
  * a *static* guarantee — no decision anywhere in the Stage A code path
    mentions star vocabulary (AST audit);
  * a *behavioural* guarantee — permuting every review status leaves the cohort
    bit-identical, which is the property the static rule exists to protect.
"""
from __future__ import annotations

import pytest

from hotspot3d.data import audit, cohort, stage
from hotspot3d.data.audit import (
    STAGE_A_AUDITED_MODULES,
    assert_no_star_based_inclusion,
    audit_module_source,
)
from hotspot3d.data.review_status import REVIEW_STATUS_STARS, star_level
from hotspot3d.data.stage import run_stage_a
from hotspot3d.utils.errors import LeakageError

pytestmark = pytest.mark.unit


# --- static guarantee -------------------------------------------------------

def test_no_stage_a_decision_keys_on_a_review_star():
    record = assert_no_star_based_inclusion()
    assert record["passed"] is True
    assert record["frozen_decision"] == "F12"
    assert set(record["modules_scanned"]) == set(STAGE_A_AUDITED_MODULES)
    assert record["n_modules_scanned"] == len(STAGE_A_AUDITED_MODULES)


def test_the_audit_actually_catches_a_star_filter():
    """The audit must fail on code that filters by star — otherwise it proves nothing."""
    offending = (
        "def build(records):\n"
        "    return [r for r in records if r['star_level'] >= 2]\n"
    )
    findings = audit_module_source(offending, "fake", "fake.py")
    assert findings, "a star-keyed comprehension must be detected"
    assert findings[0].token == "star_level"
    assert findings[0].position == "comprehension_if"


@pytest.mark.parametrize("snippet,position", [
    ("def f(r):\n    if r.review_status == 'x':\n        return 1\n    return 0\n", "if"),
    ("def f(r):\n    return 1 if r.star_level else 0\n", "conditional_expression"),
    ("def f(rs):\n    return list(filter(lambda r: r.max_star > 1, rs))\n", "call:filter"),
])
def test_audit_detects_every_decision_position(snippet, position):
    findings = audit_module_source(snippet, "fake", "fake.py")
    assert [f.position for f in findings] == [position]


def test_audit_permits_star_reporting():
    """Describing stars is required; only selecting by them is forbidden."""
    reporting = (
        "def tally(records):\n"
        "    out = {}\n"
        "    for r in records:\n"
        "        out[r.star_level] = out.get(r.star_level, 0) + 1\n"
        "    return out\n"
        "def top(levels):\n"
        "    return max(levels)\n"
    )
    assert audit_module_source(reporting, "fake", "fake.py") == []


def test_audit_raises_leakage_error_on_violation(monkeypatch, tmp_path):
    module = tmp_path / "bad_stage_module.py"
    module.write_text("def f(r):\n    return r if r['star_level'] > 0 else None\n")

    class _Fake:
        __file__ = str(module)

    monkeypatch.setattr(audit, "import_module", lambda name: _Fake())
    with pytest.raises(LeakageError, match="F12 VIOLATION"):
        assert_no_star_based_inclusion(modules=("fake.module",))


# --- the star vocabulary itself --------------------------------------------

@pytest.mark.parametrize("status,expected", [
    ("practice guideline", 4),
    ("reviewed by expert panel", 3),
    ("criteria provided, multiple submitters, no conflicts", 2),
    ("criteria provided, single submitter", 1),
    ("no assertion criteria provided", 0),
    ("NO ASSERTION CRITERIA PROVIDED", 0),
    ("  criteria provided,   single submitter ", 1),
])
def test_star_level_mapping(status, expected):
    assert star_level(status) == expected


def test_unknown_review_status_is_na_not_zero():
    """Inventing a star for an unknown wording would be inventing metadata."""
    assert star_level("wording introduced after this release") is None
    assert star_level(None) is None
    assert 0 in set(REVIEW_STATUS_STARS.values())


# --- behavioural guarantee --------------------------------------------------

def test_permuting_every_star_leaves_the_cohort_identical(clustered, make_ctx):
    original = clustered
    permuted = clustered.with_permuted_stars()

    handoffs = []
    for case in (original, permuted):
        ctx = make_ctx()
        handoffs.append(run_stage_a(ctx, gene=case.gene, source=case.source))

    for key in ("M", "N", "N_P", "N_B", "n_conflict"):
        assert handoffs[0].payload[key] == handoffs[1].payload[key], key


def test_permuting_stars_leaves_residue_classes_identical(clustered):
    def classes(case):
        from hotspot3d.data.classify import SignificancePolicy
        from hotspot3d.data.clinvar import normalize_records
        from hotspot3d.utils.config import load_config

        cfg = load_config("config/pipeline.yaml")
        evals = cohort.evaluate_records(
            normalize_records(case.records, gene=case.gene),
            policy=SignificancePolicy.from_config(cfg),
            canonical_sequence=case.sequence,
            mane_transcript=case.mane_transcript,
        )
        rows = cohort.collapse_to_residues(evals, canonical_sequence=case.sequence)
        return {r.residue_index: r.residue_class for r in rows}

    assert classes(clustered) == classes(clustered.with_permuted_stars(seed=99))


def test_config_freezes_the_no_star_filter_policy(config):
    assert config.get("clinvar.review_star_filter") is None
    assert config.get("clinvar.use_review_stars_for_inclusion") is False
    assert ("clinvar.review_star_filter", None) in stage.FROZEN_STAGE_A_ASSERTIONS
