"""The UNDERPOWERED terminal state, end to end through the orchestrator.

Workflow v2 §0 and §5.4. Stage B's own tests establish that the power certificate
fails and that Stage B declares `terminal_state = "UNDERPOWERED"`. This module
covers the half that lives in the Lead's code: that the declaration survives the
orchestrator and reaches the run.

That half is where the original defect lived. Stage B returns a complete handoff
rather than raising, so before the orchestrator hook existed, the pipeline called
Stage CD next, which refused on `qc_status = FAIL` with a `BlockedError` — and the
run was typed `BLOCKED / TECHNICAL_FAILURE`. An uninformative run reported as a
pipeline fault is the mirror of the misreport that motivated v2, and Appendix A
forbids the other mistake — reporting it as a negative — just as firmly.

**Trigger: an existing, purpose-built fixture — not a config override.**

An earlier version of this module reached UNDERPOWERED by overlaying
``permutation.B_default`` down to 200. That is invalid: ``B_default`` is FROZEN at
10000 (F7/F11), `run_pipeline` calls `assert_frozen_methodology` before any stage
runs, and the overlay was correctly rejected — a `BlockedError` at setup, not the
terminal state under test.

The permutation-resolution floor `p_res = 1/(B+1)` is fixed once `B` is fixed. At
the frozen `B = 10000`, `p_res = 1/10001 = 9.999e-05`. The rank-1 critical value is
`c_1 = q/m`. Whenever the REALIZED test-family size `m` exceeds `q*(B+1) = 500.05`,
`c_1 < p_res`, so `p_floor = max(p_res, p_comb_best) >= p_res > c_1` — the run is
UNDERPOWERED **regardless of the data**, because the floor here is set by `B` and
`m` alone, never by how the P/LP residues happen to be arranged.

`tests/fixtures/synthetic_clinvar.py` (`data-structure`, frozen, used unmodified)
already ships exactly this case: `"permutation_resolution"` — "a large U_struct
(m > q * (B + 1) = 500 ...) with a handful of extremely tight P/LP residues whose
empirical p pins at that floor — the II.4 conditions for
PERMUTATION_RESOLUTION_LIMITED." It was built for the diagnostic v1 could only
warn about; v2 §5.4 is the rule that says what such a run must terminate as. No
new fixture is constructed and no cohort parameter is chosen to force this
outcome — the case already exists, for this exact regime, at the frozen `B`.

**Why a specific `run_id` is pinned.** The realized family size `m` — how many
`U_center` residues end up with >=1 classified residue in-sphere at the SELECTED
`r_hot` — depends on which candidate radius the multi-objective optimization
picks, and near a tie that pick can move with the permutation draw, exactly as
`test_seed_determinism.py::test_clustered_conclusion_is_seed_stable_within_the_
radius_grid` already documents for the `clustered` case. At the default
`FIXED_RUN_ID`-style id the case lands at `m = 497`, three short of the
`m > 500.05` line — PASSING the certificate by 0.6%, not failing it. That is a
genuine near-miss, not a defect in the fixture or the certificate; it is why the
diagnostic's own design note says "m > 500" as an approximate target rather than
a guarantee at every seed. Sweeping six deterministic, pre-registered-shaped
run_ids (the same sanctioned mechanism the KCNA2 seed-sensitivity work used, per
METHOD_SPEC II.12 — seeds legitimately derive from run_id) found one,
`"20250101T000000Z_seedE001_seedE002"`, where the SAME frozen case, frozen `B`,
frozen radius domain and frozen cohort realize `m = 529`, clearing the line by a
comfortable margin (`p_res/c_1 ~= 1.058`, not a coin flip). Nothing about the
gene, the cohort, or the methodology changed between the two seeds — only which
of several legitimate permutation draws the run happened to use, which is
precisely what II.12 licenses varying.

The only config overlay used is the standard integration one
(`tests/integration/conftest.py`, session-scoped `config` fixture): a reduced
`robustness.N_CAP` for wall-clock reasons only. It is irrelevant here regardless,
since Stage D never runs on an underpowered chain. `permutation.B_default` is
never overlaid and is asserted equal to the frozen 10000 by the first test below.
"""
from __future__ import annotations

import json

import pytest

from hotspot3d.utils.status import OutcomeType, RunStatus

pytestmark = [pytest.mark.integration, pytest.mark.slow]

#: The `permutation_resolution` case at its default-style run_id realizes
#: m=497 (certificate PASSES by 0.6%). This id was found — by sweeping several
#: deterministic ids, the II.12-sanctioned degree of freedom — to realize m=529
#: at the SAME frozen case, cohort, B and radius domain, clearing the m>500.05
#: line decisively (p_res/c_1 ~= 1.058). See the module docstring.
UNDERPOWERED_RUN_ID = "20250101T000000Z_seedE001_seedE002"


@pytest.fixture(scope="module")
def underpowered_run(tmp_path_factory, config):
    from hotspot3d.orchestration.pipeline import SourceBundle, run_pipeline
    from tests.fixtures.synthetic_annotation import MockAnnotationSource
    from tests.fixtures.synthetic_clinvar import make_synthetic_case

    assert config.get("permutation.B_default") == 10_000, (
        "this test must exercise the frozen production B, not an overlay")

    tmp = tmp_path_factory.mktemp("underpowered")
    fixture = make_synthetic_case("permutation_resolution")
    sources = SourceBundle(variant_source=fixture.source,
                           structure_source=fixture.source,
                           annotation_source=MockAnnotationSource())
    return run_pipeline(
        gene="SYNTH", config=config, results_root=tmp / "results",
        sources=sources, synthetic=True, run_id=UNDERPOWERED_RUN_ID)


def _certificate(result):
    path = (result.run_root / "FULL_RESULTS" / "06_FINAL_HOTSPOTS"
            / "power_certificate.json")
    assert path.is_file(), "the power certificate is missing"
    return json.loads(path.read_text())


def test_the_fixture_is_genuinely_underpowered_at_the_frozen_b(underpowered_run):
    """Confirms the premise before trusting any downstream assertion.

    `B_default` was never overlaid — `config.get("permutation.B_default")` reads
    the shipped, frozen 10000 for this run. If the realized family size did not
    clear `q*(B+1)`, every assertion below would be vacuous, so this is checked
    directly against the certificate rather than assumed.
    """
    cert = _certificate(underpowered_run)
    posthoc = cert["posthoc"]
    assert posthoc["B"] == 10_000, "B was not the frozen production value"
    assert posthoc["p_res"] == pytest.approx(1.0 / 10_001, rel=1e-9)
    m = posthoc["m_test_family_size"]
    q = 0.05
    assert m > q * 10_001, (
        f"m={m} does not clear q*(B+1)={q * 10_001:g} at the frozen B — the "
        f"fixture did not reach the regime this test exercises")
    assert posthoc["binding_floor"] == "permutation_resolution"
    assert posthoc["p_floor"] == pytest.approx(posthoc["p_res"])


def test_the_run_terminates_underpowered_and_not_as_a_fault(underpowered_run):
    result = underpowered_run
    assert result.run_status is RunStatus.UNDERPOWERED, (
        f"expected UNDERPOWERED, got {result.run_status} ({result.error})")
    assert result.outcome_type is OutcomeType.UNINFORMATIVE
    assert result.exit_code == 3, "a caller reading only the exit code must not " \
                                  "mistake this for a completed analysis"
    assert result.error is None, "an uninformative run is not a technical failure"


def test_it_is_never_reported_as_a_negative(underpowered_run):
    """Workflow v2 Appendix A — the prohibited inference."""
    result = underpowered_run
    assert result.negative_result is None, (
        "an underpowered run carried a negative_result; v2 Appendix A forbids any "
        "claim about the absence of hotspots when the certificate fails")
    assert result.run_status is not RunStatus.COMPLETED_NEGATIVE

    summary = result.terminal_summary
    assert "NEGATIVE RESULT" not in summary
    for banned in ("No 3D hotspot is detectable",
                   "valid, complete scientific negative",
                   "licenses no modification"):
        assert banned not in summary, f"prohibited inference in summary: {banned!r}"


def test_it_is_never_reported_as_blocked(underpowered_run):
    """The defect the orchestrator hook fixed: this used to come back BLOCKED."""
    result = underpowered_run
    assert result.run_status is not RunStatus.BLOCKED
    assert result.outcome_type is not OutcomeType.TECHNICAL_FAILURE


def test_the_certificate_is_written_and_explains_which_floor_binds(underpowered_run):
    """§5.4 + §11 — a first-class artifact in every terminal state."""
    cert = _certificate(underpowered_run)
    assert cert["PASSES"] is False
    posthoc = cert["posthoc"]
    assert posthoc["p_floor"] > posthoc["c_1_rank1_critical_value"]
    for key in ("N_P", "N_B", "m_test_family_size", "B", "p_res", "p_comb_best",
               "p_floor"):
        assert key in posthoc, f"certificate omits {key}"

    # p_res binds here, so the remedy is a concrete B, not the "increasing B is
    # futile" wording the combinatorial branch must give instead.
    assert posthoc["binding_floor"] == "permutation_resolution", posthoc
    assert posthoc["increasing_B_is_futile"] is False
    assert posthoc.get("B_to_clear_c_1")


def test_downstream_stages_are_uninformative_not_failed(underpowered_run):
    """The chain stopped because nothing could be concluded, not because it broke."""
    fr = underpowered_run.run_root / "FULL_RESULTS"
    for stage in ("07_FOOTPRINT_RADIUS", "08_FINAL_FOOTPRINT", "09_ROBUSTNESS",
                  "10_ANNOTATION"):
        status = json.loads((fr / stage / "stage_status.json").read_text())
        assert status["status"] == "NOT_RUN", f"{stage}: {status['status']}"
        assert status["outcome_type"] != OutcomeType.TECHNICAL_FAILURE.value, (
            f"{stage} was typed a technical failure on an underpowered run")
        assert (fr / stage / "NOT_RUN.txt").is_file()
