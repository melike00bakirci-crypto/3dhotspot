"""Stage A -> Stage B across the real handoff_01 interface.

Uses the shared fixture `tests/fixtures/synthetic_clinvar.py` (owned by
`data-structure`) so that Stage B is exercised against Stage A's ACTUAL output
schema rather than against a hand-built imitation of it — which is the only way the
column allowlist, the manifest and the review-star channel are really tested.

Stage A refuses a gene overlay (it re-hashes `config/pipeline.yaml` and requires an
exact match), so these runs use the FROZEN configuration unmodified: ``B = 10000``,
``q = 0.05``, ``BH``. They are correspondingly slow.

Assertions are about *contracts and branches*, never about specific p-values, radii
or hotspot counts — those are the pipeline's to compute, and the fixture documents
its cases the same way.
"""
from __future__ import annotations

import json

import pytest

from hotspot3d.data.stage import run_stage_a
from hotspot3d.data.synthetic import make_synthetic_case
from hotspot3d.hotspot.stage import EXPECTED_OUTPUTS, STAGES, run_stage_b
from hotspot3d.orchestration.contracts import (
    FORBIDDEN_DOWNSTREAM_COLUMNS,
    HANDOFF_REQUIRED_KEYS,
)
from hotspot3d.utils.config import load_config
from hotspot3d.utils.hashing import verify_manifest
from hotspot3d.utils.io import read_tsv
from hotspot3d.utils.runctx import RunContext

pytestmark = [pytest.mark.integration, pytest.mark.slow]


def _run(tmp_path, case_name: str, run_id: str = "20250101T000000Z_ab_test"):
    cfg = load_config("config/pipeline.yaml")
    ctx = RunContext.create(gene="SYNG", config=cfg, results_root=tmp_path / "results",
                            synthetic=True, run_id=run_id)
    case = make_synthetic_case(case_name)
    handoff_01 = run_stage_a(ctx, gene="SYNG", source=case.source)
    return ctx, handoff_01, run_stage_b(ctx, upstream=handoff_01), case


# Each case is a full Stage A + Stage B execution at the frozen B = 10000 (~15 s).
# Cases are therefore executed once per module and shared by every test that
# interrogates them. No assertion is weakened; the tests read the same completed
# runs they previously each recomputed.

@pytest.fixture(scope="module")
def clustered_run(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("clustered"), "clustered")


@pytest.fixture(scope="module")
def dispersed_run(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("dispersed"), "no_significant_hotspot")


@pytest.fixture(scope="module")
def resolution_run(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("resolution"), "permutation_resolution")


@pytest.fixture(scope="module")
def low_plddt_run(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("lowplddt"), "low_plddt")


@pytest.fixture(scope="module")
def alternate_run_id_run(tmp_path_factory):
    """The same case under a DIFFERENT RUN_ID, for the seed-derivation test."""
    return _run(tmp_path_factory.mktemp("altrunid"), "clustered",
                run_id="20250101T000000Z_bbb")


def test_stage_b_consumes_a_real_handoff_01(clustered_run):
    ctx, h1, h2, case = clustered_run

    assert h1.qc_status != "FAIL"
    assert h2.run_id == h1.run_id == ctx.run_id
    assert h2.config_sha256 == h1.config_sha256 == ctx.config.sha256
    for key in HANDOFF_REQUIRED_KEYS["handoff_02"]:
        assert key in h2.payload, key
    assert h2.payload["B"] == 10_000          # the FROZEN default, unmodified
    assert h2.payload["q"] == 0.05 and h2.payload["fdr_method"] == "BH"
    assert verify_manifest(h2.manifest, ctx.run_root) == []

    for stage in STAGES:
        assert (ctx.full_results / stage / "stage_status.json").is_file()
    assert (ctx.full_results / "06_FINAL_HOTSPOTS" / "handoff_02.json").is_file()


def test_only_allowlisted_columns_reach_the_primary_path(clustered_run):
    """Stage A emits `aa`, `plddt_band` and CB columns; none may enter the analysis."""
    ctx, h1, h2, _ = clustered_run
    load_report = json.loads(
        (ctx.full_results / "12_REPRODUCIBILITY" / "provenance" /
         "provenance_06_final_hotspots.json").read_text())["inputs"]["load_report"]

    permitted = set(h1.payload["stage_b_permitted_columns"])
    for table, report in load_report.items():
        assert set(report["columns_loaded"]) <= permitted, table
        assert not set(report["columns_loaded"]) & set(FORBIDDEN_DOWNSTREAM_COLUMNS)
    # Stage A really did offer more than the allowlist, so the projection is doing work.
    assert load_report["residue_coordinates"]["columns_discarded_not_allowlisted"]


def test_review_star_channel_is_separate_from_the_primary_path(clustered_run):
    """>=1*/>=2* read stars only through the declared channel — F12."""
    ctx, h1, h2, _ = clustered_run
    channel = h1.payload["sensitivity_channel"]
    assert channel["declared"] is True
    assert channel["may_redefine_primary"] is False
    assert set(channel["columns"]) & set(FORBIDDEN_DOWNSTREAM_COLUMNS)

    path = ctx.full_results / "11_SENSITIVITY" / "sensitivity_review_status.json"
    if not path.is_file():
        # The run ended in a negative before r_hot was frozen, so no sensitivity
        # analysis is applicable. Absence must still be explicit (P3).
        assert (ctx.full_results / "11_SENSITIVITY" / "NOT_RUN.txt").is_file()
        return
    review = json.loads(path.read_text())
    assert set(review["strata"]) == {"star1", "star2"}
    assert review["may_redefine_primary"] is False
    for stratum in review["strata"].values():
        assert stratum["status"] in ("COMPLETED", "NOT_EVALUABLE", "NOT_APPLICABLE")
    if review["channel"]["available"]:
        # The written `channel` field is the whole return of open_review_channel:
        # availability and the star payload at the top level, the channel DESCRIPTOR
        # nested one level below it. This assertion read the descriptor's key off the
        # outer dict, so it could only ever have raised — it was unreachable for as
        # long as the fixture left the channel unavailable, and became reachable when
        # the clustered fixture started emitting variants_residue_level.tsv with the
        # declared `residue_index`/`max_star` columns.
        assert review["channel"]["channel"]["max_star_column"] \
            in FORBIDDEN_DOWNSTREAM_COLUMNS


def test_dispersed_cohort_reaches_an_empty_center_set(dispersed_run):
    """The fixture's `no_significant_hotspot` case must NOT produce a hotspot.

    v2 §5.4 adds a precondition: this may be reported as a negative only because the
    power certificate passed. If it had failed, the same empty center set would be an
    UNDERPOWERED run and no claim about the gene would be permitted.
    """
    ctx, h1, h2, case = dispersed_run
    assert case.expected["expected_negative_branch"] == "S_empty"
    assert h2.payload["terminal_state"] == "COMPLETED_NEGATIVE"
    assert h2.payload["n_significant_centers"] == 0
    assert h2.negative_result is not None
    assert h2.negative_result["licenses_no_method_modification"] is True
    assert h2.qc_status != "FAIL"                  # a negative is not a technical failure

    assert h2.payload["power_certificate_passed"] is True
    certificate = json.loads((ctx.full_results / "06_FINAL_HOTSPOTS" /
                              "power_certificate.json").read_text())
    assert certificate["PASSES"] is True
    assert certificate["posthoc"]["p_floor"] <= \
        certificate["posthoc"]["c_1_rank1_critical_value"]

    # Whichever branch fires, the chain stops and nothing is invented.
    assert h2.negative_result["condition"] in ("NO_ADMISSIBLE_RADIUS",
                                               "NO_SIGNIFICANT_HOTSPOT_CENTERS")
    centers = read_tsv(ctx.full_results / "06_FINAL_HOTSPOTS" /
                       "significant_hotspot_centers.tsv") \
        if (ctx.full_results / "06_FINAL_HOTSPOTS" /
            "significant_hotspot_centers.tsv").is_file() else []
    assert centers == []


def test_large_family_trips_the_permutation_resolution_diagnostic(resolution_run):
    """The fixture's `permutation_resolution` case: m so large that q/m < 1/(B+1).

    This case also pins the relationship between the two diagnostics, which measure
    the same floor against *different* family sizes and may therefore legitimately
    disagree:

      * II.4's condition C1 is evaluated pre-flight against ``m = |U_struct|`` — every
        residue is a candidate center;
      * the v2 §5.4 certificate is evaluated against the **realised** family of §5.1 —
        only centers whose sphere holds a classified residue are tested.

    The realised family is the smaller of the two, so ``c_1 = q/m`` is larger and the
    certificate can pass while C1 fires. The run then completes, carrying the
    resolution caveat; it is not underpowered, because the test it actually ran could
    have rejected. Neither diagnostic is allowed to silence the other.
    """
    ctx, h1, h2, case = resolution_run
    diagnostic = json.loads((ctx.full_results / "06_FINAL_HOTSPOTS" /
                             "permutation_resolution_diagnostic.json").read_text())

    assert diagnostic["PERMUTATION_RESOLUTION_LIMITED"] is True
    assert h2.payload["permutation_resolution_limited"] is True
    assert diagnostic["preflight"]["conditions"]["C1"] is True
    assert diagnostic["preflight"]["p_res"] == pytest.approx(1 / 10_001)
    assert diagnostic["preflight"]["m_family_size"] == h1.payload["M"]

    certificate = json.loads((ctx.full_results / "06_FINAL_HOTSPOTS" /
                              "power_certificate.json").read_text())
    posthoc = certificate["posthoc"]
    assert posthoc["m_test_family_size"] <= h1.payload["M"]      # §5.1 restriction
    assert posthoc["p_res"] == pytest.approx(1 / 10_001)
    assert posthoc["binding_floor"] == "permutation_resolution"
    assert posthoc["PASSES"] is (posthoc["p_floor"] <=
                                 posthoc["c_1_rank1_critical_value"])
    assert certificate["PASSES"] is posthoc["PASSES"]

    # Whatever the verdict, the two are consistent with each other and with the run.
    if certificate["PASSES"]:
        assert h2.payload["terminal_state"] in ("COMPLETED", "COMPLETED_NEGATIVE")
        assert h2.payload["power_certificate_passed"] is True
    else:
        assert h2.payload["terminal_state"] == "UNDERPOWERED"
        assert h2.negative_result is None
        assert h2.qc_status == "FAIL"
        assert certificate["TEST_CANNOT_REJECT"] is True
    # B and the FDR method are untouched: the diagnostic reports, it never acts.
    assert diagnostic["B_used"] == 10_000
    assert diagnostic["B_was_changed_by_the_agent"] is False
    assert diagnostic["fdr_method_was_changed_by_the_agent"] is False
    assert diagnostic["B_recommended"] > 10_000

    # Five P/LP residues also trip FT-4 (N_P < N_MIN_PCF), and the fallback is surfaced.
    assert h2.payload["fallback_radius_domain"] is True
    domain = json.loads((ctx.full_results / "04_GLOBAL_CLUSTERING" /
                         "candidate_radius_domain.json").read_text())
    assert "FT-4" in domain["triggers_fired"]
    assert h2.payload["fallback_trigger_id"] == domain["triggers_fired"][0]
    warnings = read_tsv(ctx.full_results / "06_FINAL_HOTSPOTS" /
                        "warnings_06_FINAL_HOTSPOTS.tsv")
    assert any(w["warning_code"] == "PERMUTATION_RESOLUTION_LIMITED" for w in warnings)
    assert any(e["ambiguity"].startswith("PERMUTATION_RESOLUTION_LIMITED")
               for e in h2.payload["escalations"])


def test_low_plddt_never_filters_the_primary_analysis(low_plddt_run):
    """F2: pLDDT is recorded, never a primary filter; it appears only in 11_SENSITIVITY."""
    ctx, h1, h2, case = low_plddt_run
    cohort = read_tsv(ctx.full_results / "03_STRUCTURE_QC" / "classified_cohort.tsv")
    coords = {r["residue_index"]: r for r in
              read_tsv(ctx.full_results / "03_STRUCTURE_QC" / "residue_coordinates.tsv")}

    low = [r for r in cohort if coords[r["residue_index"]]["plddt"] < 70]
    assert low, "the fixture is supposed to carry a low-confidence tail"

    universe = read_tsv(ctx.full_results / "03_STRUCTURE_QC" / "positional_universe.tsv")
    summary = json.loads((ctx.full_results / "04_GLOBAL_CLUSTERING" /
                          "ripleys_k_summary.json").read_text())
    # Every residue of U_struct entered the global analysis; nothing was filtered out.
    assert summary["global_test"]["PLP"]["n_points"] + \
        summary["global_test"]["BLB"]["n_points"] == len(cohort)
    assert len(universe) == h1.payload["M"]

    path = ctx.full_results / "11_SENSITIVITY" / "sensitivity_plddt70.json"
    if not path.is_file():
        assert (ctx.full_results / "11_SENSITIVITY" / "NOT_RUN.txt").is_file()
        return
    plddt = json.loads(path.read_text())
    assert plddt["IS_NON_REDEFINING"] is True
    if plddt["status"] == "COMPLETED":
        assert plddt["n_centers_universe"] < h1.payload["M"]     # the filter bit only here


def test_run_id_changes_the_derived_seeds_but_not_the_contract(clustered_run,
                                                                alternate_run_id_run):
    """Seeds are context-derived AND run-scoped: a new RUN_ID is a genuinely new draw."""
    a_ctx, _, a, _ = clustered_run
    b_ctx, _, b, _ = alternate_run_id_run

    shared = set(a.payload["derived_seeds"]) & set(b.payload["derived_seeds"])
    assert {"final_detection", "global_null|PLP", "global_null|BLB"} <= shared
    # Same context strings, different RUN_ID -> genuinely different draws.
    assert all(a.payload["derived_seeds"][c] != b.payload["derived_seeds"][c]
               for c in shared)
    for handoff in (a, b):
        assert handoff.payload["B"] == 10_000
        assert handoff.payload["fdr_method"] == "BH"
        for stage in STAGES:
            for rel in EXPECTED_OUTPUTS[stage]:
                root = (a_ctx if handoff is a else b_ctx).full_results
                if (root / stage / "NOT_RUN.txt").is_file() and rel.startswith("figures/"):
                    continue                       # a stopped stage writes no figures
                assert (root / stage / rel).is_file() or \
                    (root / stage / "NOT_RUN.txt").is_file(), f"{stage}/{rel}"
