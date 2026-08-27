"""Workflow v2 §2/§4 — the U_center split and the biological-plausibility diagnostics.

v2 §2 excludes low-pLDDT residues from the CANDIDATE-CENTER universe U_center while
leaving U_struct (the positional reference) and the classified cohort L untouched.
v2 §4 requires structural_coverage/cohort_absorption/r_to_domain_ratio to be computed
and reported, but F15 freezes the Pareto objective count at exactly four (Lead ruling,
config/pipeline.yaml) — so these three quantities must never move a selection result.

These tests pin:

  * U_center = {r in U_struct : pLDDT >= plddt.center_universe_min_plddt};
    |U_center| < |U_struct| when a low-confidence tail exists, |U_center| ==
    |U_struct| when it does not;
  * the positional null still places P/LP over the FULL U_struct even when the
    candidate-center row space is restricted (n_centers != M);
  * structural_coverage/cohort_absorption/r_to_domain_ratio are computed correctly
    and are reported, never admissibility constraints, never Pareto objectives;
  * MAJOR/ADVISORY warnings fire exactly at the documented thresholds;
  * selection is IDENTICAL whether or not the diagnostics are computed;
  * the v2 §6 one-page decision summary states its rows in the documented order.
"""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from synthetic_hotspot import full_synthetic_run, planted_cohort, synthetic_protein

from hotspot3d.hotspot.detection import final_detection
from hotspot3d.hotspot.positional import positional_pass
from hotspot3d.hotspot.report import _decision_summary
from hotspot3d.hotspot.scan import RadiusResult, scan_one_radius
from hotspot3d.hotspot.selection import run_selection
from hotspot3d.hotspot.stage import run_stage_b
from hotspot3d.utils.config import load_config
from hotspot3d.utils.io import read_tsv

pytestmark = pytest.mark.unit

SELECT = dict(tie_chain=["distance", "neighbor_stability", "perm_evidence_zg",
                         "smaller_radius"],
              w_k=1.0, metric="L2", near_tie_rel=0.05, near_tie_sep_steps=2, step=0.5,
              band_fraction=0.10, bw2_fraction=0.50, bw3_top_n=3, redundancy_rho=0.98)


# --- §2: U_center construction, low-level -------------------------------------

def _grid_inputs(n_side=5, spacing=6.0, seed=3):
    grid = np.arange(n_side) * spacing
    coords = np.array([[x, y, z] for x in grid for y in grid for z in grid],
                      dtype=np.float64)
    index = np.arange(1, len(coords) + 1)
    rng = np.random.default_rng(seed)
    positions = np.sort(rng.choice(len(coords), size=16, replace=False))
    y = np.array([1.0] * 8 + [0.0] * 8)
    return coords, index, positions, coords[positions], y


def test_positional_pass_decouples_candidate_count_from_population():
    """U need not be square: n_centers (rows) and M = |U_struct| (columns) differ."""
    coords, index, positions, labeled, y = _grid_inputs()
    M = len(coords)
    center_mask = np.zeros(M, dtype=bool)
    center_mask[:10] = True                     # only the first 10 residues are U_center
    n_centers = int(center_mask.sum())

    from hotspot3d.utils.geometry import cross_distances, membership_matrix
    radius = 8.0
    U = (cross_distances(coords[center_mask], coords) <= radius).astype(np.uint8)
    A = membership_matrix(coords[center_mask], labeled, radius)
    plp_positions = positions[y > 0]

    out = positional_pass(U, plp_positions, A.sum(axis=1), np.random.default_rng(1),
                          200, radius, "test|u_center_decouple")

    assert U.shape == (n_centers, M)
    assert len(out.n_universe) == n_centers          # one row per CANDIDATE center
    assert out.M == M                                # population is the FULL U_struct
    # n_U(c) can never exceed the population it is drawn from.
    assert (out.n_universe <= M).all()


def test_scan_one_radius_restricts_centers_without_touching_placement():
    """centers/center_rows/S(r) all live in U_center; the null still draws over M."""
    coords, index, positions, labeled, y = _grid_inputs()
    M = len(coords)
    center_mask = np.zeros(M, dtype=bool)
    center_mask[:15] = True
    center_coords, center_index = coords[center_mask], index[center_mask]

    out_restricted = scan_one_radius(
        6.0, coords, index, positions, labeled, y, np.random.default_rng(1),
        label_rng=np.random.default_rng(2), B=200, q=0.05, kappa=2.0,
        tie_tolerance=1e-12, sparse_cutoff=3, undefined_mcc=0.0,
        qc_h2_max_coverage=0.50, qc_h3_min_median=2, qc_h4_max_isolated=0.50,
        r_floor=3.0, r_cap=30.0, seed_context="scan|r=6",
        center_coords=center_coords, center_index=center_index)

    # Every reported center residue index is drawn from U_center, never from the
    # excluded tail of U_struct.
    allowed = set(int(i) for i in center_index)
    assert set(out_restricted.centers) <= allowed
    for row in out_restricted.center_rows:
        assert row["center_residue_index"] in allowed

    # The positional null's population is unaffected: n_universe_in_sphere can
    # exceed n_centers (it counts ALL of U_struct within r, not just U_center).
    if out_restricted.center_rows:
        assert max(r["n_universe_in_sphere"] for r in out_restricted.center_rows) <= M


def test_final_detection_center_universe_default_matches_full_universe():
    """Omitting center_coords/center_index preserves pre-v2 (centers == U_struct)."""
    coords, index, positions, labeled, y = _grid_inputs()
    det = final_detection(
        6.0, coords, index, positions, labeled, y, np.random.default_rng(1),
        B=200, q=0.05, seed_context="final|default")
    assert len(det.center_rows) == len(coords)
    assert det.counts["M_universe"] == len(coords)
    assert det.counts["n_center_universe"] == len(coords)


def test_final_detection_restricted_center_universe_shrinks_the_family():
    """A strict U_center subset must never grow the family beyond |U_center|."""
    coords, index, positions, labeled, y = _grid_inputs()
    M = len(coords)
    center_mask = np.zeros(M, dtype=bool)
    center_mask[:20] = True
    center_coords, center_index = coords[center_mask], index[center_mask]

    det = final_detection(
        6.0, coords, index, positions, labeled, y, np.random.default_rng(1),
        B=200, q=0.05, seed_context="final|restricted",
        center_coords=center_coords, center_index=center_index)

    assert len(det.center_rows) == int(center_mask.sum())
    assert det.counts["M_universe"] == M
    assert det.counts["n_center_universe"] == int(center_mask.sum())
    assert det.counts["n_in_test_family"] <= int(center_mask.sum())
    # Every published residue index is a genuine U_center member.
    allowed = set(int(i) for i in center_index)
    for row in det.center_rows:
        assert row["center_residue_index"] in allowed
    for row in det.covered_rows:                       # coverage is over the FULL universe
        assert row["covered_residue_index"] in set(int(i) for i in index)


# --- §2: integration — the |U_center| < |U_struct| / == cases -----------------

@pytest.mark.integration
def test_u_center_is_strictly_smaller_with_a_low_confidence_tail(tmp_path):
    ctx, handoff, geo = full_synthetic_run(
        tmp_path, B=200, low_confidence_fraction=0.2)
    try:
        run_stage_b(ctx, upstream=handoff)
    except Exception:
        pass          # B=200 with the default geometry is expected to be UNDERPOWERED;
                       # center_universe.json is written before that termination.
    report = ctx.full_results / "05_HOTSPOT_RADIUS" / "center_universe.json"
    assert report.is_file()
    import json
    payload = json.loads(report.read_text())
    assert payload["n_center_universe_U_center"] < payload["n_universe_U_struct"]
    assert payload["n_excluded_from_center_universe"] == (
        payload["n_universe_U_struct"] - payload["n_center_universe_U_center"])
    assert payload["min_plddt_threshold"] == 50.0
    assert payload["n_contiguous_excluded_regions"] >= 1

    exclusions = read_tsv(ctx.full_results / "05_HOTSPOT_RADIUS" /
                          "center_universe_exclusions.tsv")
    assert len(exclusions) == payload["n_excluded_from_center_universe"]
    for row in exclusions:
        assert float(row["plddt"]) < 50.0


@pytest.mark.integration
def test_u_center_equals_u_struct_without_a_low_confidence_tail(tmp_path):
    ctx, handoff, geo = full_synthetic_run(tmp_path, B=200)          # no low_confidence
    try:
        run_stage_b(ctx, upstream=handoff)
    except Exception:
        pass
    import json
    payload = json.loads((ctx.full_results / "05_HOTSPOT_RADIUS" /
                          "center_universe.json").read_text())
    assert payload["n_center_universe_U_center"] == payload["n_universe_U_struct"]
    assert payload["n_excluded_from_center_universe"] == 0
    exclusions = read_tsv(ctx.full_results / "05_HOTSPOT_RADIUS" /
                          "center_universe_exclusions.tsv")
    assert exclusions == []


# --- §4: biological-plausibility formulas --------------------------------------

def test_structural_coverage_is_exactly_the_coverage_fraction():
    coords, index, positions, labeled, y = planted_cohort_inputs()
    out = scan_one_radius(
        8.0, coords, index, positions, labeled, y, np.random.default_rng(4),
        label_rng=np.random.default_rng(5), B=300, q=0.05, kappa=2.0,
        tie_tolerance=1e-12, sparse_cutoff=3, undefined_mcc=0.0,
        qc_h2_max_coverage=0.50, qc_h3_min_median=2, qc_h4_max_isolated=0.50,
        r_floor=3.0, r_cap=30.0, seed_context="scan|r=8",
        radius_of_gyration_A=20.0)
    assert out.structural_coverage == pytest.approx(out.coverage_fraction)


def test_cohort_absorption_matches_the_largest_significant_sphere():
    coords, index, positions, labeled, y = planted_cohort_inputs()
    out = scan_one_radius(
        10.0, coords, index, positions, labeled, y, np.random.default_rng(4),
        label_rng=np.random.default_rng(5), B=300, q=0.05, kappa=2.0,
        tie_tolerance=1e-12, sparse_cutoff=3, undefined_mcc=0.0,
        qc_h2_max_coverage=0.50, qc_h3_min_median=2, qc_h4_max_isolated=0.50,
        r_floor=3.0, r_cap=30.0, seed_context="scan|r=10",
        radius_of_gyration_A=20.0)
    N = len(y)
    if out.center_rows:
        expected = max(r["n_labeled_in_sphere"] for r in out.center_rows) / N
    else:
        expected = 0.0
    assert out.cohort_absorption == pytest.approx(expected)
    assert 0.0 <= out.cohort_absorption <= 1.0


def test_r_to_domain_ratio_is_radius_over_the_supplied_gyration_radius():
    coords, index, positions, labeled, y = planted_cohort_inputs()
    out = scan_one_radius(
        7.5, coords, index, positions, labeled, y, np.random.default_rng(4),
        label_rng=np.random.default_rng(5), B=300, q=0.05, kappa=2.0,
        tie_tolerance=1e-12, sparse_cutoff=3, undefined_mcc=0.0,
        qc_h2_max_coverage=0.50, qc_h3_min_median=2, qc_h4_max_isolated=0.50,
        r_floor=3.0, r_cap=30.0, seed_context="scan|r=7.5",
        radius_of_gyration_A=15.0)
    assert out.r_to_domain_ratio == pytest.approx(7.5 / 15.0)


def test_r_to_domain_ratio_is_nan_without_a_supplied_gyration_radius():
    coords, index, positions, labeled, y = planted_cohort_inputs()
    out = scan_one_radius(
        7.5, coords, index, positions, labeled, y, np.random.default_rng(4),
        label_rng=np.random.default_rng(5), B=300, q=0.05, kappa=2.0,
        tie_tolerance=1e-12, sparse_cutoff=3, undefined_mcc=0.0,
        qc_h2_max_coverage=0.50, qc_h3_min_median=2, qc_h4_max_isolated=0.50,
        r_floor=3.0, r_cap=30.0, seed_context="scan|r=7.5")
    assert np.isnan(out.r_to_domain_ratio)


def planted_cohort_inputs():
    coords = synthetic_protein(sphere_radius=18.0, spacing=5.0)
    positions, labels = planted_cohort(coords, n_cluster_plp=16, n_scatter_plp=8, n_blb=40)
    index = np.arange(1, len(coords) + 1)
    return coords, index, positions, coords[positions], labels.astype(np.float64)


def test_config_thresholds_match_the_documented_values():
    """A drift here would silently change which warning fires where (v2 §4)."""
    cfg = load_config("config/pipeline.yaml")
    assert cfg.get("biological_plausibility.structural_coverage_penalty_threshold") == 0.25
    assert cfg.get("biological_plausibility.structural_coverage_major_threshold") == 0.50
    assert cfg.get("biological_plausibility.cohort_absorption_penalty_threshold") == 0.33
    assert cfg.get("biological_plausibility.r_to_domain_ratio_penalty_threshold") == 0.5
    assert cfg.get("biological_plausibility.enters_pareto_objectives") is False


# --- §4: never a Pareto objective, never an admissibility constraint -----------

def _result(radius, mcc, zg, fe, stab, structural_coverage=0.0, cohort_absorption=0.0,
           r_to_domain_ratio=0.0) -> RadiusResult:
    return RadiusResult(
        radius_A=radius, centers=[], center_rows=[],
        loo=SimpleNamespace(mcc=mcc), perm=SimpleNamespace(zg=zg),
        positional=SimpleNamespace(null_model="structure_aware_positional"),
        fold_enrichment=fe, fold_enrichment_defined=True,
        n_plp_in=0, n_labeled_in=0, coverage_fraction=structural_coverage,
        n_singleton_centers=0,
        structural_coverage=structural_coverage, cohort_absorption=cohort_absorption,
        r_to_domain_ratio=r_to_domain_ratio,
        median_n_labeled=3.0, qc_failures=[], diagnostics_fired=[], diagnostic_values={},
        boundary_neutral_verified=True, bh_boundary_p=None, n_in_family=10,
        neighbor_stability=stab)


def test_selection_is_identical_whether_or_not_bio_plausibility_is_computed():
    """The three diagnostics never enter Candidate.objectives() (F15 stays at four)."""
    live = [
        _result(5.0, 0.2, 1.0, 2.0, 0.1, structural_coverage=0.9, cohort_absorption=0.8,
               r_to_domain_ratio=1.5),
        _result(5.5, 0.6, 3.0, 1.0, 0.5, structural_coverage=0.4, cohort_absorption=0.2,
               r_to_domain_ratio=0.6),
        _result(6.0, 0.9, 2.0, 1.5, 0.9, structural_coverage=0.1, cohort_absorption=0.05,
               r_to_domain_ratio=0.2),
    ]
    stubbed = [
        _result(r.radius_A, r.loo.mcc, r.perm.zg, r.fold_enrichment, r.neighbor_stability,
               structural_coverage=0.0, cohort_absorption=0.0, r_to_domain_ratio=0.0)
        for r in live
    ]

    sel_live = run_selection(live, [5.0, 5.5, 6.0], **SELECT)
    sel_stub = run_selection(stubbed, [5.0, 5.5, 6.0], **SELECT)

    assert sel_live.selected == sel_stub.selected
    assert sel_live.result.pareto_members == sel_stub.result.pareto_members
    assert sel_live.result.distances == sel_stub.result.distances
    assert sel_live.result.normalized == sel_stub.result.normalized


def test_bio_plausibility_never_appears_in_the_pareto_objective_names():
    from hotspot3d.hotspot.scan import OBJECTIVE_NAMES
    assert OBJECTIVE_NAMES == ["loo_mcc", "perm_evidence_zg", "fold_enrichment",
                               "neighbor_stability"]
    assert "structural_coverage" not in OBJECTIVE_NAMES
    assert "cohort_absorption" not in OBJECTIVE_NAMES
    assert "r_to_domain_ratio" not in OBJECTIVE_NAMES


def test_bio_plausibility_never_appears_in_the_admissibility_rules():
    from hotspot3d.hotspot.scan import ADMISSIBILITY_RULES, SIGNIFICANCE_CONDITIONED_DIAGNOSTICS
    for rules in (ADMISSIBILITY_RULES, SIGNIFICANCE_CONDITIONED_DIAGNOSTICS):
        assert "structural_coverage" not in rules
        assert "cohort_absorption" not in rules
        assert "r_to_domain_ratio" not in rules


# --- §4: warnings fire exactly at the documented thresholds, integration -------

@pytest.mark.integration
@pytest.mark.slow
def test_bio_plausibility_warnings_are_self_consistent_with_their_thresholds(tmp_path):
    """Whichever way the fixture's numbers land, the warning must fire IFF the
    threshold is crossed — checked against the SELECTED radius's own reported row,
    so the test is not tied to engineering a fixture that crosses a specific bound.
    """
    ctx, handoff, geo = full_synthetic_run(tmp_path, B=10_000)
    handoff2 = run_stage_b(ctx, upstream=handoff)
    if handoff2.payload.get("terminal_state") == "UNDERPOWERED":
        pytest.skip("fixture went UNDERPOWERED before select_radius(); nothing to check")

    rows = read_tsv(ctx.full_results / "05_HOTSPOT_RADIUS" / "r_hot_scan.tsv")
    selected = [r for r in rows if r["selected"] is True]
    assert len(selected) == 1
    row = selected[0]
    for col in ("structural_coverage", "cohort_absorption", "r_to_domain_ratio"):
        assert col in row and row[col] is not None

    cfg = load_config("config/pipeline.yaml")
    major = float(cfg.get("biological_plausibility.structural_coverage_major_threshold"))
    absorb = float(cfg.get("biological_plausibility.cohort_absorption_penalty_threshold"))
    ratio = float(cfg.get("biological_plausibility.r_to_domain_ratio_penalty_threshold"))

    warnings = read_tsv(ctx.full_results / "05_HOTSPOT_RADIUS" /
                        "warnings_05_HOTSPOT_RADIUS.tsv")
    codes = {w["warning_code"] for w in warnings}

    assert ("EXCESSIVE_STRUCTURAL_COVERAGE" in codes) == (float(row["structural_coverage"]) > major)
    assert ("HIGH_COHORT_ABSORPTION" in codes) == (float(row["cohort_absorption"]) > absorb)
    assert ("HIGH_R_TO_DOMAIN_RATIO" in codes) == (float(row["r_to_domain_ratio"]) > ratio)

    # These never remove a radius from the Pareto set or change which one is selected.
    decision_path = ctx.full_results / "05_HOTSPOT_RADIUS" / "radius_decision.json"
    import json
    decision = json.loads(decision_path.read_text())
    assert decision["selected_r_hot_A"] == float(row["radius_A"])


# --- §6: the one-page decision summary, content and order ----------------------

def test_decision_summary_states_rows_in_the_v2_section_6_order():
    """terminal state; N_P; N_B; |U_struct|; |U_center|; m; null; p_floor; c_1;
    ratio; ...; selected r_hot; vacuous; n_significant — in THAT order (v2 §6)."""
    s = {
        "inputs": {"N_P": 12, "N_B": 40, "M": 500, "n_center_universe": 470},
        "detection": {"n_in_test_family": 88, "n_significant_bh": 3},
        "selection": {"r_hot": 9.5, "vacuous_pareto": False},
        "primary_null": "structure_aware_positional",
        "power_certificate": {
            "posthoc": {
                "n_center_universe_U_center": 470, "m_test_family_size": 88,
                "null_model": "structure_aware_positional",
                "p_floor": 1.2e-6, "c_1_rank1_critical_value": 5.7e-4,
                "p_floor_over_c_1": 2.1e-3, "binding_floor": "permutation_resolution",
                "PASSES": True,
            },
        },
    }
    lines = _decision_summary(s, "COMPLETED")
    table_block = next(l for l in lines if l.startswith("| quantity | value |"))
    table_lines = table_block.split("\n")
    data_lines = table_lines[2:2 + 15]   # skip header + separator
    # rsplit, not split: some labels (e.g. "|U_struct|") embed pipe characters.
    rows = [line[2:-2].rsplit(" | ", 1)[0] for line in data_lines]

    expected_order = [
        "terminal state", "N_P (P/LP residues)", "N_B (B/LB residues)",
        "|U_struct| (positional reference)",
        "|U_center| (test family before restriction)",
        "m (family actually tested)",
        "null used for the PRIMARY per-center test",
        "p_floor", "c_1 = q/m", "p_floor / c_1", "binding floor", "power certificate",
        "selected r_hot (A)", "Pareto set vacuous", "significant hotspot centers",
    ]
    assert rows == expected_order


def test_decision_summary_reads_u_center_from_inputs_when_no_certificate_yet():
    """Before any certificate exists (e.g. an early BLOCKED path), |U_center| must
    still read from inputs rather than silently falling back to |U_struct|."""
    s = {"inputs": {"N_P": 5, "N_B": 5, "M": 100, "n_center_universe": 80},
         "detection": {}, "selection": None, "power_certificate": {}}
    lines = _decision_summary(s, "BLOCKED")
    u_center_line = next(l for l in lines if "U_center" in l)
    assert "| 80 |" in u_center_line
