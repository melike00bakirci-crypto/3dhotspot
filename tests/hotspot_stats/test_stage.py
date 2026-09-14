"""Stage B end to end: contracts, barriers, negative branches and determinism."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from hotspot3d.hotspot.stage import (
    EXPECTED_OUTPUTS,
    STAGES,
    is_underpowered,
    run_stage_b,
)
from hotspot3d.orchestration.contracts import (
    HANDOFF_REQUIRED_KEYS,
    assert_handoff_shape,
)
from hotspot3d.utils.config import assert_frozen_methodology, load_config
from hotspot3d.utils.errors import BlockedError, LeakageError
from hotspot3d.utils.hashing import sha256_file, verify_manifest
from hotspot3d.utils.io import read_tsv, write_tsv
from synthetic_hotspot import full_synthetic_run

FULL = "FULL_RESULTS"
# Files that legitimately carry a wall-clock timestamp, a wall time or an environment
# fingerprint. Everything else must be byte-identical; the warning files are compared
# separately with the detection timestamp stripped, so their CONTENT is still checked.
TIMESTAMPED = {"stage_status.json", "stage_manifest.tsv", "stage_b_report.md"}


def _is_timestamped(path: Path) -> bool:
    # handoff_02.json embeds the manifest hashes of the timestamped files above, so it
    # inherits their variability; its scientific payload is compared separately.
    return (path.name in TIMESTAMPED or path.name.startswith("warnings_")
            or path.name == "handoff_02.json")


def _payload_without_timestamps(handoff) -> dict:
    payload = dict(handoff.payload)
    payload["escalations"] = [{k: v for k, v in e.items() if k != "detected_utc"}
                              for e in payload.get("escalations", [])]
    # code_version hashes the whole source tree, which other stages are still growing.
    # It pins the run; it is not an output of the computation.
    payload.pop("code_version", None)
    return payload


def _warnings_without_timestamps(root: Path) -> dict:
    out = {}
    for path in sorted(root.rglob("warnings_*.tsv")):
        rows = read_tsv(path)
        for row in rows:
            row.pop("detected_utc", None)
        out[str(path.relative_to(root))] = rows
    return out


def _positive_run(tmp_path, **kwargs):
    ctx, handoff, geo = full_synthetic_run(tmp_path, B=10_000, **kwargs)
    return ctx, run_stage_b(ctx, upstream=handoff), geo


# --- shared stage runs -------------------------------------------------------
# Executing run_stage_b is the expensive part of this module (~16 s at the frozen
# B = 10000, dominated by the global permutation null and matplotlib rendering).
# The contract assertions below all interrogate the SAME completed run, so the run
# is performed once per module and shared between them.
#
# This is a TEST-EXECUTION optimization only. No assertion is weakened, removed or
# relaxed: every test still checks exactly what it checked before, against a run
# produced by the frozen production configuration (B = 10000, q = 0.05, BH).

@pytest.fixture(scope="module")
def shared_positive(tmp_path_factory):
    """One complete Stage B run at the FROZEN default B = 10000."""
    return _positive_run(tmp_path_factory.mktemp("positive"))


@pytest.fixture(scope="module")
def determinism_pair(tmp_path_factory):
    """Two independent runs with identical inputs and identical derived seeds."""
    root = tmp_path_factory.mktemp("determinism")
    a_ctx, a_handoff, _ = _positive_run(root / "a")
    b_ctx, b_handoff, _ = _positive_run(root / "b")
    return (a_ctx, a_handoff), (b_ctx, b_handoff)


@pytest.fixture(scope="module")
def shared_underpowered(tmp_path_factory):
    """One reduced-B run that terminates UNDERPOWERED (Workflow v2 §5.4).

    B = 200 here is NOT a performance shortcut and is not interchangeable with the
    production default: at this family size the permutation resolution floor
    ``1/(B+1)`` sits above the rank-1 BH critical value ``q/m``, so no center could
    be declared significant regardless of the data. Under v1 this run was reported
    as ``NO_ADMISSIBLE_RADIUS`` — a scientific negative — which is precisely the
    misreport v2 §5.4 exists to prevent.

    Decoy 07_*/08_*/09_*/10_* directories are planted before the run, so every
    consumer of this fixture also exercises the information barrier.
    """
    ctx, handoff, geo = full_synthetic_run(tmp_path_factory.mktemp("underpowered"), B=200)
    for stage in ("07_FOOTPRINT_RADIUS", "08_FINAL_FOOTPRINT", "09_ROBUSTNESS",
                  "10_ANNOTATION"):
        d = ctx.full_results / stage
        d.mkdir(parents=True, exist_ok=True)
        (d / "tempting.tsv").write_text("schema_version\tr_fp\n1.0.0\t9.5\n")
    return ctx, run_stage_b(ctx, upstream=handoff), geo


@pytest.fixture(scope="module")
def shared_certified_negative(tmp_path_factory):
    """A CERTIFIED negative: the design could have rejected, and did not.

    The cohort is drawn uniformly over the protein — exactly the distribution the
    structure-aware positional null assumes — at the FROZEN B = 10000, so the §5.4
    certificate passes and an empty center set is a statement about the data rather
    than about the design. This is the only configuration in which
    ``COMPLETED_NEGATIVE`` is a legitimate terminal state.
    """
    ctx, handoff, geo = full_synthetic_run(
        tmp_path_factory.mktemp("certified_negative"), B=10_000, dispersed=True)
    return ctx, run_stage_b(ctx, upstream=handoff), geo



# --- happy path --------------------------------------------------------------

@pytest.mark.integration
@pytest.mark.slow
def test_end_to_end_at_the_frozen_default_b(shared_positive):
    """One full run at the FROZEN default B = 10000, start to finish."""
    ctx, handoff, _ = shared_positive
    p = handoff.payload

    assert handoff.name == "handoff_02"
    assert handoff.qc_status in ("PASS", "PASS_WITH_WARNINGS")
    assert p["B"] == 10_000 and p["q"] == 0.05 and p["fdr_method"] == "BH"
    assert p["n_significant_centers"] > 0
    assert p["hotspot_radius"] is not None
    assert_handoff_shape(handoff)
    for key in HANDOFF_REQUIRED_KEYS["handoff_02"]:
        assert key in p, key

    # Every canonical output exists in every owned stage directory.
    for stage in STAGES:
        for rel in EXPECTED_OUTPUTS[stage]:
            assert (ctx.full_results / stage / rel).is_file(), f"{stage}/{rel}"

    # Provenance: exactly our own four files, nothing else.
    prov = sorted(p.name for p in
                  (ctx.full_results / "12_REPRODUCIBILITY" / "provenance").iterdir())
    assert prov == ["provenance_04_global_clustering.json",
                    "provenance_05_hotspot_radius.json",
                    "provenance_06_final_hotspots.json",
                    "provenance_11_sensitivity.json"]

    # The handoff manifest verifies against what is actually on disk.
    assert verify_manifest(handoff.manifest, ctx.run_root) == []


@pytest.mark.integration
@pytest.mark.slow
def test_stage_writes_nothing_outside_its_owned_directories(shared_positive):
    ctx, handoff, _ = shared_positive
    written = {p.relative_to(ctx.full_results).parts[0]
               for p in ctx.full_results.rglob("*") if p.is_file()}
    allowed = set(STAGES) | {"01_INPUT_RAW", "02_CLINVAR", "03_STRUCTURE_QC",
                             "12_REPRODUCIBILITY"}
    assert written <= allowed
    for forbidden in ("07_FOOTPRINT_RADIUS", "08_FINAL_FOOTPRINT", "09_ROBUSTNESS",
                      "10_ANNOTATION", "00_RUN_SUMMARY"):
        assert not (ctx.full_results / forbidden).exists()
    assert not (ctx.run_root / "REVIEW_PACK").exists()


@pytest.mark.integration
@pytest.mark.slow
def test_the_three_residue_objects_are_separate_and_nested_on_disk(shared_positive):
    ctx, handoff, _ = shared_positive
    final = ctx.full_results / "06_FINAL_HOTSPOTS"
    centers = read_tsv(final / "significant_hotspot_centers.tsv")
    covered = read_tsv(final / "hotspot_covered_residues.tsv")
    classified = read_tsv(final / "hotspot_classified_variants.tsv")
    every = read_tsv(final / "all_residue_center_tests.tsv")

    c_idx = {r["center_residue_index"] for r in centers}
    v_idx = {r["variant_residue_index"] for r in classified}
    o_idx = {r["covered_residue_index"] for r in covered}
    assert c_idx <= o_idx and v_idx <= o_idx

    # QC rule 14: the published center file is EXACTLY the significant=TRUE subset,
    # with identical columns, so the two files can never disagree.
    assert {r["center_residue_index"] for r in every if r["significant"]} == c_idx
    assert set(centers[0].keys()) == set(every[0].keys())
    assert all(r["q_bh"] <= 0.05 for r in centers)
    assert all(r["p_emp"] <= handoff.payload["bh_boundary_p"] for r in centers)


@pytest.mark.integration
@pytest.mark.slow
def test_scan_trace_retains_every_radius_with_a_verdict(shared_positive):
    ctx, handoff, _ = shared_positive
    radius_dir = ctx.full_results / "05_HOTSPOT_RADIUS"
    rows = read_tsv(radius_dir / "r_hot_scan.tsv")
    domain = json.loads((ctx.full_results / "04_GLOBAL_CLUSTERING" /
                         "candidate_radius_domain.json").read_text())

    assert len(rows) == domain["n_grid_points"]
    assert [r["radius_A"] for r in rows] == domain["grid_A"]
    admissibility_rules = set(ctx.config.get("radius_qc.admissibility_constraints"))
    for row in rows:                                    # no gaps anywhere
        assert row["qc_status"] in ("PASS", "FAIL")
        assert row["admissible"] is (row["qc_status"] == "PASS")
        assert row["primary_null"] == "structure_aware_positional"
        assert row["secondary_null"] == "label_permutation"
        if row["qc_status"] == "PASS":
            assert row["qc_failure_reason"] is None     # written as "NA", never blank
        if row["qc_status"] == "FAIL":
            assert isinstance(row["qc_failure_reason"], str)
            # v2 §4 — only DECLARED, genuinely-undefined quantities may do this.
            assert set(row["qc_failure_reason"].split(";")) <= admissibility_rules
            assert isinstance(row["admissibility_reason"], str)
            assert row["pareto_member"] is None         # NA, not silently blank
            assert row["distance_to_ideal"] is None
        # v2 §4 — a radius that found nothing keeps its place in the optimization.
        if row["n_significant_centers"] == 0:
            assert row["admissible"] is True
            assert row["qc_h1_zero_significant_centers"] is True
    assert sum(1 for r in rows if r["selected"]) == 1

    # Provisional center sets survive for EVERY radius, selected or not.
    for row in rows:
        path = radius_dir / "radius_scan_centers" / f"r_{row['radius_A']:.1f}.tsv"
        assert path.is_file()
        assert len(read_tsv(path)) == row["n_significant_centers"]

    dominated = json.loads((radius_dir / "pareto_dominated.json").read_text())
    inadmissible = {float(k) for k in dominated["inadmissible"]}
    assert inadmissible == {r["radius_A"] for r in rows if r["qc_status"] == "FAIL"}


@pytest.mark.integration
@pytest.mark.slow
def test_selection_is_recorded_as_pareto_not_argmax_mcc(shared_positive):
    ctx, handoff, _ = shared_positive
    decision = json.loads((ctx.full_results / "05_HOTSPOT_RADIUS" /
                           "radius_decision.json").read_text())
    assert decision["selected_r_hot_A"] == handoff.payload["hotspot_radius"]
    assert decision["weights"] == {k: 1.0 for k in decision["objectives"]}
    assert decision["weights_frozen_equal"] is True
    assert decision["normalization"] == "min_max"
    assert decision["normalization_domain"] == "admissible_set"
    assert decision["utopia_point"] == {k: 1.0 for k in decision["objectives"]}
    assert decision["distance_metric"] == "L2"
    assert decision["max_mcc_radius_has_no_privileged_status"] is True
    assert decision["independent_dominance_check_passed"] is True
    # Every candidate is retained with a reason; only the selected one has none.
    assert len(decision["candidates"]) == decision["n_candidates"]
    assert sum(1 for c in decision["candidates"] if c["rejection_reason"] == "NA") == 1


@pytest.mark.integration
@pytest.mark.slow
def test_sensitivity_outputs_are_labelled_non_redefining(shared_positive):
    ctx, handoff, _ = shared_positive
    sens = ctx.full_results / "11_SENSITIVITY"
    statement = (sens / "SENSITIVITY_IS_NON_REDEFINING.txt").read_text()
    assert "NON-REDEFINING" in statement
    assert "never resolved by adopting the sensitivity result" in statement
    plddt = json.loads((sens / "sensitivity_plddt70.json").read_text())
    review = json.loads((sens / "sensitivity_review_status.json").read_text())
    assert plddt["IS_NON_REDEFINING"] is True and plddt["may_redefine_primary"] is False
    assert plddt["run_at_frozen_r_hot"] == handoff.payload["hotspot_radius"]
    assert set(review["strata"]) == {"star1", "star2"}
    for row in read_tsv(sens / "sensitivity_overlap.tsv"):
        assert row["is_non_redefining"] is True

    # The primary center set is untouched by whatever the sensitivity runs found.
    centers = {r["center_residue_index"] for r in
               read_tsv(ctx.full_results / "06_FINAL_HOTSPOTS" /
                        "significant_hotspot_centers.tsv")}
    assert centers == set(plddt["primary_centers"])


@pytest.mark.integration
@pytest.mark.slow
def test_review_stars_are_not_evaluable_without_a_declared_channel(tmp_path):
    ctx, handoff, _ = full_synthetic_run(tmp_path, B=10_000, stars=False)
    out = run_stage_b(ctx, upstream=handoff)
    review = json.loads((ctx.full_results / "11_SENSITIVITY" /
                         "sensitivity_review_status.json").read_text())
    assert review["channel"]["available"] is False
    for stratum in review["strata"].values():
        assert stratum["status"] == "NOT_EVALUABLE"
        assert "sensitivity_channel" in stratum["reason"]
    assert out.payload["sensitivity_summary"]["star1"]["status"] == "NOT_EVALUABLE"


# --- determinism -------------------------------------------------------------

@pytest.mark.integration
@pytest.mark.slow
def test_rerun_at_the_same_seeds_is_byte_identical(determinism_pair):
    """The determinism assertion: identical inputs and derived seeds, identical bytes."""
    (a_ctx, a_handoff), (b_ctx, b_handoff) = determinism_pair

    assert a_handoff.payload["hotspot_radius"] == b_handoff.payload["hotspot_radius"]
    assert a_handoff.payload["derived_seeds"] == b_handoff.payload["derived_seeds"]

    compared = 0
    for path in sorted((a_ctx.full_results).rglob("*")):
        if not path.is_file() or _is_timestamped(path):
            continue
        rel = path.relative_to(a_ctx.full_results)
        if rel.parts[0] == "12_REPRODUCIBILITY":       # provenance carries wall times
            continue
        other = b_ctx.full_results / rel
        assert other.is_file(), rel
        assert sha256_file(path) == sha256_file(other), f"non-deterministic output: {rel}"
        compared += 1
    assert compared > 40                                # the whole result set, not a token

    # The warning files carry a detection timestamp; their content must still match.
    assert (_warnings_without_timestamps(a_ctx.full_results)
            == _warnings_without_timestamps(b_ctx.full_results))
    # ... and so must the whole handoff payload, timestamps aside.
    assert _payload_without_timestamps(a_handoff) == _payload_without_timestamps(b_handoff)


@pytest.mark.integration
@pytest.mark.slow
def test_figures_are_deterministic(determinism_pair):
    (a_ctx, _), (b_ctx, _) = determinism_pair
    figures = sorted(p for p in a_ctx.full_results.rglob("figures/*") if p.is_file())
    assert len(figures) == 12                          # F1-F6, PNG + SVG
    for path in figures:
        rel = path.relative_to(a_ctx.full_results)
        assert sha256_file(path) == sha256_file(b_ctx.full_results / rel), rel


# --- input validation and barriers -------------------------------------------

@pytest.mark.unit
def test_upstream_qc_fail_blocks(tmp_path):
    ctx, handoff, _ = full_synthetic_run(tmp_path, B=200, qc_status="FAIL")
    with pytest.raises(BlockedError, match="qc_status=FAIL"):
        run_stage_b(ctx, upstream=handoff)


@pytest.mark.unit
def test_manifest_hash_mismatch_blocks(tmp_path):
    ctx, handoff, geo = full_synthetic_run(tmp_path, B=200)
    target = Path(geo["paths"]["classified_cohort"])
    target.write_text(target.read_text() + "\n")           # tamper after hashing
    with pytest.raises(BlockedError, match="manifest verification failed"):
        run_stage_b(ctx, upstream=handoff)


@pytest.mark.unit
def test_forbidden_downstream_column_aborts_as_a_leakage_event(tmp_path):
    ctx, handoff, geo = full_synthetic_run(tmp_path, B=200)
    path = Path(geo["paths"]["classified_cohort"])
    rows = read_tsv(path)
    for row in rows:
        row["review_status"] = "criteria_provided"          # a forbidden column
    write_tsv(path, rows, ["residue_index", "class", "review_status"])
    handoff.manifest[str(path.relative_to(ctx.run_root))] = sha256_file(path)
    with pytest.raises(LeakageError, match="forbidden column"):
        run_stage_b(ctx, upstream=handoff)


@pytest.mark.unit
def test_non_default_b_is_blocked_outside_synthetic_mode(tmp_path):
    ctx, handoff, _ = full_synthetic_run(tmp_path, B=200)
    object.__setattr__(ctx, "synthetic", False)
    with pytest.raises(BlockedError, match="not the frozen default"):
        run_stage_b(ctx, upstream=handoff)


@pytest.mark.unit
def test_a_non_empty_owned_directory_blocks_the_rerun(tmp_path):
    ctx, handoff, _ = full_synthetic_run(tmp_path, B=200)
    stale = ctx.full_results / "06_FINAL_HOTSPOTS"
    stale.mkdir(parents=True, exist_ok=True)
    (stale / "significant_hotspot_centers.tsv").write_text("stale\n")
    with pytest.raises(BlockedError, match="is not empty"):
        run_stage_b(ctx, upstream=handoff)


@pytest.mark.unit
def test_labelled_residue_outside_the_universe_blocks(tmp_path):
    ctx, handoff, geo = full_synthetic_run(tmp_path, B=200)
    path = Path(geo["paths"]["classified_cohort"])
    rows = read_tsv(path)
    rows.append({"residue_index": 999_999, "class": "PLP"})
    write_tsv(path, rows, ["residue_index", "class"])
    handoff.manifest[str(path.relative_to(ctx.run_root))] = sha256_file(path)
    with pytest.raises(BlockedError, match="not in U_struct"):
        run_stage_b(ctx, upstream=handoff)


@pytest.mark.unit
def test_an_empty_class_blocks(tmp_path):
    ctx, handoff, geo = full_synthetic_run(tmp_path, B=200)
    path = Path(geo["paths"]["classified_cohort"])
    rows = [r for r in read_tsv(path) if r["class"] == "PLP"]
    write_tsv(path, rows, ["residue_index", "class"])
    handoff.manifest[str(path.relative_to(ctx.run_root))] = sha256_file(path)
    with pytest.raises(BlockedError, match="a class is empty"):
        run_stage_b(ctx, upstream=handoff)


@pytest.mark.unit
def test_downstream_directories_are_never_read(shared_underpowered):
    # The fixture plants decoy 07_*/08_*/09_*/10_* files before the run.
    ctx, _, _ = shared_underpowered
    barrier = json.loads((ctx.full_results / "12_REPRODUCIBILITY" / "provenance" /
                          "provenance_06_final_hotspots.json").read_text()
                         )["information_barrier"]
    assert barrier["downstream_reads"] == 0
    assert not any("07_FOOTPRINT_RADIUS" in p or "08_FINAL_FOOTPRINT" in p or
                   "09_ROBUSTNESS" in p or "10_ANNOTATION" in p
                   for p in barrier["paths_read"])


# --- v2 §5.4: UNDERPOWERED is not a negative ---------------------------------

@pytest.mark.integration
def test_underpowered_run_is_never_reported_as_a_negative(shared_underpowered):
    """v2 §0/§5.4/Appendix A — the design could not have rejected, so nothing is claimed.

    Under v1 this exact configuration returned ``NO_ADMISSIBLE_RADIUS`` with
    ``qc_status = PASS`` and a ``negative_result`` block asserting the effect was
    absent. That is the misreport the certificate exists to prevent.
    """
    ctx, out, _ = shared_underpowered

    assert out.payload["terminal_state"] == "UNDERPOWERED"
    assert out.payload["outcome_type"] == "UNINFORMATIVE"
    assert out.negative_result is None                   # NOT a negative result
    assert out.qc_status == "FAIL"                       # every consumer blocks on this
    assert out.payload["power_certificate_passed"] is False
    assert out.payload["n_significant_centers"] == 0

    cert = out.payload["power_certificate"]
    assert cert["p_floor"] > cert["c_1_rank1_critical_value"]
    assert cert["TEST_CANNOT_REJECT"] is True
    assert cert["binding_floor"] == "permutation_resolution"
    # p_res binds, so the report must name the B that would clear c_1 (v2 §5.4).
    assert out.payload["increasing_B_is_futile"] is False
    assert 1 / (out.payload["B_to_clear_c_1"] + 1) <= cert["c_1_rank1_critical_value"]
    assert any("B = " in r for r in cert["remedies"])

    status = json.loads((ctx.full_results / "05_HOTSPOT_RADIUS" /
                         "stage_status.json").read_text())
    assert status["status"] == "UNDERPOWERED"
    assert status["outcome_type"] == "UNINFORMATIVE"
    assert "negative_result" not in status

    prose = (ctx.full_results / "05_HOTSPOT_RADIUS" / "NOT_RUN.txt").read_text()
    assert "UNINFORMATIVE" in prose
    assert "NOT a scientific negative" in prose

    warnings = read_tsv(ctx.full_results / "05_HOTSPOT_RADIUS" /
                        "warnings_05_HOTSPOT_RADIUS.tsv")
    blocking = [w for w in warnings if w["warning_code"] == "TEST_CANNOT_REJECT"]
    assert blocking and blocking[0]["severity"] == "BLOCKING"

    # The handoff is still a well-formed handoff_02: the orchestrator must be able to
    # read the terminal state rather than crash on a truncated payload.
    assert_handoff_shape(out)
    assert is_underpowered(out)
    for key in HANDOFF_REQUIRED_KEYS["handoff_02"]:
        assert key in out.payload, key
    # v2 forbids NO_SIGNIFICANT_HOTSPOTS on an underpowered run.
    for stage in STAGES:
        rows = read_tsv(ctx.full_results / stage / f"warnings_{stage}.tsv")
        assert not any(w["warning_code"] == "NO_SIGNIFICANT_HOTSPOTS" for w in rows)


@pytest.mark.integration
def test_underpowered_report_states_the_prohibited_inferences(shared_underpowered):
    """v2 Appendix A — the three forbidden sentences never appear; the fix does."""
    ctx, _, _ = shared_underpowered
    report = (ctx.full_results / "06_FINAL_HOTSPOTS" / "stage_b_report.md").read_text()

    assert "UNDERPOWERED" in report
    assert ("the test could not have rejected any center regardless of the data"
            in report)
    assert "uninformative about the presence or absence of hotspots" in report
    for forbidden in ("No 3D hotspot is detectable in this gene.\n",
                      "This is a valid, complete scientific negative.\n"):
        # They appear only inside the explicit "prohibited inferences" list, which is
        # a quotation of what may NOT be said, never an assertion.
        assert forbidden not in report


@pytest.mark.integration
def test_power_certificate_is_emitted_in_the_underpowered_terminal_state(
        shared_underpowered):
    """v2 §11 — a first-class artefact in EVERY terminal state, not only on success."""
    ctx, _, _ = shared_underpowered
    cert = json.loads((ctx.full_results / "06_FINAL_HOTSPOTS" /
                       "power_certificate.json").read_text())
    assert cert["PASSES"] is False
    assert cert["TEST_CANNOT_REJECT"] is True
    assert cert["preflight"]["null_model"] == "structure_aware_positional"
    assert cert["preflight"]["per_radius"], "the per-radius evidence must survive"
    assert all(not row["passes"] for row in cert["preflight"]["per_radius"]), \
        "pre-flight termination requires that NO radius in the domain could reject"
    assert cert["B_was_changed_by_the_agent"] is False
    assert cert["fdr_method_was_changed_by_the_agent"] is False


# --- negative results --------------------------------------------------------

@pytest.mark.integration
@pytest.mark.slow
def test_certified_negative_carries_its_power_certificate(shared_certified_negative):
    """v2 §5.4 — COMPLETED_NEGATIVE is legitimate only with a PASSED certificate."""
    ctx, out, _ = shared_certified_negative

    assert out.payload["terminal_state"] == "COMPLETED_NEGATIVE"
    assert out.qc_status != "FAIL"                       # a negative is not a failure
    assert out.payload["n_significant_centers"] == 0
    assert out.negative_result["condition"] == "NO_SIGNIFICANT_HOTSPOT_CENTERS"
    assert out.negative_result["power_certificate_passed"] is True
    assert out.negative_result["licenses_no_method_modification"] is True

    cert = out.negative_result["power_certificate"]
    assert cert["PASSES"] is True
    assert cert["p_floor"] <= cert["c_1_rank1_critical_value"]

    # The scan trace exists in full and every radius stayed in the optimization.
    rows = read_tsv(ctx.full_results / "05_HOTSPOT_RADIUS" / "r_hot_scan.tsv")
    assert rows and any(r["n_significant_centers"] == 0 for r in rows)
    assert all(r["admissible"] for r in rows if r["n_significant_centers"] == 0), \
        "v2 §4 — zero significant centers must never remove a radius"

    status = json.loads((ctx.full_results / "06_FINAL_HOTSPOTS" /
                         "stage_status.json").read_text())
    assert status["status"] == "COMPLETED_NEGATIVE"
    assert status["outcome_type"] == "SCIENTIFIC_NEGATIVE"


@pytest.mark.integration
def test_negative_run_still_writes_every_stage_status_and_manifest(shared_underpowered):
    ctx, _, _ = shared_underpowered
    for stage in STAGES:
        assert (ctx.full_results / stage / "stage_status.json").is_file()
        assert (ctx.full_results / stage / f"warnings_{stage}.tsv").is_file()
        manifest = read_tsv(ctx.full_results / stage / "stage_manifest.tsv")
        assert manifest
        for row in manifest:                             # absence is explicit (P3)
            assert row["status"] in ("CREATED", "NOT_CREATED")
            if row["status"] == "NOT_CREATED":
                assert row["reason"] is not None


@pytest.mark.integration
@pytest.mark.slow
def test_global_clustering_is_reported_but_never_terminates_discovery(
        shared_certified_negative):
    """F5: a non-significant global result sets a caveat flag and nothing more."""
    ctx, out, _ = shared_certified_negative
    summary = json.loads((ctx.full_results / "04_GLOBAL_CLUSTERING" /
                          "ripleys_k_summary.json").read_text())
    assert summary["non_significant_terminates_discovery"] is False
    assert set(summary["global_test"]) == {"PLP", "BLB"}
    assert out.payload["global_clustering_flag"] in ("significant", "non_significant")
    # The radius scan ran regardless of the global verdict.
    assert (ctx.full_results / "05_HOTSPOT_RADIUS" / "r_hot_scan.tsv").is_file()


@pytest.mark.integration
def test_pcf_peak_is_recorded_but_never_selected(shared_underpowered):
    ctx, _, _ = shared_underpowered
    peaks = json.loads((ctx.full_results / "04_GLOBAL_CLUSTERING" /
                        "pair_correlation_peaks.json").read_text())
    assert peaks["USED_AS_R_HOT"] is False
    domain = json.loads((ctx.full_results / "04_GLOBAL_CLUSTERING" /
                         "candidate_radius_domain.json").read_text())
    assert domain["pcf_peak_used_as_r_hot"] is False
    assert domain["automatic_domain_expansion"] is False


# --- Lead rulings: config-sourced constants and mandated reporting ------------

@pytest.mark.unit
def test_k_target_and_alpha_are_read_from_config(shared_underpowered):
    """Both constants are auditable config keys, not implicit agent choices."""
    ctx, _, _ = shared_underpowered
    assert ctx.config.get("permutation_resolution_diagnostic.k_target") == 1
    assert ctx.config.get("global_clustering.alpha") == 0.05

    summary = json.loads((ctx.full_results / "04_GLOBAL_CLUSTERING" /
                          "ripleys_k_summary.json").read_text())
    assert summary["alpha"] == 0.05
    assert "config:global_clustering.alpha" in summary["alpha_source"]

    path = ctx.full_results / "06_FINAL_HOTSPOTS" / "permutation_resolution_diagnostic.json"
    if path.is_file():                       # absent only if the run stopped earlier
        diagnostic = json.loads(path.read_text())
        for phase in ("preflight", "posthoc"):
            block = diagnostic[phase]
            if block.get("computed") is False:   # the run stopped before this phase
                assert block["reason"]
                continue
            assert block["k_target"] == 1
            assert block["k_target_source"] == \
                "config:permutation_resolution_diagnostic.k_target"


@pytest.mark.unit
def test_alpha_is_guarded_as_a_frozen_constant(tmp_path):
    ctx, handoff, _ = full_synthetic_run(
        tmp_path, B=200, overrides={"global_clustering": {"alpha": 0.10}})
    with pytest.raises(BlockedError, match="FROZEN methodology violation"):
        run_stage_b(ctx, upstream=handoff)


@pytest.mark.unit
def test_alpha_must_match_the_envelope_percentiles(tmp_path):
    """The plotted envelope and the global-test decision level must be one level.

    Moving the percentiles while leaving alpha at 0.05 would draw a 90% envelope next
    to a 5% decision — the second, independent guard catches exactly that.
    """
    ctx, handoff, _ = full_synthetic_run(
        tmp_path, B=200,
        overrides={"global_clustering": {"envelope_percentiles": [5.0, 95.0]}})
    with pytest.raises(BlockedError, match="envelope percentiles"):
        run_stage_b(ctx, upstream=handoff)


@pytest.mark.integration
@pytest.mark.slow
def test_negative_zg_is_explained_rather_than_left_bare(shared_positive):
    """A negative Zg must never be readable as evidence AGAINST a hotspot (II.5B)."""
    ctx, handoff, _ = shared_positive
    zg = handoff.payload["perm_evidence_zg_at_r_hot"]
    note = handoff.payload["perm_evidence_zg_direction"]
    assert "NOT evidence against a hotspot" in note
    assert "never a significance test" in note

    objectives = json.loads((ctx.full_results / "05_HOTSPOT_RADIUS" /
                             "objective_matrix.json").read_text())
    assert objectives["perm_evidence_zg_direction"] == note

    report = (ctx.full_results / "06_FINAL_HOTSPOTS" / "stage_b_report.md").read_text()
    assert "Reading `perm_evidence_zg`" in report
    assert "NOT evidence against a hotspot" in report
    if zg < 0:
        assert "is the **expected** reading when P/LP clustering is present" in report


@pytest.mark.integration
@pytest.mark.slow
def test_bh_margin_is_reported_against_the_rank_one_critical_value(shared_positive):
    """A reader can see how close the evidence came without recomputing anything."""
    ctx, handoff, _ = shared_positive
    p = handoff.payload
    m = json.loads((ctx.full_results / "06_FINAL_HOTSPOTS" /
                    "permutation_resolution_diagnostic.json").read_text())["bh_margin"]

    assert m["rank1_critical_value"] == pytest.approx(p["q"] / m["m_family"])
    assert m["p_resolution_floor"] == pytest.approx(1 / (p["B"] + 1))
    assert m["smallest_observed_p"] >= m["p_resolution_floor"] - 1e-15
    assert m["margin_ratio"] == pytest.approx(
        m["smallest_observed_p"] / m["rank1_critical_value"])
    assert p["bh_rank1_critical_value"] == m["rank1_critical_value"]
    assert p["bh_smallest_observed_p"] == m["smallest_observed_p"]
    assert p["bh_margin_ratio"] == m["margin_ratio"]

    # Ranks are ascending in p and each carries its own critical value.
    ranks = m["top_ranks"]
    assert [r["rank"] for r in ranks] == list(range(1, len(ranks) + 1))
    assert [r["p_emp"] for r in ranks] == sorted(r["p_emp"] for r in ranks)
    for row in ranks:
        assert row["bh_critical_value"] == pytest.approx(
            row["rank"] * p["q"] / m["m_family"])

    report = (ctx.full_results / "06_FINAL_HOTSPOTS" / "stage_b_report.md").read_text()
    assert "rank-1 BH critical value q/m" in report
    assert "margin ratio" in report


@pytest.mark.unit
def test_production_default_permutation_count():
    """The FROZEN production configuration, asserted independently of any fixture.

    Some tests in this module deliberately run at a reduced B through an explicit
    synthetic config overlay, because they exercise the resolution-limited branch.
    None of that can reach production: this test reads `config/pipeline.yaml`
    directly — no overlay, no fixture — and pins the shipped values. Stage B
    additionally refuses a non-default B unless the run is flagged synthetic, which
    is covered by test_non_default_b_is_blocked_outside_synthetic_mode.
    """
    cfg = load_config(Path(__file__).resolve().parents[2] / "config" / "pipeline.yaml")
    # DECISION-B-DEFAULT-0001 raised this from 10_000; B is a resolution parameter,
    # not a threshold, so no scientific criterion moved with it.
    assert cfg.get("permutation.B_default") == 100_000
    assert cfg.get("fdr.q") == 0.05
    assert cfg.get("fdr.method") == "BH"
    assert cfg.get("loo_mcc.kappa") == 2.0
    assert cfg.get("radius_selection.w_k") == 1.0
    assert cfg.get("permutation_resolution_diagnostic.auto_change_B") is False
    assert_frozen_methodology(cfg)          # the whole frozen block, unmodified
