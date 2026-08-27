"""Phase D end to end: baseline reproduction, both metric families, honesty rules."""
from __future__ import annotations

import json
import tarfile

import pytest

from synthetic_footprint import synthetic_run, two_cluster_centers

from hotspot3d.footprint.stage import run_phase_c
from hotspot3d.robustness import perturb
from hotspot3d.robustness.classify import CLASSIFICATION_FIELDS
from hotspot3d.robustness.compare import GEOMETRY_FIELDS
from hotspot3d.robustness.stage import (EXPECTED_09, assert_no_hotspot_statistics_import,
                                        run_phase_d)
from hotspot3d.utils.errors import BlockedError
from hotspot3d.utils.io import read_tsv

pytestmark = pytest.mark.unit


def _stage(ctx):
    return ctx.full_results / "09_ROBUSTNESS"


# --- preconditions -----------------------------------------------------------

def test_baseline_reproduction_is_byte_identical(phase_d_run):
    ctx, _, _, _ = phase_d_run
    payload = json.loads((_stage(ctx) / "baseline_reproduction.json").read_text())

    assert payload["performed"] is True
    assert payload["byte_identical"] is True
    assert payload["mismatched_fields"] == []
    assert payload["original"]["covered_residues_sha256"] == \
        payload["reproduced"]["covered_residues_sha256"]
    assert payload["original"]["occupancy_sha256"] == \
        payload["reproduced"]["occupancy_sha256"]


def test_code_version_drift_blocks_phase_d(tmp_path):
    ctx, upstream = synthetic_run(tmp_path, center_ids=two_cluster_centers()[:4])
    handoff_03, phase_c = run_phase_c(ctx, upstream=upstream)

    freeze_path = ctx.full_results / "08_FINAL_FOOTPRINT" / "footprint_freeze.json"
    payload = json.loads(freeze_path.read_text())
    payload["footprint_code_version"] = "0" * 64
    freeze_path.write_text(json.dumps(payload))

    with pytest.raises(BlockedError) as excinfo:
        run_phase_d(ctx, phase_c=phase_c, handoff_03=handoff_03)
    assert "code_version drift" in str(excinfo.value)


def test_frozen_artifacts_are_re_verified_at_phase_d_entry(tmp_path):
    ctx, upstream = synthetic_run(tmp_path, center_ids=two_cluster_centers()[:4])
    handoff_03, phase_c = run_phase_c(ctx, upstream=upstream)

    target = ctx.full_results / "08_FINAL_FOOTPRINT" / "footprint_residues.tsv"
    target.write_text(target.read_text() + "# tampered\n")

    with pytest.raises(BlockedError) as excinfo:
        run_phase_d(ctx, phase_c=phase_c, handoff_03=handoff_03)
    assert "BASELINE INTEGRITY FAILURE" in str(excinfo.value)


def test_the_robustness_procedure_does_not_import_hotspot_statistics():
    assert_no_hotspot_statistics_import()


# --- outputs -----------------------------------------------------------------

def test_every_expected_output_exists(phase_d_run):
    ctx, _, _, _ = phase_d_run
    for name, _description in EXPECTED_09:
        assert (_stage(ctx) / name).is_file(), name


def test_design_is_recorded_and_accounts_for_every_iteration(phase_d_run):
    ctx, _, _, handoff_04 = phase_d_run
    levels = read_tsv(_stage(ctx) / "perturbation_design.tsv")
    design = json.loads((_stage(ctx) / "subset_design.json").read_text())
    results = read_tsv(_stage(ctx) / "perturbation_results.tsv")

    assert design["recorded_before_iteration_1"] is True
    assert sum(level["n_evaluated"] for level in levels) == len(results)
    assert design["total_evaluated"] == len(results)
    assert design["total_subset_space_sum_T_k"] == \
        sum(level["total_subsets_T_k"] for level in levels)
    for level in levels:
        assert level["mode"] in ("exhaustive", "sampled")
        assert level["sampling_fraction"] is not None

    iterations = handoff_04.payload["robustness_profile"]["iterations"]
    assert iterations["planned"] == iterations["completed"] + iterations["failed"]
    assert iterations["planned"] == len(results)


def test_results_carry_all_metric_blocks(phase_d_run):
    """DECISION-STAGE-D-FIXED-RFP-0001: one fixed-r_fp geometry block per
    iteration, plus both classification families (footprint and r_hot coverage)."""
    ctx, _, _, _ = phase_d_run
    rows = read_tsv(_stage(ctx) / "perturbation_results.tsv")
    assert rows

    for row in rows:
        for field in GEOMETRY_FIELDS:                      # the sole reconstruction
            assert field in row
        for field in CLASSIFICATION_FIELDS:                # footprint family
            assert field in row
        for field in CLASSIFICATION_FIELDS:                # r_hot coverage family
            assert f"{field}_hs" in row
        assert row["status"] in ("OK", "FAILED_GEOMETRY_FAULT")
        assert "robustness_r_fp" in row


def test_both_classification_families_are_distinctly_labelled(phase_d_run):
    """F10 permits r_fp < r_hot; a footprint-only MCC would then deflate by design."""
    ctx, _, _, handoff_04 = phase_d_run
    profile = handoff_04.payload["robustness_profile"]

    assert profile["mcc_footprint_family"]["n"] > 0
    assert profile["mcc_hotspot_coverage_family"]["n"] > 0

    summary = read_tsv(_stage(ctx) / "robustness_summary.tsv")
    families = {row["reconstruction"] for row in summary}
    assert {"fixed_radius", "classification_footprint",
            "classification_hotspot_coverage", "diagnostic"} <= families


def test_every_reported_metric_is_aggregated_with_the_full_summary(phase_d_run):
    ctx, _, _, _ = phase_d_run
    rows = read_tsv(_stage(ctx) / "robustness_summary.tsv")
    metrics = {(r["metric"], r["reconstruction"]) for r in rows}

    for field in GEOMETRY_FIELDS:
        assert (field, "fixed_radius") in metrics
    for field in CLASSIFICATION_FIELDS:
        assert (field, "classification_footprint") in metrics
        assert (field, "classification_hotspot_coverage") in metrics
    # unflattering metrics are present alongside the flattering ones
    assert ("jaccard_residues", "fixed_radius") in metrics
    assert ("mcc", "classification_footprint") in metrics

    for row in rows:
        assert row["n_total"] is not None
        assert row["ci_method"] is not None


def test_center_recurrence_is_the_reappearance_rate_when_removed(phase_d_run):
    """¶51 recurrence: does the footprint reappear over a center that is GONE?"""
    ctx, phase_c, _, _ = phase_d_run
    rows = read_tsv(_stage(ctx) / "center_sensitivity.tsv")

    for row in rows:
        assert row["n_iterations_removing"] + row["n_iterations_retaining"] == \
            len(read_tsv(_stage(ctx) / "perturbation_results.tsv"))
        rate, hits, n = (row["center_recurrence_rate"],
                         row["n_recurrent_iterations"],
                         row["n_iterations_removing"])
        assert 0.0 <= rate <= 1.0
        assert rate == pytest.approx(hits / n)

    # recurrence must actually discriminate; a constant would be no evidence
    assert len({row["center_recurrence_rate"] for row in rows}) > 1 or \
        all(r["center_recurrence_rate"] < 1.0 for r in rows)


def test_the_vacuous_retained_center_form_is_not_emitted(phase_d_run):
    """A retained center is trivially covered; such a column would be constant.

    Constant-by-construction is the worst kind of reported number: it reads as a
    strong result while carrying no information. The rejection is recorded in
    RECURRENCE_DEFINITION so it cannot be quietly reinstated.
    """
    from hotspot3d.robustness.sensitivity import (CENTER_SENSITIVITY_COLUMNS,
                                                  RECURRENCE_DEFINITION)

    ctx, _, _, handoff_04 = phase_d_run
    header = read_tsv(_stage(ctx) / "center_sensitivity.tsv")[0]

    assert "retained_containment_rate" not in CENTER_SENSITIVITY_COLUMNS
    assert "retained_containment_rate" not in header
    assert not any("retained" in c and "containment" in c for c in header)

    assert "vacuous and must not be substituted" in RECURRENCE_DEFINITION
    profile = handoff_04.payload["robustness_profile"]
    assert profile["center_recurrence_definition"] == RECURRENCE_DEFINITION


def test_recurrence_denominator_zero_gives_null_not_zero():
    """"Never measured" must not be presentable as "never recurred"."""
    from hotspot3d.robustness.sensitivity import center_influence

    class _R:
        removed, residues = (), frozenset()
        row, failed = {"jaccard_residues": 1.0}, False

    rows = center_influence([_R()], (11, 22), {}, {})
    for row in rows:
        assert row["n_iterations_removing"] == 0
        assert row["center_recurrence_rate"] is None
        assert row["mean_jaccard_drop"] is None


def test_failed_iterations_count_as_non_recurrence():
    from hotspot3d.robustness.sensitivity import center_influence

    class _R:
        def __init__(self, removed, residues, failed):
            self.removed, self.residues = removed, residues
            self.failed = failed
            self.row = {"jaccard_residues": 0.0 if failed else 0.8}

    results = [_R((7,), frozenset({7}), False),      # recurred
               _R((7,), frozenset(), True)]          # failed -> non-recurrence
    row = center_influence(results, (7,), {}, {})[0]
    assert row["n_iterations_removing"] == 2
    assert row["n_recurrent_iterations"] == 1
    assert row["center_recurrence_rate"] == pytest.approx(0.5)
    assert row["n_failed_iterations_removing"] == 1


def test_center_sensitivity_reports_contributing_counts(phase_d_run):
    ctx, phase_c, _, _ = phase_d_run
    rows = read_tsv(_stage(ctx) / "center_sensitivity.tsv")
    assert len(rows) == len(phase_c.inputs.centers)

    ranks = [r["rank"] for r in rows]
    assert ranks == sorted(ranks)
    drops = [r["mean_jaccard_drop"] for r in rows]
    assert drops == sorted(drops, reverse=True)
    for row in rows:
        assert row["n_iterations_removing"] > 0
        assert row["carries_clinvar_variant"] in (True, False)


def test_footprints_are_bundled_not_scattered_across_directories(phase_d_run):
    """Inode control: one archive replaces per-iteration directories (IX.3)."""
    ctx, _, _, _ = phase_d_run
    stage = _stage(ctx)
    subdirs = {p.name for p in stage.iterdir() if p.is_dir()}
    assert subdirs == {"figures"}

    with tarfile.open(stage / "perturbation_footprints.tar.gz") as tar:
        names = tar.getnames()
    results = read_tsv(stage / "perturbation_results.tsv")
    assert len(names) == len(results)
    assert names[0] == "iter_000000.txt"
    assert all(name.endswith(".txt") for name in names)


def test_robustness_r_fp_is_invariant_across_all_iterations(phase_d_run):
    """DECISION-STAGE-D-FIXED-RFP-0001: r_fp is frozen after Stage C; every
    completed iteration reconstructs at exactly that value, never a re-selected
    one. Superseded the old rfp_diagnostics.tsv drift-reporting test, since
    there is no longer any per-iteration radius selection that could drift."""
    ctx, phase_c, _, handoff_04 = phase_d_run
    rows = read_tsv(_stage(ctx) / "perturbation_results.tsv")
    assert rows

    completed = [r for r in rows if r["status"] == "OK"]
    assert completed
    for row in completed:
        assert row["robustness_r_fp"] == pytest.approx(phase_c.r_fp)

    invariant = handoff_04.payload["robustness_r_fp_invariant"]
    assert invariant["invariant_holds"] is True
    assert invariant["n_equal_to_original"] == invariant["n_iterations_completed"]
    assert invariant["r_fp_original"] == phase_c.r_fp


# --- honesty rules -----------------------------------------------------------

def test_no_categorical_robustness_verdict_is_emitted(phase_d_run):
    ctx, _, _, handoff_04 = phase_d_run
    banned = ("NOT_ROBUST", "\"ROBUST\"", "'ROBUST'", "ROBUST/NOT_ROBUST")

    for path in sorted(_stage(ctx).rglob("*")):
        if not path.is_file() or path.suffix in (".png", ".svg", ".gz"):
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if path.name in ("stage_d_report.md", "robustness_profile.json"):
            assert "NOT_ROBUST" not in text.replace(
                "ROBUST / NOT_ROBUST verdict", "")
            continue
        for token in banned:
            assert token not in text, f"{path.name} contains {token}"

    profile = handoff_04.payload["robustness_profile"]
    assert profile["categorical_verdict_emitted"] is False
    assert "verdict" not in {k.lower() for k in profile if k != "categorical_verdict_emitted"}


def test_profile_exposes_the_keys_the_lead_reads_verbatim(phase_d_run):
    """REVIEW_PACK and gene_summary index these names directly; they are contract."""
    ctx, _, _, handoff_04 = phase_d_run
    on_disk = json.loads((_stage(ctx) / "robustness_profile.json").read_text())

    for profile in (on_disk, handoff_04.payload["robustness_profile"]):
        for key in ("median", "iqr_lo", "iqr_hi"):
            assert profile["jaccard_residues"][key] is not None, key
        assert profile["jaccard_residues"]["iqr_lo"] <= \
            profile["jaccard_residues"]["median"] <= \
            profile["jaccard_residues"]["iqr_hi"]
        assert 0.0 <= profile["preservation_freq_050"] <= 1.0
        assert 0.0 <= profile["preservation_freq_070"] <= 1.0
        # the aliases are the same numbers as the detailed block, not a recount
        assert profile["preservation_freq_050"] == \
            profile["preservation"]["P_jaccard_ge_primary"]
        assert profile["preservation_freq_070"] == \
            profile["preservation"]["P_jaccard_ge_secondary"]


def test_limited_robustness_uses_its_own_warning_code():
    """Instability is surfaced under a code that names it, never a proxy."""
    from hotspot3d.utils.status import WARNING_CODES, Severity

    assert WARNING_CODES["LIMITED_GEOMETRIC_ROBUSTNESS"] == Severity.MAJOR


def test_handoff_04_declares_the_scope_of_the_result(phase_d_run):
    _, _, _, handoff_04 = phase_d_run
    payload = handoff_04.payload

    assert payload["clinical_dataset_unmodified"] is True
    assert payload["geometric_footprint_robustness_only"] is True
    assert "not validation of hotspot discovery" in payload["statement"].lower()
    for key in ("n_S", "K_MAX", "total_subset_space", "total_evaluated",
                "per_level", "robustness_r_fp_invariant", "robustness_profile",
                "footprint_code_version"):
        assert key in payload


def test_the_cohort_is_untouched_by_perturbation(phase_d_run):
    ctx, phase_c, _, _ = phase_d_run
    from hotspot3d.utils.hashing import sha256_file

    cohort = phase_c.inputs.paths["classified_cohort.tsv"]
    assert sha256_file(cohort) == phase_c.inputs.hashes["classified_cohort.tsv"]

    universe = read_tsv(_stage(ctx) / "perturbation_universe.tsv")
    assert len(universe) == len(phase_c.inputs.centers)
    for row in universe:
        assert row["carries_clinvar_variant"] in (True, False)


def test_a_failed_iteration_keeps_its_place_with_jaccard_zero(monkeypatch,
                                                              phase_d_run):
    """Failures are outcomes with status codes, not omissions.

    DECISION-STAGE-D-FIXED-RFP-0001: an iteration's sole reconstruction is
    reconstruct_at, so a forced failure is injected there, not in build_footprint
    (which a Phase D iteration no longer calls at all — see
    test_the_robustness_procedure_never_re_selects_r_fp below)."""
    _, phase_c, _, _ = phase_d_run
    up = phase_c.inputs

    def exploding(*_args, **_kwargs):
        raise BlockedError("forced for the test")

    monkeypatch.setattr(perturb, "reconstruct_at", exploding)
    result = perturb.run_iteration(
        7, (up.centers.ids[0],), "exhaustive", "iter|0007",
        centers=up.centers, universe=up.universe, params=phase_c.params,
        original_state=phase_c.solution.selected,
        original_residues=frozenset(phase_c.solution.covered_residue_ids(up.universe)),
        r_fp_original=phase_c.r_fp, plp=up.cohort_plp, blb=up.cohort_blb,
        r_hot_value=up.r_hot.value_for_coverage_query, undefined_mcc=0.0)

    assert result.failed is True
    assert result.status == "FAILED_GEOMETRY_FAULT"
    assert result.row["jaccard_residues"] == 0.0
    assert result.row["dice_residues"] == 0.0
    assert result.row["robustness_r_fp"] is None
    assert result.row["failure_reason"] != "NA"
    # the independent _hs family still carries information; it never reconstructs
    assert result.row["mcc_hs"] is not None


def test_the_robustness_procedure_never_re_selects_r_fp():
    """DECISION-STAGE-D-FIXED-RFP-0001, the source-level pin required by it."""
    from hotspot3d.robustness.stage import assert_no_radius_optimization_in_perturb

    assert_no_radius_optimization_in_perturb()      # raises BlockedError on violation
    assert perturb.reconstruct_at is not None        # actually importable, not just named


def test_no_stage_d_iteration_invokes_footprint_radius_optimization(monkeypatch,
                                                                     tmp_path):
    """DECISION-STAGE-D-FIXED-RFP-0001, the dynamic pin: a real run, not just source.

    build_footprint (domain re-derivation -> multi-scale sweep -> QC -> Pareto ->
    distance-to-ideal selection) is spied on across a real Phase C + Phase D run.
    It legitimately fires exactly twice — once for Stage C's own FP_original
    selection, once for Phase D's mandatory baseline reproduction on the FULL,
    unperturbed center set (agent §6.10) — and must fire ZERO additional times
    regardless of how many perturbation iterations run. If a future edit ever
    reintroduces per-iteration radius selection, this test fails on the call
    count, not merely on an import-time source scan.
    """
    import hotspot3d.footprint.api as api_mod

    calls = []
    real_build_footprint = api_mod.build_footprint

    def spy(*args, **kwargs):
        calls.append(1)
        return real_build_footprint(*args, **kwargs)

    monkeypatch.setattr(api_mod, "build_footprint", spy)
    monkeypatch.setattr("hotspot3d.footprint.stage.build_footprint", spy)
    monkeypatch.setattr("hotspot3d.robustness.stage.build_footprint", spy)

    ctx, upstream = synthetic_run(tmp_path, center_ids=two_cluster_centers()[:4])
    handoff_03, phase_c = run_phase_c(ctx, upstream=upstream)
    run_phase_d(ctx, phase_c=phase_c, handoff_03=handoff_03)

    design = json.loads(
        (ctx.full_results / "09_ROBUSTNESS" / "subset_design.json").read_text())
    assert design["planned_iterations"] >= 1, "the design must actually run iterations"

    assert len(calls) == 2, (
        f"build_footprint was called {len(calls)} times for a run with "
        f"{design['planned_iterations']} planned iterations; expected exactly 2 "
        f"(Stage C's own selection + Phase D's mandatory baseline reproduction), "
        f"never once per iteration.")


def test_stage_status_and_warnings_are_always_present(phase_d_run):
    ctx, _, _, _ = phase_d_run
    status = json.loads((_stage(ctx) / "stage_status.json").read_text())
    assert status["status"] in ("COMPLETED", "COMPLETED_NEGATIVE")
    assert status["outcome_type"] in ("COMPLETED", "NOT_APPLICABLE")
    assert (_stage(ctx) / "warnings_09_ROBUSTNESS.tsv").is_file()


def test_provenance_records_the_design_and_the_seeds(phase_d_run):
    ctx, _, _, _ = phase_d_run
    path = ctx.full_results / "12_REPRODUCIBILITY" / "provenance" \
        / "provenance_09_robustness.json"
    payload = json.loads(path.read_text())

    assert payload["agent_owner"] == "footprint-robustness"
    assert payload["parameters"]["hotspot_statistics_recomputed"] is False
    assert payload["parameters"]["clinvar_records_modified"] is False
    assert payload["parameters"]["per_level"]
    assert any(command.startswith("iteration ") for command in payload["commands"])


# --- n_S = 1 -----------------------------------------------------------------

def test_single_center_is_reported_as_not_evaluable(tmp_path):
    ctx, upstream = synthetic_run(tmp_path, center_ids=two_cluster_centers()[:1])
    handoff_03, phase_c = run_phase_c(ctx, upstream=upstream)
    handoff_04 = run_phase_d(ctx, phase_c=phase_c, handoff_03=handoff_03)

    profile = handoff_04.payload["robustness_profile"]
    assert profile["evaluable"] is False
    assert "ROBUSTNESS_NOT_EVALUABLE" in profile["reason"]
    assert profile["iterations"] == {"planned": 0, "completed": 0, "failed": 0}

    warnings = read_tsv(_stage(ctx) / "warnings_09_ROBUSTNESS.tsv")
    codes = {w["warning_code"]: w["severity"] for w in warnings}
    assert codes["ROBUSTNESS_NOT_EVALUABLE"] == "MAJOR"

    status = json.loads((_stage(ctx) / "stage_status.json").read_text())
    assert status["negative_result"]["condition"] == "ROBUSTNESS_NOT_EVALUABLE"
    assert (_stage(ctx) / "NOT_RUN.txt").is_file()
    # empty-but-present tables, never missing files
    assert read_tsv(_stage(ctx) / "perturbation_results.tsv") == []

    # the contract keys are present and NULL: a consumer must not be able to read
    # "not evaluable" as a low score
    assert profile["preservation_freq_050"] is None
    assert profile["preservation_freq_070"] is None
    assert profile["jaccard_residues"]["median"] is None
    assert profile["jaccard_residues"]["iqr_lo"] is None


# --- II.11 operational cost preflight (advisory, never inferential) ----------

def test_cost_preflight_escalates_before_iteration_one(tmp_path):
    """``robustness.max_wall_seconds`` warns BEFORE the budget is spent.

    Activated in production so an unusually expensive Stage D is surfaced ahead of
    computation rather than discovered hours into it. The guard is exercised here
    with a deliberately tiny budget; the shipped 40 h value is far above any normal
    run, and that it does NOT fire normally is covered by every other test in this
    module running Phase D to completion under the shipped config.

    The properties that make it advisory-only rather than inferential:
      * it raises EscalationRequired, handing the Lead the decision, and never
        chooses for them;
      * the design is already RECORDED at the full configured ``N_CAP`` when it
        fires, so nothing has been silently reduced, truncated or resampled;
      * no iteration has run, so no perturbation result exists to be biased.
    """
    import yaml

    from hotspot3d.utils.config import load_config
    from hotspot3d.utils.errors import EscalationRequired
    from hotspot3d.utils.runctx import RunContext
    from synthetic_footprint import CONFIG_PATH, two_cluster_centers, write_upstream
    from synthetic_footprint import lattice

    overlay = tmp_path / "tiny_budget.yaml"
    overlay.write_text(yaml.safe_dump({"robustness": {"max_wall_seconds": 0.001}}),
                       encoding="utf-8")
    cfg = load_config(CONFIG_PATH, gene_overlay=overlay)
    assert cfg.get("robustness.N_CAP") == 10000, "the budget guard must not touch N_CAP"

    ctx = RunContext.create(gene="SYNTH", config=cfg, results_root=tmp_path / "results",
                            synthetic=True, run_id="20250101T000000Z_cost0001_cost0002")
    ids, coords = lattice(6, 5.5)
    upstream = write_upstream(ctx, center_ids=two_cluster_centers()[:4],
                              ids=ids, coords=coords)
    handoff_03, phase_c = run_phase_c(ctx, upstream=upstream)

    with pytest.raises(EscalationRequired) as excinfo:
        run_phase_d(ctx, phase_c=phase_c, handoff_03=handoff_03)
    assert "exceeds the configured budget" in str(excinfo.value)

    stage = _stage(ctx)
    design = json.loads((stage / "subset_design.json").read_text())
    assert design["N_CAP"] == 10000, (
        "the preflight altered the design instead of escalating — it must never "
        "silently reduce the perturbation budget")
    assert not (stage / "perturbation_results.tsv").exists(), (
        "iterations ran before the cost escalation; the guard is not preflight")


def test_cost_preflight_is_dormant_when_no_budget_is_configured(tmp_path):
    """A null budget keeps the escalation path dormant — the pre-activation behaviour."""
    import yaml

    from hotspot3d.utils.config import load_config
    from hotspot3d.robustness.params import RobustnessParams
    from synthetic_footprint import CONFIG_PATH

    overlay = tmp_path / "no_budget.yaml"
    overlay.write_text(yaml.safe_dump({"robustness": {"max_wall_seconds": None}}),
                       encoding="utf-8")
    params = RobustnessParams.from_config(load_config(CONFIG_PATH, gene_overlay=overlay))
    assert params.max_wall_seconds is None
    assert params.n_cap == 10000
