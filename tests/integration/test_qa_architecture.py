"""The QA architecture, validated against deliberately injected/synthetic failures.

This module proves the DETERMINISTIC MACHINERY behind the QA architecture upgrade
behaves correctly — the stage-contract validator, the failure taxonomy, the
INFRA_CONFIG/SCIENTIFIC_CONFIG split, the write-access manifest, the bug-brief
protocol, and the bounded repair loop. It does not spawn a live
``pipeline-qa-debugger`` agent for every case: what is under test here is the
STRUCTURE that agent's charter is built on (`hotspot3d.utils.qa_taxonomy`,
`hotspot3d.utils.stage_contract`) and that the Lead uses to accept or reject its
output — proving THAT works is what "prove routing works" means for a
policy/plumbing layer, independent of any one agent's live diagnostic skill
(already demonstrated repeatedly elsewhere in this codebase's history).

No real gene, no network. Injected failures are synthetic and local.
"""
from __future__ import annotations

import json

import pytest

from hotspot3d.utils.errors import BlockedError
from hotspot3d.utils.qa_taxonomy import (
    QA_DIRECT_WRITE_GLOBS,
    QA_SYMBOL_CARVEOUTS,
    BugBrief,
    ConfigSubclass,
    FailureClass,
    RepairTrail,
    assert_brief_acceptable_to_lead,
    classify_config_key,
    owner_of_path,
    qa_may_edit_directly,
    root_cause_signature,
)
from hotspot3d.utils.stage_contract import validate_stage_contract

pytestmark = pytest.mark.integration


def _minimal_brief(**overrides) -> BugBrief:
    base = dict(
        bug_id="QA-TEST-0001",
        failure_class=FailureClass.IMPLEMENTATION_BUG,
        stage="06_FINAL_HOTSPOTS",
        owner="hotspot-statistics",
        observed_behavior="observed",
        expected_behavior="expected",
        minimal_reproduction="repro",
        root_cause="cause",
        evidence_artifact="evidence",
        why_this_is_an_implementation_defect="why",
        scientific_methodology_change_required=False,
        files_or_components_affected=["x.py"],
        targeted_test_required="test_x",
    )
    base.update(overrides)
    return BugBrief(**base)


# --------------------------------------------------------------------------- #
# CASE 1 — schema/handoff defect: deterministic validator catches it,         #
# QA classifies IMPLEMENTATION_BUG                                            #
# --------------------------------------------------------------------------- #

def test_case1_schema_defect_caught_by_validator_and_classified(tmp_path):
    """A stage claims COMPLETED but its PRIMARY outputs and handoff are broken."""
    from hotspot3d.utils.config import load_config
    from hotspot3d.utils.runctx import RunContext

    stage_dir = tmp_path / "FULL_RESULTS" / "06_FINAL_HOTSPOTS"
    stage_dir.mkdir(parents=True)
    (stage_dir / "stage_status.json").write_text(json.dumps(
        {"status": "COMPLETED", "outcome_type": "COMPLETED"}), encoding="utf-8")
    (stage_dir / "handoff_02.json").write_text("", encoding="utf-8")  # truncated

    cfg = load_config("config/pipeline.yaml")
    ctx = RunContext.create(gene="INJECT1", config=cfg, results_root=tmp_path / "results",
                            synthetic=True, run_id="20250101T000000Z_qainjc1_qainjc2")
    ctx.run_root = tmp_path

    result = validate_stage_contract(ctx, "06_FINAL_HOTSPOTS")
    assert not result.passed, "the validator must catch a stage claiming COMPLETED " \
        "with missing PRIMARY outputs and a truncated handoff"
    assert any(v.rule == "EXPECTED_OUTPUT_MISSING" for v in result.violations)
    assert any(v.rule == "EXPECTED_OUTPUT_MALFORMED" for v in result.violations)

    # QA classifies from the validator's own evidence — never a narrative retelling.
    evidence = result.summary()
    brief = _minimal_brief(
        failure_class=FailureClass.IMPLEMENTATION_BUG,
        stage="06_FINAL_HOTSPOTS", owner="hotspot-statistics",
        evidence_artifact=evidence,
        why_this_is_an_implementation_defect=(
            "stage_status.json declares COMPLETED but the declared PRIMARY "
            "outputs do not exist and the handoff is a 0-byte file — the writer "
            "did not actually produce what it claimed to"),
    )
    assert brief.validate() == []
    assert "EXPECTED_OUTPUT_MISSING" in brief.evidence_artifact


def test_case1_a_well_formed_stage_passes_clean(tmp_path):
    """The converse: a genuinely well-formed NOT_RUN stage must not be flagged."""
    from hotspot3d.utils.config import load_config
    from hotspot3d.utils.runctx import RunContext

    stage_dir = tmp_path / "FULL_RESULTS" / "09_ROBUSTNESS"
    stage_dir.mkdir(parents=True)
    (stage_dir / "stage_status.json").write_text(json.dumps(
        {"status": "NOT_RUN", "outcome_type": "NOT_APPLICABLE"}), encoding="utf-8")
    (stage_dir / "NOT_RUN.txt").write_text("upstream terminated the chain",
                                           encoding="utf-8")

    cfg = load_config("config/pipeline.yaml")
    ctx = RunContext.create(gene="INJECT1b", config=cfg, results_root=tmp_path / "results",
                            synthetic=True, run_id="20250101T000000Z_qainjc3_qainjc4")
    ctx.run_root = tmp_path

    result = validate_stage_contract(ctx, "09_ROBUSTNESS")
    assert result.passed, result.summary()


# --------------------------------------------------------------------------- #
# CASE 2 — pure glue-code defect: the QA-allowed repair path is real          #
# --------------------------------------------------------------------------- #

def test_case2_glue_files_are_directly_editable_by_qa():
    for path in ("src/hotspot3d/orchestration/pipeline.py",
                 "src/hotspot3d/orchestration/cli.py",
                 "src/hotspot3d/reporting/manifest.py",
                 "src/hotspot3d/reporting/archives.py",
                 "src/hotspot3d/utils/io.py",
                 "src/hotspot3d/utils/hashing.py",
                 "src/hotspot3d/utils/status.py"):
        assert qa_may_edit_directly(path), f"{path} should be QA-direct-editable"


def test_case2_a_glue_defect_brief_validates_as_qa_repairable():
    brief = _minimal_brief(
        bug_id="QA-TEST-0002", failure_class=FailureClass.IMPLEMENTATION_BUG,
        stage="00_RUN_SUMMARY", owner="lead",
        observed_behavior="manifest.tsv omits a row for a file that exists on disk",
        expected_behavior="every CREATED file has a manifest row",
        root_cause="build_manifest() iterated EXPECTED but skipped an ad-hoc file "
                  "written outside the declared list",
        evidence_artifact="manifest.tsv has 40 rows; ls FULL_RESULTS/**/*.tsv "
                          "counts 41 files",
        why_this_is_an_implementation_defect="pure reporting/aggregation plumbing, "
                                             "no scientific value involved",
        files_or_components_affected=["src/hotspot3d/reporting/manifest.py"],
        fix_applied_by_qa=True,
    )
    assert qa_may_edit_directly("src/hotspot3d/reporting/manifest.py")
    assert brief.validate() == []
    assert_brief_acceptable_to_lead(brief)  # does not raise


# --------------------------------------------------------------------------- #
# CASE 3 — scientific-module defect: QA does not edit directly, routes to     #
# the owner, and a brief falsely claiming a direct QA fix is rejected         #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("path,expected_owner", [
    ("src/hotspot3d/hotspot/scan.py", "hotspot-statistics"),
    ("src/hotspot3d/hotspot/selection.py", "hotspot-statistics"),
    ("src/hotspot3d/footprint/selection.py", "footprint-robustness"),
    ("src/hotspot3d/robustness/perturb.py", "footprint-robustness"),
    ("src/hotspot3d/data/classify.py", "data-structure"),
    ("src/hotspot3d/annotation/rubric.py", "biological-annotation"),
])
def test_case3_scientific_files_route_to_their_owner_not_qa(path, expected_owner):
    assert owner_of_path(path) == expected_owner
    assert not qa_may_edit_directly(path), \
        f"{path} must NOT be QA-direct-editable — it is {expected_owner}'s territory"


def test_case3_a_scientific_defect_routes_and_a_false_qa_fix_claim_is_rejected():
    # Correctly routed: QA diagnoses, does not claim to have fixed it.
    routed = _minimal_brief(
        bug_id="QA-TEST-0003", failure_class=FailureClass.IMPLEMENTATION_BUG,
        stage="05_HOTSPOT_RADIUS", owner="hotspot-statistics",
        evidence_artifact="r_hot_scan.tsv column 'fold_enrichment' is negative for "
                          "3 radii, which the formula FE=[N_P_in/N_L_in]/pi cannot "
                          "produce (a ratio of non-negative counts over a positive "
                          "prevalence)",
        files_or_components_affected=["src/hotspot3d/hotspot/scan.py"],
        targeted_test_required="test_fold_enrichment_is_never_negative",
        fix_applied_by_qa=False,
    )
    assert routed.validate() == []
    assert_brief_acceptable_to_lead(routed)

    # Incorrectly claimed: QA says it fixed a file it is not permitted to touch.
    false_claim = _minimal_brief(
        bug_id="QA-TEST-0003b", failure_class=FailureClass.IMPLEMENTATION_BUG,
        stage="05_HOTSPOT_RADIUS", owner="hotspot-statistics",
        evidence_artifact="same as above",
        files_or_components_affected=["src/hotspot3d/hotspot/scan.py"],
        fix_applied_by_qa=True,  # <- the violation under test
    )
    problems = false_claim.validate()
    assert problems, "a QA-applied fix to a scientific-module file must be rejected"
    assert any("write manifest" in p or "not 'lead'" in p for p in problems)
    with pytest.raises(BlockedError):
        assert_brief_acceptable_to_lead(false_claim)

    # The "fix returns to QA, targeted regression PASS" half of CASE 3: once the
    # OWNER applies the fix (not QA), a follow-up brief closing the loop is valid.
    closed = _minimal_brief(
        bug_id="QA-TEST-0003c", failure_class=FailureClass.IMPLEMENTATION_BUG,
        stage="05_HOTSPOT_RADIUS", owner="hotspot-statistics",
        evidence_artifact="hotspot-statistics fixed scan.py; "
                          "test_fold_enrichment_is_never_negative now passes",
        files_or_components_affected=["src/hotspot3d/hotspot/scan.py"],
        targeted_test_required="test_fold_enrichment_is_never_negative",
        fix_applied_by_qa=False,  # QA validated, did not author the fix
        repair_attempt_number=1,
    )
    assert closed.validate() == []


# --------------------------------------------------------------------------- #
# CASE 4 / 5 — legitimate UNDERPOWERED / SCIENTIFIC_NEGATIVE: QA validates    #
# consistency and stops, never repairs methodology                            #
# --------------------------------------------------------------------------- #

def test_case4_underpowered_is_validated_not_repaired():
    # Shaped exactly like a real power_certificate.json posthoc block that FAILS.
    certificate_posthoc = {
        "N_P": 38, "N_B": 12, "m_test_family_size": 452, "B": 10000,
        "p_res": 9.999e-05, "p_comb_best": 7.1247e-04, "p_floor": 7.1247e-04,
        "c_1_rank1_critical_value": 1.1062e-04, "binding_floor": "combinatorial",
    }
    assert certificate_posthoc["p_floor"] > certificate_posthoc["c_1_rank1_critical_value"]

    brief = _minimal_brief(
        bug_id="QA-TEST-0004", failure_class=FailureClass.UNDERPOWERED,
        stage="06_FINAL_HOTSPOTS", owner="hotspot-statistics",
        observed_behavior="run terminated UNDERPOWERED with 0 significant centers",
        expected_behavior="(none — this is the correct terminal state for this design)",
        evidence_artifact=json.dumps(certificate_posthoc),
        why_this_is_an_implementation_defect="N/A — not a defect; validated as an "
                                             "internally consistent UNDERPOWERED result",
        scientific_methodology_change_required=False,
        fix_applied_by_qa=False,
    )
    assert brief.validate() == []
    # The defining property of this case: nothing about B, q, m, or any threshold
    # was touched to "fix" the result. Assert the evidence is exactly the
    # certificate's own numbers, not a narrative substitute.
    assert json.loads(brief.evidence_artifact)["p_floor"] > \
        json.loads(brief.evidence_artifact)["c_1_rank1_critical_value"]


def test_case5_scientific_negative_is_validated_not_repaired():
    certificate_posthoc = {
        "N_P": 6, "N_B": 45, "m_test_family_size": 40, "B": 10000,
        "p_res": 9.999e-05, "p_comb_best": 1.2e-07, "p_floor": 9.999e-05,
        "c_1_rank1_critical_value": 1.25e-03, "binding_floor": "permutation_resolution",
    }
    assert certificate_posthoc["p_floor"] <= certificate_posthoc["c_1_rank1_critical_value"]

    brief = _minimal_brief(
        bug_id="QA-TEST-0005", failure_class=FailureClass.SCIENTIFIC_NEGATIVE,
        stage="06_FINAL_HOTSPOTS", owner="hotspot-statistics",
        observed_behavior="0 significant centers; certificate PASSED",
        expected_behavior="(none — a valid negative)",
        evidence_artifact=json.dumps(certificate_posthoc),
        why_this_is_an_implementation_defect="N/A — validated negative, not a defect",
        fix_applied_by_qa=False,
    )
    assert brief.validate() == []
    assert json.loads(brief.evidence_artifact)["p_floor"] <= \
        json.loads(brief.evidence_artifact)["c_1_rank1_critical_value"]


# --------------------------------------------------------------------------- #
# CASE 6 — external-data failure: classified separately, QA repairs directly  #
# --------------------------------------------------------------------------- #

def test_case6_external_data_failure_is_qa_repairable_via_carveout():
    # The concrete carve-out this case is built on: retry/backoff constants
    # inside an otherwise off-limits (data-structure-owned) scientific file.
    assert qa_may_edit_directly("src/hotspot3d/data/sources.py", "CLINVAR_MAX_ATTEMPTS")
    assert qa_may_edit_directly("src/hotspot3d/data/sources.py", "CLINVAR_BACKOFF_BASE_S")
    # But the file as a whole, or an unrelated symbol in it, is NOT carved out.
    assert not qa_may_edit_directly("src/hotspot3d/data/sources.py")
    assert not qa_may_edit_directly("src/hotspot3d/data/sources.py", "SUBSTITUTION_TYPES")

    brief = _minimal_brief(
        bug_id="QA-TEST-0006", failure_class=FailureClass.EXTERNAL_DATA_FAILURE,
        stage="01_INPUT_RAW", owner="data-structure",
        observed_behavior="ClinVar download failed after 3 attempts",
        expected_behavior="download succeeds or exhausts a bounded retry budget "
                          "with a clear error",
        evidence_artifact=(
            "urllib3.exceptions.ReadTimeoutError: HTTPSConnectionPool"
            "(host='ftp.ncbi.nlm.nih.gov', port=443): Read timed out. "
            "3 attempts, all timed out after 30s — failure originates in the "
            "remote host's response time, not in hotspot3d's parsing"),
        why_this_is_an_implementation_defect=(
            "the retry BUDGET (CLINVAR_MAX_ATTEMPTS=3) is too low for this "
            "host's observed latency; raising it is an infra/retry-limit change, "
            "not a scientific one"),
        files_or_components_affected=["src/hotspot3d/data/sources.py:CLINVAR_MAX_ATTEMPTS"],
        fix_applied_by_qa=True,
    )
    assert brief.validate() == []
    assert_brief_acceptable_to_lead(brief)  # no escalation needed


# --------------------------------------------------------------------------- #
# CASE 7 — bounded repair loop stops a recurring signature                    #
# --------------------------------------------------------------------------- #

def test_case7_bounded_repair_loop_stops_recurrence(tmp_path):
    sig = root_cause_signature(
        stage="06_FINAL_HOTSPOTS", failure_class=FailureClass.IMPLEMENTATION_BUG,
        detail="center_index/universe_index confusion in detection.py")

    trail = RepairTrail(run_id="20250101T000000Z_qacase7a_qacase7b",
                        stage="06_FINAL_HOTSPOTS", max_attempts=2)
    trail.record_attempt(root_cause_signature=sig, bug_id="QA-TEST-0007a",
                         owner="hotspot-statistics", outcome="reappeared")
    assert not trail.exceeded()

    trail.record_attempt(root_cause_signature=sig, bug_id="QA-TEST-0007b",
                         owner="hotspot-statistics", outcome="reappeared")
    assert not trail.exceeded()  # exactly at the bound, not yet over it

    trail.record_attempt(root_cause_signature=sig, bug_id="QA-TEST-0007c",
                         owner="hotspot-statistics", outcome="reappeared")
    assert trail.exceeded(), "a third non-fixed attempt at the same signature " \
        "must trip the bound"
    payload = trail.as_dict()
    assert payload["QA_REPAIR_EXHAUSTED"] is True
    assert len(payload["attempts"]) == 3

    # Persists and reloads faithfully — the trail must survive across separate
    # agent invocations within one run, not just live in one process's memory.
    path = trail.write(tmp_path / "qa" / "triage_trail.json")
    reloaded = RepairTrail.read(path, run_id=trail.run_id, stage=trail.stage)
    assert reloaded.exceeded()
    assert len(reloaded.attempts) == 3


def test_case7_a_signature_marked_fixed_does_not_by_itself_trip_the_bound():
    sig = root_cause_signature(stage="08_FINAL_FOOTPRINT",
                               failure_class=FailureClass.IMPLEMENTATION_BUG,
                               detail="voxel grid off-by-one")
    trail = RepairTrail(run_id="20250101T000000Z_qacase7d_qacase7e",
                        stage="08_FINAL_FOOTPRINT", max_attempts=2)
    trail.record_attempt(root_cause_signature=sig, bug_id="QA-TEST-0007d",
                         owner="footprint-robustness", outcome="fixed")
    assert not trail.exceeded()

    # The SAME signature reappearing after being marked fixed is the case that
    # must escalate immediately (Section 8), not silently re-enter the loop.
    trail.record_attempt(root_cause_signature=sig, bug_id="QA-TEST-0007e",
                         owner="footprint-robustness", outcome="reappeared")
    assert trail.signature_reappeared(sig)


# --------------------------------------------------------------------------- #
# CASE 8 — the methodology firewall blocks a QA-attempted methodology change  #
# --------------------------------------------------------------------------- #

def test_case8_methodology_change_flag_blocks_a_qa_applied_fix():
    attempted = _minimal_brief(
        bug_id="QA-TEST-0008", failure_class=FailureClass.CONFIGURATION_PROBLEM,
        config_subclass=ConfigSubclass.SCIENTIFIC_CONFIG, config_key="fdr.q",
        stage="06_FINAL_HOTSPOTS", owner="lead",
        observed_behavior="0 significant centers at q=0.05",
        expected_behavior="(the analyst wants a positive result)",
        evidence_artifact="config_key=fdr.q value=0.05 subclass=SCIENTIFIC_CONFIG",
        root_cause="raising q to 0.10 would let 2 more centers clear BH",
        why_this_is_an_implementation_defect=(
            "N/A — this is not a defect at all, it is an attempted methodology "
            "change disguised as a fix, which is exactly what this test proves "
            "the firewall rejects"),
        scientific_methodology_change_required=True,
        fix_applied_by_qa=True,  # <- the firewall violation
    )
    problems = attempted.validate()
    assert problems
    assert any("QA must never" in p for p in problems)
    with pytest.raises(BlockedError):
        assert_brief_acceptable_to_lead(attempted)

    # The correct disposition of the SAME finding: escalate, never apply.
    escalated = _minimal_brief(
        bug_id="QA-TEST-0008b", failure_class=FailureClass.CONFIGURATION_PROBLEM,
        config_subclass=ConfigSubclass.SCIENTIFIC_CONFIG, config_key="fdr.q",
        stage="06_FINAL_HOTSPOTS", owner="lead",
        evidence_artifact="config_key=fdr.q value=0.05 subclass=SCIENTIFIC_CONFIG",
        scientific_methodology_change_required=True,
        fix_applied_by_qa=False,
    )
    assert escalated.validate() == []


@pytest.mark.parametrize("blocked_key", [
    "fdr.q", "permutation.B_default", "permutation.primary_null",
    "radius_domain.R_FLOOR", "radius_qc.significance_may_determine_admissibility",
    "clinvar.compound_significance_policy", "plddt.center_universe_min_plddt",
    "robustness.N_CAP", "robustness.max_wall_seconds", "structure.source",
    "seeding.MASTER_SEED",
])
def test_case8_every_firewall_listed_key_classifies_scientific(blocked_key):
    """Section 3's explicit list, re-checked mechanically against the manifest —
    not merely asserted in prose in the agent's own charter."""
    assert classify_config_key(blocked_key) is ConfigSubclass.SCIENTIFIC_CONFIG


# --------------------------------------------------------------------------- #
# CASE 9 — a SCIENTIFIC_CONFIG value misclassified as INFRA_CONFIG is caught  #
# --------------------------------------------------------------------------- #

def test_case9_misclassified_scientific_config_is_caught_not_repaired():
    # The worked trap from the manifest itself: structure.source reads like an
    # infrastructure "which data source" choice; it is II.13 FROZEN methodology.
    assert classify_config_key("structure.source") is ConfigSubclass.SCIENTIFIC_CONFIG

    misclassified = _minimal_brief(
        bug_id="QA-TEST-0009", failure_class=FailureClass.CONFIGURATION_PROBLEM,
        config_subclass=ConfigSubclass.INFRA_CONFIG,  # <- wrong, the violation
        config_key="structure.source",
        stage="03_STRUCTURE_QC", owner="lead",
        observed_behavior="structure.source could be pointed at a different "
                          "database to work around a retrieval failure",
        evidence_artifact="config_key=structure.source value=alphafold "
                          "subclass=INFRA_CONFIG",
        root_cause="mis-triaged as a data-source path setting",
        why_this_is_an_implementation_defect="N/A — this brief is itself the "
                                             "defect under test",
        fix_applied_by_qa=True,
    )
    problems = misclassified.validate()
    assert problems, "declaring structure.source as INFRA_CONFIG must be rejected"
    assert any("classify_config_key reports" in p for p in problems)
    with pytest.raises(BlockedError):
        assert_brief_acceptable_to_lead(misclassified)

    # Correctly classified, the same key routes to the owner instead.
    corrected = _minimal_brief(
        bug_id="QA-TEST-0009b", failure_class=FailureClass.CONFIGURATION_PROBLEM,
        config_subclass=ConfigSubclass.SCIENTIFIC_CONFIG, config_key="structure.source",
        stage="03_STRUCTURE_QC", owner="lead",
        evidence_artifact="config_key=structure.source value=alphafold "
                          "subclass=SCIENTIFIC_CONFIG",
        scientific_methodology_change_required=True,
        fix_applied_by_qa=False,
    )
    assert corrected.validate() == []


def test_case9_a_genuinely_infra_key_is_qa_editable():
    """The non-trap control case: execution.n_jobs really is INFRA_CONFIG."""
    assert classify_config_key("execution.n_jobs") is ConfigSubclass.INFRA_CONFIG
    brief = _minimal_brief(
        bug_id="QA-TEST-0009c", failure_class=FailureClass.CONFIGURATION_PROBLEM,
        config_subclass=ConfigSubclass.INFRA_CONFIG, config_key="execution.n_jobs",
        stage="00_RUN_SUMMARY", owner="lead",
        evidence_artifact="config_key=execution.n_jobs value=1 subclass=INFRA_CONFIG",
        fix_applied_by_qa=True,
    )
    assert brief.validate() == []
    assert_brief_acceptable_to_lead(brief)


# --------------------------------------------------------------------------- #
# Manifest sanity — the write-boundary and taxonomy stay well-formed          #
# --------------------------------------------------------------------------- #

def test_qa_write_manifest_globs_are_well_formed():
    assert QA_DIRECT_WRITE_GLOBS
    for pat in QA_DIRECT_WRITE_GLOBS:
        assert "\\" not in pat, f"{pat}: use forward slashes, checked cross-platform"


def test_symbol_carveouts_only_target_files_outside_the_direct_write_globs():
    """A carve-out only means something if the file ISN'T already fully
    QA-editable — otherwise the distinction between 'this symbol' and 'this
    file' collapses and the manifest stops being meaningfully restrictive."""
    for file_path in QA_SYMBOL_CARVEOUTS:
        assert not qa_may_edit_directly(file_path), (
            f"{file_path} is both a symbol carve-out AND fully QA-editable — "
            f"the carve-out is redundant, which likely means the file was "
            f"mis-added to QA_DIRECT_WRITE_GLOBS")


def test_owner_of_path_covers_every_stage_owner():
    from hotspot3d.utils.runctx import STAGE_OWNERS

    for stage, owner in STAGE_OWNERS.items():
        assert owner_of_path(f"FULL_RESULTS/{stage}/x.tsv") == owner
