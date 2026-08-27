"""End-to-end Stage E against the synthetic run and the mock source."""
from __future__ import annotations

import pytest

from hotspot3d.annotation.schema import (
    FORBIDDEN_UPSTREAM_KEYS,
    HANDOFF_05_PASSTHROUGH,
    MANDATORY_ROW_CONTEXT,
)
from hotspot3d.annotation.stage import REQUIRED_OUTPUTS, run_stage_e
from hotspot3d.orchestration.contracts import HANDOFF_REQUIRED_KEYS
from hotspot3d.utils.errors import BlockedError
from hotspot3d.utils.io import read_json, read_tsv

from tests.fixtures.synthetic_annotation import (
    CONFLICTING_STRONG_RESIDUE,
    NO_EVIDENCE_RESIDUE,
    WEAK_ONLY_RESIDUE,
    MockAnnotationSource,
    build_synthetic_run,
)

pytestmark = pytest.mark.integration


@pytest.fixture
def completed(tmp_path):
    run = build_synthetic_run(tmp_path)
    handoff = run_stage_e(
        run.ctx, handoff_02=run.handoff_02, handoff_03=run.handoff_03,
        handoff_04=run.handoff_04, source=MockAnnotationSource(),
    )
    return run, handoff


def stage_dir(run):
    return run.ctx.full_results / "10_ANNOTATION"


# --- outputs ----------------------------------------------------------------

def test_all_canonical_outputs_exist(completed):
    run, _ = completed
    for name in REQUIRED_OUTPUTS:
        assert (stage_dir(run) / name).is_file(), f"missing canonical output {name}"
    assert (stage_dir(run) / "figures").is_dir()


def test_stage_status_is_completed(completed):
    run, _ = completed
    status = read_json(stage_dir(run) / "stage_status.json")
    assert status["status"] == "COMPLETED"
    assert status["outcome_type"] == "COMPLETED"
    assert status["agent_owner"] == "biological-annotation"
    assert not (stage_dir(run) / "NOT_RUN.txt").exists()


def test_manifest_records_optional_figures_as_not_created(completed):
    run, _ = completed
    rows = read_tsv(stage_dir(run) / "stage_manifest.tsv")
    optional = [r for r in rows if r["kind"] == "optional_figure"]
    assert optional, "optional figures must be declared"
    assert all(r["status"] == "NOT_CREATED" and r["reason"] != "NA" for r in optional)
    required = [r for r in rows if r["kind"] == "required"]
    assert all(r["status"] == "CREATED" for r in required)


def test_provenance_is_written_to_the_single_owned_file(completed):
    run, _ = completed
    prov_dir = run.ctx.full_results / "12_REPRODUCIBILITY" / "provenance"
    assert [p.name for p in prov_dir.iterdir()] == ["provenance_10_annotation.json"]
    prov = read_json(prov_dir / "provenance_10_annotation.json")
    assert prov["agent_owner"] == "biological-annotation"
    assert prov["network_used"] is False
    assert prov["feeds_upstream"] is False
    assert prov["synthetic_input_mode"] is True


# --- verbatim statistics ----------------------------------------------------

def test_statistics_are_byte_identical_to_stage_b(completed):
    run, _ = completed
    source_rows = {
        r["hotspot_id"]: r
        for r in _rows(
            run.ctx.full_results / "06_FINAL_HOTSPOTS" / "hotspot_regions.tsv")
    }
    out_rows = {r["hotspot_id"]: r for r in _rows(stage_dir(run) / "hotspot_annotation.tsv")}
    assert source_rows and set(source_rows) == set(out_rows)

    for hid, src in source_rows.items():
        for field in ("n_plp", "n_blb", "fold_enrichment", "p_emp", "q_bh"):
            assert out_rows[hid][field] == src[field], (
                f"{hid}.{field} was not copied byte-for-byte"
            )


def _rows(path):
    """Raw-text rows, so the comparison is on bytes and not on parsed floats."""
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def test_plddt_is_read_not_recomputed(completed):
    run, _ = completed
    profile = {r["residue_index"]: r["plddt"] for r in _rows(
        run.ctx.full_results / "03_STRUCTURE_QC" / "plddt_profile.tsv")}
    for row in _rows(stage_dir(run) / "hotspot_annotation.tsv"):
        covered = [
            r["covered_residue_index"] for r in _rows(
                run.ctx.full_results / "06_FINAL_HOTSPOTS" / "hotspot_covered_residues.tsv")
            if r["hotspot_ids"] == row["hotspot_id"]
        ]
        expected_min = min(float(profile[c]) for c in covered)
        assert float(row["plddt_min"]) == expected_min


# --- mandatory context ------------------------------------------------------

def test_every_row_carries_all_flags_and_the_robustness_profile(completed):
    run, _ = completed
    rows = read_tsv(stage_dir(run) / "hotspot_annotation.tsv")
    assert rows
    for row in rows:
        for field in MANDATORY_ROW_CONTEXT:
            assert field in row
            assert row[field] not in (None, "", "NA"), f"{field} was dropped"


def test_three_residue_objects_are_reported_separately(completed):
    run, _ = completed
    rows = read_tsv(stage_dir(run) / "hotspot_annotation.tsv")
    for row in rows:
        assert row["n_significant_hotspot_centers"] is not None
        assert row["n_hotspot_sphere_classified_variants"] is not None
        assert row["n_hotspot_covered_residues"] is not None
    report = (stage_dir(run) / "stage_e_report.md").read_text(encoding="utf-8")
    assert "geometric test position and need not carry a ClinVar variant" in report


# --- mechanism curation -----------------------------------------------------

def test_the_three_rubric_cases_land_where_they_should(completed):
    run, _ = completed
    rows = {r["residue_index"]: r
            for r in read_tsv(stage_dir(run) / "functional_mechanism_variants.tsv")}

    mixed = rows[CONFLICTING_STRONG_RESIDUE]
    assert mixed["mechanism"] == "Mixed"
    assert "PMID:MOCK0001" in mixed["reference"] and "PMID:MOCK0002" in mixed["reference"]

    assert rows[WEAK_ONLY_RESIDUE]["mechanism"] == "Unclear"
    assert rows[NO_EVIDENCE_RESIDUE]["mechanism"] == "Not_Experimentally_Characterized"


def test_literature_log_records_inclusions_and_exclusions(completed):
    run, _ = completed
    rows = read_tsv(stage_dir(run) / "literature_search_log.tsv")
    decisions = {r["decision"] for r in rows}
    assert decisions == {"INCLUDED", "EXCLUDED"}
    for row in rows:
        assert row["query"] not in (None, "NA")
        assert row["database"] not in (None, "NA")
        assert row["search_date"] not in (None, "NA")
        assert row["hit_count"] is not None
        assert row["decision_reason"] not in (None, "NA")


def test_gaps_file_is_present_and_non_empty(completed):
    run, _ = completed
    gaps = read_tsv(stage_dir(run) / "annotation_gaps.tsv")
    assert gaps, "a missing annotation is recorded, never inferred"
    assert any(g["missing_item"] == "functional_mechanism" for g in gaps)


# --- the gate ---------------------------------------------------------------

def test_default_run_fails_the_gate_and_records_numbers(completed):
    run, handoff = completed
    posthoc = read_json(stage_dir(run) / "mechanism_spatial_posthoc.json")
    assert posthoc["status"] == "NOT_RUN"
    assert posthoc["gate"]["passed"] is False
    assert "counts_by_category_at_or_above_min_evidence" in posthoc["gate"]
    assert handoff.payload["posthoc_gate_passed"] is False
    assert handoff.payload["posthoc_gate_numbers"]["threshold_min_categories"] == 2


def test_gate_passing_run_executes_the_labelled_analysis(tmp_path):
    run = build_synthetic_run(tmp_path)
    handoff = run_stage_e(
        run.ctx, handoff_02=run.handoff_02, handoff_03=run.handoff_03,
        handoff_04=run.handoff_04, source=MockAnnotationSource(gate_passing=True),
    )
    posthoc = read_json(run.ctx.full_results / "10_ANNOTATION"
                        / "mechanism_spatial_posthoc.json")
    assert posthoc["analysis_type"] == "post_hoc"
    assert posthoc["feeds_upstream"] is False
    assert posthoc["correction"]["separate_from_stage_b"] is True
    assert handoff.payload["posthoc_gate_passed"] is True


# --- the handoff ------------------------------------------------------------

def test_handoff_shape_matches_the_contract(completed):
    _, handoff = completed
    assert handoff.name == "handoff_05"
    for key in HANDOFF_REQUIRED_KEYS["handoff_05"]:
        assert key in handoff.payload


def test_handoff_carries_no_upstream_parameter(completed):
    """The leakage audit, asserted on the artifact the Lead actually reads."""
    run, handoff = completed
    written = read_json(run.ctx.handoff_path("handoff_05"))

    def walk(node, trail=""):
        if isinstance(node, dict):
            for key, value in node.items():
                if trail.split(".")[0] in HANDOFF_05_PASSTHROUGH:
                    continue
                assert str(key).lower() not in FORBIDDEN_UPSTREAM_KEYS, (
                    f"handoff_05 exposes upstream parameter at {trail}.{key}"
                )
                walk(value, f"{trail}.{key}" if trail else str(key))
        elif isinstance(node, list):
            for i, value in enumerate(node):
                walk(value, f"{trail}[{i}]")

    walk(written)
    assert handoff.payload["terminal_stage"] is True
    assert handoff.payload["may_influence_discovery"] is False


def test_robustness_profile_is_an_exact_upstream_echo(completed):
    run, handoff = completed
    assert handoff.payload["robustness_profile"] == \
        run.handoff_04.payload["robustness_profile"]


def test_handoff_manifest_verifies(completed):
    run, handoff = completed
    handoff.verify(run.ctx)


# --- upstream negative ------------------------------------------------------

def test_upstream_negative_invents_no_regions(tmp_path):
    run = build_synthetic_run(tmp_path, negative=True)
    handoff = run_stage_e(
        run.ctx, handoff_02=run.handoff_02, handoff_03=run.handoff_03,
        handoff_04=run.handoff_04, source=MockAnnotationSource(),
    )
    sd = run.ctx.full_results / "10_ANNOTATION"

    status = read_json(sd / "stage_status.json")
    # NOT_RUN, not COMPLETED_NEGATIVE — the negative finding is the upstream stage's;
    # Stage E never became applicable and produced no negative result of its own.
    assert status["status"] == "NOT_RUN"
    assert status["outcome_type"] == "NOT_APPLICABLE"
    assert status["negative_result"]["condition"] == "NO_SIGNIFICANT_HOTSPOTS"

    assert (sd / "NOT_RUN.txt").is_file(), "no mysteriously empty folder ever exists"
    assert read_tsv(sd / "hotspot_annotation.tsv") == []
    assert handoff.payload["n_hotspots_annotated"] == 0
    assert handoff.negative_result is not None

    report = (sd / "stage_e_report.md").read_text(encoding="utf-8")
    assert "did not invent regions" in report


# --- blocking conditions ----------------------------------------------------

def test_unversioned_source_blocks_the_stage(tmp_path):
    run = build_synthetic_run(tmp_path)
    with pytest.raises(BlockedError, match="ANNOTATION_SOURCE_UNVERSIONED"):
        run_stage_e(
            run.ctx, handoff_02=run.handoff_02, handoff_03=run.handoff_03,
            handoff_04=run.handoff_04, source=MockAnnotationSource(versioned=False),
        )
    status = read_json(run.ctx.full_results / "10_ANNOTATION" / "stage_status.json")
    assert status["status"] == "BLOCKED"
    assert status["outcome_type"] == "TECHNICAL_FAILURE"


def test_missing_source_blocks_the_stage(tmp_path):
    """A synthetic run given no source must block: Stage E never fabricates
    annotation data. (A non-synthetic run with no source instead auto-
    constructs LiveAnnotationSource — see test_live_sources.py.)"""
    run = build_synthetic_run(tmp_path)
    with pytest.raises(BlockedError, match="never fabricates annotation data"):
        run_stage_e(
            run.ctx, handoff_02=run.handoff_02, handoff_03=run.handoff_03,
            handoff_04=run.handoff_04, source=None,
        )


def test_upstream_qc_fail_blocks_the_stage(tmp_path):
    run = build_synthetic_run(tmp_path)
    run.handoff_02.qc_status = "FAIL"
    with pytest.raises(BlockedError, match="qc_status=FAIL"):
        run_stage_e(
            run.ctx, handoff_02=run.handoff_02, handoff_03=run.handoff_03,
            handoff_04=run.handoff_04, source=MockAnnotationSource(),
        )


def test_tampered_upstream_file_blocks_the_stage(tmp_path):
    """Trust is never transitive: the manifest is recomputed against disk."""
    run = build_synthetic_run(tmp_path)
    target = run.ctx.full_results / "06_FINAL_HOTSPOTS" / "hotspot_regions.tsv"
    target.write_text(target.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(BlockedError, match="manifest verification failed"):
        run_stage_e(
            run.ctx, handoff_02=run.handoff_02, handoff_03=run.handoff_03,
            handoff_04=run.handoff_04, source=MockAnnotationSource(),
        )


def test_missing_robustness_profile_blocks_the_stage(tmp_path):
    run = build_synthetic_run(tmp_path)
    run.handoff_04.payload.pop("robustness_profile")
    with pytest.raises(BlockedError, match="no robustness_profile"):
        run_stage_e(
            run.ctx, handoff_02=run.handoff_02, handoff_03=run.handoff_03,
            handoff_04=run.handoff_04, source=MockAnnotationSource(),
        )


def test_missing_upstream_flag_blocks_the_stage(tmp_path):
    run = build_synthetic_run(tmp_path)
    run.handoff_02.payload.pop("fallback_radius_domain")
    with pytest.raises(BlockedError, match="do not state flag"):
        run_stage_e(
            run.ctx, handoff_02=run.handoff_02, handoff_03=run.handoff_03,
            handoff_04=run.handoff_04, source=MockAnnotationSource(),
        )


# --- nothing is written outside 10_ANNOTATION -------------------------------

def test_stage_e_writes_nowhere_else(completed):
    run, _ = completed
    allowed_prefixes = (
        "FULL_RESULTS/10_ANNOTATION/",
        "FULL_RESULTS/12_REPRODUCIBILITY/provenance/provenance_10_annotation.json",
    )
    upstream_written = []
    for path in run.ctx.run_root.rglob("*"):
        if not path.is_file():
            continue
        rel = str(path.relative_to(run.ctx.run_root)).replace("\\", "/")
        if rel.startswith(allowed_prefixes):
            continue
        # everything else must be an upstream artifact created by the fixture
        if rel not in {
            str(p.relative_to(run.ctx.run_root)).replace("\\", "/")
            for p in run.paths.values()
        }:
            upstream_written.append(rel)
    assert not upstream_written, f"Stage E wrote outside its scope: {upstream_written}"


def test_upstream_artifacts_are_unmodified(tmp_path):
    from hotspot3d.utils.hashing import sha256_file

    run = build_synthetic_run(tmp_path)
    before = {name: sha256_file(p) for name, p in run.paths.items()}
    run_stage_e(
        run.ctx, handoff_02=run.handoff_02, handoff_03=run.handoff_03,
        handoff_04=run.handoff_04, source=MockAnnotationSource(),
    )
    after = {name: sha256_file(p) for name, p in run.paths.items()}
    assert before == after, "nothing downstream may modify Stage A–D artifacts"


# --- reproducibility --------------------------------------------------------

def test_stage_is_deterministic(tmp_path):
    """Same run_id and same inputs give bit-identical outputs, seeded analysis included."""
    outputs = []
    for i in range(2):
        run = build_synthetic_run(
            tmp_path / f"run{i}", gene="MOCKG", run_id="20260101T000000Z_abcdef12_34567890"
        )
        run_stage_e(
            run.ctx, handoff_02=run.handoff_02, handoff_03=run.handoff_03,
            handoff_04=run.handoff_04, source=MockAnnotationSource(gate_passing=True),
        )
        sd = run.ctx.full_results / "10_ANNOTATION"
        outputs.append({
            name: (sd / name).read_text(encoding="utf-8")
            for name in ("hotspot_annotation.tsv", "functional_mechanism_variants.tsv",
                         "mechanism_by_hotspot.tsv", "literature_search_log.tsv",
                         "mechanism_spatial_posthoc.json")
        })
    assert outputs[0] == outputs[1]


def test_upstream_modified_mid_stage_is_blocking(tmp_path, monkeypatch):
    """If an input changes underneath Stage E, every verbatim copy is invalidated."""
    from hotspot3d.annotation import stage as stage_mod

    run = build_synthetic_run(tmp_path)
    target = run.ctx.full_results / "09_ROBUSTNESS" / "center_sensitivity.tsv"
    original = stage_mod._qc_plddt

    def tamper(rows, data):
        target.write_text(target.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        return original(rows, data)

    monkeypatch.setattr(stage_mod, "_qc_plddt", tamper)
    with pytest.raises(BlockedError, match="changed during Stage E"):
        run_stage_e(
            run.ctx, handoff_02=run.handoff_02, handoff_03=run.handoff_03,
            handoff_04=run.handoff_04, source=MockAnnotationSource(),
        )
