"""End-to-end Stage A: outputs, handoff shape, contract compliance, negatives."""
from __future__ import annotations

import json

import pytest

from hotspot3d.data.sources import (
    ClinVarClient,
    GeneResolution,
    LiveStageASource,
    UniProtClient,
    assert_release_metadata,
    assert_source_complete,
)
from hotspot3d.data.stage import EXPECTED_FILES, STAGES, run_stage_a
from hotspot3d.orchestration.contracts import (
    FORBIDDEN_DOWNSTREAM_COLUMNS,
    HANDOFF_REQUIRED_KEYS,
    QC_FAIL,
    QC_PASS,
    QC_PASS_WITH_WARNINGS,
    STAGE_B_PERMITTED_COLUMNS,
    assert_handoff_shape,
    load_handoff,
)
from hotspot3d.structure.sources import AlphaFoldClient
from hotspot3d.utils.errors import BlockedError, EscalationRequired
from hotspot3d.utils.io import read_json, read_tsv
from hotspot3d.utils.status import Status

pytestmark = pytest.mark.unit


@pytest.fixture
def run(make_ctx):
    def _run(case, **kwargs):
        ctx = make_ctx(**kwargs)
        return ctx, run_stage_a(ctx, gene=case.gene, source=case.source)

    return _run


# --- happy path -------------------------------------------------------------

def test_stage_a_completes_on_every_case(cases, run):
    for name, case in cases.items():
        ctx, handoff = run(case)
        expected = case.expected
        assert handoff.qc_status != QC_FAIL, name
        for key in ("M", "N", "N_P", "N_B", "n_conflict"):
            assert handoff.payload[key] == expected[key], f"{name}.{key}"


def test_every_canonical_output_exists(clustered, run):
    ctx, _ = run(clustered)
    for stage in STAGES:
        stage_dir = ctx.full_results / stage
        for name, _description in EXPECTED_FILES[stage]:
            assert (stage_dir / name).is_file(), f"{stage}/{name}"


def test_provenance_files_are_written_for_stage_a_only(clustered, run):
    ctx, _ = run(clustered)
    provenance = ctx.full_results / "12_REPRODUCIBILITY" / "provenance"
    names = sorted(p.name for p in provenance.glob("*.json"))
    assert names == ["provenance_01_input_raw.json", "provenance_02_clinvar.json",
                     "provenance_03_structure_qc.json"]

    payload = read_json(provenance / "provenance_02_clinvar.json")
    assert payload["derived_seeds"] is None
    assert payload["derived_seeds_reason"] == "not_applicable_deterministic_stage"
    assert payload["agent_owner"] == "data-structure"
    assert payload["information_barrier"]["downstream_directories_read"] == []
    assert payload["star_filter_audit"]["passed"] is True


def test_stage_a_writes_nothing_downstream(clustered, run):
    ctx, _ = run(clustered)
    written = {p.relative_to(ctx.full_results).parts[0]
               for p in ctx.full_results.rglob("*") if p.is_file()}
    assert written <= set(STAGES) | {"12_REPRODUCIBILITY"}


def test_run_root_is_left_partial_for_the_lead_to_rename(clustered, run):
    ctx, _ = run(clustered)
    assert ctx.run_root.name.endswith(".partial")


# --- handoff ----------------------------------------------------------------

def test_handoff_shape_and_column_contract(clustered, run):
    ctx, handoff = run(clustered)
    assert handoff.name == "handoff_01"
    assert_handoff_shape(handoff)
    for key in HANDOFF_REQUIRED_KEYS["handoff_01"]:
        assert key in handoff.payload

    assert handoff.payload["stage_b_permitted_columns"] == list(STAGE_B_PERMITTED_COLUMNS)
    assert handoff.payload["forbidden_downstream"] == list(FORBIDDEN_DOWNSTREAM_COLUMNS)
    assert handoff.payload["sensitivity_channel"]["may_redefine_primary"] is False
    assert handoff.payload["derived_seeds"] is None


def test_handoff_manifest_verifies_and_covers_the_outputs(clustered, run):
    ctx, handoff = run(clustered)
    handoff.verify(ctx)                       # recomputes every hash

    assert handoff.manifest
    covered = set(handoff.manifest)
    for stage in STAGES:
        for name, _ in EXPECTED_FILES[stage]:
            if name in ("stage_manifest.tsv", "handoff_01.json"):
                continue
            assert f"FULL_RESULTS/{stage}/{name}" in covered, name


def test_handoff_reloads_from_disk_identically(clustered, run):
    ctx, handoff = run(clustered)
    reloaded = load_handoff(ctx, "handoff_01")
    assert reloaded.run_id == handoff.run_id
    assert reloaded.qc_status == handoff.qc_status
    assert reloaded.manifest == handoff.manifest
    assert reloaded.payload["N_P"] == handoff.payload["N_P"]


def test_manifest_verification_detects_tampering(clustered, run):
    ctx, handoff = run(clustered)
    target = ctx.full_results / "03_STRUCTURE_QC" / "classified_cohort.tsv"
    target.write_text(target.read_text() + "\n", encoding="utf-8")
    with pytest.raises(BlockedError, match="manifest verification failed"):
        handoff.verify(ctx)


# --- Output Contract compliance --------------------------------------------

def test_cohort_file_exposes_no_forbidden_column(clustered, run):
    ctx, _ = run(clustered)
    path = ctx.full_results / "03_STRUCTURE_QC" / "classified_cohort.tsv"
    rows = read_tsv(path, allowed_columns=STAGE_B_PERMITTED_COLUMNS,
                    forbidden_columns=FORBIDDEN_DOWNSTREAM_COLUMNS)
    assert rows
    assert set(rows[0]) == set(STAGE_B_PERMITTED_COLUMNS) | {"schema_version"}


def test_universe_file_matches_the_declared_column_set(clustered, run):
    ctx, handoff = run(clustered)
    path = ctx.full_results / "03_STRUCTURE_QC" / "positional_universe.tsv"
    rows = read_tsv(path, forbidden_columns=FORBIDDEN_DOWNSTREAM_COLUMNS)
    declared = set(handoff.payload["stage_b_universe_columns"])
    assert set(rows[0]) == declared | {"schema_version"}


#: The write-once payload in 01_INPUT_RAW is preserved exactly as the source
#: delivered it. Stamping a schema_version column onto it would be editing the
#: evidence, so the derived-table conventions deliberately do not apply to it.
IMMUTABLE_RAW_PAYLOADS = {"clinvar_raw.tsv"}


def test_tsv_conventions_are_honoured(clustered, run):
    ctx, _ = run(clustered)
    for stage in STAGES:
        for path in (ctx.full_results / stage).glob("*.tsv"):
            raw = path.read_bytes().decode("utf-8")
            assert "\r" not in raw, path
            if path.name in IMMUTABLE_RAW_PAYLOADS:
                continue
            header = raw.splitlines()[0].split("\t")
            assert "schema_version" in header, path
            for line in raw.splitlines()[1:]:
                assert "" not in line.split("\t"), f"empty cell in {path}"


def test_raw_payload_is_preserved_unmodified(clustered, run):
    """01_INPUT_RAW holds evidence, not a Stage A product."""
    ctx, _ = run(clustered)
    header = (ctx.full_results / "01_INPUT_RAW" / "clinvar_raw.tsv"
              ).read_text().splitlines()[0].split("\t")
    assert "schema_version" not in header
    assert header[:3] == list(clustered.records[0])[:3]


def test_stage_status_is_written_for_every_stage(clustered, run):
    ctx, _ = run(clustered)
    for stage in STAGES:
        payload = read_json(ctx.full_results / stage / "stage_status.json")
        assert payload["status"] == Status.COMPLETED.value
        assert payload["agent_owner"] == "data-structure"
        assert payload["n_outputs_created"] > 0
        assert not (ctx.full_results / stage / "NOT_RUN.txt").exists()


def test_stage_manifest_lists_expected_files_with_hashes(clustered, run):
    ctx, _ = run(clustered)
    for stage in STAGES:
        rows = read_tsv(ctx.full_results / stage / "stage_manifest.tsv")
        assert {r["path"] for r in rows}
        listed = {r["path"].rsplit(f"{stage}/", 1)[-1] for r in rows}
        assert listed == {name for name, _ in EXPECTED_FILES[stage]}
        for row in rows:
            assert row["status"] in ("CREATED", "NOT_CREATED", "DEFERRED")
            if row["status"] == "CREATED":
                assert len(row["sha256"]) == 64


def test_cohort_summary_records_the_audits(clustered, run):
    ctx, _ = run(clustered)
    payload = read_json(ctx.full_results / "02_CLINVAR" / "cohort_summary.json")
    assert payload["star_filter_audit"]["passed"] is True
    assert payload["review_stars_used_for_inclusion"] is False
    assert payload["vus_or_conflicting_admitted_to_binary_classes"] is False
    assert payload["plddt_used_as_filter"] is False
    assert payload["record_conservation"]["conserved"] is True
    # The shipped ruling is whole_field, so component matching is OFF: only an
    # exact bare-label match enters the binary classes. Pinned to the frozen value
    # rather than to whichever reading happens to be configured, so a silent revert
    # to component_wise fails here as well as in the Lead audit guard.
    assert payload["significance_policy"]["component_matching"] is False


def test_cohort_summary_reports_n_b_prominently_and_the_v2_additions(clustered, run):
    """Workflow v2 §1 — N_B at the top level (not only nested), the
    compound-significance policy counterfactual cost, and the true-cause
    breakdown nested under the two reasons it still groups together.
    """
    ctx, handoff = run(clustered)
    payload = read_json(ctx.full_results / "02_CLINVAR" / "cohort_summary.json")

    assert payload["N_B"] == payload["residue_counts"]["N_B"] == handoff.payload["N_B"]
    assert payload["N_P"] == handoff.payload["N_P"]
    assert payload["N"] == handoff.payload["N"]
    assert "benign_cohort_adequacy_note" in payload

    cost = payload["compound_significance_policy_cost"]
    assert cost["configured_policy"] == "whole_field"
    assert cost["counterfactual_policy"] == "component_wise"
    assert cost["n_excluded_by_configured_but_kept_by_counterfactual"] >= 0

    detail = payload["exclusion_reason_detail_counts"]
    # Every reason key present here must itself be one of the two that cover
    # more than one true cause; every nested count must sum back to the
    # top-level exclusion_reason_counts for that same reason.
    for reason, per_cause in detail.items():
        assert reason in ("unparseable_hgvs", "non_missense_variant_type")
        assert sum(per_cause.values()) == payload["exclusion_reason_counts"][reason]


def test_raw_inputs_carry_a_sha256_sidecar(clustered, run):
    from hotspot3d.utils.hashing import sha256_file

    ctx, _ = run(clustered)
    raw = ctx.full_results / "01_INPUT_RAW"
    for name in ("clinvar_raw.tsv", "uniprot_canonical.fasta", "alphafold_model.cif"):
        digest = (raw / f"{name}.sha256").read_text().split()[0]
        assert digest == sha256_file(raw / name)

    log = read_json(raw / "retrieval_log.json")
    assert {entry["artifact"] for entry in log["retrievals"]} == {
        "clinvar_raw.tsv", "uniprot_canonical.fasta", "alphafold_model.cif"}
    assert all(entry["release_date"] for entry in log["retrievals"])


def test_structure_export_is_not_edited(clustered, run):
    ctx, _ = run(clustered)
    payload = read_json(ctx.full_results / "03_STRUCTURE_QC" / "structure_qc.json")
    assert payload["structure_export"]["model_edited"] is False
    assert payload["structure_export"]["b_factor_column"] == "pLDDT"
    assert payload["structure_edited"] is False
    assert payload["sequence_trimmed"] is False


def test_report_is_written_with_the_reporting_protocol_headings(clustered, run):
    ctx, _ = run(clustered)
    text = (ctx.full_results / "03_STRUCTURE_QC" / "stage_a_report.md").read_text()
    for heading in ("STATUS:", "INPUTS USED:", "METHODS EXECUTED:", "OUTPUTS GENERATED:",
                    "QC RESULTS:", "SCIENTIFIC DECISIONS:", "WARNINGS:",
                    "FAILED OR REJECTED ANALYSES:", "UNRESOLVED ISSUES:", "HANDOFF:"):
        assert heading in text, heading


# --- conflicts, unmapped residues, negatives --------------------------------

def test_conflict_case_warns_and_writes_the_conflict_table(cases, run):
    ctx, handoff = run(cases["conflict"])
    rows = read_tsv(ctx.full_results / "02_CLINVAR" / "residue_class_conflicts.tsv")
    assert len(rows) == cases["conflict"].expected["n_conflict"] == handoff.payload["n_conflict"]
    for row in rows:
        assert row["exclusion_reason"] == "residue_class_conflict"
        assert row["n_records_plp"] > 0 and row["n_records_blb"] > 0

    warnings = read_tsv(ctx.full_results / "02_CLINVAR" / "warnings_02_clinvar.tsv")
    codes = {w["warning_code"] for w in warnings}
    assert "RESIDUE_CLASS_CONFLICT" in codes


def test_unusable_ca_residue_is_excluded_and_reported(clustered, run):
    victim = clustered.expected["plp_residues"][0]
    ctx, handoff = run(clustered.with_unusable_ca(victim))

    assert handoff.payload["M"] == clustered.expected["M"] - 1
    assert handoff.payload["N"] == clustered.expected["N"] - 1
    assert handoff.payload["n_unmapped"] == 1

    report = read_json(ctx.full_results / "03_STRUCTURE_QC" / "mapping_report.json")
    assert [r["residue_index"] for r in report["unmapped_residues"]] == [victim]
    assert report["unmapped_residues"][0]["reason"] == "no_ca_coordinate"

    excluded = read_tsv(ctx.full_results / "02_CLINVAR" / "variants_excluded_from_primary.tsv")
    assert any(r["exclusion_reason"] == "no_ca_coordinate" and r["residue_index"] == victim
               for r in excluded)


@pytest.mark.parametrize("dropped", ["PLP", "BLB"])
def test_empty_class_blocks_the_handoff_but_keeps_the_evidence(cases, run, dropped):
    """Case (B): the INPUT cannot support the analysis, so nothing was tested.

    An empty class is BLOCKING in the frozen enum (IX.7) and blocks the A->B
    handoff. Reporting PASS here would assert "we looked and found nothing" when
    in fact we could not look. ``outcome_type`` stays SCIENTIFIC_NEGATIVE because
    the cause is ClinVar's content rather than a fault in the code.
    """
    case = cases["sparse"].without_class(dropped)
    ctx, handoff = run(case)

    assert handoff.qc_status == QC_FAIL
    assert handoff.negative_result["condition"] == "INSUFFICIENT_CLASSIFIED_RESIDUES"
    assert handoff.negative_result["outcome_type"] == "SCIENTIFIC_NEGATIVE"
    assert handoff.negative_result["nothing_was_tested"] is True

    warnings = read_tsv(ctx.full_results / "02_CLINVAR" / "warnings_02_clinvar.tsv")
    blocking = {w["warning_code"] for w in warnings if w["severity"] == "BLOCKING"}
    assert "INSUFFICIENT_CLASSIFIED_RESIDUES" in blocking, (
        "the frozen BLOCKING severity must not be softened")

    status = read_json(ctx.full_results / "02_CLINVAR" / "stage_status.json")
    assert status["status"] == Status.BLOCKED.value
    assert status["outcome_type"] == "SCIENTIFIC_NEGATIVE"

    # NOT_RUN.txt must name the cause without implying anything was tested.
    not_run = (ctx.full_results / "02_CLINVAR" / "NOT_RUN.txt").read_text()
    assert "No spatial analysis was attempted" in not_run
    assert "nothing was tested" in not_run

    # the evidence is still complete: nothing was relaxed to rescue the result
    assert (ctx.full_results / "03_STRUCTURE_QC" / "positional_universe.tsv").is_file()
    assert handoff.payload["M"] == case.expected["M"]


def test_a_blocked_handoff_stops_downstream(cases, run):
    """qc_status FAIL is the interlock: no consumer may proceed against it."""
    ctx, handoff = run(cases["sparse"].without_class("BLB"))
    with pytest.raises(BlockedError, match="qc_status=FAIL"):
        handoff.verify(ctx)


def test_a_viable_cohort_with_a_warning_reports_pass_with_warnings(clustered, run):
    """A warning degrades the verdict on its own strength, with no outcome override.

    The cohort here is viable, so there is no negative result in play: this pins
    the ordinary PASS_WITH_WARNINGS path that the insufficient-cohort branch must
    not be confused with.
    """
    case = clustered.with_unusable_ca(clustered.expected["plp_residues"][0])
    ctx, handoff = run(case)

    assert handoff.negative_result is None
    assert handoff.qc_status == QC_PASS_WITH_WARNINGS
    assert handoff.payload["n_unmapped"] == 1
    assert handoff.payload["M"] == clustered.expected["M"] - 1

    codes = {w["warning_code"] for w in read_tsv(
        ctx.full_results / "03_STRUCTURE_QC" / "warnings_03_structure_qc.tsv")}
    assert "STRUCTURAL_CONFIDENCE_SENSITIVE" in codes

    status = read_json(ctx.full_results / "03_STRUCTURE_QC" / "stage_status.json")
    assert status["status"] == Status.COMPLETED.value


def test_a_clean_viable_cohort_reports_pass(cases, run):
    ctx, handoff = run(cases["clustered"])
    assert handoff.qc_status == QC_PASS
    assert handoff.negative_result is None
    handoff.verify(ctx)


def test_multi_fragment_entry_blocks_and_still_writes_status(clustered, make_ctx):
    case = clustered.with_fragments(2)
    ctx = make_ctx()
    with pytest.raises(BlockedError, match="MULTI_FRAGMENT_AFDB_ENTRY"):
        run_stage_a(ctx, gene=case.gene, source=case.source)

    for stage in STAGES:
        payload = read_json(ctx.full_results / stage / "stage_status.json")
        assert payload["status"] == Status.BLOCKED.value
        assert payload["outcome_type"] == "TECHNICAL_FAILURE"
        assert (ctx.full_results / stage / "NOT_RUN.txt").is_file()


# --- preconditions and sources ---------------------------------------------

def test_rerunning_into_a_populated_directory_is_blocked(clustered, make_ctx):
    ctx = make_ctx()
    run_stage_a(ctx, gene=clustered.gene, source=clustered.source)
    with pytest.raises(BlockedError, match="already contains"):
        run_stage_a(ctx, gene=clustered.gene, source=clustered.source)


def test_gene_mismatch_with_the_run_context_is_blocked(clustered, make_ctx):
    ctx = make_ctx(gene="OTHERGENE")
    with pytest.raises(BlockedError, match="one gene"):
        run_stage_a(ctx, gene=clustered.gene, source=clustered.source)


def test_synthetic_run_without_a_source_is_blocked(make_ctx):
    ctx = make_ctx()
    with pytest.raises(BlockedError, match="never fabricates"):
        run_stage_a(ctx, gene="SYNGENE", source=None)


def test_live_source_is_enabled_and_gated_by_allow_network():
    """Supersedes the implementation-phase test that refused construction.

    The bundle now exists; the gate is ``allow_network``, not the constructor, so
    a run with the gate closed refuses the retrieval rather than the object.
    Composition, delegation and the config wiring are covered in
    ``test_live_source.py``.
    """
    source = LiveStageASource(allow_network=False)
    assert_source_complete(source)
    with pytest.raises(BlockedError, match="allow_network=FALSE"):
        source.fetch_missense("SYNGENE")


def test_incomplete_source_is_blocked():
    class Partial:
        def fetch_missense(self, gene):
            return []

    with pytest.raises(BlockedError, match="does not implement"):
        assert_source_complete(Partial())


def test_undated_payload_is_blocked():
    with pytest.raises(BlockedError, match="undated or unversioned"):
        assert_release_metadata({"source": "x", "query": "y"}, "ClinVar")


# --- live clients: request shaping and parsing, without a network -----------

def test_clinvar_request_carries_no_star_filter():
    request = ClinVarClient().build_request("BRAF")
    assert request["review_star_filter"] is None
    assert request["url"].endswith("variant_summary.txt.gz")


def test_network_calls_are_gated():
    with pytest.raises(BlockedError, match="allow_network=FALSE"):
        ClinVarClient(allow_network=False)._download("https://example.invalid/x")
    with pytest.raises(BlockedError, match="allow_network=FALSE"):
        UniProtClient(allow_network=False)._get_json("https://example.invalid/x", {})
    with pytest.raises(BlockedError, match="allow_network=FALSE"):
        AlphaFoldClient(allow_network=False)._get("https://example.invalid/x", "text/plain")


def test_uniprot_entry_parsing_and_mane_selection():
    entry = {
        "primaryAccession": "P15056", "uniProtkbId": "BRAF_HUMAN",
        "entryAudit": {"entryVersion": "245"},
        "sequence": {"value": "MAALSGGGGG"},
        "uniProtKBCrossReferences": [{"database": "MANE-Select", "id": "NM_004333.6"}],
        "proteinDescription": {"recommendedName": {"fullName": {"value": "BRAF kinase"}}},
    }
    resolution = UniProtClient.parse_entry(entry, "BRAF")
    assert isinstance(resolution, GeneResolution)
    assert resolution.mane_transcript == "NM_004333.6"
    assert resolution.mane_fallback_used is False
    assert resolution.aa_at(2) == "A"
    assert resolution.aa_at(0) is None and resolution.aa_at(999) is None


def test_ambiguous_gene_resolution_escalates():
    payload = {"results": [{"primaryAccession": "P1"}, {"primaryAccession": "P2"}]}
    with pytest.raises(EscalationRequired):
        UniProtClient.parse_search_payload(payload, "AMBIG")
    with pytest.raises(EscalationRequired):
        UniProtClient.parse_search_payload({"results": []}, "MISSING")


def test_afdb_prediction_parsing_counts_fragments():
    meta = AlphaFoldClient.parse_prediction_payload(
        [{"entryId": "AF-P15056-F1", "latestVersion": 4, "modelCreatedDate": "2022-06-01",
          "cifUrl": "https://example/AF.cif"}], "P15056")
    assert meta["n_fragments"] == 1
    assert meta["model_version"] == "4"

    two = AlphaFoldClient.parse_prediction_payload(
        [{"entryId": "AF-X-F1", "latestVersion": 4}, {"entryId": "AF-X-F2"}], "X")
    assert two["n_fragments"] == 2


def test_afdb_missing_entry_escalates():
    with pytest.raises(EscalationRequired):
        AlphaFoldClient.parse_prediction_payload([], "NOSUCH")


# --- determinism ------------------------------------------------------------

def test_two_runs_produce_identical_scientific_outputs(clustered, make_ctx):
    payloads = []
    tables = []
    for _ in range(2):
        ctx = make_ctx()
        handoff = run_stage_a(ctx, gene=clustered.gene, source=clustered.source)
        payloads.append({k: handoff.payload[k]
                         for k in ("M", "N", "N_P", "N_B", "n_conflict", "D_max")})
        tables.append(
            (ctx.full_results / "03_STRUCTURE_QC" / "classified_cohort.tsv").read_text())
    assert payloads[0] == payloads[1]
    assert tables[0] == tables[1]
