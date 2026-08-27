"""Structural guardrails: mandatory row context, write ownership, no back-propagation."""
from __future__ import annotations

import pytest

from hotspot3d.annotation.schema import (
    MANDATORY_ROW_CONTEXT,
    UPSTREAM_FLAGS,
    assert_mandatory_context,
    assert_no_upstream_influence,
    assert_write_target,
)
from hotspot3d.annotation.sources import (
    NetworkAnnotationSource,
    SourceVersion,
    validate_source,
)
from hotspot3d.utils.errors import BlockedError, LeakageError

pytestmark = pytest.mark.unit


def good_row(**overrides) -> dict:
    row = {
        "hotspot_id": "HS1",
        "global_clustering_flag": "significant",
        "permutation_resolution_limited": False,
        "fallback_radius_domain": False,
        "domain_boundary_warning": False,
        "robustness_profile_summary": "preservation_rate=0.82;mean_jaccard=0.79",
        "per_center_influence_summary": "104=0.10;108=0.15",
        "center_recurrence_rate": "104=0.95;108=0.92",
    }
    row.update(overrides)
    return row


# --- robustness and flags travel with every claim ---------------------------

def test_complete_row_passes():
    assert_mandatory_context([good_row()]) is None


@pytest.mark.parametrize("field", MANDATORY_ROW_CONTEXT)
def test_missing_mandatory_field_cannot_be_written(field):
    row = good_row()
    del row[field]
    with pytest.raises(BlockedError, match="missing mandatory context"):
        assert_mandatory_context([row])


@pytest.mark.parametrize("field", MANDATORY_ROW_CONTEXT)
def test_na_mandatory_field_cannot_be_written(field):
    """'NA' and 'absent' must not be indistinguishable for a flag."""
    with pytest.raises(BlockedError, match="empty/NA mandatory context"):
        assert_mandatory_context([good_row(**{field: "NA"})])


def test_all_four_upstream_flags_are_mandatory():
    assert set(UPSTREAM_FLAGS) <= set(MANDATORY_ROW_CONTEXT)


def test_false_flag_is_a_legal_value():
    """FALSE is information; only absence is forbidden."""
    assert_mandatory_context([good_row(domain_boundary_warning=False)])


# --- write ownership --------------------------------------------------------

def test_annotation_directory_is_writable(tmp_path):
    target = tmp_path / "FULL_RESULTS" / "10_ANNOTATION" / "hotspot_annotation.tsv"
    assert assert_write_target(tmp_path, target) == target.resolve()


def test_own_provenance_file_is_writable(tmp_path):
    target = (tmp_path / "FULL_RESULTS" / "12_REPRODUCIBILITY" / "provenance"
              / "provenance_10_annotation.json")
    assert assert_write_target(tmp_path, target) == target.resolve()


@pytest.mark.parametrize(
    "stage_dir",
    ["01_INPUT_RAW", "02_CLINVAR", "03_STRUCTURE_QC", "04_GLOBAL_CLUSTERING",
     "05_HOTSPOT_RADIUS", "06_FINAL_HOTSPOTS", "07_FOOTPRINT_RADIUS",
     "08_FINAL_FOOTPRINT", "09_ROBUSTNESS", "11_SENSITIVITY"],
)
def test_upstream_stage_directories_are_refused(tmp_path, stage_dir):
    target = tmp_path / "FULL_RESULTS" / stage_dir / "sneaky.tsv"
    with pytest.raises(LeakageError, match="OWNERSHIP violation"):
        assert_write_target(tmp_path, target)


@pytest.mark.parametrize(
    "rel",
    ["FULL_RESULTS/00_RUN_SUMMARY/run_summary.json",
     "FULL_RESULTS/12_REPRODUCIBILITY/manifest.tsv",
     "REVIEW_PACK/index.html", "ARCHIVES/bundle.tar.gz", "manifest.tsv"],
)
def test_lead_owned_paths_are_refused(tmp_path, rel):
    with pytest.raises(LeakageError, match="OWNERSHIP violation"):
        assert_write_target(tmp_path, tmp_path / rel)


def test_writing_outside_the_run_root_is_refused(tmp_path):
    with pytest.raises(LeakageError, match="outside the run root"):
        assert_write_target(tmp_path / "run", tmp_path / "elsewhere" / "x.tsv")


# --- no back-propagation ----------------------------------------------------

def minimal_payload(**extra) -> dict:
    payload = {
        "annotation_sources": ["dbNSFP"],
        "annotation_source_versions": {},
        "n_hotspots_annotated": 2,
        "mechanism_counts": {"GOF": 1},
        "posthoc_gate_passed": False,
        "posthoc_gate_numbers": {"passed": False},
        "robustness_profile": {"preservation_rate": 0.82},
        "global_clustering_flag": "significant",
        "interpretation_caveats": ["overlay only"],
        "upstream_handoff_hashes": {},
        "residue_objects": {},
        "terminal_stage": True,
        "may_influence_discovery": False,
    }
    payload.update(extra)
    return payload


def test_clean_payload_passes():
    assert assert_no_upstream_influence(minimal_payload()) is None


@pytest.mark.parametrize(
    "key", ["r_hot", "r_fp", "hotspot_radius", "footprint_radius",
            "recommended_radius", "preservation_threshold"],
)
def test_top_level_upstream_parameter_is_refused(key):
    with pytest.raises(BlockedError, match="undeclared key"):
        assert_no_upstream_influence(minimal_payload(**{key: 8.5}))


def test_nested_upstream_parameter_is_refused():
    """A radius suggestion buried inside an allowed key is still back-propagation."""
    payload = minimal_payload(
        interpretation_caveats=[{"note": "widen", "r_hot": 12.0}]
    )
    with pytest.raises(LeakageError, match="BACK-PROPAGATION violation"):
        assert_no_upstream_influence(payload)


def test_deeply_nested_upstream_parameter_is_refused():
    payload = minimal_payload(
        posthoc_gate_numbers={"detail": {"suggested_radius": 9.0}}
    )
    with pytest.raises(LeakageError, match="BACK-PROPAGATION violation"):
        assert_no_upstream_influence(payload)


def test_upstream_passthrough_is_exempt():
    """Echoing Stage D's own profile back to the Lead is not authoring a parameter."""
    payload = minimal_payload(
        robustness_profile={"footprint_radius": 6.0, "preservation_rate": 0.82}
    )
    assert assert_no_upstream_influence(payload) is None


# --- source version pinning -------------------------------------------------

class _Bare:
    pass


def test_incomplete_source_surface_is_blocking():
    with pytest.raises(BlockedError, match="does not implement"):
        validate_source(_Bare())


def test_unversioned_source_is_blocking(tmp_path):
    from tests.fixtures.synthetic_annotation import MockAnnotationSource
    with pytest.raises(BlockedError, match="ANNOTATION_SOURCE_UNVERSIONED"):
        validate_source(MockAnnotationSource(versioned=False))


@pytest.mark.parametrize("bad", ["", "  ", "NA", "unknown", "latest"])
def test_placeholder_versions_do_not_count_as_pinned(bad):
    assert not SourceVersion("dbNSFP", bad).is_versioned()


def test_pinned_version_counts():
    assert SourceVersion("dbNSFP", "4.9a").is_versioned()


# --- the network gate -------------------------------------------------------

def test_network_calls_are_gated_when_network_is_off():
    source = NetworkAnnotationSource(
        allow_network=False, versions={"UniProt": SourceVersion("UniProt", "2026_01")}
    )
    for call in (
        lambda: source.domains("P00000"),
        lambda: source.functional_sites("P00000"),
        lambda: source.conservation("P00000", [1, 2]),
        lambda: source.disease_associations("GENE"),
        lambda: source.literature_mechanisms("GENE"),
        lambda: source.secondary_structure("model.cif"),
    ):
        with pytest.raises(BlockedError, match="allow_network is FALSE"):
            call()
    assert source.call_log == [], "no request may be made when the gate is closed"


def test_no_live_retrieval_is_implemented():
    """The gate is the only implemented part; there is no code path to a real call."""
    source = NetworkAnnotationSource(
        allow_network=True, versions={"UniProt": SourceVersion("UniProt", "2026_01")}
    )
    with pytest.raises(NotImplementedError, match="not implemented in this build"):
        source.domains("P00000")
