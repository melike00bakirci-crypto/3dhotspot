"""The fixture contract other agents build their integration tests against.

If any of these fail, a teammate's suite breaks, so they are asserted here rather
than left to discover downstream.
"""
from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path

import pytest

from hotspot3d.data.synthetic import (
    CASE_NAMES,
    SyntheticCase,
    SyntheticStageASource,
    all_synthetic_cases,
    make_synthetic_case,
    make_synthetic_source,
)
from hotspot3d.data.sources import assert_source_complete

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_both_documented_import_paths_resolve_to_the_same_objects():
    package = importlib.import_module("hotspot3d.data.synthetic")
    mirror = importlib.import_module("tests.fixtures.synthetic_clinvar")
    assert mirror.make_synthetic_case is package.make_synthetic_case
    assert mirror.SyntheticCase is package.SyntheticCase
    assert set(mirror.CASE_NAMES) == set(package.CASE_NAMES)


@pytest.mark.parametrize("name", CASE_NAMES)
def test_case_exposes_records_structure_expected(name):
    case = make_synthetic_case(name)
    assert isinstance(case, SyntheticCase)
    assert isinstance(case.records, list) and case.records
    assert all(isinstance(r, dict) for r in case.records)
    assert case.structure.n_residues == len(case.sequence)
    assert isinstance(case.expected, dict) and case.expected["case"] == name

    # positional construction, as written in the build plan
    rebuilt = SyntheticCase(case.records, case.structure, case.expected)
    assert rebuilt.records is case.records


@pytest.mark.parametrize("name", CASE_NAMES)
def test_case_source_satisfies_the_stage_a_provider_interface(name):
    source = make_synthetic_source(name)
    assert isinstance(source, SyntheticStageASource)
    assert assert_source_complete(source) is None

    resolution = source.resolve_gene("SYNGENE")
    assert resolution.sequence_length == len(source.case.sequence)
    assert source.fetch_model(resolution.uniprot_acc) is source.case.structure
    assert len(source.fetch_missense("SYNGENE")) == len(source.case.records)

    for metadata in (source.release_metadata(), source.sequence_metadata()):
        assert metadata["source"] and metadata["release_date"] and metadata["query"]


EIGHT_CASES = (
    "clustered", "no_cluster", "sparse", "conflict", "low_plddt",
    "no_significant_hotspot", "excessive_coverage", "permutation_resolution",
)


def test_all_eight_cases_are_exposed():
    assert CASE_NAMES == EIGHT_CASES
    assert set(all_synthetic_cases()) == set(EIGHT_CASES)


def test_the_required_cases_exist_with_their_advertised_properties():
    built = all_synthetic_cases()
    assert set(built) == set(CASE_NAMES)

    assert built["clustered"].expected["has_spatial_cluster"] is True
    assert built["no_cluster"].expected["has_spatial_cluster"] is False
    assert built["sparse"].expected["N"] <= 6
    assert built["conflict"].expected["n_conflict"] > 0
    assert built["low_plddt"].expected["has_low_confidence_tail"] is True
    assert built["low_plddt"].expected["n_low_confidence_residues"] > 20


def test_every_case_declares_what_it_is_for():
    """The integration test asserts intent, so intent must be recorded."""
    for name, case in all_synthetic_cases().items():
        expected = case.expected
        assert expected["purpose"], name
        assert expected["exercises"], name
        assert isinstance(expected["reaches_stage_cd"], bool), name
        # no expected p-values, radii or hotspot counts: those are the pipeline's job
        forbidden = {"r_hot", "r_fp", "p_value", "p_values", "n_significant",
                     "footprint_radius", "hotspot_radius", "q_value"}
        assert not forbidden & set(expected), name


@pytest.mark.parametrize("name,branch", [
    ("no_significant_hotspot", "S_empty"),
    ("excessive_coverage", "all_footprint_candidates_excessive_coverage"),
    ("no_cluster", "no_global_clustering_or_S_empty"),
    ("sparse", "underpowered_cohort"),
])
def test_negative_branch_cases_declare_their_branch(name, branch):
    assert make_synthetic_case(name).expected["expected_negative_branch"] == branch


def test_permutation_resolution_case_makes_the_floor_reachable():
    """m must exceed q * (B + 1) or the k=1 BH threshold cannot sit below 1/(B+1)."""
    from hotspot3d.utils.config import load_config

    config = load_config(REPO_ROOT / "config" / "pipeline.yaml")
    q = float(config.get("fdr.q"))
    B = int(config.get("permutation.B_default"))

    expected = make_synthetic_case("permutation_resolution").expected
    m = expected["M"]
    assert m > q * (B + 1), "U_struct too small to drive the BH threshold below 1/(B+1)"
    assert q / m < 1.0 / (B + 1)
    assert expected["N_P"] <= 6, "the pinned-p condition needs very few P/LP residues"


def test_excessive_coverage_case_is_small_and_broadly_occupied():
    expected = make_synthetic_case("excessive_coverage").expected
    assert expected["M"] < 100, "QC-F1 needs a small U_struct denominator"
    assert expected["N_P"] / expected["M"] > 0.25, "P/LP must occupy much of the protein"
    assert expected["n_plp_sequence_segments"] >= 2


def test_dispersed_case_is_more_dispersed_than_the_random_one():
    """`no_significant_hotspot` must be genuinely non-clustered, not merely small."""
    import numpy as np

    from hotspot3d.utils.geometry import pairwise_distances

    def mean_nearest_neighbour(name):
        case = make_synthetic_case(name)
        by_index = case.structure.by_index()
        coords = np.array([by_index[i].ca.xyz for i in case.expected["plp_residues"]])
        distances = pairwise_distances(coords)
        np.fill_diagonal(distances, np.inf)
        return distances.min(axis=1).mean()

    assert (mean_nearest_neighbour("no_significant_hotspot")
            > mean_nearest_neighbour("no_cluster")
            > mean_nearest_neighbour("clustered"))


def test_records_are_clinvar_variant_summary_shaped():
    record = make_synthetic_case("clustered").records[0]
    for column in ("VariationID", "Type", "Name", "GeneSymbol",
                   "ClinicalSignificance", "ReviewStatus", "Assembly"):
        assert column in record
    assert record["Name"].startswith("NM_")


def test_every_review_star_level_including_zero_is_represented():
    from hotspot3d.data.review_status import star_level

    levels = {star_level(r["ReviewStatus"]) for r in make_synthetic_case("clustered").records}
    assert {0, 1, 2, 3, 4} <= levels
    assert None in levels, "an unrecognised wording must be present too"


def test_fixtures_import_without_a_circular_import_in_either_order():
    for statement in ("import hotspot3d.structure, hotspot3d.data",
                      "import hotspot3d.data, hotspot3d.structure",
                      "from hotspot3d.data import run_stage_a",
                      "from hotspot3d.structure import map_cohort"):
        result = subprocess.run([sys.executable, "-c", statement],
                                capture_output=True, text=True, cwd=REPO_ROOT)
        assert result.returncode == 0, f"{statement}: {result.stderr}"
