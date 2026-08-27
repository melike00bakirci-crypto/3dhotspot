"""Structure QC, pLDDT recording, U_struct/L construction and mmCIF round-tripping."""
from __future__ import annotations

import math

import pytest

from hotspot3d.data.synthetic import (
    CA_SPACING_A,
    CASE_NAMES,
    CASE_SPECS,
    MIN_MST_EDGE_FOR_FOOTPRINT_A,
    consecutive_ca_distances,
    make_synthetic_case,
)
from hotspot3d.structure.cif import parse_mmcif, render_mmcif, split_cif_row
from hotspot3d.structure.mapping import (
    CLASSIFIED_COHORT_COLUMNS,
    POSITIONAL_UNIVERSE_COLUMNS,
    RESIDUE_COORDINATE_COLUMNS,
    map_cohort,
    positional_universe,
    residue_coordinates,
)
from hotspot3d.structure.plddt import (
    band_for,
    band_labels,
    low_confidence_regions,
    plddt_profile,
    summarize,
)
from hotspot3d.structure.qc import align_numbering, assert_numbering_usable, structure_qc
from hotspot3d.structure.sources import assert_model_admissible
from hotspot3d.utils.errors import BlockedError, EscalationRequired

pytestmark = pytest.mark.unit

BANDS = [50.0, 70.0, 90.0]


# --- the synthetic scaffold itself -----------------------------------------

def test_backbone_spacing_is_exactly_38_angstrom(cases):
    for name, case in cases.items():
        distances = consecutive_ca_distances(case.structure)
        assert len(distances) == len(case.sequence) - 1
        for d in distances:
            assert math.isclose(d, CA_SPACING_A, rel_tol=1e-9), name


def test_fixtures_are_deterministic():
    a, b = make_synthetic_case("clustered"), make_synthetic_case("clustered")
    assert a.records == b.records
    assert a.expected == b.expected
    assert [r.ca.xyz for r in a.structure.residues] == [r.ca.xyz for r in b.structure.residues]


def test_clustered_case_is_spatially_tighter_than_the_random_one(cases):
    """The fixture must actually contain the cluster it claims to contain.

    `clustered` and `no_cluster` share a scaffold, cohort sizes and record
    structure, so this is a controlled comparison rather than two unrelated
    proteins.
    """
    import numpy as np

    from hotspot3d.utils.geometry import pairwise_distances

    def mean_pairwise(case):
        by_index = case.structure.by_index()
        coords = np.array([by_index[i].ca.xyz for i in case.expected["plp_residues"]])
        d = pairwise_distances(coords)
        return d[np.triu_indices_from(d, k=1)].mean()

    assert mean_pairwise(cases["clustered"]) < mean_pairwise(cases["no_cluster"])
    assert cases["clustered"].expected["has_spatial_cluster"] is True
    assert cases["no_cluster"].expected["has_spatial_cluster"] is False
    assert (cases["clustered"].expected["sequence_length"]
            == cases["no_cluster"].expected["sequence_length"])
    assert cases["clustered"].expected["N_P"] == cases["no_cluster"].expected["N_P"]


# --- the geometric precondition every downstream stage depends on -----------

def test_a_3d_hotspot_draws_from_several_sequence_segments(cases):
    """A cluster confined to one sequence run is a 1D artefact, not a 3D hotspot."""
    for name in ("clustered", "conflict", "low_plddt", "excessive_coverage",
                 "permutation_resolution"):
        expected = cases[name].expected
        assert expected["n_plp_sequence_segments"] >= 2, name


@pytest.mark.parametrize("name", CASE_NAMES)
def test_max_mst_edge_keeps_the_footprint_domain_non_empty(cases, name):
    """METHOD_SPEC II.9: rho_max = 1.25 * (max MST edge / 2) must clear rho_min = 3.0.

    A cluster of sequence-adjacent residues has every MST edge at 3.8 A, giving
    rho_max = 2.375 and an empty domain, so Stage C terminates before it can
    exercise anything. This is the fixture-side precondition for that not to
    happen.
    """
    expected = cases[name].expected
    max_edge = expected["max_mst_edge_plp_A"]
    assert max_edge > MIN_MST_EDGE_FOR_FOOTPRINT_A, name

    rho_all = max_edge / 2.0
    rho_max = min(1.25 * rho_all, 0.25 * expected["D_max_A"], 20.0)
    assert rho_max > 3.0, f"{name}: rho_min 3.0 would exceed rho_max {rho_max:.3f}"


def test_helix_cores_pack_in_the_real_ca_contact_range(clustered):
    """Inter-helix CA-CA contacts must land where real contacts do (~5-12 A).

    Helix *cores* only: a loop is chain-bonded to the helix it runs into, so the
    junction is 3.8 A by definition and says nothing about packing.
    """
    import numpy as np

    from hotspot3d.data.synthetic import fold_layout
    from hotspot3d.utils.geometry import cross_distances

    layout = fold_layout(clustered.expected["sequence_length"])
    by_index = clustered.structure.by_index()

    def core(helix):
        first, last = helix
        return np.array([by_index[i].ca.xyz for i in range(first, last + 1)])

    closest = cross_distances(core(layout.helices[0]), core(layout.helices[1])).min()
    assert MIN_MST_EDGE_FOR_FOOTPRINT_A < closest < 12.0, closest


def test_center_neighbourhoods_keep_a_long_mst_edge(cases):
    """The invariant Stage C actually depends on, checked without running Stage B.

    Stage A cannot compute the significant-center set — that is Stage B's output
    and reading it would breach the information barrier. What can be checked is
    that every *plausible* center set stays multi-segment: for a range of sphere
    radii and enrichment thresholds, the max MST edge over residues whose
    neighbourhood contains P/LP stays above 4.8 A. Denser center sets have
    shorter MST edges, so sweeping the thresholds is the pessimistic check.
    """
    import numpy as np

    from hotspot3d.data.synthetic import backbone_coordinates, max_mst_edge
    from hotspot3d.utils.geometry import cross_distances

    for name in ("clustered", "conflict", "low_plddt", "excessive_coverage",
                 "permutation_resolution"):
        case = cases[name]
        spec = CASE_SPECS[name]
        coords = backbone_coordinates(spec.length, spec.seed, n_cols=spec.n_cols,
                                      pitch=spec.pitch)
        plp = np.asarray(case.expected["plp_residues"]) - 1
        distances = cross_distances(coords, coords[plp])

        for radius in (5.0, 6.0, 7.0, 8.0):
            n_near = (distances <= radius).sum(axis=1)
            for threshold in (1, 2):
                centers = [i + 1 for i in range(len(coords)) if n_near[i] >= threshold]
                if len(centers) < 2:
                    continue
                edge = max_mst_edge(coords, centers)
                assert edge > MIN_MST_EDGE_FOR_FOOTPRINT_A, (
                    f"{name}: centers within {radius} A of >={threshold} P/LP give "
                    f"max MST edge {edge:.2f} A, so rho_max would fall below rho_min")


# --- structural QC ----------------------------------------------------------

def test_qc_passes_on_a_well_formed_model(clustered):
    result = structure_qc(clustered.structure, canonical_sequence=clustered.sequence)
    assert result.passed, [c.as_dict() for c in result.failures]
    assert result.as_dict()["verdict"] == "PASS"
    assert result.payload["chains"] == ["A"]
    assert result.payload["n_residues_modelled"] == len(clustered.sequence)
    assert result.payload["structure_edited"] is False
    assert result.payload["numbering_repaired"] is False


def test_every_check_is_recorded_pass_or_fail(clustered):
    result = structure_qc(clustered.structure, canonical_sequence=clustered.sequence)
    ids = [c.check_id for c in result.checks]
    assert len(ids) == len(set(ids)) >= 12
    assert all("passed" in c.as_dict() for c in result.checks)


def test_numbering_offset_is_reported_and_never_applied(clustered):
    result = structure_qc(clustered.structure, canonical_sequence=clustered.sequence)
    assert result.alignment.offset == 0
    assert result.alignment.identity == 1.0
    assert result.alignment.as_dict()["offset_applied"] is False
    assert result.alignment.as_dict()["policy"] == "reported, never silently applied"
    assert assert_numbering_usable(result) is None


def test_a_shifted_sequence_escalates_rather_than_being_repaired(clustered):
    shifted = "GGGGG" + clustered.sequence
    result = structure_qc(clustered.structure, canonical_sequence=shifted)
    assert result.alignment.offset == 5
    with pytest.raises(EscalationRequired) as excinfo:
        assert_numbering_usable(result)
    assert "never applied" in excinfo.value.recommendation


def test_alignment_offset_direction():
    alignment = align_numbering("ACDEF", "XXACDEF")
    assert alignment.offset == 2
    assert alignment.longest_block == 5


def test_multi_fragment_entry_is_blocked_not_stitched(clustered):
    case = clustered.with_fragments(3)
    with pytest.raises(BlockedError, match="MULTI_FRAGMENT_AFDB_ENTRY"):
        assert_model_admissible(case.structure)


def test_undeclared_fragment_count_escalates(clustered):
    from dataclasses import replace

    model = replace(clustered.structure, n_fragments=None)
    with pytest.raises(EscalationRequired):
        assert_model_admissible(model)


def test_missing_entry_escalates():
    with pytest.raises(EscalationRequired):
        assert_model_admissible(None)


# --- pLDDT (F2) -------------------------------------------------------------

def test_band_labels_and_assignment():
    assert band_labels(BANDS) == ["<50", "50-70", "70-90", ">90"]
    assert band_for(10.0, BANDS) == "<50"
    assert band_for(49.999, BANDS) == "<50"
    assert band_for(50.0, BANDS) == "50-70"
    assert band_for(70.0, BANDS) == "70-90"
    assert band_for(90.0, BANDS) == ">90"
    assert band_for(None, BANDS) == "NA"


def test_plddt_is_recorded_for_every_modelled_residue(clustered):
    profile = plddt_profile(clustered.structure, BANDS)
    assert len(profile) == len(clustered.sequence)
    assert all(row["plddt"] is not None for row in profile)
    assert summarize(profile, BANDS).n_with_plddt == len(profile)


def test_low_confidence_regions_are_recorded_and_never_exclude(cases):
    case = cases["low_plddt"]
    profile = plddt_profile(case.structure, BANDS)
    regions = low_confidence_regions(profile, [50.0, 70.0])

    assert regions, "the low_plddt fixture must produce low-confidence regions"
    assert all(r["excluded"] is False for r in regions)

    covered = {i for r in regions if r["threshold"] == 70.0
               for i in range(r["start_residue"], r["end_residue"] + 1)}
    assert covered == set(case.expected["low_confidence_residues"])


def test_the_low_confidence_tail_does_not_shrink_the_cohort(cases, make_ctx):
    from hotspot3d.data.stage import run_stage_a

    case = cases["low_plddt"]
    handoff = run_stage_a(make_ctx(), gene=case.gene, source=case.source)
    assert handoff.payload["N"] == case.expected["N"]
    assert handoff.payload["N_P"] == case.expected["N_P"]
    assert handoff.payload["N_B"] == case.expected["N_B"]
    assert handoff.payload["M"] == case.expected["M"], "U_struct is not pLDDT-filtered"


def test_low_confidence_regions_are_maximal_runs():
    profile = [{"residue_index": i, "plddt": p, "plddt_band": "NA", "aa": "A",
                "chain_id": "A", "ca_usable": True}
               for i, p in enumerate([95, 40, 42, 41, 95, 30, 95], start=1)]
    regions = low_confidence_regions(profile, [50.0])
    assert [(r["start_residue"], r["end_residue"]) for r in regions] == [(2, 4), (6, 6)]
    assert [r["length"] for r in regions] == [3, 1]


# --- U_struct and L ---------------------------------------------------------

def test_universe_and_cohort_column_sets_are_frozen(clustered):
    rows = residue_coordinates(clustered.structure, BANDS)
    assert set(rows[0]) == set(RESIDUE_COORDINATE_COLUMNS)
    assert all(row["unused_under_F1"] is True for row in rows)

    universe = positional_universe(rows)
    assert set(universe[0]) == set(POSITIONAL_UNIVERSE_COLUMNS)


def test_cohort_columns_are_exactly_the_stage_b_allowlist():
    from hotspot3d.orchestration.contracts import (
        FORBIDDEN_DOWNSTREAM_COLUMNS,
        STAGE_B_PERMITTED_COLUMNS,
    )

    assert set(CLASSIFIED_COHORT_COLUMNS) == set(STAGE_B_PERMITTED_COLUMNS)
    assert not set(CLASSIFIED_COHORT_COLUMNS) & set(FORBIDDEN_DOWNSTREAM_COLUMNS)


def test_glycine_has_no_cb_and_reports_na(clustered):
    rows = {r["residue_index"]: r for r in residue_coordinates(clustered.structure, BANDS)}
    glycines = [i for i, r in rows.items() if r["aa"] == "G"]
    assert glycines, "the fixture must contain glycine"
    for index in glycines:
        assert rows[index]["x_cb"] is None
        assert rows[index]["x_sc"] is None
        assert rows[index]["ca_usable"] is True


def test_l_is_a_subset_of_u_struct(cases, config):
    from hotspot3d.data.classify import SignificancePolicy
    from hotspot3d.data.clinvar import normalize_records
    from hotspot3d.data.cohort import collapse_to_residues, evaluate_records

    policy = SignificancePolicy.from_config(config)
    for name, case in cases.items():
        evals = evaluate_records(normalize_records(case.records, gene=case.gene),
                                 policy=policy, canonical_sequence=case.sequence,
                                 mane_transcript=case.mane_transcript)
        rows = collapse_to_residues(evals, canonical_sequence=case.sequence)
        coords = residue_coordinates(case.structure, BANDS)
        mapping = map_cohort(rows, coords)

        universe = {r["residue_index"] for r in mapping.universe}
        assert {r["residue_index"] for r in mapping.cohort} <= universe, name
        assert mapping.report["cohort_subset_of_universe"] is True
        assert mapping.report["position_correspondence_verified"] is True
        assert mapping.M == case.expected["M"], name
        assert mapping.N == case.expected["N"], name


def test_a_residue_without_a_usable_ca_leaves_u_struct_and_is_listed(clustered):
    victim = clustered.expected["plp_residues"][0]
    case = clustered.with_unusable_ca(victim)
    coords = residue_coordinates(case.structure, BANDS)
    universe = positional_universe(coords)

    assert victim not in {r["residue_index"] for r in universe}
    assert len(universe) == len(case.sequence) - 1
    row = next(r for r in coords if r["residue_index"] == victim)
    assert row["ca_usable"] is False


# --- mmCIF ------------------------------------------------------------------

def test_cif_round_trip_preserves_coordinates_and_plddt(clustered):
    text = render_mmcif(clustered.structure)
    parsed = parse_mmcif(text, accession=clustered.uniprot_acc)

    assert parsed.n_residues == clustered.structure.n_residues
    original = clustered.structure.by_index()
    for index, residue in parsed.by_index().items():
        assert residue.aa == original[index].aa
        assert residue.plddt == pytest.approx(original[index].plddt, abs=0.01)
        for axis, value in zip("xyz", residue.ca.xyz):
            assert value == pytest.approx(getattr(original[index].ca, axis), abs=0.001)


def test_cif_row_tokenizer_handles_quoting():
    assert split_cif_row("ATOM 1 C 'CA 1' . ALA") == ["ATOM", "1", "C", "CA 1", ".", "ALA"]


def test_payload_without_atom_site_is_blocked():
    with pytest.raises(BlockedError, match="STRUCTURE_MAPPING_FAILURE"):
        parse_mmcif("data_X\n#\n_entry.id X\n#\n")
