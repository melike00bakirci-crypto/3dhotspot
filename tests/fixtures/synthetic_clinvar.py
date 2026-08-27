"""Synthetic Stage A fixtures — the path named in the build plan.

The implementation lives in :mod:`hotspot3d.data.synthetic` so that it is
importable from anywhere the package is importable, without depending on how
pytest happens to arrange ``sys.path``. This module re-exports it verbatim;
both import paths return the same objects.

    from tests.fixtures.synthetic_clinvar import make_synthetic_case
    from hotspot3d.data.synthetic import make_synthetic_case      # identical

    case = make_synthetic_case("clustered")
    case.records     # list[dict]  — ClinVar variant_summary-shaped rows
    case.structure   # StructureModel — deterministic 3.8 A CA scaffold
    case.expected    # dict — ground truth (N, N_P, N_B, n_conflict, ...)
    case.source      # a StageASource ready for run_stage_a(..., source=...)

Eight cases::

    Stage A behaviour        clustered, no_cluster, sparse, conflict, low_plddt
    downstream branches      no_significant_hotspot, excessive_coverage,
                             permutation_resolution

``case.expected`` records what each case is FOR — ``purpose``, ``exercises``,
``expected_negative_branch``, ``reaches_stage_cd`` — so an integration test can
assert intent rather than incidental numbers. It deliberately contains no
expected p-values, radii or hotspot counts: those are the pipeline's to compute.
"""
from __future__ import annotations

from hotspot3d.data.synthetic import (  # noqa: F401  (re-export)
    ALT_TRANSCRIPT,
    BLB_LABELS,
    BUNDLE_COLS,
    BUNDLE_PITCH_A,
    CA_SPACING_A,
    CASE_NAMES,
    CASE_SPECS,
    CLINVAR_RELEASE,
    GENE,
    HELIX_LEN,
    MANE_TRANSCRIPT,
    MIN_MST_EDGE_FOR_FOOTPRINT_A,
    PLP_LABELS,
    REVIEW_STATUSES,
    UNIPROT_ACC,
    FoldLayout,
    SyntheticCase,
    SyntheticStageASource,
    all_synthetic_cases,
    backbone_coordinates,
    build_structure,
    consecutive_ca_distances,
    fold_layout,
    make_synthetic_case,
    make_synthetic_source,
    max_mst_edge,
    sequence_segments,
    synthetic_sequence,
)

__all__ = [
    "make_synthetic_case", "make_synthetic_source", "all_synthetic_cases",
    "SyntheticCase", "SyntheticStageASource", "CASE_NAMES", "CASE_SPECS",
    "synthetic_sequence", "backbone_coordinates", "build_structure",
    "consecutive_ca_distances", "fold_layout", "FoldLayout", "max_mst_edge",
    "sequence_segments", "GENE", "UNIPROT_ACC", "MANE_TRANSCRIPT",
    "ALT_TRANSCRIPT", "CLINVAR_RELEASE", "PLP_LABELS", "BLB_LABELS",
    "REVIEW_STATUSES", "CA_SPACING_A", "HELIX_LEN", "BUNDLE_PITCH_A",
    "BUNDLE_COLS", "MIN_MST_EDGE_FOR_FOOTPRINT_A",
]
