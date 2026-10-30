"""METHOD_SPEC II.7 / Workflow v2 §5.1-§5.3 — final detection at the selected ``r_hot``.

**Every** residue of ``U_center`` is a candidate sphere center — not only residues
that carry a ClinVar variant (F9). A significant hotspot center is a *geometric test
position*; the number of significant centers carrying no ClinVar variant at all is
counted and published, because conflating "center" with "variant" is the single
easiest way to misread this pipeline.

**[v2 §2] U_center vs U_struct.** ``U_center`` is a candidate-center-eligibility
restriction of ``U_struct`` (residues with ``pLDDT >= plddt.center_universe_min_plddt``
are coordinate noise as sphere centers and are excluded). It does **not** revise F9's
core point — that ANY residue, not only cohort residues, may be a center — it only
narrows which residues of ``U_struct`` are tested. The positional-null PLACEMENT
population remains the full, unrestricted ``U_struct`` in every case (F9/§3
unchanged); only the row space of the membership matrix — which residues get a
center-test row at all — is restricted.

**[v2 §5.2] The primary per-center test is the structure-aware positional null** —
the same null as the global analysis of §3. The label-permutation null is computed
and published **beside** it as a SECONDARY analysis of a different hypothesis:

  positional : are P/LP residues placed more tightly than uniform placement over
               ``U_struct``, conditioned on ``N_P``?
  label      : are P/LP residues more clustered than B/LB residues, GIVEN the set of
               variant-bearing positions?

The two are never mixed. ``significant``, ``q_bh``, the three residue objects and
the hotspot regions all come from the primary null alone; the secondary columns
carry the ``_label_null`` suffix and decide nothing.

Three objects are emitted, **separately named and never conflated**, with
deliberately DISTINCT key column names so that an accidental join fails loudly:

  ``SIGNIFICANT_HOTSPOT_CENTERS``          key ``center_residue_index``
  ``HOTSPOT_SPHERE_CLASSIFIED_VARIANTS``   key ``variant_residue_index``
  ``HOTSPOT_COVERED_RESIDUES``             key ``covered_residue_index``

Their nesting (``centers subset-of covered``; ``classified subset-of covered & L``)
is verified in code, not assumed.

Hotspot regions are the connected components of ``S`` under ``d(c,c') <= 2*r_hot``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..utils.geometry import (
    covered_mask,
    cross_distances,
    membership_matrix,
    neighbours_within,
)
from ..utils.geometry import connected_components
from .fdr import FDRResult, benjamini_hochberg
from .permutation import PermutationOutcome, permutation_pass
from .positional import PositionalOutcome, positional_pass

ALL_CENTER_TEST_COLUMNS = [
    "schema_version", "center_residue_index", "x_ca", "y_ca", "z_ca", "radius_A",
    "n_labeled_in_sphere", "n_plp_in_sphere", "n_blb_in_sphere", "expected_plp",
    "null_mean_plp", "null_sd_plp", "z_score", "p_emp", "q_bh", "q_by",
    "significant", "significant_by", "in_test_family", "exclusion_reason",
    "carries_clinvar_variant", "hotspot_id", "B", "q", "fdr_method",
    # --- v2 §5.2: which null produced the decision, and its exact tail --------
    "null_model", "n_universe_in_sphere", "p_exact_hypergeom", "p_comb_floor",
    # --- v2 §5.2: the SECONDARY label-permutation null, reported side by side --
    "p_emp_label_null", "q_bh_label_null", "significant_label_null",
    "in_test_family_label_null", "expected_plp_label_null",
]

CLASSIFIED_VARIANT_COLUMNS = [
    "schema_version", "variant_residue_index", "class", "hotspot_ids",
    "n_covering_centers", "x_ca", "y_ca", "z_ca",
]

COVERED_RESIDUE_COLUMNS = [
    "schema_version", "covered_residue_index", "hotspot_ids", "is_classified",
    "class", "is_significant_center", "n_covering_centers", "x_ca", "y_ca", "z_ca",
]

HOTSPOT_REGION_COLUMNS = [
    "schema_version", "hotspot_id", "n_centers", "center_residue_indices",
    "residue_index_min", "residue_index_max", "n_covered_residues",
    "n_classified_residues", "n_plp", "n_blb", "fold_enrichment", "min_p_emp",
    "min_q_bh", "centroid_x", "centroid_y", "centroid_z",
    "max_center_pairwise_distance_A", "radius_A",
]

BH_TABLE_COLUMNS = [
    "schema_version", "rank", "center_residue_index", "p_emp", "bh_critical_value",
    "q_bh", "q_by", "reject_bh", "reject_by", "is_boundary", "m", "q",
]


@dataclass
class DetectionResult:
    radius_A: float
    perm: PositionalOutcome              # PRIMARY null (v2 §5.2)
    fdr: FDRResult
    family_positions: np.ndarray
    significant_positions: np.ndarray
    center_rows: list[dict]
    significant_rows: list[dict]
    classified_rows: list[dict]
    covered_rows: list[dict]
    region_rows: list[dict]
    bh_rows: list[dict]
    hotspot_of_center: dict[int, str]
    n_centers_without_variant: int
    boundary_p: float | None
    counts: dict = field(default_factory=dict)
    #: SECONDARY label-permutation analysis (v2 §5.2). ``None`` when not computed
    #: (the sensitivity re-runs do not need it). It NEVER feeds primary inference.
    label_perm: PermutationOutcome | None = None
    label_fdr: FDRResult | None = None
    secondary: dict = field(default_factory=dict)

    @property
    def n_significant(self) -> int:
        return len(self.significant_positions)

    @property
    def positional(self) -> PositionalOutcome:
        """Explicit alias — ``perm`` is the PRIMARY positional outcome (v2 §5.2)."""
        return self.perm


def _hotspot_ids(covering: np.ndarray, sig_positions: np.ndarray,
                 hotspot_of_position: dict[int, str]) -> str:
    ids = sorted({hotspot_of_position[int(sig_positions[j])]
                  for j in np.nonzero(covering)[0]})
    return ";".join(ids) if ids else "NA"


def final_detection(radius: float, universe_coords: np.ndarray, universe_index: np.ndarray,
                    labeled_positions: np.ndarray, labeled_coords: np.ndarray,
                    y: np.ndarray, rng: np.random.Generator, *, B: int, q: float,
                    seed_context: str, chunk: int = 2000,
                    label_rng: np.random.Generator | None = None,
                    label_seed_context: str | None = None,
                    center_coords: np.ndarray | None = None,
                    center_index: np.ndarray | None = None) -> DetectionResult:
    """The complete II.7 procedure at the frozen ``r_hot``, under the PRIMARY null.

    ``label_rng`` enables the SECONDARY label-permutation analysis of v2 §5.2. It is
    drawn from its own seed context so the two nulls never share a random stream,
    and it never influences ``significant``, ``q_bh`` or any published object.

    ``center_coords``/``center_index`` restrict the CANDIDATE sphere centers to v2
    §2's ``U_center`` (a subset of ``universe_coords``/``universe_index``); the
    universe arguments remain the FULL positional reference ``U_struct`` — the
    placement population and the coverage/region targets are unaffected. Defaults
    to the universe, preserving pre-v2 (centers == U_struct, F9) behaviour for
    direct library callers (e.g. the sensitivity re-runs, which restrict the
    universe wholesale and need no separate center restriction).
    """
    universe_coords = np.asarray(universe_coords, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    M = len(universe_coords)
    N = len(y)
    pi = float(y.sum()) / N if N else 0.0

    if center_coords is None:
        center_coords, center_index = universe_coords, universe_index
    center_coords = np.asarray(center_coords, dtype=np.float64)
    n_centers = len(center_coords)

    A = membership_matrix(center_coords, labeled_coords, radius)
    U = (cross_distances(center_coords, universe_coords) <= radius).astype(np.uint8)
    plp_positions = np.asarray(labeled_positions, dtype=np.int64)[y > 0]
    perm = positional_pass(U, plp_positions, A.sum(axis=1), rng, B, radius, seed_context)

    family = np.nonzero(perm.in_family)[0]
    fdr = benjamini_hochberg(perm.p_emp[family], q)

    q_bh = np.full(n_centers, np.nan)
    q_by = np.full(n_centers, np.nan)
    q_bh[family] = fdr.q_bh
    q_by[family] = fdr.q_by
    significant = np.zeros(n_centers, dtype=bool)
    significant_by = np.zeros(n_centers, dtype=bool)
    significant[family] = fdr.reject_bh
    significant_by[family] = fdr.reject_by
    sig_positions = np.nonzero(significant)[0]

    # --- SECONDARY: label permutation, reported beside the primary (v2 §5.2) ---
    label_perm: PermutationOutcome | None = None
    label_fdr: FDRResult | None = None
    label_q_bh = np.full(n_centers, np.nan)
    label_significant = np.zeros(n_centers, dtype=bool)
    if label_rng is not None:
        label_perm = permutation_pass(
            A, y, label_rng, B, radius,
            label_seed_context or f"{seed_context}|secondary_label_null", chunk=chunk)
        label_family = np.nonzero(label_perm.in_family)[0]
        label_fdr = benjamini_hochberg(label_perm.p_emp[label_family], q)
        label_q_bh[label_family] = label_fdr.q_bh
        label_significant[label_family] = label_fdr.reject_bh

    # --- hotspot regions: connected components of S under d(c,c') <= 2 r ------
    # ``regions``/``sig_positions`` live in CENTER space; residue-index values are
    # used as the join key with UNIVERSE space wherever the two meet below.
    hotspot_of_position: dict[int, str] = {}
    regions: list[list[int]] = []
    if len(sig_positions):
        s_coords = center_coords[sig_positions]
        adjacency = neighbours_within(s_coords, 2.0 * radius, exclude_self=True)
        for comp in connected_components(adjacency):
            regions.append([int(sig_positions[i]) for i in comp])
        regions.sort(key=lambda comp: int(center_index[comp[0]]))
        for n, comp in enumerate(regions, start=1):
            for pos in comp:
                hotspot_of_position[pos] = f"H{n}"

    # Residue-index sets: the canonical join key between CENTER space and UNIVERSE
    # space, which diverge once U_center is a strict subset of U_struct (v2 §2).
    # ``labeled_positions`` are positions WITHIN universe_coords/universe_index.
    labeled_residue_index_set = {int(universe_index[p]) for p in labeled_positions}
    significant_residue_index_set = {int(center_index[p]) for p in sig_positions}
    class_of_position = {int(p): ("PLP" if y[i] == 1 else "BLB")
                         for i, p in enumerate(labeled_positions)}

    z = perm.z_center
    center_rows: list[dict] = []
    for c in range(n_centers):
        in_family = bool(perm.in_family[c])
        residue_idx = int(center_index[c])
        center_rows.append({
            "center_residue_index": residue_idx,
            "x_ca": float(center_coords[c, 0]), "y_ca": float(center_coords[c, 1]),
            "z_ca": float(center_coords[c, 2]), "radius_A": float(radius),
            "n_labeled_in_sphere": int(perm.n_labeled[c]),
            "n_plp_in_sphere": int(perm.n_plp_obs[c]),
            "n_blb_in_sphere": int(perm.n_labeled[c] - perm.n_plp_obs[c]),
            # Expectation under the PRIMARY (positional) null: N_P * n_U(c) / M.
            "expected_plp": float(perm.expected_plp[c]),
            "null_mean_plp": float(perm.null_mean[c]),
            "null_sd_plp": float(perm.null_sd[c]),
            "z_score": float(z[c]),
            "p_emp": float(perm.p_emp[c]),
            "q_bh": (float(q_bh[c]) if in_family else None),
            "q_by": (float(q_by[c]) if in_family else None),
            "significant": bool(significant[c]),
            "significant_by": bool(significant_by[c]),
            "in_test_family": in_family,
            "exclusion_reason": perm.exclusion_reason[c],
            "carries_clinvar_variant": bool(residue_idx in labeled_residue_index_set),
            "hotspot_id": hotspot_of_position.get(c, "NA"),
            "B": int(B), "q": float(q), "fdr_method": "BH",
            "null_model": perm.null_model,
            "n_universe_in_sphere": int(perm.n_universe[c]),
            "p_exact_hypergeom": float(perm.p_exact[c]),
            "p_comb_floor": float(perm.p_comb_floor[c]),
            "p_emp_label_null": (float(label_perm.p_emp[c]) if label_perm else None),
            "q_bh_label_null": (float(label_q_bh[c])
                                if (label_perm is not None and
                                    bool(label_perm.in_family[c])) else None),
            "significant_label_null": (bool(label_significant[c]) if label_perm else None),
            "in_test_family_label_null": (bool(label_perm.in_family[c])
                                          if label_perm else None),
            "expected_plp_label_null": (float(pi * perm.n_labeled[c])
                                        if label_perm else None),
        })

    significant_rows = [center_rows[c] for c in sig_positions]
    n_without_variant = int(sum(1 for r in significant_rows
                                if not r["carries_clinvar_variant"]))

    # --- the two covered-residue objects (over the FULL universe U_struct) -----
    classified_rows: list[dict] = []
    covered_rows: list[dict] = []
    if len(sig_positions):
        s_coords = center_coords[sig_positions]
        cover = cross_distances(universe_coords, s_coords) <= radius     # (M, |S|)
        for c in np.nonzero(cover.any(axis=1))[0]:
            residue_idx = int(universe_index[c])
            is_classified = residue_idx in labeled_residue_index_set
            ids = _hotspot_ids(cover[c], sig_positions, hotspot_of_position)
            covered_rows.append({
                "covered_residue_index": residue_idx,
                "hotspot_ids": ids, "is_classified": is_classified,
                "class": class_of_position.get(int(c), "NA"),
                "is_significant_center": bool(residue_idx in significant_residue_index_set),
                "n_covering_centers": int(cover[c].sum()),
                "x_ca": float(universe_coords[c, 0]),
                "y_ca": float(universe_coords[c, 1]),
                "z_ca": float(universe_coords[c, 2]),
            })
            if is_classified:
                classified_rows.append({
                    "variant_residue_index": residue_idx,
                    "class": class_of_position[int(c)], "hotspot_ids": ids,
                    "n_covering_centers": int(cover[c].sum()),
                    "x_ca": float(universe_coords[c, 0]),
                    "y_ca": float(universe_coords[c, 1]),
                    "z_ca": float(universe_coords[c, 2]),
                })

    # --- region summary rows ---------------------------------------------------
    region_rows: list[dict] = []
    for n, comp in enumerate(regions, start=1):
        hid = f"H{n}"
        comp_coords = center_coords[comp]
        comp_cover = covered_mask(comp_coords, universe_coords, radius)
        comp_labeled = covered_mask(comp_coords, np.asarray(labeled_coords), radius)
        n_lab = int(comp_labeled.sum())
        n_plp = int(y[comp_labeled].sum())
        d = cross_distances(comp_coords, comp_coords)
        region_rows.append({
            "hotspot_id": hid, "n_centers": len(comp),
            "center_residue_indices": ";".join(str(int(center_index[p])) for p in comp),
            "residue_index_min": int(min(center_index[p] for p in comp)),
            "residue_index_max": int(max(center_index[p] for p in comp)),
            "n_covered_residues": int(comp_cover.sum()),
            "n_classified_residues": n_lab,
            "n_plp": n_plp, "n_blb": n_lab - n_plp,
            "fold_enrichment": ((n_plp / n_lab) / pi if (n_lab > 0 and pi > 0) else None),
            "min_p_emp": float(min(perm.p_emp[p] for p in comp)),
            "min_q_bh": float(min(q_bh[p] for p in comp)),
            "centroid_x": float(comp_coords[:, 0].mean()),
            "centroid_y": float(comp_coords[:, 1].mean()),
            "centroid_z": float(comp_coords[:, 2].mean()),
            "max_center_pairwise_distance_A": float(d.max()),
            "radius_A": float(radius),
        })

    # --- the BH table ----------------------------------------------------------
    bh_rows: list[dict] = []
    if len(family):
        order = np.argsort(fdr.p, kind="stable")
        m = len(family)
        for rank, j in enumerate(order, start=1):
            c = int(family[j])
            bh_rows.append({
                "rank": rank,
                "center_residue_index": int(center_index[c]),
                "p_emp": float(fdr.p[j]),
                "bh_critical_value": float(rank * q / m),
                "q_bh": float(fdr.q_bh[j]), "q_by": float(fdr.q_by[j]),
                "reject_bh": bool(fdr.reject_bh[j]), "reject_by": bool(fdr.reject_by[j]),
                "is_boundary": bool(fdr.boundary_p is not None and
                                    fdr.p[j] == fdr.boundary_p),
                "m": m, "q": float(q),
            })

    result = DetectionResult(
        radius_A=float(radius), perm=perm, fdr=fdr, family_positions=family,
        significant_positions=sig_positions, center_rows=center_rows,
        significant_rows=significant_rows, classified_rows=classified_rows,
        covered_rows=covered_rows, region_rows=region_rows, bh_rows=bh_rows,
        hotspot_of_center={int(center_index[p]): hid
                           for p, hid in hotspot_of_position.items()},
        n_centers_without_variant=n_without_variant, boundary_p=fdr.boundary_p,
        label_perm=label_perm, label_fdr=label_fdr,
        secondary=_secondary_summary(label_perm, label_fdr, center_index,
                                     label_significant, q),
        counts={
            "M_universe": M, "n_center_universe": n_centers, "N_labeled": N,
            "N_P": int(y.sum()), "N_B": int(N - y.sum()), "pi": pi,
            "primary_null": perm.null_model,
            "secondary_null": ("label_permutation" if label_perm is not None
                               else "not_computed"),
            "n_in_test_family": int(len(family)),
            "n_excluded_from_family": int(n_centers - len(family)),
            "n_excluded_no_labeled": int(sum(
                1 for r in perm.exclusion_reason if r == "no_labeled_residue_in_sphere")),
            "n_excluded_zero_variance": int(sum(
                1 for r in perm.exclusion_reason if r == "zero_permutation_variance")),
            "n_significant_bh": int(len(sig_positions)),
            "n_significant_by": int(significant_by.sum()),
            "n_regions": len(regions),
            "n_covered_residues": len(covered_rows),
            "n_classified_in_hotspots": len(classified_rows),
            "n_centers_without_variant": n_without_variant,
        },
    )
    verify_object_nesting(result)
    return result


def _secondary_summary(label_perm: PermutationOutcome | None, label_fdr: FDRResult | None,
                       center_index: np.ndarray, label_significant: np.ndarray,
                       q: float) -> dict:
    """v2 §5.2 — the secondary null reported side by side, never mixed in.

    Emitted with an explicit statement of the hypothesis it addresses, so no reader
    can mistake it for the primary result or aggregate the two.
    """
    if label_perm is None or label_fdr is None:
        return {
            "computed": False,
            "reason": "the secondary null is not computed for restricted re-runs",
        }
    sig = np.nonzero(label_significant)[0]
    return {
        "computed": True,
        "null_model": "label_permutation",
        "hypothesis": ("are P/LP residues more clustered than B/LB residues, GIVEN the "
                       "set of variant-bearing positions?"),
        "primary_hypothesis_for_contrast": (
            "are P/LP residues placed more tightly than uniform without-replacement "
            "placement over U_struct, conditioned on N_P?"),
        "m_test_family_size": int(label_perm.in_family.sum()),
        "n_significant_bh": int(len(sig)),
        "significant_center_residue_indices": sorted(int(center_index[c]) for c in sig),
        "smallest_p_emp": (float(label_perm.p_emp[label_perm.in_family].min())
                           if label_perm.in_family.any() else None),
        "bh_boundary_p": label_fdr.boundary_p,
        "q": float(q),
        "IS_SECONDARY_NEVER_PRIMARY": True,
        "note": ("v2 §5.2 — the two nulls answer DIFFERENT questions and are never "
                 "mixed across the stages of one run. Nothing in this block feeds "
                 "significance, the three residue objects, the hotspot regions or any "
                 "downstream stage."),
    }


def verify_object_nesting(d: DetectionResult) -> None:
    """QC rule 15 — the three objects are disjointly defined and correctly nested."""
    centers = {r["center_residue_index"] for r in d.significant_rows}
    covered = {r["covered_residue_index"] for r in d.covered_rows}
    classified = {r["variant_residue_index"] for r in d.classified_rows}
    if not centers <= covered:
        raise AssertionError(
            f"nesting violation: {len(centers - covered)} significant center(s) are not "
            f"in HOTSPOT_COVERED_RESIDUES; a center is always within r of itself.")
    if not classified <= covered:
        raise AssertionError(
            "nesting violation: HOTSPOT_SPHERE_CLASSIFIED_VARIANTS is not a subset of "
            "HOTSPOT_COVERED_RESIDUES.")
    covered_classified = {r["covered_residue_index"] for r in d.covered_rows
                          if r["is_classified"]}
    if classified != covered_classified:
        raise AssertionError(
            "nesting violation: the classified-variant object disagrees with the "
            "is_classified flag of the covered-residue object.")


def significant_subset_matches(all_rows: list[dict], sig_rows: list[dict]) -> bool:
    """QC rule 14 — the published center file is EXACTLY the significant=TRUE subset."""
    expected = [r for r in all_rows if r["significant"]]
    if len(expected) != len(sig_rows):
        return False
    return all(a == b for a, b in zip(expected, sig_rows))
