"""METHOD_SPEC II.5 / v2 §4 — the radius scan and its admissibility constraints.

Every grid radius is analysed **independently**: it gets its own seed derived from
its own context string (``scan|r=7.5``), its own membership matrices, its own null
draws and its own complete metric row. No state is carried from one radius to the
next, so radii could be evaluated in any order — or in parallel — and the result
would be bit-identical.

Four Pareto objectives, all maximizing, all kept separate:

  A  ``loo_mcc``            selection objective only (II.5A, F7/F13)
  B  ``perm_evidence_zg``   standardized evidence from T(r) (II.5B) — SECONDARY null
  C  ``fold_enrichment``    effect size, deliberately kept apart from significance (II.5C)
  D  ``neighbor_stability`` mean Jaccard of S(r) against S(r +/- step) (II.5D)

The secondary classification metrics (sens, spec, acc, ppv, npv, F1, balanced
accuracy) are computed and reported but are **not** objectives.

**[v2 §4] Significance is never an admissibility filter.** v1 treated ``|S(r)| = 0``
(QC-H1) as an admissibility failure, and QC-H2/QC-H4 are computed over ``S(r)`` too.
A radius that found nothing was therefore struck from the Pareto set, so
significance chose the radius and the radius chose significance. That circularity is
what v2 §4 forbids. Here:

  * a radius may be inadmissible ONLY when a declared quantity is genuinely
    undefined (``radius_qc.admissibility_constraints``);
  * QC-H1/QC-H2/QC-H4 are computed, recorded and reported as DIAGNOSTICS and never
    touch admissibility (``radius_qc.significance_conditioned_diagnostics``);
  * a radius with no significant center simply scores 0 on the enrichment
    objective — it stays a legitimate Pareto candidate.

**[v2 §5.2] The per-center test is the structure-aware positional null.** The
label-permutation pass is still run at every radius, because ``T(r)`` and its
``Zg`` are defined on it (II.5B), but it is a SECONDARY analysis of a different
hypothesis and never decides significance.

Inadmissible radii are retained in full — raw values, QC verdict and reason —
because a scan that silently drops radii cannot be audited.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..utils.geometry import covered_mask, cross_distances, jaccard, membership_matrix
from .fdr import benjamini_hochberg
from .loo import LooOutcome, kappa_is_boundary_neutral, leave_one_out
from .permutation import PermutationOutcome, permutation_pass
from .positional import PositionalOutcome, positional_pass

# The single decision trace (Output Contract IX, extended by v2 §4). The frozen
# prefix is unchanged; v2's per-radius admissibility record is appended.
R_HOT_SCAN_COLUMNS = [
    "schema_version", "radius_A", "n_labeled_in_universe", "loo_mcc", "loo_sens",
    "loo_spec", "loo_ppv", "loo_npv", "loo_f1", "loo_balacc",
    "loo_n_zero_neighbour", "loo_prop_zero_neighbour", "loo_n_sparse", "loo_n_tie",
    "perm_evidence_raw", "perm_evidence_normalized", "fold_enrichment",
    "fold_enrichment_normalized", "neighbor_stability", "neighbor_stability_normalized",
    "loo_mcc_normalized", "n_significant_centers", "coverage_fraction",
    "n_singleton_centers", "qc_status", "qc_failure_reason", "pareto_member",
    "distance_to_ideal", "selected", "search_domain_source", "boundary_warning",
    # --- v2 §4: admissibility recorded explicitly, rule by rule ---------------
    "admissible", "admissibility_rule", "admissibility_reason",
    "significance_diagnostics_fired",
    "qc_h1_zero_significant_centers", "qc_h2_coverage_exceeds_max",
    "qc_h3_median_n_labeled", "qc_h3_below_threshold",
    "qc_h4_isolated_center_fraction", "qc_h4_exceeds_max",
    "fold_enrichment_defined", "neighbor_stability_defined", "n_in_test_family",
    # --- v2 §5.2 / §5.4: which null decided, and could it have decided at all --
    "primary_null", "secondary_null", "power_p_res", "power_p_comb_best",
    "power_p_floor", "power_c_1", "power_certificate_passes", "power_binding_floor",
    # --- v2 §4: biological-plausibility DIAGNOSTICS — reported, NEVER Pareto
    # objectives and NEVER admissibility constraints (F15 keeps exactly four
    # objectives; see config/pipeline.yaml `biological_plausibility` for the ruling).
    "structural_coverage", "cohort_absorption", "r_to_domain_ratio",
]

SCAN_CENTER_COLUMNS = [
    "center_residue_index", "x_ca", "y_ca", "z_ca", "n_labeled_in_sphere",
    "n_universe_in_sphere", "n_plp_in_sphere", "n_blb_in_sphere", "null_mean_plp",
    "null_sd_plp", "z_score", "p_emp", "p_exact_hypergeom", "q_bh", "in_test_family",
    "exclusion_reason", "null_model",
]

#: v2 §4 — the ONLY permitted admissibility constraints: a stated quantity is
#: genuinely undefined. Each is declared in config with its rule.
ADMISSIBILITY_RULES = {
    "FE_UNDEF": ("Fold enrichment is undefined: S(r) is non-empty but its spheres "
                 "contain no classified residue, so the ratio has an empty "
                 "denominator (II.5C). NOT triggered by |S(r)| = 0 — that case scores "
                 "zero on the objective and stays admissible (v2 §4)."),
    "QC_H5_DOMAIN": ("r lies outside [R_FLOOR, r_cap] — an assertion; it must never "
                     "fire for an on-grid radius."),
    "QC_H3_UNDEFINED_MEDIAN": ("median_{c in U_struct} n_L(c) is undefined because "
                               "U_struct is empty. The v1 THRESHOLD form of QC-H3 is "
                               "demoted to a recorded diagnostic (v2 §4): thin "
                               "evidence is not an undefined quantity."),
}

#: v2 §4 — computed, recorded and reported, but NEVER admissibility. Each of these
#: is conditioned on S(r), so using it as a filter would let significance choose the
#: radius and the radius choose significance.
SIGNIFICANCE_CONDITIONED_DIAGNOSTICS = {
    "QC_H1": "|S(r)| = 0 — no center survives FDR at this radius.",
    "QC_H2": "coverage(r) > radius_qc.QC_H2_max_coverage — the spheres swallow the protein.",
    "QC_H4": ("fraction of centers in S(r) with no other center within 2r exceeds "
              "radius_qc.QC_H4_max_isolated_center_fraction — S(r) is scattered."),
}

#: Recorded, non-binding: the v1 threshold form of QC-H3.
UNCONDITIONED_DIAGNOSTICS = {
    "QC_H3_THIN_EVIDENCE": ("median_{c in U_struct} n_L(c) < "
                            "radius_qc.QC_H3_min_median_n_labeled — evidence is thin. "
                            "Reported; never an admissibility filter (v2 §4)."),
}

QC_RULES = {**ADMISSIBILITY_RULES, **SIGNIFICANCE_CONDITIONED_DIAGNOSTICS,
            **UNCONDITIONED_DIAGNOSTICS}

OBJECTIVE_NAMES = ["loo_mcc", "perm_evidence_zg", "fold_enrichment", "neighbor_stability"]


@dataclass
class RadiusResult:
    """Everything computed at ONE radius, before any cross-radius comparison."""

    radius_A: float
    centers: list[int]                    # provisional S(r): residue indices, sorted
    center_rows: list[dict]
    loo: LooOutcome
    perm: PermutationOutcome              # SECONDARY label null (Zg only)
    positional: PositionalOutcome         # PRIMARY null — decides S(r)
    fold_enrichment: float
    fold_enrichment_defined: bool
    n_plp_in: int
    n_labeled_in: int
    coverage_fraction: float
    n_singleton_centers: int
    median_n_labeled: float
    qc_failures: list[str]                # admissibility ONLY (v2 §4)
    diagnostics_fired: list[str]          # recorded, never binding
    diagnostic_values: dict
    boundary_neutral_verified: bool
    bh_boundary_p: float | None
    n_in_family: int
    neighbor_stability: float = 0.0       # filled by the cross-radius pass (II.5D)
    #: False when S(r) and BOTH neighbours are empty: the 0.0 above is then a
    #: sentinel for an undefined quantity, not a measured dissimilarity. Reported
    #: so a reader can tell the two apart; it does NOT change admissibility, which
    #: v2 explicitly forbids conditioning on significance (QC_H1 is a diagnostic).
    neighbor_stability_defined: bool = True
    certificate: object | None = None     # per-radius power certificate (v2 §5.4)
    # --- v2 §4: biological-plausibility DIAGNOSTICS (never Pareto objectives) --
    structural_coverage: float = 0.0      # == coverage_fraction; U_struct fraction covered
    cohort_absorption: float = 0.0        # fraction of L in the single largest sphere
    r_to_domain_ratio: float = 0.0        # radius / radius of gyration of U_struct
    extra: dict = field(default_factory=dict)

    @property
    def admissible(self) -> bool:
        return not self.qc_failures

    @property
    def qc_failure_reason(self) -> str:
        return ";".join(self.qc_failures) if self.qc_failures else "NA"

    @property
    def admissibility_reason(self) -> str:
        if not self.qc_failures:
            return "NA"
        return " | ".join(ADMISSIBILITY_RULES.get(rule, rule) for rule in self.qc_failures)

    @property
    def diagnostics_reason(self) -> str:
        return ";".join(self.diagnostics_fired) if self.diagnostics_fired else "NA"

    def objectives(self) -> dict[str, float]:
        return {
            "loo_mcc": float(self.loo.mcc),
            "perm_evidence_zg": float(self.perm.zg),
            # v2 §4: a radius with no significant center scores ZERO here; it is not
            # struck from the Pareto set for having found nothing.
            "fold_enrichment": float(self.fold_enrichment),
            "neighbor_stability": float(self.neighbor_stability),
        }


def scan_one_radius(radius: float, universe_coords: np.ndarray, universe_index: np.ndarray,
                    labeled_positions: np.ndarray, labeled_coords: np.ndarray,
                    y: np.ndarray, rng: np.random.Generator, *, B: int, q: float,
                    kappa: float, tie_tolerance: float, sparse_cutoff: int,
                    undefined_mcc: float, qc_h2_max_coverage: float,
                    qc_h3_min_median: float, qc_h4_max_isolated: float,
                    r_floor: float, r_cap: float, seed_context: str,
                    label_rng: np.random.Generator | None = None,
                    chunk: int = 2000,
                    center_coords: np.ndarray | None = None,
                    center_index: np.ndarray | None = None,
                    radius_of_gyration_A: float | None = None) -> RadiusResult:
    """One completely independent radius evaluation (II.5 A-D plus v2 §4).

    ``rng`` drives the PRIMARY positional null; ``label_rng`` drives the SECONDARY
    label null and must come from its own derived seed context so the two never share
    a random stream (v2 §5.5). Stage B always supplies both; a direct library caller
    that omits ``label_rng`` gets the primary generator, which is still deterministic
    but couples the two streams and is not how the stage runs.

    ``center_coords``/``center_index`` restrict the CANDIDATE centers to v2 §2's
    ``U_center`` (a subset of ``universe_coords``/``universe_index``); the universe
    arguments remain the FULL positional reference ``U_struct``, used for the
    positional-null population and for coverage. Defaulting to the universe
    preserves the pre-v2 (centers == U_struct) behaviour for direct library callers.
    """
    universe_coords = np.asarray(universe_coords, dtype=np.float64)
    labeled_coords = np.asarray(labeled_coords, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if center_coords is None:
        center_coords, center_index = universe_coords, universe_index
    center_coords = np.asarray(center_coords, dtype=np.float64)
    n_centers, N = len(center_coords), len(y)
    n_p = float(y.sum())
    pi = n_p / N

    # --- PRIMARY test: structure-aware positional null (v2 §5.2) ---------------
    # U is (n_centers candidates x M = |U_struct|): centers may be a strict subset
    # of the FULL positional reference the placement is drawn over.
    U = (cross_distances(center_coords, universe_coords) <= radius).astype(np.uint8)
    A = membership_matrix(center_coords, labeled_coords, radius)
    n_labeled = A.sum(axis=1)
    plp_positions = np.asarray(labeled_positions, dtype=np.int64)[y > 0]
    positional = positional_pass(U, plp_positions, n_labeled, rng, B, radius,
                                 seed_context)

    # --- SECONDARY: label permutation, for T(r)/Zg (II.5B) only ---------------
    perm = permutation_pass(A, y, label_rng if label_rng is not None else rng, B,
                            radius, f"{seed_context}|secondary_label_null", chunk=chunk)

    # --- provisional S(r): the same BH machinery as the final detection --------
    family = np.nonzero(positional.in_family)[0]
    fdr = benjamini_hochberg(positional.p_emp[family], q)
    q_bh_full = np.ones(n_centers, dtype=np.float64)
    q_bh_full[family] = fdr.q_bh
    significant = np.zeros(n_centers, dtype=bool)
    significant[family] = fdr.reject_bh
    sel = np.nonzero(significant)[0]

    z_center = positional.z_center
    center_rows = [
        {
            "center_residue_index": int(center_index[c]),
            "x_ca": float(center_coords[c, 0]), "y_ca": float(center_coords[c, 1]),
            "z_ca": float(center_coords[c, 2]),
            "n_labeled_in_sphere": int(positional.n_labeled[c]),
            "n_universe_in_sphere": int(positional.n_universe[c]),
            "n_plp_in_sphere": int(positional.n_plp_obs[c]),
            "n_blb_in_sphere": int(positional.n_labeled[c] - positional.n_plp_obs[c]),
            "null_mean_plp": float(positional.null_mean[c]),
            "null_sd_plp": float(positional.null_sd[c]),
            "z_score": float(z_center[c]), "p_emp": float(positional.p_emp[c]),
            "p_exact_hypergeom": float(positional.p_exact[c]),
            "q_bh": float(q_bh_full[c]), "in_test_family": bool(positional.in_family[c]),
            "exclusion_reason": positional.exclusion_reason[c],
            "null_model": positional.null_model,
        }
        for c in sel
    ]

    # --- II.5A LOO-MCC (selection objective only) ------------------------------
    loo = leave_one_out(labeled_coords, y.astype(np.int64), radius, kappa,
                        tie_tolerance, sparse_cutoff, undefined_mcc)
    boundary_neutral = kappa_is_boundary_neutral(loo.n_plp, loo.n_labeled,
                                                 loo.pi_minus, kappa, tie_tolerance)

    # --- II.5C fold enrichment + coverage --------------------------------------
    if len(sel):
        s_coords = center_coords[sel]
        covered_u = covered_mask(s_coords, universe_coords, radius)
        covered_l = covered_mask(s_coords, labeled_coords, radius)
        coverage = float(covered_u.mean())
        n_labeled_in = int(covered_l.sum())
        n_plp_in = int(y[covered_l].sum())
        d_centers = cross_distances(s_coords, s_coords)
        np.fill_diagonal(d_centers, np.inf)
        n_singleton = int((d_centers.min(axis=1) > 2.0 * radius).sum())
        isolated_fraction = n_singleton / len(sel)
    else:
        coverage, n_labeled_in, n_plp_in = 0.0, 0, 0
        n_singleton, isolated_fraction = 0, 0.0

    if n_labeled_in > 0 and pi > 0:
        fold_enrichment, fe_defined = (n_plp_in / n_labeled_in) / pi, True
    elif len(sel) == 0:
        # v2 §4 — no significant center is not an undefined quantity. The radius
        # scores zero on the enrichment objective and remains a Pareto candidate.
        fold_enrichment, fe_defined = 0.0, False
    else:
        fold_enrichment, fe_defined = 0.0, False

    # --- v2 §4 biological-plausibility DIAGNOSTICS (never objectives, F15) -----
    # structural_coverage IS coverage as defined above (fraction of U_struct within
    # r of any significant center); cohort_absorption is the fraction of the FULL
    # classified cohort L held by the single largest significant sphere.
    structural_coverage = coverage
    cohort_absorption = (max(row["n_labeled_in_sphere"] for row in center_rows) / N
                         if center_rows and N else 0.0)
    r_to_domain_ratio = (float(radius) / radius_of_gyration_A
                         if radius_of_gyration_A else float("nan"))

    # --- v2 §4 admissibility: ONLY genuinely undefined quantities ---------------
    median_n_labeled = (float(np.median(positional.n_labeled)) if n_centers
                        else float("nan"))
    failures: list[str] = []
    if len(sel) and n_labeled_in == 0:
        failures.append("FE_UNDEF")
    if not (r_floor <= radius <= r_cap + 1e-9):
        failures.append("QC_H5_DOMAIN")
    if n_centers == 0 or not np.isfinite(median_n_labeled):
        failures.append("QC_H3_UNDEFINED_MEDIAN")

    # --- v2 §4 diagnostics: recorded, reported, never binding ------------------
    diagnostics: list[str] = []
    if len(sel) == 0:
        diagnostics.append("QC_H1")
    if coverage > qc_h2_max_coverage:
        diagnostics.append("QC_H2")
    if median_n_labeled < qc_h3_min_median:
        diagnostics.append("QC_H3_THIN_EVIDENCE")
    if len(sel) and isolated_fraction > qc_h4_max_isolated:
        diagnostics.append("QC_H4")

    return RadiusResult(
        radius_A=float(radius),
        centers=[int(center_index[c]) for c in sel],
        center_rows=center_rows, loo=loo, perm=perm, positional=positional,
        fold_enrichment=fold_enrichment, fold_enrichment_defined=fe_defined,
        n_plp_in=n_plp_in, n_labeled_in=n_labeled_in,
        coverage_fraction=coverage, n_singleton_centers=n_singleton,
        structural_coverage=structural_coverage, cohort_absorption=cohort_absorption,
        r_to_domain_ratio=r_to_domain_ratio,
        median_n_labeled=median_n_labeled, qc_failures=failures,
        diagnostics_fired=diagnostics,
        diagnostic_values={
            "QC_H1_zero_significant_centers": bool(len(sel) == 0),
            "QC_H2_coverage_fraction": coverage,
            "QC_H2_max_coverage": qc_h2_max_coverage,
            "QC_H2_exceeds_max": bool(coverage > qc_h2_max_coverage),
            "QC_H3_median_n_labeled": median_n_labeled,
            "QC_H3_min_median_n_labeled": qc_h3_min_median,
            "QC_H3_below_threshold": bool(median_n_labeled < qc_h3_min_median),
            "QC_H4_isolated_center_fraction": isolated_fraction,
            "QC_H4_max_isolated_center_fraction": qc_h4_max_isolated,
            "QC_H4_exceeds_max": bool(len(sel) and isolated_fraction > qc_h4_max_isolated),
        },
        boundary_neutral_verified=boundary_neutral, bh_boundary_p=fdr.boundary_p,
        n_in_family=int(positional.in_family.sum()),
        extra={
            "pi_prevalence": pi, "isolated_center_fraction": isolated_fraction,
            "n_plp_in_hotspot_spheres": n_plp_in,
            "n_labeled_in_hotspot_spheres": n_labeled_in,
            "T_obs": perm.t_obs, "T_null_mean": perm.t_null_mean,
            "T_null_sd": perm.t_null_sd, "p_ratio_descriptive_only": perm.p_ratio,
            "n_in_test_family": int(positional.in_family.sum()),
            "bh_boundary_p": fdr.boundary_p, "seed_context": seed_context,
            "primary_null": positional.null_model,
            "secondary_null": "label_permutation",
            "smallest_p_emp_in_family": (float(positional.p_emp[family].min())
                                         if len(family) else None),
            "smallest_p_exact_in_family": (float(positional.p_exact[family].min())
                                           if len(family) else None),
        },
    )


def attach_neighbor_stability(results: list[RadiusResult], step: float) -> None:
    """II.5D — mean Jaccard of S(r) against S(r +/- step), on CENTER SETS.

    Computed as a cross-radius pass *after* every radius has been evaluated
    independently, so no radius influences another's primary statistics. Both sets
    empty -> 0.0 by the recorded convention in ``utils.geometry.jaccard``.
    """
    by_radius = {round(r.radius_A, 6): r for r in results}
    for r in results:
        sets = []
        for offset in (-step, step):
            neighbour = by_radius.get(round(r.radius_A + offset, 6))
            if neighbour is not None:
                sets.append(jaccard(set(r.centers), set(neighbour.centers)))
        r.neighbor_stability = float(np.mean(sets)) if sets else 0.0
        # Defined only if at least one comparison had a non-empty side. All-empty ->
        # jaccard(empty, empty) = 0.0 by convention, which is a sentinel, not a
        # measurement; flag it so the scan table does not present the two alike.
        r.neighbor_stability_defined = bool(sets) and any(
            r.centers or by_radius[round(r.radius_A + o, 6)].centers
            for o in (-step, step)
            if by_radius.get(round(r.radius_A + o, 6)) is not None)
