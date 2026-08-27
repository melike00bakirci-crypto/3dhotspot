"""Workflow v2 §5.4 — the POWER CERTIFICATE.

The mandatory precondition for any negative result, and the single check whose
absence let the KCNA2 run ``20260817T155124Z`` be misreported as a scientific
negative (v2 Appendix B).

Before the per-center test is run, and again after it, compute the **minimum
p-value the test could possibly return** and compare it with the hardest bar the
correction imposes:

    p_res        = 1 / (B + 1)                    permutation-resolution floor
    p_comb_best  = smallest attainable tail under THE NULL ACTUALLY USED
    p_floor      = max(p_res, p_comb_best)
    c_1          = q / m                          rank-1 BH critical value

``p_floor > c_1`` means **no center can be declared significant under any
realization of the data**. The run terminates ``UNDERPOWERED`` with a BLOCKING
``TEST_CANNOT_REJECT``; it is not a negative result and licenses no statement about
the gene (v2 Appendix A).

The two floors are remediated differently, and conflating them is how an
unremediable design gets recommended more permutations:

* ``p_res`` binds  -> report the ``B`` that would clear ``c_1``; re-run at that ``B``.
* ``p_comb_best`` binds -> **increasing B is futile**. The remedies are a larger
  benign cohort, a different null (§5.2), or a smaller test family (§2, §5.1).

The certificate never changes ``B``, never changes the FDR method, never widens a
domain and never converts an empty result into a positive one. It decides only
whether the run is *allowed to make a claim*.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..utils.geometry import cross_distances
from .positional import hypergeom_sf

LABEL_NULL = "label_permutation"
POSITIONAL_NULL = "structure_aware_positional"

BINDING_PERMUTATION = "permutation_resolution"
BINDING_COMBINATORIAL = "combinatorial"
BINDING_BOTH = "both"

CERTIFICATE_SCHEMA = "1.0.0"

UNINFORMATIVE_STATEMENT = (
    "Under the configuration executed, the test could not have rejected any center "
    "regardless of the data; the analysis is uninformative about the presence or "
    "absence of hotspots in this gene."
)

COMBINATORIAL_REMEDIES = [
    "a larger benign cohort (N_B sets the resolution of any P/LP-vs-B/LB contrast)",
    "a different null model (v2 §5.2 — the structure-aware positional null has a far "
    "lower combinatorial floor than the label-permutation null at the same cohort size)",
    "a smaller test family (v2 §2/§5.1 — m sets the bar c_1 = q/m directly)",
]


class UnderpoweredResult(Exception):
    """Terminal state ``UNDERPOWERED`` (v2 §0): the design could not have rejected.

    Deliberately NOT a subclass of the pipeline's negative-result or blocked-error
    types. It is neither a finding about the gene nor a fault in the code, and
    typing it as either is precisely the misreport v2 exists to prevent.
    """

    def __init__(self, certificate: "PowerCertificate", stage: str):
        self.certificate = certificate
        self.stage = stage
        super().__init__(certificate.headline())


# --- combinatorial floors ----------------------------------------------------

def label_null_floor(N_P: int, N_B: int, n_L: int) -> float:
    """Smallest p attainable by a sphere holding ``n_L`` classified residues.

    Label permutation over the cohort ``L``: the extreme configuration is
    ``k_max = min(n_L, N_P)`` pathogenic residues in the sphere, whose tail is
    ``C(N_P, n_L) / C(N_P + N_B, n_L)`` when ``n_L <= N_P``.
    """
    N = int(N_P) + int(N_B)
    n_L = int(n_L)
    if n_L <= 0 or N <= 0:
        return 1.0
    return hypergeom_sf(min(n_L, int(N_P)), N, int(N_P), n_L)


def positional_null_floor(M: int, n_U: int, N_P: int) -> float:
    """Smallest p attainable by a sphere covering ``n_U`` residues of ``U_struct``."""
    n_U = int(n_U)
    if n_U <= 0 or M <= 0 or N_P <= 0:
        return 1.0
    return hypergeom_sf(min(int(N_P), n_U), int(M), n_U, int(N_P))


def best_attainable_floor(null_model: str, *, N_P: int, N_B: int, M: int,
                          occupancies_labeled: np.ndarray,
                          occupancies_universe: np.ndarray) -> tuple[float, dict]:
    """``(p_comb_best, detail)`` — the minimum over the occupancies that OCCUR.

    v2 §5.4 asks for "the combinatorial floor evaluated at the most favourable
    sphere occupancy attainable at the selected radius" and, in the same paragraph,
    for "the minimum p-value the test could possibly return". Those coincide for the
    label null, whose floor decreases monotonically in ``n_L``; they do **not**
    coincide for the positional null, whose floor is minimised at ``n_U = N_P`` and
    rises again for larger spheres. The minimum over realised occupancies is
    therefore taken, and the value at the largest occupancy — the reading v2 spells
    out — is recorded alongside so both are auditable.
    """
    occ_l = np.asarray(occupancies_labeled, dtype=np.int64)
    occ_u = np.asarray(occupancies_universe, dtype=np.int64)
    if null_model == LABEL_NULL:
        values = sorted({int(v) for v in occ_l if v > 0})
        floors = {v: label_null_floor(N_P, N_B, v) for v in values}
        largest = max(values) if values else 0
        driver = "n_L"
    elif null_model == POSITIONAL_NULL:
        values = sorted({int(v) for v in occ_u if v > 0})
        floors = {v: positional_null_floor(M, v, N_P) for v in values}
        largest = max(values) if values else 0
        driver = "n_U"
    else:                                                   # pragma: no cover
        raise ValueError(f"unknown null model {null_model!r}")

    if not floors:
        return 1.0, {"reason": "no sphere holds any residue", "occupancy_driver": driver}

    best_occ = min(floors, key=lambda v: (floors[v], v))
    return floors[best_occ], {
        "occupancy_driver": driver,
        "best_occupancy": int(best_occ),
        "largest_occupancy": int(largest),
        "floor_at_largest_occupancy": float(floors[largest]),
        "n_distinct_occupancies": len(floors),
        "note": ("The floor is evaluated at every occupancy that actually occurs; the "
                 "minimum is the smallest p the test could return anywhere in the "
                 "family. For the positional null the floor is NOT monotone in n_U, "
                 "so the largest sphere is not automatically the most favourable."),
    }


def b_required_to_clear(c_1: float) -> int | None:
    """Smallest ``B`` with ``1/(B+1) <= c_1``. ``None`` when ``c_1 <= 0``."""
    if c_1 <= 0:
        return None
    b = int(math.ceil(1.0 / c_1)) - 1
    while b > 0 and 1.0 / (b + 1) > c_1:
        b += 1
    return max(b, 1)


# --- the certificate ---------------------------------------------------------

@dataclass
class PowerCertificate:
    """v2 §5.4 — a first-class Stage B artefact, emitted in EVERY terminal state."""

    phase: str                       # preflight | posthoc
    null_model: str
    B: int
    q: float
    m: int
    N_P: int
    N_B: int
    n_universe: int                  # |U_struct|
    n_center_universe: int           # |U_center|
    p_res: float
    p_comb_best: float
    p_floor: float
    c_1: float
    ratio: float
    binding_floor: str
    passes: bool
    B_to_clear_c_1: int | None
    increasing_B_is_futile: bool
    remedies: list[str]
    radius_A: float | None = None
    detail: dict = field(default_factory=dict)

    def headline(self) -> str:
        verdict = "PASSES" if self.passes else "FAILS — TEST_CANNOT_REJECT"
        return (f"power certificate ({self.phase}, null = {self.null_model}, "
                f"r = {self.radius_A if self.radius_A is not None else 'NA'}): "
                f"p_floor = {self.p_floor:.6g} vs c_1 = q/m = {self.c_1:.6g} "
                f"(ratio {self.ratio:.6g}, binding floor {self.binding_floor}, "
                f"m = {self.m}, N_P = {self.N_P}, N_B = {self.N_B}, B = {self.B}) "
                f"-> {verdict}")

    def as_json(self) -> dict:
        return {
            "schema_version": CERTIFICATE_SCHEMA,
            "phase": self.phase,
            "null_model": self.null_model,
            "null_model_source": "config:permutation.primary_null",
            "radius_A": self.radius_A,
            "N_P": self.N_P, "N_B": self.N_B,
            "n_universe_U_struct": self.n_universe,
            "n_center_universe_U_center": self.n_center_universe,
            "m_test_family_size": self.m,
            "B": self.B, "q": self.q,
            "p_res": self.p_res,
            "p_comb_best": self.p_comb_best,
            "p_floor": self.p_floor,
            "c_1_rank1_critical_value": self.c_1,
            "p_floor_over_c_1": self.ratio,
            "binding_floor": self.binding_floor,
            "PASSES": self.passes,
            "TEST_CANNOT_REJECT": not self.passes,
            "B_to_clear_c_1": self.B_to_clear_c_1,
            "increasing_B_is_futile": self.increasing_B_is_futile,
            "remedies": self.remedies,
            "headline": self.headline(),
            "interpretation": (
                "p_floor <= c_1: the design could return a p-value small enough for "
                "the rank-1 BH critical value, so an empty rejection set would be a "
                "statement about the data."
                if self.passes else UNINFORMATIVE_STATEMENT),
            "prohibited_inferences_if_failed": [
                "No 3D hotspot is detectable in this gene.",
                "This is a valid, complete scientific negative.",
                "This result licenses no modification of the method.",
            ],
            "definitions": {
                "p_res": "1 / (B + 1) — the permutation-resolution floor; remediable by B.",
                "p_comb_best": ("the smallest tail attainable under the null actually "
                                "used, given the cohort composition and the sphere "
                                "occupancies that occur; NOT remediable by B."),
                "c_1": "q / m — the rank-1 Benjamini-Hochberg critical value.",
            },
            **self.detail,
        }


def certify(phase: str, *, null_model: str, B: int, q: float, m: int, N_P: int,
            N_B: int, n_universe: int, n_center_universe: int,
            occupancies_labeled: np.ndarray, occupancies_universe: np.ndarray,
            radius_A: float | None = None, extra: dict | None = None
            ) -> PowerCertificate:
    """Evaluate v2 §5.4 for one (phase, radius, family) triple."""
    p_res = 1.0 / (B + 1)
    p_comb_best, comb_detail = best_attainable_floor(
        null_model, N_P=N_P, N_B=N_B, M=n_universe,
        occupancies_labeled=occupancies_labeled,
        occupancies_universe=occupancies_universe)
    p_floor = max(p_res, p_comb_best)

    if m <= 0:
        # No hypothesis is tested at all, so nothing can be rejected. Unreachable
        # while both classes are non-empty (a classified residue is its own center
        # and always sees itself), but it must never read as a vacuous PASS: an
        # infinite c_1 would make `p_floor <= c_1` true while the test does not exist.
        return PowerCertificate(
            phase=phase, null_model=null_model, B=int(B), q=float(q), m=0,
            N_P=int(N_P), N_B=int(N_B), n_universe=int(n_universe),
            n_center_universe=int(n_center_universe), p_res=float(p_res),
            p_comb_best=float(p_comb_best), p_floor=float(p_floor),
            c_1=float("inf"), ratio=float("inf"), binding_floor="empty_test_family",
            passes=False, B_to_clear_c_1=None, increasing_B_is_futile=True,
            remedies=["The test family is empty: no center's sphere contains a "
                      "classified residue, so no hypothesis exists to reject. Enlarging "
                      "B cannot help; the cohort or the radius domain must change."],
            radius_A=radius_A,
            detail={"combinatorial_floor_detail": comb_detail, **(extra or {})})

    c_1 = q / m

    if p_comb_best > p_res:
        binding = BINDING_COMBINATORIAL
    elif p_comb_best == p_res:
        binding = BINDING_BOTH
    else:
        binding = BINDING_PERMUTATION
    futile = binding in (BINDING_COMBINATORIAL, BINDING_BOTH)

    passes = bool(p_floor <= c_1)
    remedies: list[str] = []
    if not passes:
        if futile:
            remedies.extend(COMBINATORIAL_REMEDIES)
            remedies.append(
                "Increasing B is FUTILE: the binding floor is combinatorial, not the "
                "permutation resolution.")
        if binding in (BINDING_PERMUTATION, BINDING_BOTH):
            need = b_required_to_clear(c_1)
            remedies.append(
                f"Re-run at B = {need} under a NEW RUN_ID: that is the smallest B with "
                f"1/(B+1) <= c_1 = {c_1:.6g}. B is never raised inside a run.")

    return PowerCertificate(
        phase=phase, null_model=null_model, B=int(B), q=float(q), m=int(m),
        N_P=int(N_P), N_B=int(N_B), n_universe=int(n_universe),
        n_center_universe=int(n_center_universe), p_res=float(p_res),
        p_comb_best=float(p_comb_best), p_floor=float(p_floor), c_1=float(c_1),
        ratio=float(p_floor / c_1),
        binding_floor=binding, passes=passes,
        B_to_clear_c_1=b_required_to_clear(c_1),
        increasing_B_is_futile=futile, remedies=remedies, radius_A=radius_A,
        detail={"combinatorial_floor_detail": comb_detail, **(extra or {})},
    )


# --- pre-flight over the candidate domain ------------------------------------

@dataclass
class RadiusOccupancy:
    """Sphere occupancies at one candidate radius — no permutations spent."""

    radius_A: float
    n_labeled: np.ndarray
    n_universe: np.ndarray
    family: np.ndarray               # bool mask: n_L >= 1 (v2 §5.1)

    @property
    def m(self) -> int:
        return int(self.family.sum())


def domain_occupancies(universe_coords: np.ndarray, labeled_positions: np.ndarray,
                       radii, center_coords: np.ndarray | None = None
                       ) -> list[RadiusOccupancy]:
    """Occupancy counts for every candidate radius, from one distance matrix.

    Cheap by construction — this is what makes a pre-flight certificate affordable
    BEFORE any permutation compute is spent (v2 §5.4).

    ``center_coords`` restricts the CANDIDATE rows to v2 §2's ``U_center`` (a
    subset of ``universe_coords``); ``universe_coords`` itself remains the FULL
    positional reference ``U_struct`` used for ``n_universe`` and for locating the
    labelled positions. Defaults to ``universe_coords`` when omitted.
    """
    coords = np.asarray(universe_coords, dtype=np.float64)
    centers = (np.asarray(center_coords, dtype=np.float64)
              if center_coords is not None else coords)
    d = cross_distances(centers, coords)
    d_lab = d[:, np.asarray(labeled_positions, dtype=np.int64)]
    out = []
    for r in radii:
        n_u = (d <= float(r)).sum(axis=1).astype(np.int64)
        n_l = (d_lab <= float(r)).sum(axis=1).astype(np.int64)
        out.append(RadiusOccupancy(radius_A=float(r), n_labeled=n_l, n_universe=n_u,
                                   family=n_l >= 1))
    return out


def preflight(null_model: str, *, B: int, q: float, N_P: int, N_B: int,
              n_universe: int, n_center_universe: int,
              occupancies: list[RadiusOccupancy]) -> tuple[PowerCertificate, list[dict]]:
    """Per-radius certificates over the candidate domain (v2 §5.4 pre-flight).

    Returns ``(representative, per_radius_rows)``. The representative is the
    **most favourable** radius — the one whose ``p_floor / c_1`` is smallest. The
    run may only be terminated pre-flight when *no* radius in the domain could
    reject, so the certificate that decides is the best case, not the worst: an
    early termination must mean "no configuration in this domain could have
    worked", never "some configuration could not".
    """
    certificates = []
    for occ in occupancies:
        cert = certify(
            "preflight", null_model=null_model, B=B, q=q, m=occ.m, N_P=N_P, N_B=N_B,
            n_universe=n_universe, n_center_universe=n_center_universe,
            occupancies_labeled=occ.n_labeled[occ.family],
            occupancies_universe=occ.n_universe[occ.family],
            radius_A=occ.radius_A)
        certificates.append(cert)

    if not certificates:                                    # pragma: no cover
        raise ValueError("pre-flight requires at least one candidate radius")

    best = min(certificates, key=lambda c: (c.ratio, c.radius_A))
    rows = [{"radius_A": c.radius_A, "m": c.m, "p_res": c.p_res,
             "p_comb_best": c.p_comb_best, "p_floor": c.p_floor, "c_1": c.c_1,
             "p_floor_over_c_1": c.ratio, "binding_floor": c.binding_floor,
             "passes": c.passes} for c in certificates]
    best.detail["per_radius"] = rows
    best.detail["n_radii_evaluated"] = len(rows)
    best.detail["n_radii_passing"] = int(sum(1 for c in certificates if c.passes))
    best.detail["representative_radius_rule"] = (
        "the radius with the smallest p_floor / c_1 over the candidate domain — the "
        "most favourable configuration the domain offers. Pre-flight termination "
        "requires that even this radius cannot reject.")
    return best, rows
