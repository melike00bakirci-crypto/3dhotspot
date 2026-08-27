"""Workflow v2 §3/§5.2 — the structure-aware POSITIONAL null (the PRIMARY test).

The per-center test uses **the same null as the global analysis of §3**: the P/LP
point set is redrawn uniformly **without replacement** over ``U_struct``,
conditioned on the observed ``N_P``, and the observed count of P/LP residues inside
the sphere is compared with its permutation distribution.

    H0 : the N_P pathogenic residues are placed uniformly at random over U_struct.
    n_P^(b)(c) = |{ drawn positions } ∩ { residues within r of c }|
    p_emp(c)   = (1 + #{b : n_P^(b)(c) >= n_P^obs(c)}) / (1 + B)

This replaces the v1 label-permutation null, which is retained by
:mod:`hotspot3d.hotspot.permutation` as an explicitly SECONDARY analysis testing a
different hypothesis. v2 §5.2 forbids mixing the two across the stages of one run: a
workflow that declares global clustering under one null and then fails to detect it
under another has produced an inconsistency, not a negative result.

**Why this matters numerically.** Under the positional null the in-sphere count is
*exactly* hypergeometric — ``n_P(c) ~ Hypergeom(M = |U_struct|, K = n_U(c),
n = N_P)`` — so the smallest p-value the test can attain is the hypergeometric tail
at the most extreme attainable occupancy, which is typically many orders of
magnitude below the corresponding label-permutation floor
``C(N_P, n_L) / C(N_P + N_B, n_L)``. That difference is the whole of the KCNA2
finding recorded in v2 Appendix B, and it is why §5.4's power certificate must be
evaluated **under the null actually used**.

The exact tail is computed here alongside the empirical p-value and published as a
diagnostic. It is **not** the decision statistic: BH operates on ``p_emp``, whose
resolution floor ``1/(B+1)`` is exactly what §5.4's first floor governs. Publishing
both makes the resolution floor visible instead of implicit.

Bit-reproducibility: all B draws are taken in one call before any chunking, so the
chunk size is invisible in the results; every accumulated count is an exact integer.

**[NEW — Workflow v2 §2] Candidate centers vs. the placement population.** The set
of positions ``H0`` places ``N_P`` residues over is ALWAYS the full ``U_struct``
(the positional reference is unchanged). The set of centers the test is actually
RUN AT is ``U_center`` — residues of ``U_struct`` with ``pLDDT >=
plddt.center_universe_min_plddt`` — which may be a strict subset. The two are
therefore decoupled here: the membership matrix ``U`` need not be square, and its
row count (candidate centers) and column count (the population ``M = |U_struct|``
draws are taken from) are tracked separately. The null formula above is otherwise
untouched.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import lgamma

import numpy as np

from ..spatial.ripley import positional_null_draws

NULL_MODEL = "structure_aware_positional"

#: Draw-chunk sizing. Chosen so the gathered (M, chunk, N_P) uint8 block stays a few
#: MB; it is a memory device only and cannot influence a single bit of the result.
_TARGET_GATHER_ELEMENTS = 4_000_000


@dataclass
class PositionalOutcome:
    """Everything one radius needs from one positional-null pass."""

    radius_A: float
    B: int
    seed_context: str
    null_model: str
    n_universe: np.ndarray       # n_U(c) — residues of U_struct within r of c
    n_labeled: np.ndarray        # n_L(c) — classified residues within r of c
    n_plp_obs: np.ndarray        # n_P^obs(c)
    null_mean: np.ndarray        # mean_b n_P^(b)(c)
    null_sd: np.ndarray          # sd_b n_P^(b)(c) -> sigma_c
    ge_count: np.ndarray         # #{b : n_P^(b)(c) >= n_P^obs(c)}
    p_emp: np.ndarray            # (1 + ge) / (1 + B) — the DECISION statistic
    p_exact: np.ndarray          # exact hypergeometric tail — diagnostic only
    p_comb_floor: np.ndarray     # smallest p this center could ever return (§5.4)
    in_family: np.ndarray        # bool: n_L(c) > 0 AND sigma_c > 0
    exclusion_reason: list[str]
    N_P: int
    M: int

    @property
    def z_center(self) -> np.ndarray:
        out = np.zeros_like(self.n_plp_obs, dtype=np.float64)
        ok = self.null_sd > 0
        out[ok] = (self.n_plp_obs[ok] - self.null_mean[ok]) / self.null_sd[ok]
        return out

    @property
    def expected_plp(self) -> np.ndarray:
        """``N_P * n_U(c) / M`` — the analytic null mean, for cross-checking."""
        return self.N_P * self.n_universe.astype(np.float64) / max(self.M, 1)


# --- exact hypergeometric tails ---------------------------------------------
# Implemented locally rather than imported so that the floor in the power
# certificate (§5.4) has no optional-dependency branch: the number that decides
# whether a run may be reported as a negative must always be computable.

def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -np.inf
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def hypergeom_sf(k: int, M: int, K: int, n: int) -> float:
    """``P(X >= k)`` for ``X ~ Hypergeom(M, K, n)`` — n draws without replacement.

    ``M`` is the population, ``K`` the number of successes in it, ``n`` the draw
    size. Summed in the log domain over the (always short) support, so the result
    is accurate down to the smallest floors this pipeline reports.
    """
    if M <= 0 or n <= 0 or K <= 0:
        return 1.0 if k <= 0 else 0.0
    k_min, k_max = max(0, n - (M - K)), min(K, n)
    if k <= k_min:
        return 1.0                       # every realization satisfies X >= k
    if k > k_max:
        return 0.0                       # no realization can reach k
    logs = np.array([_log_comb(K, i) + _log_comb(M - K, n - i) - _log_comb(M, n)
                     for i in range(int(k), k_max + 1)], dtype=np.float64)
    top = float(logs.max())
    return float(np.exp(top) * np.exp(logs - top).sum())


def positional_p_exact(n_plp_obs: np.ndarray, n_universe: np.ndarray, M: int,
                       N_P: int) -> np.ndarray:
    """Exact positional-null tail per center. Cached by ``(n_U, n_obs)``."""
    cache: dict[tuple[int, int], float] = {}
    out = np.ones(len(n_plp_obs), dtype=np.float64)
    for c in range(len(n_plp_obs)):
        key = (int(n_universe[c]), int(n_plp_obs[c]))
        if key not in cache:
            cache[key] = hypergeom_sf(key[1], M, key[0], N_P)
        out[c] = cache[key]
    return out


def positional_combinatorial_floor(n_universe: np.ndarray, M: int,
                                   N_P: int) -> np.ndarray:
    """Smallest p each center could EVER return under the positional null (§5.4).

    Attained at the maximal occupancy ``k_max = min(N_P, n_U(c))``. Note this is
    **not** monotone in ``n_U``: it decreases up to ``n_U = N_P`` and increases
    afterwards, so "the largest sphere" is not automatically the most favourable
    configuration — the minimum must be taken over the occupancies that actually
    occur.
    """
    cache: dict[int, float] = {}
    out = np.ones(len(n_universe), dtype=np.float64)
    for c in range(len(n_universe)):
        n_u = int(n_universe[c])
        if n_u not in cache:
            cache[n_u] = hypergeom_sf(min(N_P, n_u), M, n_u, N_P)
        out[c] = cache[n_u]
    return out


# --- the pass ----------------------------------------------------------------

def positional_pass(U: np.ndarray, plp_positions: np.ndarray, n_labeled: np.ndarray,
                    rng: np.random.Generator, B: int, radius_A: float,
                    seed_context: str, chunk: int | None = None) -> PositionalOutcome:
    """One structure-aware positional null at one radius (v2 §5.2).

    ``U`` is the ``(n_centers candidate rows x M = |U_struct| population columns)``
    binary membership matrix — the geometry, held FIXED. Rows are the CANDIDATE
    centers under test (v2 §2 restricts these to ``U_center``); columns are the
    FULL positional reference ``U_struct``, unchanged. ``plp_positions`` are
    positions of the observed P/LP residues WITHIN the column space (``U_struct``);
    only their *positions* are randomized, over the full population, exactly as
    before — ``n_centers`` and ``M`` are decoupled so ``U`` need not be square.
    """
    U = np.asarray(U)
    if U.dtype != np.uint8:
        U = (U > 0).astype(np.uint8)
    n_centers, M = U.shape[0], U.shape[1]
    plp_positions = np.asarray(plp_positions, dtype=np.int64)
    N_P = int(len(plp_positions))

    n_universe = U.sum(axis=1, dtype=np.int64)
    n_plp_obs = U[:, plp_positions].sum(axis=1, dtype=np.int64).astype(np.float64)
    n_labeled = np.asarray(n_labeled, dtype=np.float64)

    # Every draw taken BEFORE chunking: the chunk size is a memory device only.
    # The draw pool is the POPULATION M = |U_struct|, never the candidate-center count.
    draws = positional_null_draws(rng, M, N_P, B)

    if chunk is None:
        chunk = max(1, min(B, _TARGET_GATHER_ELEMENTS // max(n_centers * max(N_P, 1), 1)))

    total = np.zeros(n_centers, dtype=np.float64)
    total_sq = np.zeros(n_centers, dtype=np.float64)
    ge_count = np.zeros(n_centers, dtype=np.int64)

    done = 0
    while done < B:
        size = min(chunk, B - done)
        idx = draws[done:done + size]                      # (size, N_P)
        counts = U[:, idx].sum(axis=2, dtype=np.int32)     # (n_centers, size) — exact integers
        counts = counts.astype(np.float64)
        total += counts.sum(axis=1)
        total_sq += (counts * counts).sum(axis=1)
        ge_count += (counts >= n_plp_obs[:, None]).sum(axis=1)
        done += size

    mean = total / B
    var = np.maximum(total_sq / B - mean ** 2, 0.0) * (B / (B - 1) if B > 1 else 1.0)
    sd = np.sqrt(var)
    p_emp = (1.0 + ge_count) / (1.0 + B)

    reasons: list[str] = []
    in_family = np.zeros(n_centers, dtype=bool)
    for c in range(n_centers):
        if n_labeled[c] == 0:
            # v2 §5.1 — the family is the centers whose sphere holds >= 1 classified
            # residue. m sets the significance bar and is reported with every q.
            reasons.append("no_labeled_residue_in_sphere")
        elif sd[c] == 0:
            reasons.append("zero_permutation_variance")
        else:
            reasons.append("NA")
            in_family[c] = True

    return PositionalOutcome(
        radius_A=float(radius_A), B=int(B), seed_context=seed_context,
        null_model=NULL_MODEL, n_universe=n_universe, n_labeled=n_labeled,
        n_plp_obs=n_plp_obs, null_mean=mean, null_sd=sd, ge_count=ge_count,
        p_emp=p_emp, p_exact=positional_p_exact(n_plp_obs, n_universe, M, N_P),
        p_comb_floor=positional_combinatorial_floor(n_universe, M, N_P),
        in_family=in_family, exclusion_reason=reasons, N_P=N_P, M=M,
    )
