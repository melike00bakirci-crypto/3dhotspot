"""The label-permutation null — the SECONDARY analysis (v2 §5.2) and II.5B evidence.

**This is no longer the primary per-center test.** Workflow v2 §5.2 requires the
primary test to use the same structure-aware positional null as the global analysis
of §3; that null lives in :mod:`hotspot3d.hotspot.positional`. Label permutation
tests a *different* hypothesis — "are P/LP residues more clustered than B/LB
residues, GIVEN the set of variant-bearing positions?" — and is permitted only as a
secondary, reported analysis, never as the primary test and never without the §5.4
power certificate. The two are reported side by side and never mixed.

It is retained for two uses, both of which are explicitly not primary inference:

  * the radius-selection objective ``Zg`` derived from ``T(r)`` (II.5B), which is one
    of four equally weighted objectives and never a significance test;
  * the secondary per-center analysis published beside the primary result at
    ``r_hot``.

Its combinatorial floor ``C(N_P, n_L) / C(N_P + N_B, n_L)`` is typically orders of
magnitude coarser than the positional one at the same cohort size — that gap is the
whole of the KCNA2 finding in v2 Appendix B — so any use of this null must carry the
§5.4 certificate computed under it.

The null holds the **geometry fixed** and permutes only the P/LP vs B/LB *labels*
across the labelled cohort ``L``, so ``N_P`` and ``N_B`` are preserved exactly in
every replicate (F4). One pass over the permutations yields, simultaneously:

  * per-center pathogenic counts ``n_P(c)`` for every candidate center, giving the
    add-one empirical p-value of II.3, and
  * the radius-level ratio statistic ``T(r) = sum_{c : n_L(c) >= 1} n_P(c)/n_L(c)``
    and its standardized evidence ``Zg(r)`` (II.5B).

Both statistics therefore come from the *same* permutation draw at a given radius,
which is both cheaper and scientifically tidier than two independent nulls.

Implementation note (the II.3 matrix trick): with the binary membership matrix
``A`` (centers x labelled residues), ``n_P`` for ALL centers under one permutation
is the single matrix-vector product ``A @ y_perm``; a whole chunk of permutations
is one matrix-matrix product. There is no Python loop over permutations.

Bit-reproducibility: ``A`` and ``Y`` are 0/1 matrices, so every partial sum in
``A @ Y`` is an exact integer far below ``2**53``. Any BLAS blocking or threading
order therefore yields bit-identical results. Reductions whose operands are *not*
exact integers (the ratio statistic, means, sds) are computed with NumPy's own
deterministic pairwise summation rather than BLAS.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# n_P is compared to the observed count with ">=", the conservative direction
# required by the add-one (Phipson-Smyth) estimator: p is never zero.
_DEFAULT_CHUNK = 2000


@dataclass
class PermutationOutcome:
    """Everything one radius needs from one permutation pass."""

    radius_A: float
    B: int
    seed_context: str
    n_labeled: np.ndarray        # n_L(c), (M,)
    n_plp_obs: np.ndarray        # n_P^obs(c), (M,)
    null_mean: np.ndarray        # mean_b n_P^(b)(c)
    null_sd: np.ndarray          # sd_b n_P^(b)(c)  -> sigma_c
    ge_count: np.ndarray         # #{b : n_P^(b)(c) >= n_P^obs(c)}
    p_emp: np.ndarray            # (1 + ge) / (1 + B)
    in_family: np.ndarray        # bool: n_L(c) > 0 AND sigma_c > 0
    exclusion_reason: list[str]
    t_obs: float
    t_null_mean: float
    t_null_sd: float
    zg: float
    p_ratio: float               # descriptive only — it saturates at 1/(B+1)

    @property
    def z_center(self) -> np.ndarray:
        out = np.zeros_like(self.n_plp_obs, dtype=np.float64)
        ok = self.null_sd > 0
        out[ok] = (self.n_plp_obs[ok] - self.null_mean[ok]) / self.null_sd[ok]
        return out


def permute_labels(rng: np.random.Generator, y: np.ndarray, n: int) -> np.ndarray:
    """``n`` independent label permutations as an ``(N, n)`` int8 matrix.

    ``Generator.permuted`` shuffles each row independently, so every column of the
    returned matrix is a permutation of ``y`` — ``N_P`` and ``N_B`` are preserved
    exactly, by construction rather than by assertion.

    All ``B`` permutations are drawn in ONE call, before any chunking, so the chunk
    size used for the matrix products cannot influence the random stream and is
    therefore invisible in the results.
    """
    tiled = np.tile(np.asarray(y, dtype=np.int8), (n, 1))
    return rng.permuted(tiled, axis=1).T


def permutation_pass(A: np.ndarray, y: np.ndarray, rng: np.random.Generator, B: int,
                     radius_A: float, seed_context: str,
                     chunk: int = _DEFAULT_CHUNK) -> PermutationOutcome:
    """One label-permutation null at one radius — SECONDARY analysis plus II.5B.

    ``A`` is the ``(M centers x N labelled residues)`` binary membership matrix.
    The per-center ``p_emp`` returned here is never the primary decision statistic
    (v2 §5.2); primary significance comes from
    :func:`hotspot3d.hotspot.positional.positional_pass`.
    """
    A = np.asarray(A, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    M = A.shape[0]

    n_labeled = A.sum(axis=1)
    n_plp_obs = A @ y

    # II.5B weights: centers with no labelled residue contribute nothing to T(r).
    w = np.zeros(M, dtype=np.float64)
    informative = n_labeled >= 1
    w[informative] = 1.0 / n_labeled[informative]
    t_obs = float(np.sum(w * n_plp_obs))

    total = np.zeros(M, dtype=np.float64)
    total_sq = np.zeros(M, dtype=np.float64)
    ge_count = np.zeros(M, dtype=np.int64)
    t_values = np.empty(B, dtype=np.float64)

    # Draw every permutation first: chunking below is a memory device only.
    Y_all = permute_labels(rng, y, B)                    # (N, B) int8

    done = 0
    while done < B:
        size = min(chunk, B - done)
        Y = Y_all[:, done:done + size].astype(np.float64)
        NP = A @ Y                                       # (M, size) — exact integers
        total += NP.sum(axis=1)
        total_sq += (NP * NP).sum(axis=1)
        ge_count += (NP >= n_plp_obs[:, None]).sum(axis=1)
        t_values[done:done + size] = (w[:, None] * NP).sum(axis=0)
        done += size

    mean = total / B
    var = np.maximum(total_sq / B - mean ** 2, 0.0) * (B / (B - 1) if B > 1 else 1.0)
    sd = np.sqrt(var)

    p_emp = (1.0 + ge_count) / (1.0 + B)

    reasons: list[str] = []
    in_family = np.zeros(M, dtype=bool)
    for c in range(M):
        if n_labeled[c] == 0:
            reasons.append("no_labeled_residue_in_sphere")
        elif sd[c] == 0:
            reasons.append("zero_permutation_variance")
        else:
            reasons.append("NA")
            in_family[c] = True

    t_mean = float(t_values.mean())
    t_sd = float(t_values.std(ddof=1)) if B > 1 else 0.0

    # T(r) mixes exact integer counts with the inexact weights 1/n_L, so a statistic
    # that is mathematically CONSTANT under permutation (it is, whenever the spheres
    # partition the cohort into disjoint blocks) still shows float noise. Treating that
    # noise as evidence would be indefensible, so a numerically degenerate T is
    # recorded as such: Zg = 0 and the descriptive p = 1.
    scale = max(1.0, abs(t_obs))
    degenerate = t_sd <= 1e-12 * scale
    zg = 0.0 if degenerate else float((t_obs - t_mean) / t_sd)
    if degenerate:
        p_ratio = 1.0
    else:
        p_ratio = float((1 + int((t_values >= t_obs - 1e-12 * scale).sum())) / (1 + B))

    return PermutationOutcome(
        radius_A=float(radius_A), B=int(B), seed_context=seed_context,
        n_labeled=n_labeled, n_plp_obs=n_plp_obs, null_mean=mean, null_sd=sd,
        ge_count=ge_count, p_emp=p_emp, in_family=in_family,
        exclusion_reason=reasons, t_obs=t_obs, t_null_mean=t_mean, t_null_sd=t_sd,
        zg=zg, p_ratio=p_ratio,
    )
