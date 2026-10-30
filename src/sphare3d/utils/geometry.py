"""Shared geometric primitives — the SINGLE implementation of every distance.

Placed in ``utils`` (Lead-owned) precisely so that no two stages can implement
the same formula differently. F1 is frozen: residue representation is CA and
distance is CA-CA Euclidean. Every stage that needs a distance calls this module.

Cross-package *imports* are allowed; cross-package *edits* are not.
"""
from __future__ import annotations

import numpy as np


def pairwise_distances(coords: np.ndarray) -> np.ndarray:
    """Full symmetric CA-CA Euclidean distance matrix (F1)."""
    coords = np.asarray(coords, dtype=np.float64)
    diff = coords[:, None, :] - coords[None, :, :]
    return np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))


def cross_distances(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Euclidean distances between every row of ``a`` and every row of ``b``."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    diff = a[:, None, :] - b[None, :, :]
    return np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))


def protein_diameter(coords: np.ndarray) -> float:
    """``D_max`` — maximum pairwise CA-CA distance in U_struct (METHOD_SPEC II.0)."""
    if len(coords) < 2:
        return 0.0
    return float(pairwise_distances(coords).max())


def r_max_for_estimators(d_max: float) -> int:
    """``floor(D_max / 2)`` — Ripley/PCF estimates are unreliable beyond this (II.1)."""
    return int(np.floor(d_max / 2.0))


def membership_matrix(centers: np.ndarray, points: np.ndarray, radius: float) -> np.ndarray:
    """Binary center x point sphere-membership matrix ``A`` (II.3 implementation note).

    With ``A`` in hand, ``n_P`` for ALL centers under one label permutation is a
    single matrix-vector product ``A @ y_perm`` — the whole scan is a handful of
    GFLOPs rather than a nested loop.
    """
    return (cross_distances(centers, points) <= radius).astype(np.float64)


def neighbours_within(coords: np.ndarray, radius: float,
                      exclude_self: bool = True) -> np.ndarray:
    """Boolean adjacency for ``d <= radius`` over one coordinate set."""
    adj = pairwise_distances(coords) <= radius
    if exclude_self:
        np.fill_diagonal(adj, False)
    return adj


def covered_mask(centers: np.ndarray, points: np.ndarray, radius: float) -> np.ndarray:
    """Boolean mask of points within ``radius`` of at least one center."""
    if len(centers) == 0:
        return np.zeros(len(points), dtype=bool)
    return (cross_distances(points, centers) <= radius).any(axis=1)


def connected_components(adjacency: np.ndarray) -> list[list[int]]:
    """Connected components of a boolean adjacency matrix (deterministic order)."""
    n = adjacency.shape[0]
    seen = np.zeros(n, dtype=bool)
    components: list[list[int]] = []
    for start in range(n):
        if seen[start]:
            continue
        stack, comp = [start], []
        seen[start] = True
        while stack:
            node = stack.pop()
            comp.append(node)
            for nbr in np.nonzero(adjacency[node] & ~seen)[0]:
                seen[nbr] = True
                stack.append(int(nbr))
        components.append(sorted(comp))
    return sorted(components, key=lambda c: c[0])


def jaccard(a: set, b: set) -> float:
    """Jaccard index. Both empty -> 0.0 by recorded convention (II.5D)."""
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def dice(a: set, b: set) -> float:
    """Dice coefficient. Both empty -> 0.0 by recorded convention."""
    if not a and not b:
        return 0.0
    return 2 * len(a & b) / (len(a) + len(b))


def mst_edges(coords: np.ndarray) -> list[tuple[int, int, float]]:
    """Euclidean minimum spanning tree edges as ``(i, j, weight)`` (II.9).

    Deterministic Prim's algorithm; ties broken by the smallest node index, so
    the edge list is reproducible for any input ordering of identical geometry.
    """
    n = len(coords)
    if n < 2:
        return []
    dist = pairwise_distances(coords)
    in_tree = np.zeros(n, dtype=bool)
    in_tree[0] = True
    best = dist[0].copy()
    parent = np.zeros(n, dtype=int)
    edges: list[tuple[int, int, float]] = []
    for _ in range(n - 1):
        masked = np.where(in_tree, np.inf, best)
        j = int(np.argmin(masked))
        i = int(parent[j])
        edges.append((min(i, j), max(i, j), float(dist[i, j])))
        in_tree[j] = True
        closer = (~in_tree) & (dist[j] < best)
        best = np.where(closer, dist[j], best)
        parent = np.where(closer, j, parent)
    return sorted(edges, key=lambda e: (e[2], e[0], e[1]))


def centroid(coords: np.ndarray) -> np.ndarray:
    return np.asarray(coords, dtype=np.float64).mean(axis=0)


def hausdorff95(a: np.ndarray, b: np.ndarray) -> float:
    """Symmetric 95th-percentile Hausdorff distance (II.11 comparison metric)."""
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    d = cross_distances(a, b)
    return float(max(np.percentile(d.min(axis=1), 95),
                     np.percentile(d.min(axis=0), 95)))


def mcc_from_confusion(tp: int, fp: int, tn: int, fn: int,
                       undefined_value: float = 0.0) -> float:
    """Matthews correlation coefficient.

    A zero denominator yields ``undefined_value`` (config: ``loo_mcc.undefined_mcc_value``,
    default 0.0) by the recorded convention "indistinguishable from chance".
    """
    denom = float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denom <= 0:
        return undefined_value
    return float((tp * tn - fp * fn) / np.sqrt(denom))


def classification_metrics(tp: int, fp: int, tn: int, fn: int,
                           undefined_mcc: float = 0.0) -> dict[str, float]:
    """The full ¶49 metric family from one confusion matrix. Single implementation."""
    def _safe(num: float, den: float) -> float:
        return float(num / den) if den > 0 else float("nan")

    sens = _safe(tp, tp + fn)
    spec = _safe(tn, tn + fp)
    ppv = _safe(tp, tp + fp)
    npv = _safe(tn, tn + fn)
    f1 = _safe(2 * tp, 2 * tp + fp + fn)
    return {
        "tp": int(tp), "fp": int(fp), "tn": int(tn), "fn": int(fn),
        "mcc": mcc_from_confusion(tp, fp, tn, fn, undefined_mcc),
        "sens": sens, "spec": spec,
        "acc": _safe(tp + tn, tp + tn + fp + fn),
        "ppv": ppv, "npv": npv, "f1": f1,
        "balacc": float(np.nanmean([sens, spec])),
    }
