"""¶49 classification metrics, reported as TWO distinctly labelled families.

The cohort ``L`` is **unchanged** in every iteration: perturbation removes
geometric center positions, never ClinVar records and never class labels. Only the
*predicted-positive region* moves.

Both families are always reported:

* **footprint family** — predicted positive iff the residue lies in ``FP_iter``;
* **``_hs`` family** — predicted positive iff the residue lies within ``r_hot`` of
  ``S_iter``.

The second family exists because F10 permits ``r_fp < r_hot``. When it does, a
footprint-only confusion matrix deflates MCC by construction, and reporting only
that number would understate classification performance for a purely geometric
reason. Computing residues within ``r_hot`` of ``S_iter`` is a **pure geometric
coverage query** — it is explicitly permitted and is emphatically not hotspot
re-testing: no permutation, no FDR, no significance is recomputed anywhere.
"""
from __future__ import annotations

import numpy as np

from ..footprint.api import CenterSet, Universe
from ..utils.geometry import classification_metrics, covered_mask

CLASSIFICATION_FIELDS = ["tp", "fp", "tn", "fn", "mcc", "sens", "spec", "acc",
                         "ppv", "npv", "f1", "balacc"]


def confusion(predicted_positive: frozenset, plp: frozenset,
              blb: frozenset) -> tuple[int, int, int, int]:
    tp = len(plp & predicted_positive)
    fp = len(blb & predicted_positive)
    fn = len(plp) - tp
    tn = len(blb) - fp
    return tp, fp, tn, fn


def metrics_for(predicted_positive: frozenset, plp: frozenset, blb: frozenset,
                undefined_mcc: float, suffix: str = "") -> dict:
    tp, fp, tn, fn = confusion(predicted_positive, plp, blb)
    values = classification_metrics(tp, fp, tn, fn, undefined_mcc=undefined_mcc)
    return {f"{k}{suffix}": values[k] for k in CLASSIFICATION_FIELDS}


def hotspot_coverage_family(iter_centers: CenterSet, universe: Universe,
                            r_hot_value: float, plp: frozenset, blb: frozenset,
                            undefined_mcc: float) -> dict:
    """Residues within ``r_hot`` of ``S_iter`` — a geometric query, nothing more."""
    if len(iter_centers) == 0:
        covered = frozenset()
    else:
        mask = covered_mask(iter_centers.coords, universe.coords, float(r_hot_value))
        covered = frozenset(int(universe.ids[i]) for i in np.nonzero(mask)[0])
    return metrics_for(covered, plp, blb, undefined_mcc, suffix="_hs")


def failed_row(suffix: str = "") -> dict:
    """A failed iteration keeps an explicit, empty confusion — never a blank."""
    return {f"{k}{suffix}": None for k in CLASSIFICATION_FIELDS}
