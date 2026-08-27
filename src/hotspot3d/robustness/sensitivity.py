"""¶48 center-loss sensitivity — which centers the footprint actually depends on.

The influence of a center is the **mean drop in covered-residue Jaccard across all
evaluated iterations that removed it**, with the number of contributing iterations
reported alongside. The count matters: an influence estimated from three
iterations and one estimated from three thousand are not the same claim, and a
ranked table without ``n`` invites reading noise as structure.

Baseline Jaccard is 1.0 by definition (``FP_original`` against itself), so the drop
is ``1 - J_iter``. Failed iterations contribute ``J = 0``, hence a drop of 1.0, and
stay in the denominator.
"""
from __future__ import annotations

import numpy as np

CENTER_SENSITIVITY_COLUMNS = [
    "rank", "center_residue_index", "region_id", "carries_clinvar_variant",
    "n_iterations_removing", "n_iterations_retaining",
    "mean_jaccard_when_removed", "mean_jaccard_drop", "sd_jaccard_drop",
    "median_jaccard_when_removed", "min_jaccard_when_removed",
    "n_failed_iterations_removing",
    "center_recurrence_rate", "n_recurrent_iterations",
]

# The retained-center form of this quantity — "is a RETAINED center covered by
# FP_iter?" — is identically 1.0 by construction, because a retained center is an
# element of S_iter, so its distance to S_iter is 0 and it lies inside FP_iter at
# every radius; see QC-F5, which asserts exactly that containment and BLOCKS as a
# computational fault if it ever fails. It was proposed, measured (24/24 retained
# vs 8/16 removed on the four-center fixture) and rejected: it would have reported
# the most fragile center in the run as perfectly recurrent. Do not reinstate it.
#
#: FROZEN semantics. Recorded here, in stage_d_report.md and in handoff_04 so that
#: nobody later "fixes" this back to the vacuous form.
RECURRENCE_DEFINITION = (
    "center_recurrence_rate(c) = #(iterations that REMOVED c and whose reconstructed "
    "footprint still covers c) / #(iterations that removed c). This is the para-51 "
    "quantity 'the frequency with which individual hotspots reappear', ported to the "
    "frozen footprint-only design: the hotspot reappears in the reconstruction even "
    "though its own center is gone. Failed iterations stay in the denominator and "
    "count as non-recurrence. A zero denominator gives NULL, never 0.0. "
    "NOTE FOR ANY FUTURE READER: the complementary form — 'a RETAINED center is "
    "covered by FP_iter' — is vacuous and must not be substituted. A retained center "
    "is an element of S_iter, so its distance to S_iter is 0 and it lies inside "
    "FP_iter at every radius; QC-F5 asserts exactly that containment and would BLOCK "
    "as a computational fault if it ever failed. That form is identically 1.0 except "
    "where an iteration fails, at which point it measures exposure to reconstruction "
    "failure — which n_failed_iterations_removing already carries. A constant that "
    "reads as strong evidence while carrying none is the worst kind of reported "
    "number, so it is not emitted under any name."
)


def center_influence(results, center_ids: tuple[int, ...],
                     regions: dict[int, str],
                     carries: dict[int, bool]) -> list[dict]:
    per_center: dict[int, list[float]] = {int(c): [] for c in center_ids}
    failures: dict[int, int] = {int(c): 0 for c in center_ids}
    recurrent: dict[int, int] = {int(c): 0 for c in center_ids}
    retained: dict[int, int] = {int(c): 0 for c in center_ids}

    for result in results:
        jaccard = result.row.get("jaccard_residues")
        value = 0.0 if jaccard is None else float(jaccard)
        removed = {int(c) for c in result.removed}
        for center in center_ids:
            center = int(center)
            if center in removed:
                per_center[center].append(value)
                if result.failed:
                    failures[center] += 1
                # ¶51 recurrence: does the footprint reappear over this position
                # even though the center that anchored it is gone?
                if center in result.residues:
                    recurrent[center] += 1
            else:
                retained[center] += 1

    def _rate(numerator: int, denominator: int) -> float | None:
        """NULL on a zero denominator — never 0.0.

        A center that no evaluated iteration happened to remove has no recurrence
        evidence at all; reporting 0.0 would present "never measured" as "never
        recurred", which is the strongest possible claim from no data.
        """
        return float(numerator / denominator) if denominator else None

    rows = []
    for center in center_ids:
        values = np.asarray(per_center[int(center)], dtype=np.float64)
        common = {
            "center_residue_index": int(center),
            "region_id": regions.get(int(center), "NA"),
            "carries_clinvar_variant": bool(carries.get(int(center), False)),
            "n_iterations_retaining": retained[int(center)],
            "center_recurrence_rate": _rate(recurrent[int(center)], values.size),
            "n_recurrent_iterations": recurrent[int(center)],
        }
        if values.size:
            drops = 1.0 - values
            rows.append({
                **common,
                "n_iterations_removing": int(values.size),
                "mean_jaccard_when_removed": float(values.mean()),
                "mean_jaccard_drop": float(drops.mean()),
                "sd_jaccard_drop": float(drops.std(ddof=1)) if values.size > 1 else 0.0,
                "median_jaccard_when_removed": float(np.median(values)),
                "min_jaccard_when_removed": float(values.min()),
                "n_failed_iterations_removing": failures[int(center)],
            })
        else:
            rows.append({
                **common,
                "n_iterations_removing": 0,
                "mean_jaccard_when_removed": None,
                "mean_jaccard_drop": None,
                "sd_jaccard_drop": None,
                "median_jaccard_when_removed": None,
                "min_jaccard_when_removed": None,
                "n_failed_iterations_removing": 0,
            })

    rows.sort(key=lambda r: (-(r["mean_jaccard_drop"] if r["mean_jaccard_drop"]
                               is not None else -1.0),
                             r["center_residue_index"]))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows
