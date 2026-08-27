"""Stage B input validation and the information barriers, enforced in code.

Two barriers are implemented here rather than merely promised in prose:

**The column allowlist.** Only ``stage_b_permitted_columns`` enter the primary path.
Any ``forbidden_downstream`` column present in a table aborts the stage as a
LEAKAGE EVENT. Columns that are neither permitted nor forbidden (Stage A also emits
``aa``, ``plddt_band`` and the CB / side-chain columns that F1 marks unused) are
**discarded at load time** and the discard is logged, so the set of columns that
actually reached the analysis is auditable rather than assumed.

**No downstream reads.** ``07_FOOTPRINT_RADIUS``, ``08_FINAL_FOOTPRINT``,
``09_ROBUSTNESS`` and ``10_ANNOTATION`` are never opened. Every read goes through
:class:`ReadGuard`, which refuses those paths and records everything it did read, so
the final assertion is made against a real ledger.

Every failure in this module is BLOCKING. Nothing is ever defaulted.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..orchestration.contracts import (
    FORBIDDEN_DOWNSTREAM_COLUMNS,
    STAGE_B_PERMITTED_COLUMNS,
    Handoff,
)
from ..utils.errors import BlockedError, LeakageError
from ..utils.io import read_tsv
from ..utils.runctx import RunContext

FORBIDDEN_STAGE_DIRS = ("07_FOOTPRINT_RADIUS", "08_FINAL_FOOTPRINT", "09_ROBUSTNESS",
                        "10_ANNOTATION")
OWNED_STAGE_DIRS = ("04_GLOBAL_CLUSTERING", "05_HOTSPOT_RADIUS", "06_FINAL_HOTSPOTS",
                    "11_SENSITIVITY")
UPSTREAM_STAGE_DIRS = ("01_INPUT_RAW", "02_CLINVAR", "03_STRUCTURE_QC")


class ReadGuard:
    """Records every path Stage B opens and refuses anything downstream of Stage B."""

    def __init__(self) -> None:
        self.paths: list[str] = []

    def check(self, path: str | Path) -> Path:
        path = Path(path)
        parts = set(path.parts)
        hit = sorted(parts & set(FORBIDDEN_STAGE_DIRS))
        if hit:
            raise LeakageError(
                f"INFORMATION BARRIER violation: Stage B attempted to read {hit} at "
                f"{path}. Stage C/D/E output is never an input to primary discovery.")
        self.paths.append(str(path))
        return path

    def assert_clean(self) -> None:
        offending = [p for p in self.paths
                     if any(d in Path(p).parts for d in FORBIDDEN_STAGE_DIRS)]
        if offending:
            raise LeakageError(
                f"INFORMATION BARRIER violation detected in the read ledger: {offending}")


@dataclass
class StageBInputs:
    """The only data that reaches the primary path, plus its provenance."""

    universe_index: np.ndarray            # residue indices of U_struct, ascending
    universe_coords: np.ndarray           # (M, 3) CA coordinates
    universe_plddt: np.ndarray            # recorded, NEVER a primary filter (F2)
    # --- [NEW — Workflow v2 §2] the candidate-CENTER universe U_center ---------
    # A residue with pLDDT < plddt.center_universe_min_plddt is coordinate noise as
    # a sphere CENTER and is excluded from U_center. It is UNCHANGED in every other
    # role: still in U_struct (the positional reference, F9) and still in L if it
    # carries a variant (F2 is unaffected). ``center_universe_mask`` is aligned with
    # ``universe_index``/``universe_coords`` — True where the residue is eligible.
    center_universe_mask: np.ndarray
    center_universe_min_plddt: float
    labeled_index: np.ndarray             # residue indices of L, ascending
    labeled_positions: np.ndarray         # positions of L within U_struct
    labeled_coords: np.ndarray
    y: np.ndarray                         # 1 = P/LP, 0 = B/LB
    allowlist_used: list[str]
    load_report: dict = field(default_factory=dict)
    paths: dict = field(default_factory=dict)

    @property
    def M(self) -> int:
        return len(self.universe_index)

    @property
    def N(self) -> int:
        return len(self.y)

    @property
    def N_P(self) -> int:
        return int(self.y.sum())

    @property
    def N_B(self) -> int:
        return int(self.N - self.N_P)

    @property
    def center_index(self) -> np.ndarray:
        """Residue indices of U_center — ONLY the eligible candidate centers."""
        return self.universe_index[self.center_universe_mask]

    @property
    def center_coords(self) -> np.ndarray:
        return self.universe_coords[self.center_universe_mask]

    @property
    def n_center_universe(self) -> int:
        return int(self.center_universe_mask.sum())


def read_allowlisted(path: Path, required: list[str], guard: ReadGuard,
                     allowlist: list[str] = STAGE_B_PERMITTED_COLUMNS,
                     forbidden: list[str] = FORBIDDEN_DOWNSTREAM_COLUMNS
                     ) -> tuple[list[dict], dict]:
    """Load a Stage A table under the Stage B column allowlist.

    A forbidden column raises :class:`LeakageError` (via ``utils.io.read_tsv``). Every
    column outside the allowlist is dropped before the data reaches the analysis, and
    the projection is asserted in code — a non-allowlisted key can never appear in a
    returned row.
    """
    guard.check(path)
    rows = read_tsv(path, forbidden_columns=forbidden)
    header = sorted(rows[0].keys()) if rows else []
    missing = [c for c in required if c not in header]
    if missing:
        raise BlockedError(
            f"BLOCKED — {Path(path).name} is missing required column(s) {missing}. "
            f"Required columns are never reconstructed or defaulted.")

    loaded = [c for c in header if c in allowlist]
    discarded = [c for c in header if c not in allowlist and c != "schema_version"]
    projected = [{c: r.get(c) for c in loaded} for r in rows]
    for row in projected:                      # asserted, not assumed
        extra = set(row) - set(allowlist)
        if extra:
            raise LeakageError(
                f"LEAKAGE EVENT: column(s) {sorted(extra)} survived the allowlist "
                f"projection of {Path(path).name}. Stage aborted.")

    return projected, {
        "file": str(path), "n_rows": len(projected),
        "header_seen": header, "columns_loaded": loaded,
        "columns_discarded_not_allowlisted": discarded,
        "allowlist_asserted": sorted(allowlist),
        "forbidden_checked": sorted(forbidden),
    }


def assert_owned_dirs_empty(ctx: RunContext, guard: ReadGuard,
                            stages: tuple[str, ...] = OWNED_STAGE_DIRS) -> None:
    """Stage B never overwrites a previous Stage B result for the same RUN_ID."""
    for stage in stages:
        path = ctx.full_results / stage
        if not path.exists():
            continue
        existing = [p for p in path.rglob("*") if p.is_file()]
        if existing:
            raise BlockedError(
                f"BLOCKED — {stage} is not empty for run_id {ctx.run_id} "
                f"({len(existing)} file(s), e.g. {existing[0].name}). Stage B never "
                f"overwrites a previous result; the Lead declares a new RUN_ID.")


def assert_no_downstream_dirs_touched(ctx: RunContext, guard: ReadGuard) -> dict:
    """Final assertion (agent section 6.7): 07_*, 08_*, 09_*, 10_* were never opened."""
    guard.assert_clean()
    return {
        "forbidden_stage_dirs": list(FORBIDDEN_STAGE_DIRS),
        "n_paths_read": len(guard.paths),
        "paths_read": sorted(set(guard.paths)),
        "downstream_reads": 0,
        "assertion": "no path under 07_*, 08_*, 09_* or 10_* was opened by Stage B",
    }


def load_inputs(ctx: RunContext, upstream: Handoff, guard: ReadGuard) -> StageBInputs:
    """Verify ``handoff_01`` and load ``U_struct``, ``L`` and the coordinates."""
    qc_dir = ctx.full_results / "03_STRUCTURE_QC"

    coords_rows, coords_report = read_allowlisted(
        qc_dir / "residue_coordinates.tsv",
        ["residue_index", "x_ca", "y_ca", "z_ca", "ca_usable"], guard)
    universe_rows, universe_report = read_allowlisted(
        qc_dir / "positional_universe.tsv", ["residue_index"], guard)
    cohort_rows, cohort_report = read_allowlisted(
        qc_dir / "classified_cohort.tsv", ["residue_index", "class"], guard)

    coord_by_index = {int(r["residue_index"]): r for r in coords_rows}
    universe_index = np.array(sorted(int(r["residue_index"]) for r in universe_rows),
                              dtype=np.int64)
    if len(universe_index) == 0:
        raise BlockedError("BLOCKED — U_struct is empty; there is nothing to test.")
    if len(universe_index) != len(universe_rows):
        raise BlockedError("BLOCKED — positional_universe.tsv contains duplicate residues.")

    missing_coords = [int(i) for i in universe_index if i not in coord_by_index]
    if missing_coords:
        raise BlockedError(
            f"BLOCKED — {len(missing_coords)} residue(s) of U_struct have no row in "
            f"residue_coordinates.tsv (e.g. {missing_coords[:5]}).")

    unusable = [int(i) for i in universe_index
                if coord_by_index[int(i)].get("ca_usable") is not True]
    if unusable:
        raise BlockedError(
            f"BLOCKED — {len(unusable)} residue(s) of U_struct are flagged "
            f"ca_usable != TRUE (e.g. {unusable[:5]}). U_struct is defined as the set "
            f"of residues WITH a usable CA.")

    universe_coords = np.array(
        [[float(coord_by_index[int(i)]["x_ca"]), float(coord_by_index[int(i)]["y_ca"]),
          float(coord_by_index[int(i)]["z_ca"])] for i in universe_index],
        dtype=np.float64)
    if not np.isfinite(universe_coords).all():
        raise BlockedError("BLOCKED — non-finite CA coordinate in U_struct.")
    universe_plddt = np.array(
        [_as_float(coord_by_index[int(i)].get("plddt")) for i in universe_index],
        dtype=np.float64)

    position_of = {int(idx): pos for pos, idx in enumerate(universe_index)}
    labeled_index, y = [], []
    for row in sorted(cohort_rows, key=lambda r: int(r["residue_index"])):
        idx = int(row["residue_index"])
        label = str(row["class"]).strip().upper()
        if label not in ("PLP", "BLB"):
            raise BlockedError(
                f"BLOCKED — classified_cohort.tsv residue {idx} carries class "
                f"{row['class']!r}; Stage B accepts only the binary classes PLP and "
                f"BLB. Every other disposition stays outside the primary comparison (F12).")
        if idx not in position_of:
            raise BlockedError(
                f"BLOCKED — labelled residue {idx} is not in U_struct. Every residue of "
                f"L must be present in U_struct with a usable CA.")
        labeled_index.append(idx)
        y.append(1 if label == "PLP" else 0)

    if len(labeled_index) != len(set(labeled_index)):
        raise BlockedError("BLOCKED — classified_cohort.tsv contains duplicate residues.")

    labeled_index_arr = np.array(labeled_index, dtype=np.int64)
    y_arr = np.array(y, dtype=np.int64)
    if y_arr.sum() == 0 or (len(y_arr) - y_arr.sum()) == 0:
        raise BlockedError(
            f"BLOCKED — a class is empty (N_P = {int(y_arr.sum())}, "
            f"N_B = {int(len(y_arr) - y_arr.sum())}). N_P = 0 leaves the primary "
            f"positional null nothing to place; N_B = 0 leaves the LOO prevalence and "
            f"the secondary label-permutation null undefined. v2 §1 also requires N_B "
            f"to be reported and weighed before any negative result is issued.")

    positions = np.array([position_of[i] for i in labeled_index], dtype=np.int64)

    # --- v2 §2: U_center = {r in U_struct : plddt(r) >= center_universe_min_plddt}
    # NaN pLDDT (never observed in practice — AlphaFold always scores every residue)
    # compares False under `>=` and is therefore conservatively EXCLUDED: an
    # unscored residue is exactly the "coordinate noise" case this rule targets.
    min_plddt = float(ctx.config.get("plddt.center_universe_min_plddt"))
    excludes_low_confidence = bool(
        ctx.config.get("plddt.center_universe_excludes_low_confidence"))
    center_universe_mask = (universe_plddt >= min_plddt if excludes_low_confidence
                            else np.ones(len(universe_index), dtype=bool))

    return StageBInputs(
        universe_index=universe_index, universe_coords=universe_coords,
        universe_plddt=universe_plddt,
        center_universe_mask=center_universe_mask,
        center_universe_min_plddt=min_plddt,
        labeled_index=labeled_index_arr,
        labeled_positions=positions, labeled_coords=universe_coords[positions],
        y=y_arr, allowlist_used=sorted(STAGE_B_PERMITTED_COLUMNS),
        load_report={"residue_coordinates": coords_report,
                     "positional_universe": universe_report,
                     "classified_cohort": cohort_report},
        paths={"residue_coordinates": str(qc_dir / "residue_coordinates.tsv"),
               "positional_universe": str(qc_dir / "positional_universe.tsv"),
               "classified_cohort": str(qc_dir / "classified_cohort.tsv")},
    )


def _as_float(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def load_review_channel(ctx: RunContext, upstream: Handoff, guard: ReadGuard) -> dict:
    """The SEPARATE, explicitly non-redefining channel for ClinVar review stars.

    Review metadata is a ``forbidden_downstream`` column for the Stage B primary path
    (F12: stars are never an inclusion criterion). The >= 1* / >= 2* sensitivity
    analyses may read it only through a channel that ``handoff_01`` declares
    explicitly. If no channel is declared, the analysis is NOT_EVALUABLE — the
    columns are never fetched from the primary tables as a workaround.
    """
    declaration = (upstream.payload.get("sensitivity_channel")
                   or upstream.payload.get("sensitivity_review_channel"))
    if not declaration or declaration.get("declared") is False:
        return {"available": False,
                "reason": ("handoff_01 declares no sensitivity_channel; the "
                           "review-status sensitivity analysis is NOT_EVALUABLE. Star "
                           "columns are forbidden in the Stage B primary path and are "
                           "never read from the primary tables as a substitute."),
                "stars": {}}
    if declaration.get("may_redefine_primary"):
        raise LeakageError(
            "LEAKAGE EVENT: handoff_01 declares a sensitivity channel with "
            "may_redefine_primary=TRUE. The review-status channel is non-redefining by "
            "construction; Stage B refuses to open it on those terms.")

    rel = declaration.get("file")
    declared_columns = declaration.get("columns") or []
    key = declaration.get("key_column", "residue_index")
    column = declaration.get("max_star_column") or (
        "max_star" if "max_star" in declared_columns or not declared_columns
        else declared_columns[-1])
    path = ctx.run_root / rel if rel else None
    if path is None or not Path(path).is_file():
        return {"available": False,
                "reason": f"declared review channel file not found: {rel}", "stars": {}}

    guard.check(path)
    rows = read_tsv(path)                       # separate channel: NO primary allowlist
    if rows and (key not in rows[0] or column not in rows[0]):
        return {"available": False,
                "reason": (f"declared review channel {rel} lacks column(s) "
                           f"{key!r}/{column!r}"), "stars": {}}
    stars = {}
    for r in rows:
        try:
            stars[int(r[key])] = float(r[column])
        except (TypeError, ValueError):
            continue
    return {"available": True, "reason": "NA", "stars": stars,
            "channel": {"file": rel, "key_column": key, "max_star_column": column,
                        "declared_columns": declared_columns},
            "note": ("Read through the declared non-redefining channel only. This data "
                     "never re-enters the primary path.")}
