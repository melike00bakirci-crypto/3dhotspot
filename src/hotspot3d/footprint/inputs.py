"""Upstream loading and BLOCKING input validation (agent §6).

Everything this stage consumes is read-only and hash-recorded. Two structural
guarantees are implemented here rather than merely promised:

* ``10_ANNOTATION/`` is never opened — :func:`_guard_path` refuses the path, and
  reads go through :meth:`RunContext.assert_may_read` which excludes it from this
  agent's read scope entirely;
* ``r_hot`` is loaded into a wrapper that raises if it is ever compared with, or
  substituted for, a footprint radius (F10).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..orchestration.contracts import FORBIDDEN_DOWNSTREAM_COLUMNS, Handoff, QC_FAIL
from ..utils.errors import BlockedError, LeakageError
from ..utils.hashing import sha256_file
from ..utils.io import read_tsv, read_tsv_columns
from ..utils.runctx import RunContext
from .api import CenterSet, Universe

AGENT = "footprint-robustness"
FORBIDDEN_STAGE = "10_ANNOTATION"

#: Canonical upstream files (agent §5) and the stage directories they may live in.
UPSTREAM_FILES: dict[str, tuple[str, ...]] = {
    "significant_hotspot_centers.tsv": ("06_FINAL_HOTSPOTS",),
    "hotspot_covered_residues.tsv": ("06_FINAL_HOTSPOTS",),
    "hotspot_classified_variants.tsv": ("06_FINAL_HOTSPOTS",),
    "residue_coordinates.tsv": ("03_STRUCTURE_QC",),
    "positional_universe.tsv": ("03_STRUCTURE_QC",),
    "classified_cohort.tsv": ("02_CLINVAR", "03_STRUCTURE_QC"),
    "structures/structure_with_plddt.cif": ("03_STRUCTURE_QC",),
}

#: Accepted column aliases, so a purely cosmetic upstream naming choice does not
#: block a run. Absence of the underlying quantity still blocks.
_CENTER_ID_ALIASES = ("center_residue_index", "residue_index", "center_residue")
_REGION_ALIASES = ("region_id", "hotspot_region_id", "region")
_CARRIES_ALIASES = ("carries_clinvar_variant", "carries_variant")
_CLASS_ALIASES = ("class", "residue_class", "primary_class")

_PLP = "PLP"
_BLB = "BLB"
_CLASS_MAP = {
    "PLP": _PLP, "P/LP": _PLP, "PATHOGENIC": _PLP, "P": _PLP, "LP": _PLP,
    "PATHOGENIC/LIKELY PATHOGENIC": _PLP, "LIKELY PATHOGENIC": _PLP,
    "BLB": _BLB, "B/LB": _BLB, "BENIGN": _BLB, "B": _BLB, "LB": _BLB,
    "BENIGN/LIKELY BENIGN": _BLB, "LIKELY BENIGN": _BLB,
}


class ReadOnlyHotspotRadius:
    """``r_hot`` as recorded context that cannot leak into the ``r_fp`` decision.

    F10 says ``r_fp`` is neither set equal to nor bounded below by ``r_hot``. This
    wrapper makes the prohibited operations raise instead of merely being frowned
    upon: the value can be recorded and used for the explicitly permitted
    ``_hs`` coverage query, and nothing else.
    """

    __slots__ = ("_value",)

    def __init__(self, value: float):
        object.__setattr__(self, "_value", float(value))

    @property
    def value_for_coverage_query(self) -> float:
        """The ONLY sanctioned use: residues within ``r_hot`` of ``S_iter`` (¶49)."""
        return self._value

    def _blocked(self, *_args, **_kwargs):
        raise BlockedError(
            "F10 violation: r_hot was coerced into arithmetic or comparison in a "
            "footprint context. r_hot and r_fp are separate analytical quantities; "
            "r_fp is never set equal to, bounded by, or compared against r_hot. The "
            "only sanctioned use is value_for_coverage_query."
        )

    __float__ = _blocked
    __int__ = _blocked
    __lt__ = _blocked
    __le__ = _blocked
    __gt__ = _blocked
    __ge__ = _blocked
    __add__ = _blocked
    __sub__ = _blocked
    __mul__ = _blocked
    __truediv__ = _blocked

    def __repr__(self) -> str:
        return f"ReadOnlyHotspotRadius(recorded={self._value!r})"


@dataclass
class UpstreamInputs:
    centers: CenterSet
    universe: Universe
    center_regions: dict[int, str]
    center_carries_variant: dict[int, bool]
    cohort_plp: frozenset
    cohort_blb: frozenset
    r_hot: ReadOnlyHotspotRadius
    paths: dict[str, Path]
    hashes: dict[str, str]
    n_significant_centers: int
    residue_names: dict[int, str]


def _guard_path(path: Path) -> Path:
    if FORBIDDEN_STAGE in Path(path).parts:
        raise LeakageError(
            f"INFORMATION BARRIER: {AGENT} attempted to open {path}. "
            f"{FORBIDDEN_STAGE}/ is never read in either phase — biological and "
            f"mechanism information plays no part in any footprint or robustness "
            f"decision."
        )
    return path


def resolve(ctx: RunContext, name: str) -> Path:
    """Locate a canonical upstream file inside this agent's read scope."""
    for stage in UPSTREAM_FILES[name]:
        ctx.assert_may_read(AGENT, stage)
        candidate = _guard_path(ctx.full_results / stage / name)
        if candidate.is_file():
            return candidate
    raise BlockedError(
        f"BLOCKED — required upstream input {name!r} not found in "
        f"{list(UPSTREAM_FILES[name])}. Missing upstream artifacts are never "
        f"substituted, defaulted or reconstructed."
    )


def assert_stage_dirs_empty(ctx: RunContext, stages: tuple[str, ...]) -> None:
    """§6.5 — the owned stage directories must be empty for this RUN_ID."""
    for stage in stages:
        ctx.assert_may_write(AGENT, stage)
        existing = [p for p in (ctx.full_results / stage).rglob("*") if p.is_file()]
        if existing:
            raise BlockedError(
                f"BLOCKED — {stage}/ is not empty for RUN_ID {ctx.run_id}: "
                f"{[str(p.name) for p in existing[:5]]}. A stage never overwrites a "
                f"previous run; the Lead declares a new RUN_ID."
            )


def load(ctx: RunContext, upstream: Handoff) -> UpstreamInputs:
    """Validate and load every upstream artifact. Every failure here is BLOCKING."""
    # 1-2. handoff integrity, QC verdict, non-empty center set
    upstream.verify(ctx)
    if upstream.qc_status == QC_FAIL:
        raise BlockedError("BLOCKED — handoff_02 reports qc_status=FAIL.")
    n_centers = int(upstream.payload.get("n_significant_centers", -1))
    if n_centers < 0:
        raise BlockedError("BLOCKED — handoff_02 does not declare "
                           "n_significant_centers.")
    if n_centers == 0:
        raise BlockedError(
            "BLOCKED — handoff_02 reports n_significant_centers = 0. The chain "
            "already ended upstream in a valid scientific negative; there is no "
            "footprint to construct and nothing here to rescue."
        )
    if "hotspot_radius" not in upstream.payload:
        raise BlockedError("BLOCKED — handoff_02 does not declare hotspot_radius.")

    paths = {name: resolve(ctx, name) for name in UPSTREAM_FILES}
    hashes = {name: sha256_file(path) for name, path in paths.items()}

    # 3-4. coordinates, universe, cohort
    universe, coord_by_id, residue_names = _load_universe(paths)
    centers, regions, carries = _load_centers(paths, universe, coord_by_id)
    plp, blb = _load_cohort(paths, universe)

    return UpstreamInputs(
        centers=centers, universe=universe, center_regions=regions,
        center_carries_variant=carries, cohort_plp=plp, cohort_blb=blb,
        r_hot=ReadOnlyHotspotRadius(upstream.payload["hotspot_radius"]),
        paths=paths, hashes=hashes, n_significant_centers=n_centers,
        residue_names=residue_names,
    )


def _load_universe(paths: dict[str, Path]
                   ) -> tuple[Universe, dict[int, np.ndarray], dict[int, str]]:
    coords_rows = read_tsv_columns(
        paths["residue_coordinates.tsv"],
        ["residue_index", "x_ca", "y_ca", "z_ca", "ca_usable",
         "residue_aa3", "residue_aa", "aa3"])
    universe_rows = read_tsv_columns(paths["positional_universe.tsv"],
                                     ["residue_index"])

    coord_by_id: dict[int, np.ndarray] = {}
    usable: dict[int, bool] = {}
    names: dict[int, str] = {}
    for row in coords_rows:
        rid = row.get("residue_index")
        if rid is None:
            continue
        xyz = [row.get("x_ca"), row.get("y_ca"), row.get("z_ca")]
        if any(v is None for v in xyz):
            continue
        coord_by_id[int(rid)] = np.array([float(v) for v in xyz], dtype=np.float64)
        flag = row.get("ca_usable")
        usable[int(rid)] = True if flag is None else bool(flag)
        aa = row.get("residue_aa3") or row.get("residue_aa") or row.get("aa3")
        names[int(rid)] = str(aa)[:3].upper() if aa else "UNK"

    ids = sorted({int(r["residue_index"]) for r in universe_rows
                  if r.get("residue_index") is not None})
    if not ids:
        raise BlockedError("BLOCKED — positional_universe.tsv declares no residues; "
                           "|U_struct| = 0 and coverage would be undefined.")
    missing = [i for i in ids if i not in coord_by_id]
    if missing:
        raise BlockedError(
            f"BLOCKED — {len(missing)} residues of U_struct have no CA coordinate "
            f"(first: {missing[:10]}). U_struct is by definition the set of residues "
            f"with a usable CA."
        )
    unusable = [i for i in ids if not usable.get(i, True)]
    if unusable:
        raise BlockedError(
            f"BLOCKED — U_struct contains residues flagged ca_usable=FALSE "
            f"(first: {unusable[:10]}); the upstream universe is inconsistent."
        )
    coords = np.vstack([coord_by_id[i] for i in ids])
    if not np.all(np.isfinite(coords)):
        raise BlockedError("BLOCKED — non-finite CA coordinate in U_struct.")
    return (Universe(ids=tuple(ids), coords=coords), coord_by_id,
            {i: names.get(i, "UNK") for i in ids})


def _load_centers(paths: dict[str, Path], universe: Universe,
                  coord_by_id: dict[int, np.ndarray]
                  ) -> tuple[CenterSet, dict[int, str], dict[int, bool]]:
    rows = read_tsv(paths["significant_hotspot_centers.tsv"],
                    forbidden_columns=FORBIDDEN_DOWNSTREAM_COLUMNS)
    if not rows:
        raise BlockedError(
            "BLOCKED — significant_hotspot_centers.tsv is empty while handoff_02 "
            "declares a non-zero center count."
        )
    id_col = _first_present(rows[0], _CENTER_ID_ALIASES)
    if id_col is None:
        raise BlockedError(
            f"BLOCKED — significant_hotspot_centers.tsv has no center identifier "
            f"column (accepted: {list(_CENTER_ID_ALIASES)})."
        )
    carries_col = _first_present(rows[0], _CARRIES_ALIASES)
    if carries_col is None:
        raise BlockedError(
            "BLOCKED — significant_hotspot_centers.tsv does not carry "
            "'carries_clinvar_variant'. QC-D17 requires it populated: a center is a "
            "geometric position and need not carry a variant, and the perturbation "
            "universe must state which do. This column is produced by "
            "hotspot-statistics and is never derived or defaulted here."
        )
    region_col = _first_present(rows[0], _REGION_ALIASES)

    ids, regions, carries = [], {}, {}
    for row in rows:
        rid = row.get(id_col)
        if rid is None:
            raise BlockedError("BLOCKED — a significant center has a null identifier.")
        rid = int(rid)
        ids.append(rid)
        regions[rid] = str(row.get(region_col)) if region_col else "NA"
        flag = row.get(carries_col)
        if flag is None:
            raise BlockedError(
                f"BLOCKED — carries_clinvar_variant is null for center {rid}."
            )
        carries[rid] = bool(flag)

    if len(set(ids)) != len(ids):
        raise BlockedError("BLOCKED — duplicate significant hotspot centers.")
    ids = sorted(ids)

    index_of = {rid: i for i, rid in enumerate(universe.ids)}
    outside = [i for i in ids if i not in index_of]
    if outside:
        raise BlockedError(
            f"BLOCKED — significant centers {outside[:10]} are not in U_struct. "
            f"Every center must have a usable CA (agent §6.3)."
        )
    coords = np.vstack([coord_by_id[i] for i in ids])
    if not np.all(np.isfinite(coords)):
        raise BlockedError("BLOCKED — a significant center has a non-finite CA.")

    centers = CenterSet(ids=tuple(ids), coords=coords,
                        universe_rows=tuple(index_of[i] for i in ids))
    return centers, regions, carries


def _load_cohort(paths: dict[str, Path],
                 universe: Universe) -> tuple[frozenset, frozenset]:
    """The classified cohort ``L``. Read-only, never modified, never re-derived."""
    rows = read_tsv(paths["classified_cohort.tsv"])
    if not rows:
        raise BlockedError("BLOCKED — classified_cohort.tsv is empty; the cohort L "
                           "is required for the ¶49 classification metrics.")
    class_col = _first_present(rows[0], _CLASS_ALIASES)
    if class_col is None:
        raise BlockedError(
            f"BLOCKED — classified_cohort.tsv has no class column "
            f"(accepted: {list(_CLASS_ALIASES)})."
        )
    if class_col in FORBIDDEN_DOWNSTREAM_COLUMNS:       # defensive, never reachable
        raise LeakageError(f"forbidden column {class_col!r} requested")

    universe_ids = set(universe.ids)
    plp, blb = set(), set()
    for row in rows:
        rid = row.get("residue_index")
        raw = row.get(class_col)
        if rid is None or raw is None:
            continue
        key = str(raw).strip().upper()
        if key not in _CLASS_MAP:
            raise BlockedError(
                f"BLOCKED — classified_cohort.tsv contains class {raw!r}, which is "
                f"not one of the two primary binary classes. Class labels are "
                f"upstream-owned and are never reinterpreted here."
            )
        target = plp if _CLASS_MAP[key] == _PLP else blb
        if int(rid) in universe_ids:
            target.add(int(rid))
    if not plp and not blb:
        raise BlockedError("BLOCKED — the classified cohort L is empty within "
                           "U_struct.")
    overlap = plp & blb
    if overlap:
        raise BlockedError(
            f"BLOCKED — residues {sorted(overlap)[:10]} appear in both cohort "
            f"classes. Residue class conflicts are resolved upstream (F3), never here."
        )
    return frozenset(plp), frozenset(blb)


def _first_present(row: dict, aliases: tuple[str, ...]) -> str | None:
    for name in aliases:
        if name in row:
            return name
    return None
