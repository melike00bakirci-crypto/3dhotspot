"""Stage A — ClinVar ingestion, classification and cohort construction.

Public surface (agent ``data-structure``)::

    from hotspot3d.data import run_stage_a
    handoff = run_stage_a(ctx, gene="BRAF", source=provider)

``source`` satisfies :class:`hotspot3d.data.sources.StageASource`. Synthetic
providers come from :mod:`hotspot3d.data.synthetic`.

Re-exports are resolved lazily (PEP 562). ``hotspot3d.structure`` reads the
shared amino-acid vocabulary from :mod:`hotspot3d.data.hgvs` while
:mod:`hotspot3d.data.stage` reads coordinates from ``hotspot3d.structure``;
importing either package eagerly here would close that loop.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:                                   # pragma: no cover
    from .sources import GeneResolution, SequenceSource, StageASource, VariantSource
    from .stage import run_stage_a

_LAZY = {
    "run_stage_a": ".stage",
    "StageASource": ".sources",
    "VariantSource": ".sources",
    "SequenceSource": ".sources",
    "GeneResolution": ".sources",
}

__all__ = sorted(_LAZY)


def __getattr__(name: str) -> Any:
    module = _LAZY.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from importlib import import_module

    return getattr(import_module(module, __name__), name)


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(_LAZY))
