"""Annotation source interface, version pinning and the network gate.

Stage E is the only stage that legitimately needs external biological and literature
data, and only *after* discovery and robustness are complete and frozen (agent §14.5).
This module keeps that privilege narrow:

  * :class:`AnnotationSource` is the whole surface a source may expose. Stage E talks
    to nothing else, so swapping a live source for a mock changes no pipeline code.
  * Every source must declare a pinned release/version. An unversioned source may not
    be used — ``ANNOTATION_SOURCE_UNVERSIONED`` is BLOCKING (agent §6.4).
  * :class:`NetworkAnnotationSource` places an explicit ``allow_network`` gate in
    front of every retrieval. In this build the gate is the only implemented part;
    no retrieval is performed.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Protocol, runtime_checkable

from ..utils.errors import BlockedError
from .rubric import MechanismRecord

#: Methods that constitute the annotation surface, in report order.
SOURCE_METHODS = (
    "secondary_structure",
    "domains",
    "functional_sites",
    "conservation",
    "disease_associations",
    "literature_mechanisms",
)


@dataclass(frozen=True)
class SourceVersion:
    """A pinned external source. ``version`` is mandatory and never defaulted."""

    name: str
    version: str
    url: str = "NA"
    accession: str = "NA"
    retrieved_utc: str = "NA"

    def is_versioned(self) -> bool:
        text = (self.version or "").strip()
        return bool(text) and text.upper() not in {"NA", "NONE", "UNKNOWN", "LATEST"}

    def as_row(self) -> dict:
        return {
            "source_name": self.name,
            "version": self.version,
            "url": self.url,
            "accession": self.accession,
            "retrieved_utc": self.retrieved_utc,
        }


@runtime_checkable
class AnnotationSource(Protocol):
    """The complete annotation surface Stage E consumes.

    Implementations are free to be live, cached or synthetic; Stage E cannot tell
    the difference and must not try to.
    """

    def versions(self) -> dict[str, SourceVersion]:
        """Pinned release/version for every source this provider draws on."""
        ...

    def secondary_structure(self, structure: Any) -> dict[int, str]:
        """Residue index -> DSSP secondary-structure code."""
        ...

    def domains(self, acc: str) -> list[dict]:
        """InterPro/Pfam/UniProt domains: ``domain_id, name, start, end, source``."""
        ...

    def functional_sites(self, acc: str) -> list[dict]:
        """Functional regions, ligand-binding and protein-interaction sites."""
        ...

    def conservation(self, acc: str, residues: Iterable[int]) -> dict[int, dict]:
        """Residue index -> pinned dbNSFP metrics (GERP++ RS, phyloP100way).

        Genomic conservation mapped through the MANE transcript to a residue, which
        is what it is — recorded descriptively, never used to redefine a region.
        """
        ...

    def disease_associations(self, gene: str) -> list[dict]:
        """Previously reported disease associations, each with a citation."""
        ...

    def literature_mechanisms(self, gene: str) -> list[MechanismRecord]:
        """Curated functional-mechanism records, including excluded candidates.

        Excluded records (``included=False``) are returned deliberately: the search
        log must record per-record inclusion *and* exclusion decisions (agent §9.7).
        """
        ...


# --- validation -------------------------------------------------------------


def validate_source(source: Any) -> dict[str, SourceVersion]:
    """Check the surface and the version pinning. Every failure is BLOCKING."""
    missing = [m for m in ("versions",) + SOURCE_METHODS if not callable(getattr(source, m, None))]
    if missing:
        raise BlockedError(
            f"BLOCKED — annotation source {type(source).__name__} does not implement "
            f"the required surface: missing {missing}."
        )

    versions = source.versions()
    if not isinstance(versions, dict) or not versions:
        raise BlockedError(
            "BLOCKED — annotation source declared no versions. Every annotation "
            "source must carry a recorded release/version (agent §6.4)."
        )

    unversioned = sorted(
        name for name, ver in versions.items()
        if not isinstance(ver, SourceVersion) or not ver.is_versioned()
    )
    if unversioned:
        raise BlockedError(
            f"ANNOTATION_SOURCE_UNVERSIONED — BLOCKING. Unversioned annotation "
            f"source(s) {unversioned} may not be used. Pin a release and re-run; "
            f"an unpinned source makes the annotation irreproducible."
        )
    return versions


# --- the network gate -------------------------------------------------------

_GATE_MESSAGE = (
    "BLOCKED — network retrieval for {what!r} was requested but "
    "execution.allow_network is FALSE. Stage E holds network privileges only after "
    "discovery and robustness are frozen, and only when the run explicitly enables "
    "them. No request was made."
)


class NetworkAnnotationSource:
    """Base class for live sources. The gate is in front of every retrieval.

    Every method checks ``allow_network`` *before* anything else and raises
    :class:`BlockedError` when it is off. In this build no retrieval is implemented
    beyond the gate: reaching :meth:`_retrieve` raises, so there is no code path
    that performs a network call. A live subclass overrides :meth:`_retrieve`.
    """

    def __init__(self, *, allow_network: bool, versions: dict[str, SourceVersion] | None = None):
        self.allow_network = bool(allow_network)
        self._versions = dict(versions or {})
        self.call_log: list[dict] = []

    def versions(self) -> dict[str, SourceVersion]:
        return dict(self._versions)

    # -- gate ---------------------------------------------------------------
    def _gate(self, what: str, **detail) -> None:
        if not self.allow_network:
            raise BlockedError(_GATE_MESSAGE.format(what=what))
        self.call_log.append({"what": what, **detail})

    def _retrieve(self, what: str, **detail) -> Any:
        """Live retrieval. Not implemented in this build — never executed."""
        raise NotImplementedError(
            f"Live retrieval of {what!r} is not implemented in this build. "
            f"Stage E currently runs against a supplied AnnotationSource; "
            f"implementing this method is the only place a network call may appear."
        )

    # -- gated surface ------------------------------------------------------
    def secondary_structure(self, structure: Any) -> dict[int, str]:
        self._gate("secondary_structure")
        return self._retrieve("secondary_structure", structure=structure)

    def domains(self, acc: str) -> list[dict]:
        self._gate("domains", acc=acc)
        return self._retrieve("domains", acc=acc)

    def functional_sites(self, acc: str) -> list[dict]:
        self._gate("functional_sites", acc=acc)
        return self._retrieve("functional_sites", acc=acc)

    def conservation(self, acc: str, residues: Iterable[int]) -> dict[int, dict]:
        self._gate("conservation", acc=acc)
        return self._retrieve("conservation", acc=acc, residues=list(residues))

    def disease_associations(self, gene: str) -> list[dict]:
        self._gate("disease_associations", gene=gene)
        return self._retrieve("disease_associations", gene=gene)

    def literature_mechanisms(self, gene: str) -> list[MechanismRecord]:
        self._gate("literature_mechanisms", gene=gene)
        return self._retrieve("literature_mechanisms", gene=gene)
