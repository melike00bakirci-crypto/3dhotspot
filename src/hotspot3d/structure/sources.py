"""Structure acquisition interface and the AlphaFold DB client — METHOD_SPEC II.13.

FROZEN acquisition policy:
  * the model comes from the AlphaFold Database under the pinned
    ``structure.model_template`` in ``config/pipeline.yaml`` — that config value
    is the single authority for which version is analysed, and the version
    actually retrieved is recorded alongside AFDB's current latest. The template
    is deliberately NOT resolved from the API response: letting AFDB decide would
    change the structural input between runs with nobody deciding. Moving the pin
    is a pre-registered revision under a new RUN_ID (v4 -> v6 on 2026-08-17, after
    AFDB withdrew v4 and v5);
  * a **multi-fragment** AFDB entry is BLOCKED and escalated — fragments are
    never stitched;
  * a **missing** AFDB entry is escalated;
  * homology models, ortholog structures and experimental substitutes are
    **never** used.

The model itself is held verbatim. Stage A extracts coordinates and pLDDT from
it; it never edits it.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Protocol, Sequence, runtime_checkable

from ..data.hgvs import AA3_NONSTANDARD, AA3_TO_1
from ..data.sources import OLIGOMERIC_STATE_OLIGOMERIC
from ..utils.errors import BlockedError, EscalationRequired

AFDB_PREDICTION_URL = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
AFDB_CIF_URL = "https://alphafold.ebi.ac.uk/files/{model_id}.cif"
#: Fallback only — every live client is constructed with the template from
#: ``structure.model_template``, which is the authority. Kept in step with the
#: config so the default cannot silently request a version AFDB has withdrawn
#: (it read model_v4 after the 2026-08-17 revision to v6, and v4 now 404s).
AFDB_MODEL_TEMPLATE = "AF-{acc}-F1-model_v6"


def _version_of(model_id: str, default: str = "NA") -> str:
    """The version encoded in an AFDB model id, e.g. 'AF-P16389-F1-model_v6' -> '6'."""
    import re

    match = re.search(r"-model_v(\d+)$", str(model_id or ""))
    return match.group(1) if match else default


def _fragment_number(entry_id: str) -> int:
    """The fragment index in an AFDB entry id, e.g. 'AF-Q8WZ42-F7' -> 7.

    Fragments of one accession are distinguished by this suffix; an unparseable
    id sorts last rather than silently claiming to be fragment 1.
    """
    import re

    match = re.search(r"-F(\d+)$", str(entry_id or ""))
    return int(match.group(1)) if match else 1 << 30


def _accession_of(entry_id: str) -> str | None:
    """The accession in an AFDB entry id: ``AF-<accession>-F<n>``.

    The accession part is matched greedily so an ISOFORM keeps its suffix —
    ``AF-P16389-2-F1`` -> ``P16389-2``, a different molecule from
    ``AF-P16389-F1`` -> ``P16389``. This is what lets an isoform be recognised
    even in a payload that does not spell out ``uniprotAccession``.
    """
    import re

    match = re.match(r"^AF-(.+)-F\d+$", str(entry_id or "").strip())
    return match.group(1) if match else None

BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})


@dataclass(frozen=True)
class Atom:
    """One atom record. ``b_factor`` carries pLDDT in an AlphaFold model."""

    name: str
    x: float
    y: float
    z: float
    occupancy: float = 1.0
    b_factor: float = 0.0
    element: str = "C"
    alt_loc: str = "."

    @property
    def xyz(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def is_finite(self) -> bool:
        return all(math.isfinite(v) for v in (self.x, self.y, self.z))


@dataclass
class StructureResidue:
    """One modelled residue, with whatever atoms the model actually provides."""

    residue_index: int
    aa3: str
    chain_id: str
    atoms: dict[str, Atom] = field(default_factory=dict)

    # -- identity ------------------------------------------------------------
    @property
    def aa(self) -> str:
        """1-letter code; ``X`` for anything outside the standard 20."""
        upper = self.aa3.strip().upper()
        return AA3_TO_1.get(upper, AA3_NONSTANDARD.get(upper, "X"))

    @property
    def is_standard(self) -> bool:
        return self.aa3.strip().upper() in AA3_TO_1

    # -- representation (F1: CA only; CB / SC emitted but unused) ------------
    @property
    def ca(self) -> Atom | None:
        return self.atoms.get("CA")

    @property
    def cb(self) -> Atom | None:
        return self.atoms.get("CB")

    @property
    def ca_usable(self) -> bool:
        """Present, non-zero occupancy, finite coordinates (U_struct membership)."""
        atom = self.ca
        return bool(atom) and atom.occupancy > 0.0 and atom.is_finite

    @property
    def plddt(self) -> float | None:
        """Per-residue pLDDT = the CA B-factor of an AlphaFold model."""
        atom = self.ca
        return atom.b_factor if atom is not None else None

    def sidechain_centroid(self) -> tuple[float, float, float] | None:
        """Mean of side-chain heavy atoms. Emitted for future flexibility only.

        Marked ``unused_under_F1`` in every output: the frozen representation is
        CA and the frozen distance is CA-CA Euclidean.
        """
        heavy = [a for name, a in sorted(self.atoms.items())
                 if name not in BACKBONE_ATOMS and a.element.upper() != "H" and a.is_finite]
        if not heavy:
            return None
        n = float(len(heavy))
        return (sum(a.x for a in heavy) / n,
                sum(a.y for a in heavy) / n,
                sum(a.z for a in heavy) / n)


@dataclass
class StructureModel:
    """A retrieved structural model, held as delivered.

    ``n_fragments`` must be declared by the source. ``None`` means "the source
    could not tell us", which is escalated rather than assumed to be 1 — an
    undeclared fragment count is exactly the situation METHOD_SPEC II.13 blocks.
    """

    accession: str
    entry_id: str
    model_version: str
    db_version: str
    residues: list[StructureResidue] = field(default_factory=list)
    n_fragments: int | None = None
    raw_cif_text: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    # -- derived, descriptive ------------------------------------------------
    @property
    def chain_ids(self) -> list[str]:
        return sorted({r.chain_id for r in self.residues})

    @property
    def residue_indices(self) -> list[int]:
        return [r.residue_index for r in self.residues]

    @property
    def sequence(self) -> str:
        return "".join(r.aa for r in sorted(self.residues, key=lambda r: r.residue_index))

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    @property
    def n_atoms(self) -> int:
        return sum(len(r.atoms) for r in self.residues)

    def by_index(self) -> dict[int, StructureResidue]:
        return {r.residue_index: r for r in self.residues}

    def model_id(self) -> str:
        return self.entry_id or AFDB_MODEL_TEMPLATE.format(acc=self.accession)


# --- provider interfaces ----------------------------------------------------

@runtime_checkable
class StructureSource(Protocol):
    """Anything that can deliver a structural model for a UniProt accession."""

    def fetch_model(self, uniprot_acc: str) -> StructureModel: ...

    def model_metadata(self) -> dict: ...


def assert_model_admissible(model: StructureModel) -> None:
    """Apply the FROZEN II.13 fallback/refusal rules to a retrieved model.

    Raises :class:`BlockedError` for a multi-fragment entry (the caller records
    ``MULTI_FRAGMENT_AFDB_ENTRY``) and :class:`EscalationRequired` when the
    source did not declare enough to decide.
    """
    if model is None:
        raise EscalationRequired(
            ambiguity="No AlphaFold DB entry was returned for the resolved accession.",
            options=["Escalate and stop", "Substitute an experimental structure",
                     "Use a homology model or an ortholog"],
            consequences=[
                "Stage A stops; no cohort is produced and nothing is claimed.",
                "Requires an explicit Lead authorization plus a numbering "
                "reconciliation; config.structure.allow_experimental_substitute is FALSE.",
                "Prohibited by METHOD_SPEC II.13 — coordinates would not correspond "
                "to the analysed sequence.",
            ],
            recommendation="Escalate and stop. Homology and ortholog models are never used.",
        )
    if model.n_fragments is None:
        raise EscalationRequired(
            ambiguity=(
                f"The structure source did not declare a fragment count for "
                f"{model.accession}. Multi-fragment entries must be BLOCKED, and an "
                f"undeclared count cannot be assumed to be 1."
            ),
            options=["Escalate and stop",
                     "Assume a single fragment and continue"],
            consequences=[
                "Stage A stops until the source declares the fragment count.",
                "Risks silently analysing one fragment of a fragmented entry as if "
                "it were the whole protein — coordinates and numbering would be wrong.",
            ],
            recommendation="Escalate and stop; require the source to declare n_fragments.",
        )
    if model.n_fragments > 1:
        raise BlockedError(
            f"MULTI_FRAGMENT_AFDB_ENTRY — AFDB serves {model.n_fragments} fragments for "
            f"{model.accession}. Fragment stitching is prohibited (METHOD_SPEC II.13); "
            f"this is BLOCKED and escalated to the Lead."
        )
    if not model.residues:
        raise BlockedError(
            f"STRUCTURE_MAPPING_FAILURE — model {model.model_id()} contains no residues."
        )


def assert_oligomeric_assembly_permitted(oligomeric_state: str, *,
                                         allow_monomer: bool) -> None:
    """Workflow v2 §2 — an obligate oligomer analysed on a monomer model is
    missing every inter-subunit neighbourhood a real variant sits in.

    A no-op unless UniProt's own subunit annotation gave a *confident*
    ``oligomeric`` reading (:func:`hotspot3d.data.sources.classify_oligomeric_state`).
    ``monomer`` and ``unknown`` states both proceed silently here — this pipeline
    has no multimer-fetching capability at all, so there is nothing this check
    could do differently for either of them; the caller records the state either
    way (mirroring how :func:`assert_model_admissible` never raises on n_fragments
    <= 1).

    Raises :class:`BlockedError` for ``OLIGOMERIC_ASSEMBLY_UNAVAILABLE`` when the
    monomer exception is not authorized. When it IS authorized (the current
    config default), this function does not raise; the caller is responsible for
    recording ``OLIGOMERIC_ASSEMBLY_MONOMER_DECLARED`` and propagating
    ``intra_subunit_only`` downstream.
    """
    if oligomeric_state != OLIGOMERIC_STATE_OLIGOMERIC:
        return
    if not allow_monomer:
        raise BlockedError(
            "OLIGOMERIC_ASSEMBLY_UNAVAILABLE — UniProt's subunit annotation indicates "
            "the protein assembles as a homo- or hetero-oligomer. This pipeline has no "
            "multimer-fetching capability, and "
            "structure.allow_monomer_for_obligate_oligomer is FALSE, so proceeding on "
            "the monomeric AlphaFold model is not authorized under this configuration. "
            "The run is BLOCKED rather than silently missing every inter-subunit "
            "neighbourhood a real variant may occupy (METHOD_SPEC/Workflow v2 §2)."
        )


# --- live client (request shaping + parsing; network gated) ------------------

@dataclass
class AlphaFoldClient:
    """AlphaFold DB client. Every network call is gated on ``allow_network``.

    Request shaping and response parsing are pure and unit-tested; only
    :meth:`_get` touches the network, and it refuses to unless explicitly enabled.
    """

    allow_network: bool = False
    timeout: int = 120
    model_template: str = AFDB_MODEL_TEMPLATE
    _last_metadata: dict = field(default_factory=dict, init=False, repr=False)

    # -- request shaping -----------------------------------------------------
    def build_prediction_request(self, uniprot_acc: str) -> dict:
        return {
            "url": AFDB_PREDICTION_URL.format(acc=uniprot_acc),
            "method": "GET",
            "headers": {"Accept": "application/json"},
            "purpose": "enumerate AFDB predictions to determine the fragment count",
        }

    def build_cif_request(self, uniprot_acc: str, model_id: str | None = None) -> dict:
        mid = model_id or self.model_template.format(acc=uniprot_acc)
        return {
            "url": AFDB_CIF_URL.format(model_id=mid),
            "method": "GET",
            "headers": {"Accept": "text/plain"},
            "model_id": mid,
            "purpose": "retrieve the mmCIF model verbatim",
        }

    # -- response parsing (pure) --------------------------------------------
    @staticmethod
    def parse_prediction_payload(payload: Sequence[dict], uniprot_acc: str) -> dict:
        """Read the AFDB prediction list into the metadata Stage A must record."""
        if not payload:
            raise EscalationRequired(
                ambiguity=f"AlphaFold DB has no prediction for {uniprot_acc}.",
                options=["Escalate and stop", "Use a non-AFDB structure"],
                consequences=["Stage A stops and nothing is claimed.",
                              "Prohibited by METHOD_SPEC II.13."],
                recommendation="Escalate and stop.",
            )
        # AFDB returns ONE ENTRY PER ISOFORM, not per fragment. For P16389 the
        # payload carries AF-P16389-F1 (acc P16389, residues 1-499) *and*
        # AF-P16389-2-F1 (acc P16389-2, residues 1-356) — two isoforms, each a
        # single-fragment model. Counting payload entries therefore reported two
        # fragments for a 499-residue protein that AFDB cannot fragment at all
        # (it fragments only above ~2700 residues), and the II.13 multi-fragment
        # rule fired on a condition that did not exist. That false positive would
        # have blocked any protein with more than one AFDB isoform.
        #
        # Fragmentation is expressed by the -F<n> suffix on entries sharing ONE
        # accession. An isoform accession ("P16389-2") is a different molecule,
        # never a fragment of the canonical one, so it is excluded here. Selecting
        # the canonical isoform introduces no new decision: resolve_gene already
        # returns the canonical accession under the frozen MANE_SELECT /
        # UniProt-canonical transcript policy.
        # An entry is set aside only when it is demonstrably a DIFFERENT molecule.
        # uniprotAccession is the primary evidence; when it is absent the entryId
        # carries the same fact (AF-P16389-2-F1 -> P16389-2), and reading it there
        # matters: treating an unlabelled entry as canonical would reinstate the
        # isoform-as-fragment miscount on any payload that omitted the field.
        # Only when NEITHER names an accession is the entry taken to be the
        # requested one, because inferring exclusion from nothing at all would
        # silently drop the canonical entry.
        def _accession(entry: dict) -> str:
            declared = str(entry.get("uniprotAccession") or "").strip()
            return declared or (_accession_of(entry.get("entryId")) or "")

        canonical = [e for e in payload
                     if not _accession(e) or _accession(e) == str(uniprot_acc)]
        returned = [_accession(e) for e in payload]
        if not canonical:
            raise EscalationRequired(
                ambiguity=(f"AlphaFold DB returned {len(payload)} prediction(s) for "
                           f"{uniprot_acc} but none for that accession itself; the "
                           f"payload carries only {sorted(set(returned))}."),
                options=["Escalate and stop", "Use a non-AFDB structure"],
                consequences=["Stage A stops and nothing is claimed.",
                              "Prohibited by METHOD_SPEC II.13."],
                recommendation="Escalate and stop.",
            )
        # Deterministic fragment order, so the entry chosen never depends on the
        # order AFDB happened to serve.
        canonical.sort(key=lambda e: _fragment_number(e.get("entryId")))
        first = canonical[0]
        excluded = sorted({a for a in returned if a and a != str(uniprot_acc)})
        return {
            "uniprot_acc": uniprot_acc,
            # Fragments OF THE CANONICAL ACCESSION. The II.13 rule below is
            # unchanged and still blocks when this is genuinely > 1.
            "n_fragments": len(canonical),
            "entry_id": first.get("entryId"),
            "model_version": str(first.get("latestVersion", "NA")),
            "db_version": str(first.get("modelCreatedDate", "NA")),
            "cif_url": first.get("cifUrl"),
            "uniprot_sequence_length": first.get("uniprotEnd"),
            "model_entity_id": first.get("modelEntityId"),
            # Provenance: what AFDB offered and what was set aside, so the choice
            # of isoform is auditable rather than implicit.
            "afdb_accessions_returned": returned,
            "afdb_isoform_accessions_excluded": excluded,
            "afdb_fragment_entry_ids": [e.get("entryId") for e in canonical],
        }

    # -- network (gated) -----------------------------------------------------
    def _get(self, url: str, accept: str) -> bytes:
        if not self.allow_network:
            raise BlockedError(
                f"BLOCKED — network retrieval is disabled for this run "
                f"(allow_network=FALSE); refusing to fetch {url}. Stage A never "
                f"substitutes a cached, hand-made or previous-run file for a "
                f"retrieval it did not perform."
            )
        import requests                                       # pragma: no cover

        response = requests.get(url, headers={"Accept": accept}, timeout=self.timeout)
        response.raise_for_status()
        return response.content

    def fetch_model(self, uniprot_acc: str) -> StructureModel:  # pragma: no cover
        from .cif import parse_mmcif

        predictions = self._get(self.build_prediction_request(uniprot_acc)["url"],
                                "application/json")
        import json

        meta = self.parse_prediction_payload(json.loads(predictions), uniprot_acc)
        # The CONFIG TEMPLATE builds the URL, not the entryId from the payload.
        # Passing meta["entry_id"] here made build_cif_request ignore
        # model_template entirely (it consults the template only when model_id is
        # None), so the frozen pin never reached the wire: the request became
        # .../AF-<acc>-F1.cif, which AFDB does not serve (404). The pin is the
        # decision about WHICH model version is analysed; it has to govern the
        # retrieval or it is decorative.
        request = self.build_cif_request(uniprot_acc)
        cif_text = self._get(request["url"], "text/plain").decode("utf-8")

        model = parse_mmcif(cif_text, accession=uniprot_acc)
        model.entry_id = meta.get("entry_id") or model.entry_id
        # Record the version RETRIEVED, which is the pinned one, and keep AFDB's
        # current latest beside it. These agree only while the pin is current; a
        # pin deliberately lags once AFDB publishes a newer version, and at that
        # point recording latestVersion would misreport the analysed structure.
        meta["afdb_latest_version"] = meta["model_version"]
        meta["requested_model_id"] = request["model_id"]
        meta["model_version"] = _version_of(request["model_id"],
                                            default=meta["model_version"])
        model.model_version = meta["model_version"]
        model.db_version = meta["db_version"]
        model.n_fragments = meta["n_fragments"]
        model.metadata.update(meta)
        self._last_metadata = {**meta, "cif_request": request,
                               "retrieved_from": request["url"]}
        return model

    def model_metadata(self) -> dict:
        return {"source": "alphafold_db", "allow_network": self.allow_network,
                "model_template": self.model_template, **self._last_metadata}
