"""Minimal, deterministic mmCIF reader/writer for AlphaFold models.

Only the ``_atom_site`` loop (coordinates, occupancy, B-factor) and the
``_ma_qa_metric_local`` loop (per-residue pLDDT) are interpreted; everything else
in a retrieved file is carried through untouched as ``raw_cif_text``. Stage A
reads structures, it does not edit them.

In an AlphaFold model the atom B-factor *is* pLDDT, so
``structures/structure_with_plddt.cif`` is the model verbatim whenever the source
supplied its own text; the equality is verified and recorded rather than assumed.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from ..utils.errors import BlockedError
from .sources import Atom, StructureModel, StructureResidue

ATOM_SITE_PREFIX = "_atom_site."
QA_METRIC_PREFIX = "_ma_qa_metric_local."
NULL_TOKENS = frozenset({".", "?"})


# --- tokenizing -------------------------------------------------------------

def split_cif_row(line: str) -> list[str]:
    """Split one mmCIF loop row, honouring single and double quoting."""
    tokens: list[str] = []
    buf: list[str] = []
    quote: str | None = None
    for char in line.strip():
        if quote is not None:
            if char == quote:
                quote = None
                tokens.append("".join(buf))
                buf = []
            else:
                buf.append(char)
        elif char in ("'", '"'):
            quote = char
        elif char.isspace():
            if buf:
                tokens.append("".join(buf))
                buf = []
        else:
            buf.append(char)
    if buf:
        tokens.append("".join(buf))
    return tokens


def _as_float(token: str, default: float = 0.0) -> float:
    if token in NULL_TOKENS:
        return default
    try:
        return float(token)
    except ValueError:
        return float("nan")


def _as_int(token: str, default: int | None = None) -> int | None:
    if token in NULL_TOKENS:
        return default
    try:
        return int(float(token))
    except ValueError:
        return default


@dataclass
class _Loop:
    columns: list[str]
    rows: list[list[str]]

    def dicts(self):
        for row in self.rows:
            yield dict(zip(self.columns, row))


def _read_loop(lines: list[str], start: int, prefix: str) -> tuple[_Loop | None, int]:
    """Read the ``loop_`` block beginning at ``start`` if it carries ``prefix``."""
    i = start + 1
    columns: list[str] = []
    while i < len(lines) and lines[i].strip().startswith("_"):
        columns.append(lines[i].strip())
        i += 1
    if not columns or not columns[0].startswith(prefix):
        return None, i

    rows: list[list[str]] = []
    while i < len(lines):
        text = lines[i].strip()
        if not text or text.startswith("#") or text.startswith("loop_") or text.startswith("_"):
            break
        if text.startswith("data_"):
            break
        row = split_cif_row(lines[i])
        if len(row) == len(columns):
            rows.append(row)
        i += 1
    return _Loop(columns=columns, rows=rows), i


def parse_mmcif(text: str, accession: str = "NA") -> StructureModel:
    """Parse an mmCIF into a :class:`StructureModel`, preserving the raw text."""
    lines = text.splitlines()
    entry_id = "NA"
    atom_loop: _Loop | None = None
    qa_loop: _Loop | None = None

    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        if stripped.startswith("data_"):
            entry_id = stripped[len("data_"):].strip() or entry_id
            i += 1
            continue
        if stripped == "loop_":
            for prefix, slot in ((ATOM_SITE_PREFIX, "atom"), (QA_METRIC_PREFIX, "qa")):
                loop, nxt = _read_loop(lines, i, prefix)
                if loop is not None:
                    if slot == "atom":
                        atom_loop = loop
                    else:
                        qa_loop = loop
                    i = nxt
                    break
            else:
                _, i = _read_loop(lines, i, "_")
            continue
        i += 1

    if atom_loop is None:
        raise BlockedError(
            "STRUCTURE_MAPPING_FAILURE — mmCIF payload has no _atom_site loop; "
            "the retrieved file is not a usable coordinate model."
        )

    residues: dict[tuple[str, int], StructureResidue] = {}
    for row in atom_loop.dicts():
        if row.get("_atom_site.group_PDB", "ATOM") not in ("ATOM", "HETATM"):
            continue
        index = _as_int(row.get("_atom_site.label_seq_id", "."))
        if index is None:
            index = _as_int(row.get("_atom_site.auth_seq_id", "."))
        if index is None:
            continue
        chain = (row.get("_atom_site.auth_asym_id")
                 or row.get("_atom_site.label_asym_id") or "A")
        aa3 = (row.get("_atom_site.label_comp_id")
               or row.get("_atom_site.auth_comp_id") or "UNK")
        name = (row.get("_atom_site.label_atom_id")
                or row.get("_atom_site.auth_atom_id") or "").strip('"')

        residue = residues.setdefault((chain, index),
                                      StructureResidue(residue_index=index, aa3=aa3,
                                                       chain_id=chain))
        residue.atoms[name] = Atom(
            name=name,
            x=_as_float(row.get("_atom_site.Cartn_x", ".")),
            y=_as_float(row.get("_atom_site.Cartn_y", ".")),
            z=_as_float(row.get("_atom_site.Cartn_z", ".")),
            occupancy=_as_float(row.get("_atom_site.occupancy", "."), 1.0),
            b_factor=_as_float(row.get("_atom_site.B_iso_or_equiv", "."), 0.0),
            element=(row.get("_atom_site.type_symbol") or "C").strip(),
            alt_loc=(row.get("_atom_site.label_alt_id") or ".").strip(),
        )

    qa_plddt: dict[int, float] = {}
    if qa_loop is not None:
        for row in qa_loop.dicts():
            index = _as_int(row.get("_ma_qa_metric_local.label_seq_id", "."))
            value = _as_float(row.get("_ma_qa_metric_local.metric_value", "."), float("nan"))
            if index is not None and math.isfinite(value):
                qa_plddt[index] = value

    ordered = [residues[key] for key in sorted(residues, key=lambda k: (k[0], k[1]))]
    return StructureModel(
        accession=accession, entry_id=entry_id, model_version="NA", db_version="NA",
        residues=ordered, raw_cif_text=text,
        metadata={"n_atom_site_rows": len(atom_loop.rows),
                  "ma_qa_metric_local_present": qa_loop is not None,
                  "ma_qa_metric_local_plddt": qa_plddt},
    )


# --- writing ----------------------------------------------------------------

_ATOM_SITE_COLUMNS = (
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
    "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
    "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
    "B_iso_or_equiv", "auth_seq_id", "auth_comp_id", "auth_asym_id",
    "auth_atom_id", "pdbx_PDB_model_num",
)


def render_mmcif(model: StructureModel) -> str:
    """Render a model as mmCIF with ``B_iso_or_equiv`` = pLDDT.

    Used only when the source provided no original text (synthetic runs). When
    original text exists it is written verbatim — the model is never rewritten.
    """
    out = [f"data_{model.model_id()}", "#",
           "_entry.id   " + model.model_id(), "#",
           "loop_"]
    out += [f"_atom_site.{col}" for col in _ATOM_SITE_COLUMNS]

    serial = 0
    for residue in sorted(model.residues, key=lambda r: (r.chain_id, r.residue_index)):
        plddt = residue.plddt
        for name, atom in sorted(residue.atoms.items()):
            serial += 1
            b_iso = plddt if plddt is not None else atom.b_factor
            out.append(
                f"ATOM {serial} {atom.element} {name} . {residue.aa3} "
                f"{residue.chain_id} 1 {residue.residue_index} . "
                f"{atom.x:.3f} {atom.y:.3f} {atom.z:.3f} {atom.occupancy:.2f} "
                f"{b_iso:.2f} {residue.residue_index} {residue.aa3} "
                f"{residue.chain_id} {name} 1"
            )
    out.append("#")
    return "\n".join(out) + "\n"


def write_structure_with_plddt(model: StructureModel, path: str | Path) -> dict:
    """Write ``structures/structure_with_plddt.cif``.

    Returns the record of *how* it was written, so a reader can tell a verbatim
    copy from a rendered one without re-deriving it.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    verbatim = model.raw_cif_text is not None
    text = model.raw_cif_text if verbatim else render_mmcif(model)
    path.write_text(text, encoding="utf-8", newline="\n")
    return {
        "path": str(path),
        "mode": "verbatim_source_model" if verbatim else "rendered_from_parsed_model",
        "b_factor_column": "pLDDT",
        "model_edited": False,
        "n_residues": model.n_residues,
        "n_atoms": model.n_atoms,
    }
