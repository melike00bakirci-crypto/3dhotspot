"""Structural exports for 06_FINAL_HOTSPOTS/structures/ (methodology paragraph 57).

CA-only models, written deterministically as text so two runs produce byte-identical
files.

Residue names are written as ``UNK``. That is not a placeholder oversight: the
amino-acid identity column is **outside** ``stage_b_permitted_columns``, so Stage B
never loads it. The structural export therefore carries geometry and statistics, and
the residue identity is restored downstream by agents that are permitted to know it.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PDB_MAX_RESSEQ = 9999


@dataclass
class CaRecord:
    residue_index: int
    x: float
    y: float
    z: float
    bfactor: float = 0.0
    occupancy: float = 1.00
    chain: str = "A"
    resname: str = "UNK"


def _pdb_atom_line(serial: int, rec: CaRecord) -> str:
    return (
        f"ATOM  {serial:5d}  CA  {rec.resname:>3} {rec.chain:1}{rec.residue_index:4d}"
        f"    {rec.x:8.3f}{rec.y:8.3f}{rec.z:8.3f}"
        f"{rec.occupancy:6.2f}{rec.bfactor:6.2f}           C"   # element right-justified in 77-78
    )


def write_ca_pdb(path: str | Path, records: list[CaRecord], title: str,
                 remarks: list[str] | None = None) -> int:
    """Write a CA-only PDB. Returns the number of residues omitted by the format.

    The fixed-column PDB format cannot hold a residue number above 9999. Rather than
    silently wrapping such a number to a wrong residue, those residues are omitted
    here and listed in a REMARK; the mmCIF companion file carries the complete model.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    keep = [r for r in records if r.residue_index <= PDB_MAX_RESSEQ]
    omitted = len(records) - len(keep)

    lines = [f"TITLE     {title[:60]}"]
    for text in (remarks or []):
        lines.append(f"REMARK   1 {text[:69]}")
    if omitted:
        lines.append(
            f"REMARK   1 {omitted} residue(s) with index > {PDB_MAX_RESSEQ} omitted; "
            f"see the mmCIF file")
    for serial, rec in enumerate(keep, start=1):
        lines.append(_pdb_atom_line(serial, rec))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return omitted


_CIF_HEADER = """data_hotspot3d
#
_entry.id   hotspot3d
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num"""


def write_ca_cif(path: str | Path, records: list[CaRecord], comments: list[str] | None = None
                 ) -> Path:
    """Write a minimal, valid mmCIF ``atom_site`` loop (no residue-number ceiling)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# {c}" for c in (comments or [])]
    lines.append(_CIF_HEADER)
    for serial, rec in enumerate(records, start=1):
        lines.append(
            f"ATOM {serial} C CA {rec.resname} {rec.chain} {rec.residue_index} "
            f"{rec.x:.3f} {rec.y:.3f} {rec.z:.3f} {rec.occupancy:.2f} "
            f"{rec.bfactor:.2f} {rec.residue_index} {rec.chain} 1")
    lines.append("#")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


_PALETTE = ["red", "orange", "yellow", "green", "cyan", "blue", "magenta", "salmon"]


def write_chimerax_script(path: str | Path, hotspots: dict[str, list[int]], radius: float,
                          model_file: str = "final_hotspots.pdb") -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# ChimeraX visualization of the Stage B hotspot centers.",
        "# Run from this directory:  chimerax view_hotspots.cxc",
        f"# r_hot = {radius:g} A. Residue names are UNK by design: Stage B never loads",
        "# the amino-acid identity column (outside stage_b_permitted_columns).",
        f"open {model_file}",
        "style ball",
        "color #1 gray(70)",
    ]
    for i, (hid, residues) in enumerate(sorted(hotspots.items())):
        if not residues:
            continue
        spec = ",".join(str(r) for r in sorted(residues))
        colour = _PALETTE[i % len(_PALETTE)]
        lines.append(f"color /A:{spec} {colour}   # {hid}")
    lines.append("view")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def write_pymol_script(path: str | Path, hotspots: dict[str, list[int]], radius: float,
                       model_file: str = "final_hotspots.pdb") -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# PyMOL visualization of the Stage B hotspot centers.",
        "# Run from this directory:  pymol view_hotspots.pml",
        f"# r_hot = {radius:g} A. Residue names are UNK by design (see the .cxc script).",
        f"load {model_file}, hotspots",
        "hide everything",
        "show spheres, hotspots",
        "set sphere_scale, 0.4",
        "color grey70, hotspots",
    ]
    for i, (hid, residues) in enumerate(sorted(hotspots.items())):
        if not residues:
            continue
        spec = "+".join(str(r) for r in sorted(residues))
        colour = _PALETTE[i % len(_PALETTE)]
        lines.append(f"select {hid}, hotspots and resi {spec}")
        lines.append(f"color {colour}, {hid}")
    lines.append("deselect")
    lines.append("orient")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path
