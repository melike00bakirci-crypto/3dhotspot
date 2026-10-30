"""Protein HGVS parsing and amino-acid vocabulary — Stage A, METHOD_SPEC II.13.

Pure functions only: no I/O, no configuration, no decisions. A record either
yields a *missense* protein change (reference AA, 1-based position, alternate AA)
or it does not, and the reason it does not is returned rather than swallowed —
``unparseable_hgvs`` is a recorded exclusion, never a silent drop.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

# --- amino-acid vocabulary --------------------------------------------------

AA3_TO_1: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
# Recorded as non-standard: present in some models/records, never a missense
# reference or alternate under the frozen 20-letter canonical alphabet.
AA3_NONSTANDARD: dict[str, str] = {
    "SEC": "U", "PYL": "O", "MSE": "M", "ASX": "B", "GLX": "Z",
    "XAA": "X", "UNK": "X", "XLE": "J",
}
STANDARD_AA1: frozenset[str] = frozenset(AA3_TO_1.values())
TERMINATION_TOKENS: frozenset[str] = frozenset({"TER", "*", "TRM"})

# Change suffixes that make a record something other than a simple substitution.
_NON_SUBSTITUTION_TOKENS = (
    "fs", "del", "ins", "dup", "ext", "delins", "=/", "?", "[", "]",
)

# ``NM_004333.6(BRAF):c.1799T>A (p.Val600Glu)``
_NAME_RE = re.compile(
    r"^\s*(?P<transcript>[A-Za-z]{2}_[0-9]+(?:\.[0-9]+)?)"
    r"(?:\((?P<gene>[^)]*)\))?"
    r"(?::(?P<hgvs_c>[^\s]+))?"
    r"(?:\s*\((?P<hgvs_p>p\.[^)]+(?:\)[^)]*)?)\))?\s*$"
)
# fallback: pull any ``(p....)`` group out of a free-form Name
_P_IN_NAME_RE = re.compile(r"\((?P<hgvs_p>p\.[^()]+)\)")
_TRANSCRIPT_RE = re.compile(r"^\s*(?P<transcript>[A-Za-z]{2}_[0-9]+(?:\.[0-9]+)?)")

_HGVS_P_RE = re.compile(
    r"^p\.\(?"
    r"(?P<ref>[A-Za-z]{3}|[A-Z])"
    r"(?P<pos>[0-9]+)"
    r"(?P<alt>[A-Za-z]{3}|[A-Z]|\*|=)"
    r"\)?$"
)


def strip_version(accession: str | None) -> str:
    """``NM_004333.6`` -> ``NM_004333``. Versions are recorded, never compared."""
    text = str(accession or "").strip()
    return text.split(".", 1)[0]


def aa_one_letter(token: str | None) -> str | None:
    """Normalize a 1- or 3-letter amino-acid token to its 1-letter code.

    Returns ``None`` for anything that is not a standard proteinogenic residue,
    including termination and non-standard residues, which the caller records.
    """
    text = str(token or "").strip()
    upper = text.upper()
    one = {1: upper, 3: AA3_TO_1.get(upper)}.get(len(upper))
    return one if (one in STANDARD_AA1) else None


def is_termination(token: str | None) -> bool:
    return str(token or "").strip().upper() in TERMINATION_TOKENS


# --- parsed protein change --------------------------------------------------

@dataclass(frozen=True)
class ProteinChange:
    """One parsed ``p.`` change. ``is_missense`` is the only inclusion signal."""

    raw: str
    ref_aa: str | None = None
    position: int | None = None
    alt_aa: str | None = None
    is_missense: bool = False
    reason: str | None = None          # why it is not a usable missense change

    @property
    def short(self) -> str:
        parts = (self.ref_aa, self.position, self.alt_aa)
        return "".join(str(p) for p in parts) if all(p is not None for p in parts) else "NA"


def parse_protein_hgvs(text: str | None) -> ProteinChange:
    """Parse a protein-level HGVS string into a :class:`ProteinChange`.

    Accepts ``p.Val600Glu``, ``p.V600E`` and the parenthesised ``p.(Val600Glu)``
    predicted form. Everything else — synonymous, nonsense, stop-loss, frameshift,
    indel, extension, unparseable — comes back with ``is_missense=False`` and a
    machine-readable ``reason``.
    """
    raw = str(text or "").strip()
    if not raw:
        return ProteinChange(raw="NA", reason="empty_protein_hgvs")
    if not raw.startswith("p."):
        return ProteinChange(raw=raw, reason="not_protein_level_hgvs")

    body = raw[2:]
    lowered = body.lower()
    for token in _NON_SUBSTITUTION_TOKENS:
        if token in lowered:
            return ProteinChange(raw=raw, reason=f"not_a_substitution:{token}")

    match = _HGVS_P_RE.match(raw)
    if match is None:
        return ProteinChange(raw=raw, reason="hgvs_p_pattern_unmatched")

    ref_token, alt_token = match.group("ref"), match.group("alt")
    position = int(match.group("pos"))

    if is_termination(ref_token):
        return ProteinChange(raw=raw, position=position, reason="reference_is_termination")
    ref_aa = aa_one_letter(ref_token)
    if ref_aa is None:
        return ProteinChange(raw=raw, position=position,
                             reason=f"non_standard_reference_aa:{ref_token}")

    if alt_token == "=":
        return ProteinChange(raw=raw, ref_aa=ref_aa, position=position, alt_aa=ref_aa,
                             reason="synonymous")
    if is_termination(alt_token):
        return ProteinChange(raw=raw, ref_aa=ref_aa, position=position,
                             reason="nonsense")
    alt_aa = aa_one_letter(alt_token)
    if alt_aa is None:
        return ProteinChange(raw=raw, ref_aa=ref_aa, position=position,
                             reason=f"non_standard_alternate_aa:{alt_token}")
    if alt_aa == ref_aa:
        return ProteinChange(raw=raw, ref_aa=ref_aa, position=position, alt_aa=alt_aa,
                             reason="synonymous")
    if position < 1:
        return ProteinChange(raw=raw, ref_aa=ref_aa, position=position, alt_aa=alt_aa,
                             reason="non_positive_position")

    return ProteinChange(raw=raw, ref_aa=ref_aa, position=position, alt_aa=alt_aa,
                         is_missense=True)


@dataclass(frozen=True)
class ParsedName:
    """Decomposition of the ClinVar ``Name`` field."""

    transcript: str | None = None
    gene: str | None = None
    hgvs_c: str | None = None
    hgvs_p: str | None = None


def parse_clinvar_name(name: str | None) -> ParsedName:
    """Split ``NM_004333.6(BRAF):c.1799T>A (p.Val600Glu)`` into its parts."""
    text = str(name or "").strip()
    if not text:
        return ParsedName()

    match = _NAME_RE.match(text)
    if match is not None:
        return ParsedName(transcript=match.group("transcript"), gene=match.group("gene"),
                          hgvs_c=match.group("hgvs_c"), hgvs_p=match.group("hgvs_p"))

    transcript_match = _TRANSCRIPT_RE.match(text)
    protein_match = _P_IN_NAME_RE.search(text)
    return ParsedName(
        transcript=transcript_match.group("transcript") if transcript_match else None,
        gene=None,
        hgvs_c=None,
        hgvs_p=protein_match.group("hgvs_p") if protein_match else None,
    )
