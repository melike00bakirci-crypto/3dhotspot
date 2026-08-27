"""ClinVar record normalization and the ``variant_summary`` parser.

Stage A ingests **every** missense record for the gene, at every review star
level (F12 / Decision 1). This module turns a raw ``variant_summary`` row into a
normalized record and records, for every row that does not survive, the reason
code it failed on. Nothing is dropped silently and nothing is dropped by star.
"""
from __future__ import annotations

import gzip
import io
from dataclasses import dataclass, field
from typing import Iterable, Iterator, Sequence

from .hgvs import ProteinChange, parse_clinvar_name, parse_protein_hgvs, strip_version
from .review_status import review_status_label, star_level

# Columns of the weekly ``variant_summary.txt.gz`` release that Stage A reads.
# The full file carries ~40 columns; the rest are carried through untouched in
# the immutable 01_INPUT_RAW payload and are never interpreted here.
VARIANT_SUMMARY_FIELDS: dict[str, tuple[str, ...]] = {
    "variation_id": ("VariationID",),
    "allele_id": ("AlleleID", "#AlleleID"),
    "gene_symbol": ("GeneSymbol",),
    "name": ("Name",),
    "variant_type": ("Type",),
    "clinical_significance": ("ClinicalSignificance",),
    "clinsig_simple": ("ClinSigSimple",),
    "review_status": ("ReviewStatus",),
    "condition_text": ("PhenotypeList",),
    "last_evaluated": ("LastEvaluated",),
    "n_submitters": ("NumberSubmitters",),
    "assembly": ("Assembly",),
    "chromosome": ("Chromosome",),
    "position_vcf": ("PositionVCF", "Start"),
    "ref_allele": ("ReferenceAlleleVCF", "ReferenceAllele"),
    "alt_allele": ("AlternateAlleleVCF", "AlternateAllele"),
    "rcv": ("RCVaccession",),
    "origin": ("Origin",),
    "submitter_categories": ("SubmitterCategories",),
}

#: Free-text fields held by Stage A and forbidden to Stage B's primary path.
PASSTHROUGH_TEXT_FIELDS = ("condition_text", "name", "rcv", "origin")

#: ClinVar ``Type`` values that can carry a small, point-level protein-coding
#: change. [REVISED — Workflow v2 §1, Lead ruling 2026-08-18, Finding B] added
#: "deletion", "duplication", "insertion" — real frameshift/indel variants are
#: routinely typed by ClinVar as one of these rather than as "Indel", so
#: without them such records never reached stratum 1 at all: not "excluded with
#: a reason" but silently absent from n_retrieved before any accounting layer
#: saw them. hgvs.parse_protein_hgvs() already resolves the records this
#: promotes into their correct true-cause reason (frameshift / in-frame indel);
#: no further change was needed there. CNV, fusion, translocation, inversion
#: and microsatellite stay excluded — none is a point-level protein-coding
#: change a residue-level HGVS parse can represent.
SUBSTITUTION_TYPES = frozenset({
    "single nucleotide variant", "snv", "variation", "indel", "missense variant",
    "deletion", "duplication", "insertion",
})


@dataclass
class ClinVarRecord:
    """One normalized ClinVar record. Immutable evidence; never edited in place."""

    record_id: str
    variation_id: str = "NA"
    allele_id: str = "NA"
    gene_symbol: str = "NA"
    name: str = "NA"
    variant_type: str = "NA"
    clinical_significance: str = "NA"
    review_status: str = "NA"
    star_level: int | None = None
    transcript: str | None = None
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    condition_text: str = "NA"
    last_evaluated: str = "NA"
    n_submitters: int | None = None
    assembly: str = "NA"
    chromosome: str = "NA"
    position_vcf: str = "NA"
    ref_allele: str = "NA"
    alt_allele: str = "NA"
    rcv: str = "NA"
    origin: str = "NA"
    change: ProteinChange = field(default_factory=lambda: parse_protein_hgvs(None))
    raw: dict = field(default_factory=dict)

    # -- derived, purely descriptive ----------------------------------------
    @property
    def residue_index(self) -> int | None:
        return self.change.position

    @property
    def aa_ref_hgvs(self) -> str | None:
        return self.change.ref_aa

    @property
    def aa_alt_hgvs(self) -> str | None:
        return self.change.alt_aa

    @property
    def is_missense(self) -> bool:
        return self.change.is_missense


def _first_present(row: dict, keys: Sequence[str]) -> str | None:
    for key in keys:
        value = row.get(key)
        text = str(value).strip() if value is not None else ""
        if text and text.lower() not in ("na", "-", "nan", "none"):
            return text
    return None


def _as_int(value) -> int | None:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def normalize_row(row: dict, index: int, gene: str | None = None) -> ClinVarRecord:
    """Normalize one raw row (``variant_summary`` names or already-normalized).

    Accepting both shapes is what lets the synthetic fixtures and the live
    ``ClinVarClient`` feed byte-identical code paths.
    """
    def pick(field_name: str) -> str | None:
        direct = row.get(field_name)
        if direct is not None and str(direct).strip():
            return str(direct).strip()
        return _first_present(row, VARIANT_SUMMARY_FIELDS.get(field_name, ()))

    variation_id = pick("variation_id") or "NA"
    allele_id = pick("allele_id") or "NA"
    name = pick("name") or "NA"

    parsed_name = parse_clinvar_name(name)
    transcript = pick("transcript") or parsed_name.transcript
    hgvs_c = pick("hgvs_c") or parsed_name.hgvs_c
    hgvs_p = pick("hgvs_p") or parsed_name.hgvs_p

    review = pick("review_status") or "NA"
    record_id = f"CV{variation_id}" if variation_id != "NA" else f"REC{index:06d}"

    return ClinVarRecord(
        record_id=record_id,
        variation_id=variation_id,
        allele_id=allele_id,
        gene_symbol=pick("gene_symbol") or (gene or "NA"),
        name=name,
        variant_type=pick("variant_type") or "NA",
        clinical_significance=pick("clinical_significance") or "NA",
        review_status=review_status_label(review),
        star_level=star_level(review),
        transcript=transcript,
        hgvs_c=hgvs_c,
        hgvs_p=hgvs_p,
        condition_text=pick("condition_text") or "NA",
        last_evaluated=pick("last_evaluated") or "NA",
        n_submitters=_as_int(pick("n_submitters")),
        assembly=pick("assembly") or "NA",
        chromosome=pick("chromosome") or "NA",
        position_vcf=pick("position_vcf") or "NA",
        ref_allele=pick("ref_allele") or "NA",
        alt_allele=pick("alt_allele") or "NA",
        rcv=pick("rcv") or "NA",
        origin=pick("origin") or "NA",
        change=parse_protein_hgvs(hgvs_p),
        raw=dict(row),
    )


def normalize_records(rows: Iterable[dict], gene: str | None = None) -> list[ClinVarRecord]:
    """Normalize a payload and guarantee globally unique ``record_id`` values.

    Record identity carries the whole four-stratum conservation proof, so a
    duplicated ClinVar VariationID is disambiguated deterministically rather than
    collapsed — collapsing happens later, at residue level, and only there.
    """
    records: list[ClinVarRecord] = []
    seen: dict[str, int] = {}
    for index, row in enumerate(rows):
        record = normalize_row(row, index, gene=gene)
        count = seen.get(record.record_id, 0)
        seen[record.record_id] = count + 1
        if count:
            record.record_id = f"{record.record_id}-{count + 1}"
        records.append(record)
    return records


# --- ``variant_summary.txt.gz`` parsing -------------------------------------

def iter_variant_summary(handle: Iterable[str]) -> Iterator[dict]:
    """Stream a tab-delimited ``variant_summary`` into dicts (header-driven)."""
    header: list[str] | None = None
    for line in handle:
        text = line.rstrip("\n").rstrip("\r")
        if not text:
            continue
        if header is None:
            header = text.lstrip("#").split("\t")
            header[0] = header[0].strip()
            continue
        values = text.split("\t")
        yield dict(zip(header, values))


def parse_variant_summary(payload: bytes | str, gene: str,
                          assembly: str = "GRCh38") -> list[dict]:
    """Select the gene's substitution records from a ``variant_summary`` payload.

    The only selection applied here is *gene*, *assembly* and *record type*.
    There is deliberately no significance filter and no review-status filter:
    strata 1 and 3 must be able to account for every retrieved record.
    """
    if isinstance(payload, bytes):
        head = payload[:2]
        raw = gzip.decompress(payload) if head == b"\x1f\x8b" else payload
        text = raw.decode("utf-8", errors="replace")
    else:
        text = payload

    wanted_gene = gene.strip().upper()
    wanted_assembly = assembly.strip().upper()
    out: list[dict] = []
    for row in iter_variant_summary(io.StringIO(text)):
        symbols = {s.strip().upper() for s in str(row.get("GeneSymbol", "")).split(";")}
        if wanted_gene not in symbols:
            continue
        row_assembly = str(row.get("Assembly", "")).strip().upper()
        if row_assembly and row_assembly != wanted_assembly:
            continue
        if str(row.get("Type", "")).strip().lower() not in SUBSTITUTION_TYPES:
            continue
        out.append(row)
    return out


def raw_tsv_text(rows: Sequence[dict]) -> str:
    """Serialize a raw payload for ``01_INPUT_RAW/clinvar_raw.tsv``.

    Column order follows first appearance in the payload, so the immutable copy
    reflects what the source delivered rather than what Stage A prefers. Used
    only when the source cannot hand over the original bytes.
    """
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(
            str(row.get(col, "")).replace("\t", " ").replace("\n", " ")
            for col in columns
        ))
    return "\n".join(lines) + "\n"


def transcript_matches(record_transcript: str | None, mane_transcript: str | None) -> bool:
    """Accession comparison ignoring version (versions are recorded, not compared)."""
    left = strip_version(record_transcript)
    right = strip_version(mane_transcript)
    return bool(left) and bool(right) and left == right
