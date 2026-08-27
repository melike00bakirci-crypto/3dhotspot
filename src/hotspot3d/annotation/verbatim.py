"""Verbatim copying of upstream statistics — agent §7.2, §9.1.

Every statistic in ``hotspot_annotation.tsv`` is *copied* from ``06_FINAL_HOTSPOTS/``.
Recomputation is prohibited, and a suspect number is escalated, never "corrected".

The guarantee here is structural rather than procedural: :func:`load_verbatim_table`
keeps the **raw text of every cell**, and :func:`verbatim_copy` returns that raw text.
An upstream p-value never becomes a Python float inside Stage E, so there is no
representation in which it *could* be recomputed, reformatted or rounded. The
:class:`VerbatimAudit` then re-reads what was actually written and compares it
byte-for-byte against the source file, recording the result in QC.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..utils.errors import BlockedError, EscalationRequired

#: Statistics that must be carried across byte-for-byte (agent §3, lead item 1).
#: These are Stage E's LOGICAL names. Stage B emits some of them under different
#: column names; see :data:`SOURCE_COLUMN_ALIASES`. Renaming a column is not
#: recomputation — the bytes of the value are unchanged — so the verbatim guarantee
#: is preserved. Z is deliberately absent: METHOD_SPEC II.7 defines the per-region
#: statistics as member centers, covered intervals, classified-variant content, FE,
#: min/median p_emp and min FDR-hat. Z is a PER-CENTER statistic and has no
#: region-level definition, so requiring it here would have forced Stage E either to
#: aggregate (i.e. author a statistic) or to block on a column that correctly does
#: not exist.
COPIED_STATISTICS = (
    "residue_start",
    "residue_end",
    "n_plp",
    "n_blb",
    "fold_enrichment",
    "p_emp",
    "q_bh",
)

#: Stage E logical name -> the column Stage B actually emits in hotspot_regions.tsv.
#: Declared in ONE place so a B-side rename fails here loudly and is repaired at a
#: single site rather than being guessed at each call.
SOURCE_COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "residue_start": ("residue_start", "residue_index_min"),
    "residue_end": ("residue_end", "residue_index_max"),
    "p_emp": ("p_emp", "min_p_emp"),
    "q_bh": ("q_bh", "min_q_bh"),
}


def resolve_source_column(logical: str, available: tuple[str, ...] | list[str]) -> str:
    """Map a logical statistic name onto the column the source actually provides."""
    for candidate in SOURCE_COLUMN_ALIASES.get(logical, (logical,)):
        if candidate in available:
            return candidate
    return logical


@dataclass(frozen=True)
class VerbatimTable:
    """A source table held as raw text, keyed by one column."""

    path: Path
    key_column: str
    columns: tuple[str, ...]
    rows: dict[str, dict[str, str]]

    def __contains__(self, key: object) -> bool:
        return key in self.rows

    def keys(self) -> list[str]:
        """Keys in file order — Stage E never re-sorts upstream output."""
        return list(self.rows)


def load_verbatim_table(path: str | Path, key_column: str) -> VerbatimTable:
    """Load a canonical TSV **without parsing any value**.

    Unlike :func:`hotspot3d.utils.io.read_tsv`, no cell is converted to int/float/bool.
    That is the point: values pass through Stage E as the exact bytes upstream wrote.
    """
    path = Path(path)
    if not path.is_file():
        raise BlockedError(f"BLOCKED — required upstream input not found: {path}")

    raw = path.read_text(encoding="utf-8").splitlines()
    if not raw:
        raise BlockedError(f"BLOCKED — empty file (no header): {path}")

    header = raw[0].split("\t")
    if key_column not in header:
        raise BlockedError(
            f"BLOCKED — {path.name} has no key column {key_column!r}; "
            f"columns are {header}."
        )

    rows: dict[str, dict[str, str]] = {}
    for lineno, line in enumerate(raw[1:], start=2):
        if not line:
            continue
        values = line.split("\t")
        if len(values) != len(header):
            raise BlockedError(
                f"BLOCKED — {path.name} line {lineno} has {len(values)} fields but the "
                f"header declares {len(header)}. Malformed upstream input is never "
                f"repaired downstream."
            )
        record = dict(zip(header, values))
        key = record[key_column]
        if key in rows:
            raise BlockedError(
                f"BLOCKED — duplicate {key_column}={key!r} in {path.name}. Stage E "
                f"cannot decide which upstream row is authoritative; escalate."
            )
        rows[key] = record

    return VerbatimTable(
        path=path, key_column=key_column, columns=tuple(header), rows=rows
    )


def verbatim_copy(source_table: VerbatimTable, key: str, field_name: str) -> str:
    """Return the raw upstream cell for ``key``/``field_name``.

    Never parses, never reformats, never fills a default. A missing key or column is
    an escalation, not a gap to be papered over: Stage E is annotating a region that
    upstream says exists, so an absent statistic means the two disagree about what
    was produced.
    """
    if key not in source_table.rows:
        raise EscalationRequired(
            ambiguity=(
                f"Hotspot {key!r} is being annotated but has no row in "
                f"{source_table.path.name}."
            ),
            options=[
                "Stop and have the Lead re-verify Stage B output integrity.",
                "Drop the hotspot from the annotation.",
            ],
            consequences=[
                "Stage E stops; nothing is published from an inconsistent input set.",
                "A statistically significant region silently disappears from the "
                "report — an undetectable omission.",
            ],
            recommendation=(
                "Stop. A hotspot present in one Stage B artifact and absent from "
                "another is an upstream integrity problem that Stage E must not mask."
            ),
        )

    record = source_table.rows[key]
    # Resolve the logical statistic name onto the column the source actually emits.
    # This is a rename, not a computation: the returned bytes are the upstream cell.
    field_name = resolve_source_column(field_name, tuple(record))
    if field_name not in record:
        raise EscalationRequired(
            ambiguity=(
                f"Statistic {field_name!r} is required for hotspot {key!r} but "
                f"{source_table.path.name} does not provide that column "
                f"(has {list(record)})."
            ),
            options=[
                "Stop and have Stage B emit the missing statistic.",
                "Recompute the statistic in Stage E.",
            ],
            consequences=[
                "Stage E stops; the contract between B and E is repaired at source.",
                "PROHIBITED — Stage E recomputing any upstream statistic breaks the "
                "verbatim-copy guarantee and the audit that depends on it.",
            ],
            recommendation=(
                "Stop. Recomputation is prohibited by agent §16; the column must come "
                "from Stage B."
            ),
        )
    return record[field_name]


@dataclass
class VerbatimAudit:
    """Records every copy, then proves the written file matches the source.

    QC rule §9.1 requires an *automated* byte-for-byte comparison whose result is
    recorded — not an assertion that copying was intended.
    """

    source: VerbatimTable
    entries: list[tuple[str, str, str]] = field(default_factory=list)

    def copy(self, key: str, field_name: str) -> str:
        value = verbatim_copy(self.source, key, field_name)
        self.entries.append((key, field_name, value))
        return value

    def copy_all(self, key: str, fields=COPIED_STATISTICS) -> dict[str, str]:
        return {f: self.copy(key, f) for f in fields}

    # -- post-write verification -------------------------------------------
    def verify_written(self, written_path: str | Path, key_column: str = "hotspot_id") -> dict:
        """Re-read the emitted TSV and compare every copied cell to the source.

        Returns a QC record. Any mismatch raises :class:`EscalationRequired` — Stage E
        never "corrects" a discrepancy it discovers in its own output.
        """
        written = load_verbatim_table(written_path, key_column)
        mismatches: list[dict] = []

        for key, field_name, expected in self.entries:
            actual = written.rows.get(key, {}).get(field_name)
            if actual != expected:
                mismatches.append(
                    {
                        "hotspot_id": key,
                        "field": field_name,
                        "source_value": expected,
                        "written_value": "<ABSENT>" if actual is None else actual,
                    }
                )

        record = {
            "check": "verbatim_copy_byte_for_byte",
            "source_file": str(self.source.path.name),
            "written_file": str(Path(written_path).name),
            "n_cells_compared": len(self.entries),
            "n_mismatches": len(mismatches),
            "result": "PASS" if not mismatches else "FAIL",
            "mismatches": mismatches,
        }
        if mismatches:
            raise EscalationRequired(
                ambiguity=(
                    f"{len(mismatches)} copied statistic(s) in "
                    f"{Path(written_path).name} do not match "
                    f"{self.source.path.name} byte-for-byte: {mismatches[:3]}"
                ),
                options=[
                    "Stop; the Lead re-verifies Stage B artifacts and Stage E code.",
                    "Overwrite the annotation with the Stage B values.",
                ],
                consequences=[
                    "No annotation is published until the discrepancy is understood.",
                    "PROHIBITED — silently rewriting hides whichever of the two "
                    "stages is actually wrong.",
                ],
                recommendation=(
                    "Stop and escalate. A verbatim-copy mismatch is a declared "
                    "BLOCKING failure (agent §11) and is never self-resolved."
                ),
            )
        return record
