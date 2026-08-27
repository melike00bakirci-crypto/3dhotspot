"""ClinVar review-status vocabulary — RECORDED METADATA ONLY (F12).

Review status and its star level are derived here so that they can be *reported*:
``review_star_distribution.tsv``, the ``star_levels``/``max_star`` columns, and the
explicitly separate >=1*/>=2* sensitivity channel of METHOD_SPEC II.7.

They are never an inclusion criterion. Nothing in this module makes a decision:
every function is a total lookup, so there is no branch anywhere in Stage A that
can key on a star. :mod:`hotspot3d.data.audit` asserts that mechanically.
"""
from __future__ import annotations

# Canonical ClinVar review-status -> star level. Historical and current wordings
# are both present; ClinVar renamed "interpretations" to "classifications" in 2024
# and both spellings still occur across weekly releases.
REVIEW_STATUS_STARS: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, multiple submitters": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, conflicting classifications": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no interpretation for the single variant": 0,
    "no classification for the single variant": 0,
    "no classifications from unflagged records": 0,
    "no interpretations from unflagged records": 0,
    "flagged submission": 0,
}

STAR_LEVELS: tuple[int, ...] = (0, 1, 2, 3, 4)
UNKNOWN_REVIEW_STATUS = "NA"


def normalize_review_status(text: str | None) -> str:
    """Lowercase, collapse internal whitespace. Total: ``None`` -> ``''``."""
    return " ".join(str(text or "").strip().lower().split())


def star_level(review_status: str | None) -> int | None:
    """Star level for a review status, or ``None`` when the wording is unknown.

    An unknown wording is recorded as ``NA`` in the distribution table. It is
    never mapped to 0, because inventing a star would be inventing metadata.
    """
    return REVIEW_STATUS_STARS.get(normalize_review_status(review_status))


def review_status_label(text: str | None) -> str:
    """The review status as it will appear in output tables (verbatim, or ``NA``)."""
    return str(text or "").strip() or UNKNOWN_REVIEW_STATUS


# --- presentation helpers ---------------------------------------------------
# These take a bare value rather than reaching into a record, so no caller ever
# needs a branch whose test mentions star vocabulary. That keeps the F12 AST
# audit in hotspot3d.data.audit sharp instead of riddled with exemptions.

def star_or_na(value: int | None) -> object:
    """Render a star level for a table cell. Total: ``None`` -> ``'NA'``."""
    return {None: UNKNOWN_REVIEW_STATUS}.get(value, value)


def star_sort_key(value: int | None) -> int:
    """Sort ordinal for a star level; unknown sorts first. Presentation only."""
    return {None: -1}.get(value, value)
