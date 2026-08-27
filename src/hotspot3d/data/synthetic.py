"""Deterministic synthetic Stage A fixtures — the team's integration substrate.

Canonical import path::

    from hotspot3d.data.synthetic import make_synthetic_case, make_synthetic_source
    case = make_synthetic_case("clustered")     # -> SyntheticCase(records, structure, expected)

``tests/fixtures/synthetic_clinvar.py`` re-exports these under the path named in
the build plan; both are the same objects.

Everything here is a pure function of the case name: coordinates come from a
fixed-seed PCG64 generator, records are generated in a fixed order, and no clock,
filesystem or network is touched. Two calls return equal data, on any machine.

The scaffold is an explicitly synthetic **antiparallel helix bundle**, not a
predicted fold, but it reproduces the geometry the analysis depends on:

  * consecutive CA-CA distances are exactly 3.8 A everywhere, inside helices and
    across loops alike;
  * helices pack ~11.5 A apart, so residues on different helices sit ~7-12 A
    apart — the range real CA contacts occupy;
  * a spatial neighbourhood therefore draws from SEVERAL SEQUENCE SEGMENTS, which
    is the entire reason this analysis is done in 3D rather than on the sequence.

That last property is load-bearing downstream. A cluster built as one contiguous
sequence run has every MST edge at 3.8 A, so ``rho_all = 1.9`` and
``rho_max = 1.25 * rho_all = 2.375`` falls below the frozen ``rho_min = 3.0`` —
the II.9 footprint domain comes out empty and Stage C terminates before it can
exercise anything. Multi-segment clusters put the maximum MST edge comfortably
above 4.8 A, which is the condition every case intended to reach Stage C/D must
meet. ``expected["max_mst_edge_plp_A"]`` records it per case.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Iterable, Sequence

import numpy as np

from ..structure.sources import Atom, StructureModel, StructureResidue
from ..utils.geometry import cross_distances, mst_edges
from .hgvs import AA3_TO_1
from .sources import GeneResolution

CA_SPACING_A = 3.8
CB_BOND_A = 1.53
CG_BOND_A = 2.5
MASTER_FIXTURE_SEED = 20250101

# --- helix-bundle scaffold geometry -----------------------------------------
#: Alpha-helical rise per residue and turn per residue.
HELIX_RISE_A = 1.5
HELIX_TURN_DEG = 100.0
#: Radius chosen so that the helical step is EXACTLY CA_SPACING_A:
#:     chord^2 + rise^2 = 3.8^2,  chord = 2 R sin(turn/2)
HELIX_RADIUS_A = (math.sqrt(CA_SPACING_A ** 2 - HELIX_RISE_A ** 2)
                  / (2.0 * math.sin(math.radians(HELIX_TURN_DEG) / 2.0)))
#: Residues per helix, and axis-to-axis spacing in the bundle.
HELIX_LEN = 18
BUNDLE_PITCH_A = 11.5
BUNDLE_COLS = 3
#: Phase offset per helix, so facing residues are not all in register.
HELIX_PHASE_STEP_DEG = 37.0
#: Direction jitter, applied before each step is renormalized back to exactly
#: 3.8 A. Breaks the exact pairwise-distance ties a perfect helix would produce;
#: real structures have none, and downstream tie-breaking should not be exercised
#: by an artefact of the fixture.
JITTER_A = 0.30
#: A case meant to reach Stage C/D needs max MST edge > this, so that
#: rho_all > 2.4 and rho_max = 1.25 * rho_all clears the frozen rho_min = 3.0.
MIN_MST_EDGE_FOR_FOOTPRINT_A = 4.8

GENE = "SYNGENE"
UNIPROT_ACC = "Q9SYN1"
MANE_TRANSCRIPT = "NM_900001.4"
ALT_TRANSCRIPT = "NM_900777.2"
CLINVAR_RELEASE = "2025-01-06"
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA1_TO_3 = {one: three.capitalize() for three, one in AA3_TO_1.items()}

#: Cycled so that every star level, including 0, is represented in every case.
REVIEW_STATUSES: tuple[str, ...] = (
    "no assertion criteria provided",                        # 0 stars
    "criteria provided, single submitter",                   # 1 star
    "criteria provided, multiple submitters, no conflicts",  # 2 stars
    "reviewed by expert panel",                              # 3 stars
    "no classification provided",                            # 0 stars
    "practice guideline",                                    # 4 stars
)
PLP_LABELS = ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")
BLB_LABELS = ("Benign", "Likely benign", "Benign/Likely benign")

# CASE_NAMES is derived from CASE_SPECS below, so the two can never drift apart.


# --- sequence and geometry --------------------------------------------------

def synthetic_sequence(length: int) -> str:
    """Deterministic canonical sequence; includes glycines (no CB) on purpose."""
    return "".join(AA_ALPHABET[(i * 7 + 3) % len(AA_ALPHABET)] for i in range(length))


def _perpendicular(direction: np.ndarray) -> np.ndarray:
    """A deterministic unit vector orthogonal to ``direction``."""
    axis = np.zeros(3)
    axis[int(np.argmin(np.abs(direction)))] = 1.0
    perpendicular = np.cross(direction, axis)
    return perpendicular / np.linalg.norm(perpendicular)


def _loop_points(start: np.ndarray, end: np.ndarray) -> list[np.ndarray]:
    """Loop residues from ``start`` to ``end``, every step exactly 3.8 A.

    A planar zig-zag: each step advances ``D/n`` along the start-end axis and
    flips a fixed perpendicular offset, so its length is
    ``sqrt(axial^2 + perp^2) = 3.8`` by construction. ``n`` is even, so the
    perpendicular offsets cancel and the walk lands exactly on ``end``.
    """
    delta = end - start
    span = float(np.linalg.norm(delta))
    n = max(2, int(math.ceil(span / CA_SPACING_A)))
    n += n % 2                                    # even: the zig-zag closes on ``end``

    axial = span / n
    perpendicular_length = math.sqrt(max(CA_SPACING_A ** 2 - axial ** 2, 0.0))
    along = delta / span
    across = _perpendicular(along)
    return [start + along * (k * axial) + across * (perpendicular_length * (k % 2))
            for k in range(1, n)]


@dataclass(frozen=True)
class FoldLayout:
    """Nominal helix-bundle trace plus the sequence extent of each helix."""

    coords: np.ndarray
    helices: tuple[tuple[int, int], ...]          # 1-based inclusive residue ranges

    def segment_of(self) -> np.ndarray:
        """Segment id per residue (0-based index); loops join the helix before them."""
        ids = np.zeros(len(self.coords), dtype=int)
        for segment, (first, last) in enumerate(self.helices):
            ids[first - 1:last] = segment
            if segment + 1 < len(self.helices):
                ids[last:self.helices[segment + 1][0] - 1] = segment
        return ids

    @property
    def n_helices(self) -> int:
        return len(self.helices)


def fold_layout(length: int, *, helix_len: int = HELIX_LEN, n_cols: int = BUNDLE_COLS,
                pitch: float = BUNDLE_PITCH_A) -> FoldLayout:
    """Nominal antiparallel helix bundle: every consecutive step is exactly 3.8 A."""
    height = (helix_len - 1) * HELIX_RISE_A
    turn = math.radians(HELIX_TURN_DEG)
    points: list[np.ndarray] = []
    helices: list[tuple[int, int]] = []
    h = 0

    while len(points) < length:
        row, col = divmod(h, n_cols)
        if row % 2:                               # serpentine: consecutive helices adjacent
            col = n_cols - 1 - col
        centre_x, centre_y = col * pitch, row * pitch
        direction = 1 if h % 2 == 0 else -1
        phase = math.radians(HELIX_PHASE_STEP_DEG * h)

        helix = []
        for j in range(helix_len):
            angle = phase + direction * turn * j
            z = j * HELIX_RISE_A if direction > 0 else height - j * HELIX_RISE_A
            helix.append(np.array([centre_x + HELIX_RADIUS_A * math.cos(angle),
                                   centre_y + HELIX_RADIUS_A * math.sin(angle), z]))

        if points:
            points.extend(_loop_points(points[-1], helix[0]))
        first = len(points) + 1
        points.extend(helix)
        helices.append((first, len(points)))
        h += 1

    coords = np.asarray(points[:length], dtype=np.float64)
    helices = tuple((a, min(b, length)) for a, b in helices if a <= length)
    return FoldLayout(coords=coords, helices=helices)


def backbone_coordinates(length: int, seed: int, *, helix_len: int = HELIX_LEN,
                         n_cols: int = BUNDLE_COLS,
                         pitch: float = BUNDLE_PITCH_A) -> np.ndarray:
    """Jittered CA trace with exactly ``CA_SPACING_A`` between consecutive residues.

    Each step re-targets the nominal fold before it is renormalized, so the
    jitter perturbs the trace without letting it drift away from the bundle.
    """
    nominal = fold_layout(length, helix_len=helix_len, n_cols=n_cols, pitch=pitch).coords
    rng = np.random.Generator(np.random.PCG64(seed))
    jitter = rng.uniform(-JITTER_A, JITTER_A, size=nominal.shape)

    coords = np.zeros_like(nominal)
    coords[0] = nominal[0]
    for i in range(1, len(nominal)):
        step = (nominal[i] - coords[i - 1]) + jitter[i]
        norm = float(np.linalg.norm(step))
        if norm == 0.0:                                       # pragma: no cover
            step, norm = np.array([1.0, 0.0, 0.0]), 1.0
        coords[i] = coords[i - 1] + CA_SPACING_A * step / norm
    return coords


# --- geometric diagnostics recorded in ``expected`` -------------------------

def max_mst_edge(coords: np.ndarray, residues: Sequence[int]) -> float:
    """Longest edge of the Euclidean MST over ``residues`` (1-based indices).

    ``rho_all = max_edge / 2`` in METHOD_SPEC II.9, so this is what decides
    whether the footprint domain ``[3.0, 1.25 * rho_all]`` is non-empty.
    """
    if len(residues) < 2:
        return 0.0
    points = np.asarray([coords[i - 1] for i in sorted(residues)])
    return max(edge[2] for edge in mst_edges(points))


def sequence_segments(residues: Sequence[int]) -> list[tuple[int, int]]:
    """Contiguous runs of residue indices — how many sequence stretches a set spans."""
    ordered = sorted(residues)
    if not ordered:
        return []
    runs, start, previous = [], ordered[0], ordered[0]
    for index in ordered[1:]:
        if index != previous + 1:
            runs.append((start, previous))
            start = index
        previous = index
    runs.append((start, previous))
    return runs


def _plddt_values(length: int, seed: int, low_tail_from: int | None) -> list[float]:
    """High-confidence body, optional long low-confidence tail (case ``low_plddt``)."""
    rng = np.random.Generator(np.random.PCG64(seed))
    values = []
    for i in range(length):
        if low_tail_from is not None and i >= low_tail_from:
            span = max(length - low_tail_from - 1, 1)
            decay = (i - low_tail_from) / span
            base = 68.0 - 38.0 * decay
            values.append(round(float(base + rng.uniform(-2.5, 2.5)), 2))
        else:
            values.append(round(float(rng.uniform(82.0, 97.0)), 2))
    return values


def build_structure(sequence: str, *, seed: int, plddt: Sequence[float],
                    accession: str = UNIPROT_ACC, coords: np.ndarray | None = None,
                    n_cols: int = BUNDLE_COLS,
                    pitch: float = BUNDLE_PITCH_A) -> StructureModel:
    """Assemble a :class:`StructureModel` on the CA scaffold.

    CB and a single CG stand-in are placed at deterministic offsets so that the
    CB and side-chain-centroid columns exercise their real code path. They are
    marked ``unused_under_F1`` on every emitted row; glycine carries neither, so
    the ``NA`` path is exercised too.
    """
    if coords is None:
        coords = backbone_coordinates(len(sequence), seed, n_cols=n_cols, pitch=pitch)
    rng = np.random.Generator(np.random.PCG64(seed + 1))
    directions = rng.normal(size=(len(sequence), 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)

    residues: list[StructureResidue] = []
    for i, aa in enumerate(sequence):
        index = i + 1
        confidence = float(plddt[i])
        x, y, z = coords[i]
        atoms = {"CA": Atom("CA", float(x), float(y), float(z), 1.0, confidence, "C")}
        if aa != "G":
            dx, dy, dz = directions[i]
            atoms["CB"] = Atom("CB", float(x + CB_BOND_A * dx), float(y + CB_BOND_A * dy),
                               float(z + CB_BOND_A * dz), 1.0, confidence, "C")
            atoms["CG"] = Atom("CG", float(x + CG_BOND_A * dx), float(y + CG_BOND_A * dy),
                               float(z + CG_BOND_A * dz), 1.0, confidence, "C")
        residues.append(StructureResidue(residue_index=index,
                                         aa3=AA1_TO_3[aa].upper(),
                                         chain_id="A", atoms=atoms))

    return StructureModel(
        accession=accession, entry_id=f"AF-{accession}-F1-model_v4",
        model_version="4", db_version="2024-06-01", residues=residues,
        n_fragments=1, raw_cif_text=None,
        metadata={"source": "synthetic_fixture", "scaffold": "antiparallel_helix_bundle",
                  "ca_spacing_A": CA_SPACING_A, "helix_len": HELIX_LEN,
                  "bundle_pitch_A": pitch, "bundle_cols": n_cols},
    )


# --- record construction ----------------------------------------------------

def _alt_aa(ref: str, salt: int) -> str:
    alt = AA_ALPHABET[(AA_ALPHABET.index(ref) + 5 + salt) % len(AA_ALPHABET)]
    return alt if alt != ref else AA_ALPHABET[(AA_ALPHABET.index(ref) + 1) % 20]


def _record(counter: list[int], *, position: int, ref_aa: str, significance: str,
            review_status: str, transcript: str = MANE_TRANSCRIPT,
            hgvs_p: str | None = None, salt: int = 0,
            condition: str = "Synthetic condition") -> dict:
    counter[0] += 1
    n = counter[0]
    alt_aa = _alt_aa(ref_aa, salt)
    protein = hgvs_p or f"p.{AA1_TO_3[ref_aa]}{position}{AA1_TO_3[alt_aa]}"
    alt_base = "ACGT"[(AA_ALPHABET.index(alt_aa) + salt) % 4]
    coding = f"c.{position * 3 - 2}A>{alt_base if alt_base != 'A' else 'G'}"
    return {
        "VariationID": str(900000 + n),
        "AlleleID": str(800000 + n),
        "Type": "single nucleotide variant",
        "Name": f"{transcript}({GENE}):{coding} ({protein})",
        "GeneSymbol": GENE,
        "ClinicalSignificance": significance,
        "ClinSigSimple": "1" if significance in PLP_LABELS else "0",
        "LastEvaluated": "2024-06-01",
        "ReviewStatus": review_status,
        "NumberSubmitters": str(1 + (n % 4)),
        "PhenotypeList": condition,
        "Assembly": "GRCh38",
        "Chromosome": "1",
        "Start": str(1000000 + position * 3),
        "PositionVCF": str(1000000 + position * 3),
        "ReferenceAllele": "A",
        "AlternateAllele": "G",
        "ReferenceAlleleVCF": "A",
        "AlternateAlleleVCF": "G",
        "RCVaccession": f"RCV{9000000 + n}",
        "Origin": "germline",
        "SubmitterCategories": "1",
    }


def _noise_records(counter: list[int], sequence: str, *, plp_anchor: int,
                   spare: Sequence[int]) -> list[tuple[dict, str]]:
    """Records that must NOT enter the binary classes, each with its reason.

    One of each exclusion path, so every branch of the strata is exercised in
    every case and the expected counts are derived rather than hand-tallied.
    """
    def aa(pos: int) -> str:
        return sequence[pos - 1]

    def wrong_aa(pos: int) -> str:
        """A reference AA guaranteed to disagree with the canonical sequence."""
        return AA_ALPHABET[(AA_ALPHABET.index(aa(pos)) + 9) % len(AA_ALPHABET)]

    a, b, c, d, e, f = spare[:6]

    out: list[tuple[dict, str]] = [
        # significance categories that are preserved but never binary
        (_record(counter, position=b, ref_aa=aa(b), significance="Uncertain significance",
                 review_status=REVIEW_STATUSES[1]), "vus"),
        (_record(counter, position=plp_anchor, ref_aa=aa(plp_anchor),
                 significance="Uncertain significance",
                 review_status=REVIEW_STATUSES[3], salt=3), "vus"),
        (_record(counter, position=c, ref_aa=aa(c),
                 significance="Conflicting classifications of pathogenicity",
                 review_status="criteria provided, conflicting classifications"),
         "conflicting"),
        (_record(counter, position=d, ref_aa=aa(d), significance="drug response",
                 review_status=REVIEW_STATUSES[0]), "non_binary_significance"),
        # reference-AA mismatch on the MANE transcript
        (_record(counter, position=a, ref_aa=wrong_aa(a), significance="Pathogenic",
                 review_status=REVIEW_STATUSES[2]), "ref_aa_mismatch"),
        # non-MANE transcript, reference AA inconsistent -> not remappable
        (_record(counter, position=e, ref_aa=wrong_aa(e), significance="Likely pathogenic",
                 review_status=REVIEW_STATUSES[1], transcript=ALT_TRANSCRIPT),
         "non_canonical_transcript"),
        # protein changes that are not missense — all three parse CLEANLY into a
        # non-missense consequence (nonsense / synonymous / frameshift), so all
        # three are "non_missense_variant_type", never "unparseable_hgvs"
        # (Workflow v2 §1, Lead ruling 2026-08-18).
        (_record(counter, position=f, ref_aa=aa(f), significance="Pathogenic",
                 review_status=REVIEW_STATUSES[2],
                 hgvs_p=f"p.{AA1_TO_3[aa(f)]}{f}Ter"), "non_missense_variant_type"),
        (_record(counter, position=f, ref_aa=aa(f), significance="Likely benign",
                 review_status=REVIEW_STATUSES[0],
                 hgvs_p=f"p.{AA1_TO_3[aa(f)]}{f}="), "non_missense_variant_type"),
        (_record(counter, position=f, ref_aa=aa(f), significance="Pathogenic",
                 review_status=REVIEW_STATUSES[1],
                 hgvs_p=f"p.{AA1_TO_3[aa(f)]}{f}fs"), "non_missense_variant_type"),
    ]
    return out


def _class_records(counter: list[int], sequence: str, positions: Sequence[int],
                   labels: Sequence[str], offset: int = 0) -> list[dict]:
    """One or two records per residue, cycling labels and review statuses."""
    rows: list[dict] = []
    for k, position in enumerate(positions):
        ref = sequence[position - 1]
        rows.append(_record(counter, position=position, ref_aa=ref,
                            significance=labels[(k + offset) % len(labels)],
                            review_status=REVIEW_STATUSES[(k + offset)
                                                          % len(REVIEW_STATUSES)]))
        if k % 3 == 0:
            rows.append(_record(counter, position=position, ref_aa=ref,
                                significance=labels[(k + offset + 1) % len(labels)],
                                review_status=REVIEW_STATUSES[(k + offset + 2)
                                                              % len(REVIEW_STATUSES)],
                                salt=1))
    return rows


# --- the case object --------------------------------------------------------

@dataclass
class SyntheticCase:
    """One deterministic Stage A scenario: records, structure, expected truth."""

    records: list[dict]
    structure: StructureModel
    expected: dict
    name: str = "unnamed"
    gene: str = GENE
    uniprot_acc: str = UNIPROT_ACC
    mane_transcript: str = MANE_TRANSCRIPT
    sequence: str = ""
    clinvar_release: str = CLINVAR_RELEASE
    # Workflow v2 §2 — override for the three oligomeric-state test branches.
    # Defaults to "unknown" (no subunit comment at all), matching a fixture that
    # says nothing about the biological assembly.
    oligomeric_state: str = "unknown"
    oligomeric_state_evidence: str | None = None
    oligomeric_state_matched_phrase: str | None = None

    # -- provider ------------------------------------------------------------
    @property
    def source(self) -> "SyntheticStageASource":
        return SyntheticStageASource(self)

    def resolution(self) -> GeneResolution:
        return GeneResolution(
            gene=self.gene, uniprot_acc=self.uniprot_acc,
            uniprot_entry_version="12", canonical_isoform=f"{self.uniprot_acc}-1",
            sequence=self.sequence, mane_transcript=self.mane_transcript,
            mane_source="MANE_SELECT", mane_fallback_used=False,
            entry_name="SYN1_HUMAN", protein_name="Synthetic fixture protein",
            oligomeric_state=self.oligomeric_state,
            oligomeric_state_evidence=self.oligomeric_state_evidence,
            oligomeric_state_matched_phrase=self.oligomeric_state_matched_phrase,
        )

    def with_oligomeric_state(self, state: str, *, evidence: str | None = None,
                              matched_phrase: str | None = None) -> "SyntheticCase":
        """Workflow v2 §2 test fixture — override the UniProt subunit reading."""
        return replace(self, oligomeric_state=state, oligomeric_state_evidence=evidence,
                       oligomeric_state_matched_phrase=matched_phrase,
                       name=f"{self.name}_oligomeric_{state}")

    # -- deterministic variations used by tests ------------------------------
    def without_class(self, klass: str) -> "SyntheticCase":
        """Drop every record of one binary class — the empty-class negative branch."""
        labels = PLP_LABELS if klass == "PLP" else BLB_LABELS
        kept = [r for r in self.records if r["ClinicalSignificance"] not in labels]
        expected = dict(self.expected)
        expected.update({"N_P": 0 if klass == "PLP" else self.expected["N_P"],
                         "N_B": 0 if klass == "BLB" else self.expected["N_B"],
                         "empty_class": klass})
        return replace(self, records=kept, expected=expected,
                       name=f"{self.name}_no_{klass.lower()}")

    def with_permuted_stars(self, seed: int = 7) -> "SyntheticCase":
        """Reassign every review status at random — the F12 invariance fixture.

        The cohort must come out bit-identical: stars are metadata, so permuting
        them may not move a single record between strata.
        """
        rng = np.random.Generator(np.random.PCG64(seed))
        shuffled = []
        for record in self.records:
            row = dict(record)
            row["ReviewStatus"] = REVIEW_STATUSES[int(rng.integers(len(REVIEW_STATUSES)))]
            shuffled.append(row)
        return replace(self, records=shuffled, name=f"{self.name}_stars_permuted")

    def with_unusable_ca(self, residue_index: int) -> "SyntheticCase":
        """Zero the occupancy of one CA — the ``no_ca_coordinate`` branch."""
        residues = []
        for residue in self.structure.residues:
            if residue.residue_index != residue_index:
                residues.append(residue)
                continue
            atoms = dict(residue.atoms)
            ca = atoms["CA"]
            atoms["CA"] = Atom(ca.name, ca.x, ca.y, ca.z, 0.0, ca.b_factor, ca.element)
            residues.append(StructureResidue(residue.residue_index, residue.aa3,
                                             residue.chain_id, atoms))
        structure = replace(self.structure, residues=residues)
        expected = dict(self.expected)
        expected["unusable_ca_residue"] = residue_index
        return replace(self, structure=structure, expected=expected,
                       name=f"{self.name}_unusable_ca")

    def with_fragments(self, n_fragments: int) -> "SyntheticCase":
        """Declare a multi-fragment AFDB entry — the BLOCKED branch of II.13."""
        return replace(self, structure=replace(self.structure, n_fragments=n_fragments),
                       name=f"{self.name}_fragments{n_fragments}")


class SyntheticStageASource:
    """A :class:`hotspot3d.data.sources.StageASource` backed by a fixture.

    Satisfies the variant, sequence and structure interfaces at once, so
    ``run_stage_a(ctx, gene=..., source=case.source)`` exercises exactly the code
    path a live run would, minus the network.
    """

    def __init__(self, case: SyntheticCase) -> None:
        self.case = case
        self.calls: list[str] = []

    # -- variants ------------------------------------------------------------
    def fetch_missense(self, gene: str) -> list[dict]:
        self.calls.append(f"fetch_missense({gene})")
        return [dict(record) for record in self.case.records]

    def release_metadata(self) -> dict:
        return {
            "source": "synthetic_clinvar_fixture",
            "release_date": self.case.clinvar_release,
            "query": (f"GeneSymbol=={self.case.gene} AND Assembly==GRCh38 AND "
                      f"Type==missense AND review_star_filter=NONE"),
            "url": "NA",
            "assembly": "GRCh38",
            "review_star_filter": None,
            "use_review_stars_for_inclusion": False,
            "synthetic_case": self.case.name,
            "http_status": "NA",
            "retries": 0,
        }

    # -- sequence ------------------------------------------------------------
    def resolve_gene(self, gene: str) -> GeneResolution:
        self.calls.append(f"resolve_gene({gene})")
        return self.case.resolution()

    def sequence_metadata(self) -> dict:
        return {
            "source": "synthetic_uniprot_fixture",
            "release_date": "2024-06-12",
            "query": f"gene_exact:{self.case.gene} AND organism_id:9606 AND reviewed:true",
            "url": "NA",
            "entry_version": "12",
            "synthetic_case": self.case.name,
        }

    # -- structure -----------------------------------------------------------
    def fetch_model(self, uniprot_acc: str) -> StructureModel:
        self.calls.append(f"fetch_model({uniprot_acc})")
        return self.case.structure

    def model_metadata(self) -> dict:
        return {
            "source": "synthetic_alphafold_fixture",
            "release_date": self.case.structure.db_version,
            "query": f"AF-{self.case.uniprot_acc}-F1-model_v4",
            "url": "NA",
            "model_version": self.case.structure.model_version,
            "n_fragments": self.case.structure.n_fragments,
            "synthetic_case": self.case.name,
        }

    def raw_payload(self, kind: str) -> bytes | None:
        return None


# --- case builders ----------------------------------------------------------

def _spatially_nearest(coords: np.ndarray, center_index: int, k: int,
                       exclude: Iterable[int] = ()) -> list[int]:
    """The ``k`` residues nearest a center, by CA-CA Euclidean distance (1-based)."""
    distances = cross_distances(coords[center_index - 1][None, :], coords)[0]
    order = np.argsort(distances, kind="stable")
    blocked = set(exclude)
    picked: list[int] = []
    for position in order:
        residue = int(position) + 1
        if residue in blocked:
            continue
        picked.append(residue)
        if len(picked) == k:
            break
    return sorted(picked)


def _spread(length: int, k: int, seed: int, exclude: Iterable[int] = ()) -> list[int]:
    """``k`` residues drawn without replacement from the whole protein."""
    rng = np.random.Generator(np.random.PCG64(seed))
    blocked = set(exclude)
    pool = [i for i in range(1, length + 1) if i not in blocked]
    chosen = rng.choice(np.asarray(pool), size=k, replace=False)
    return sorted(int(c) for c in chosen)


def _quotas(k: int, parts: int) -> list[int]:
    """Split ``k`` into ``parts`` as evenly as possible, front-loaded."""
    base, extra = divmod(k, parts)
    return [base + (1 if i < extra else 0) for i in range(parts)]


def _cluster_across_segments(coords: np.ndarray, segment_of: np.ndarray,
                             center_index: int, k: int, n_segments: int) -> list[int]:
    """A tight 3D neighbourhood drawn from ``n_segments`` distinct sequence stretches.

    This is what makes a synthetic hotspot a *3D* hotspot. Residues are ranked by
    CA-CA distance to the center, the ``n_segments`` nearest helices are chosen,
    and each contributes its share of the nearest residues. The resulting MST
    therefore carries at least one inter-helix edge (~7-12 A), which is what keeps
    ``rho_all`` above 2.4 A and the METHOD_SPEC II.9 footprint domain non-empty.
    """
    distances = cross_distances(coords[center_index - 1][None, :], coords)[0]
    order = np.argsort(distances, kind="stable")

    nearest_by_segment: dict[int, float] = {}
    for position in order:
        nearest_by_segment.setdefault(int(segment_of[position]), float(distances[position]))
    chosen = sorted(nearest_by_segment,
                    key=lambda s: (nearest_by_segment[s], s))[:n_segments]

    picked: list[int] = []
    for segment, quota in zip(chosen, _quotas(k, len(chosen))):
        members = [int(p) + 1 for p in order if int(segment_of[p]) == segment]
        picked.extend(members[:quota])
    return sorted(picked)


def _dispersed(coords: np.ndarray, k: int, seed: int) -> list[int]:
    """``k`` maximally separated residues (deterministic farthest-point sampling).

    Stronger than a uniform draw: a uniform sample can still leave a locally
    enriched neighbourhood by chance, whereas farthest-point sampling actively
    empties every neighbourhood. That is what "genuinely non-clustered" requires.
    """
    rng = np.random.Generator(np.random.PCG64(seed))
    start = int(rng.integers(len(coords)))
    chosen = [start]
    nearest = cross_distances(coords[start][None, :], coords)[0]
    for _ in range(k - 1):
        nxt = int(np.argmax(nearest))
        chosen.append(nxt)
        nearest = np.minimum(nearest, cross_distances(coords[nxt][None, :], coords)[0])
    return sorted(i + 1 for i in chosen)


def _outside_shell(coords: np.ndarray, center_index: int, k: int,
                   min_distance: float, seed: int,
                   exclude: Iterable[int] = ()) -> list[int]:
    """``k`` residues drawn from OUTSIDE a sphere around the cluster centre.

    Keeping B/LB clearly out of the pathogenic core is what makes a positive
    control decisive rather than marginal: a sphere over the core then contains
    P/LP and almost nothing else, so its enrichment is unambiguous instead of
    sitting on the BH boundary where a seed change can flip it.
    """
    distances = cross_distances(coords[center_index - 1][None, :], coords)[0]
    blocked = set(exclude)
    pool = [i + 1 for i in range(len(coords))
            if distances[i] >= min_distance and (i + 1) not in blocked]
    if len(pool) < k:
        raise ValueError(
            f"only {len(pool)} residues lie beyond {min_distance} A of residue "
            f"{center_index}; cannot place {k} B/LB outside the core")
    rng = np.random.Generator(np.random.PCG64(seed))
    chosen = rng.choice(np.asarray(pool), size=k, replace=False)
    return sorted(int(c) for c in chosen)


def _from_helices(layout: FoldLayout, helix_ids: Sequence[int], k: int) -> list[int]:
    """``k`` residues taken as contiguous runs from the middle of chosen helices."""
    picked: list[int] = []
    for helix_id, quota in zip(helix_ids, _quotas(k, len(helix_ids))):
        first, last = layout.helices[helix_id]
        span = last - first + 1
        start = first + max((span - quota) // 2, 0)
        picked.extend(range(start, min(start + quota, last + 1)))
    return sorted(picked)


#: How each case selects its P/LP residues.
SELECT_CLUSTER = "cluster_across_segments"     # tight 3D hotspot from N helices
SELECT_SPREAD = "uniform_spread"               # uniform random draw
SELECT_DISPERSED = "farthest_point"            # maximally dispersed
SELECT_HELICES = "whole_helices"               # contiguous runs from named helices


@dataclass(frozen=True)
class _CaseSpec:
    length: int
    n_plp: int
    n_blb: int
    n_conflict: int
    selection: str
    low_tail_from: int | None
    seed: int
    cluster_center: int | None = None
    n_segments: int = 3
    helix_ids: tuple[int, ...] = ()
    blb_helix_ids: tuple[int, ...] = ()
    #: When set, B/LB residues are drawn only from beyond this distance of the
    #: cluster centre, keeping the pathogenic core clean.
    blb_min_distance_A: float | None = None
    n_cols: int = BUNDLE_COLS
    pitch: float = BUNDLE_PITCH_A
    purpose: str = ""
    reaches_stage_cd: bool = False
    expected_negative_branch: str | None = None
    exercises: tuple[str, ...] = ()


CASE_SPECS: dict[str, _CaseSpec] = {
    # --- Stage A behaviour ---------------------------------------------------
    # A POSITIVE CONTROL MUST BE DECISIVE, not merely correct. The earlier
    # 180-residue / 24 P/LP version put per-center p-values right at the BH
    # boundary and left only two admissible radii, so min-max normalization over
    # the admissible set was degenerate (both candidates at distance 1.0) and the
    # result — r_hot AND |S| — flipped with the run_id-derived seed. That is a
    # fixture fault, not a pipeline fault: the methodology behaved correctly on
    # ambiguous data. The cohort is now large enough, the core dense enough and
    # B/LB far enough outside it that the answer is unambiguous.
    "clustered": _CaseSpec(
        length=260, n_plp=44, n_blb=90, n_conflict=0,
        selection=SELECT_CLUSTER, low_tail_from=None,
        seed=MASTER_FIXTURE_SEED + 1, cluster_center=130, n_segments=4,
        blb_min_distance_A=18.0,
        purpose=("the happy path: a decisive 3D hotspot drawn from four separate "
                 "sequence segments packed together in the bundle, with B/LB held "
                 "outside an 18 A shell so enrichment is unambiguous, so the whole "
                 "chain A -> B -> C -> D -> E can run and its result is stable "
                 "across run_ids"),
        reaches_stage_cd=True,
        exercises=("global_clustering", "r_hot_selection", "significant_centers",
                   "footprint_domain", "robustness"),
    ),
    "no_cluster": _CaseSpec(
        length=260, n_plp=44, n_blb=90, n_conflict=0,
        selection=SELECT_SPREAD, low_tail_from=None,
        seed=MASTER_FIXTURE_SEED + 2,
        purpose=("matched control for `clustered`: identical scaffold, identical "
                 "protein length and identical cohort sizes, with P/LP placed at "
                 "random instead of in a 3D cluster. The two differ in ONE "
                 "variable — the spatial arrangement of the labels — which is the "
                 "only way the comparison means anything"),
        expected_negative_branch="no_global_clustering_or_S_empty",
        exercises=("global_clustering_non_significant",),
    ),
    "sparse": _CaseSpec(
        length=120, n_plp=3, n_blb=2, n_conflict=0,
        selection=SELECT_SPREAD, low_tail_from=None,
        seed=MASTER_FIXTURE_SEED + 3,
        purpose="very few classified residues; exercises the underpowered path",
        expected_negative_branch="underpowered_cohort",
        exercises=("sparse_evidence_reporting",),
    ),
    "conflict": _CaseSpec(
        length=150, n_plp=12, n_blb=14, n_conflict=4,
        selection=SELECT_CLUSTER, low_tail_from=None,
        seed=MASTER_FIXTURE_SEED + 4, cluster_center=70, n_segments=3,
        purpose="residues carrying both P/LP and B/LB evidence (F3)",
        exercises=("residue_class_conflict",),
    ),
    "low_plddt": _CaseSpec(
        length=200, n_plp=18, n_blb=20, n_conflict=0,
        selection=SELECT_CLUSTER, low_tail_from=140,
        seed=MASTER_FIXTURE_SEED + 5, cluster_center=60, n_segments=3,
        purpose="a long low-confidence tail that must be recorded and never filter (F2)",
        exercises=("plddt_recorded_not_filtered", "plddt70_sensitivity"),
    ),

    # --- downstream negative and diagnostic branches -------------------------
    "no_significant_hotspot": _CaseSpec(
        length=180, n_plp=20, n_blb=40, n_conflict=0,
        selection=SELECT_DISPERSED, low_tail_from=None,
        seed=MASTER_FIXTURE_SEED + 6,
        purpose=("P/LP present in reasonable number but MAXIMALLY dispersed by "
                 "farthest-point sampling, so no residue center survives BH-FDR "
                 "at q=0.05 and S comes out empty"),
        expected_negative_branch="S_empty",
        exercises=("no_significant_hotspots", "valid_negative_at_stage_b"),
    ),
    "excessive_coverage": _CaseSpec(
        length=66, n_plp=22, n_blb=14, n_conflict=0,
        selection=SELECT_HELICES, low_tail_from=None,
        seed=MASTER_FIXTURE_SEED + 7, helix_ids=(0, 1), blb_helix_ids=(2,),
        n_cols=2,
        purpose=("a small compact protein (small D_max) whose P/LP residues occupy "
                 "two of its three helices, so significant centers are spread "
                 "through it and any usable rho already covers most of U_struct — "
                 "QC-F1 (coverage > 0.50) should reject every candidate"),
        reaches_stage_cd=True,
        expected_negative_branch="all_footprint_candidates_excessive_coverage",
        exercises=("footprint_qc_f1", "coverage_constraint_failure"),
    ),
    "permutation_resolution": _CaseSpec(
        length=620, n_plp=6, n_blb=45, n_conflict=0,
        selection=SELECT_CLUSTER, low_tail_from=None,
        seed=MASTER_FIXTURE_SEED + 8, cluster_center=300, n_segments=3,
        n_cols=6,
        purpose=("a large U_struct (m > q * (B + 1) = 500, so the k=1 BH threshold "
                 "q/m sits BELOW the resolution floor 1/(B+1)) with a handful of "
                 "extremely tight P/LP residues whose empirical p pins at that "
                 "floor — the II.4 conditions for PERMUTATION_RESOLUTION_LIMITED"),
        reaches_stage_cd=True,
        exercises=("permutation_resolution_diagnostic", "b_ladder_recommendation"),
    ),
}

CASE_NAMES: tuple[str, ...] = tuple(CASE_SPECS)


def _select_plp(spec: _CaseSpec, coords: np.ndarray, layout: FoldLayout) -> list[int]:
    """Place the P/LP residues according to the case's declared strategy."""
    if spec.selection == SELECT_CLUSTER:
        return _cluster_across_segments(coords, layout.segment_of(),
                                        spec.cluster_center, spec.n_plp,
                                        spec.n_segments)
    if spec.selection == SELECT_DISPERSED:
        return _dispersed(coords, spec.n_plp, spec.seed + 11)
    if spec.selection == SELECT_HELICES:
        return _from_helices(layout, spec.helix_ids, spec.n_plp)
    return _spread(spec.length, spec.n_plp, spec.seed + 11)


def make_synthetic_case(name: str) -> SyntheticCase:
    """Build one named case. Deterministic: equal inputs, equal outputs, always.

    Cases: ``clustered``, ``no_cluster``, ``sparse``, ``conflict``, ``low_plddt``,
    ``no_significant_hotspot``, ``excessive_coverage``, ``permutation_resolution``.
    """
    if name not in CASE_SPECS:
        raise KeyError(f"unknown synthetic case {name!r}; available: {list(CASE_SPECS)}")
    spec = CASE_SPECS[name]

    sequence = synthetic_sequence(spec.length)
    layout = fold_layout(spec.length, n_cols=spec.n_cols, pitch=spec.pitch)
    coords = backbone_coordinates(spec.length, spec.seed, n_cols=spec.n_cols,
                                  pitch=spec.pitch)
    plddt = _plddt_values(spec.length, spec.seed + 100, spec.low_tail_from)
    structure = build_structure(sequence, seed=spec.seed, plddt=plddt, coords=coords,
                                n_cols=spec.n_cols, pitch=spec.pitch)

    # -- choose the cohort residues -----------------------------------------
    plp = _select_plp(spec, coords, layout)
    if spec.blb_helix_ids:
        blb = [r for r in _from_helices(layout, spec.blb_helix_ids, spec.n_blb)
               if r not in set(plp)]
    elif spec.blb_min_distance_A is not None:
        blb = _outside_shell(coords, spec.cluster_center, spec.n_blb,
                             spec.blb_min_distance_A, spec.seed + 21, exclude=plp)
    else:
        blb = _spread(spec.length, spec.n_blb, spec.seed + 21, exclude=plp)
    conflict = _spread(spec.length, spec.n_conflict, spec.seed + 31,
                       exclude=set(plp) | set(blb)) if spec.n_conflict else []
    spare = _spread(spec.length, 6, spec.seed + 41,
                    exclude=set(plp) | set(blb) | set(conflict))

    # -- records -------------------------------------------------------------
    counter = [0]
    records: list[dict] = []
    records += _class_records(counter, sequence, plp, PLP_LABELS)
    records += _class_records(counter, sequence, blb, BLB_LABELS, offset=1)
    for k, position in enumerate(conflict):
        ref = sequence[position - 1]
        records.append(_record(counter, position=position, ref_aa=ref,
                               significance=PLP_LABELS[k % len(PLP_LABELS)],
                               review_status=REVIEW_STATUSES[k % len(REVIEW_STATUSES)]))
        records.append(_record(counter, position=position, ref_aa=ref,
                               significance=BLB_LABELS[k % len(BLB_LABELS)],
                               review_status=REVIEW_STATUSES[(k + 2)
                                                             % len(REVIEW_STATUSES)],
                               salt=2))

    # a P/LP record on a non-MANE transcript whose position and reference AA both
    # agree with the canonical sequence: remapped and INCLUDED (METHOD_SPEC II.13)
    anchor = plp[0]
    records.append(_record(counter, position=anchor, ref_aa=sequence[anchor - 1],
                           significance="Pathogenic", review_status=REVIEW_STATUSES[2],
                           transcript=ALT_TRANSCRIPT, salt=4))
    # a P/LP record whose review status is an unrecognised wording: star level NA
    records.append(_record(counter, position=anchor, ref_aa=sequence[anchor - 1],
                           significance="Likely pathogenic",
                           review_status="wording introduced after this release",
                           salt=5))

    noise = _noise_records(counter, sequence, plp_anchor=anchor, spare=spare)
    records += [row for row, _ in noise]

    # -- expected truth ------------------------------------------------------
    reason_counts: dict[str, int] = {}
    for _, reason in noise:
        reason_counts[reason] = reason_counts.get(reason, 0) + 1
    if conflict:
        reason_counts["residue_class_conflict"] = 2 * len(conflict)

    low_confidence = [i + 1 for i, value in enumerate(plddt) if value < 70.0]
    plp_segments = sequence_segments(plp)
    expected = {
        "case": name,
        # -- what the case is FOR (assert intent, not incidental numbers) -----
        "purpose": spec.purpose,
        "selection_strategy": spec.selection,
        "exercises": list(spec.exercises),
        "reaches_stage_cd": spec.reaches_stage_cd,
        "expected_negative_branch": spec.expected_negative_branch,
        # -- Stage A ground truth ---------------------------------------------
        "sequence_length": spec.length,
        "n_records": len(records),
        "M": spec.length,
        "N": len(plp) + len(blb),
        "N_P": len(plp),
        "N_B": len(blb),
        "n_conflict": len(conflict),
        "plp_residues": list(plp),
        "blb_residues": list(blb),
        "conflict_residues": list(conflict),
        "spare_residues": list(spare),
        "n_residues_stratum2": len(plp) + len(blb) + len(conflict),
        "exclusion_reason_counts": dict(sorted(reason_counts.items())),
        "n_excluded_records": sum(reason_counts.values()),
        # -- geometry the downstream stages depend on -------------------------
        "has_spatial_cluster": spec.selection in (SELECT_CLUSTER, SELECT_HELICES),
        "cluster_center_residue": spec.cluster_center,
        "n_plp_sequence_segments": len(plp_segments),
        "plp_sequence_segments": plp_segments,
        "max_mst_edge_plp_A": round(max_mst_edge(coords, plp), 6),
        "min_mst_edge_for_footprint_A": MIN_MST_EDGE_FOR_FOOTPRINT_A,
        "D_max_A": round(float(cross_distances(coords, coords).max()), 6),
        "n_helices": layout.n_helices,
        "bundle_cols": spec.n_cols,
        "bundle_pitch_A": spec.pitch,
        # -- pLDDT (recorded, never a filter) ---------------------------------
        "n_low_confidence_residues": len(low_confidence),
        "low_confidence_residues": low_confidence,
        "has_low_confidence_tail": spec.low_tail_from is not None,
        "low_confidence_tail_starts_at": (None if spec.low_tail_from is None
                                          else spec.low_tail_from + 1),
        "plddt_used_as_filter": False,
        "review_stars_used_for_inclusion": False,
        "ca_spacing_A": CA_SPACING_A,
    }

    return SyntheticCase(records=records, structure=structure, expected=expected,
                         name=name, sequence=sequence)


def make_synthetic_source(name: str) -> SyntheticStageASource:
    """Shorthand for ``make_synthetic_case(name).source``."""
    return make_synthetic_case(name).source


def all_synthetic_cases() -> dict[str, SyntheticCase]:
    """Every named case, built once — convenient for parametrized tests."""
    return {name: make_synthetic_case(name) for name in CASE_NAMES}


def consecutive_ca_distances(model: StructureModel) -> list[float]:
    """CA-CA distances along the chain; used to assert the 3.8 A backbone."""
    ordered = sorted(model.residues, key=lambda r: r.residue_index)
    points = [r.ca for r in ordered]
    return [math.dist(a.xyz, b.xyz) for a, b in zip(points, points[1:])]
