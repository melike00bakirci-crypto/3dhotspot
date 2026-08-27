"""Deterministic synthetic inputs for Stage E.

Two things live here:

* :class:`MockAnnotationSource` — a complete :class:`AnnotationSource` returning fixed
  synthetic records with explicit **fake** version strings. It deliberately contains
  the awkward cases the rubric exists for: a variant with two conflicting STRONG
  sources, a WEAK-only variant, a variant with no evidence at all, and (by default) a
  mechanism distribution that fails the post hoc gate.

* :func:`build_synthetic_run` — a full synthetic upstream run root (Stages A–D
  artifacts plus handoffs 02/03/04) so ``run_stage_e`` can be executed end to end
  without any network access and without any other stage being implemented.

Nothing here is real biological data. Version strings begin with ``MOCK-`` so a
synthetic run can never be mistaken for a pinned live one.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from hotspot3d.annotation.inputs import UPSTREAM_PATHS
from hotspot3d.annotation.rubric import Evidence, Mechanism, MechanismRecord
from hotspot3d.annotation.sources import SourceVersion
from hotspot3d.orchestration.contracts import Handoff
from hotspot3d.utils.config import load_config
from hotspot3d.utils.hashing import manifest_for
from hotspot3d.utils.io import write_json, write_tsv
from hotspot3d.utils.runctx import RunContext

SEARCH_DATE = "2026-01-15"
DB = "PubMed"

#: Resolved from this file, so tests do not depend on the working directory.
DEFAULT_CONFIG = Path(__file__).resolve().parents[2] / "config" / "pipeline.yaml"

#: Fake, explicit versions. Every source Stage E touches must be pinned.
MOCK_VERSIONS: dict[str, SourceVersion] = {
    "DSSP": SourceVersion("DSSP", "MOCK-4.4.0", "https://example.invalid/dssp",
                          "NA", "2026-01-15T00:00:00Z"),
    "InterPro": SourceVersion("InterPro", "MOCK-98.0", "https://example.invalid/interpro",
                              "IPR000000", "2026-01-15T00:00:00Z"),
    "UniProt": SourceVersion("UniProt", "MOCK-2026_01", "https://example.invalid/uniprot",
                             "P00000", "2026-01-15T00:00:00Z"),
    "dbNSFP": SourceVersion("dbNSFP", "MOCK-4.9a", "https://example.invalid/dbnsfp",
                            "NA", "2026-01-15T00:00:00Z"),
    "PubMed": SourceVersion("PubMed", "MOCK-2026-01-15", "https://example.invalid/pubmed",
                            "NA", "2026-01-15T00:00:00Z"),
}

# --- protein / hotspot layout ----------------------------------------------

N_RESIDUES = 200
HOTSPOTS: dict[str, tuple[int, int]] = {"HS1": (100, 115), "HS2": (150, 160)}
CENTERS: dict[str, list[int]] = {"HS1": [104, 108], "HS2": [153, 157]}

#: The three interesting rubric cases, all inside HS1.
CONFLICTING_STRONG_RESIDUE = 101      # -> Mixed
WEAK_ONLY_RESIDUE = 105               # -> Unclear
NO_EVIDENCE_RESIDUE = 110             # -> Not_Experimentally_Characterized

#: Residues used to build a gate-passing distribution.
GOF_RESIDUES = [102, 103, 104, 106, 107, 108, 109, 111, 112, 113]
LOF_RESIDUES = list(range(150, 160))


def coords_for(residue: int) -> tuple[float, float, float]:
    """Deterministic synthetic Cα coordinates."""
    return (residue * 3.0, float(residue % 7), float(residue % 5))


def plddt_for(residue: int) -> float:
    return round(90.0 - (residue % 20), 1)


# --- mechanism records ------------------------------------------------------

def _record(residue: int, mech: Mechanism, evidence: Evidence, ref: str, **kw) -> MechanismRecord:
    return MechanismRecord(
        residue_index=residue,
        aa_change=kw.pop("aa_change", f"p.Xaa{residue}Tyr"),
        mechanism=mech,
        reference=ref,
        experimental_system=kw.pop("experimental_system", "HEK293 transient expression"),
        assay_type=kw.pop("assay_type", "whole-cell patch clamp"),
        principal_finding=kw.pop("principal_finding", f"{mech.value} effect reported"),
        evidence_strength=evidence,
        query=kw.pop("query", f"MOCK[gene] AND residue {residue}"),
        database=kw.pop("database", DB),
        search_date=kw.pop("search_date", SEARCH_DATE),
        hit_count=kw.pop("hit_count", 3),
        included=kw.pop("included", True),
        exclusion_reason=kw.pop("exclusion_reason", "NA"),
    )


def special_records() -> list[MechanismRecord]:
    """The three cases the rubric must get right, plus one excluded record."""
    return [
        # Two conflicting STRONG sources -> Mixed, both references retained.
        _record(CONFLICTING_STRONG_RESIDUE, Mechanism.GOF, Evidence.STRONG,
                "PMID:MOCK0001", aa_change="p.Arg101Gln",
                principal_finding="Increased current density vs wild type"),
        _record(CONFLICTING_STRONG_RESIDUE, Mechanism.LOF, Evidence.STRONG,
                "PMID:MOCK0002", aa_change="p.Arg101Gln",
                experimental_system="Xenopus oocytes",
                principal_finding="Reduced surface expression and current"),
        # WEAK only -> Unclear.
        _record(WEAK_ONLY_RESIDUE, Mechanism.LOF, Evidence.WEAK, "PMID:MOCK0003",
                aa_change="p.Gly105Ser", assay_type="in silico prediction",
                principal_finding="Computational prediction reported in source"),
        # No evidence -> Not_Experimentally_Characterized.
        _record(NO_EVIDENCE_RESIDUE, Mechanism.LOF, Evidence.NONE, "PMID:MOCK0004",
                aa_change="p.Val110Met", assay_type="NA",
                principal_finding="Variant listed; no functional assay performed"),
        # An excluded candidate — exclusions must leave a trace in the search log.
        _record(175, Mechanism.GOF, Evidence.MODERATE, "PMID:MOCK0005",
                aa_change="p.Ala175Thr", included=False,
                exclusion_reason="Assay performed on a different gene product"),
    ]


def make_mechanism_records(
    *, n_gof: int = 0, n_lof: int = 0, evidence: Evidence = Evidence.MODERATE,
    include_special: bool = True,
) -> list[MechanismRecord]:
    """Build a mechanism set with an exactly controlled category distribution.

    Used to drive the post hoc gate boundary tests (9 vs 10 variants, 1 vs 2
    categories) without touching the rubric.
    """
    records = special_records() if include_special else []
    for residue in GOF_RESIDUES[:n_gof]:
        records.append(_record(residue, Mechanism.GOF, evidence,
                               f"PMID:MOCKG{residue}"))
    for residue in LOF_RESIDUES[:n_lof]:
        records.append(_record(residue, Mechanism.LOF, evidence,
                               f"PMID:MOCKL{residue}"))
    return records


# --- the mock source --------------------------------------------------------

class MockAnnotationSource:
    """Deterministic, offline, version-pinned annotation source.

    ``gate_passing=False`` (the default) yields too few characterized variants for the
    post hoc mechanism-spatial analysis, which is the common real-world case and the
    one the gate exists to handle.
    """

    def __init__(
        self, *, gate_passing: bool = False, versioned: bool = True,
        records: list[MechanismRecord] | None = None,
        with_domains: bool = True, with_conservation: bool = True,
    ):
        self.gate_passing = gate_passing
        self.versioned = versioned
        self.with_domains = with_domains
        self.with_conservation = with_conservation
        self._records = records
        self.calls: list[str] = []

    # -- version pinning ----------------------------------------------------
    def versions(self) -> dict[str, SourceVersion]:
        if not self.versioned:
            # An unversioned source. Using it is BLOCKING.
            return {
                name: SourceVersion(name, "", sv.url, sv.accession, sv.retrieved_utc)
                for name, sv in MOCK_VERSIONS.items()
            }
        return dict(MOCK_VERSIONS)

    # -- annotation surface -------------------------------------------------
    def secondary_structure(self, structure: Any) -> dict[int, str]:
        self.calls.append("secondary_structure")
        codes = "HHHHEEEECCCC"
        return {i: codes[i % len(codes)] for i in range(1, N_RESIDUES + 1)}

    def domains(self, acc: str) -> list[dict]:
        self.calls.append("domains")
        if not self.with_domains:
            return []
        return [
            {"domain_id": "IPR000001", "name": "MOCK_kinase_domain", "start": 90,
             "end": 130, "source": "InterPro", "overlap_residues": "NA"},
            {"domain_id": "IPR000002", "name": "MOCK_regulatory_domain", "start": 300,
             "end": 340, "source": "InterPro", "overlap_residues": "NA"},
        ]

    def functional_sites(self, acc: str) -> list[dict]:
        self.calls.append("functional_sites")
        return [
            {"site_id": "SITE1", "site_type": "ligand_binding", "residue_start": 103,
             "residue_end": 107, "description": "MOCK ATP-binding pocket",
             "source": "UniProt"},
            {"site_id": "SITE2", "site_type": "protein_interaction",
             "residue_start": 152, "residue_end": 156,
             "description": "MOCK partner-binding interface", "source": "UniProt"},
            {"site_id": "SITE3", "site_type": "functional_region", "residue_start": 100,
             "residue_end": 116, "description": "MOCK activation loop",
             "source": "UniProt"},
        ]

    def conservation(self, acc: str, residues: Iterable[int]) -> dict[int, dict]:
        self.calls.append("conservation")
        if not self.with_conservation:
            return {}
        return {
            int(r): {
                "GERP++_RS": round(4.0 + (int(r) % 5) * 0.25, 3),
                "phyloP100way": round(2.0 + (int(r) % 7) * 0.1, 3),
            }
            for r in residues
        }

    def disease_associations(self, gene: str) -> list[dict]:
        self.calls.append("disease_associations")
        # HS2 deliberately has no association — a valid, informative negative.
        return [
            {"disease": "MOCK developmental disorder", "residue_start": 98,
             "residue_end": 120, "reference": "PMID:MOCK9001", "source": "PubMed"},
        ]

    def literature_mechanisms(self, gene: str) -> list[MechanismRecord]:
        self.calls.append("literature_mechanisms")
        if self._records is not None:
            return list(self._records)
        if self.gate_passing:
            return make_mechanism_records(n_gof=10, n_lof=10)
        return special_records()


# --- synthetic upstream run -------------------------------------------------

@dataclass
class SyntheticRun:
    ctx: RunContext
    handoff_02: Handoff
    handoff_03: Handoff
    handoff_04: Handoff
    paths: dict[str, Path]


def upstream_path(full_results: Path, name: str) -> Path:
    """Resolve a synthetic artifact through Stage E's declared read contract.

    The fixture never hard-codes an upstream filename: if the Lead re-points a path
    in ``UPSTREAM_PATHS``, the synthetic run follows automatically instead of
    silently writing a file Stage E no longer looks for.
    """
    stage, rel = UPSTREAM_PATHS[name]
    return full_results / stage / rel


def build_synthetic_run(
    tmp_path: Path, *, gene: str = "MOCKG", negative: bool = False,
    config_path: str | Path | None = None,
    flags: dict | None = None, run_id: str | None = None,
) -> SyntheticRun:
    """Materialize a complete synthetic Stage A–D run root plus handoffs 02/03/04.

    ``negative=True`` produces a valid upstream scientific negative: no hotspot
    regions and ``n_significant_centers = 0``.
    """
    cfg = load_config(config_path or DEFAULT_CONFIG)
    ctx = RunContext.create(gene, cfg, results_root=tmp_path / "results", synthetic=True,
                            run_id=run_id)
    fr = ctx.full_results
    paths: dict[str, Path] = {}

    hotspots = {} if negative else HOTSPOTS
    residues_by_hotspot = {
        hid: list(range(start, end + 1)) for hid, (start, end) in hotspots.items()
    }

    # --- 06_FINAL_HOTSPOTS ---
    region_rows = []
    for i, (hid, (start, end)) in enumerate(hotspots.items()):
        region_rows.append({
            "hotspot_id": hid, "residue_start": start, "residue_end": end,
            "n_plp": 12 + i, "n_blb": 2 + i,
            "fold_enrichment": round(3.25 + i * 0.5, 6),
            "p_emp": 0.0001 + i * 0.0002, "q_bh": 0.002 + i * 0.001,
            "Z": round(4.12 - i * 0.3, 6),
        })
    paths["hotspot_regions"] = write_tsv(
        upstream_path(fr, "hotspot_regions"), region_rows,
        ["hotspot_id", "residue_start", "residue_end", "n_plp", "n_blb",
         "fold_enrichment", "p_emp", "q_bh", "Z"],
    )

    paths["significant_hotspot_centers"] = write_tsv(
        upstream_path(fr, "significant_hotspot_centers"),
        [{"hotspot_id": hid, "center_residue_index": c, "p_emp": 0.0001, "q_bh": 0.002}
         for hid in hotspots for c in CENTERS[hid]],
        ["hotspot_id", "center_residue_index", "p_emp", "q_bh"],
    )
    # NOTE: the real Stage B output (hotspot-statistics) names this column
    # "hotspot_ids" (plural, ";"-separated) on BOTH these two files, because a
    # variant or a covered residue can legitimately fall inside more than one
    # overlapping hotspot's sphere -- unlike significant_hotspot_centers.tsv
    # above, where "hotspot_id" is genuinely singular. This fixture matches
    # that real schema so _group_by_hotspot's real behavior is what's tested.
    paths["hotspot_sphere_classified_variants"] = write_tsv(
        upstream_path(fr, "hotspot_sphere_classified_variants"),
        [{"hotspot_ids": hid, "variant_residue_index": r,
          "class": "P/LP" if r % 3 else "B/LB"}
         for hid, rs in residues_by_hotspot.items() for r in rs[::2]],
        ["hotspot_ids", "variant_residue_index", "class"],
    )
    paths["hotspot_covered_residues"] = write_tsv(
        upstream_path(fr, "hotspot_covered_residues"),
        [{"hotspot_ids": hid, "covered_residue_index": r}
         for hid, rs in residues_by_hotspot.items() for r in rs],
        ["hotspot_ids", "covered_residue_index"],
    )
    paths["all_residue_center_tests"] = write_tsv(
        upstream_path(fr, "all_residue_center_tests"),
        [{"residue_index": r, "p_emp": 0.5, "q_bh": 0.9}
         for r in range(1, N_RESIDUES + 1)],
        ["residue_index", "p_emp", "q_bh"],
    )

    # --- 03_STRUCTURE_QC ---
    paths["plddt_profile"] = write_tsv(
        upstream_path(fr, "plddt_profile"),
        [{"residue_index": r, "plddt": plddt_for(r)} for r in range(1, N_RESIDUES + 1)],
        ["residue_index", "plddt"],
    )
    cif = upstream_path(fr, "structure_with_plddt")
    cif.parent.mkdir(parents=True, exist_ok=True)
    cif.write_text("data_MOCK\n# synthetic placeholder structure\n", encoding="utf-8")
    paths["structure_with_plddt"] = cif

    # --- 02_CLINVAR ---
    cohort_rows = []
    for r in range(1, N_RESIDUES + 1):
        x, y, z = coords_for(r)
        cohort_rows.append({
            "residue_index": r, "class": "P/LP" if r % 3 else "B/LB",
            "x_ca": x, "y_ca": y, "z_ca": z, "plddt": plddt_for(r), "ca_usable": True,
        })
    paths["classified_cohort"] = write_tsv(
        upstream_path(fr, "classified_cohort"), cohort_rows,
        ["residue_index", "class", "x_ca", "y_ca", "z_ca", "plddt", "ca_usable"],
    )
    paths["variants_residue_level"] = write_tsv(
        upstream_path(fr, "variants_residue_level"),
        [{"residue_index": r, "class": "P/LP" if r % 3 else "B/LB", "n_records": 1 + r % 3}
         for r in range(1, N_RESIDUES + 1)],
        ["residue_index", "class", "n_records"],
    )
    paths["review_star_distribution"] = write_tsv(
        upstream_path(fr, "review_star_distribution"),
        [{"residue_index": r, "max_star": r % 4} for r in range(1, N_RESIDUES + 1)],
        ["residue_index", "max_star"],
    )

    # --- 08_FINAL_FOOTPRINT / 09_ROBUSTNESS ---
    paths["footprint_residues"] = write_tsv(
        upstream_path(fr, "footprint_residues"),
        [{"residue_index": r} for rs in residues_by_hotspot.values() for r in rs],
        ["residue_index"],
    )
    profile = {
        "preservation_rate": 0.82, "mean_jaccard": 0.79, "sd_jaccard": 0.08,
        "n_iterations": 120, "median_dice": 0.86,
        "geometric_footprint_robustness_only": True,
    }
    paths["robustness_profile"] = write_json(
        upstream_path(fr, "robustness_profile"), profile)
    # Distinct per-center values, so a test that collapsed centers would be visible.
    all_centers = [c for hid in hotspots for c in CENTERS[hid]]
    paths["center_sensitivity"] = write_tsv(
        upstream_path(fr, "center_sensitivity"),
        [{"center_residue_index": c, "mean_jaccard_drop": round(0.10 + i * 0.07, 3),
          "center_recurrence_rate": round(0.95 - i * 0.04, 3)}
         for i, c in enumerate(all_centers)],
        ["center_residue_index", "mean_jaccard_drop", "center_recurrence_rate"],
    )

    # --- 11_SENSITIVITY ---
    paths["sensitivity_overlap"] = write_tsv(
        upstream_path(fr, "sensitivity_overlap"),
        [{"hotspot_id": hid, "overlap_ge1star": 0.91, "overlap_ge2star": 0.74}
         for hid in hotspots],
        ["hotspot_id", "overlap_ge1star", "overlap_ge2star"],
    )

    # --- handoffs ---
    base_flags = {
        "global_clustering_flag": "significant",
        "permutation_resolution_limited": False,
        "fallback_radius_domain": False,
        "domain_boundary_warning": False,
    }
    base_flags.update(flags or {})

    upstream_files = list(paths.values())
    manifest = manifest_for(upstream_files, root=ctx.run_root)

    h02 = Handoff(
        name="handoff_02", run_id=ctx.run_id, config_sha256=cfg.sha256,
        qc_status="PASS", manifest=manifest,
        payload={
            "hotspot_radius": 8.5, "search_domain_source": "pcf",
            "q": 0.05, "fdr_method": "BH", "B": 10000,
            "n_significant_centers": 0 if negative else sum(
                len(CENTERS[h]) for h in hotspots),
            "n_centers_without_variant": 1, "bh_boundary_p": 0.004,
            "b_recommended": 10000, "near_tie_flag": False,
            "loo_sparse_proportion": 0.05, "loo_zero_neighbour_proportion": 0.01,
            "sensitivity_summary": {"ge1star_overlap": 0.91, "ge2star_overlap": 0.74},
            "code_version": "MOCK", "uniprot_acc": "P00000",
            **base_flags,
        },
        negative_result=(
            {"condition": "NO_SIGNIFICANT_HOTSPOTS",
             "detail": "No center reached BH-FDR significance."} if negative else None
        ),
    )
    h03 = Handoff(
        name="handoff_03", run_id=ctx.run_id, config_sha256=cfg.sha256,
        qc_status="PASS", manifest=manifest,
        payload={
            "hotspot_radius": 8.5, "footprint_radius": 6.0, "n_components": 2,
            "volume_A3": 12000.0, "surface_area_A2": 3400.0,
            "n_footprint_residues": sum(len(v) for v in residues_by_hotspot.values()),
            "coverage": 0.14, "voxel_h": 0.5, "near_tie_flag": False,
            "domain_boundary_warning": base_flags["domain_boundary_warning"],
            "code_version": "MOCK",
        },
    )
    h04 = Handoff(
        name="handoff_04", run_id=ctx.run_id, config_sha256=cfg.sha256,
        qc_status="PASS", manifest=manifest,
        payload={
            "footprint_code_version": "MOCK", "robustness_profile": profile,
            "n_S": sum(len(CENTERS[h]) for h in hotspots), "K_MAX": 2,
            "total_subset_space": 6, "total_evaluated": 6,
            "per_level": [{"k": 1, "mode": "exhaustive", "n_evaluated": 4}],
            "robustness_r_fp_invariant": {"boundary": "none"},
            "clinical_dataset_unmodified": True,
            "geometric_footprint_robustness_only": True,
        },
    )
    return SyntheticRun(ctx=ctx, handoff_02=h02, handoff_03=h03, handoff_04=h04,
                        paths=paths)
