"""Targeted seed-determinism and seed-sensitivity tests — Phase 6, Root Cause 1.

Two distinct properties are asserted here and they must not be confused. Doing so
is a live hazard, not a hypothetical one: an earlier version of this module
required one identical ``r_hot`` from every seed, which is a determinism demand
made of a stochastic quantity, and the clustered fixture duly "failed" it while
every actual determinism test passed.

**Determinism** — same input, same effective config, same run_id must produce
byte-identical decision-bearing output: in the same process, in a fresh process,
under a different PYTHONHASHSEED, and under a different BLAS thread count. This is
a correctness property of the code, it is absolute, and all four forms hold.

**Seed sensitivity** — a different run_id is a different *legitimate* experiment.
Seeds derive from run_id (METHOD_SPEC II.12), so the permutation stream changes,
and ``perm_evidence_zg`` changes with it. Where two adjacent candidate radii are
near-tied in the II.6 distance-to-ideal ranking, that is enough to flip which one
is selected. The frozen rule is behaving correctly; what must be checked is that
no seed reaches a QUALITATIVELY different conclusion — a different hotspot, a
vanished hotspot, a non-adjacent radius, a fallback domain.

Nothing here changes any frozen rule: no threshold, null model, FDR procedure or
selection rule is touched, and the one tolerance used is read from the frozen
candidate grid spacing rather than invented. These tests only pin behaviour that
already exists.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

pytestmark = pytest.mark.integration

#: Decision-bearing Stage B artifacts. If any byte of these differs between two
#: runs at the same seed, the pipeline is nondeterministic.
DECISION_BEARING = [
    "FULL_RESULTS/04_GLOBAL_CLUSTERING/candidate_radius_domain.json",
    "FULL_RESULTS/05_HOTSPOT_RADIUS/r_hot_scan.tsv",
    "FULL_RESULTS/05_HOTSPOT_RADIUS/radius_decision.json",
    "FULL_RESULTS/05_HOTSPOT_RADIUS/pareto_front.json",
    "FULL_RESULTS/06_FINAL_HOTSPOTS/all_residue_center_tests.tsv",
    "FULL_RESULTS/06_FINAL_HOTSPOTS/significant_hotspot_centers.tsv",
    "FULL_RESULTS/06_FINAL_HOTSPOTS/hotspot_regions.tsv",
    "FULL_RESULTS/06_FINAL_HOTSPOTS/bh_fdr_table.tsv",
    "FULL_RESULTS/06_FINAL_HOTSPOTS/permutation_resolution_diagnostic.json",
]

FIXED_RUN_ID = "20250101T000000Z_det00001_det00002"

#: Distinct run_ids for the seed-robustness assertion. Each derives a completely
#: different permutation stream via II.12.
ROBUSTNESS_RUN_IDS = [
    "20250101T000000Z_seed0001_aaaaaaaa",
    "20250101T000000Z_seed0002_bbbbbbbb",
    "20250101T000000Z_seed0003_cccccccc",
    "20250101T000000Z_seed0004_dddddddd",
    "20250101T000000Z_seed0005_eeeeeeee",
    "20250101T000000Z_seed0006_ffffffff",
    "20250101T000000Z_seed0007_11111111",
    "20250101T000000Z_seed0008_22222222",
]


# --------------------------------------------------------------------------- #
# helpers                                                                     #
# --------------------------------------------------------------------------- #

@dataclass
class StageBRun:
    """What this module needs from a run: where it wrote, and how it ended."""

    run_root: Path
    status: str


def _run_stage_b(run_id: str, case: str = "clustered", root: Path | None = None):
    """Run **A -> B only** on a synthetic fixture and return a :class:`StageBRun`.

    Every artifact this module examines is a Stage B output — ``DECISION_BEARING``
    is entirely 04/05/06, and ``_summary`` reads ``handoff_02.json`` — so Stages
    C, D and E are executed for nothing. This helper used to call ``run_pipeline``,
    which runs the whole chain: each of the module's ~15 runs paid for the footprint
    optimization and the Stage D perturbation design, then ended ``BLOCKED`` at
    Stage E anyway because no annotation source is supplied here. Running only what
    is asserted on takes the module from ~61 minutes to a few.

    This is a TEST-EXECUTION change and nothing else. Stage A and Stage B are
    invoked through their own public entry points, in the order and with the
    arguments ``run_pipeline`` uses, under the same ``RunContext`` and the same
    frozen-methodology assertion. No stage is stubbed, no assertion is relaxed and
    production execution is untouched: ``run_pipeline`` still runs every stage.
    """
    sys.path.insert(0, str(REPO_ROOT))
    from hotspot3d.data.stage import run_stage_a
    from hotspot3d.hotspot.stage import run_stage_b as _stage_b
    from hotspot3d.reporting.manifest import ALL_STAGES
    from hotspot3d.utils.config import assert_frozen_methodology
    from hotspot3d.utils.errors import BlockedError, LeakageError, NegativeResult
    from hotspot3d.utils.runctx import RunContext
    from tests.fixtures.synthetic_clinvar import make_synthetic_case
    from tests.integration._config import load_integration_config

    root = root or Path(tempfile.mkdtemp())
    fixture = make_synthetic_case(case)
    cfg = load_integration_config(root / "config")
    assert_frozen_methodology(cfg)             # as run_pipeline does, before any stage

    ctx = RunContext.create(gene="SYNTH", config=cfg, results_root=root / "results",
                            synthetic=True, run_id=run_id)
    for stage in ALL_STAGES:                   # P3 — every stage directory exists
        ctx.stage_dir(stage)

    try:
        handoff_01 = run_stage_a(ctx, gene="SYNTH", source=fixture.source)
        _stage_b(ctx, upstream=handoff_01)
    except NegativeResult as nr:
        return StageBRun(run_root=ctx.run_root, status=f"NEGATIVE:{nr.condition}")
    except (BlockedError, LeakageError) as exc:
        return StageBRun(run_root=ctx.run_root, status=f"{type(exc).__name__}")
    return StageBRun(run_root=ctx.run_root, status="COMPLETED")


def _digests(run_root: Path) -> dict[str, str]:
    out = {}
    for rel in DECISION_BEARING:
        path = run_root / rel
        out[rel] = (hashlib.sha256(path.read_bytes()).hexdigest()
                    if path.is_file() else "ABSENT")
    return out


def _summary(result) -> dict:
    """The decision-bearing scalars, read from the handoff (never recomputed)."""
    path = result.run_root / "FULL_RESULTS/06_FINAL_HOTSPOTS/handoff_02.json"
    if not path.is_file():
        return {"present": False, "run_status": result.status}
    h2 = json.loads(path.read_text())
    centers_file = (result.run_root
                    / "FULL_RESULTS/06_FINAL_HOTSPOTS/significant_hotspot_centers.tsv")
    centers = []
    if centers_file.is_file():
        rows = centers_file.read_text().splitlines()
        if rows:
            header = rows[0].split("\t")
            idx = header.index("center_residue_index")
            centers = sorted(int(r.split("\t")[idx]) for r in rows[1:] if r)
    return {
        "present": True,
        "run_status": result.status,
        "r_hot": h2.get("hotspot_radius"),
        "n_significant_centers": h2.get("n_significant_centers"),
        "centers": centers,
        "bh_boundary_p": h2.get("bh_boundary_p"),
        "permutation_resolution_limited": h2.get("permutation_resolution_limited"),
        "b_recommended": h2.get("b_recommended"),
        "search_domain_source": h2.get("search_domain_source"),
        "fallback_radius_domain": h2.get("fallback_radius_domain"),
        "n_hotspot_regions": h2.get("n_hotspot_regions"),
    }


# --------------------------------------------------------------------------- #
# 1. determinism at a FIXED seed                                              #
# --------------------------------------------------------------------------- #

def test_same_seed_same_process_is_byte_identical():
    a = _run_stage_b(FIXED_RUN_ID)
    b = _run_stage_b(FIXED_RUN_ID)
    da, db = _digests(a.run_root), _digests(b.run_root)
    differing = [k for k in da if da[k] != db[k]]
    assert not differing, f"nondeterministic within one process: {differing}"
    assert _summary(a) == _summary(b)


def test_same_seed_fresh_process_is_byte_identical():
    """A fresh interpreter re-derives every seed from (MASTER_SEED, run_id, context).

    Catches any dependence on module-level RNG state carried across calls.
    """
    script = (
        "import sys, json; sys.path.insert(0, %r)\n"
        "sys.argv = ['x']\n"
        "from tests.integration.test_seed_determinism import "
        "_run_stage_b, _digests, _summary, FIXED_RUN_ID\n"
        "r = _run_stage_b(FIXED_RUN_ID)\n"
        "print('@@' + json.dumps({'d': _digests(r.run_root), 's': _summary(r)}))\n"
        % str(REPO_ROOT)
    )
    outs = []
    for hashseed in ("0", "424242"):
        env = {**os.environ, "PYTHONHASHSEED": hashseed,
               "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
        proc = subprocess.run([sys.executable, "-c", script], capture_output=True,
                              text=True, cwd=REPO_ROOT, env=env, timeout=3600)
        assert proc.returncode == 0, proc.stderr[-3000:]
        line = [l for l in proc.stdout.splitlines() if l.startswith("@@")][-1]
        outs.append(json.loads(line[2:]))
    assert outs[0]["d"] == outs[1]["d"], "PYTHONHASHSEED changed the output"
    assert outs[0]["s"] == outs[1]["s"]


def test_same_seed_different_thread_count_is_byte_identical():
    """BLAS thread count must not alter a single bit.

    The permutation null's hot path is a product of 0/1 matrices whose partial
    sums are exact integers far below 2**53, so no reordering can change it.
    This test is what makes that claim checkable rather than asserted.
    """
    script = (
        "import sys, json; sys.path.insert(0, %r)\n"
        "from tests.integration.test_seed_determinism import "
        "_run_stage_b, _digests, _summary, FIXED_RUN_ID\n"
        "r = _run_stage_b(FIXED_RUN_ID)\n"
        "print('@@' + json.dumps({'d': _digests(r.run_root), 's': _summary(r)}))\n"
        % str(REPO_ROOT)
    )
    outs = []
    for threads in ("1", "4"):
        env = {**os.environ, "OMP_NUM_THREADS": threads,
               "OPENBLAS_NUM_THREADS": threads, "MKL_NUM_THREADS": threads,
               "PYTHONHASHSEED": "0"}
        proc = subprocess.run([sys.executable, "-c", script], capture_output=True,
                              text=True, cwd=REPO_ROOT, env=env, timeout=3600)
        assert proc.returncode == 0, proc.stderr[-3000:]
        line = [l for l in proc.stdout.splitlines() if l.startswith("@@")][-1]
        outs.append(json.loads(line[2:]))
    assert outs[0]["d"] == outs[1]["d"], "thread count changed the output"
    assert outs[0]["s"] == outs[1]["s"]


# --------------------------------------------------------------------------- #
# 2. seed ROBUSTNESS of the scientific conclusion                             #
# --------------------------------------------------------------------------- #

@pytest.mark.slow
def test_clustered_conclusion_is_seed_stable_within_the_radius_grid():
    """Changing the scientific RNG seed may move r_hot by one grid step, no more.

    **This is a seed-SENSITIVITY assertion, not a determinism one.** The two are
    different properties and this module tests both; conflating them was the fault
    in the assertion this replaces.

    Determinism is pinned by the three tests above: at a FIXED run_id the
    decision-bearing bytes are identical in-process, in a fresh interpreter, under
    a different PYTHONHASHSEED and under a different BLAS thread count. All four
    hold. Nothing here weakens them.

    Seed sensitivity is a different question. Seeds derive from ``run_id``
    (METHOD_SPEC II.12), so a different run_id draws a different permutation
    stream, and ``perm_evidence_zg`` — one of the four equally-weighted objectives
    in the II.6 selection — moves with it. When two adjacent candidate radii are
    near-tied in distance-to-ideal, that movement can flip which one wins. That is
    the frozen methodology behaving correctly on stochastic evidence, not a defect,
    and demanding one identical r_hot from every seed would be demanding that a
    permutation test stop being a permutation test.

    What must hold instead is that no seed reaches a QUALITATIVELY different
    conclusion. This asserts exactly that:

      1. every seed reaches Stage B and finds a hotspot;
      2. the selected radii span at most ONE step of the candidate grid — the
         tolerance is read from ``radius_domain.step_hot_A``, the frozen grid
         spacing, so no new scientific threshold is introduced here;
      3. the core (intersection across every seed) is non-empty, and every seed
         agrees on the qualitative descriptors — region count, domain source,
         fallback flag, permutation-resolution flag — so no seed reaches a
         differently-shaped conclusion, only a slightly different edge on the
         same one.

    A fixture whose |S| flipped between 0 and 3, which is what the original
    assertion was written against, still fails 1 and 3.

    **Why this no longer demands identical center sets.** An earlier version of
    this test demanded byte-identical center sets from seeds that agreed on
    r_hot. Once §5.2's primary test became the structure-aware positional null
    (Workflow v2), that stopped holding: eight seeds at the frozen B and the SAME
    selected r_hot=5.5 A produced center counts ranging 38-43, not one shared
    count. Inspecting the actual sets showed why this is not a regression: the
    INTERSECTION across all eight is a 38-residue core, and every seed's "extra"
    residues sit immediately adjacent to that core's own boundary — exactly what
    a permutation test is expected to do at a decision edge. Demanding byte
    identity would be demanding the test stop being stochastic; the core-non-empty
    plus qualitative-descriptor-agreement checks below are what actually
    distinguishes that from a genuine disagreement about where the hotspot is.
    """
    from hotspot3d.utils.config import load_config
    from tests.integration._config import PIPELINE_CONFIG

    step = float(load_config(PIPELINE_CONFIG).get("radius_domain.step_hot_A"))
    summaries = {rid: _summary(_run_stage_b(rid)) for rid in ROBUSTNESS_RUN_IDS}

    detail = "\n".join(
        f"  {rid}: status={s['run_status']} r_hot={s['r_hot']} "
        f"n_sig={s['n_significant_centers']} centers={s.get('centers')}"
        for rid, s in summaries.items())

    # 1. a conclusion at all, and a positive one — the fixture is the happy path.
    assert all(s["present"] for s in summaries.values()), (
        f"the clustered fixture failed to reach Stage B under some seeds:\n{detail}")
    assert all(s["n_significant_centers"] > 0 for s in summaries.values()), (
        f"some seed found no hotspot in the clustered fixture:\n{detail}")

    # 2. r_hot varies by at most one step of the frozen candidate grid.
    radii = sorted({float(s["r_hot"]) for s in summaries.values()})
    spread = radii[-1] - radii[0]
    assert spread <= step + 1e-9, (
        f"r_hot spans {spread:g} A across seeds, more than one {step:g} A step of "
        f"the candidate grid — the selection is not merely breaking a near-tie "
        f"between adjacent candidates:\n{detail}")

    # 3. a non-empty core every seed agrees on, and no qualitative contradiction —
    #    including, via the loop below, n_hotspot_regions itself.
    #
    #    An earlier version of this check also tried to verify "every seed stays
    #    within the same regions" by grouping each seed's centers into maximal
    #    RESIDUE-INDEX-contiguous runs. That is the wrong tool for this fixture on
    #    its face, not merely imprecise: hotspot regions are 3D spatial connected
    #    components (`d(c,c') <= 2*r_hot`), and this synthetic structure is a
    #    packed helix bundle deliberately built so ONE spatially-contiguous
    #    hotspot is assembled from several residue-index-distant sequence
    #    segments (see data/synthetic.py's module docstring — "a spatial
    #    neighbourhood therefore draws from SEVERAL SEQUENCE SEGMENTS, which is
    #    the entire reason this analysis is done in 3D rather than on the
    #    sequence"). Every seed here reports n_hotspot_regions=1 for a center set
    #    spanning three residue-index-separated runs — checking residue-index
    #    contiguity would have flagged that as 3 disjoint regions, which is a
    #    wrong reading of the data, not a stricter one. n_hotspot_regions IS
    #    Stage B's own already-validated spatial computation; re-deriving a
    #    competing, sequence-based notion of "region" here would test this
    #    module's own invented concept instead of the actual methodology.
    all_sets = [set(s["centers"]) for s in summaries.values()]
    core = set.intersection(*all_sets)
    assert core, f"the seeds share no common significant center:\n{detail}"
    for key in ("n_hotspot_regions", "search_domain_source",
                "fallback_radius_domain", "permutation_resolution_limited"):
        values = {s[key] for s in summaries.values()}
        assert len(values) == 1, (
            f"seeds disagree on {key} ({sorted(values, key=str)}), which is a "
            f"qualitative difference and not a near-tie:\n{detail}")
    assert summaries[ROBUSTNESS_RUN_IDS[0]]["fallback_radius_domain"] is False, (
        f"the candidate domain fell back to the envelope:\n{detail}")


@pytest.mark.slow
def test_radius_scan_is_not_degenerate():
    """At least three admissible radii, so selection is genuinely discriminating.

    With only two admissible candidates, min-max normalization over the
    admissible set collapses every objective to 0 or 1 and both candidates land
    at distance_to_ideal = 1.0 — selection then rests entirely on the tie chain,
    which is exactly how a seed-driven flip in r_hot arises. The frozen rule is
    correct; a fixture that starves it of candidates is not a fair test of it.
    """
    from hotspot3d.utils.io import read_tsv
    result = _run_stage_b(FIXED_RUN_ID)
    scan = result.run_root / "FULL_RESULTS/05_HOTSPOT_RADIUS/r_hot_scan.tsv"
    assert scan.is_file(), "no radius scan was produced"
    rows = read_tsv(scan)
    admissible = [r for r in rows
                  if str(r.get("qc_status")).upper() in ("PASS", "ADMISSIBLE")]
    distances = {r.get("distance_to_ideal") for r in admissible}
    assert len(admissible) >= 3, (
        f"only {len(admissible)} admissible radii of {len(rows)} scanned — "
        f"normalization is degenerate and the tie chain decides r_hot")
    assert len(distances) > 1, (
        f"every admissible candidate has the same distance-to-ideal {distances} — "
        f"the objectives are not discriminating on this fixture")
