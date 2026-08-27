"""Stage D robustness iteration parallelism: bit-identical results at any
worker count (II.12).

Scope, deliberately bounded per the authorized task: a SMALL fixed set of 20
real iterations, run serially and with 2 and 4 workers. Never spawns more than
4 workers, never runs a scaling benchmark, and never runs the full
N_CAP=10,000 robustness set -- the SubsetDesign object built here is the real,
unmodified production design (N_CAP/K_MAX/sampling policy all untouched), only
its first 20 already-computed subsets are actually EXECUTED.
"""
from __future__ import annotations

import os

import pytest

from synthetic_footprint import synthetic_run, two_cluster_centers

from hotspot3d.footprint.stage import run_phase_c
from hotspot3d.robustness.design import build_design
from hotspot3d.robustness.parallel import run_iterations
from hotspot3d.robustness.params import RobustnessParams
from hotspot3d.robustness.stage import _aggregate_all
from hotspot3d.utils.seeds import ctx_iteration

pytestmark = pytest.mark.unit

N_TEST_ITERATIONS = 20


@pytest.fixture(scope="module")
def parallel_fixture(tmp_path_factory):
    """One real Phase C run (six centers, two merging clusters) plus the real,
    unmodified production subset design -- built once and shared by every test
    below so Stage C only runs once."""
    tmp = tmp_path_factory.mktemp("parallel")
    ctx, upstream = synthetic_run(tmp, center_ids=two_cluster_centers())
    handoff_03, phase_c = run_phase_c(ctx, upstream=upstream)
    up = phase_c.inputs
    params = RobustnessParams.from_config(ctx.config)

    design = build_design(len(up.centers), params.n_cap, ctx.seeds, up.centers.ids)
    assert design.total_evaluated >= N_TEST_ITERATIONS, (
        "fixture must produce at least N_TEST_ITERATIONS real subsets to slice "
        "from without touching N_CAP")
    assert design.N_CAP == params.n_cap == 10000, "N_CAP must be the real, unaltered value"

    mode_by_k = {level.k: level.mode for level in design.levels}
    all_tasks = [(index, removed, mode_by_k[len(removed)], ctx_iteration(index))
                for index, removed in enumerate(design.subsets)]
    tasks = all_tasks[:N_TEST_ITERATIONS]

    original_state = phase_c.solution.selected
    original_residues = frozenset(phase_c.solution.covered_residue_ids(up.universe))
    common = dict(
        centers=up.centers, universe=up.universe, params=phase_c.params,
        original_state=original_state, original_residues=original_residues,
        r_fp_original=phase_c.r_fp, plp=up.cohort_plp, blb=up.cohort_blb,
        r_hot_value=up.r_hot.value_for_coverage_query,
        undefined_mcc=params.undefined_mcc)
    return ctx, tasks, common


@pytest.fixture(scope="module")
def three_runs(parallel_fixture):
    """Run the SAME 20 real tasks serially, with 2 workers, and with 4 workers."""
    _ctx, tasks, common = parallel_fixture
    serial = run_iterations(tasks, n_workers=1, worker_nice=10, **common)
    two_workers = run_iterations(tasks, n_workers=2, worker_nice=10, **common)
    four_workers = run_iterations(tasks, n_workers=4, worker_nice=10, **common)
    return serial, two_workers, four_workers


def test_worker_count_never_exceeds_four_in_this_module():
    """Structural guard on the test file itself, per the authorized task scope."""
    import re
    import sys
    from pathlib import Path

    source = Path(sys.modules[__name__].__file__).read_text(encoding="utf-8")
    for match in re.finditer(r"n_workers\s*=\s*(\d+)", source):
        assert int(match.group(1)) <= 4, (
            f"this test module must not spawn more than 4 workers "
            f"(found n_workers={match.group(1)})")


def test_serial_path_uses_no_process_pool(monkeypatch, parallel_fixture):
    """n_workers <= 1 must never construct a ProcessPoolExecutor at all."""
    import hotspot3d.robustness.parallel as parallel_mod

    _ctx, tasks, common = parallel_fixture

    def exploding(*_args, **_kwargs):
        raise AssertionError("ProcessPoolExecutor must not be constructed for n_workers=1")

    monkeypatch.setattr(parallel_mod, "ProcessPoolExecutor", exploding)
    result = run_iterations(tasks[:3], n_workers=1, worker_nice=10, **common)
    assert len(result) == 3


def test_effective_workers_never_exceeds_task_count(monkeypatch, parallel_fixture):
    """A run with fewer iterations than configured workers must not over-spawn."""
    import hotspot3d.robustness.parallel as parallel_mod

    _ctx, tasks, common = parallel_fixture
    captured = {}
    real_executor = parallel_mod.ProcessPoolExecutor

    class _Spy(real_executor):
        def __init__(self, *args, **kwargs):
            captured["max_workers"] = kwargs.get("max_workers")
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(parallel_mod, "ProcessPoolExecutor", _Spy)
    result = run_iterations(tasks[:3], n_workers=4, worker_nice=10, **common)
    assert len(result) == 3
    assert captured["max_workers"] == 3, (
        "must cap DOWN to the number of tasks, never spawn idle workers — and "
        "this must never be confused with auto-scaling UP to the node's CPU count")


def test_serial_and_2_and_4_workers_produce_identical_per_iteration_rows(three_runs):
    serial, two_workers, four_workers = three_runs
    assert len(serial) == len(two_workers) == len(four_workers) == N_TEST_ITERATIONS

    for i, (s, w2, w4) in enumerate(zip(serial, two_workers, four_workers)):
        assert s.index == w2.index == w4.index == i
        assert s.k == w2.k == w4.k
        assert s.removed == w2.removed == w4.removed
        assert s.status == w2.status == w4.status
        assert s.outcome_type == w2.outcome_type == w4.outcome_type
        assert s.robustness_r_fp == w2.robustness_r_fp == w4.robustness_r_fp
        assert s.residues == w2.residues == w4.residues
        # the full per-iteration metric row: every geometry + classification
        # field, bit-for-bit, not merely a spot check
        assert s.row == w2.row == w4.row, f"iteration {i}: row differs across worker counts"


def test_serial_and_2_and_4_workers_produce_identical_aggregate_summaries(
        parallel_fixture, three_runs):
    ctx, tasks, _common = parallel_fixture
    params = RobustnessParams.from_config(ctx.config)
    serial, two_workers, four_workers = three_runs

    agg_serial = _aggregate_all(ctx, serial, len(tasks), params)
    agg_two = _aggregate_all(ctx, two_workers, len(tasks), params)
    agg_four = _aggregate_all(ctx, four_workers, len(tasks), params)

    assert agg_serial == agg_two == agg_four, (
        "aggregate summaries (n, mean, SD, CI, CV, median, IQR per metric) must "
        "be identical regardless of worker count")


def test_at_least_one_real_iteration_actually_completed(three_runs):
    """A guard against the comparison being vacuously true over all-failed rows."""
    serial, _two, _four = three_runs
    assert any(r.status == "OK" for r in serial), (
        "the fixture must produce at least one COMPLETED iteration, or the "
        "identical-results assertions above would be vacuous")


def test_workers_pin_blas_threads_and_run_niced(monkeypatch):
    """Ground-truth check, not merely env vars: threadpoolctl.threadpool_info()
    on the ACTUAL loaded BLAS library inside a real worker process, plus the
    worker's own niceness -- and that the PARENT process's niceness is
    untouched (only workers are niced, never the caller)."""
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor

    from hotspot3d.robustness.parallel import _init_worker

    # simulate a caller that already had many BLAS threads configured, to
    # prove the worker's pinning genuinely overrides it rather than merely
    # happening to already be 1
    monkeypatch.setenv("OMP_NUM_THREADS", "8")
    monkeypatch.setenv("MKL_NUM_THREADS", "8")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "8")
    parent_nice_before = os.nice(0)

    mp_ctx = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=2, mp_context=mp_ctx,
                             initializer=_init_worker,
                             initargs=({}, 10)) as pool:
        reports = list(pool.map(_thread_pinning_probe, range(2)))

    assert os.nice(0) == parent_nice_before, "the PARENT process must never be niced"
    assert len(reports) == 2
    for report in reports:
        assert report["OMP_NUM_THREADS"] == "1"
        assert report["MKL_NUM_THREADS"] == "1"
        assert report["OPENBLAS_NUM_THREADS"] == "1"
        # the ground truth: the ACTUAL loaded BLAS runtime, not just env vars
        assert all(n == 1 for n in report["threadpool_info"].values()), report
        assert report["nice"] > parent_nice_before, (
            "the worker must be MORE niced than the parent, i.e. lower priority")


def _thread_pinning_probe(_):
    """Runs INSIDE the worker process (must be a top-level, picklable function)."""
    import os

    import threadpoolctl

    info = threadpoolctl.threadpool_info()
    limits = {d.get("user_api", d.get("internal_api")): d.get("num_threads")
             for d in info}
    return {
        "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
        "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS"),
        "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
        "threadpool_info": limits,
        "nice": os.nice(0),
    }
