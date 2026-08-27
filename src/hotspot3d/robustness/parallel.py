"""Process-pool parallelism for Stage D robustness iterations (II.11/II.12).

Iterations are independent by construction: each is a pure function of
``(index, removed-subset, mode, seed_context)`` plus the SAME read-only inputs
(``centers``, ``universe``, ``params``, ``original_state``, ``original_residues``,
``r_fp_original``, ``plp``, ``blb``, ``r_hot_value``, ``undefined_mcc`` — see
:func:`hotspot3d.robustness.perturb.run_iteration`). No iteration consumes a
random-number generator: the one seeded step, subset assignment, already ran
single-threaded in :func:`hotspot3d.robustness.design.build_design` before this
module is ever invoked (design.py's own module docstring: "seeded from its
context string, never from execution order"). No iteration mutates shared
state. Running iterations out of order or on different workers therefore
cannot change any individual iteration's own result.

**BIT-IDENTICAL guarantee.** This module distributes PRE-COMPUTED
``(index, removed, mode, seed_context)`` tuples to workers and collects
results via ``ProcessPoolExecutor.map()``, which returns results in INPUT
order regardless of completion order — worker count and execution order
therefore never affect which result lands at which position in the output
list, which is what makes ``perturbation_results.tsv`` and every aggregate
built from it identical across worker counts. Verified directly (not merely
argued) by ``tests/robustness/test_robustness_parallel.py``: 20 real iterations
run serially and with 2 and 4 workers, asserting every per-iteration row AND
every aggregate summary are identical.

**Process isolation, not threads.** This sidesteps the class of global-state
leakage this codebase has hit before (a module-level rcParams mutation in one
call bleeding into another call in the SAME process — see
``tests/integration/test_figure_determinism.py``). Each worker is a fresh OS
process with its own copy of every module's global state. The one thing that
IS sent to workers — ``centers``/``universe``/``params``/``original_state``/etc.
— is sent ONCE per worker at startup, read-only, and never mutated by
``run_iteration`` (a pure function); this is worker-LOCAL read caching, not
state shared ACROSS workers, and nothing is ever sent back except each
iteration's own result.

**Machine policy (permanent, not situational).** This node is shared with
other users and has no scheduler. The worker count NEVER auto-scales to the
node's CPU count — it is exactly ``execution.n_jobs`` from config (default
16), capped only DOWN to the number of iterations actually planned, so a
small run does not spawn idle worker processes. Workers run niced
(``execution.worker_nice``, default 10, validated positive by
:class:`hotspot3d.robustness.params.RobustnessParams`) so they yield CPU to
other work on the node automatically — the default is safe to run at any
time without checking node load first. ``execution.n_jobs`` and
``execution.worker_nice`` are both INFRA_CONFIG (worker/thread counts —
:mod:`hotspot3d.utils.qa_taxonomy`'s textbook example), never a scientific
decision; changing them changes wall-clock time only, never a result.

**Thread oversubscription guard.** Each worker pins its own BLAS/OpenMP
thread pools to 1 via both environment variables (``OMP_NUM_THREADS`` etc.,
set before this worker's first task) and ``threadpoolctl`` (which reconfigures
whatever BLAS/OpenMP runtime is ALREADY loaded via its C API, so it works
regardless of whether the pool's start method already imported numpy/scipy
before the environment variables could take effect). Without this, N worker
processes each spawning their own M-thread BLAS pool would oversubscribe the
node by a factor of N*M against the caller's own thread budget.
"""
from __future__ import annotations

import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor

from ..footprint.api import CenterSet, Universe
from ..footprint.params import FootprintParams
from .perturb import IterationResult, run_iteration

#: Populated ONCE per WORKER PROCESS by _init_worker, before that worker
#: accepts its first task, and never mutated afterward. This is worker-LOCAL
#: state (one independent copy per OS process), never shared ACROSS workers —
#: see the module docstring.
_worker_inputs: dict | None = None

_THREAD_ENV_VARS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")


def _init_worker(inputs: dict, nice_value: int) -> None:
    """Runs once per worker process, before it accepts any task."""
    for var in _THREAD_ENV_VARS:
        os.environ[var] = "1"
    try:
        import threadpoolctl
        threadpoolctl.threadpool_limits(1)
    except ImportError:                     # pragma: no cover - always present here
        pass
    if nice_value:
        try:
            os.nice(nice_value)
        except OSError:                     # pragma: no cover - platform-dependent
            pass
    global _worker_inputs
    _worker_inputs = inputs


def _run_one(task: tuple[int, tuple[int, ...], str, str]) -> IterationResult:
    """Executed IN the worker process; reads only this worker's own copy."""
    index, removed, mode, seed_context = task
    s = _worker_inputs
    return run_iteration(
        index, removed, mode, seed_context,
        centers=s["centers"], universe=s["universe"], params=s["params"],
        original_state=s["original_state"], original_residues=s["original_residues"],
        r_fp_original=s["r_fp_original"], plp=s["plp"], blb=s["blb"],
        r_hot_value=s["r_hot_value"], undefined_mcc=s["undefined_mcc"])


def run_iterations(tasks, *, n_workers: int, worker_nice: int, centers: CenterSet,
                   universe: Universe, params: FootprintParams, original_state,
                   original_residues: frozenset, r_fp_original: float,
                   plp: frozenset, blb: frozenset, r_hot_value: float,
                   undefined_mcc: float) -> list[IterationResult]:
    """Run every ``(index, removed, mode, seed_context)`` task in ``tasks``.

    ``n_workers <= 1`` (or fewer than 2 tasks) runs serially in THIS process —
    no pool at all. This is both the simplest path and the reference the
    parallel paths are verified against (identical results at any worker
    count is the whole point, not merely faster ones).

    Never auto-scales to the node's CPU count: the effective worker count is
    ``min(n_workers, len(tasks))``, capped DOWN only to avoid spawning idle
    worker processes for a run with fewer iterations than workers.
    """
    tasks = list(tasks)
    if n_workers <= 1 or len(tasks) <= 1:
        return [
            run_iteration(index, removed, mode, seed_context,
                          centers=centers, universe=universe, params=params,
                          original_state=original_state,
                          original_residues=original_residues,
                          r_fp_original=r_fp_original, plp=plp, blb=blb,
                          r_hot_value=r_hot_value, undefined_mcc=undefined_mcc)
            for index, removed, mode, seed_context in tasks
        ]

    effective_workers = min(n_workers, len(tasks))
    inputs = dict(centers=centers, universe=universe, params=params,
                  original_state=original_state, original_residues=original_residues,
                  r_fp_original=r_fp_original, plp=plp, blb=blb,
                  r_hot_value=r_hot_value, undefined_mcc=undefined_mcc)

    # spawn, not the platform default (fork on Linux): a forked child can
    # inherit an already-initialized, multi-threaded BLAS pool from the
    # parent in a broken state (the classic fork-after-threads hazard).
    # spawn gives every worker a genuinely fresh interpreter, so the thread
    # pinning above reliably precedes this worker's first BLAS call.
    mp_ctx = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=effective_workers, mp_context=mp_ctx,
                             initializer=_init_worker,
                             initargs=(inputs, worker_nice)) as pool:
        return list(pool.map(_run_one, tasks))
