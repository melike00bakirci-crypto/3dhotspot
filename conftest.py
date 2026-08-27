"""Repository-root pytest configuration — TEST EXECUTION ONLY.

Nothing here touches scientific methodology. No formula, threshold, seed,
permutation count or expected output is affected by this file; it governs how the
test *process* is run, not what the code computes.

Two problems are solved here, both observed in practice on this machine (256 cores):

1. **BLAS/OpenMP oversubscription.** NumPy's OpenBLAS opens one thread per core by
   default, so a single pytest process was holding ~127 threads. Five concurrent
   suites meant ~635 threads contending for 256 cores and suites that should take
   minutes ran for over 90 minutes. The caps below are set BEFORE NumPy is imported
   — this file is the earliest conftest pytest loads, which is the only point at
   which OpenBLAS reads them.

   The cap is deliberately conservative and FIXED, not tuned per machine, so test
   timings are comparable across runs and machines.

   Results are unaffected: the hot matrix products in the permutation null are
   products of 0/1 matrices whose partial sums are exact integers well below 2**53,
   so no reordering or threading change can alter a single bit of the result. This
   was verified by running the same stage at 1 thread and uncapped and comparing
   r_hot, |S| and the permutation-resolution flag.

2. **Concurrent full suites.** A full-suite invocation takes an advisory lock, so a
   second full suite fails fast with a clear message instead of silently competing
   for CPU. Targeted runs (``pytest tests/integration``, ``pytest -k ...``) are NOT
   locked, so per-agent iteration stays unblocked.

Override the caps deliberately with ``HOTSPOT3D_TEST_THREADS=<n>`` when
benchmarking; the default is what CI and release validation use.
"""
from __future__ import annotations

import os

# --- 1. thread caps: MUST run before numpy/scipy are imported ----------------
TEST_THREAD_CAP = os.environ.get("HOTSPOT3D_TEST_THREADS", "1")

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, TEST_THREAD_CAP)

import fcntl                                                        # noqa: E402
import sys                                                          # noqa: E402
from pathlib import Path                                            # noqa: E402

import pytest                                                       # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent
LOCK_PATH = REPO_ROOT / ".pytest-full-suite.lock"

_lock_handle = None


def _is_full_suite(config) -> bool:
    """True when the invocation collects the whole ``testpaths`` root.

    Uses ``config.args`` — the resolved collection targets — rather than the raw
    command line, so that ``pytest`` and ``pytest tests/`` are both recognised as
    full suites while ``pytest tests/integration`` is not. Reading the raw argv
    would misclassify option VALUES (``-p no:cacheprovider``) as test paths.
    """
    testpaths = config.getini("testpaths") or []
    if not testpaths:
        return False
    roots = {(Path(config.rootpath) / p).resolve() for p in testpaths}
    targets = {(Path(config.rootpath) / str(a).split("::")[0]).resolve()
               for a in config.args}
    return bool(targets) and targets == roots


def pytest_configure(config):
    global _lock_handle
    config.addinivalue_line(
        "markers", "slow: long-running scientific validation; excluded by -m 'not slow'")

    if not _is_full_suite(config):
        return

    _lock_handle = open(LOCK_PATH, "w")
    try:
        fcntl.flock(_lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        _lock_handle.close()
        _lock_handle = None
        raise pytest.UsageError(
            f"ANOTHER FULL PYTEST SUITE IS ALREADY RUNNING (lock: {LOCK_PATH}).\n"
            f"Concurrent full suites oversubscribe the CPU and make every run slower "
            f"and its timings meaningless. Wait for the running suite, or scope this "
            f"run to a path (e.g. `pytest tests/hotspot_stats`), which is not locked."
        )
    _lock_handle.write(f"pid={os.getpid()}\n")
    _lock_handle.flush()


def pytest_unconfigure(config):
    global _lock_handle
    if _lock_handle is not None:
        fcntl.flock(_lock_handle.fileno(), fcntl.LOCK_UN)
        _lock_handle.close()
        _lock_handle = None


#: Files whose content defines an authoritative run. Exclusion-based ON PURPOSE:
#: an inclusion list silently omits whatever nobody thought to name, and its
#: silence is indistinguishable from a clean result. This list previously missed
#: ./conftest.py, ./pyproject.toml (whose `filterwarnings = ["error::RuntimeWarning"]`
#: can flip pass/fail suite-wide) and the frozen methodology document itself.
FINGERPRINT_SUFFIXES = (".py", ".yaml", ".yml", ".toml", ".md", ".cfg", ".ini")
FINGERPRINT_EXCLUDE = (".venv", "__pycache__", ".pytest_cache", ".git",
                       "results", "logs", "htmlcov")


def authoritative_files(root: Path | None = None) -> list[Path]:
    """Every file whose content could change what a run computes or collects."""
    root = root or REPO_ROOT
    out = []
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix not in FINGERPRINT_SUFFIXES:
            continue
        if any(part in FINGERPRINT_EXCLUDE for part in path.relative_to(root).parts):
            continue
        out.append(path)
    return sorted(out)


def tree_fingerprint(root: Path | None = None) -> tuple[str, int]:
    """(sha256 over all authoritative file contents + relative paths, file count)."""
    import hashlib
    root = root or REPO_ROOT
    digest = hashlib.sha256()
    files = authoritative_files(root)
    for path in files:
        digest.update(str(path.relative_to(root)).encode())
        digest.update(hashlib.sha256(path.read_bytes()).digest())
    return digest.hexdigest(), len(files)


def pytest_report_header(config):
    import numpy
    # Report the ACTUAL environment, not the value this file intended to set.
    # os.environ.setdefault does not override a caller-supplied cap, so echoing
    # TEST_THREAD_CAP misreported the real thread count on every run.
    actual = {v: os.environ.get(v, "<unset>") for v in
              ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
    caps = ",".join(sorted(set(actual.values())))
    fingerprint, n_files = tree_fingerprint()
    return (
        f"hotspot3d: thread caps actual={caps} "
        f"(requested={TEST_THREAD_CAP}) detail={actual}\n"
        f"hotspot3d: numpy={numpy.__version__} python={sys.version.split()[0]}\n"
        f"hotspot3d: full-suite lock="
        f"{'HELD' if _lock_handle else 'not taken (targeted run)'}\n"
        f"hotspot3d: tree fingerprint={fingerprint} ({n_files} authoritative files)"
    )
