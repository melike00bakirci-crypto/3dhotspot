#!/usr/bin/env bash
# Authoritative full-suite runner.
#
# Guarantees, in order:
#   1. exactly one full suite runs (advisory lock, enforced in conftest.py too);
#   2. conservative, EXPLICIT thread caps (no reliance on a default);
#   3. a BEFORE fingerprint over every authoritative file — exclusion-based, so a
#      file nobody thought to name cannot slip through;
#   4. complete unbuffered output to a persistent log, never piped through
#      tail/head/grep;
#   5. an AFTER fingerprint, and a hard failure if the tree moved during the run.
#
# Usage: scripts/run_definitive_suite.sh [pytest args...]
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

LOG_DIR="logs"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG="${LOG_DIR}/definitive_suite_${STAMP}.log"
LATEST="${LOG_DIR}/definitive_suite.log"
FP_BEFORE="${LOG_DIR}/fingerprint_before_${STAMP}.txt"
FP_AFTER="${LOG_DIR}/fingerprint_after_${STAMP}.txt"
mkdir -p "$LOG_DIR"

# --- conservative thread caps, set EXPLICITLY so they are never defaulted -----
THREADS="${HOTSPOT3D_TEST_THREADS:-2}"
export HOTSPOT3D_TEST_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       MKL_NUM_THREADS="$THREADS" NUMEXPR_NUM_THREADS="$THREADS" \
       VECLIB_MAXIMUM_THREADS="$THREADS"

fingerprint() {
  python - "$1" <<'PY'
import sys
sys.path.insert(0, ".")
from conftest import authoritative_files, tree_fingerprint, REPO_ROOT
digest, n = tree_fingerprint()
with open(sys.argv[1], "w") as fh:
    fh.write(f"{digest}  {n} files\n")
    for p in authoritative_files():
        fh.write(f"{p.relative_to(REPO_ROOT)}\n")
print(f"{digest}  ({n} authoritative files)")
PY
}

echo "=== refusing to start if another full suite holds the lock ==="
if pgrep -f "python -u -m pytest tests/" >/dev/null 2>&1; then
  echo "ABORT: another full suite is already running." >&2
  exit 3
fi

echo "=== BEFORE fingerprint ==="
BEFORE="$(fingerprint "$FP_BEFORE")" || { echo "ABORT: fingerprint failed" >&2; exit 4; }
echo "$BEFORE"

echo "=== running (threads=${THREADS}, unbuffered, complete log) ==="
echo "log: $LOG"
python -u -m pytest tests/ -p no:cacheprovider -v -rA --tb=short --durations=25 \
  "$@" > "$LOG" 2>&1
EXIT=$?
echo "EXIT=${EXIT}" >> "$LOG"
ln -sf "$(basename "$LOG")" "$LATEST"

echo "=== AFTER fingerprint ==="
AFTER="$(fingerprint "$FP_AFTER")" || { echo "ABORT: fingerprint failed" >&2; exit 4; }
echo "$AFTER"

BEFORE_HASH="${BEFORE%% *}"
AFTER_HASH="${AFTER%% *}"
{
  echo
  echo "TREE_FINGERPRINT_BEFORE=${BEFORE_HASH}"
  echo "TREE_FINGERPRINT_AFTER=${AFTER_HASH}"
} >> "$LOG"

if [[ "$BEFORE_HASH" != "$AFTER_HASH" ]]; then
  echo "FREEZE VIOLATION: the tree changed during the run." >&2
  echo "  before: $BEFORE_HASH" >&2
  echo "  after : $AFTER_HASH" >&2
  diff <(tail -n +2 "$FP_BEFORE") <(tail -n +2 "$FP_AFTER") >&2 || true
  echo "FREEZE_VIOLATION=TRUE" >> "$LOG"
  exit 5
fi

echo "FREEZE_VIOLATION=FALSE" >> "$LOG"
echo "tree unchanged during the run — result is authoritative"
exit "$EXIT"
