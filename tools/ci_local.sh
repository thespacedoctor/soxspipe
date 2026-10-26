#!/usr/bin/env bash
# APPROXIMATE THE CI TEST JOB LOCALLY: ONE COMMAND, EVERY GATE, ONE SUMMARY.
#
# THIS IS NOT CI. IT RUNS THE SAME GATES WITH THE SAME THREAD-COUNT SETTINGS AND PYTHONHASHSEED,
# BUT IT DOES NOT REPRODUCE:
#   - THE READ-ONLY soxspipe PACKAGE DIRECTORY (CI RUNS chmod -R a-w soxspipe)
#   - NETWORK ISOLATION (CI RUNS THE SUITE UNDER `unshare --net`; HERE ONLY conftest.py BLOCKS SOCKETS)
#   - THE RUNNER USER ID AND ITS CLEAN HOME DIRECTORY
#   - -n auto (CI USES ONE WORKER PER PHYSICAL CORE; HERE -n 4, THE LOCAL CAP SET FOR MEMORY)
# A GREEN RUN HERE IS A STRONG SIGNAL, NOT A GUARANTEE THAT CI PASSES.
#
# IT DELETES .coverage AND .coverage.* IN THE REPOSITORY ROOT BEFORE THE TESTS RUN, SO DO NOT
# KEEP COVERAGE DATA THERE THAT YOU WANT TO KEEP.
#
# USAGE (FROM ANYWHERE INSIDE THE CHECKOUT):
#   PYTHON=~/anaconda/envs/soxspipe/bin/python tools/ci_local.sh
#
# PYTHON IS THE INTERPRETER EVERY TOOL RUNS THROUGH (pytest, ruff, diff-cover, coverage), SO THEY
# ALL COME FROM ONE ENVIRONMENT. IT DEFAULTS TO THE python ON PATH, WHICH MUST BE PYTHON 3.12 WITH
# THE `tests` EXTRA INSTALLED.
#
# SETTINGS (ENVIRONMENT VARIABLES):
#   CI_LOCAL_BASE      BRANCH THE CHANGE IS JUDGED AGAINST. DEFAULT origin/develop
#   CI_LOCAL_WORKERS   pytest-xdist WORKERS. DEFAULT 4, THE LOCAL MEMORY CAP
#   CI_LOCAL_REPORTS   WHERE THE COVERAGE REPORT GOES. DEFAULT A NEW DIRECTORY UNDER $TMPDIR
#
# EVERY STEP RUNS EVEN IF AN EARLIER ONE FAILED. THE EXIT STATUS IS NONZERO IF ANY STEP FAILED.
# THE LINT RATCHET AND diff-cover JUDGE THE WORKING TREE, UNCOMMITTED AND UNTRACKED FILES INCLUDED.

set -euo pipefail

PYTHON="${PYTHON:-python}"
BASE="${CI_LOCAL_BASE:-origin/develop}"
WORKERS="${CI_LOCAL_WORKERS:-4}"
DIFF_COVER_FLOOR=80

# THE ENVIRONMENT tests.yml SETS ON THE PYTEST STEP: SINGLE-THREADED NUMERICS AND A FIXED HASH SEED
CI_ENVIRONMENT=(
    BLAS_NUM_THREADS=1
    MKL_NUM_THREADS=1
    NUMEXPR_NUM_THREADS=1
    OPENBLAS_NUM_THREADS=1
    OMP_NUM_THREADS=1
    VECLIB_MAXIMUM_THREADS=1
    PYTHONHASHSEED=0
)

cd "$(git rev-parse --show-toplevel)"
REPORTS="${CI_LOCAL_REPORTS:-$(mktemp -d "${TMPDIR:-/tmp}/soxspipe-ci-local.XXXXXX")}"
mkdir -p "$REPORTS"
COVERAGE_XML="$REPORTS/coverage.xml"

SUMMARY=()
FAILED=0

# RUN ONE STEP AND RECORD ITS RESULT. A FAILURE MUST NOT STOP THE SCRIPT, SO THE COMMAND SITS IN AN if
run_step() {
    local name="$1"
    shift
    echo
    echo "=== ${name} ==="
    if "$@"; then
        SUMMARY+=("PASS  ${name}")
    else
        SUMMARY+=("FAIL  ${name}")
        FAILED=1
    fi
}

# LEFTOVER STATEMENT-ONLY DATA MIXED WITH BRANCH DATA RAISES "Can't combine statement coverage data with branch data"
rm -f .coverage .coverage.*

run_step "lint ratchet (changed lines and hard rules)" \
    "$PYTHON" tools/lint_ratchet.py --worktree "$BASE"
run_step "changelog duplicates" \
    "$PYTHON" tools/check_changelog.py CHANGES.md
run_step "ruff format" \
    "$PYTHON" -m ruff format --check .
run_step "tests (-n ${WORKERS})" \
    env "${CI_ENVIRONMENT[@]}" "$PYTHON" -m pytest tests/unit tests/integration -m "not slow" -n "$WORKERS" \
    --cov=soxspipe --cov=tools --cov-branch --cov-report=term-missing:skip-covered \
    --cov-report="xml:${COVERAGE_XML}"
run_step "coverage floor" \
    "$PYTHON" -m coverage report
run_step "diff-cover >= ${DIFF_COVER_FLOOR}% against ${BASE}" \
    "$PYTHON" -m diff_cover.diff_cover_tool "$COVERAGE_XML" \
    --compare-branch="$BASE" --include-untracked --fail-under="$DIFF_COVER_FLOOR"

echo
echo "=== summary (reports in ${REPORTS}) ==="
printf '%s\n' "${SUMMARY[@]}"

exit "$FAILED"
