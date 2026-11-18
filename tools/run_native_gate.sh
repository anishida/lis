#!/usr/bin/env bash
# Native test gate to run INSIDE a real LIS checkout with CA-MFAMG PR #94.
# This script never commits, pushes or modifies GitHub state.
# Unlike the local mock tests, this one *requires* genuine LIS FGMRES.
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
  printf 'Usage: %s /path/to/lis-checkout [output-directory]\n' "$0" >&2
  exit 2
fi
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIS_ROOT="$(cd "$1" && pwd)"
REPORT_DIR="${2:-$LIS_ROOT/native-condense-gate}"
mkdir -p "$REPORT_DIR"
REPORT_DIR="$(cd "$REPORT_DIR" && pwd)"
if [[ ! -f "$LIS_ROOT/include/lis.h" || ! -f "$LIS_ROOT/src/precon/lis_precon_mfamg.c" ]]; then
  echo 'ERROR: Not a LIS checkout with MFAMG sources.' >&2
  exit 2
fi
# Refuse a misleading Ready result for a CA branch missing upstream changes.
if ! git -C "$LIS_ROOT" rev-parse --verify upstream/develop >/dev/null 2>&1; then
  echo 'ERROR: upstream/develop ref is missing; fetch latest upstream/develop first.' >&2
  exit 2
fi
if ! git -C "$LIS_ROOT" merge-base --is-ancestor upstream/develop HEAD; then
  echo 'ERROR: checkout HEAD does not contain upstream/develop; merge or rebase first.' >&2
  exit 2
fi
if ! grep -q 'lis_solver_set_mfamg_contact_partition' "$LIS_ROOT/include/lis.h"; then
  echo 'ERROR: This gate needs the CA-MFAMG branch/API from PR #94.' >&2
  exit 2
fi
export OMP_NUM_THREADS=1
python3 "$PROJECT_ROOT/tools/install_into_lis.py" "$LIS_ROOT" --regenerate \
  > "$REPORT_DIR/install.log" 2>&1
(
  cd "$LIS_ROOT"
  ./configure --enable-test --disable-mpi --disable-fortran --disable-f90 --disable-saamg > "$REPORT_DIR/configure.log" 2>&1
  make -j2 > "$REPORT_DIR/build.log" 2>&1
  make -C test testcondense_core testcomposite_core testcondense_lis \
      > "$REPORT_DIR/build-tests.log" 2>&1
  ./test/testcondense_core > "$REPORT_DIR/core.log" 2>&1
  ./test/testcomposite_core > "$REPORT_DIR/composite.log" 2>&1
  for precon in none mfamg ca-mfamg ilutp mfamg_ilutp ca-mfamg_ilutp; do
    ./test/testcondense_lis "$precon" > "$REPORT_DIR/lis-$precon.log" 2>&1
    # A mock/unsupported-configuration SKIP must not masquerade as native PASS.
    if ! grep -Fq "LIS updated Schur: $precon PASS" "$REPORT_DIR/lis-$precon.log"; then
      echo "ERROR: missing native two-step proof for $precon" >&2
      exit 1
    fi
  done
  make check > "$REPORT_DIR/make-check.log" 2>&1
  grep -q 'PASS: testca_mfamg_profiles.sh' "$REPORT_DIR/make-check.log" || {
    echo 'ERROR: six-profile regression missing from native Automake test suite' >&2
    exit 1
  }
  grep -q 'PASS: testcondense_differential.sh' "$REPORT_DIR/make-check.log" || {
    echo 'ERROR: stdlib differential regression missing from native Automake suite' >&2
    exit 1
  }
  grep -q 'PASS: testcondense_oom' "$REPORT_DIR/make-check.log" || {
    echo 'ERROR: fault-injected condensation regression missing from native Automake suite' >&2
    exit 1
  }
  grep -q 'PASS: testca_mfamg_ilutp.sh' "$REPORT_DIR/make-check.log" || {
    echo 'ERROR: contact-aware composite native Automake test missing or not PASS' >&2
    exit 1
  }
)
echo "Native LIS gate PASSED; logs in $REPORT_DIR"
