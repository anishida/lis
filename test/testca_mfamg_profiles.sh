#!/usr/bin/env bash
# Native end-to-end profiles, invoked by Automake's make check.
# No external solver is used. Never treat unsupported configurations as PASS.
set -euo pipefail
if [[ "${enable_mpi:-no}" == "yes" || "${enable_complex:-no}" == "yes" ]]; then
  echo 'SKIP: bridge supports serial real-double LIS only'
  exit 77
fi
export OMP_NUM_THREADS=1
if [[ -n "${top_builddir:-}" ]]; then
  binary="${top_builddir}/test/testcondense_lis"
else
  binary="$(dirname "$0")/testcondense_lis"
fi
if [[ ! -x "$binary" ]]; then
  echo "FAIL: native test executable missing: $binary" >&2
  exit 1
fi
for profile in none mfamg ca-mfamg ilutp mfamg_ilutp ca-mfamg_ilutp; do
  echo "NATIVE_PROFILE: $profile"
  "$binary" "$profile"
done
echo 'NATIVE_PROFILES: all six passed'
