#!/usr/bin/env bash
# Execute the REAL native FGMRES + CA-MFAMG + native ILUTP regression.
# Automake considers exit 77 a genuine SKIP in unsupported configurations.
set -euo pipefail
if [[ "${enable_mpi:-no}" == "yes" || "${enable_complex:-no}" == "yes" ]]; then
  echo "SKIP: CA-MFAMG + ILUTP condensation adapter requires serial real-double LIS"
  exit 77
fi
export OMP_NUM_THREADS=1
if [[ -n "${top_builddir:-}" ]]; then
  binary="${top_builddir}/test/testcondense_lis"
else
  binary="$(dirname "$0")/testcondense_lis"
fi
if [[ ! -x "$binary" ]]; then
  echo "FAIL: native CA+ILUTP test binary is missing: $binary" >&2
  exit 1
fi
"$binary" ca-mfamg_ilutp
