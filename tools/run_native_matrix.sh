#!/usr/bin/env bash
# Stage 10 multi-config native LIS Ready gate. Requires real local LIS checkout.
# Read-only with respect to GitHub. Builds outside the source checkout.
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 /path/to/integrated-lis-checkout [out-dir]" >&2
  exit 2
fi
lis="$(cd "$1" && pwd)"
out="${2:-$PWD/ca-mfamg-native-matrix-logs}"
mkdir -p "$out"; out="$(cd "$out" && pwd)"
for need in "$lis/configure" "$lis/include/lis_condense_lis.h" \
            "$lis/src/precon/lis_precon_mfamg_ilutp.c" \
            "$lis/test/testcondense_oom.c" "$lis/test/testca_mfamg_profiles.sh"; do
  if [[ ! -f "$need" ]]; then
    echo "NOT READY: missing integrated source: $need" >&2
    exit 2
  fi
done
if ! git -C "$lis" rev-parse --verify upstream/develop >/dev/null 2>&1 ||
   ! git -C "$lis" merge-base --is-ancestor upstream/develop HEAD; then
  echo 'NOT READY: updated upstream/develop is not an ancestor of HEAD' >&2
  exit 2
fi
if ! command -v autoreconf >/dev/null || ! command -v gcc >/dev/null ||
   ! command -v make >/dev/null; then
  echo 'NOT READY: missing compiler/Autotools' >&2
  exit 2
fi
# Build artifacts always outside the source tree; no git commit/push.
(cd "$lis" && autoreconf -fi > "$out/autoreconf.log" 2>&1)
common=(--enable-test --enable-oom-tests --disable-fortran --disable-f90 --disable-saamg)
configurations=(serial openmp int64 sanitized)
if command -v mpicc >/dev/null && command -v mpiexec >/dev/null; then
  configurations+=(mpi)
else
  echo 'MPI toolchain absent; MPI remains UNTESTED' > "$out/mpi.SKIPPED"
fi
for mode in "${configurations[@]}"; do
  build="$out/build-$mode"
  mkdir -p "$build"
  flags=("${common[@]}")
  cc=gcc
  cflags='-O2'
  ldflags=''
  case "$mode" in
    serial)    flags+=(--disable-mpi --disable-omp);;
    openmp)    flags+=(--disable-mpi --enable-omp);;
    int64)     flags+=(--disable-mpi --disable-omp --enable-longlong);;
    sanitized) flags+=(--disable-mpi --disable-omp)
               cflags='-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined'
               ldflags='-fsanitize=address,undefined';;
    mpi)       flags+=(--enable-mpi --disable-omp)
               cc=mpicc;;
  esac
  echo "BEGIN $mode"
  (cd "$build" && \
     CC="$cc" CFLAGS="$cflags" LDFLAGS="$ldflags" \
       "$lis/configure" "${flags[@]}" > configure.log 2>&1 && \
     make -j2 > build.log 2>&1 && \
     OMP_NUM_THREADS=1 make -j2 check > check.log 2>&1) || {
        echo "NOT READY: $mode build or check failed. See $build" >&2
        exit 1
     }
  # No false green: every serial/OMP/64/sanitized build must RUN all profiles
  # and injected-OOM regression, not merely compile their binaries.
  grep -Eq 'PASS: testcondense_differential.sh( |$)' "$build/check.log" || {
    echo "NOT READY: $mode independent differential regression missing" >&2; exit 1; }
  if [[ "$mode" != mpi ]]; then
    grep -Eq 'PASS: testca_mfamg_profiles.sh' "$build/check.log" || {
      echo "NOT READY: $mode six-profile test not PASS" >&2; exit 1; }
    grep -Eq 'PASS: testca_mfamg_ilutp.sh' "$build/check.log" || {
      echo "NOT READY: $mode CA+ILUTP test not PASS" >&2; exit 1; }
    grep -Eq 'PASS: testcondense_oom' "$build/check.log" || {
      echo "NOT READY: $mode injected-OOM test not PASS" >&2; exit 1; }
  else
    grep -Eq 'SKIP: testca_mfamg_profiles.sh' "$build/check.log" || {
      echo 'NOT READY: MPI bridge must explicitly SKIP' >&2; exit 1; }
  fi
  echo "PASS $mode"
done
if [[ -e "$out/mpi.SKIPPED" ]]; then
  echo 'PARTIAL PASS: MPI not tested; do not claim all configurations Ready.'
else
  echo 'NATIVE MATRIX PASSED: serial/OpenMP/LIS_INT64/sanitizers/MPI checks.'
fi
