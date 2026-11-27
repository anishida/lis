#!/usr/bin/env bash
# Run the original RC32 LU test drivers against one configured LIS build.
# The historical test_inner_basis.c is intentionally excluded (patch 04).
set -euo pipefail
src=$(realpath "${1:?configured static LIS build tree}")
mode=${2:?serial or omp}
out=$(realpath -m "${3:?absolute results directory}")
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
cd "$out"
cc=${CC:-gcc}
flags=()
case "$mode" in serial) ;; omp) flags+=(-fopenmp);; *) exit 2;; esac
# Respect instrumented build flags if passed by the CI caller.
read -r -a cflags <<< "${RC32_TEST_CFLAGS:-}"
read -r -a ldflags <<< "${RC32_TEST_LDFLAGS:-}"
for driver in reuse_test integration_test; do
  "$cc" -O2 -g "${flags[@]}" "${cflags[@]}" -I"$src/include" "$here/tests/$driver.c" "$src/src/.libs/liblis.a" "${ldflags[@]}" -lm -o "$driver"
done
python3 "$here/generate_fixtures.py" fixtures
export OMP_DYNAMIC=FALSE
# Deliberately disable experimental parallel numeric factorization;
# the recommended RC32 path uses serial numeric LU and serial apply.
unset LIS_EXPERIMENT_ILUTP_PARALLEL_NUMERIC LIS_EXPERIMENT_ILUTP_RUNS LIS_ACTUAL_RHS
export LIS_EXPERIMENT_ILUTP_SYMBOLIC=1 LIS_EXPERIMENT_ILUTP_PAIRS=1 LIS_EXPERIMENT_ILUTP_GROUPS=1 LIS_EXPERIMENT_ILUTP_GROUP8=1 LIS_EXPERIMENT_ILUTP_SERIAL_APPLY=1 LIS_EXPERIMENT_ILUTP_DOT4=1 LIS_EXPERIMENT_ILUTP_MD=1 LIS_EXPERIMENT_ILUTP_REUSE=1
opts=(-i fgmres -p ilut -ilutp 1 -ilutp_threshold 0 -ilutp_full 1 -ilutp_full_limit 45000000 -tol 1e-10 -maxiter 40 -restart 35 -print none)
printf 'build,threads,test,case,status\n' > summary.csv
threads=(1)
if [[ $mode == omp ]]; then threads+=(4); fi
failed=0
run_case() {
  local kind=$1 id=$2 exe=$3 arg=$4; shift 4
  local log="${mode}_t${OMP_NUM_THREADS}_${kind}_${id}.log" status=PASS
  if ! timeout 45 "./$exe" "$arg" "$@" "${opts[@]}" > "$log" 2>&1; then
    status=FAIL; failed=1
    cat "$log"
  fi
  printf '%s,%s,%s,%s,%s\n' "$mode" "$OMP_NUM_THREADS" "$kind" "$id" "$status" >> summary.csv
}
for nt in "${threads[@]}"; do
  export OMP_NUM_THREADS=$nt
  for i in $(seq -w 0 11); do
    files=(fixtures/case${i}_{0,1,2}.csr)
    run_case direct "$i" reuse_test 0 "${files[@]}"
    run_case fresh "$i" reuse_test 1 "${files[@]}"
    run_case sequence "$i" integration_test 1 "${files[@]}"
    run_case rebuild "$i" integration_test 0 "${files[@]}"
  done
  for n in 1 2 3 4; do
    NEGATIVE=$n run_case negative "$n" reuse_test 0 fixtures/case08_{0,1}.csr
  done
  for i in $(seq -w 2 11); do
    run_case changed_pattern "$i" integration_test 0 "fixtures/case${i}_0.csr" "fixtures/case${i}_changed.csr"
  done
done
cat summary.csv
# Reject missing cases and any driver failure.
expected=$((62 * ${#threads[@]}))
test "$(($(wc -l < summary.csv)-1))" -eq "$expected"
test "$failed" -eq 0
