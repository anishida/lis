#!/usr/bin/env bash
set -euo pipefail
base=${1:?directory with test binaries and fixtures}; out=${2:?output folder}; mkdir -p "$out"; cd "$base"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export LIS_EXPERIMENT_ILUTP_SYMBOLIC=1 LIS_EXPERIMENT_ILUTP_PAIRS=1 LIS_EXPERIMENT_ILUTP_GROUPS=1 LIS_EXPERIMENT_ILUTP_GROUP8=1 LIS_EXPERIMENT_ILUTP_SERIAL_APPLY=1 LIS_EXPERIMENT_ILUTP_DOT4=1 LIS_EXPERIMENT_ILUTP_MD=1 LIS_EXPERIMENT_ILUTP_REUSE=1
opts=(-i fgmres -p ilut -ilutp 1 -ilutp_threshold 0 -ilutp_full 1 -ilutp_full_limit 45000000 -tol 1e-10 -maxiter 40 -restart 35 -print none)
printf 'build,test,case,status\n' > "$out/random_summary.csv"
for mode in omp serial; do
 for i in $(seq -w 0 11); do
  files=(fixtures/case${i}_{0,1,2}.csr)
  for kind in direct fresh sequence rebuild; do
   case $kind in direct) exe=reuse_test; arg=0;; fresh) exe=reuse_test; arg=1;; sequence) exe=integration_test; arg=1;; rebuild) exe=integration_test; arg=0;; esac
   if ./${exe}_${mode} "$arg" "${files[@]}" "${opts[@]}" > "$out/${mode}_${kind}_${i}.log" 2>&1; then result=PASS; else result=FAIL; fi
   printf '%s,%s,%s,%s\n' "$mode" "$kind" "$i" "$result" >> "$out/random_summary.csv"
  done
 done
 for n in 1 2 3 4; do
  if NEGATIVE=$n ./reuse_test_${mode} 0 fixtures/case08_{0,1}.csr "${opts[@]}" > "$out/${mode}_negative_${n}.log" 2>&1; then result=PASS; else result=FAIL; fi
  printf '%s,negative,%s,%s\n' "$mode" "$n" "$result" >> "$out/random_summary.csv"
 done
 for i in $(seq -w 2 11); do
  if ./integration_test_${mode} 0 fixtures/case${i}_0.csr fixtures/case${i}_changed.csr "${opts[@]}" > "$out/${mode}_changed_${i}.log" 2>&1; then result=PASS; else result=FAIL; fi
  printf '%s,changed_pattern,%s,%s\n' "$mode" "$i" "$result" >> "$out/random_summary.csv"
 done
done
