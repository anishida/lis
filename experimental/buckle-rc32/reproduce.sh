#!/usr/bin/env bash
# Use after applying the consolidated patch and building the original LIS 2.1.25 source.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
tree=$(realpath "${1:?configured and built LIS tree}")
data=$(realpath "${2:?private CSR directory (A,B,S,d.npy etc.)}")
rhs=$(realpath "${3:?private Lanczos RHS binary file}")
out=$(realpath -m "${4:-$PWD/reproduction}"); mkdir -p "$out"
cc=${CC:-gcc}; lib="$tree/src/.libs/liblis.a"
for name in native_full integration_test reuse_test; do
 "$cc" -O3 -fopenmp -I"$tree/include" "$here/tests/$name.c" "$lib" -lm -o "$out/$name"
done
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
export LIS_EXPERIMENT_ILUTP_SYMBOLIC=1 LIS_EXPERIMENT_ILUTP_PAIRS=1 LIS_EXPERIMENT_ILUTP_GROUPS=1 LIS_EXPERIMENT_ILUTP_GROUP8=1 LIS_EXPERIMENT_ILUTP_SERIAL_APPLY=1 LIS_EXPERIMENT_ILUTP_DOT4=1 LIS_EXPERIMENT_ILUTP_MD=1 LIS_EXPERIMENT_ILUTP_REUSE=1
unset LIS_EXPERIMENT_ILUTP_PARALLEL_NUMERIC LIS_EXPERIMENT_ILUTP_RUNS
export LIS_BENCH_LINEAR='-i fgmres -p ilut -ilutp 1 -ilutp_threshold 0 -ilutp_full 1 -ilutp_full_limit 45000000 -tol 1e-10 -maxiter 40 -restart 35 -print none'
/usr/bin/time -f 'TIME elapsed=%e rss_KiB=%M exit=%x' "$out/native_full" "$data/A.csr" "$data/B.csr" "$out/native.bin" -e gli -etransform buckling -shift .00107 -ss 52 -nev 6 -etol 1e-8 -emaxiter 100 -omp_num_threads 4 -eprint out > "$out/native.log" 2>&1
python "$here/certify_portable.py" "$out/native.bin" "$data" "$here/reference_eigenvalues.csv" "$out/native_certified.csv"
export LIS_ACTUAL_RHS="$rhs"
opts=(-i fgmres -p ilut -ilutp 1 -ilutp_threshold 0 -ilutp_full 1 -ilutp_full_limit 45000000 -omp_num_threads 4 -print none -tol 1e-10 -maxiter 40 -restart 35)
for cache in on off; do
 if [[ $cache == off ]]; then unset LIS_EXPERIMENT_ILUTP_REUSE; fi
 /usr/bin/time -f 'TIME elapsed=%e rss_KiB=%M exit=%x' "$out/integration_test" 1 "$data/S0.0009.csr" "$data/S0.00105.csr" "$data/S.csr" "${opts[@]}" > "$out/sequence_$cache.log" 2>&1
done
