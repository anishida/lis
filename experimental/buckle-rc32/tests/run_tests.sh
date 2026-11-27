#!/usr/bin/env bash
# Run against a configured/built LIS 2.1.25 tree with patches 01..04 applied.
set -euo pipefail
src=$(realpath "${1:?path to configured and built LIS tree}")
here=$(cd "$(dirname "$0")" && pwd)
out=${2:-"$here/rerun"}
mkdir -p "$out"; out=$(realpath "$out")
cc=${CC:-gcc}
lib="$src/src/.libs/liblis.a"
export OMP_NUM_THREADS=1 LIS_MFAMG_SPARSE_TERMINAL=1
"$cc" -O3 -fopenmp -I"$src/include" "$here/test_inner_basis.c" "$lib" -Wl,--wrap=lis_precon_create -lm -o "$out/test_inner_basis"
for transform in buckling shift_invert; do
 for method in gli gii; do
  for pre in mfamg mfamg_ilutp; do
   for b in 1 2; do
    if [[ $method == gii && $b == 2 ]]; then continue; fi
    name="${method}_${pre}_b${b}"; if [[ $transform == shift_invert ]]; then name+='_si'; fi
    timeout 30 "$out/test_inner_basis" "$b" "$pre" "$method" "$transform" > "$out/$name.log" 2>&1
    grep -q '^PASS$' "$out/$name.log"
    echo "$name PASS"
   done
  done
 done
done
for ss in 3 4; do
 TEST_SS=$ss timeout 30 "$out/test_inner_basis" 1 mfamg gli > "$out/refine_ss$ss.log" 2>&1
 grep -q '^PASS$' "$out/refine_ss$ss.log"
 test "$(grep -c '^OBS' "$out/refine_ss$ss.log")" -ge 6
 echo "refine_ss$ss PASS"
done
for name in testblockgram testblockorthogonal testblockrank testprojectedeigen testesolverselection testgesolvepair testeigenbreakdown; do
 "$cc" -O3 -fopenmp -I"$src/include" "$src/test/$name.c" "$lib" -lm -o "$out/$name"
 timeout 30 "$out/$name" > "$out/$name.log" 2>&1
 echo "$name PASS"
done
