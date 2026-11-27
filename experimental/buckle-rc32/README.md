# LIS BUCKLE RC32 (stages 32–34)

Recommended experimental successor to RC31 on the original uploaded LIS 2.1.25. No changes to GitHub. Not a port to current develop. The consolidated patch is against pristine LIS; do not apply the incremental patches on top of it. It retains RC31's experimental ABI change from the inner near-nullspace prototype; rebuild callers. The same restrictions as RC31 apply (see RC31_background.md).

## Build

```
unzip lis-2.1.25.zip
cd lis-2.1.25
patch -p1 < /path/to/LIS_2.1.25_RC32_consolidated.patch
env -u TARGET ./configure --enable-omp --disable-shared
make -j4
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 make -j4 check
```
For the serial build omit `--enable-omp`. Restore executable permissions after ZIP extraction if needed (see RC31_background.md). The new consolidated patch was applied to copies of all 12 pristine affected files, and the results compared byte-for-byte with the RC32 tested source files. The full RC32 builds used RC31 plus incremental patch 20; RC33 was built separately. Standard regression: 125 PASS, 1 SKIP, 0 failures in each build; 124 additional checks across both builds (small random matrices, changed pattern, invalid data, reuse, rebuild, fresh factors). Logs in results/. MPI, wide indices, complex, production OOM and full sanitizer coverage were not tested.

## Recommended changes

Patch 20 shares the numeric row elimination kernel between initial and reused LU, and enables the experimental cache in serial builds. No meaningful time improvement measured over RC31. For RC32 use the RC31 experimental environment switches and solver parameters; the existing reproduce.sh is inherited from RC31 (requires private input data). For the full GLI benchmark, `-restart 2` is an optional lower-memory alternative to `-restart 35`; both returned identical binary eigenvectors on this model, but this is not a universal recommendation for incomplete preconditioners. Restart 1 also met certification but produced different vectors. All timings are process wall time, not repeated statistics except where multiple individual runs are listed.

## Rejected experiment

Patch 21 is NOT included in the consolidated patch. It tries direct full-LU correction in FGMRES and accepts only if the physical residual meets tolerance. In two full runs it accepted 29 of 52 solves, rejected 23, and was slower than RC32 (13.49/13.59 s vs 12.01/11.85 s). It is provided for inspection only, not recommended. No guarantee outside tested configurations.

## Files

summary.csv: timings, status, certification and peak RSS; vector_hashes.txt: hashes of private binary vectors (vectors not included). results/: standard and targeted regressions and original benchmark logs; generate_fixtures.py and run_random.sh regenerate synthetic test matrices and run tests. Tests and reproduction scripts inherited from RC31 require a built LIS tree and, for full benchmarks, private matrices. To build small test binaries for run_random.sh use the test drivers in tests/ and link static liblis.a, with -fopenmp for the OpenMP build and without it for the serial build. No private model matrix/RHS/vector files in this archive. The benchmark driver emits an out-of-range mode diagnostic after retrieving the six requested modes; solve returned SUCCESS and independent certificate verified all six. The test driver itself returns zero regardless of solver status, so the certificate and RESULT status must be checked.
