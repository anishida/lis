# Reproducible ILUT Matrix Market review benchmark

This change extends the existing optional **`test/testhalfilut --fea-matrix`**
harness. It does not add a new solver or change default ILUT behavior.

## Build and run

The existing test executable is built with the LIS half-enabled test profile:

```sh
./configure --enable-test --enable-half --enable-omp --disable-mpi
make -j2
OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE \
  LIS_FEA_BASELINE_THREADS=1 LIS_FEA_RCM_THREADS=8 \
  ./test/testhalfilut --fea-matrix test/testmat0.mtx 8 8 0.005
```

The original CLI `--fea-matrix MATRIX BASE_RATE RCM_RATE [DROP]` and output
variant names `baseline` and `rcm_scaled` are preserved. The optional
environment variables allow per-variant **OpenMP factorization and
application** thread counts; if unset, each variant uses the process's
existing OpenMP thread count. In a non-OpenMP build, only a thread count of
1 is accepted.

The harness constructs a deterministic known vector `x_true` and sets
`b = A x_true`. Both variants operate on the original imported matrix,
not a reordered or row-scaled surrogate. Reported **true relative
residual** is `||Ax-b||_2 / ||b||_2` on that original system.
`rcm_scaled` enables the opt-in ILUT RCM ordering and row equilibration
already provided by upstream. The test is limited to single-process
real-valued CSR matrices.

## Interpreting the results

CSV contains the legacy columns and additionally:

- `setup_wall_s`, `solve_wall_s`, `total_wall_s`: **wall-clock**
  measurements, unlike the existing CPU-time counters. These timings are
  approximate and can vary with load and OpenMP scheduling.
- `threads`: actual thread budget selected for each variant.
- `target_1e-8_met`, `solver_error`, `status`: distinguish a completed
  solver run from satisfying the configured `1e-8` tolerance.

**Both variant rows are printed even if either fails**. Nonconvergence
returns a nonzero process status when its true residual exceeds the legacy
`1e-7` smoke acceptance band, but a residual between `1e-8` and
`1e-7` is still explicitly reported as `not_converged` rather than
misrepresented as meeting the target. File and matrix-import failures
are reported separately.

`fair_fill=yes` means the stored L/U **entry-count proxy** differs by at
most 5%; it does **not** mean equal peak RSS. If factor statistics are
missing, fill matching is reported as `unknown`, never guessed. For fair
comparisons, tune `BASE_RATE` and `RCM_RATE` independently and measure
multiple runs. A green CI smoke check is not a performance claim.

No private contact or KKT data are bundled, and a synthetic or standard
structural matrix is not evidence of contact-KKT performance.
