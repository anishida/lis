# MFAMG: sparse terminal fallback for stalled hierarchies

MFAMG normally reduces the hierarchy to a small coarse matrix and uses a
dense factorization at the last level. With multiple supplied near-nullspace
vectors, some matrices stop coarsening before the dense terminal limit.
For example, a 12-vector elasticity/contact basis sets the normal terminal
limit to 96 unknowns; a coarse matrix of thousands of rows cannot safely be
passed to that dense path.

## Default behavior

The default is unchanged. If the hierarchy does not reach its dense
terminal limit, setup fails with a diagnostic reporting the number of levels,
the terminal matrix size, and the limit. It no longer returns an unexplained
invalid-argument status from the dense terminal allocation.

## Experimental opt-in: serial sparse terminal

Set the following **before launching** the process:

```bat
set "LIS_MFAMG_SPARSE_TERMINAL=1"
set "LIS_MFAMG_DIAG=1"
set "OMP_NUM_THREADS=8"
set "OMP_DYNAMIC=FALSE"
set "CCX_LIS_NL_CONTACT_OPTIONS=-i fgmres -restart 120 -p mfamg -mfamg_block_size 4 -mfamg_fine_agg_factor 2 -mfamg_deep_sibling_merge 1 -mfamg_vcycle_smooth_steps 2 -scale none -omp_num_threads 8 -maxiter 1000"
```

When enabled, **serial CSR** multivector hierarchies stop after two consecutive
transfers each reduce the dimension by less than 10%. If the remaining
terminal matrix still exceeds the dense limit, MFAMG uses the existing
Weighted Jacobi or Block Jacobi smoother on that sparse terminal matrix,
with eight passes by default. To control the strength and cost of the approximate
terminal correction, you may **opt in** to 1–32 passes using
`LIS_MFAMG_SPARSE_TERMINAL_PASSES=16`, for example. The setting is read
once during preconditioner setup; an invalid value fails explicitly with
`LIS_ERR_ILL_ARG` if the sparse terminal is activated. It has no effect
on normal dense-terminal hierarchies and does not enable sparse terminal
automatically. The default remains eight passes, with no change in results
when the new variable is unset. The fallback avoids the O(n²) storage and
O(n³) factorization of a large dense terminal solve.

This is an **approximate** terminal correction. It does not preserve the
spectral properties of an exact coarse solve and can converge slowly or fail
for difficult contact matrices. An outer **FGMRES** solver is recommended;
always verify its true residual and the nonlinear Newton convergence.
The sparse-terminal fallback is **not implemented for multi-rank MPI**.
OpenMP inside a serial MPI rank remains independently configurable via
`LIS_MFAMG_OMP_MIN_ROWS` and `-omp_num_threads`.

Without the opt-in, existing MFAMG numerical paths and their coarse terminal
factorizations are unchanged.

## Diagnostics

The following messages distinguish hierarchy stagnation from Krylov failure:

- `MFAMG COARSENING STAGNATION` — consecutive weak reductions under the opt-in.
- `MFAMG SPARSE TERMINAL` — sparse terminal is active; reports final size,
  smoother, and the configured (default eight) terminal passes.
- `MFAMG HIERARCHY STALLED` — normal dense terminal is too small and the
  fallback has not been enabled.
- `LIS_MFAMG_DIAG_TRANSFER` — per-level rank and dimension statistics when
  `LIS_MFAMG_DIAG=1` is enabled.

The regression in `test/testmfamgvcycle.c` uses 36 disconnected 16-row
positive-definite chains and 12 near-nullspace vectors. It verifies that a
large sparse terminal is used instead of a forbidden dense factor and that
repeated corrections are finite, reduce the defect, and match bitwise across
the exercised OpenMP thread counts. It also checks that the default eight-pass
correction remains reproducible, and that an explicit 16-pass correction
strictly improves the manufactured SPD defect. A second test with a weaker
diagonal (2.25 instead of 4.0 and off-diagonal -1) requires at least a 1%
reduction in squared relative defect when increasing from 8 to 16 passes.
Both tests check that invalid pass counts fail closed. This synthetic
improvement is not evidence of better FGMRES or contact/Newton convergence.

## Limitations and follow-up

The fallback addresses a setup failure; it is **not** a replacement for
effective coarse-space construction. More robust algebraic aggregation,
true sparse coarse solves, and validation on exported CalculiX matrices
remain follow-up work. Do not infer any full-model speedup from the
stand-alone OpenMP smoother benchmarks.
