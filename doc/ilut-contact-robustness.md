# Optional ILUT reordering and row equilibration for FEA

The LIS ILUT preconditioner accepts two **opt-in** solver-local options:

- `-ilut_rcm 1` computes deterministic reverse Cuthill–McKee ordering on the **undirected sparsity graph** of the input CSR matrix. It permutes both rows and columns before factorization, and inversely permutes the triangular-solve result.
- `-ilut_row_scale 1` scales each reordered row by the reciprocal of its largest absolute coefficient before factorization. The corresponding scaling is applied to the right-hand side (and to the output of adjoint solves).

Both defaults remain zero. The factorization and direct `off` behavior are unchanged when neither option is requested. Use `-ilut_drop` and `-ilut_rate` as readable aliases for the existing `-iluc_drop` and `-iluc_rate` options (the same numerical parameter slots).

## Example

```sh
./test/testhalfilut
```

The regression compares ILUT normal and adjoint residuals under RCM, equilibration and their combination. For a custom solve:

```
-i fgmres -p ilut -scale none -ilut_rcm 1 -ilut_row_scale 1 \
-ilut_drop 0.01 -ilut_rate 10 -ilu_pivot_tol 1e-12
```

## Applicability and caveats

This first implementation supports **serial CSR in ordinary double precision**, including nonsymmetric sparse matrices. It explicitly rejects distributed, USER/OPERATOR and quad-precision configurations when the new transform options are requested. This is neither ddPQ (diagonal matching) nor distributed RAS/Schur. Matrix entries must be finite; a zero row or singular matrix may still produce a factorization failure.

A read-only internal diagnostic `lis_precon_ilut_factor_stats` reports the stored L/U nonzero counts and minimum/maximum actual pivot magnitudes (the stored D values are reciprocals). The regression prints these statistics to support controlled fill and stability comparisons.

RCM is a fill/bandwidth heuristic. It does **not** guarantee lower Krylov iteration counts or greater stability for contact KKT saddle-point systems. Report true relative residual, pivot failures, factor nonzeros, setup time and process memory when benchmarking. Preserve the default ILUT path as a control. Near-singular rows should use `-ilu_pivot_tol` as an explicit policy rather than silently replacing pivots.

Hypre's ILU framework additionally offers local RCM, dual-threshold ILUT, ddPQ matching, RAS-ILUT and two-level Schur correction; these are distinct algorithms and should not be conflated with this opt-in serial preprocessing step: https://hypre.readthedocs.io/en/latest/solvers-ilu.html
