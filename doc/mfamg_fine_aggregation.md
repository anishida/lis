# Optional first-level multivector MFAMG aggregates

By default, MFAMG retains its historical fine-level aggregate target of
**two local unknowns**, even if multiple near-nullspace modes are supplied.
That can leave many fine-level aggregates with local rank close to the
number of unknowns, making the first coarse space unnecessarily large.

## Option

Use the per-solver option:

```text
-i gmres -p mfamg -mfamg_fine_agg_factor 2
```

Accepted values are `0` (default, unchanged behavior) and `2`, `3`, `4`.
When enabled and `near_nullspace_dim > 1`, the first-level target size
becomes `min(local_n, factor * mode_count)`. The multiplication is
capped using division before evaluating the product to prevent
`LIS_INT` overflow.

MFAMG reuses its existing connected aggregate builder, local QR/rank
determination, Galerkin setup and physical convergence tests. The
option changes **only the first transfer**; the default transfer
and deep-level aggregation remain unchanged. For one near-nullspace
mode, the target stays at two.

The option is stored per solver, not in process-global environment
variables, and participates in transactional option-profile replacement
and change detection. Invalid values are rejected.

## Validation and limitations

The `testmfamgfineagg` regression uses a deterministic 240x240 SPD
tridiagonal matrix. It compares **the actual first coarse-level
unknown counts** for the default (`factor=0`), `factor=2`, and
`factor=3`, with both **two and three globally independent
near-nullspace modes**. Within blocks of eight fine unknowns, the
chosen modes are locally constant (rank one), permitting construction
of the default hierarchy; globally they remain independent.
The test also checks physical preconditioner defects, option-profile
change detection and rejection of invalid values.
Build and run via `make check`.

The regression exercises a one-rank numerical problem. MPI builds
are compiled and checked by the CI suite; parallel aggregation with
multiple ranks has not yet been independently proven to improve
convergence. No nonlinear-contact or solver-performance gain is claimed.

This is an **opt-in** capability. MFAMG is unchanged unless the
new option is explicitly specified.

## Degrees of freedom and numerical tradeoffs

An "aggregate size" here is a count of **local algebraic unknowns
(matrix rows / DOFs)**, not physical mesh nodes. In multiphysics or
multi-DOF finite-element discretizations, the two numbers may differ.

Increasing the aggregate target is not guaranteed to improve wall
time, iteration count, or interpolation quality. On unstructured
multidimensional graphs, a greedy connected aggregate can become
elongated, and its near-nullspace approximation properties can
worsen. The new option therefore remains disabled by default.

The numerical regression measures a first-level hierarchy-size
reduction on one structured SPD example. It is **not** a benchmark
of 3D FE matrices or contact problems.
