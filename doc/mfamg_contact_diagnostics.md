# Opt-in MFAMG transfer diagnostics for contact / elasticity

Set `LIS_MFAMG_DIAG=1` in the environment of the application using LIS,
for example:

```sh
LIS_MFAMG_DIAG=1 ./test/testmfamg
```

This logs **one line per successfully constructed multivector tentative
transfer** and **one line per supplied near-nullspace mode**. Output goes to
`stderr`, and the default behavior is unchanged. It does not change
aggregation, rank selection, coarse solves, convergence or exit status.

## Log records

- `LIS_MFAMG_DIAG_TRANSFER n=... aggregates=... coarse_n=... ratio=...`
  includes histogram bins `singles`, `pairs`, `size3to8`,
  `size9plus`, `rank0`, `rank1`, `rank2`, `rank3`,
  `rank4plus`, and `poor_coarsening`.
  The ratio is `coarse_n / n`. The heuristic flag is 1 if the ratio
  exceeds 0.9. It is not a solver failure criterion.
- `LIS_MFAMG_DIAG_MODE n=... mode=... input_norm2=... relerr=...`
  reports the relative Euclidean error
  `||B_m - P Bc_m||_2 / ||B_m||_2`, where `Bc_m` is the tentative
  coarse basis constructed by the existing local QR procedure.
  The zero-input case uses zero error if its reconstruction is also zero.

These are **per-transfer, per-MPI-rank local** measurements, not
global reductions. In MPI, records from ranks can interleave. For serial
CalculiX, `n` identifies the input dimension of each transfer.
Rank bins report *algebraic* local ranks, not physical node blocks.

This diagnostic currently covers **multivector** MFAMG
(`near_nullspace_dim > 1`). Scalar near-nullspace uses a separate
prolongation path, which is not instrumented here.

## Interpretation

- A large `relerr` for a relevant rigid-body mode indicates that the
  tentative transfer does not span that supplied vector. This says
  nothing by itself about whether the mode is an exact null vector of
  the current contact tangent matrix.
- A large singleton count or `ratio > 0.9` indicates weak coarsening,
  but should be examined together with iteration histories.
- For contact, preserve the identity of the structural bodies before
  activating penalty springs; the matrix graph alone may merge them.
  The diagnostic cannot infer bodies or DOF-to-node maps.

## Limitations and cost

The optional projection loop costs `O(n * modes * local_rank)` floating
operations in addition to setup, and prints `modes + 1` lines per transfer.
It uses already allocated arrays and performs **no new allocations**.
Do not enable this diagnostic for production throughput without assessing
setup overhead and log volume.

This is intentionally separate from the experimental
MFAMG terminal block-Schur PR #201. It does not claim to solve the
CalculiX contact convergence problems by itself.

The Linux CI exercises the diagnostics using `testmfamg` and verifies
that both transfer and mode records are emitted.
