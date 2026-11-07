# Opt-in generic MFAMG near-nullspace transfer diagnostics

Set `LIS_MFAMG_DIAG=1` in the environment of any application using LIS:

```sh
LIS_MFAMG_DIAG=1 ./test/testmfamg
```

The diagnostic is **FEA-agnostic**: it consumes algebraic aggregates and
user-supplied near-nullspace vectors. It makes no assumptions about the
exporting FEA program, node numbering, element type, dimension or number of
DOFs per node. Rigid-body modes in 3D contact are a motivating example,
not a required input type.

The flag emits one `LIS_MFAMG_DIAG_TRANSFER` record for each constructed
multivector tentative transfer and one `LIS_MFAMG_DIAG_MODE` record per
input near-nullspace vector. Output goes to `stderr`, with no effect on
numerical setup, convergence decisions or exit status.

## Record fields

- Both kinds of record include `level`, `mpi_rank`, and local input
  dimension `n`. Level 0 is the finest transfer and higher levels are
  deeper hierarchy transfers. Stand-alone calls to the low-level builder
  use `level=-1` (no hierarchy context). Rank is zero in serial.
- Transfer records also provide `aggregates`, `coarse_n`,
  `ratio=coarse_n/n`, `modes`, aggregate size bins
  `singles/pairs/size3to8/size9plus`, algebraic rank bins
  `rank0/rank1/rank2/rank3/rank4plus`, and an advisory
  `poor_coarsening=1` when `ratio > 0.9`.
- Mode records report `mode` (zero-based vector index), squared
  Euclidean input norm `input_norm2`, and
  `relerr=||B_m-P Bc_m||_2/||B_m||_2` on this transfer.
  An exactly zero input reconstructed as zero has error zero.

The transfer QR is the real measured projection. These are local
Euclidean errors, not cumulative hierarchy errors or energy-norm
estimates. Large RBM errors signal loss of supplied modes, not necessarily
exact null vectors of the assembled operator. Singleton-dominated
aggregation and high `ratio` are warnings, not failure conditions.

## Applicability, runtime and limits

Diagnostics cover **multivector** MFAMG only (`near_nullspace_dim > 1`);
the scalar near-nullspace path is unchanged. Under MPI, all metrics are
rank-local; lines emitted from different processes may interleave. Use the
`mpi_rank` and `level` fields for log aggregation. The builder itself
requires no FE metadata.

When enabled, computing all per-mode projections costs roughly
`O(n * modes * local_rank)` and emits `modes + 1` records per transfer.
The diagnostic reuses allocated arrays and makes no new allocations.
Leave it disabled for production throughput unless diagnostic overhead
has been measured.

The smoke test checks the log schema, histogram consistency and
near-exact recovery of a manufactured two-mode basis. It is independent
of any FEA adapter and of the experimental terminal Block-Schur work.
