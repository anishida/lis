# CA-MFAMG: contact-aware algebraic aggregation (opt-in)

CA-MFAMG extends the existing matrix-free-friendly algebraic MFAMG hierarchy
without changing the behavior of `-p mfamg` or requiring mesh/FE metadata.

The application supplies **one partition label per locally owned algebraic
row (DOF)**. Equal labels permit MFAMG's normal strength-based aggregation;
different labels are *hard boundaries* when constructing every fine and
coarse aggregate. Unlike dropping contact terms or splitting the matrix,
this does **not** remove cross-partition/contact couplings from `A`,
`P^H A P`, the residual, or Krylov convergence checks. Every coarse
basis unknown inherits the partition label of its fine aggregate.

## API and usage

```c
/* A is the assembled LIS matrix; labels cover A->n locally owned rows.
 * The example assumes the application has already constructed the
 * near-nullspace vectors (e.g., appropriate rigid-body modes).
 */
LIS_INT *labels; /* application-supplied local DOF IDs, length A->n */
LIS_SOLVER solver; /* previously created */
/* ... create solver, A, b, x and modes ... */
lis_solver_set_near_nullspace(solver, mode_count, modes);
lis_solver_set_mfamg_contact_partition(solver, A->n, labels);
lis_solver_set_option(
    "-i fgmres -p ca-mfamg -scale none -print none", solver);
lis_solve(A, b, x, solver);
```

The setter copies labels into solver-owned memory. Its `nlocal` argument
**must match `A->n`** at setup. `(0, NULL)` clears the labels and
invalidates a persistent sequence preconditioner. The alias
`-p ca-mfamg` selects the same built-in MFAMG implementation and sets
`-mfamg_contact_aware 1`; `-p mfamg -mfamg_contact_aware 1` is equivalent.
Without the flag, the original MFAMG implementation remains unchanged.

Labels are equivalence classes, **not** penalty weights or numerical
thresholds. For a two-body mechanical system with an interface,
using one label for all DOFs of each body prevents aggregates that merge
DOFs across bodies. It does not prevent appropriate coupling through
the Galerkin operator. More detailed segmentation (e.g., separating
interior/interface or active/inactive contact regions) is possible,
but labels should be assigned consistently across DOFs that may
legitimately share algebraic aggregates.

When the active contact set or partition membership changes, call the
setter again before the next solve. The setter invalidates cached
sequence-preconditioner state even when the matrix object is unchanged.
A change in matrix coefficients is handled by MFAMG's existing
rebuild/refresh infrastructure; CA-MFAMG opts for a full rebuild instead
of the scalar numerical-only fast path to preserve the partition hierarchy.

## Constraints and limitations

* Existing MFAMG numerical restrictions still apply: double precision,
  `-scale none`, compatible near-nullspace modes and a CSR or capable
  `LIS_MATRIX_USER` operator; distributed MFAMG requires CSR.
* `-mfamg_block_size 1` is required in this first version. Contiguous
  block-Jacobi smoothing can cross partition boundaries, and is rejected.
  Weighted Jacobi and normal V-cycles are retained.
* The optional deep serial sibling merge is disabled **only in contact-aware
  mode**, because a merge of sibling aggregates could cross boundaries.
* CA-MFAMG allows a weak local-edge fallback **within the same partition**
  when strength-based pairing stalls. A bounded dense terminal of at most
  64 coarse DOFs is used for contact mode (the ordinary MFAMG default is 8).
  If many disconnected partitions prevent the hierarchy from reaching this
  size, setup reports a failure rather than allocating an unbounded dense LU.
* The partition is defined **per MPI rank's owned rows**; aggregates remain
  rank-owned as in MFAMG. No cross-rank contact-interface detection or
  cross-rank aggregation is introduced.
* This implementation is a **contact-boundary-aware hierarchy**, not a
  frictional/inequality contact active-set solver, block constraint
  preconditioner, or specialized contact Schur complement. It does not
  infer active contact nodes, normals, or penalty multipliers from matrix
  coefficients. A future block/interface smoother could improve performance.

## Regression

`test/testmfamgfineagg.c` also runs a deterministic 240-DOF SPD contact-
boundary proxy, with class boundaries that cut strong matrix edges.
It checks the alias/option profile, rejection of missing metadata,
solver ownership of copied labels, first coarse dimension, homogeneous
partition labels on **every coarsening level** for scalar and two-mode spaces, and the physical
preconditioner defect. The existing MFAMG tests continue unchanged.

These assertions establish boundary preservation and basic algebraic
correctness, **not** improved performance or convergence for 3D nonlinear
contact FEA. Benchmark against unmodified MFAMG and ILUTP on representative
real contact snapshots before selecting a production default.

Changing the public `LIS_SOLVER_STRUCT` size extends the ABI; applications
using the struct must be rebuilt against matching LIS headers/library.
