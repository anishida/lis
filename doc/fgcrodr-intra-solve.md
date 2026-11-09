# FGCRO-DR: experimental within-solve harmonic Ritz recycling

This is an opt-in LIS FGCRO-DR feature, independent of MFAMG and offline near-nullspace vectors.

## Behavior

With the opt-in switch, FGCRO-DR activates the most recent harmonic
Ritz recycle subspace **after an Arnoldi restart within the same solve**.
It builds/updates the subspace from right-preconditioned search
directions, checks operator identity and structure/value generation,
then projects the physical residual using the saved relation A*Z=C
before proceeding with the current recycle dimension and new B stride.
True residual must still pass convergence criteria.

This is **not** an AMG mode generator and does not automatically
construct the offline six near-nullspace modes in MFAMG experiments.

## Opt-in configuration

Use normal LIS options, e.g.:

    -i fgcrodr -p jacobi -fgcrodr_k 1 -restart 4 -maxiter 300 -tol 1e-5

Per-solver LIS options are explicitly configured with `lis_solver_set_option`
(e.g. `-fgcrodr_intra 1 -fgcrodr_projected_relative 1`).
They are disabled by default; no process-global environment flags
control these paths.

| Solver option | Enabled value | Effect |
| --- | --- | --- |
| `-fgcrodr_intra` | 1 | Activate valid recycling on subsequent Arnoldi restarts |
| `-fgcrodr_intra_trace` | 1 | Print restart diagnostics, including projected residual norm before and after |
| `-fgcrodr_projected_relative` | 1 | Use relative projected-residual semantics; true residual remains authoritative |

All three accept only 0 or 1. They are stored in the individual solver's
private FGCRO-DR state. Under MPI, the algorithm checks that all ranks
agree on the mode switches before entering the Krylov core; mismatches
return a deterministic error rather than risking divergent collectives.
The settings must be configured coherently on all participating ranks.

The standalone test sets projected-relative normalization **both**
with intra-solve recycling ON and OFF, so it does **not** isolate that
normalization switch's impact. Convergence requires the true residual.
When tracing is enabled, the implementation checks that an orthogonal
projection does not materially increase the residual norm.

Recycled rank must remain smaller than restart size. Without the new
intra-solve switch, historical cross-solve FGCRO-DR remains the default.
For recycle dimension zero, the existing FGMRES path is used. Nonzero
recycle dimension is not implemented for the quad-precision FGCRO-DR
core.

## Independent public regressions

A new standalone test, test/testfgcrodrintra.c, runs automatically in
make check. It checks:

1. Restarted Arnoldi, harmonic Ritz build/update and actual reuse,
   using a deterministic 120-row nonsymmetric CSR system.
2. OFF and ON paths with identical operator, rhs, restart and tolerance.
3. True relative residual < 1e-5; no rejected directions; bounded
   orthogonality and A*Z-C errors below 1e-8.
4. A synthetic sequence A1 -> A2 with the *same* persistent solver but
   a distinct matrix and shifted diagonal, including the new operator
   generation in the saved recycling state.
5. The above synthetic operator sequence under two MPI ranks, with
   process-local row assembly.

PR #49 independently adds a MAXITER regression (test/testfgcrodrmaxiter.c)
for the active-recycle workspace cleanup path. Additional CI runs all
three targeted regression executables with AddressSanitizer,
LeakSanitizer and UndefinedBehaviorSanitizer enabled.

These checks cover specific source paths. They cannot guarantee that
every memory or numerical problem has been excluded.

## Non-goals and further validation

- No MFAMG v7/v8 source code, CalculiX integration or private snapshot
  file is included.
- The same-matrix value-update regression changes the diagonal,
  explicitly marks the value generation, reuses one solver and checks
  that the updated AZ=C relation and true residual are valid. The
  distinct-object A1 -> A2 case exercises a separate invalidation path.
- The manufactured A1 -> A2 sequence is **not** a CalculiX
  changing-contact Newton sequence.
- The feature is not a demonstrated
  general performance improvement; real application benchmarking,
  robustness review and further integration testing remain necessary.
- Full contact replay and offline mode research stay in
  [draft PR #47](https://github.com/3rav/lis/pull/47).

## B-workspace allocation bound

When intra-solve recycling is enabled, workspace B reserves capacity
for the largest projected recycle rank. The requested rank can exceed
restart-1, but no cycle can retain more than restart-1 directions.
The allocator therefore caps prospective rows to restart-1, while
still preserving enough space for an already active cross-solve rank.

The standalone test also uses a deliberately oversized
`-fgcrodr_k 1000000` with restart=2. It checks that the solver
still converges, maintains actual rank 1, uses recycling, and passes
the true residual / A*Z=C checks, without reserving a million rows.
This is a resource-safety regression, not a performance benchmark.

## Changing restart with a retained recycle basis

A caller may reduce the restart between two solves using the same
persistent solver. A recycle basis with dimension greater than or
equal to the new restart cannot be used as an Arnoldi augmentation:
at least one Krylov direction must remain available.

The FGCRO-DR entry path now drops such a stale-dimensional recycle
basis before it runs the first Arnoldi cycle and records a recycle
rejection. The next cycle may rebuild up to restart-1 directions.
The standalone test checks restart=4 with a rank-2 retained basis
followed by restart=2 with rank-1 rebuilding, true residual and
A*Z=C checks. This does not change the default behavior for valid
recycle dimensions.

The MPI regression intentionally applies conflicting per-rank intra-solve switches and verifies that the solve refuses to enter the numerical core rather than deadlocking.
