# Native Krylov solver extensions

The serial, real-double symmetric eigensolver supports optional controls
for active Krylov subspace capacity, reusable allocation, and additional
restart seeds. The original `lis_krylov_symmetric()` interface and its
existing public option, problem, and statistics structures remain unchanged.

These extensions are opt-in. They are not intended as default performance
optimizations, and their effectiveness depends on the operator, eigenvalue
distribution, and cost of shifted solves.

## Extended configuration

Include `lis_krylov.h`, initialize the new structures, and pass them to
`lis_krylov_symmetric_ex()`:

~~~c
lis_krylov_options options;
lis_krylov_stats stats;
lis_krylov_ext_options ext;
lis_krylov_ext_stats extra;

lis_krylov_options_default(&options);
lis_krylov_ext_options_default(&ext);
lis_krylov_ext_stats_default(&extra);

options.nev = 4;
options.ncv = 8;                    /* initial active dimension */
options.maxiter = 600;
options.tol = 1e-9;

ext.adaptive_ncv = 1;
ext.ncv_max = 24;                   /* allocated capacity */
ext.restart_directions = 3;

status = lis_krylov_symmetric_ex(
    &problem, &options, &ext, NULL,
    eigenvalues, eigenvectors, residuals, &stats, &extra);
~~~

Passing `NULL` for `ext` keeps the original fixed-capacity,
single-direction behavior. Passing `NULL` for the workspace uses
per-solve allocation. Passing `NULL` for `extra` suppresses extended
diagnostics. The original `lis_krylov_symmetric()` remains supported.

### Structure versions

`lis_krylov_ext_options` and `lis_krylov_ext_stats` each contain
`struct_size` and `version`. Initialize both through their respective
default functions before passing them to LIS.

Version 1 intentionally accepts only **exactly** the version-1 structure
size and version. A mismatched input *or output* structure returns
`LIS_KRYLOV_BAD_ARGUMENT`; no unknown trailing fields are read.
This is a conservative validation policy, not a promise that future
versions can be passed unchanged to older libraries. Supporting larger
structures or newer versions will require an explicit, reviewed contract.

Both typedefs also have named structure tags, allowing declarations
such as `struct lis_krylov_ext_options;` in client headers.

## Adaptive subspace sizing

When `ext.adaptive_ncv = 1`, `options.ncv` denotes the initial active
dimension. `ext.ncv_max` is the reserved maximum capacity. A zero
maximum selects `min(n, 2*ncv)`. After repeated restarts without
additional converged eigenpairs, the solver may grow the active
dimension while retaining its current basis.

**The maximum storage is reserved from the beginning**, even while
the active subspace is small. Adaptive sizing does not save the maximum
workspace memory compared with fixed `ncv=ncv_max`; a different
convergence trajectory may increase operator work. Equal-capacity
fixed-`ncv` runs are essential performance baselines.

### Spectrum-aware restart retention (experimental)

When `ext.adaptive_ncv = 1`, the solver can retain some additional Ritz
vectors across a full-subspace restart if the sorted Ritz values are
tightly clustered at the usual retention boundary. This is distinct
from increasing the active `ncv`: it avoids splitting a near-degenerate
group of Ritz values while reserving at least one fresh direction for
expansion. The number retained is bounded by the existing wanted count,
the active subspace dimension and a cap of roughly three quarters of
the reserved capacity.

For `ext.adaptive_ncv = 0` or the original
`lis_krylov_symmetric()` entry point, the pre-existing half-space
restart policy and scratch allocation are unchanged. Adaptive
workspaces reserve up to three quarters of the maximum `ncv` in
temporary restart storage rather than the historical half-space
scratch, so memory use **increases** when spectrum-aware retention is
enabled. The same buffer is reused across solves; no per-restart
allocation is introduced.

When combined with the separate `BE` (BothEnds) selector, adaptive
restarts use the **original fixed Ritz-retention boundary**. Interleaved
high/low Ritz ordering does not have adjacent gaps that correspond to
neighboring eigenvalues; applying a single-boundary gap heuristic to
that ordering would be mathematically misleading. `BE:ASC` and
`BE:DESC` use the same guarded fallback because suffixes change only
output order. A future two-flank spectral policy requires its own
numerical evidence and tests.

The heuristic is *not* guaranteed to improve every spectrum. Compare
to equal-capacity fixed-size solves and assess operator applications,
restarts, convergence and peak memory. The existing strict residual
acceptance rules remain unchanged. This is not an implicitly restarted
Lanczos implementation and does not introduce QR shifts.

## Reusable allocation

Use a workspace when multiple solves share a compatible solver shape:

~~~c
lis_krylov_workspace *workspace = NULL;

status = lis_krylov_workspace_create(&problem, &options, &ext,
                                     &workspace);
if (status == LIS_KRYLOV_SUCCESS) {
    /* Prepare the inverse operator associated with problem.sigma. */
    status = lis_krylov_symmetric_ex(
        &problem, &options, &ext, workspace,
        eigenvalues, eigenvectors, residuals, &stats, &extra);

    lis_krylov_workspace_destroy(workspace);
}
~~~

The workspace retains buffers, not eigenvectors, previous Krylov
iterations, or external matrix factorizations. Every solve resets
numerical state, locked-pair counts, and extended statistics.

The workspace shape is fixed at creation: dimension, transform, number
of requested eigenpairs, initial `ncv`, reserved maximum capacity,
and hard-locking setting must match on subsequent solves. The caller
supplies the extended configuration again on each solve so run-time
policies (including restart directions) may change while the reserved
shape remains compatible. An incompatible shape returns
`LIS_KRYLOV_BAD_ARGUMENT`. There is no implicit workspace resizing.

Applications are responsible for ensuring that the supplied
shift-invert callback corresponds to **each current** `problem.sigma`.
The workspace cannot verify a factorization managed by external code.
Do not use the same workspace concurrently in multiple threads.
`lis_krylov_workspace_destroy(NULL)` is safe. A failed solve does
not invalidate the allocated workspace, provided the application
restores valid callbacks and configuration before reusing it.

## Additional restart directions

`ext.restart_directions = 1` selects the existing scalar restart
behavior. Larger values attempt more deterministic, metric-orthogonal
search directions after restarts. Dependent or invalid candidates can
be rejected, especially with a semidefinite metric, and the actual
number of accepted directions may be lower than requested.

This is **restart enrichment**, not a full block Lanczos/Arnoldi
method: all matvec and inverse callbacks remain scalar. More
directions can increase operator work without improving convergence.
The feature requires independent numerical evidence before it should
be recommended for production workloads.

## Extended diagnostics

The output `lis_krylov_ext_stats` includes:

| Field | Meaning |
|---|---|
| `ncv_initial` | Initial active subspace size |
| `ncv_final` | Final active subspace size |
| `ncv_growths` | Number of growth events |
| `basis_vectors` | Accepted basis-vector count |
| `projected_solves` | Dense projected eigensolver invocations |
| `extra_seed_attempts` | Additional restart candidates tried |
| `extra_seed_accepted` | Additional candidates accepted |

These counts complement, rather than replace, the original
`lis_krylov_stats` fields.

## Building and testing

The extension tests use only the scalar operator callbacks and the
existing C99 solver implementation. For example:

~~~sh
cc -O2 -std=c99 -Iinclude src/esolver/lis_krylov.c \
   test/test_krylov_unified.c -lm -o test_krylov_unified
./test_krylov_unified

cc -O2 -std=c99 -Iinclude src/esolver/lis_krylov.c \
   test/test_krylov_adaptive_ex.c -lm -o test_krylov_adaptive_ex
./test_krylov_adaptive_ex

cc -O2 -std=c99 -Iinclude src/esolver/lis_krylov.c \
   test/test_krylov_api_contract.c -lm -o test_krylov_api_contract
./test_krylov_api_contract
~~~

Build and run the final test with `-m32` on a toolchain providing
32-bit development libraries to exercise the size-overflow rejection
path. On 64-bit targets it still checks version and structure-size
errors, but skips that platform-specific overflow case.

The two standalone benchmark drivers compare adaptive sizing with
fixed initial, intermediate, and maximum capacities and compare
repeated one-shot solves with workspace reuse. These are synthetic
diagnostics, not representative end-to-end sparse-eigensolver timings.
