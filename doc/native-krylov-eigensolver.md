# Native restarted symmetric subspace eigensolver

This module implements an iterative, thick-restarted Rayleigh–Ritz subspace
method directly in LIS. The implementation is self-contained; its only
runtime dependencies are the standard C library and libm. It supports
symmetric standard and generalized eigenproblems through arbitrary storage
formats or user-defined matrix-vector operators.

## Two interfaces

The standard LIS interface can be selected using `-e krylov` in `lis_esolve`
(standard) or `lis_gesolve` (generalized). It accepts `LIS_MATRIX` and
`LIS_VECTOR`, including `LIS_MATRIX_USER` for serial matrices. For example:

~~~c
lis_esolver_create(&esolver);
lis_esolver_set_option(
    "-e krylov -ss 4 -nev 4 -kwhich LM -etol 1e-9 -emaxiter 500",
    esolver);
status = lis_gesolve(A, B, x, &lambda0, esolver);
~~~

- `-ss` is the maximum number of stored results; `-nev` (if supplied) is
  the number requested (at most `-ss`).
- `-kncv` sets the internal Krylov subspace capacity; zero selects
  `min(n, 3*nev+8)`.
- `-kwhich` accepts `LM`, `SM`, `LA`, `SA`; ordering refers to the
  **transformed** Ritz spectrum.
- `-kresidual_scale` supplies an application-based estimate of `||A||`.
  Zero selects an automatically sampled estimate.
- `-kinner_tol` (zero by default) enables explicit verification of each
  shifted inverse with two additional matrix-vector operations.
- `-ksemidefinite 1` opts in to a positive-semidefinite B metric in the
  generalized shift-invert mode. Starting and replenishment directions
  are mapped through `(A-sigma B)^-1 B` before orthogonalization.
- `-etransform shift_invert` uses the positive-definite B metric by default, and
  `-etransform buckling` uses the positive-definite A metric, with nonzero
  shift required by the latter. If B is present and transform is AUTO,
  shift-invert with the B metric is used. Standard mode uses A directly.
- Existing LIS linear-solver options configure the internal shifted solve;
  preconditioner setup is reused across Krylov iterations.

The LIS adapter explicitly rejects MPI-distributed matrices, complex or
non-double scalar types, and unsupported precision settings. A build with
MPI enabled may use it with one rank; it is **not** a distributed eigensolver.
Existing Fortran applications can select the method through LIS's standard
options interface; no new dedicated Fortran API is provided.

The independent C callback API is declared in `lis_krylov.h`:

~~~c
lis_krylov_problem problem = {
    .n = n, .transform = LIS_KRYLOV_SHIFT_INVERT_B,
    .sigma = sigma, .ctx = ctx,
    .apply_a = apply_a, .apply_b = apply_b,
    .solve = solve_prepared_shifted_system
};
lis_krylov_options options;
lis_krylov_stats stats;
lis_krylov_options_default(&options);
options.nev = nev;
status = lis_krylov_symmetric(&problem, &options,
                               eigenvalues, eigenvectors, residuals, &stats);
~~~

Callbacks receive contiguous **serial**, real-double arrays. `solve`
implements a fixed inverse action of `(A - sigma B)` and must reuse any
prepared solver/factorization state; the eigensolver does not assemble
or factorize matrices. It requires sufficiently accurate, nearly symmetric
inverse actions. For iterative solves, set the inner tolerance well below
`options.tol` (typically at least one order of magnitude tighter) and
consider enabling `options.inner_tol` to detect inaccurate solutions.

## Transformations and metric normalization

| Transform | Inner product | Operator | Original eigenvalue |
|---|---|---|---|
| `STANDARD` | Euclidean | `A` | `lambda = mu` |
| `SHIFT_INVERT_B` | `B` (positive definite, or PSD when enabled) | `(A-sigma B)^-1 B` | `lambda = sigma + 1/mu` |
| `SHIFT_INVERT_A` | `A` (positive definite) | `(A-sigma B)^-1 A` | `lambda = sigma*mu/(mu-1)` |

For PSD B, the operator `(A-sigma B)` must remain nonsingular, and
only finite eigenpairs with nonzero B-norm are sought. The solver rejects
null B-metric directions and may return `NO_CONVERGENCE` with partial
results if fewer finite eigenpairs exist than requested. This mode does
not implement a general singular-pencil eigensolver.

Vectors are normalized in the **chosen metric** (I, B or A). Generalized
operators must satisfy the corresponding symmetry and definiteness
requirements. In the last mode, B may be indefinite and sigma must be
nonzero. The result order follows `-kwhich` on the transformed eigenvalues.

## Convergence and partial results

The solver reports a backward-error indicator per returned pair:

`||A*x - lambda*B*x|| / max(||A*x|| + |lambda|*||B*x||,
                            residual_scale*||x||)`.

The scale-free fallback is a sampled estimate of `||A||`; for near-zero
modes and highly nonuniform matrix scaling, providing an externally
computed operator norm estimate is recommended. The residual is a
numerical accuracy indicator, **not proof** that all requested eigenvalues
were found or that multiplicities are complete.

At `LIS_KRYLOV_NO_CONVERGENCE`, the leading `stats.nreturned` approximate
pairs remain in `which` order. `stats.nconv` counts pairs with residual at
most the requested tolerance, but converged pairs are **not compacted** to
the front. Each entry of `residuals[]` must be inspected independently.
`stats.max_residual` is the largest residual among *all* returned pairs.
Errors other than nonconvergence invalidate outputs (`nreturned=0`).
Starting vectors can be supplied via `options.initial_vector` or chosen
with `options.random_seed` for reproducible restarts.

The implementation uses adaptive DGKS metric reorthogonalization, reuses
cached metric/operator products, updates the projected matrix incrementally,
retains Ritz vectors on restart, injects reproducible independent directions
to resolve repeated eigenvalues, and skips expansion of already converged
Ritz vectors. Optional persistent hard locking (`-khardlock 1`, or
`options.hard_locking=1`) freezes fully verified eigenpairs in a separate
metric-orthogonal subspace and purges the active Ritz space after locking.
Hard locking is off by default because it can require extra solves and
additional storage proportional to `n*nev`. No global process state is used.

## Build and validation

~~~sh
autoreconf -fiv
mkdir build && cd build
../configure
make -j2
make check
~~~

The test suite includes analytic standard and generalized problems, shifted
inverse transforms, indefinite B, 1e9 scaling, zero and near-zero
values, selection `SM`/`LA`/`SA`, repeated values, non-diagonal sparse
problems, partial convergence, callback failures and inexact solves.
The PSD regression covers a coupled matrix with massless degrees of freedom,
a start in `ker(B)`, a scaled pencil, and a rank-deficient B metric.
The sparse regression uses an analytic relative eigenvalue tolerance.
The stress regression checks clustered and fourfold repeated generalized
values, B-orthogonal invariant subspaces, and repeated restarts with
small search subspaces. Additional large sparse performance evaluation
is recorded separately from the implementation.

## Inertia-based interval counting

The native C callback API exposes `lis_krylov_count_interval` as an
**optional** interval-counting helper. For a symmetric pencil `A-lambda B`
with **positive-definite B**, a caller-provided inertia callback reports the
negative and zero pivot counts of `A-shift*B`, typically using a stable
symmetric indefinite factorization. The count in the open interval
`(lower, upper)` equals the difference in negative pivot counts. Endpoints
must be nonsingular; a reported zero pivot returns
`LIS_KRYLOV_SINGULAR_SHIFT`. Callback failures and inconsistent counts
have explicit error codes.

The eigensolver neither constructs the factorization nor uses the callback
to enumerate the eigenvectors in the interval. Inertia counting does not
by itself prove that Ritz extraction is complete. General semidefinite B
and complex/Hermitian pencils are outside this helper's contract.

## Limitations

- Native C99 real-double implementation, serial/local operators only.
- Complete eigenvalue coverage and multiplicities cannot be inferred solely
  from small residuals; use multiple starting vectors or external checks.
- No distributed MPI, complex/quad arithmetic, or dedicated new Fortran
  ABI is implemented in this version.
- No performance claim for large sparse matrices without timed benchmarks
  that include the cost of the application's inverse operation.
