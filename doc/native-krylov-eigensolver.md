# Native restarted symmetric subspace eigensolver

This module implements an iterative, thick-restarted Rayleigh–Ritz subspace
method directly in LIS. The implementation is self-contained; its only
runtime dependencies are the standard C library and libm. It supports
symmetric standard and generalized eigenproblems through arbitrary storage
formats or user-defined matrix-vector operators.

## Why an additional eigensolver?

LIS already provides subspace iteration, Lanczos, Arnoldi, and FEAST
eigensolvers. This module is an **additional explicit option**, not a
replacement for any of them. Its intended contribution is a standalone,
callback-driven symmetric Rayleigh–Ritz workflow with user-managed
shifted inverse actions, thick restart, explicit original-problem
residuals, optional positive-semidefinite metric handling, optional
persistent locking, and a separate inertia-based interval helper.
These features are grouped behind a new algorithm selector rather than
silently changing the convergence and storage semantics of existing
eigensolvers.

No claim of superior speed or robustness relative to the existing LIS
methods is made here. A controlled benchmark against those methods,
including the application's shifted linear-solve/factorization costs,
is still required to establish a performance advantage.

## Scope of the serial implementation

This is a deliberately bounded first implementation: dense Rayleigh–Ritz
projections and callbacks operate on complete local vectors, and the
adapter explicitly rejects distributed inputs instead of silently
performing an incorrect rank-local eigensolve. An MPI implementation
would require coordinated global inner products, restart operations,
and distributed result ownership; a one-rank execution inside an
MPI-enabled build is not distributed support.

The standalone C callback interface uses `int` indices and `double`
scalars intentionally. The LIS adapter checks dimensions before
conversion from `LIS_INT` and rejects unsupported sizes; successful
`--enable-longlong` builds do **not** mean this interface can address
arbitrarily large 64-bit dimensions. The restarted method requires
`nev < n` because `ncv > nev`; requesting the full spectrum should
use another eigensolver.

## OpenMP shared-memory execution

When LIS is configured with `--enable-omp`, the native Krylov core
parallelizes long compensated dot products and stable norms with fixed
independent chunks, as well as vector scaling, metric projections,
Ritz vector construction, and restart basis copies. Small vectors
remain on the low-overhead serial path; builds without OpenMP retain
the standard C99 implementation. The chunk layout and summation order
are fixed independently of OpenMP scheduling; tiny roundoff differences
between one and multiple threads can still occur.

The number of threads is selected by the OpenMP runtime, for example
`OMP_NUM_THREADS=4`. Matrix-vector products and shifted linear solves
remain delegated to LIS or application callbacks, whose threading
behavior is independent of the Krylov vector kernels. Dense projected
eigendecomposition is intentionally serial at the modest subspace sizes
used here. No MPI-distributed execution is provided.

The `test_krylov_openmp` regression runs a 65,536-dimensional
generalized problem with a nontrivial restart, both with and without
persistent hard locking. OpenMP builds exercise 1, 2, 4 and 8 threads
and check the analytic eigenvalues, original-problem residuals,
B-metric orthogonality, and agreement with the single-thread solution.
This validates correctness, not a guaranteed parallel speedup.

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
- `-kwhich` accepts `LM`, `SM`, `LA`, `SA`, and `BE`. All modes
  select Ritz values from the **transformed** operator spectrum. `BE`
  alternates the largest and smallest algebraic Ritz values, starting
  from the largest, so odd `nev` favors the upper end.
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

An independent validation workflow was run in the contributor fork
rather than added to this feature diff. The verified run
https://github.com/3rav/lis/actions/runs/41603936979 reports 89/89
`make check` tests and 8/8 focused Krylov tests in GCC long-long,
GCC long-long with OpenMP, and Clang configurations, as well as
standalone numerical AddressSanitizer and UndefinedBehaviorSanitizer
checks. These are **not** upstream PR Checks, nor a substitute for the
maintainer's platform and configuration testing.

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

The eigensolver does not construct the factorization. An optional
`lis_krylov_solve_interval` wrapper combines the independently counted
interval with a caller-prepared midpoint shift-invert operator and requests
exactly that many eigenpairs. SUCCESS requires all expected eigenpairs to
converge with true residuals and lie inside the open interval; otherwise
NO_CONVERGENCE or a more specific error is returned. The midpoint must not
be an eigenvalue (the shifted operator must be invertible), and the output
capacity and search subspace must be sufficient. Boundary and midpoint
singularities require a different shift or interval decomposition.

Inertia supplies an independent **exact count only under the assumed
symmetric pencil, SPD B, and reliable inertia factorization**. The wrapper
checks that the number of converged, in-window approximate pairs matches
that count. These are useful numerical consistency checks; at finite
precision, the reported residuals and eigenvalues are approximations,
and their agreement with the inertia count is **not a rigorous interval
certificate**, especially for modes near interval boundaries or clusters
narrower than the requested residual tolerance. More stringent guarantees
would require appropriate eigenvalue-enclosure/error bounds in addition
to inertia. General semidefinite B and complex/Hermitian pencils are
outside this helper's contract.

## Limitations

- Native C99 real-double implementation, serial/local operators only.
- Complete eigenvalue coverage and multiplicities cannot be inferred solely
  from small residuals; use multiple starting vectors or external checks.
- No distributed MPI, complex/quad arithmetic, or dedicated new Fortran
  ABI is implemented in this version.
- No performance claim for large sparse matrices without timed benchmarks
  that include the cost of the application's inverse operation.

## Thresholded Jacobi refinement of the projected problem

The native thick-restarted Krylov solver diagonalizes a relatively small
symmetric projected matrix during each Rayleigh–Ritz step. This implementation
uses an internal cyclic Jacobi eigendecomposition without requiring LAPACK.

For the first three Jacobi sweeps, the projection solver skips
off-diagonal pivots no larger than 5% of the **largest off-diagonal
entry measured at the beginning of that sweep**. The 5% threshold is
a heuristic, not an error bound or a universally optimal tuning value.
Later sweeps revert to the original machine-epsilon pivot cutoff, and
the stopping test remains unchanged. However, a changed rotation sequence
can affect roundoff, eigenvector bases in degenerate eigenspaces, and
the number of Krylov iterations. Public ABI, callbacks, and requested
eigenvalue ordering are unchanged.

A standalone regression (`test_krylov_projected`) exercises a symmetric,
orthogonally transformed matrix with known clustered and repeated
eigenvalues, verifying the returned eigenpairs, original-problem residuals
and orthogonality. It also exercises hard locking.

A small standalone timing/accuracy experiment is available without
installing LIS:

~~~sh
cc -O3 -std=c99 -Iinclude test/bench_krylov_projected.c \
   -lm -o bench_krylov_projected
./bench_krylov_projected
~~~

The benchmark directly exercises the internal Jacobi kernel at four
projected dimensions with **different repeat counts** (24: 1600,
48: 480, 72: 120, 96: 40). It performs an untimed warm-up, seven timed
trials per dimension and reports the median, minimum and maximum process
CPU time as well as the **median microseconds per eigendecomposition**.
Compare variants at the *same dimension*, not total times between
different dimensions. Each trial checks the reconstruction residual and
eigenvector orthogonality of the resulting matrix decomposition.

Kernel-only results are not end-to-end Krylov timings. A baseline and an
optimized build should be compared on the same hardware, preferably with
independent runs, and uncertainty reported. The new threshold can change
the number of sweeps; any speedup is workload-dependent. No new
dependencies, user-facing solver options or CI workflow files are
introduced.

## BothEnds selection and independent output sorting

The native solver now accepts `options.which = "BE"`. It interleaves
**high/low extrema of the transformed Ritz spectrum**, starting at the
high end; for `nev = 5`, the target ordering is
`largest, smallest, 2nd-largest, 2nd-smallest, 3rd-largest`. Ties
are handled deterministically, without a special ordering guarantee
among equal Ritz values. This
behavior is deliberately **not** defined as both ends of the original
generalized pencil: after shift-invert or buckling transformation the
corresponding original eigenvalues may lie near the shift instead.

An independent `":ASC"` or `":DESC"` suffix can be appended to **any**
of the five native `which` selectors, for example:

~~~c
options.which = "BE:ASC";  /* choose transformed ends, report ascending lambda */
options.which = "LM:DESC"; /* choose largest |mu|, report descending lambda */
~~~

The suffix affects only the **published original eigenvalues** (and
moves their residuals and eigenvectors together). It does not modify
Ritz candidate selection, restarts, convergence tolerances or which
values are requested. It also sorts available partial results when
the solver returns `LIS_KRYLOV_NO_CONVERGENCE`. Omitting the suffix
preserves the original output ordering for `LM`, `SM`, `LA` and
`SA`. BothEnds interleaving is applied to returned results even
after hard locking.

`-kwhich BE` is supported through the existing LIS eigensolver
command-line options without changing the existing three-character
`krylov_which` field. **Output sort suffixes are native-API only**
(`lis_krylov_symmetric()` or `lis_krylov_symmetric_ex()`). The
LIS `-kwhich` option remains two characters due to the public
eigensolver structure's fixed binary layout; `-kwhich BE:ASC` is
therefore not a valid CLI option. No public structure layout changes
or added dependencies are introduced.

The regression `test_krylov_both_ends` validates correct selected and
published eigenvalues, physical residuals, metric orthogonality,
standard / B-metric shift-invert / A-metric shift-invert modes,
hard locking, both native entry points and reusable workspaces. The
LIS integration regression also tests `-kwhich BE`.

## Optional extended Krylov API

[The native Krylov extension guide](native-krylov-extensions.md)
documents adaptive active subspace sizing, reusable allocation,
optional multi-direction restart seeding, validation rules, and
extended diagnostics.

The original `lis_krylov_symmetric()` and existing public structures
remain unchanged. Optional controls use versioned extension structures
and a single `lis_krylov_symmetric_ex()` entry point. The workspace
retains allocation only, not external shifted factorizations. The
caller must prepare the inverse callback for each spectral shift.
Extended restart seeding is not a full block Krylov algorithm.
