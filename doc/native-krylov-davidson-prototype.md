# Experimental diagonal Davidson eigensolver (standard symmetric)

This is an **independent, self-contained C99 proof of concept**, inspired by the
Davidson-family algorithm available in Spectra. No Spectra source code is copied,
and **no public LIS interfaces or production solver implementation are changed**.

The prototype handles standard real symmetric problems `A x = lambda x` with
the largest algebraic eigenvalues requested. Fixed `DV_NMAX=72` and
`DV_MMAX=28` exist only in the test harness.

## Algorithm

The method starts with several coordinate vectors chosen from the largest
entries of `diag(A)`, orthogonalizes them twice (modified Gram-Schmidt), and
calculates both the basis `V` and its images `A V`.

Each outer iteration solves a small symmetric Rayleigh–Ritz problem:

```text
T = V^T A V
T y = theta y
x = V y
r = A x - theta x
```

The original-problem physical residual, not just a projected indicator,
determines convergence.

For the first unconverged wanted Ritz pair, a **one-step diagonally
preconditioned approximation to the Jacobi–Davidson correction** is
constructed. Writing `M = regularized(diag(A) - theta I)`, the direction is:

```text
u       = M^(-1) r
w       = M^(-1) x
epsilon = (x^T u) / (x^T w)
t       = -u + epsilon*w
```

When `x^T w` can be evaluated reliably, this enforces `x^T t = 0`
(up to roundoff). Crucially, the inverse-preconditioned `M^(-1)x`
direction is **not** generally parallel to `x`, so this is not the
redundant Euclidean projection `t <- t - x*(x^T t)`. A deterministic
unit test checks both the formula and that the resulting correction
differs from the simply projected Davidson direction.

The prototype still **does not solve the full Jacobi–Davidson projected
correction equation**. It makes one inexpensive diagonal approximate step.
Since `M` can be indefinite, `x^T M^(-1)x` may vanish by cancellation.
Such cases fall back to ordinary regularized Davidson with removal of the
component parallel to `x`; a dedicated test exercises an exactly
cancelling diagonal. Tiny diagonal denominators are regularized, and
nonfinite or linearly dependent directions trigger the existing fallback
to the physical residual or a deterministic independent coordinate vector.

When the basis reaches capacity, a **thick restart** retains the leading Ritz
vectors and their cached `A V` images. It does not reapply `A` to rotated
vectors. Neither preconditioners nor matrix factorizations are exposed via
new public callbacks.

## Tests and accuracy

`test/test_krylov_davidson_prototype.c` exercises **36 deterministic cases**:
three symmetric spectral shapes (well separated, moderately clustered and
tightly clustered), four maximum subspace sizes (12, 16, 20, 24), and
three different initial basis sizes (3, 4, 5). Checks cover:

- Accurate top three eigenvalues relative to independently established
  references, together with physical residuals below `2e-8`.
- Orthonormal converged eigenvectors.
- Finite regularized corrections when `diag(A)_i == theta`.
- Verified oblique diagonal correction formula and nonredundant search
  direction relative to an ordinary projected Davidson step.
- Safe fall back when `x^T M^(-1)x` cancels for an indefinite diagonal.
- Compiler and sanitizer validation.

Two additional **equal-diagonal tridiagonal Toeplitz** cases (n=32 and
n=40) deliberately make the diagonal preconditioner weak. Their three
largest eigenvalues have the exact analytic form
`1 + 2*b*cos(j*pi/(n+1))`. These tests independently recompute physical
residuals directly from matrix coefficients and require at least one
full-subspace restart across the two cases. This exercises the thick
restart and guards against a test suite that converges only because the
initial coordinate vectors are near exact eigenvectors.

`test/bench_krylov_davidson_compare.c` additionally constructs the **same
three matrices** for the native symmetric Krylov solver and this prototype,
then independently recomputes `A*x - lambda*x` from the dense matrix before
accepting either solver's output.

Original experimental `apply_a` counts before this correction change (n=64, nev=3, ncv=24):

| Spectrum | Native Krylov | Diagonal Davidson |
| --- | ---: | ---: |
| Separated | 42 | 8 |
| Clustered | 54 | 12 |
| Tightly clustered | 43 | 13 |

These numbers **are not production performance results**. Davidson starts
from multiple coordinate directions informed by `diag(A)`, while native
Krylov starts from one supplied direction. The results do not separate the
benefit of diagonal-informed initialization from that of the Davidson
correction. Counting `A` calls also omits the cost of obtaining a diagonal,
small-matrix eigendecompositions, orthogonalization and memory traffic.

## Running standalone

```sh
cc -std=c99 -O2 test/test_krylov_davidson_prototype.c \
  -lm -o /tmp/test_davidson
/tmp/test_davidson

cc -std=c99 -O2 -Iinclude src/esolver/lis_krylov.c \
  test/bench_krylov_davidson_compare.c -lm -o /tmp/davidson_compare
/tmp/davidson_compare
```

The permanent numerical test is also registered with `make check`. The
tracked `test/Makefile.in` is synchronized so full builds do not require
regenerating Autotools outputs.

## Not production-ready

- No generalized metric B, shift-invert, inexact inner solves, semidefinite
  matrices, interval extraction, locking or production-safe restart policy.
- Requires access to `diag(A)`, which the existing LIS callback API does
  not currently expose. Acquiring this cheaply or receiving a caller-owned
  diagonal would require an explicit, versioned extension contract.
- A diagonal may be a poor preconditioner; small denominators require
  regularization, and eigenvalue multiplicity or a strong off-diagonal
  coupling can cause convergence failures.
- Hard-coded demonstration dimensions and allocation; no production
  error/status, cancellation, MPI or scalable sparse-matrix support.
- No real cost-per-solve benchmark or comparison with equal starting
  subspaces. Additional tests on sparse and clustered problems are needed.

The method now reports how often the oblique diagonal correction was used
and how often it fell back to ordinary Davidson. Any claimed improvement
must compare the old and new correction on the **same** matrices, seeds,
tolerances, and subspace sizes, including total operator calls and
restarts; fewer calls alone do not establish end-to-end speedups.

The next engineering milestone would be a **separate, versioned opt-in
preconditioner interface**, with explicit memory ownership, callback-error
semantics, a well-defined generalized correction equation, and a controlled
benchmark where the initial search space and all costs are comparable.
