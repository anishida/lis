# Stage 6: experimental matrix-free Cayley transform (native LIS Krylov)

This is a **test-only C99 proof of concept** for real symmetric generalized
eigenproblems A x = lambda B x, where B is symmetric positive definite.
It does **not** change the production Krylov solver, public API/ABI, default
behavior, dependencies, or MPI support.

## Operator and spectral mapping

For real shifts sigma != tau and nonsingular A - sigma B:

    C = (A - sigma B)^(-1) (A - tau B)
      = I + (sigma - tau) (A - sigma B)^(-1) B.

Original eigenvalue lambda maps to

    mu = (lambda - tau) / (lambda - sigma)
    lambda = sigma + (sigma - tau) / (mu - 1).

This is **not** the orthogonal Cayley transform of a skew-symmetric
generator. For symmetric A and SPD B, C is self-adjoint in the B metric.
It is generally **not** Euclidean symmetric.

Our fixture uses B = D^2 with positive diagonal D, so the Euclidean
symmetric operator accepted by LIS_KRYLOV_STANDARD is

    C_y = D C D^(-1).

It is applied matrix-free through the application's apply_a, apply_b
and prepared solve(A - sigma B) callbacks. The generalized eigenvectors
are recovered as x = D^(-1) y, and original eigenvalues via the map above.

## Regression and reference

The test constructs A = D T D, B = D^2, with T symmetric tridiagonal
Toeplitz: diagonal 4, offdiagonal 0.25. The exact spectrum is

    lambda_j = 4 + 0.5 cos(j pi / (n + 1)).

Two shift pairs target opposite spectral ends:
(sigma, tau) = (4.65, 3.2) and (3.35, 4.8).

The permanent test/test_krylov_cayley_prototype.c checks:

- Analytic eigenvectors, eigenvalue map and affine shift-invert identity.
- Symmetry of the explicitly transformed Euclidean operator.
- Rejection of equal/nonfinite shifts, an exactly singular shift and
  propagation of callback errors.
- Native Krylov convergence on the actual matrix-free Cayley action.
- Independent original-pencil residuals, analytic eigenvalues, and
  B-metric eigenvector orthogonality.
- Native SHIFT_INVERT_B results for the same problem and factorization.

## Running

Operator-only checks require just a C99 compiler and libm:

    cc -std=c99 -O2 -DCAYLEY_ACTION_ONLY test/test_krylov_cayley_prototype.c -lm -o /tmp/cayley_action
    /tmp/cayley_action

Full checks link the existing native Krylov core:

    cc -std=c99 -O2 -Iinclude src/esolver/lis_krylov.c test/test_krylov_cayley_prototype.c -lm -o /tmp/cayley_krylov
    /tmp/cayley_krylov

The full test is registered in make check, with the tracked Makefile.in
updated along with Makefile.am.

## Scope and limitations

This is **not** a production Cayley solver mode. The similarity
B^(1/2) C B^(-1/2) is inexpensive in this fixture only because B is
positive diagonal. For arbitrary SPD B, a square-root application is
not supplied by the existing LIS callback contract. Passing the raw
unsymmetrized C into STANDARD mode would be incorrect.

For the SPD generalized problem, Cayley is an *affine* transform of
SHIFT_INVERT_B: C = I + (sigma - tau) OP_SI. Consequently it uses the
same invariant subspaces and does not automatically accelerate convergence.
A production option would first need a reason to change Ritz selection,
a versioned opt-in contract, behavior near poles, verified physical
residuals, and inexact shifted-solve handling.

The prototype's unpivoted tridiagonal factorization is valid only for
the tested shifts outside the spectrum. It is not a general sparse
factorization, nor is a callback count a runtime benchmark.
