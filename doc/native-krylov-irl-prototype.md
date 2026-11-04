# Experimental implicitly restarted Lanczos (IRL) reference

This directory contains an **independent C99 algorithm prototype**, not a replacement for
the production native Krylov eigensolver. It is inspired by the mathematical approach
of Spectra's symmetric implicitly restarted Lanczos routine; no Spectra source code
was copied.

The experiment covers **real, standard symmetric eigenproblems** (B = I), selects the
largest algebraic eigenvalues, and makes no modifications to LIS's public headers,
callbacks, ABI, workspace layout, or existing solver behavior.

## Mathematical kernel

The symmetric Lanczos relation after m steps is

```text
A V_m = V_m T_m + beta_m v_(m+1) e_m^T
```

where T_m is symmetric tridiagonal. After computing its Ritz values, the smallest
m-k values become *implicit* QR shifts. Each shifted QR step applies an orthogonal
similarity to the projected matrix, and the products of the small Q factors are
accumulated into Q_total:

```text
T_m - mu I = Q R
T_m' = Q^T T_m Q
V_m' = V_m Q_total
```

After m-k shifts, the old residual couples only to the last retained direction.
The next Lanczos seed is assembled as

```text
f_k = T'_(k+1,k) * V'_(:,k+1)
      + beta_m * (e_m^T Q_total e_k) * v_(m+1)
```

and normalized before continuing Lanczos from step k. The kernel validates that
the unwanted couplings are numerically negligible, that the retained basis remains
orthonormal, and that the transformed projected matrix remains tridiagonal.

For simplicity and auditability, the prototype explicitly constructs dense Q
factors using Givens QR. **It is not a performant bulge-chasing implementation**:
QR work scales far worse than the O(m^2) or better small-matrix machinery one
would use in a production IRL implementation. The experiment is therefore a
check of mathematical feasibility and numerical behavior, not a shipping algorithm.

## Convergence and testing

The implementation uses two-pass orthogonalization and a Lanczos Ritz-residual
estimate `abs(beta_m * y_last)` only to decide whether to request a physical
verification. **It never accepts a Ritz value without evaluating the actual
original-problem residual.**

The permanent test `test_krylov_irl_prototype.c` exercises **36 combinations**:
three deterministic matrix spectra (well-separated, moderately clustered,
tightly clustered), four projected dimensions (16, 20, 24, 28), and three retained
dimensions (5, 8, 11). It checks known reference eigenvalues, physical residuals,
orthogonality, and the QR restart structure, with GCC, Clang, OpenMP and
AddressSanitizer/UndefinedBehaviorSanitizer validation.

`bench_krylov_irl_compare.c` compares the standalone prototype with the original
native eigensolver on the **same three matrices**, validating returned physical
eigenpairs and reporting callback counts. Preliminary counts for n=64, nev=3,
m=24, k=9 were:

| Spectrum | Native `apply_a` | Experimental IRL `apply_a` |
| --- | ---: | ---: |
| separated | 235 | 74 |
| clustered | 224 | 57 |
| tightly clustered | 170 | 42 |

These counts are only **exploratory**: the solvers have different initial vectors,
restart strategies and bookkeeping, so this is not an apples-to-apples runtime
comparison. The prototype has costly dense QR transforms and may run **slower**
despite fewer A applications.

## Build standalone

From the project root:

```sh
cc -std=c99 -O2 test/test_krylov_irl_prototype.c -lm -o /tmp/irl_prototype
/tmp/irl_prototype

cc -std=c99 -O2 -Iinclude src/esolver/lis_krylov.c \
  test/bench_krylov_irl_compare.c -lm -o /tmp/irl_compare
/tmp/irl_compare
```

The standalone test is also registered with `make check`; the generated
`test/Makefile.in` is committed for projects that do not regenerate Autotools.

## Not supported / reasons not to merge as a production solver yet

- No generalized metrics, semidefinite B, shift-invert or inexact inverse solves.
- No automatic sizing, hard locking, reuse, multi-direction restarts or interval solver.
- Fixed demonstration bounds (`MAXN=80`, `MAXM=32`), not scalable storage.
- Breakdown currently aborts that prototype trial rather than performing the
  full existing solver's careful restart/error-recovery behavior.
- The dense QR implementation should be replaced with tridiagonal implicit
  bulge chasing, followed by tests for repeated/near-multiple eigenvalues,
  loss of orthogonality, end-to-end runtime and large sparse operators.
- Existing public ABI and numerical behavior remain entirely unchanged.

The next engineering milestone is to separate the small projected-QR kernel
from the fixed test harness, implement properly bounded scratch/workspace,
and compare an internal opt-in path against the existing solver with equal
starting vectors and identical callback/termination contracts.
