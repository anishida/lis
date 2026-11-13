# ILUT finite-factor and inverse-response regression

This test-only contribution extends the existing `test/testilupivot.c`. **No solver option, runtime guard, ILUT factorization algorithm, or default behavior is changed.**

## Why the extra check matters

The factor-statistics API already reports counts of stored L/U entries and the minimum/maximum absolute pivot magnitudes. Finite nonzero pivots alone do **not** prove that a triangular solve is numerically well-behaved. A complementary diagnostic on a known right-hand side is

```text
G(b) = ||A||_infinity * ||M^-1 b||_infinity / ||b||_infinity
```

The two deterministic four-by-four CSR regressions use the **existing** ILUT creation and `lis_psolve` path:

- A well-scaled positive diagonal with bounded expected inverse response.
- A diagonally ill-scaled positive matrix with all **finite** inverse pivots but a normalized single-vector response exceeding `1e15`.

They run without loading any private FEM data or requiring a network connection. The test checks finite pivot statistics, finite output, and the expected separation between a mild and an extreme response, and prints machine-readable `LIS_ILUT_NUMERICAL_PROBE` lines.

```sh
./test/testilupivot
```

The test is part of the existing LIS test executable and is exercised in CI. The cases are intentionally diagonal so they remain stable across serial/OpenMP ILUT factor partitioning.

**Interpretation:** `G(b)` is a one-RHS inverse amplification probe; it is *not* a condition-number estimator, a backward-error bound, a production-safe cutoff, or a guarantee about other right-hand sides. No fallback or diagonal shift is automatically applied. The much more complex RAS and Schur experiments remain outside this PR.
