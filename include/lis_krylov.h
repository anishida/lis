/* Native restarted symmetric subspace eigensolver.
 * Copyright (C) 2026, contributors to LIS.
 * Redistribution and use are governed by the LIS project BSD license.
 */
#ifndef LIS_KRYLOV_H
#define LIS_KRYLOV_H

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*lis_krylov_matvec)(void *ctx, int n,
                                 const double *x, double *y);

/* Symmetric standard / generalized eigenproblem A x = lambda B x.
 * All vectors are real double, local (NOT MPI distributed), length n.
 * STANDARD: B=I, Euclidean inner product, OP=A.
 * SHIFT_INVERT_B: B positive definite by default; optional PSD B.
 *                 OP=(A-sigma B)^-1 B; B metric.
 * SHIFT_INVERT_A: A positive definite; OP=(A-sigma B)^-1 A; A metric.
 * In both shifted modes, 'solve' is an application-provided inverse
 * operation on the prepared, fixed shifted matrix. Callbacks return 0
 * on success. Input and output buffers do not alias.
 */
enum {
    LIS_KRYLOV_STANDARD = 0,
    LIS_KRYLOV_SHIFT_INVERT_B = 1,
    LIS_KRYLOV_SHIFT_INVERT_A = 2
};

typedef struct {
    int n;
    int transform;
    double sigma;
    void *ctx;
    lis_krylov_matvec apply_a;
    lis_krylov_matvec apply_b;
    lis_krylov_matvec solve;
} lis_krylov_problem;

typedef struct {
    int nev;                       /* requested eigenpairs, 1 <= nev < n */
    int ncv;                       /* 0: choose min(n,3*nev+8); else >nev */
    int maxiter;                   /* 0: 300 expansion/restart cycles */
    double tol;                    /* 0: 1e-10 backward residual */
    const char *which;             /* LM (default), SM, LA, SA */
    double residual_scale;         /* 0: estimate ||A||; >0: user norm scale */
    double inner_tol;              /* 0: skip inverse-action verification */
    unsigned random_seed;          /* 0: stable default seed */
    const double *initial_vector;  /* NULL: generated deterministic start */
    int allow_semidefinite_metric;  /* 1: opt-in PSD B, mode SHIFT_INVERT_B */
} lis_krylov_options;

typedef struct {
    int nconv;                     /* number of returned residuals <= tol */
    int nreturned;                 /* number of ordered approximate pairs */
    int iterations;
    int restarts;
    int op_count;
    int metric_count;
    int solve_count;
    double max_residual;           /* worst residual among ALL returned pairs */
} lis_krylov_stats;

enum {
    LIS_KRYLOV_SUCCESS = 0,
    LIS_KRYLOV_BAD_ARGUMENT = 1,
    LIS_KRYLOV_NO_MEMORY = 2,
    LIS_KRYLOV_CALLBACK_ERROR = 3,
    LIS_KRYLOV_BREAKDOWN = 4,
    LIS_KRYLOV_NO_CONVERGENCE = 5,
    LIS_KRYLOV_INEXACT_SOLVE = 6,
    LIS_KRYLOV_SINGULAR_SHIFT = 7
};

/* Callback returning inertia of A - shift B from a symmetric indefinite
 * factorization: negative and (numerically) zero pivot counts.
 * The caller guarantees symmetric A, SPD B, and valid pivoting/inertia.
 * A nonzero callback return indicates a failed factorization/operation. */
typedef int (*lis_krylov_inertia_callback)(void *ctx, int n, double shift,
                                          int *negative, int *zero);

/* Count eigenvalues in the open interval (lower, upper) using the inertia
 * difference at two NONSINGULAR shifts. This does not compute eigenvectors
 * or guarantee complete spectral extraction. No factorization is built
 * inside LIS. Returns SINGULAR_SHIFT if a boundary has zero pivots. */
int lis_krylov_count_interval(int n, void *ctx,
                              lis_krylov_inertia_callback inertia,
                              double lower, double upper, int *count);

void lis_krylov_options_default(lis_krylov_options *options);

/* Return the leading nreturned approximate eigenpairs, ordered by 'which'.
 * eigenvalues[nev], eigenvectors[n*nev] (column-major), residuals[nev]
 * are caller-owned. At NO_CONVERGENCE the approximations remain available
 * in their ORIGINAL RITZ ORDER; inspect each residual separately rather
 * than assuming the first nconv are converged. Other errors invalidate
 * all results (nreturned=0).
 *
 * Each eigenvector has unit norm in the selected metric (I, B, or A).
 * Residual is ||A x-lambda B x||/max(||A x||+|lambda| ||B x||,
 *                                    residual_scale*||x||).
 * For residual_scale=0 a sampled estimate of ||A|| is used. This is a
 * backward-error indicator, not a proof of complete eigenvalue coverage.
 *
 * PSD metric mode requires B >= 0 (not necessarily invertible), a
 * nonsingular shifted operator, and finite eigenpairs in the generated
 * invariant subspace. Starting and restart directions are mapped through
 * (A-sigma B)^-1 B; the zero metric directions are rejected. A rank
 * smaller than nev can yield NO_CONVERGENCE with valid partial results.
 *
 * This API supports real double and SERIAL operators only. Applications
 * with MPI-distributed matrices must not pass rank-local vectors to it.
 */
int lis_krylov_symmetric(const lis_krylov_problem *problem,
                         const lis_krylov_options *options,
                         double *eigenvalues, double *eigenvectors,
                         double *residuals, lis_krylov_stats *stats);
#ifdef __cplusplus
}
#endif
#endif
