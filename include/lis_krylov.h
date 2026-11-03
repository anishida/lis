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
 * SHIFT_INVERT_B: B positive definite; OP=(A-sigma B)^-1 B; B metric.
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
    LIS_KRYLOV_INEXACT_SOLVE = 6
};

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
