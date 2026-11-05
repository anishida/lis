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
    const char *which;             /* LM (default), SM, LA, SA, BE;
                                    * optional :ASC/:DESC sorts output
                                    * ORIGINAL eigenvalues independently.
                                    * BE alternates extreme transformed
                                    * Ritz values, high end first. */
    double residual_scale;         /* 0: estimate ||A||; >0: user norm scale */
    double inner_tol;              /* 0: skip inverse-action verification */
    unsigned random_seed;          /* 0: stable default seed */
    const double *initial_vector;  /* NULL: generated deterministic start */
    int allow_semidefinite_metric; /* 1: opt-in PSD B, mode SHIFT_INVERT_B */
    int hard_locking;              /* 1: freeze and deflate converged pairs */
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

/* Enumerate a symmetric open interval with a user-prepared midpoint shift.
 * Requires a symmetric generalized pencil with SPD B, an inertia callback
 * at the nonsingular interval boundaries, and p->sigma equal to
 * lower+0.5*(upper-lower). The caller supplies a nonsingular shifted solve
 * at that sigma. options->nev is the output capacity; ncv must exceed the
 * inertia count. LM shift-invert mode is used regardless of options->which.
 * SUCCESS certifies that the number of converged, in-window pairs equals
 * the independently supplied inertia count; otherwise returns a status
 * (including NO_CONVERGENCE with partial results). Exact midpoint poles
 * cannot be inverted and must be handled by the application. */
int lis_krylov_solve_interval(const lis_krylov_problem *problem,
                              const lis_krylov_options *options,
                              lis_krylov_inertia_callback inertia,
                              double lower, double upper,
                              double *values, double *vectors,
                              double *residuals, lis_krylov_stats *stats,
                              int *expected_count);

void lis_krylov_options_default(lis_krylov_options *options);

/* Experimental, versioned extension API.
 * The ORIGINAL lis_krylov_options and lis_krylov_stats binary layouts
 * are unchanged. New controls live in separate structs and entry points.
 * Set both struct_size and version by calling the default initializer.
 */
#define LIS_KRYLOV_EXT_VERSION 1u

typedef struct lis_krylov_ext_options {
    unsigned int struct_size;       /* sizeof(lis_krylov_ext_options) */
    unsigned int version;           /* LIS_KRYLOV_EXT_VERSION */
    int adaptive_ncv;               /* opt-in active subspace growth */
    int ncv_max;                    /* 0: min(n,2*ncv) when adaptive */
    int restart_directions;         /* 1: original scalar-seed restart */
} lis_krylov_ext_options;

typedef struct lis_krylov_ext_stats {
    unsigned int struct_size;       /* sizeof(lis_krylov_ext_stats) */
    unsigned int version;           /* LIS_KRYLOV_EXT_VERSION */
    int ncv_initial;
    int ncv_final;
    int ncv_growths;
    int basis_vectors;
    int projected_solves;
    int extra_seed_attempts;        /* candidates tried after legacy seed */
    int extra_seed_accepted;        /* extra candidates actually accepted */
} lis_krylov_ext_stats;

void lis_krylov_ext_options_default(lis_krylov_ext_options *options);
void lis_krylov_ext_stats_default(lis_krylov_ext_stats *stats);

/* Opaque, reusable allocation for repeated solves with compatible
 * dimensions, transform, nev, initial ncv, hard locking and capacity.
 * If adaptive sizing is enabled at creation, the maximum workspace
 * is allocated ONCE and reused as the active subspace grows.
 * Operators, context and sigma are supplied anew for every call.
 * The caller MUST prepare the correct shifted inverse for each sigma.
 * A workspace is NOT safe for concurrent calls.
 */
typedef struct lis_krylov_workspace lis_krylov_workspace;
int lis_krylov_workspace_create(const lis_krylov_problem *problem,
                                const lis_krylov_options *options,
                                const lis_krylov_ext_options *ext_options,
                                lis_krylov_workspace **workspace);
void lis_krylov_workspace_destroy(lis_krylov_workspace *workspace);

/* Single composable entry point: optional sizing, restart enrichment,
 * diagnostics and optional workspace. NULL ext_options selects defaults.
 * NULL workspace uses per-call allocations. NULL ext_stats suppresses
 * extended counters. Existing lis_krylov_symmetric() is unchanged.
 *
 * No shared Ritz state survives between calls: reuse concerns allocation
 * only. A failed solve does not invalidate the workspace. This is not
 * a full block Lanczos algorithm; callbacks remain scalar. */
int lis_krylov_symmetric_ex(const lis_krylov_problem *problem,
                            const lis_krylov_options *options,
                            const lis_krylov_ext_options *ext_options,
                            lis_krylov_workspace *workspace,
                            double *eigenvalues, double *eigenvectors,
                            double *residuals, lis_krylov_stats *stats,
                            lis_krylov_ext_stats *ext_stats);


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
