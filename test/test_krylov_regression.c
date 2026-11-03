/* Numerical regressions for scale, partial results and callback failures. */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 48
#define NEV 4

typedef struct {
    double a[N], b[N];
    double sigma;
    int failed_callback;
    int inverses;
    double inverse_error;
} context;

static int matrix_a(void *opaque, int n, const double *x, double *y)
{
    context *c = (context *)opaque;
    int i;
    if (c->failed_callback == 1) return -1;
    for (i = 0; i < n; ++i) y[i] = c->a[i] * x[i];
    return 0;
}

static int matrix_b(void *opaque, int n, const double *x, double *y)
{
    context *c = (context *)opaque;
    int i;
    if (c->failed_callback == 2) return -1;
    for (i = 0; i < n; ++i) y[i] = c->b[i] * x[i];
    return 0;
}

static int shifted_solve(void *opaque, int n, const double *x, double *y)
{
    context *c = (context *)opaque;
    int i;
    if (c->failed_callback == 3) return -1;
    ++c->inverses;
    for (i = 0; i < n; ++i)
        y[i] = (1.0 + c->inverse_error) * x[i] /
               (c->a[i] - c->sigma * c->b[i]);
    return 0;
}

static int run_case(const char *label, context *c,
                    int transform, const char *which,
                    double sigma, const double *expected, int nev,
                    double relative_eigen_tolerance)
{
    lis_krylov_problem p;
    lis_krylov_options opt;
    lis_krylov_stats stats;
    double eigen[NEV], vectors[N * NEV], residual[NEV];
    int err, i, j;

    memset(&p, 0, sizeof(p));
    lis_krylov_options_default(&opt);
    p.n = N;
    p.transform = transform;
    p.sigma = sigma;
    p.ctx = c;
    p.apply_a = matrix_a;
    p.apply_b = matrix_b;
    p.solve = shifted_solve;
    opt.nev = nev;
    opt.ncv = 32;
    opt.tol = 1.0e-9;
    opt.maxiter = 600;
    opt.which = which;
    c->sigma = sigma;
    err = lis_krylov_symmetric(&p, &opt, eigen, vectors, residual, &stats);
    if (err || stats.nconv != nev || stats.nreturned != nev) {
        fprintf(stderr, "FAIL %-18s: rc=%d conv=%d returned=%d iter=%d worst=%g\n",
                label, err, stats.nconv, stats.nreturned,
                stats.iterations, stats.max_residual);
        return 1;
    }
    for (i = 0; i < nev; ++i) {
        double best = HUGE_VAL;
        for (j = 0; j < nev; ++j) {
            double delta = fabs(eigen[j] - expected[i]);
            if (delta < best) best = delta;
        }
        if (best > relative_eigen_tolerance *
                   fmax(1.0, fabs(expected[i])) ||
            residual[i] > opt.tol) {
            fprintf(stderr,
                    "FAIL %-18s: mode=%d eig=%g target=%g delta=%g residual=%g\n",
                    label, i, eigen[i], expected[i], best, residual[i]);
            return 1;
        }
    }
    printf("PASS %-18s: iter=%d restart=%d solve=%d worst=%.3g\n",
           label, stats.iterations, stats.restarts,
           stats.solve_count, stats.max_residual);
    return 0;
}

static int failures(void)
{
    context c = {0};
    lis_krylov_problem p;
    lis_krylov_options opt;
    lis_krylov_stats stats;
    double eigen[NEV], vectors[N * NEV], residual[NEV];
    double zero[N] = {0};
    int i, err;

    for (i = 0; i < N; ++i) {
        c.a[i] = (double)(i + 2);
        c.b[i] = 1.0;
    }
    memset(&p, 0, sizeof(p));
    p.n = N;
    p.transform = LIS_KRYLOV_SHIFT_INVERT_B;
    p.ctx = &c;
    p.apply_a = matrix_a;
    p.apply_b = matrix_b;
    p.solve = shifted_solve;
    lis_krylov_options_default(&opt);
    opt.nev = 3;
    opt.ncv = 24;
    opt.initial_vector = zero;
    err = lis_krylov_symmetric(&p, &opt, eigen, vectors, residual, &stats);
    if (err != LIS_KRYLOV_BAD_ARGUMENT) return 1;
    opt.initial_vector = NULL;

    for (i = 1; i <= 3; ++i) {
        c.failed_callback = i;
        err = lis_krylov_symmetric(&p, &opt, eigen, vectors, residual, &stats);
        if (err != LIS_KRYLOV_CALLBACK_ERROR || stats.nreturned != 0) {
            fprintf(stderr, "FAIL callback %d err=%d\n", i, err);
            return 1;
        }
    }
    c.failed_callback = 0;
    c.inverse_error = 0.01;
    opt.inner_tol = 1.0e-6;
    err = lis_krylov_symmetric(&p, &opt, eigen, vectors, residual, &stats);
    if (err != LIS_KRYLOV_INEXACT_SOLVE || stats.nreturned != 0) {
        fprintf(stderr, "FAIL inaccurate inverse: got %d\n", err);
        return 1;
    }
    c.inverse_error = 0.0;
    opt.inner_tol = 1.0e-7;
    opt.maxiter = 1;
    err = lis_krylov_symmetric(&p, &opt, eigen, vectors, residual, &stats);
    if (err != LIS_KRYLOV_NO_CONVERGENCE ||
        stats.nreturned != 1 || stats.nconv > stats.nreturned ||
        !isfinite(residual[0]) || stats.max_residual != residual[0]) {
        fprintf(stderr, "FAIL partial results: err=%d conv=%d returned=%d worst=%g\n",
                err, stats.nconv, stats.nreturned, stats.max_residual);
        return 1;
    }
    printf("PASS callback errors, inexact inverse and partial results\n");
    return 0;
}

int main(void)
{
    context c = {0};
    double expected[NEV];
    int i;

    for (i = 0; i < N; ++i) {
        c.a[i] = 1.0e9 * (double)(i + 2);
        c.b[i] = 1.0;
    }
    for (i = 0; i < NEV; ++i) expected[i] = 1.0e9 * (i + 2);
    if (run_case("scaled 1e9", &c, LIS_KRYLOV_SHIFT_INVERT_B,
                 "LM", 0.0, expected, NEV, 1.0e-8)) return 1;

    for (i = 0; i < N; ++i) {
        c.a[i] = (double)(i + 1);
        c.b[i] = 1.0;
    }
    c.a[0] = 0.0;
    c.a[1] = 0.0;
    c.a[2] = 1.0e-5;
    expected[0] = 0.0;
    expected[1] = 0.0;
    expected[2] = 1.0e-5;
    if (run_case("zero eigenvalues", &c, LIS_KRYLOV_SHIFT_INVERT_B,
                 "LM", -0.01, expected, 3, 1.0e-7)) return 1;

    /* Test all selection rules with full-window iterations. */
    for (i = 0; i < N; ++i) {
        c.a[i] = (i & 1) ? -(double)(i + 1) : (double)(i + 1);
        c.b[i] = 1.0;
    }
    expected[0] = 1.0;
    expected[1] = 3.0;
    expected[2] = -2.0;
    if (run_case("SM", &c, LIS_KRYLOV_STANDARD,
                 "SM", 0.0, expected, 3, 1.0e-8)) return 1;
    expected[0] = 47.0;
    expected[1] = 45.0;
    expected[2] = 43.0;
    if (run_case("LA", &c, LIS_KRYLOV_STANDARD,
                 "LA", 0.0, expected, 3, 1.0e-8)) return 1;
    expected[0] = -48.0;
    expected[1] = -46.0;
    expected[2] = -44.0;
    if (run_case("SA", &c, LIS_KRYLOV_STANDARD,
                 "SA", 0.0, expected, 3, 1.0e-8)) return 1;

    /* Double eigenvalues: returned eigenvectors must be independent. */
    for (i = 0; i < N; ++i) {
        c.a[i] = i < 4 ? 100.0 : (double)(N - i);
        c.b[i] = 1.0;
    }
    for (i = 0; i < NEV; ++i) expected[i] = 100.0;
    if (run_case("multiplicity x4", &c, LIS_KRYLOV_STANDARD,
                 "LA", 0.0, expected, NEV, 1.0e-7)) return 1;
    if (failures()) return 1;
    return 0;
}
