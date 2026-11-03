/* Regression for OpenMP Krylov vector kernels, DGKS, restart and locking. */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define N 65536
#define NEV 3
#define NCV 8

typedef struct { int solves; } test_context;

static double eigenvalue(int i)
{
    return i < NEV ? 0.25 * (i + 1) : 15.0 + 0.01 * (i % 17);
}

static double mass(int i)
{
    return 1.0 + 0.125 * (i % 5);
}

static int apply_a(void *unused, int n, const double *x, double *y)
{
    int i;
    (void)unused;
    if (n != N) return 1;
    for (i = 0; i < n; ++i) y[i] = mass(i) * eigenvalue(i) * x[i];
    return 0;
}

static int apply_b(void *unused, int n, const double *x, double *y)
{
    int i;
    (void)unused;
    if (n != N) return 1;
    for (i = 0; i < n; ++i) y[i] = mass(i) * x[i];
    return 0;
}

static int solve_shifted(void *state, int n, const double *x, double *y)
{
    test_context *ctx = (test_context *)state;
    int i;
    if (n != N) return 1;
    ++ctx->solves;
    for (i = 0; i < n; ++i)
        y[i] = x[i] / (mass(i) * eigenvalue(i));
    return 0;
}

static int run_case(int threads, int hard_locking, double *reference)
{
    test_context ctx = {0};
    lis_krylov_problem p = {0};
    lis_krylov_options options;
    lis_krylov_stats stats;
    double values[NEV], residuals[NEV], *vectors;
    int status, i, j, k;

#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#else
    (void)threads;
#endif
    vectors = (double *)calloc((size_t)N * NEV, sizeof(double));
    if (!vectors) return 1;
    lis_krylov_options_default(&options);
    options.nev = NEV;
    options.ncv = NCV;
    options.maxiter = 150;
    options.tol = 1e-9;
    options.random_seed = 43;
    options.hard_locking = hard_locking;

    p.n = N;
    p.ctx = &ctx;
    p.transform = LIS_KRYLOV_SHIFT_INVERT_B;
    p.sigma = 0.0;
    p.apply_a = apply_a;
    p.apply_b = apply_b;
    p.solve = solve_shifted;

    status = lis_krylov_symmetric(&p, &options, values,
                                  vectors, residuals, &stats);
    if (status != LIS_KRYLOV_SUCCESS || stats.nconv != NEV ||
        stats.nreturned != NEV || stats.restarts < 1 ||
        stats.solve_count != ctx.solves) {
        fprintf(stderr,
                "FAIL threads=%d hardlock=%d status=%d nconv=%d "
                "restarts=%d solves=%d/%d\n",
                threads, hard_locking, status, stats.nconv, stats.restarts,
                stats.solve_count, ctx.solves);
        free(vectors);
        return 1;
    }
    for (i = 0; i < NEV; ++i) {
        const double expected = 0.25 * (i + 1);
        if (fabs(values[i] - expected) > 1e-7 ||
            residuals[i] > options.tol ||
            (threads != 1 && reference && fabs(values[i] - reference[i]) > 1e-8)) {
            fprintf(stderr, "FAIL threads=%d pair=%d value=%.17g "
                    "residual=%.3g\n", threads, i, values[i], residuals[i]);
            free(vectors);
            return 1;
        }
        for (j = 0; j <= i; ++j) {
            double product = 0.0;
            for (k = 0; k < N; ++k)
                product += mass(k) * vectors[(size_t)i * N + k] *
                                     vectors[(size_t)j * N + k];
            if (fabs(product - (i == j ? 1.0 : 0.0)) > 1e-6) {
                fprintf(stderr, "FAIL threads=%d metric dot %d/%d: %.5g\n",
                        threads, i, j, product);
                free(vectors);
                return 1;
            }
        }
    }
    if (threads == 1 && reference)
        for (i = 0; i < NEV; ++i) reference[i] = values[i];
    printf("PASS threads=%d hardlock=%d n=%d restarts=%d "
           "solves=%d maxres=%.3g\n",
           threads, hard_locking, N, stats.restarts,
           ctx.solves, stats.max_residual);
    free(vectors);
    return 0;
}

int main(void)
{
#ifdef _OPENMP
    const int threads[] = {1, 2, 4, 8};
    double ref_plain[NEV] = {0}, ref_locked[NEV] = {0};
    int i;
    for (i = 0; i < 4; ++i)
        if (run_case(threads[i], 0, ref_plain) ||
            run_case(threads[i], 1, ref_locked))
            return 1;
#else
    if (run_case(1, 0, NULL) || run_case(1, 1, NULL)) return 1;
#endif
    return 0;
}
