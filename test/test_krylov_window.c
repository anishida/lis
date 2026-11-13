/* Count-checked symmetric spectral window with analytic reference. */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 64

typedef struct {
    double metric[N];
    double eigenvalues[N];
    double sigma;
    int solves;
    int fail_solve;
} test_context;

static int apply_a(void *v, int n, const double *x, double *y)
{
    test_context *c = (test_context *)v;
    int i;
    if (n != N) return 1;
    for (i = 0; i < n; ++i)
        y[i] = c->metric[i] * c->eigenvalues[i] * x[i];
    return 0;
}

static int apply_b(void *v, int n, const double *x, double *y)
{
    test_context *c = (test_context *)v;
    int i;
    if (n != N) return 1;
    for (i = 0; i < n; ++i) y[i] = c->metric[i] * x[i];
    return 0;
}

static int solve_shifted(void *v, int n, const double *x, double *y)
{
    test_context *c = (test_context *)v;
    int i;
    if (n != N || c->fail_solve) return 1;
    ++c->solves;
    for (i = 0; i < n; ++i) {
        double pivot = c->metric[i] * (c->eigenvalues[i] - c->sigma);
        if (pivot == 0.0) return 1;
        y[i] = x[i] / pivot;
    }
    return 0;
}

static int inertia(void *v, int n, double shift, int *negative, int *zero)
{
    test_context *c = (test_context *)v;
    int i;
    if (n != N) return 1;
    *negative = *zero = 0;
    for (i = 0; i < n; ++i) {
        double pivot = c->metric[i] * (c->eigenvalues[i] - shift);
        if (pivot < 0.0) ++*negative;
        else if (pivot == 0.0) ++*zero;
    }
    return 0;
}

int main(void)
{
    test_context c = {{0}, {0}, 3.0, 0, 0};
    lis_krylov_problem problem = {0};
    lis_krylov_options options;
    lis_krylov_stats stats;
    double values[4], vectors[N*4], residuals[4];
    int i, rc, expected, found_two, found_four;

    for (i = 0; i < N; ++i) {
        c.metric[i] = 1.0 + 0.01 * i;
        c.eigenvalues[i] = i == 0 ? 1.0 :
                           i == 1 ? 2.0 :
                           i == 2 ? 4.0 :
                           i == 3 ? 5.0 : 10.0 + i;
    }
    problem.n = N;
    problem.transform = LIS_KRYLOV_SHIFT_INVERT_B;
    problem.ctx = &c;
    problem.sigma = c.sigma;
    problem.apply_a = apply_a;
    problem.apply_b = apply_b;
    problem.solve = solve_shifted;
    lis_krylov_options_default(&options);
    options.nev = 4;
    options.ncv = 16;
    options.maxiter = 300;
    options.tol = 1e-9;
    options.random_seed = 19;

    rc = lis_krylov_solve_interval(&problem, &options, inertia, 1.5, 4.5,
                                   values, vectors, residuals, &stats,
                                   &expected);
    if (rc != LIS_KRYLOV_SUCCESS || expected != 2 ||
        stats.nreturned != 2 || stats.nconv != 2) {
        fprintf(stderr,
                "FAIL window rc=%d expected=%d conv=%d returned=%d\n",
                rc, expected, stats.nconv, stats.nreturned);
        return 1;
    }
    found_two = found_four = 0;
    for (i = 0; i < 2; ++i) {
        if (fabs(values[i] - 2.0) < 1e-7) found_two = 1;
        if (fabs(values[i] - 4.0) < 1e-7) found_four = 1;
        if (residuals[i] > options.tol) return 1;
    }
    if (!found_two || !found_four) return 1;
    puts("PASS count-checked window extraction and residuals");

    options.nev = 1;
    rc = lis_krylov_solve_interval(&problem, &options, inertia, 1.5, 4.5,
                                   values, vectors, residuals, &stats,
                                   &expected);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT || expected != 2) return 1;
    options.nev = 4;

    problem.sigma = 2.9;
    rc = lis_krylov_solve_interval(&problem, &options, inertia, 1.5, 4.5,
                                   values, vectors, residuals, &stats,
                                   &expected);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT) return 1;

    /* The symmetric window (2,4.5) has midpoint 3.25, while the lower
     * endpoint is an exact eigenvalue and must be rejected by inertia. */
    problem.sigma = c.sigma = 3.25;
    rc = lis_krylov_solve_interval(&problem, &options, inertia, 2.0, 4.5,
                                   values, vectors, residuals, &stats,
                                   &expected);
    if (rc != LIS_KRYLOV_SINGULAR_SHIFT) return 1;

    /* An empty interval is a valid, certified result. */
    problem.sigma = c.sigma = 7.5;
    rc = lis_krylov_solve_interval(&problem, &options, inertia, 7.0, 8.0,
                                   values, vectors, residuals, &stats,
                                   &expected);
    if (rc != LIS_KRYLOV_SUCCESS || expected != 0 ||
        stats.nreturned != 0) return 1;

    puts("PASS capacity, midpoint, boundary and empty-window handling");
    return 0;
}
