/* Regression for opt-in adaptive, operator-only Krylov restarts.
 * This file deliberately uses no LIS matrices or application-specific API.
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 64
#define NEV 4

static int diagonal(void *ctx, int n, const double *x, double *y)
{
    int i;
    (void)ctx;
    if (n != N) return 1;
    for (i = 0; i < n; ++i) y[i] = (double)(i + 1) * x[i];
    return 0;
}

static int run_case(int adaptive, int ncv_max, int expected_growth)
{
    lis_krylov_problem p = {0};
    lis_krylov_options o;
    lis_krylov_stats stats;
    double values[NEV], vectors[N*NEV], residuals[NEV];
    int i, j, rc, converged = 0;

    p.n = N;
    p.transform = LIS_KRYLOV_STANDARD;
    p.apply_a = diagonal;
    lis_krylov_options_default(&o);
    o.nev = NEV;
    o.ncv = 6;
    o.ncv_max = ncv_max;
    o.adaptive_ncv = adaptive;
    o.which = "LA";
    o.tol = 1e-12;
    o.maxiter = 200;

    rc = lis_krylov_symmetric(&p, &o, values, vectors,
                              residuals, &stats);
    if (rc != LIS_KRYLOV_SUCCESS &&
        rc != LIS_KRYLOV_NO_CONVERGENCE) {
        fprintf(stderr, "FAIL adaptive=%d solver rc=%d\n", adaptive, rc);
        return 1;
    }
    if (stats.ncv_initial != 6 ||
        stats.ncv_final < 6 ||
        stats.ncv_final > (ncv_max ? ncv_max : 12) ||
        !!stats.ncv_growths != !!expected_growth ||
        stats.ncv_growths > 0 != (stats.ncv_final > 6) ||
        stats.projected_solves != stats.iterations ||
        stats.basis_vectors != stats.op_count ||
        stats.nreturned < 1 || stats.nreturned > NEV) {
        fprintf(stderr, "FAIL adaptive=%d accounting: ncv=%d..%d "
                "growths=%d projections=%d iterations=%d basis=%d op=%d "
                "returned=%d\n", adaptive, stats.ncv_initial,
                stats.ncv_final, stats.ncv_growths,
                stats.projected_solves, stats.iterations,
                stats.basis_vectors, stats.op_count, stats.nreturned);
        return 1;
    }
    for (j = 0; j < stats.nreturned; ++j) {
        double norm = 0.0, numerator = 0.0, scale = 0.0;
        for (i = 0; i < N; ++i) {
            double x = vectors[(size_t)j*N + i];
            double r = (i + 1 - values[j])*x;
            norm += x*x;
            numerator += r*r;
            scale += ((i+1) + fabs(values[j])) * fabs(x);
        }
        /* A NO_CONVERGENCE return legitimately exposes partial Ritz pairs.
         * Demand a small true residual only for pairs marked converged. */
        if (fabs(norm - 1.0) > 1e-6 ||
            (residuals[j] <= o.tol &&
             sqrt(numerator) > 1e-8 * fmax(scale, 1.0)) ||
            !(residuals[j] >= 0.0) || !isfinite(residuals[j])) {
            fprintf(stderr, "FAIL adaptive=%d eigenpair %d "
                    "lambda=%.15g residual=%.3g\n",
                    adaptive, j, values[j], residuals[j]);
            return 1;
        }
        if (residuals[j] <= o.tol) ++converged;
    }
    if (converged != stats.nconv ||
        (rc == LIS_KRYLOV_SUCCESS && stats.nconv != NEV)) {
        fprintf(stderr, "FAIL adaptive=%d convergence count\n", adaptive);
        return 1;
    }
    printf("PASS adaptive=%d rc=%d ncv=%d..%d growths=%d "
           "restarts=%d operator_calls=%d converged=%d\n",
           adaptive, rc, stats.ncv_initial, stats.ncv_final,
           stats.ncv_growths, stats.restarts, stats.op_count, stats.nconv);
    return 0;
}

static int identity(void *ctx, int n, const double *x, double *y)
{
    int i;
    (void)ctx;
    if (n != N) return 1;
    for (i = 0; i < n; ++i) y[i] = x[i];
    return 0;
}

static int inverse(void *ctx, int n, const double *x, double *y)
{
    int i;
    (void)ctx;
    if (n != N) return 1;
    for (i = 0; i < n; ++i) y[i] = x[i] / (double)(i+1);
    return 0;
}

static int run_shifted_hardlock(void)
{
    lis_krylov_problem p = {0};
    lis_krylov_options o;
    lis_krylov_stats stats;
    double values[NEV], vectors[N*NEV], residuals[NEV];
    int rc, i;

    p.n = N;
    p.transform = LIS_KRYLOV_SHIFT_INVERT_B;
    p.apply_a = diagonal;
    p.apply_b = identity;
    p.solve = inverse;
    lis_krylov_options_default(&o);
    o.nev = NEV;
    o.ncv = 6;
    o.ncv_max = 20;
    o.adaptive_ncv = 1;
    o.hard_locking = 1;
    o.tol = 1e-10;
    o.maxiter = 300;
    rc = lis_krylov_symmetric(&p, &o, values, vectors,
                              residuals, &stats);
    if ((rc != LIS_KRYLOV_SUCCESS &&
         rc != LIS_KRYLOV_NO_CONVERGENCE) ||
        stats.op_count != stats.solve_count || stats.metric_count < 1 ||
        stats.ncv_initial != 6 || stats.ncv_final > 20) {
        fprintf(stderr, "FAIL shifted hardlock rc=%d solves=%d ncv=%d\n",
                rc, stats.solve_count, stats.ncv_final);
        return 1;
    }
    if (rc == LIS_KRYLOV_SUCCESS) {
        for (i = 0; i < NEV; ++i)
            if (fabs(values[i] - (double)(i+1)) > 1e-6 ||
                residuals[i] > o.tol)
                return 1;
    }
    printf("PASS shifted hardlock rc=%d ncv=%d growths=%d solves=%d\n",
           rc, stats.ncv_final, stats.ncv_growths, stats.solve_count);
    return 0;
}

int main(void)
{
    lis_krylov_problem p = {0};
    lis_krylov_options o;
    lis_krylov_stats stats;
    double values[NEV], vectors[N*NEV], residuals[NEV];

    if (run_case(0, 0, 0) ||
        run_case(1, 24, 1) ||
        run_case(1, 6, 0) ||
        run_shifted_hardlock())
        return 1;

    p.n = N;
    p.apply_a = diagonal;
    lis_krylov_options_default(&o);
    o.nev = NEV;
    o.ncv = 6;
    o.adaptive_ncv = 1;
    o.ncv_max = 5;
    if (lis_krylov_symmetric(&p, &o, values, vectors, residuals,
                             &stats) != LIS_KRYLOV_BAD_ARGUMENT)
        return 1;
    o.ncv_max = N+1;
    if (lis_krylov_symmetric(&p, &o, values, vectors, residuals,
                             &stats) != LIS_KRYLOV_BAD_ARGUMENT)
        return 1;
    o.ncv_max = 0;
    o.adaptive_ncv = 2;
    if (lis_krylov_symmetric(&p, &o, values, vectors, residuals,
                             &stats) != LIS_KRYLOV_BAD_ARGUMENT)
        return 1;
    o.adaptive_ncv = 0;
    o.ncv_max = 12;
    if (lis_krylov_symmetric(&p, &o, values, vectors, residuals,
                             &stats) != LIS_KRYLOV_BAD_ARGUMENT)
        return 1;
    puts("PASS invalid adaptive controls");
    return 0;
}
