/* Standalone, application-neutral benchmark for a 1-D FEM spring chain.
 * Real application benchmarks should replace the three callbacks below.
 *
 * cc -O3 -Iinclude src/esolver/lis_krylov.c \
 *    test/bench_krylov_adaptive.c -lm -o bench_krylov_adaptive
 * ./bench_krylov_adaptive [n] [nev] [initial_ncv] [ncv_max]
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
    int n;
    double *pivot, *upper, *tmp;
    unsigned long solves;
} chain;

static int apply_k(void *context, int n, const double *x, double *y)
{
    chain *c = (chain *)context;
    int i;
    if (n != c->n) return 1;
    for (i = 0; i < n; ++i)
        y[i] = 2.0*x[i] - (i ? x[i-1] : 0.0) -
               (i+1 < n ? x[i+1] : 0.0);
    return 0;
}

static int apply_m(void *context, int n, const double *x, double *y)
{
    chain *c = (chain *)context;
    int i;
    if (n != c->n) return 1;
    for (i = 0; i < n; ++i) y[i] = x[i];
    return 0;
}

/* Exact tridiagonal LDL^T factorization is prepared once, before timing. */
static int factor(chain *c)
{
    int i;
    c->pivot = (double *)malloc((size_t)c->n * sizeof(double));
    c->upper = (double *)malloc((size_t)c->n * sizeof(double));
    c->tmp = (double *)malloc((size_t)c->n * sizeof(double));
    if (!c->pivot || !c->upper || !c->tmp) return 1;
    c->pivot[0] = 2.0;
    for (i = 1; i < c->n; ++i) {
        c->upper[i-1] = -1.0 / c->pivot[i-1];
        c->pivot[i] = 2.0 + c->upper[i-1];
    }
    return 0;
}

static int inverse_k(void *context, int n, const double *x, double *y)
{
    chain *c = (chain *)context;
    int i;
    if (n != c->n) return 1;
    c->tmp[0] = x[0];
    for (i = 1; i < n; ++i)
        c->tmp[i] = x[i] + c->tmp[i-1] / c->pivot[i-1];
    y[n-1] = c->tmp[n-1] / c->pivot[n-1];
    for (i = n-2; i >= 0; --i)
        y[i] = c->tmp[i] / c->pivot[i] - c->upper[i]*y[i+1];
    ++c->solves;
    return 0;
}

static int measure(chain *c, int nev, int ncv, int cap, int adaptive,
                   const char *label)
{
    lis_krylov_problem p = {0};
    lis_krylov_options o;
    lis_krylov_stats s;
    double *values, *vectors, *residuals;
    clock_t start, finish;
    int rc, i, correct = 1;
    const double pi = 3.14159265358979323846;
    double elapsed;

    values = (double *)malloc((size_t)nev*sizeof(double));
    vectors = (double *)malloc((size_t)c->n*(size_t)nev*sizeof(double));
    residuals = (double *)malloc((size_t)nev*sizeof(double));
    if (!values || !vectors || !residuals) {
        free(values); free(vectors); free(residuals);
        return 1;
    }
    p.n = c->n;
    p.transform = LIS_KRYLOV_SHIFT_INVERT_B;
    p.ctx = c;
    p.apply_a = apply_k;
    p.apply_b = apply_m;
    p.solve = inverse_k;
    lis_krylov_options_default(&o);
    o.nev = nev;
    o.ncv = ncv;
    o.adaptive_ncv = adaptive;
    o.ncv_max = adaptive ? cap : 0;
    o.maxiter = 300;
    o.tol = 1e-9;

    c->solves = 0;
    start = clock();
    rc = lis_krylov_symmetric(&p, &o, values, vectors, residuals, &s);
    finish = clock();
    elapsed = (double)(finish-start) / CLOCKS_PER_SEC;

    if (rc == LIS_KRYLOV_SUCCESS) {
        for (i = 0; i < nev; ++i) {
            double exact = 2.0 - 2.0*cos((double)(i+1)*pi/(c->n+1));
            if (fabs(values[i]-exact) > 1e-5*fmax(1.0, fabs(exact)) ||
                residuals[i] > o.tol) correct = 0;
        }
    } else if (rc != LIS_KRYLOV_NO_CONVERGENCE) {
        correct = 0;
    }
    if (c->solves != (unsigned long)s.solve_count) correct = 0;

    printf("%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.6f,%.3e,%d\n",
           label, c->n, nev, rc,
           s.nconv, s.ncv_initial, s.ncv_final, s.ncv_growths,
           s.iterations, s.restarts, s.op_count, s.metric_count,
           s.solve_count, s.projected_solves,
           elapsed, s.max_residual, correct);
    free(values); free(vectors); free(residuals);
    return correct ? 0 : 1;
}

int main(int argc, char **argv)
{
    chain c = {0};
    int n = argc > 1 ? atoi(argv[1]) : 120;
    int nev = argc > 2 ? atoi(argv[2]) : 4;
    int ncv = argc > 3 ? atoi(argv[3]) : 8;
    int cap = argc > 4 ? atoi(argv[4]) : 24;
    int rc;
    if (n < 3 || nev < 1 || nev >= n ||
        ncv <= nev || ncv > n || cap < ncv || cap > n) {
        fprintf(stderr, "Usage: %s [n>=3] [1<=nev<n] "
                "[nev<ncv<=n] [ncv<=ncv_max<=n]\n", argv[0]);
        return 2;
    }
    c.n = n;
    if (factor(&c)) {
        fprintf(stderr, "Out of memory during factorization\n");
        free(c.pivot); free(c.upper); free(c.tmp);
        return 1;
    }
    puts("mode,n,nev,status,nconv,ncv_initial,ncv_final,ncv_growths,"
         "iterations,restarts,op_count,metric_count,solve_count,"
         "projected_solves,cpu_seconds,max_residual,valid");
    /* Equal-memory comparison: adaptive reserves storage for cap, so a
     * fixed-ncv=cap solve is an essential reference at that capacity. */
    rc = measure(&c, nev, ncv, cap, 0, "fixed_initial");
    if (cap > ncv + 1) {
        int mid = ncv + (cap - ncv) / 2;
        rc |= measure(&c, nev, mid, cap, 0, "fixed_mid");
    }
    rc |= measure(&c, nev, cap, cap, 0, "fixed_cap");
    rc |= measure(&c, nev, ncv, cap, 1, "adaptive");
    free(c.pivot); free(c.upper); free(c.tmp);
    return rc;
}
