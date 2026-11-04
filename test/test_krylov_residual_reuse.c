/* Regression for reusing verified physical Ritz residuals.
 * Instrument application callbacks, check physical eigenpairs, and emit
 * deterministic records for baseline-vs-optimized comparison in CI.
 * The small hard-locking callback budget below is deliberately asserted
 * in the installed make check regression, not only by a CI script.
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 64
#define NEV 4

typedef struct {
    double a[N], b[N], sigma;
    long calls_a, calls_b, calls_solve;
    long fail_a_at;
} context;

static int apply_a(void *v, int n, const double *x, double *y)
{
    context *c = (context *)v;
    int k;
    if (n != N) return 1;
    ++c->calls_a;
    if (c->fail_a_at && c->calls_a == c->fail_a_at) return 1;
    for (k = 0; k < n; ++k) y[k] = c->a[k] * x[k];
    return 0;
}
static int apply_b(void *v, int n, const double *x, double *y)
{
    context *c = (context *)v;
    int k;
    if (n != N) return 1;
    ++c->calls_b;
    for (k = 0; k < n; ++k) y[k] = c->b[k] * x[k];
    return 0;
}
static int shifted_solve(void *v, int n, const double *x, double *y)
{
    context *c = (context *)v;
    int k;
    if (n != N) return 1;
    ++c->calls_solve;
    for (k = 0; k < n; ++k) {
        double d = c->a[k] - c->sigma * c->b[k];
        if (fabs(d) < 1e-13) return 1;
        y[k] = x[k] / d;
    }
    return 0;
}

static void setup(context *c, lis_krylov_problem *p,
                  lis_krylov_options *opt, int transform,
                  int locking, int fixed_scale, int long_run, int psd)
{
    int k;
    memset(c, 0, sizeof(*c));
    memset(p, 0, sizeof(*p));
    for (k = 0; k < N; ++k) {
        c->a[k] = 1.0 + 0.027 * k;
        c->b[k] = psd && k >= 48 ? 0.0 : 1.0 + 0.003 * (k % 7);
    }
    c->sigma = transform == LIS_KRYLOV_SHIFT_INVERT_B ? 1.08 : 0.41;
    p->n = N;
    p->ctx = c;
    p->transform = transform;
    p->sigma = transform == LIS_KRYLOV_STANDARD ? 0.0 : c->sigma;
    p->apply_a = apply_a;
    p->apply_b = apply_b;
    p->solve = shifted_solve;
    lis_krylov_options_default(opt);
    opt->nev = NEV;
    opt->ncv = long_run ? 24 : 10;
    opt->maxiter = long_run ? 250 : 24;
    opt->tol = 1e-9;
    opt->which = transform == LIS_KRYLOV_STANDARD ? "LA" : "LM";
    opt->random_seed = 3;
    opt->hard_locking = locking;
    opt->residual_scale = fixed_scale ? 4.0 : 0.0;
    opt->allow_semidefinite_metric = psd;
}

/* Use the true A and B (or I in standard mode), without callbacks. */
static int verify(const context *c, int transform,
                  const double *v, double lambda, double reported,
                  double scale, int converged)
{
    double nr2 = 0.0, na2 = 0.0, nb2 = 0.0, nx2 = 0.0, denom, direct;
    int k;
    for (k = 0; k < N; ++k) {
        double ax = c->a[k] * v[k];
        double bx = (transform == LIS_KRYLOV_STANDARD ?
                     1.0 : c->b[k]) * v[k];
        double r = ax - lambda * bx;
        nr2 += r * r;
        na2 += ax * ax;
        nb2 += bx * bx;
        nx2 += v[k] * v[k];
    }
    denom = sqrt(na2) + fabs(lambda) * sqrt(nb2);
    if (scale * sqrt(nx2) > denom) denom = scale * sqrt(nx2);
    direct = denom > 0.0 ? sqrt(nr2) / denom : sqrt(nr2);
    if (!isfinite(direct) || !isfinite(reported) ||
        direct > 20.0 * reported + 5e-10 ||
        (converged && direct > 2e-7)) {
        fprintf(stderr, "FAIL direct residual %.17g reported %.17g\n",
                direct, reported);
        return 1;
    }
    return 0;
}

static int run(int transform, int locking, int fixed_scale,
               int long_run, int api, int psd)
{
    context c;
    lis_krylov_problem p;
    lis_krylov_options opt;
    lis_krylov_stats st;
    lis_krylov_ext_options ext;
    lis_krylov_ext_stats extra;
    lis_krylov_workspace *w = NULL;
    double values[NEV] = {0}, residuals[NEV] = {0};
    double vectors[N*NEV] = {0};
    int rc, k, i;
    setup(&c, &p, &opt, transform, locking, fixed_scale, long_run, psd);
    lis_krylov_ext_options_default(&ext);
    lis_krylov_ext_stats_default(&extra);
    if (api == 2) {
        rc = lis_krylov_workspace_create(&p, &opt, &ext, &w);
        if (rc != LIS_KRYLOV_SUCCESS || !w) {
            fprintf(stderr, "FAIL workspace_create: %d\n", rc);
            return 1;
        }
    }
    if (api == 0)
        rc = lis_krylov_symmetric(&p, &opt, values, vectors, residuals, &st);
    else
        rc = lis_krylov_symmetric_ex(&p, &opt, &ext, w, values, vectors,
                                     residuals, &st, &extra);
    lis_krylov_workspace_destroy(w);
    if (rc != LIS_KRYLOV_SUCCESS && rc != LIS_KRYLOV_NO_CONVERGENCE) {
        fprintf(stderr, "FAIL solver rc %d transform=%d lock=%d psd=%d\n",
                rc, transform, locking, psd);
        return 1;
    }
    if (st.nreturned < 0 || st.nreturned > NEV ||
        st.nconv < 0 || st.nconv > st.nreturned) return 1;
    for (k = 0; k < st.nreturned; ++k) {
        if (verify(&c, transform, vectors + (size_t)k*N,
                   values[k], residuals[k], opt.residual_scale,
                   residuals[k] <= opt.tol)) return 1;
    }
    /* Preserve all public stats, actual callbacks and all vector
     * coordinates in the reference-comparison records. */
    printf("CASE %d %d %d %d %d %d %d %d %d %d %d %d %d %ld %ld %ld",
           transform, locking, fixed_scale, long_run, api, psd, rc,
           st.nreturned, st.nconv, st.iterations, st.restarts,
           st.op_count, st.metric_count,
           c.calls_a, c.calls_b, c.calls_solve);
    for (k = 0; k < st.nreturned; ++k)
        printf(" %.17g %.17g", values[k], residuals[k]);
    printf("\nVECTOR");
    for (k = 0; k < st.nreturned; ++k)
        for (i = 0; i < N; ++i)
            printf(" %.17g", vectors[(size_t)k*N+i]);
    printf("\n");

    /* Deterministic operator-call regression: these thresholds are
     * below the unoptimized develop counts but above optimized counts.
     * They are intentionally generous across compilers. Do not hide
     * this check in a temporary workflow: make check enforces it.
     *
     * Baseline / optimized apply_a, apply_b on GCC:
     * standard 1228/0 -> 762/0; B-invert 186/394 -> 95/303;
     * PSD B-invert 164/354 -> 84/274; A-invert 1071/528 -> 809/266.
     */
#ifndef KRYLOV_BASELINE_COMPARISON
    if (api == 0 && locking && fixed_scale && long_run) {
        long max_a = 0, max_b = 0;
        if (transform == LIS_KRYLOV_STANDARD) {
            max_a = 995; max_b = 0;
        } else if (transform == LIS_KRYLOV_SHIFT_INVERT_B) {
            max_a = psd ? 125 : 150;
            max_b = psd ? 320 : 350;
        } else {
            max_a = 960; max_b = 390;
        }
        if (c.calls_a > max_a || c.calls_b > max_b) {
            fprintf(stderr,
                    "FAIL residual cache callback budget: A=%ld/%ld B=%ld/%ld\n",
                    c.calls_a, max_a, c.calls_b, max_b);
            return 1;
        }
    }
#endif
    return 0;
}

/* A-metric shift-invert with B=0 gives OP=I and mapped Ritz values at
 * the nonrepresentable pole mu=1. No partial eigenpair may be exposed. */
static int pole_case(int api)
{
    context c;
    lis_krylov_problem p;
    lis_krylov_options opt;
    lis_krylov_stats st;
    lis_krylov_ext_options ext;
    lis_krylov_ext_stats extra;
    lis_krylov_workspace *w = NULL;
    double values[NEV] = {0}, residuals[NEV] = {0};
    double vectors[N*NEV] = {0};
    int rc, k;
    setup(&c, &p, &opt, LIS_KRYLOV_SHIFT_INVERT_A, 1, 1, 0, 0);
    for (k = 0; k < N; ++k) c.b[k] = 0.0;
    c.sigma = p.sigma = -1.0;
    opt.maxiter = 16;
    opt.ncv = 16;
    opt.nev = 3;
    lis_krylov_ext_options_default(&ext);
    lis_krylov_ext_stats_default(&extra);
    if (api == 2) {
        if (lis_krylov_workspace_create(&p, &opt, &ext, &w) ||
            !w) return 1;
    }
    if (api == 0)
        rc = lis_krylov_symmetric(&p, &opt, values, vectors, residuals, &st);
    else
        rc = lis_krylov_symmetric_ex(&p, &opt, &ext, w,
                                     values, vectors, residuals,
                                     &st, &extra);
    lis_krylov_workspace_destroy(w);
    if (rc != LIS_KRYLOV_NO_CONVERGENCE ||
        st.nreturned != 0 || st.nconv != 0) {
        fprintf(stderr, "FAIL hard-locking nonrepresentable Ritz: %d %d %d\n",
                rc, st.nreturned, st.nconv);
        return 1;
    }
    printf("POLE %d %d %d %d\n", api, rc, st.nreturned, st.nconv);
    return 0;
}

/* The failure path must invalidate results and leave reusable storage
 * available for a fresh solve (resetting callback failure first). */
static int callback_failure_case(void)
{
    context c;
    lis_krylov_problem p;
    lis_krylov_options opt;
    lis_krylov_stats st;
    lis_krylov_ext_options ext;
    lis_krylov_ext_stats extra;
    lis_krylov_workspace *w = NULL;
    double values[NEV] = {0}, residuals[NEV] = {0};
    double vectors[N*NEV] = {0};
    int rc;
    setup(&c, &p, &opt, LIS_KRYLOV_SHIFT_INVERT_B, 1, 1, 1, 0);
    lis_krylov_ext_options_default(&ext);
    lis_krylov_ext_stats_default(&extra);
    if (lis_krylov_workspace_create(&p, &opt, &ext, &w) ||
        !w) return 1;
    c.fail_a_at = 4;
    rc = lis_krylov_symmetric_ex(&p, &opt, &ext, w,
                                 values, vectors, residuals, &st, &extra);
    if (rc != LIS_KRYLOV_CALLBACK_ERROR || st.nreturned != 0) {
        fprintf(stderr, "FAIL callback error: %d nreturned=%d\n",
                rc, st.nreturned);
        lis_krylov_workspace_destroy(w);
        return 1;
    }
    c.fail_a_at = 0;
    rc = lis_krylov_symmetric_ex(&p, &opt, &ext, w,
                                 values, vectors, residuals, &st, &extra);
    lis_krylov_workspace_destroy(w);
    if (rc != LIS_KRYLOV_SUCCESS || st.nconv != NEV) {
        fprintf(stderr, "FAIL workspace reuse after error: %d %d\n",
                rc, st.nconv);
        return 1;
    }
    puts("PASS callback failure and workspace reuse");
    return 0;
}

int main(void)
{
    int t, locking, scale, long_run, api, psd;
    for (t = 0; t <= 2; ++t)
        for (locking = 0; locking <= 1; ++locking)
            for (scale = 0; scale <= 1; ++scale)
                for (long_run = 0; long_run <= 1; ++long_run)
                    for (api = 0; api <= 2; ++api)
                        for (psd = 0; psd <= (t == LIS_KRYLOV_SHIFT_INVERT_B); ++psd)
                            if (run(t, locking, scale, long_run, api, psd)) return 1;
    for (api = 0; api <= 2; ++api)
        if (pole_case(api)) return 1;
    if (callback_failure_case()) return 1;
    puts("PASS verified Ritz residual reuse regression");
    return 0;
}
