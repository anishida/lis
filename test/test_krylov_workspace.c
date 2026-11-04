/* Regression: storage-independent Krylov workspace reused across shifts,
 * callback failures and repeated solves, with one-shot equivalence.
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 48
#define NEV 3

typedef struct {
    double sigma;
    int solves;
    int fail_solve;
} context;

static int apply_a(void *opaque, int n, const double *x, double *y)
{
    int i;
    (void)opaque;
    if (n != N) return -1;
    for (i = 0; i < n; ++i)
        y[i] = (1.0 + 0.25*(double)i)*x[i];
    return 0;
}

static int apply_b(void *opaque, int n, const double *x, double *y)
{
    int i;
    (void)opaque;
    if (n != N) return -1;
    for (i = 0; i < n; ++i) y[i] = x[i];
    return 0;
}

static int inverse(void *opaque, int n, const double *x, double *y)
{
    context *c = (context *)opaque;
    int i;
    if (n != N || c->fail_solve) return -1;
    for (i = 0; i < n; ++i)
        y[i] = x[i] / (1.0 + 0.25*(double)i - c->sigma);
    ++c->solves;
    return 0;
}

static int run(int hard_locking)
{
    context c = {0};
    lis_krylov_problem p = {0};
    lis_krylov_options o;
    lis_krylov_workspace *workspace = NULL;
    lis_krylov_stats reuse_stats, fresh_stats;
    double reuse_values[NEV], reuse_vectors[N*NEV], reuse_residuals[NEV];
    double fresh_values[NEV], fresh_vectors[N*NEV], fresh_residuals[NEV];
    int i, j, k, rc, calls_before;

    p.n = N;
    p.transform = LIS_KRYLOV_SHIFT_INVERT_B;
    p.ctx = &c;
    p.apply_a = apply_a;
    p.apply_b = apply_b;
    p.solve = inverse;
    lis_krylov_options_default(&o);
    o.nev = NEV;
    o.ncv = 16;
    o.maxiter = 400;
    o.tol = 1e-8;
    o.hard_locking = hard_locking;

    if (lis_krylov_workspace_create(&p, &o, &workspace) ||
        !workspace) return 1;

    for (k = 0; k < 3; ++k) {
        c.sigma = 0.0 + 0.15*(double)k;
        p.sigma = c.sigma;
        calls_before = c.solves;
        rc = lis_krylov_symmetric_reuse(workspace, &p, &o,
               reuse_values, reuse_vectors, reuse_residuals, &reuse_stats);
        if (rc != LIS_KRYLOV_SUCCESS || reuse_stats.nconv != NEV ||
            c.solves - calls_before != reuse_stats.solve_count) {
            fprintf(stderr,"FAIL reusable shift=%g lock=%d rc=%d conv=%d\n",
                    c.sigma, hard_locking, rc, reuse_stats.nconv);
            lis_krylov_workspace_destroy(workspace);
            return 1;
        }

        rc = lis_krylov_symmetric(&p, &o, fresh_values, fresh_vectors,
                                  fresh_residuals, &fresh_stats);
        if (rc != LIS_KRYLOV_SUCCESS ||
            reuse_stats.nconv != fresh_stats.nconv ||
            reuse_stats.restarts != fresh_stats.restarts ||
            reuse_stats.solve_count != fresh_stats.solve_count) {
            fprintf(stderr, "FAIL stats parity lock=%d shift=%g\n",
                    hard_locking, c.sigma);
            lis_krylov_workspace_destroy(workspace);
            return 1;
        }
        for (j = 0; j < NEV; ++j) {
            double dot = 0.0;
            if (fabs(reuse_values[j] - fresh_values[j]) > 1e-9 ||
                reuse_residuals[j] > o.tol ||
                fabs(reuse_values[j] - (1.0 + 0.25*j)) > 1e-6) {
                fprintf(stderr,"FAIL eigenpair lock=%d shift=%g j=%d "
                               "lambda=%g res=%g\n", hard_locking, c.sigma,
                        j, reuse_values[j], reuse_residuals[j]);
                lis_krylov_workspace_destroy(workspace);
                return 1;
            }
            for (i = 0; i < N; ++i)
                dot += reuse_vectors[(size_t)j*N+i] *
                       fresh_vectors[(size_t)j*N+i];
            if (fabs(fabs(dot) - 1.0) > 1e-5) {
                fprintf(stderr, "FAIL vectors parity lock=%d j=%d dot=%g\n",
                        hard_locking, j, dot);
                lis_krylov_workspace_destroy(workspace);
                return 1;
            }
        }
    }

    /* Invalid dimensions must not modify or destroy the reusable buffers. */
    o.ncv = 15;
    if (lis_krylov_symmetric_reuse(workspace, &p, &o, reuse_values,
         reuse_vectors, reuse_residuals, &reuse_stats) !=
         LIS_KRYLOV_BAD_ARGUMENT) return 1;
    o.ncv = 16;

    c.fail_solve = 1;
    if (lis_krylov_symmetric_reuse(workspace, &p, &o, reuse_values,
          reuse_vectors, reuse_residuals, &reuse_stats) !=
          LIS_KRYLOV_CALLBACK_ERROR ||
        reuse_stats.nreturned != 0) return 1;
    c.fail_solve = 0;
    if (lis_krylov_symmetric_reuse(workspace, &p, &o, reuse_values,
          reuse_vectors, reuse_residuals, &reuse_stats) !=
          LIS_KRYLOV_SUCCESS) return 1;

    lis_krylov_workspace_destroy(workspace);
    printf("PASS reusable Krylov workspace hard_locking=%d "
           "different shifts and callback recovery\n", hard_locking);
    return 0;
}

int main(void)
{
    lis_krylov_workspace *workspace = (lis_krylov_workspace *)1;
    lis_krylov_options options;
    lis_krylov_problem problem = {0};
    int rc;
    lis_krylov_workspace_destroy(NULL);
    rc = lis_krylov_workspace_create(NULL, NULL, &workspace);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT || workspace) return 1;
    lis_krylov_options_default(&options);
    problem.n = N;
    options.nev = NEV;
    options.ncv = NEV;
    rc = lis_krylov_workspace_create(&problem, &options, &workspace);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT || workspace) return 1;
    return run(0) || run(1);
}
