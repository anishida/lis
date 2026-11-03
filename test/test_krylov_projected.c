/* Dense projected Rayleigh-Ritz regression through the public API.
 * A = H * diag(lambda) * H, H = I - 2*u*u^T, ||u|| = 1.
 * Tests large projected subspaces with known repeated/clustered modes.
 */
#include "lis_krylov.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 80
#define NEV 4

typedef struct {
    double u[N];
    double lambda[N];
} fixture;

static int apply_a(void *opaque, int n, const double *x, double *y)
{
    fixture *ctx = (fixture *)opaque;
    double dot = 0.0, other = 0.0;
    int i;
    if (n != N) return -1;
    for (i = 0; i < n; ++i) dot += ctx->u[i] * x[i];
    for (i = 0; i < n; ++i)
        y[i] = ctx->lambda[i] * (x[i] - 2.0 * ctx->u[i] * dot);
    for (i = 0; i < n; ++i) other += ctx->u[i] * y[i];
    for (i = 0; i < n; ++i) y[i] -= 2.0 * ctx->u[i] * other;
    return 0;
}

static void make_fixture(fixture *ctx, int clustered)
{
    int i;
    double square = 0.0;
    for (i = 0; i < N; ++i) {
        ctx->u[i] = sin((i + 1) * 0.43) + cos((i + 1) * 0.19);
        square += ctx->u[i] * ctx->u[i];
        ctx->lambda[i] = i < NEV ?
                         100.0 + (clustered ? i * 1.0e-6 : 0.0) :
                         85.0 - (i - NEV) * 0.35;
    }
    for (i = 0; i < N; ++i)
        ctx->u[i] /= sqrt(square);
}

static int run_case(int clustered, int hard_locking, int ncv)
{
    fixture ctx = {{0}, {0}};
    lis_krylov_options o;
    lis_krylov_problem p = {0};
    lis_krylov_stats stats = {0};
    double eigenvalues[NEV], eigenvectors[N*NEV], residuals[NEV];
    double a_x[N];
    int rc, i, j, k;
    make_fixture(&ctx, clustered);
    p.n = N;
    p.transform = LIS_KRYLOV_STANDARD;
    p.ctx = &ctx;
    p.apply_a = apply_a;
    lis_krylov_options_default(&o);
    o.nev = NEV;
    o.ncv = ncv;
    o.maxiter = 500;
    o.tol = 1.0e-9;
    o.which = "LA";
    o.hard_locking = hard_locking;
    rc = lis_krylov_symmetric(&p, &o, eigenvalues, eigenvectors,
                              residuals, &stats);
    if (rc != LIS_KRYLOV_SUCCESS || stats.nconv != NEV) {
        fprintf(stderr, "FAIL projected Ritz: rc=%d ncv=%d "
                "cluster=%d hardlock=%d nconv=%d restarts=%d maxres=%.3g\n",
                rc, ncv, clustered, hard_locking,
                stats.nconv, stats.restarts, stats.max_residual);
        return 1;
    }
    for (j = 0; j < NEV; ++j) {
        double norm = 0.0, rnorm = 0.0, anorm = 0.0;
        double expected = 100.0 + (clustered ? (NEV-1-j) * 1e-6 : 0.0);
        if (fabs(eigenvalues[j] - expected) > 2.0e-5 ||
            residuals[j] > o.tol) {
            fprintf(stderr, "FAIL projected eigenvalue j=%d "
                    "got=%.15g expected=%.15g residual=%.3g\n",
                    j, eigenvalues[j], expected, residuals[j]);
            return 1;
        }
        apply_a(&ctx, N, eigenvectors + (size_t)j*N, a_x);
        for (i = 0; i < N; ++i) {
            double x = eigenvectors[(size_t)j*N+i];
            double r = a_x[i] - eigenvalues[j]*x;
            norm += x*x;
            rnorm += r*r;
            anorm += a_x[i]*a_x[i];
        }
        if (fabs(norm-1.0) > 1.0e-6 ||
            sqrt(rnorm) > 2e-8*fmax(1.0, sqrt(anorm))) {
            fprintf(stderr,"FAIL projected residual j=%d "
                           "norm=%.15g residual=%.3g\n",
                    j, norm, sqrt(rnorm));
            return 1;
        }
        for (k = 0; k < j; ++k) {
            double dot = 0.0;
            for (i = 0; i < N; ++i)
                dot += eigenvectors[(size_t)k*N+i] *
                       eigenvectors[(size_t)j*N+i];
            if (fabs(dot) > 1e-6) {
                fprintf(stderr, "FAIL projected orthogonality %d/%d %g\n",
                        j, k, dot);
                return 1;
            }
        }
    }
    printf("PASS projected Ritz cluster=%d hardlock=%d ncv=%d "
           "restarts=%d maxres=%.3g\n", clustered, hard_locking,
           ncv, stats.restarts, stats.max_residual);
    return 0;
}

int main(void)
{
    return run_case(0, 0, 32) ||
           run_case(1, 0, 48) ||
           run_case(0, 1, 32);
}
