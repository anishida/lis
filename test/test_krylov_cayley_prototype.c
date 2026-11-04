/* Stage 6: opt-in, matrix-free Cayley transform experiment for LIS Krylov.
 * Copyright (C) 2026, contributors to LIS.
 * Redistribution and use are governed by the LIS project BSD license.
 *
 * The B^(1/2) similarity below is fixture-specific (diagonal SPD B).
 * This test changes neither the production solver nor its public ABI.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#ifndef CAYLEY_ACTION_ONLY
#include "lis_krylov.h"
#endif

#define CY_NMAX 48
#define CY_NEV 3
#define CY_PI 3.14159265358979323846264338327950288

typedef struct {
    int n, a_calls, b_calls, solve_calls, cayley_calls, fail;
    double alpha, beta, sigma, tau;
    double d[CY_NMAX], pivot[CY_NMAX], upper[CY_NMAX];
} cy_problem;

#ifndef CAYLEY_ACTION_ONLY
static double cy_norm(int n, const double *x)
{
    int i;
    double s = 0.0;
    for (i = 0; i < n; ++i) s += x[i] * x[i];
    return sqrt(s);
}
#endif

/* Application-owned tridiagonal shifted factorization, prepared ONCE.
 * A = D T D, B = D^2, with T symmetric Toeplitz. The tested shifts
 * lie outside T's spectrum, so nonpivoted LU is suitable here only. */
static int cy_prepare(cy_problem *p, int n, double sigma, double tau)
{
    int i;
    double threshold;
    if (!p || n < 5 || n > CY_NMAX || !isfinite(sigma) ||
        !isfinite(tau) ||
        fabs(sigma - tau) <= 64.0 * DBL_EPSILON *
                             fmax(1.0, fmax(fabs(sigma), fabs(tau))))
        return -1;
    memset(p, 0, sizeof(*p));
    p->n = n;
    p->alpha = 4.0;
    p->beta = 0.25;
    p->sigma = sigma;
    p->tau = tau;
    for (i = 0; i < n; ++i) p->d[i] = 1.0 + 0.0125 * i;
    threshold = 64.0 * DBL_EPSILON *
                (fabs(p->alpha) + 2.0 * fabs(p->beta) + fabs(sigma) + 1.0) *
                p->d[n - 1] * p->d[n - 1];
    p->pivot[0] = (p->alpha - sigma) * p->d[0] * p->d[0];
    for (i = 0; i < n; ++i) {
        if (!isfinite(p->pivot[i]) || fabs(p->pivot[i]) <= threshold)
            return -2; /* singular or unsafe shift */
        if (i + 1 < n) {
            double off = p->beta * p->d[i] * p->d[i + 1];
            p->upper[i] = off / p->pivot[i];
            p->pivot[i + 1] = (p->alpha - sigma) *
                p->d[i + 1] * p->d[i + 1] - off * p->upper[i];
        }
    }
    return 0;
}

static int cy_apply_a(void *v, int n, const double *x, double *y)
{
    cy_problem *p = (cy_problem *)v;
    int i;
    if (!p || p->fail || n != p->n) return -1;
    ++p->a_calls;
    for (i = 0; i < n; ++i) {
        double z = p->alpha * p->d[i] * p->d[i] * x[i];
        if (i) z += p->beta * p->d[i - 1] * p->d[i] * x[i - 1];
        if (i + 1 < n) z += p->beta * p->d[i] * p->d[i + 1] * x[i + 1];
        y[i] = z;
    }
    return 0;
}

static int cy_apply_b(void *v, int n, const double *x, double *y)
{
    cy_problem *p = (cy_problem *)v;
    int i;
    if (!p || p->fail || n != p->n) return -1;
    ++p->b_calls;
    for (i = 0; i < n; ++i) y[i] = p->d[i] * p->d[i] * x[i];
    return 0;
}

static int cy_solve(void *v, int n, const double *rhs, double *x)
{
    cy_problem *p = (cy_problem *)v;
    double f[CY_NMAX];
    int i;
    if (!p || p->fail || n != p->n) return -1;
    ++p->solve_calls;
    f[0] = rhs[0] / p->pivot[0];
    for (i = 1; i < n; ++i) {
        double off = p->beta * p->d[i - 1] * p->d[i];
        f[i] = (rhs[i] - off * f[i - 1]) / p->pivot[i];
    }
    x[n - 1] = f[n - 1];
    for (i = n - 2; i >= 0; --i) x[i] = f[i] - p->upper[i] * x[i + 1];
    for (i = 0; i < n; ++i) if (!isfinite(x[i])) return -1;
    return 0;
}

/* Apply C_y = D (A-sigma B)^-1 (A-tau B) D^-1 to Euclidean y.
 * The original C_x is B-self-adjoint but generally NOT Euclidean
 * symmetric, so passing C_x directly as a STANDARD operator would
 * silently violate LIS's symmetric eigensolver contract. */
static int cy_apply_cayley(void *v, int n, const double *y, double *z)
{
    cy_problem *p = (cy_problem *)v;
    double x[CY_NMAX], ax[CY_NMAX], bx[CY_NMAX];
    double rhs[CY_NMAX] = {0}, solved[CY_NMAX] = {0};
    int i;
    if (!p || p->fail || n != p->n) return -1;
    ++p->cayley_calls;
    for (i = 0; i < n; ++i) x[i] = y[i] / p->d[i];
    if (cy_apply_a(v, n, x, ax) || cy_apply_b(v, n, x, bx)) return -1;
    for (i = 0; i < n; ++i) rhs[i] = ax[i] - p->tau * bx[i];
    if (cy_solve(v, n, rhs, solved)) return -1;
    for (i = 0; i < n; ++i) {
        z[i] = p->d[i] * solved[i];
        if (!isfinite(z[i])) return -1;
    }
    return 0;
}

static double cy_lambda(double sigma, double tau, double mu)
{
    if (!isfinite(mu) || fabs(mu - 1.0) <=
        128.0 * DBL_EPSILON * fmax(1.0, fabs(mu))) return NAN;
    return sigma + (sigma - tau) / (mu - 1.0);
}

/* Independent operator identities, analytic eigenvectors, symmetry,
 * and rejection of degenerate or singular shift pairs. */
static int cy_check_action(int n, double sigma, double tau)
{
    cy_problem p;
    double y[CY_NMAX], z[CY_NMAX], u[CY_NMAX], v[CY_NMAX];
    double cu[CY_NMAX], cv[CY_NMAX], bx[CY_NMAX], x[CY_NMAX];
    double shifted[CY_NMAX], alternative[CY_NMAX];
    int i, j;
    if (cy_prepare(&p, n, sigma, tau)) return 1;
    for (i = 0; i < n; ++i) {
        u[i] = sin(0.31 * (i + 1)) + 0.01 * i;
        v[i] = cos(0.27 * (i + 1)) - 0.01 * i;
    }
    if (cy_apply_cayley(&p, n, u, cu) ||
        cy_apply_cayley(&p, n, v, cv)) return 2;
    {
        double uv = 0.0, vu = 0.0, scale = 1.0;
        for (i = 0; i < n; ++i) {
            uv += u[i] * cv[i];
            vu += cu[i] * v[i];
            scale += fabs(u[i] * cv[i]) + fabs(cu[i] * v[i]);
        }
        if (fabs(uv - vu) > 2e-12 * scale) return 3;
    }
    for (i = 0; i < n; ++i) x[i] = u[i] / p.d[i];
    if (cy_apply_b(&p, n, x, bx) || cy_solve(&p, n, bx, shifted))
        return 4;
    for (i = 0; i < n; ++i)
        alternative[i] = u[i] + (sigma - tau) * p.d[i] * shifted[i];
    for (i = 0; i < n; ++i)
        if (fabs(alternative[i] - cu[i]) >
            8e-13 * (1.0 + fabs(cu[i]))) return 5;
    for (j = 1; j <= n; j += n - 1) {
        double lam = p.alpha + 2.0 * p.beta * cos(j * CY_PI / (n + 1));
        double mu = (lam - tau) / (lam - sigma);
        double err = 0.0;
        for (i = 0; i < n; ++i)
            y[i] = sqrt(2.0 / (n + 1)) *
                   sin((i + 1) * j * CY_PI / (n + 1));
        if (cy_apply_cayley(&p, n, y, z)) return 6;
        for (i = 0; i < n; ++i)
            err = fmax(err, fabs(z[i] - mu * y[i]));
        if (err > 3e-12 * (1.0 + fabs(mu)) ||
            fabs(cy_lambda(sigma, tau, mu) - lam) > 2e-12)
            return 7;
    }
    if (cy_prepare(&p, n, sigma, sigma) != -1) return 8;
    if (cy_prepare(&p, n, sigma, NAN) != -1) return 9;
    if (cy_prepare(&p, 41, 4.0, 3.2) != -2) return 10;
    if (cy_prepare(&p, n, sigma, tau)) return 11;
    p.fail = 1;
    if (!cy_apply_cayley(&p, n, u, cu)) return 12;
    return 0;
}

#ifndef CAYLEY_ACTION_ONLY
/* Verify original pencil residual and B-orthogonality, independently
 * of the transformed Ritz residual reported by the native solver. */
static int cy_check_pairs(cy_problem *p, const double *lambda,
                          const double *vectors, int is_euclidean,
                          int from_upper_edge)
{
    double x[CY_NMAX], ax[CY_NMAX], bx[CY_NMAX], r[CY_NMAX];
    int k, i, j;
    for (k = 0; k < CY_NEV; ++k) {
        double expected, relative;
        int index = from_upper_edge ? k + 1 : p->n - k;
        expected = p->alpha + 2.0 * p->beta *
                   cos(index * CY_PI / (p->n + 1));
        if (!isfinite(lambda[k]) || fabs(lambda[k] - expected) > 2e-5)
            return 20 + k;
        for (i = 0; i < p->n; ++i)
            x[i] = vectors[(size_t)k * p->n + i] /
                   (is_euclidean ? p->d[i] : 1.0);
        if (cy_apply_a(p, p->n, x, ax) || cy_apply_b(p, p->n, x, bx))
            return 30 + k;
        for (i = 0; i < p->n; ++i) r[i] = ax[i] - lambda[k] * bx[i];
        relative = cy_norm(p->n, r) /
                   fmax(DBL_MIN, cy_norm(p->n, ax) +
                        fabs(lambda[k]) * cy_norm(p->n, bx));
        if (!isfinite(relative) || relative > 8e-8) return 40 + k;
        for (j = 0; j <= k; ++j) {
            double dot = 0.0;
            for (i = 0; i < p->n; ++i) {
                double vi = vectors[(size_t)j * p->n + i];
                double vk = vectors[(size_t)k * p->n + i];
                dot += is_euclidean ? vi * vk :
                       vi * vk * p->d[i] * p->d[i];
            }
            if (fabs(dot - (j == k ? 1.0 : 0.0)) > 3e-6)
                return 50 + k;
        }
    }
    return 0;
}

static int cy_run_native(int n, double sigma, double tau)
{
    cy_problem p;
    lis_krylov_problem problem;
    lis_krylov_options options;
    lis_krylov_stats stats;
    double mu[CY_NEV], lambda[CY_NEV], residual[CY_NEV];
    double native_lambda[CY_NEV], native_res[CY_NEV];
    double v[CY_NMAX * CY_NEV], w[CY_NMAX * CY_NEV];
    int i, rc, expected_upper = sigma > 4.5;
    if (cy_prepare(&p, n, sigma, tau)) return 60;
    memset(&problem, 0, sizeof problem);
    problem.n = n;
    problem.ctx = &p;
    problem.transform = LIS_KRYLOV_STANDARD;
    problem.apply_a = cy_apply_cayley;
    lis_krylov_options_default(&options);
    options.nev = CY_NEV;
    options.ncv = 24;
    options.maxiter = 700;
    options.tol = 1e-8;
    options.which = "LM";
    options.residual_scale = 1.0;
    rc = lis_krylov_symmetric(&problem, &options,
                              mu, v, residual, &stats);
    if (rc != LIS_KRYLOV_SUCCESS || stats.nconv != CY_NEV) {
        fprintf(stderr, "Cayley Krylov failed rc=%d nconv=%d\n",
                rc, stats.nconv);
        return 61;
    }
    for (i = 0; i < CY_NEV; ++i) {
        lambda[i] = cy_lambda(sigma, tau, mu[i]);
        if (!isfinite(lambda[i])) return 62;
    }
    rc = cy_check_pairs(&p, lambda, v, 1, expected_upper);
    if (rc) return rc;
    printf("CAYLEY n=%d sigma=%.3f tau=%.3f op=%d solves=%d residual=%g\n",
           n, sigma, tau, stats.op_count, p.solve_calls,
           stats.max_residual);

    /* Same original pencil via existing native SHIFT_INVERT_B backend. */
    problem.transform = LIS_KRYLOV_SHIFT_INVERT_B;
    problem.sigma = sigma;
    problem.apply_a = cy_apply_a;
    problem.apply_b = cy_apply_b;
    problem.solve = cy_solve;
    rc = lis_krylov_symmetric(&problem, &options,
                              native_lambda, w, native_res, &stats);
    if (rc != LIS_KRYLOV_SUCCESS || stats.nconv != CY_NEV) {
        fprintf(stderr, "Shift-invert reference failed rc=%d nconv=%d\n",
                rc, stats.nconv);
        return 63;
    }
    rc = cy_check_pairs(&p, native_lambda, w, 0, expected_upper);
    if (rc) return rc;
    for (i = 0; i < CY_NEV; ++i)
        if (fabs(native_lambda[i] - lambda[i]) > 3e-6) return 64;
    printf("SHIFT_INVERT_B n=%d sigma=%.3f op=%d solves=%d residual=%g\n",
           n, sigma, stats.op_count, stats.solve_count,
           stats.max_residual);
    return 0;
}
#endif

int main(void)
{
    int rc;
    rc = cy_check_action(40, 4.65, 3.20);
    if (!rc) rc = cy_check_action(40, 3.35, 4.80);
#ifndef CAYLEY_ACTION_ONLY
    if (!rc) rc = cy_run_native(40, 4.65, 3.20);
    if (!rc) rc = cy_run_native(40, 3.35, 4.80);
#endif
    if (rc) {
        fprintf(stderr, "FAIL Cayley stage 6: %d\n", rc);
        return 1;
    }
    puts("PASS Cayley transform stage 6");
    return 0;
}
