/* Native thick-restarted symmetric Rayleigh-Ritz subspace eigensolver.
 * Copyright (C) 2026, contributors to LIS.
 * Redistribution and use are governed by the LIS project BSD license.
 */
#include "lis_krylov.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* A direction is considered linearly dependent only when double-pass
 * orthogonalization removes essentially all of its input metric norm. */
#define KRYLOV_ORTH_REJECT 1.0e-12
#define KRYLOV_DGKS_ETA 0.717
#define KRYLOV_DIRECTION_REJECT (128.0 * DBL_EPSILON)
#define KRYLOV_MAX_SWEEPS 80
#define KRYLOV_JACOBI_EPS (16.0 * DBL_EPSILON)

struct krylov_workspace {
    const lis_krylov_problem *problem;
    const lis_krylov_options *options;
    lis_krylov_stats *stats;
    int n;
    int ncv;
    double norm_a;
    double *basis;
    double *metric_basis;
    double *images;
    double *projected;
    double *dense_work;
    double *ritz;
    double *eval;
    double *scratch;
    double *x;
    double *metric_x;
    double *rhs;
    double *check;
    int *order;
    int *locked;
};

static double compensated_dot(int n, const double *a, const double *b)
{
    double sum = 0.0;
    double compensation = 0.0;
    int i;
    for (i = 0; i < n; ++i) {
        double value = a[i] * b[i] - compensation;
        double next = sum + value;
        compensation = (next - sum) - value;
        sum = next;
    }
    return sum;
}

static double stable_norm(int n, const double *x)
{
    double scale = 0.0;
    double square_sum = 1.0;
    int i;
    for (i = 0; i < n; ++i) {
        double a = fabs(x[i]);
        if (!isfinite(a)) return INFINITY;
        if (a == 0.0) continue;
        if (scale < a) {
            double ratio = scale / a;
            square_sum = 1.0 + square_sum * ratio * ratio;
            scale = a;
        } else {
            double ratio = a / scale;
            square_sum += ratio * ratio;
        }
    }
    return scale == 0.0 ? 0.0 : scale * sqrt(square_sum);
}

static int all_finite(int n, const double *x)
{
    int i;
    for (i = 0; i < n; ++i)
        if (!isfinite(x[i])) return 0;
    return 1;
}

static int apply_metric(struct krylov_workspace *s,
                        const double *x, double *y)
{
    const lis_krylov_problem *p = s->problem;
    lis_krylov_matvec fn;
    if (p->transform == LIS_KRYLOV_STANDARD) {
        memcpy(y, x, (size_t)s->n * sizeof(double));
    } else {
        fn = p->transform == LIS_KRYLOV_SHIFT_INVERT_A ?
             p->apply_a : p->apply_b;
        ++s->stats->metric_count;
        if (fn(p->ctx, s->n, x, y)) return -1;
    }
    return all_finite(s->n, y) ? 0 : -2;
}

/* Validate an approximate shifted inverse only when explicitly requested.
 * This costs two extra matrix-vector products per solve. */
static int verify_inverse(struct krylov_workspace *s,
                          const double *rhs, const double *solution)
{
    const lis_krylov_problem *p = s->problem;
    double numerator, denominator;
    int i;

    if (s->options->inner_tol == 0.0) return 0;
    if (p->apply_a(p->ctx, s->n, solution, s->check)) return -1;
    if (p->apply_b(p->ctx, s->n, solution, s->rhs)) return -1;
    for (i = 0; i < s->n; ++i)
        s->check[i] -= p->sigma * s->rhs[i] + rhs[i];
    numerator = stable_norm(s->n, s->check);
    denominator = stable_norm(s->n, rhs);
    if (!isfinite(numerator) || !isfinite(denominator)) return -2;
    if (numerator > s->options->inner_tol * denominator)
        return -3;
    return 0;
}

/* The metric image is already computed during orthogonalization. Reusing
 * it avoids an unnecessary matrix-vector product per Krylov extension. */
static int apply_operator(struct krylov_workspace *s,
                          const double *x, const double *metric_x,
                          double *y)
{
    const lis_krylov_problem *p = s->problem;
    int err;
    ++s->stats->op_count;
    if (p->transform == LIS_KRYLOV_STANDARD) {
        err = p->apply_a(p->ctx, s->n, x, y);
    } else {
        ++s->stats->solve_count;
        err = p->solve(p->ctx, s->n, metric_x, y);
        if (!err && all_finite(s->n, y))
            err = verify_inverse(s, metric_x, y);
    }
    if (err) return err == -3 ? -3 : -1;
    return all_finite(s->n, y) ? 0 : -2;
}

static int orthogonalize(struct krylov_workspace *s,
                         int count, double *x, double *metric_x)
{
    double initial, square, norm;
    int pass, j, i, err;

    norm = stable_norm(s->n, x);
    if (!isfinite(norm)) return -2;
    if (norm == 0.0) return 1;
    for (i = 0; i < s->n; ++i) x[i] /= norm;

    err = apply_metric(s, x, metric_x);
    if (err) return err;
    initial = compensated_dot(s->n, x, metric_x);
    if (!isfinite(initial)) return -2;
    if (initial <= 0.0) {
        if (s->options->allow_semidefinite_metric && initial == 0.0)
            return 1;
        return -2;
    }

    /* DGKS selective reorthogonalization in the chosen metric.  The
     * norm-loss test compares squared metric norms to avoid square roots.
     * Recompute B*x after each pass; reusing a heavily cancelled metric
     * image can give an inaccurate norm on semidefinite B. */
    {
        double previous_square = initial;
        for (pass = 0; pass < 3; ++pass) {
            for (j = 0; j < count; ++j) {
                const double *v = s->basis + (size_t)j * s->n;
                const double *bv = s->metric_basis + (size_t)j * s->n;
                double coefficient = compensated_dot(s->n, bv, x);
                if (!isfinite(coefficient)) return -2;
                for (i = 0; i < s->n; ++i)
                    x[i] -= coefficient * v[i];
            }
            err = apply_metric(s, x, metric_x);
            if (err) return err;
            square = compensated_dot(s->n, x, metric_x);
            if (!isfinite(square)) return -2;
            /* Tiny negative roundoff is dependence only for a PSD B. */
            if (square <= initial * KRYLOV_ORTH_REJECT *
                          KRYLOV_ORTH_REJECT) {
                if (square >= 0.0 ||
                    (s->options->allow_semidefinite_metric &&
                     square >= -64.0 * DBL_EPSILON * initial))
                    return 1;
                return -2;
            }
            if (square >= KRYLOV_DGKS_ETA * KRYLOV_DGKS_ETA *
                          previous_square)
                break;
            previous_square = square;
        }
    }
    {
        double inverse_norm = 1.0 / sqrt(square);
        for (i = 0; i < s->n; ++i) {
            x[i] *= inverse_norm;
            metric_x[i] *= inverse_norm;
        }
    }
    for (j = 0; j < count; ++j) {
        const double *bv = s->metric_basis + (size_t)j * s->n;
        if (fabs(compensated_dot(s->n, bv, x)) > 1.0e-8)
            return 1;
    }
    return 0;
}

static int append_vector(struct krylov_workspace *s, int count, double *x);

/* For a positive-semidefinite metric B, the generalized problem has
 * finite modes only in the nonzero part of the transformed spectrum.
 * Start and replenishment directions in range((A-sigma B)^-1 B) so that
 * components of ker(B) are not introduced as arbitrary Ritz directions.
 * This still requires a nonsingular shifted operator and a nonnegative
 * B-inner product on the generated Krylov subspace. */
static int append_range_seed(struct krylov_workspace *s,
                             int count, double *x)
{
    int err;
    if (s->options->allow_semidefinite_metric) {
        err = apply_metric(s, x, s->metric_x);
        if (err) return err;
        if (stable_norm(s->n, s->metric_x) == 0.0) return 1;
        /* scratch has at least n elements, and does not overlap the
         * scratch vectors used by optional inner-solve certification. */
        err = apply_operator(s, x, s->metric_x, s->scratch);
        if (err) return err;
        memcpy(x, s->scratch, (size_t)s->n * sizeof(double));
    }
    return append_vector(s, count, x);
}

static int append_vector(struct krylov_workspace *s,
                         int count, double *x)
{
    int err = orthogonalize(s, count, x, s->metric_x);
    if (err) return err;
    memcpy(s->basis + (size_t)count * s->n, x,
           (size_t)s->n * sizeof(double));
    if (s->metric_basis != s->basis)
        memcpy(s->metric_basis + (size_t)count * s->n, s->metric_x,
               (size_t)s->n * sizeof(double));
    return apply_operator(s, x, s->metric_x,
                          s->images + (size_t)count * s->n);
}

static void seed_vector(int n, double *x, unsigned seed)
{
    unsigned state = seed ^ 0x6d2b79f5u;
    int i;
    for (i = 0; i < n; ++i) {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        x[i] = (double)(state & 0xffffu) / 32768.0 - 1.0;
    }
}

/* Dense cyclic Jacobi eigendecomposition of symmetric projected matrix.
 * The transformation is stable for modest restarted subspace sizes. */
static int diagonalize(int m, double *a, double *q,
                        double *eigenvalues, int stride)
{
    int sweep, p, r, i;
    for (i = 0; i < m; ++i)
        for (p = 0; p < m; ++p)
            q[(size_t)i * stride + p] = (i == p);

    for (sweep = 0; sweep < KRYLOV_MAX_SWEEPS; ++sweep) {
        double largest = 0.0;
        double diagonal_norm = 0.0;
        for (p = 0; p < m; ++p) {
            diagonal_norm += fabs(a[(size_t)p * stride + p]);
            for (r = p + 1; r < m; ++r) {
                double off = fabs(a[(size_t)p * stride + r]);
                if (off > largest) largest = off;
            }
        }
        if (largest <= KRYLOV_JACOBI_EPS *
                       (diagonal_norm + DBL_MIN)) break;

        for (p = 0; p < m; ++p) {
            for (r = p + 1; r < m; ++r) {
                double apr = a[(size_t)p * stride + r];
                double app = a[(size_t)p * stride + p];
                double arr = a[(size_t)r * stride + r];
                double tau, angle, cosine, sine;
                if (fabs(apr) <= DBL_EPSILON *
                                (fabs(app) + fabs(arr) + DBL_MIN))
                    continue;
                tau = (arr - app) / (2.0 * apr);
                angle = copysign(1.0, tau) /
                        (fabs(tau) + hypot(1.0, tau));
                cosine = 1.0 / sqrt(1.0 + angle * angle);
                sine = angle * cosine;
                a[(size_t)p * stride + p] = app - angle * apr;
                a[(size_t)r * stride + r] = arr + angle * apr;
                a[(size_t)p * stride + r] = 0.0;
                a[(size_t)r * stride + p] = 0.0;
                for (i = 0; i < m; ++i) {
                    double vip, vir;
                    if (i != p && i != r) {
                        vip = a[(size_t)i * stride + p];
                        vir = a[(size_t)i * stride + r];
                        a[(size_t)i * stride + p] =
                            a[(size_t)p * stride + i] =
                            cosine * vip - sine * vir;
                        a[(size_t)i * stride + r] =
                            a[(size_t)r * stride + i] =
                            sine * vip + cosine * vir;
                    }
                    vip = q[(size_t)i * stride + p];
                    vir = q[(size_t)i * stride + r];
                    q[(size_t)i * stride + p] =
                        cosine * vip - sine * vir;
                    q[(size_t)i * stride + r] =
                        sine * vip + cosine * vir;
                }
            }
        }
    }
    if (sweep == KRYLOV_MAX_SWEEPS) return -1;
    for (i = 0; i < m; ++i)
        eigenvalues[i] = a[(size_t)i * stride + i];
    return 0;
}

static int earlier(double a, double b, const char *which)
{
    if (which[0] == 'L' && which[1] == 'M') return fabs(a) > fabs(b);
    if (which[0] == 'S' && which[1] == 'M') return fabs(a) < fabs(b);
    if (which[0] == 'L' && which[1] == 'A') return a > b;
    return a < b;
}

static void sort_ritz(int m, const double *values, int *order,
                      const char *which)
{
    int i, j;
    for (i = 0; i < m; ++i) order[i] = i;
    for (i = 1; i < m; ++i) {
        int key = order[i];
        j = i - 1;
        while (j >= 0 && earlier(values[key], values[order[j]], which)) {
            order[j + 1] = order[j];
            --j;
        }
        order[j + 1] = key;
    }
}

static void rotate_column(struct krylov_workspace *s,
                          const double *basis, int m,
                          int idx, double *output)
{
    int i, j;
    for (i = 0; i < s->n; ++i) output[i] = 0.0;
    for (j = 0; j < m; ++j) {
        double c = s->ritz[(size_t)j * s->ncv + idx];
        const double *col = basis + (size_t)j * s->n;
        for (i = 0; i < s->n; ++i) output[i] += c * col[i];
    }
}

static int eigenpair_residual(struct krylov_workspace *s,
                              int m, int index, double mu,
                              double *lambda, double *relative)
{
    const lis_krylov_problem *p = s->problem;
    double norm_a, norm_b, norm_r, norm_x, denom;
    int i;

    rotate_column(s, s->basis, m, index, s->x);
    if (p->transform == LIS_KRYLOV_STANDARD) {
        *lambda = mu;
    } else if (p->transform == LIS_KRYLOV_SHIFT_INVERT_B) {
        if (fabs(mu) < DBL_MIN) return 1;
        *lambda = p->sigma + 1.0 / mu;
    } else {
        /* A Ritz pole too close to 1 cannot be mapped reliably:
         * the corresponding generalized value may be infinite. */
        if (fabs(mu - 1.0) <= 64.0 * DBL_EPSILON *
                                fmax(1.0, fabs(mu))) return 1;
        *lambda = p->sigma * mu / (mu - 1.0);
    }
    if (!isfinite(*lambda)) return 1;

    if (p->apply_a(p->ctx, s->n, s->x, s->metric_x)) return -1;
    norm_a = stable_norm(s->n, s->metric_x);
    if (p->transform == LIS_KRYLOV_STANDARD) {
        memcpy(s->rhs, s->x, (size_t)s->n * sizeof(double));
    } else if (p->apply_b(p->ctx, s->n, s->x, s->rhs)) {
        return -1;
    }
    norm_b = stable_norm(s->n, s->rhs);
    norm_x = stable_norm(s->n, s->x);
    if (!isfinite(norm_a) || !isfinite(norm_b) ||
        !isfinite(norm_x)) return -2;

    if (s->options->residual_scale == 0.0 && norm_x > 0.0 &&
        norm_a / norm_x > s->norm_a)
        s->norm_a = norm_a / norm_x;

    for (i = 0; i < s->n; ++i)
        s->metric_x[i] -= *lambda * s->rhs[i];
    norm_r = stable_norm(s->n, s->metric_x);
    denom = norm_a + fabs(*lambda) * norm_b;
    if (s->norm_a * norm_x > denom)
        denom = s->norm_a * norm_x;
    if (!isfinite(denom) || !isfinite(norm_r)) return -2;
    *relative = denom > 0.0 ? norm_r / denom : norm_r;
    return isfinite(*relative) ? 0 : -2;
}

/* Update the projected matrix only for new columns. Projection entries
 * for retained Ritz vectors are restored on restart from the Ritz values. */
static void update_projected(struct krylov_workspace *s,
                             int previous, int m)
{
    int i, j;
    for (j = previous; j < m; ++j) {
        for (i = 0; i <= j; ++i) {
            double left = compensated_dot(s->n,
                               s->metric_basis + (size_t)i * s->n,
                               s->images + (size_t)j * s->n);
            double right = compensated_dot(s->n,
                               s->metric_basis + (size_t)j * s->n,
                               s->images + (size_t)i * s->n);
            double a = 0.5 * (left + right);
            s->projected[(size_t)i * s->ncv + j] = a;
            s->projected[(size_t)j * s->ncv + i] = a;
        }
    }
}

static void restart_basis(struct krylov_workspace *s,
                          int m, int retain)
{
    double *sets[3];
    int count = 0, k, j, i;
    sets[count++] = s->basis;
    if (s->metric_basis != s->basis)
        sets[count++] = s->metric_basis;
    sets[count++] = s->images;
    for (k = 0; k < count; ++k) {
        for (j = 0; j < retain; ++j)
            rotate_column(s, sets[k], m, s->order[j],
                          s->scratch + (size_t)j * s->n);
        for (j = 0; j < retain; ++j)
            for (i = 0; i < s->n; ++i)
                sets[k][(size_t)j * s->n + i] =
                    s->scratch[(size_t)j * s->n + i];
    }
    for (i = 0; i < retain; ++i) {
        for (j = 0; j < retain; ++j)
            s->projected[(size_t)i * s->ncv + j] =
                (i == j ? s->eval[s->order[i]] : 0.0);
    }
}

void lis_krylov_options_default(lis_krylov_options *o)
{
    if (!o) return;
    memset(o, 0, sizeof(*o));
    o->tol = 1.0e-10;
    o->maxiter = 300;
    o->which = "LM";
}

static void release_workspace(struct krylov_workspace *s)
{
    free(s->basis);
    if (s->metric_basis != s->basis) free(s->metric_basis);
    free(s->images);
    free(s->projected);
    free(s->dense_work);
    free(s->ritz);
    free(s->eval);
    free(s->scratch);
    free(s->x);
    free(s->metric_x);
    free(s->rhs);
    free(s->check);
    free(s->order);
    free(s->locked);
}

int lis_krylov_symmetric(const lis_krylov_problem *p,
                         const lis_krylov_options *o,
                         double *values, double *vectors, double *residuals,
                         lis_krylov_stats *stats)
{
    struct krylov_workspace s = {0};
    size_t basis_size, dense_size, scratch_size;
    int m = 0, previous = 0, iter, ncv, nev, maxiter, rc;
    int i, j;
    double tol;
    unsigned seed;
    const char *which;

    if (stats) memset(stats, 0, sizeof(*stats));
    if (!p || !o || !values || !vectors || !residuals || !stats ||
        p->n < 3 || o->nev < 1 || o->nev >= p->n || !p->apply_a ||
        p->transform < LIS_KRYLOV_STANDARD ||
        p->transform > LIS_KRYLOV_SHIFT_INVERT_A ||
        (p->transform != LIS_KRYLOV_STANDARD &&
         (!p->apply_b || !p->solve)) ||
        (p->transform == LIS_KRYLOV_SHIFT_INVERT_A && p->sigma == 0.0) ||
        !isfinite(p->sigma) || !isfinite(o->tol) || o->tol < 0.0 ||
        (o->allow_semidefinite_metric != 0 &&
         o->allow_semidefinite_metric != 1) ||
        (o->allow_semidefinite_metric &&
         p->transform != LIS_KRYLOV_SHIFT_INVERT_B) ||
        !isfinite(o->inner_tol) || o->inner_tol < 0.0 ||
        !isfinite(o->residual_scale) || o->residual_scale < 0.0 ||
        o->maxiter < 0)
        return LIS_KRYLOV_BAD_ARGUMENT;
    which = o->which ? o->which : "LM";
    if (strcmp(which, "LM") && strcmp(which, "SM") &&
        strcmp(which, "LA") && strcmp(which, "SA"))
        return LIS_KRYLOV_BAD_ARGUMENT;

    nev = o->nev;
    maxiter = o->maxiter ? o->maxiter : 300;
    tol = o->tol ? o->tol : 1.0e-10;
    ncv = o->ncv ? o->ncv :
          (nev < (p->n - 8) / 3 ? 3 * nev + 8 : p->n);
    if (ncv <= nev || ncv > p->n || p->n > INT_MAX / 8 ||
        ncv > INT_MAX / 8) return LIS_KRYLOV_BAD_ARGUMENT;
    basis_size = (size_t)p->n * (size_t)ncv;
    dense_size = (size_t)ncv * (size_t)ncv;
    scratch_size = (size_t)p->n * (size_t)
                   (nev > ncv / 2 ? nev : ncv / 2);
    if (basis_size > SIZE_MAX / sizeof(double) ||
        dense_size > SIZE_MAX / sizeof(double) ||
        scratch_size > SIZE_MAX / sizeof(double))
        return LIS_KRYLOV_BAD_ARGUMENT;

    s.n = p->n;
    s.ncv = ncv;
    s.problem = p;
    s.options = o;
    s.stats = stats;
    s.norm_a = o->residual_scale;
    seed = o->random_seed ? o->random_seed : 1u;
    rc = LIS_KRYLOV_NO_CONVERGENCE;
#define ALLOC(PTR, COUNT, TYPE)                                      \
    do {                                                             \
        (PTR) = (TYPE *)calloc((COUNT), sizeof(TYPE));               \
        if (!(PTR)) { rc = LIS_KRYLOV_NO_MEMORY; goto done; }        \
    } while (0)
    ALLOC(s.basis, basis_size, double);
    if (p->transform == LIS_KRYLOV_STANDARD) {
        s.metric_basis = s.basis;
    } else {
        ALLOC(s.metric_basis, basis_size, double);
    }
    ALLOC(s.images, basis_size, double);
    ALLOC(s.projected, dense_size, double);
    ALLOC(s.dense_work, dense_size, double);
    ALLOC(s.ritz, dense_size, double);
    ALLOC(s.eval, ncv, double);
    ALLOC(s.scratch, scratch_size, double);
    ALLOC(s.x, s.n, double);
    ALLOC(s.metric_x, s.n, double);
    ALLOC(s.rhs, s.n, double);
    ALLOC(s.check, s.n, double);
    ALLOC(s.order, ncv, int);
    ALLOC(s.locked, ncv, int);
#undef ALLOC

    if (o->initial_vector) {
        memcpy(s.x, o->initial_vector, (size_t)s.n * sizeof(double));
    } else {
        seed_vector(s.n, s.x, seed);
    }
    if (o->residual_scale == 0.0) {
        double norm_x = stable_norm(s.n, s.x);
        double norm_a;
        if (!isfinite(norm_x) || norm_x == 0.0) {
            rc = LIS_KRYLOV_BAD_ARGUMENT;
            goto done;
        }
        if (p->apply_a(p->ctx, s.n, s.x, s.rhs)) {
            rc = LIS_KRYLOV_CALLBACK_ERROR;
            goto done;
        }
        norm_a = stable_norm(s.n, s.rhs);
        if (!isfinite(norm_a)) {
            rc = LIS_KRYLOV_BREAKDOWN;
            goto done;
        }
        s.norm_a = norm_a / norm_x;
    }

    /* In semidefinite mode an explicit start may lie in ker(B).  Try
     * reproducible range-projected replacement directions instead of
     * reporting spurious breakdown or dividing by a zero metric norm. */
    {
        int attempt;
        for (attempt = 0; attempt <
             (o->allow_semidefinite_metric ? 16 : 1); ++attempt) {
            if (attempt > 0)
                seed_vector(s.n, s.x, seed + (unsigned)attempt);
            i = append_range_seed(&s, m, s.x);
            if (i == 0) break;
            if (i < 0) {
                rc = i == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                     i == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                     LIS_KRYLOV_BREAKDOWN;
                goto done;
            }
        }
        if (i != 0) {
            rc = o->allow_semidefinite_metric ?
                 LIS_KRYLOV_NO_CONVERGENCE : LIS_KRYLOV_BREAKDOWN;
            goto done;
        }
    }
    m = 1;

    for (iter = 0; iter < maxiter; ++iter) {
        int old_m = m, extensions = 0;
        stats->iterations = iter + 1;
        update_projected(&s, previous, m);
        memcpy(s.dense_work, s.projected, dense_size * sizeof(double));
        if (diagonalize(m, s.dense_work, s.ritz, s.eval, ncv)) {
            rc = LIS_KRYLOV_BREAKDOWN;
            goto done;
        }
        sort_ritz(m, s.eval, s.order, which);
        stats->nreturned = 0;
        stats->nconv = 0;
        stats->max_residual = 0.0;
        memset(s.locked, 0, (size_t)m * sizeof(int));
        for (j = 0; j < m && stats->nreturned < nev; ++j) {
            int idx = s.order[j];
            double lambda = NAN, rel = INFINITY;
            int err = eigenpair_residual(&s, m, idx, s.eval[idx],
                                         &lambda, &rel);
            if (err < 0) {
                rc = err == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                     LIS_KRYLOV_BREAKDOWN;
                goto done;
            }
            /* Do not expose a nonrepresentable original eigenvalue
             * as a valid partial result. */
            if (err == 1) continue;
            s.locked[j] = (rel <= tol);
            if (s.locked[j]) ++stats->nconv;
            if (rel > stats->max_residual) stats->max_residual = rel;
            values[stats->nreturned] = lambda;
            residuals[stats->nreturned] = rel;
            memcpy(vectors + (size_t)stats->nreturned * s.n, s.x,
                   (size_t)s.n * sizeof(double));
            ++stats->nreturned;
        }
        if (stats->nconv == nev) {
            rc = LIS_KRYLOV_SUCCESS;
            goto done;
        }
        if (m >= ncv) {
            int keep = m / 2;
            if (keep < nev) keep = nev;
            if (keep >= m) keep = m - 1;
            restart_basis(&s, m, keep);
            ++stats->restarts;
            m = keep;
            previous = keep;
            /* A degeneracy cannot be resolved by residual expansion alone
             * when the retained subspace already approximates an invariant
             * eigenspace. Inject a reproducible independent search direction
             * after restarting to discover other vectors of a cluster. */
            for (i = 0; i < 4 && m < ncv; ++i) {
                int result;
                seed_vector(s.n, s.x,
                            seed + (unsigned)iter * 997u + (unsigned)i + 17u);
                result = append_range_seed(&s, m, s.x);
                if (result < 0) {
                    rc = result == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                         result == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                         LIS_KRYLOV_BREAKDOWN;
                    goto done;
                }
                if (result == 0) { ++m; break; }
            }
            continue;
        }

        for (j = 0; j < old_m && j < nev && m < ncv; ++j) {
            int idx = s.order[j];
            double mu = s.eval[idx];
            double norm_w, norm_v, norm_r, threshold;
            if (s.locked[j]) continue;
            rotate_column(&s, s.images, old_m, idx, s.x);
            rotate_column(&s, s.basis, old_m, idx, s.metric_x);
            norm_w = stable_norm(s.n, s.x);
            norm_v = stable_norm(s.n, s.metric_x);
            for (i = 0; i < s.n; ++i)
                s.x[i] -= mu * s.metric_x[i];
            norm_r = stable_norm(s.n, s.x);
            threshold = KRYLOV_DIRECTION_REJECT *
                        (norm_w + fabs(mu) * norm_v);
            if (norm_r <= threshold) continue;
            i = append_vector(&s, m, s.x);
            if (i < 0) {
                rc = i == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                     i == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                     LIS_KRYLOV_BREAKDOWN;
                goto done;
            }
            if (i == 0) { ++m; ++extensions; }
        }
        previous = old_m;
        if (!extensions && m < ncv) {
            int trial;
            for (trial = 0; trial < 6 && !extensions; ++trial) {
                seed_vector(s.n, s.x,
                            seed + (unsigned)iter * 7u + (unsigned)trial + 1u);
                i = append_range_seed(&s, m, s.x);
                if (i < 0) {
                    rc = i == -1 ? LIS_KRYLOV_CALLBACK_ERROR :
                         i == -3 ? LIS_KRYLOV_INEXACT_SOLVE :
                         LIS_KRYLOV_BREAKDOWN;
                    goto done;
                }
                if (i == 0) { ++m; ++extensions; }
            }
            if (!extensions) {
                /* A semidefinite metric can have rank smaller than nev.
                 * Preserve finite Ritz approximations and report that
                 * the requested number of pairs was not obtained. */
                rc = o->allow_semidefinite_metric ?
                     LIS_KRYLOV_NO_CONVERGENCE : LIS_KRYLOV_BREAKDOWN;
                goto done;
            }
        }
    }
    rc = LIS_KRYLOV_NO_CONVERGENCE;

done:
    if (rc != LIS_KRYLOV_SUCCESS && rc != LIS_KRYLOV_NO_CONVERGENCE) {
        stats->nreturned = 0;
        stats->nconv = 0;
    }
    release_workspace(&s);
    return rc;
}

/* For SPD B, Sylvester's inertia law converts symmetric indefinite
 * factorizations of A-shift*B into an eigenvalue count. The callback owns
 * all factorization and pivot-sign decisions; the eigensolver has no direct
 * factorization dependency. The interval is open and its endpoints must
 * not coincide with eigenvalues (i.e. no zero inertia pivots). */
int lis_krylov_count_interval(int n, void *ctx,
                              lis_krylov_inertia_callback inertia,
                              double lower, double upper, int *count)
{
    int neg_lo = 0, neg_hi = 0, zero_lo = 0, zero_hi = 0;
    if (count) *count = 0;
    if (n <= 0 || !inertia || !count ||
        !isfinite(lower) || !isfinite(upper) || !(lower < upper))
        return LIS_KRYLOV_BAD_ARGUMENT;
    if (inertia(ctx, n, lower, &neg_lo, &zero_lo) ||
        inertia(ctx, n, upper, &neg_hi, &zero_hi))
        return LIS_KRYLOV_CALLBACK_ERROR;
    if (neg_lo < 0 || neg_lo > n || neg_hi < 0 || neg_hi > n ||
        zero_lo < 0 || zero_lo > n - neg_lo ||
        zero_hi < 0 || zero_hi > n - neg_hi ||
        neg_hi < neg_lo)
        return LIS_KRYLOV_BREAKDOWN;
    if (zero_lo || zero_hi) return LIS_KRYLOV_SINGULAR_SHIFT;
    *count = neg_hi - neg_lo;
    return LIS_KRYLOV_SUCCESS;
}
